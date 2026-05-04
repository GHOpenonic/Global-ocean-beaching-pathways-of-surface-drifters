#!/usr/bin/env python3
"""
Watershed clustering stability analysis — three-regime version.

This script:
0. Import dependencies, define helper functions
1. Initialize Dask
2. Set params / parse CLI args
3. Load and prepare trajectory data
4. Build reusable ocean/coast-detour features once
5. Precompute trajectory-history feature matrices for each sample count
6. Run the full HDBSCAN hyperparameter grid with Dask (all regimes pooled)
7. For each regime, build a consensus watershed assignment from that regime's runs
8. Save tables and diagnostic figures per regime

Expected input:
    A parquet file containing at least these columns:
        id, time, lat, lon

Noise handling:
    HDBSCAN label -1 is retained as ungrouped, with cluster_probability = 0.

Primary outputs (per regime subdirectory):
    run_summary.csv
    cluster_count_summary.json
    consensus_cluster_df.parquet
    consensus_watershed_certainty.parquet
    consensus_watershed_match_details.parquet
    labels_matrix.npy
    probabilities_matrix.npy
    coassociation_matrix.npy
    figs/*.png
"""

"""
0. Import dependencies, define helper functions
"""

import argparse
import itertools
import json
import logging
import math
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import geopandas as gpd
import geodatasets
import shapely
import scipy.sparse
import scipy.sparse.csgraph

import hdbscan
from sklearn.cluster import AgglomerativeClustering

try:
    import cartopy.crs as ccrs
except Exception:
    ccrs = None

try:
    from dask.distributed import Client, LocalCluster, as_completed
except Exception:
    Client = None
    LocalCluster = None
    as_completed = None

try:
    from scalene import scalene_profiler
except Exception:
    scalene_profiler = None


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    force=True
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def parse_list(value, dtype=float):
    """Parse comma-separated CLI list values."""
    return [dtype(item.strip()) for item in value.split(",") if item.strip() != ""]


def read_input_table(path, columns):
    """Read the trajectory table from either parquet or CSV."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path, columns=columns)
    if suffix == ".csv":
        return pd.read_csv(path, usecols=columns)
    raise ValueError(f"Unsupported input file format for {path}. Use .parquet or .csv")


def ensure_required_columns(df, required_columns):
    """Ensure required columns are regular dataframe columns (handles index cols)."""
    df = df.copy()
    missing = [c for c in required_columns if c not in df.columns]
    if not missing:
        return df
    df = df.reset_index()
    missing = [c for c in required_columns if c not in df.columns]
    if missing:
        raise KeyError(
            f"Missing required columns after reset_index(): {missing}. "
            f"Available columns are {list(df.columns)}"
        )
    return df


def try_write_parquet(df, path):
    """Write parquet when an engine is available; otherwise skip with a warning."""
    try:
        df.to_parquet(path, index=False)
        return True
    except Exception as exc:
        logger.warning("Skipping parquet output for %s: %s", path, exc)
        return False


def normalize_lon(lon):
    """Normalize longitudes to [-180, 180)."""
    return ((lon + 180) % 360) - 180


def approximate_ocean_distance(lat1, lon1, lat2, lon2):
    """Approximate distance in degrees, accounting for longitude convergence."""
    dlat = lat1 - lat2
    dlon = abs(lon1 - lon2)
    dlon = min(dlon, 360 - dlon)
    mean_lat = np.deg2rad((lat1 + lat2) / 2)
    return np.sqrt(dlat**2 + (dlon * np.cos(mean_lat))**2)


# ---------------------------------------------------------------------------
# Ocean-path graph helpers
# ---------------------------------------------------------------------------

def nearest_main_ocean_cell(lat, lon, ocean_path_data, max_radius=40):
    """Find the nearest water grid cell in the largest connected ocean component."""
    lons = ocean_path_data["lons"]
    lats = ocean_path_data["lats"]
    water_mask = ocean_path_data["water_mask"]
    node_ids = ocean_path_data["node_ids"]
    component_labels = ocean_path_data["component_labels"]
    main_component = ocean_path_data["main_component"]

    lon = normalize_lon(lon)
    row_center = int(np.argmin(np.abs(lats - lat)))
    col_center = int(np.argmin(np.abs(lons - lon)))

    best = None
    best_distance = np.inf

    for radius in range(max_radius + 1):
        for row_idx in range(max(0, row_center - radius),
                             min(len(lats), row_center + radius + 1)):
            for col_candidate in range(col_center - radius, col_center + radius + 1):
                col_idx = col_candidate % len(lons)
                node_id = node_ids[row_idx, col_idx]
                if (
                    node_id >= 0
                    and water_mask[row_idx, col_idx]
                    and component_labels[node_id] == main_component
                ):
                    dlon = abs(lons[col_idx] - lon)
                    dlon = min(dlon, 360 - dlon)
                    distance = abs(lats[row_idx] - lat) + dlon
                    if distance < best_distance:
                        best_distance = distance
                        best = (node_id, (row_idx, col_idx))
        if best is not None:
            return best

    raise ValueError(f"No nearby main-ocean cell found for lat={lat}, lon={lon}")


def build_ocean_path_data(coast_grid_resolution):
    """Build a coarse water graph and precompute anchor-to-ocean-cell distances once."""
    logger.info("Building ocean path data")

    world = gpd.read_file(geodatasets.get_path("naturalearth.land"))
    land = world.geometry.union_all()

    lons = np.arange(-180 + coast_grid_resolution / 2, 180, coast_grid_resolution)
    lats = np.arange(-80 + coast_grid_resolution / 2, 80, coast_grid_resolution)

    lon_grid, lat_grid = np.meshgrid(lons, lats)
    points = shapely.points(lon_grid.ravel(), lat_grid.ravel())

    is_land = shapely.contains(land, points).reshape(lon_grid.shape)
    water_mask = ~is_land

    node_ids = np.full(water_mask.shape, -1, dtype=np.int64)
    node_ids[water_mask] = np.arange(water_mask.sum())

    rows, cols, data = [], [], []
    for dr, dc in [(0, 1), (1, 0), (1, 1), (1, -1)]:
        valid = water_mask & np.roll(np.roll(water_mask, -dr, axis=0), -dc, axis=1)
        if dr > 0:
            valid[-dr:, :] = False
        src = node_ids[valid]
        tgt = np.roll(np.roll(node_ids, -dr, axis=0), -dc, axis=1)[valid]
        latv = lat_grid[valid]
        step_lat = abs(dr) * coast_grid_resolution
        step_lon = abs(dc) * coast_grid_resolution * np.cos(np.deg2rad(latv))
        weight = np.sqrt(step_lat**2 + step_lon**2)
        rows.extend(src.tolist());  cols.extend(tgt.tolist());  data.extend(weight.tolist())
        rows.extend(tgt.tolist());  cols.extend(src.tolist());  data.extend(weight.tolist())

    ocean_graph = scipy.sparse.coo_matrix(
        (data, (rows, cols)),
        shape=(water_mask.sum(), water_mask.sum())
    ).tocsr()

    _, component_labels = scipy.sparse.csgraph.connected_components(ocean_graph, directed=False)
    main_component = np.bincount(component_labels).argmax()

    ocean_path_data = {
        "world": world,
        "lons": lons, "lats": lats,
        "water_mask": water_mask, "node_ids": node_ids,
        "ocean_graph": ocean_graph,
        "component_labels": component_labels,
        "main_component": main_component,
    }

    anchor_nodes, anchor_coords = [], []
    for anchor_lat in [-60, -30, 0, 30, 60]:
        for anchor_lon in [-150, -90, -30, 30, 90, 150]:
            anchor_node, (ri, ci) = nearest_main_ocean_cell(anchor_lat, anchor_lon, ocean_path_data)
            if anchor_node not in anchor_nodes:
                anchor_nodes.append(anchor_node)
                anchor_coords.append((lats[ri], lons[ci]))

    ocean_path_data["anchor_nodes"] = np.array(anchor_nodes)
    ocean_path_data["anchor_coords"] = np.array(anchor_coords)
    ocean_path_data["anchor_distances"] = scipy.sparse.csgraph.dijkstra(
        ocean_graph, directed=False, indices=ocean_path_data["anchor_nodes"]
    )

    logger.info("Ocean path data complete: %d water cells, %d anchors",
                water_mask.sum(), len(anchor_nodes))
    return ocean_path_data


def build_coast_detour_features(last_points_input, ocean_path_data):
    """Build coastline-informed detour features for final beaching points."""
    logger.info("Building coast detour features")
    anchor_distances = ocean_path_data["anchor_distances"]
    anchor_coords = ocean_path_data["anchor_coords"]
    features = []
    for lat, lon in zip(last_points_input["lat"].to_numpy(),
                        last_points_input["lon"].to_numpy()):
        ocean_node, _ = nearest_main_ocean_cell(lat, lon, ocean_path_data)
        water_dist = anchor_distances[:, ocean_node]
        straight = np.array([
            approximate_ocean_distance(lat, lon, alat, alon)
            for alat, alon in anchor_coords
        ])
        features.append(water_dist - straight)
    return np.asarray(features, dtype=np.float32)


# ---------------------------------------------------------------------------
# Trajectory feature helpers
# ---------------------------------------------------------------------------

def sample_single_trajectory(group, trajectory_sample_count):
    """Sample evenly spaced pre-beaching lat/lon points for one trajectory."""
    prebeach = group.iloc[:-1] if len(group) > 1 else group.iloc[:1]
    sample_idx = np.linspace(0, len(prebeach) - 1, trajectory_sample_count).round().astype(int)
    sampled = prebeach.iloc[sample_idx][["lat", "lon"]].to_numpy(dtype=np.float32)
    tid = group["id"].iloc[0] if "id" in group.columns else group.name
    row = {"id": tid}
    for i, (slat, slon) in enumerate(sampled):
        row[f"traj_lat_{i}"] = slat
        row[f"traj_lon_{i}"] = slon
    return pd.Series(row)


def build_trajectory_history_features(trajectory_sorted, ordered_ids, trajectory_sample_count):
    """Build an id-aligned trajectory-history feature matrix."""
    logger.info("Building trajectory-history features for sample_count=%d",
                trajectory_sample_count)
    samples = (
        trajectory_sorted
        .groupby("id", group_keys=False, sort=False)
        .apply(lambda g: sample_single_trajectory(g, trajectory_sample_count))
        .reset_index(drop=True)
    )
    feat_cols = [
        c for i in range(trajectory_sample_count)
        for c in (f"traj_lat_{i}", f"traj_lon_{i}")
    ]
    samples = samples.set_index("id").loc[ordered_ids]
    return samples[feat_cols].to_numpy(dtype=np.float32)


def build_cluster_input(
    beach_lat, beach_lon, start_lat, start_lon,
    trajectory_history_features, coast_detour_features,
    start_weight, trajectory_history_weight,
    coast_detour_weight, trajectory_sample_count
):
    """Build the numerical HDBSCAN input matrix for a single hyperparameter run."""
    th_scale = trajectory_history_weight / np.sqrt(trajectory_sample_count)
    return np.column_stack((
        beach_lat, beach_lon,
        start_weight * start_lat,
        start_weight * start_lon,
        th_scale * trajectory_history_features,
        coast_detour_weight * coast_detour_features,
    )).astype(np.float32)


# ---------------------------------------------------------------------------
# Single-experiment runner
# ---------------------------------------------------------------------------

def make_run_id(params):
    """Stable readable run id."""
    return (
        f"thw_{params['trajectory_history_weight']:.3f}_"
        f"tsc_{params['trajectory_sample_count']:03d}_"
        f"mcs_{params['min_cluster_size']:03d}_"
        f"eps_{params['cluster_selection_epsilon']:.3f}"
    ).replace(".", "p")


def run_single_experiment(
    params, ids,
    beach_lat, beach_lon, start_lat, start_lon,
    trajectory_history_features, coast_detour_features,
    start_weight, coast_detour_weight
):
    """Run one HDBSCAN experiment and return compact in-memory results."""
    run_id = make_run_id(params)
    t0 = time.perf_counter()

    cluster_input = build_cluster_input(
        beach_lat=beach_lat, beach_lon=beach_lon,
        start_lat=start_lat, start_lon=start_lon,
        trajectory_history_features=trajectory_history_features,
        coast_detour_features=coast_detour_features,
        start_weight=start_weight,
        trajectory_history_weight=params["trajectory_history_weight"],
        coast_detour_weight=coast_detour_weight,
        trajectory_sample_count=params["trajectory_sample_count"],
    )

    hdb = hdbscan.HDBSCAN(
        min_cluster_size=params["min_cluster_size"],
        cluster_selection_epsilon=params["cluster_selection_epsilon"],
        cluster_selection_method="leaf",
        core_dist_n_jobs=1,
    )
    hdb.fit(cluster_input)

    labels = hdb.labels_.astype(np.int32)
    probabilities = hdb.probabilities_.astype(np.float32)
    probabilities[labels == -1] = 0.0

    n_watersheds = int(np.sum(np.unique(labels) != -1))
    n_noise = int(np.sum(labels == -1))
    n_clustered = int(np.sum(labels != -1))

    elapsed_minutes = (time.perf_counter() - t0) / 60
    summary = {
        "run_id": run_id,
        "trajectory_history_weight": params["trajectory_history_weight"],
        "trajectory_sample_count": params["trajectory_sample_count"],
        "min_cluster_size": params["min_cluster_size"],
        "cluster_selection_epsilon": params["cluster_selection_epsilon"],
        "n_watersheds": n_watersheds,
        "n_clustered": n_clustered,
        "n_noise": n_noise,
        "mean_cluster_probability": (
            float(probabilities[labels != -1].mean()) if n_clustered > 0 else np.nan
        ),
        "elapsed_minutes": elapsed_minutes,
    }
    return {"summary": summary, "labels": labels, "probabilities": probabilities}


# ---------------------------------------------------------------------------
# Consensus helpers
# ---------------------------------------------------------------------------

def relabel_clusters_by_size(labels):
    """Renumber clusters by descending size, preserving -1 as noise."""
    labels = np.asarray(labels, dtype=np.int32)
    non_noise = labels[labels != -1]
    if len(non_noise) == 0:
        return labels.copy()
    counts = pd.Series(non_noise).value_counts().sort_values(ascending=False)
    mapping = {old: new for new, old in enumerate(counts.index)}
    mapping[-1] = -1
    return np.array([mapping[l] for l in labels], dtype=np.int32)


def build_coassociation_matrix(labels_matrix):
    """Entry (i,j) = fraction of runs where i and j share a non-noise cluster."""
    logger.info("Building coassociation matrix")
    n_traj, n_runs = labels_matrix.shape
    coassoc = np.zeros((n_traj, n_traj), dtype=np.float32)
    for r in range(n_runs):
        rl = labels_matrix[:, r]
        for cl in np.unique(rl):
            if cl == -1:
                continue
            members = np.flatnonzero(rl == cl)
            coassoc[np.ix_(members, members)] += 1.0
    coassoc /= n_runs
    np.fill_diagonal(coassoc, 1.0)
    return coassoc


def build_consensus_cluster_df(
    ids, run_summary, labels_matrix, probabilities_matrix, outdir,
    noise_grouped_threshold=0.5,
    min_consensus_cluster_size=30,
    min_consensus_existence_probability=0.67,
    jaccard_threshold=0.5,
):
    """
    Build consensus cluster assignment with noise filtering and pruning.

    Steps:
      1. Trajectories grouped in < noise_grouped_threshold fraction of runs
         become consensus noise.
      2. Agglomerative clustering on the coassociation of the remaining subset.
      3. Prune any consensus cluster with fewer than min_consensus_cluster_size
         members or whose watershed existence probability (Jaccard-based) is
         below min_consensus_existence_probability back into noise.
      4. Relabel by descending size.
    """
    logger.info("Building consensus cluster assignment")

    n_traj, n_runs = labels_matrix.shape
    grouped_fraction = (labels_matrix != -1).mean(axis=1).astype(np.float32)

    clustered_mask = grouped_fraction >= noise_grouped_threshold
    clustered_indices = np.flatnonzero(clustered_mask)
    n_clustered = int(clustered_mask.sum())

    logger.info(
        "Consensus noise filter: %d clustered, %d noise (threshold %.2f)",
        n_clustered, n_traj - n_clustered, noise_grouped_threshold,
    )

    mean_n_ws = float(run_summary["n_watersheds"].mean())
    floor_mean_n_ws = max(1, int(math.floor(mean_n_ws)))

    coassociation_full = build_coassociation_matrix(labels_matrix)
    coassoc_sub = coassociation_full[np.ix_(clustered_indices, clustered_indices)]
    dist_sub = 1.0 - coassoc_sub
    np.fill_diagonal(dist_sub, 0.0)

    n_consensus = min(floor_mean_n_ws, n_clustered)

    agg = AgglomerativeClustering(n_clusters=n_consensus, metric="precomputed",
                                  linkage="average")
    subset_labels = agg.fit_predict(dist_sub)

    consensus_raw = np.full(n_traj, -1, dtype=np.int32)
    consensus_raw[clustered_indices] = subset_labels

    # --- Quick Jaccard existence probability per candidate cluster ----------
    # (needed for the pruning step; we recompute the full version later)
    def _quick_existence(candidate_label, raw_labels, lm, jt):
        members = np.flatnonzero(raw_labels == candidate_label)
        member_set = set(members.tolist())
        n_present = 0
        for ri in range(lm.shape[1]):
            rl = lm[:, ri]
            best_j = 0.0
            for cl in np.unique(rl):
                if cl == -1:
                    continue
                run_members = set(np.flatnonzero(rl == cl).tolist())
                inter = len(member_set & run_members)
                union = len(member_set | run_members)
                if union > 0:
                    best_j = max(best_j, inter / union)
            if best_j >= jt:
                n_present += 1
        return n_present / max(1, lm.shape[1])

    # Prune small or low-existence clusters
    for cand in np.unique(consensus_raw):
        if cand == -1:
            continue
        members = np.flatnonzero(consensus_raw == cand)
        if len(members) < min_consensus_cluster_size:
            logger.info("Pruning consensus cluster %d: only %d members (< %d)",
                        cand, len(members), min_consensus_cluster_size)
            consensus_raw[members] = -1
            continue
        exist_prob = _quick_existence(cand, consensus_raw, labels_matrix,
                                      jaccard_threshold)
        if exist_prob < min_consensus_existence_probability:
            logger.info(
                "Pruning consensus cluster %d: existence_prob=%.3f (< %.3f)",
                cand, exist_prob, min_consensus_existence_probability,
            )
            consensus_raw[members] = -1

    consensus_labels = relabel_clusters_by_size(consensus_raw)

    mean_cluster_prob = probabilities_matrix.mean(axis=1).astype(np.float32)

    assignment_certainty = np.zeros(n_traj, dtype=np.float32)
    for cl in np.unique(consensus_labels):
        if cl == -1:
            continue
        members = np.flatnonzero(consensus_labels == cl)
        if len(members) > 1:
            block = coassociation_full[np.ix_(members, members)]
            assignment_certainty[members] = (
                (block.sum(axis=1) - 1.0) / (len(members) - 1)
            ).astype(np.float32)
        else:
            assignment_certainty[members] = grouped_fraction[members]

    cluster_df = pd.DataFrame({
        "trajectory_index": np.arange(n_traj, dtype=np.int32),
        "id": ids,
        "HDBSCAN": consensus_labels,
        "fraction_runs_grouped": grouped_fraction,
        "mean_cluster_probability": mean_cluster_prob,
        "assignment_certainty": assignment_certainty,
    })

    n_final = int(np.sum(np.unique(consensus_labels) != -1))
    cluster_count_summary = {
        "mean_n_watersheds": mean_n_ws,
        "floor_mean_n_watersheds": floor_mean_n_ws,
        "n_consensus_clusters_before_prune": n_consensus,
        "n_consensus_clusters": n_final,
        "median_n_watersheds": float(run_summary["n_watersheds"].median()),
        "std_n_watersheds": float(run_summary["n_watersheds"].std(ddof=0)),
        "n_runs": int(n_runs),
        "n_trajectories": int(n_traj),
        "n_consensus_noise": int((consensus_labels == -1).sum()),
        "noise_grouped_threshold": noise_grouped_threshold,
        "min_consensus_cluster_size": min_consensus_cluster_size,
        "min_consensus_existence_probability": min_consensus_existence_probability,
    }

    outdir = Path(outdir)
    try_write_parquet(cluster_df.sort_values(["HDBSCAN", "id"]),
                      outdir / "consensus_cluster_df.parquet")
    cluster_df.sort_values(["HDBSCAN", "id"]).to_csv(
        outdir / "consensus_cluster_df.csv", index=False)
    np.save(outdir / "labels_matrix.npy", labels_matrix)
    np.save(outdir / "probabilities_matrix.npy", probabilities_matrix)
    np.save(outdir / "coassociation_matrix.npy", coassociation_full)
    with open(outdir / "cluster_count_summary.json", "w") as f:
        json.dump(cluster_count_summary, f, indent=2)

    logger.info(
        "Consensus complete: mean_n_ws=%.3f, n_consensus=%d, n_noise=%d",
        mean_n_ws, n_final, int((consensus_labels == -1).sum()),
    )
    return cluster_df, coassociation_full, cluster_count_summary


def build_consensus_match_details(cluster_df, run_summary, labels_matrix,
                                  jaccard_threshold):
    """Match each consensus watershed to the best-fitting cluster in every run."""
    logger.info("Matching consensus watersheds back to individual runs")

    consensus_labels = cluster_df["HDBSCAN"].to_numpy(dtype=np.int32)
    consensus_sizes = {
        k: v for k, v in
        cluster_df.groupby("HDBSCAN")["id"].size().to_dict().items()
        if k != -1
    }

    records = []
    for run_idx, run_row in enumerate(run_summary.itertuples(index=False)):
        tmp = pd.DataFrame({
            "consensus_HDBSCAN": consensus_labels,
            "run_HDBSCAN": labels_matrix[:, run_idx],
        })
        tmp = tmp[tmp["run_HDBSCAN"] != -1]

        if len(tmp) == 0:
            for cl, cs in consensus_sizes.items():
                records.append({
                    "run_id": run_row.run_id,
                    "consensus_HDBSCAN": cl,
                    "best_run_HDBSCAN": -1,
                    "consensus_cluster_size": cs,
                    "run_cluster_size": 0,
                    "overlap_size": 0,
                    "member_recall": 0.0,
                    "member_precision": 0.0,
                    "jaccard": 0.0,
                    "present_at_threshold": 0,
                })
            continue

        overlap = (
            tmp.groupby(["consensus_HDBSCAN", "run_HDBSCAN"])
            .size().rename("overlap_size").reset_index()
        )
        run_sizes = tmp.groupby("run_HDBSCAN").size().rename("run_cluster_size").reset_index()
        overlap = overlap.merge(run_sizes, on="run_HDBSCAN", how="left")
        overlap["consensus_cluster_size"] = overlap["consensus_HDBSCAN"].map(consensus_sizes)
        # skip rows for consensus noise that leaked in
        overlap = overlap.dropna(subset=["consensus_cluster_size"])
        overlap["member_recall"] = overlap["overlap_size"] / overlap["consensus_cluster_size"]
        overlap["member_precision"] = overlap["overlap_size"] / overlap["run_cluster_size"]
        overlap["jaccard"] = overlap["overlap_size"] / (
            overlap["consensus_cluster_size"] + overlap["run_cluster_size"]
            - overlap["overlap_size"]
        )

        best = (
            overlap
            .sort_values(["consensus_HDBSCAN", "jaccard", "member_recall",
                          "member_precision"],
                         ascending=[True, False, False, False])
            .drop_duplicates(subset=["consensus_HDBSCAN"])
            .set_index("consensus_HDBSCAN")
        )

        for cl, cs in consensus_sizes.items():
            if cl in best.index:
                m = best.loc[cl]
                j = float(m["jaccard"])
                records.append({
                    "run_id": run_row.run_id,
                    "consensus_HDBSCAN": int(cl),
                    "best_run_HDBSCAN": int(m["run_HDBSCAN"]),
                    "consensus_cluster_size": int(cs),
                    "run_cluster_size": int(m["run_cluster_size"]),
                    "overlap_size": int(m["overlap_size"]),
                    "member_recall": float(m["member_recall"]),
                    "member_precision": float(m["member_precision"]),
                    "jaccard": j,
                    "present_at_threshold": int(j >= jaccard_threshold),
                })
            else:
                records.append({
                    "run_id": run_row.run_id,
                    "consensus_HDBSCAN": int(cl),
                    "best_run_HDBSCAN": -1,
                    "consensus_cluster_size": int(cs),
                    "run_cluster_size": 0,
                    "overlap_size": 0,
                    "member_recall": 0.0,
                    "member_precision": 0.0,
                    "jaccard": 0.0,
                    "present_at_threshold": 0,
                })

    return pd.DataFrame.from_records(records)


def build_watershed_certainty_summary(cluster_df, coassociation, match_details,
                                      outdir):
    """Summarize reproducibility of each consensus watershed."""
    logger.info("Building watershed certainty summary")
    records = []
    for cl, group in cluster_df.groupby("HDBSCAN"):
        if cl == -1:
            continue
        members = group["trajectory_index"].to_numpy(dtype=np.int32)
        cs = len(members)
        if cs > 1:
            block = coassociation[np.ix_(members, members)]
            tri = block[np.triu_indices(cs, k=1)]
            mean_wc = float(tri.mean())
            min_wc = float(tri.min())
        else:
            mean_wc = float(group["fraction_runs_grouped"].iloc[0])
            min_wc = mean_wc
        cm = match_details[match_details["consensus_HDBSCAN"] == cl]
        records.append({   "HDBSCAN": int(cl),
            "n_trajectories": int(cs),
            "mean_assignment_certainty": float(group["assignment_certainty"].mean()),
            "mean_fraction_runs_grouped": float(group["fraction_runs_grouped"].mean()),
            "mean_cluster_probability": float(group["mean_cluster_probability"].mean()),
            "mean_within_cluster_coassociation": mean_wc,
            "min_within_cluster_coassociation": min_wc,
            "mean_best_match_jaccard": float(cm["jaccard"].mean()),
            "mean_best_match_member_recall": float(cm["member_recall"].mean()),
            "mean_best_match_member_precision": float(cm["member_precision"].mean()),
            "watershed_existence_probability": float(cm["present_at_threshold"].mean()),
        })

    watershed_summary = pd.DataFrame.from_records(records).sort_values("HDBSCAN").reset_index(drop=True)
    outdir = Path(outdir)
    try_write_parquet(watershed_summary, outdir / "consensus_watershed_certainty.parquet")
    watershed_summary.to_csv(outdir / "consensus_watershed_certainty.csv", index=False)
    try_write_parquet(match_details, outdir / "consensus_watershed_match_details.parquet")
    match_details.to_csv(outdir / "consensus_watershed_match_details.csv", index=False)
    return watershed_summary


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def stride_sample_frame(df, max_points):
    """Downsample large plotting tables with a simple stride."""
    if max_points is None or len(df) <= max_points:
        return df
    step = max(1, int(math.ceil(len(df) / max_points)))
    return df.iloc[::step]


def add_land_polygons(ax, world):
    """Add Natural Earth land polygons to a Cartopy axis."""
    ax.add_geometries(
        world.geometry, crs=ccrs.PlateCarree(),
        facecolor="gainsboro", edgecolor="none", zorder=0,
    )


def make_consensus_map_figures(cluster_df, watershed_summary, last_points_input,
                               trajectory_sorted, world, figdir):
    """Produce map-based figures for the consensus watershed assignment."""
    if ccrs is None:
        logger.warning("Cartopy unavailable; skipping consensus map figures")
        return
    if trajectory_sorted is None or world is None:
        logger.warning("Trajectory data or land polygons missing; skipping maps")
        return

    cluster_assignments = cluster_df.set_index("id")

    last_pts = last_points_input.copy()
    last_pts["HDBSCAN"] = last_pts["id"].map(cluster_assignments["HDBSCAN"])
    last_pts["fraction_runs_grouped"] = last_pts["id"].map(
        cluster_assignments["fraction_runs_grouped"])

    trajectory_sorted = trajectory_sorted.copy()
    trajectory_sorted["HDBSCAN"] = trajectory_sorted["id"].map(cluster_assignments["HDBSCAN"])
    trajectory_sorted["fraction_runs_grouped"] = trajectory_sorted["id"].map(
        cluster_assignments["fraction_runs_grouped"])
    trajectory_groups = trajectory_sorted.groupby("HDBSCAN", sort=False)

    noise_mask = last_pts["HDBSCAN"] == -1
    cluster_labels = sorted(
        [c for c in watershed_summary["HDBSCAN"].astype(int).tolist() if c != -1]
    )

    cmap = plt.get_cmap("tab20", max(len(cluster_labels), 1))
    cluster_colors = {cid: cmap(idx) for idx, cid in enumerate(cluster_labels)}
    ws_lookup = watershed_summary.set_index("HDBSCAN")

    # --- Beaching clusters with centroid labels ---
    fig = plt.figure(figsize=(22, 12))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    add_land_polygons(ax, world)
    ax.coastlines()

    noise_pts = last_pts.loc[noise_mask]
    if len(noise_pts) > 0:
        ax.scatter(noise_pts["lon"], noise_pts["lat"], marker="x", color="black",
                   s=20, label="Noise", transform=ccrs.PlateCarree())

    clustered_pts = last_pts.loc[~noise_mask]
    for cid in cluster_labels:
        cpts = clustered_pts.loc[clustered_pts["HDBSCAN"] == cid]
        if len(cpts) == 0:
            continue
        ax.scatter(cpts["lon"], cpts["lat"], marker="o", color=cluster_colors[cid],
                   alpha=0.7, s=20, transform=ccrs.PlateCarree())
        ax.text(cpts["lon"].mean(), cpts["lat"].mean(), str(int(cid)),
                fontsize=18, ha="center", va="center", color="black",
                transform=ccrs.PlateCarree())

    plt.title("Consensus watershed beaching clusters with centroid labels")
    fig.savefig(figdir / "consensus_beaching_clusters_with_centroids.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    # --- Trajectories by consensus watershed ---
    fig = plt.figure(figsize=(22, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -70, 70], crs=ccrs.PlateCarree())
    ax.set_ylim(-70, 70)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax.set_yticks([-70, -35, 0, 35, 70], crs=ccrs.PlateCarree())

    for cid in cluster_labels:
        if cid not in trajectory_groups.groups:
            continue
        cd = trajectory_groups.get_group(cid)
        cd = stride_sample_frame(cd, max_points=150000)
        if len(cd) == 0:
            continue
        ax.scatter(cd["lon"], cd["lat"], marker="o", color=cluster_colors[cid],
                   alpha=0.7, s=0.001, transform=ccrs.PlateCarree())

    ax.coastlines(color="black")
    add_land_polygons(ax, world)
    ax.set_aspect("auto")
    plt.title("Trajectories by consensus watershed")
    fig.savefig(figdir / "consensus_trajectories_by_cluster.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    # --- Subpanel per cluster ---
    n_clusters = len(cluster_labels)
    ncols = min(4, max(1, n_clusters))
    nrows = int(math.ceil(n_clusters / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows),
                             subplot_kw={"projection": ccrs.PlateCarree()})
    axes = np.atleast_1d(axes).ravel()

    for ax, cid in zip(axes, cluster_labels):
        if cid not in trajectory_groups.groups:
            ax.set_visible(False)
            continue
        cd = trajectory_groups.get_group(cid)
        cd_plot = stride_sample_frame(cd, max_points=80000)
        grouped = cd.groupby("id", sort=False)
        starts = grouped.first()
        ends = grouped.last()

        ax.set_extent([-180, 180, -70, 70], crs=ccrs.PlateCarree())
        ax.set_ylim(-70, 70)
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        ax.set_yticks([-70, -35, 0, 35, 70], crs=ccrs.PlateCarree())

        if len(cd_plot) > 0:
            ax.scatter(cd_plot["lon"], cd_plot["lat"], marker="o",
                       color=cluster_colors[cid], alpha=0.15, s=1,
                       transform=ccrs.PlateCarree())
        if len(starts) > 0:
            ax.scatter(starts["lon"], starts["lat"], color="black", s=15,
                       transform=ccrs.PlateCarree())
        if len(ends) > 0:
            ax.scatter(ends["lon"], ends["lat"], color="red", s=15,
                       transform=ccrs.PlateCarree())

        ax.coastlines(color="black", linewidth=0.7)
        add_land_polygons(ax, world)
        ax.set_aspect("auto")

        if cid in ws_lookup.index:
            ep = ws_lookup.loc[cid, "watershed_existence_probability"]
            title = f"Cluster {cid} | n={len(starts)} | p={ep:.2f}"
        else:
            title = f"Cluster {cid} | n={len(starts)}"
        ax.set_title(title, fontsize=12)

    for ax in axes[n_clusters:]:
        ax.set_visible(False)

    fig.suptitle("Consensus watershed trajectories, starts, and beaching endpoints",
                 fontsize=18)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(figdir / "consensus_cluster_trajectories_subpanels.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_figures(run_summary, cluster_df, watershed_summary, outdir,
                 last_points_input=None, trajectory_sorted=None, world=None):
    """Produce diagnostic figures."""
    figdir = Path(outdir) / "figs"
    figdir.mkdir(parents=True, exist_ok=True)
    logger.info("Producing diagnostic figures in %s", figdir)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(cluster_df.loc[cluster_df["HDBSCAN"] != -1, "assignment_certainty"],
            bins=np.linspace(0, 1, 21), edgecolor="black")
    ax.set_xlabel("Consensus assignment certainty")
    ax.set_ylabel("Trajectory count")
    ax.set_title("Trajectory certainty in the consensus watershed assignment")
    fig.savefig(figdir / "trajectory_assignment_certainty_histogram.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(watershed_summary["HDBSCAN"].astype(str),
           watershed_summary["watershed_existence_probability"])
    ax.set_xlabel("Consensus watershed")
    ax.set_ylabel("Watershed existence probability")
    ax.set_title("Probability that each consensus watershed reappears across the sweep")
    ax.set_ylim(0, 1)
    fig.savefig(figdir / "consensus_watershed_existence_probability.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(run_summary["n_watersheds"], bins=30, edgecolor="black")
    ax.set_xlabel("Number of watersheds")
    ax.set_ylabel("Experiment count")
    ax.set_title("Number of HDBSCAN watersheds across hyperparameter grid")
    fig.savefig(figdir / "n_watersheds_histogram.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    for param in ["trajectory_history_weight", "trajectory_sample_count",
                  "min_cluster_size", "cluster_selection_epsilon"]:
        grouped = (
            run_summary.groupby(param, as_index=False)
            .agg(mean_n_watersheds=("n_watersheds", "mean"),
                 std_n_watersheds=("n_watersheds", "std"),
                 mean_cluster_probability=("mean_cluster_probability", "mean"))
            .sort_values(param)
        )
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.errorbar(grouped[param], grouped["mean_n_watersheds"],
                    yerr=grouped["std_n_watersheds"], marker="o", capsize=4)
        ax.set_xlabel(param)
        ax.set_ylabel("Mean number of watersheds")
        ax.set_title(f"Watershed count sensitivity to {param}")
        fig.savefig(figdir / f"n_watersheds_by_{param}.png",
                    dpi=200, bbox_inches="tight")
        plt.close(fig)

    if last_points_input is not None and trajectory_sorted is not None and world is not None:
        make_consensus_map_figures(
            cluster_df=cluster_df,
            watershed_summary=watershed_summary,
            last_points_input=last_points_input,
            trajectory_sorted=trajectory_sorted,
            world=world,
            figdir=figdir,
        )


# ---------------------------------------------------------------------------
# Three-regime consensus pipeline
# ---------------------------------------------------------------------------

def run_regime_consensus(
    regime_name, regime_run_ids, all_summaries, all_labels, all_probabilities,
    ids, outdir, noise_grouped_threshold, min_consensus_cluster_size,
    min_consensus_existence_probability, jaccard_threshold
):
    """Build consensus outputs for a single regime subset of runs."""
    regime_dir = Path(outdir) / regime_name
    regime_dir.mkdir(parents=True, exist_ok=True)

    mask = all_summaries["run_id"].isin(regime_run_ids)
    regime_summary = all_summaries.loc[mask].reset_index(drop=True)

    if len(regime_summary) == 0:
        logger.warning("Regime %s has no matching runs; skipping.", regime_name)
        return None, None, None, None

    ordered = regime_summary["run_id"].tolist()
    labels_matrix = np.column_stack(
        [all_labels[rid] for rid in ordered]).astype(np.int32)
    probs_matrix = np.column_stack(
        [all_probabilities[rid] for rid in ordered]).astype(np.float32)

    regime_summary.to_csv(regime_dir / "run_summary.csv", index=False)

    cluster_df, coassociation, cluster_count_summary = build_consensus_cluster_df(
        ids=ids,
        run_summary=regime_summary,
        labels_matrix=labels_matrix,
        probabilities_matrix=probs_matrix,
        outdir=regime_dir,
        noise_grouped_threshold=noise_grouped_threshold,
        min_consensus_cluster_size=min_consensus_cluster_size,
        min_consensus_existence_probability=min_consensus_existence_probability,
        jaccard_threshold=jaccard_threshold,
    )

    match_details = build_consensus_match_details(
        cluster_df=cluster_df,
        run_summary=regime_summary,
        labels_matrix=labels_matrix,
        jaccard_threshold=jaccard_threshold,
    )

    watershed_summary = build_watershed_certainty_summary(
        cluster_df=cluster_df,
        coassociation=coassociation,
        match_details=match_details,
        outdir=regime_dir,
    )

    return regime_dir, cluster_df, regime_summary, watershed_summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run watershed clustering stability analysis (three-regime)."
    )

    parser.add_argument("--beach-parquet", required=True,
                        help="Input beach trajectory file (.parquet or .csv).")
    parser.add_argument("--outdir", required=True, help="Output directory.")

    parser.add_argument("--id-column", default="id")
    parser.add_argument("--time-column", default="time")
    parser.add_argument("--lat-column", default="lat")
    parser.add_argument("--lon-column", default="lon")

    # Shared across regimes
    parser.add_argument("--trajectory-history-weights", default="0.2,0.1,0.3")
    parser.add_argument("--trajectory-sample-counts", default="10,50,100")

    # Regime 1:
    parser.add_argument("--regime1-min-cluster-sizes", default="35,36,37,38,34,33,32")
    parser.add_argument("--regime1-cluster-selection-epsilons",
                        default="7.5,8,8.5,9,7,6.5,6")

    # Regime 2:
    parser.add_argument("--regime2-min-cluster-sizes", default="40,39,38,37,41,42,43")
    parser.add_argument("--regime2-cluster-selection-epsilons",
                        default="10,10.5,11,11.5,9.5,9,8.5")

    # Regime 3: 
    parser.add_argument("--regime3-min-cluster-sizes", default="45,44,43,42,46,47,48")
    parser.add_argument("--regime3-cluster-selection-epsilons",
                        default="10,10.5,11,11.5,9.5,9,8.5")

    parser.add_argument("--start-weight", type=float, default=0.0)
    parser.add_argument("--coast-detour-weight", type=float, default=0.1)
    parser.add_argument("--coast-grid-resolution", type=float, default=0.5)
    parser.add_argument("--watershed-existence-jaccard-threshold", type=float,
                        default=0.5)
    parser.add_argument("--noise-grouped-threshold", type=float, default=0.5)
    parser.add_argument("--min-consensus-cluster-size", type=int, default=30)
    parser.add_argument("--min-consensus-existence-probability", type=float,
                        default=0.67)

    parser.add_argument("--n-workers", type=int, default=None)
    parser.add_argument("--threads-per-worker", type=int, default=1)
    parser.add_argument("--memory-safety-factor", type=float, default=0.90)

    parser.add_argument("--scalene-profile", action="store_true")
    parser.add_argument("--no-figures", action="store_true")

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    t0 = time.perf_counter()

    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    slurm_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    slurm_mem_mb = int(os.environ.get("SLURM_MEM_PER_NODE", "0"))

    logger.info("Output directory: %s", outdir)

    if args.scalene_profile and scalene_profiler is not None:
        scalene_profiler.start()

    # --- Dask init ---
    if args.n_workers is None:
        n_workers = max(1, slurm_cpus // max(1, args.threads_per_worker))
    else:
        n_workers = args.n_workers

    if slurm_mem_mb > 0:
        mem_gb = slurm_mem_mb / 1024
        worker_mem = f"{args.memory_safety_factor * mem_gb / n_workers:.1f}GB"
    else:
        worker_mem = "auto"

    use_dask = Client is not None and LocalCluster is not None and n_workers > 1
    cluster = None
    client = None

    if use_dask:
        logger.info("Initializing Dask: %d workers, %s memory each",
                    n_workers, worker_mem)
        cluster = LocalCluster(
            n_workers=n_workers,
            threads_per_worker=args.threads_per_worker,
            memory_limit=worker_mem,
            dashboard_address=None,
            local_directory="/tmp",
        )
        client = Client(cluster)
        logger.info("%s", client)
    else:
        logger.warning("Running experiments serially (no Dask)")

    # --- Build regime parameter grids ---
    trajectory_history_weights = parse_list(args.trajectory_history_weights, float)
    trajectory_sample_counts = parse_list(args.trajectory_sample_counts, int)

    regimes = {
        "regime1": {
            "min_cluster_sizes": parse_list(args.regime1_min_cluster_sizes, int),
            "cluster_selection_epsilons": parse_list(
                args.regime1_cluster_selection_epsilons, float),
        },
        "regime2": {
            "min_cluster_sizes": parse_list(args.regime2_min_cluster_sizes, int),
            "cluster_selection_epsilons": parse_list(
                args.regime2_cluster_selection_epsilons, float),
        },
        "regime3": {
            "min_cluster_sizes": parse_list(args.regime3_min_cluster_sizes, int),
            "cluster_selection_epsilons": parse_list(
                args.regime3_cluster_selection_epsilons, float),
        },
    }

    # Build all param records, tagging each with its regime
    all_param_records = []
    regime_run_id_map = {name: [] for name in regimes}

    for regime_name, regime_vals in regimes.items():
        for thw, tsc, mcs, eps in itertools.product(
            trajectory_history_weights,
            trajectory_sample_counts,
            regime_vals["min_cluster_sizes"],
            regime_vals["cluster_selection_epsilons"],
        ):
            params = {
                "trajectory_history_weight": float(thw),
                "trajectory_sample_count": int(tsc),
                "min_cluster_size": int(mcs),
                "cluster_selection_epsilon": float(eps),
                "regime": regime_name,
            }
            params["run_id"] = make_run_id(params)
            all_param_records.append(params)
            regime_run_id_map[regime_name].append(params["run_id"])

    # Deduplicate (in case regime boundaries share a param combo)
    seen = set()
    deduped = []
    for p in all_param_records:
        if p["run_id"] not in seen:
            seen.add(p["run_id"])
            deduped.append(p)
    all_param_records = deduped

    logger.info("Total unique experiments across all regimes: %d",
                len(all_param_records))

    with open(outdir / "hyperparameter_grid.json", "w") as f:
        json.dump(all_param_records, f, indent=2)

    # --- Load trajectory data ---
    logger.info("Loading beach parquet: %s", args.beach_parquet)

    col_map = {
        args.id_column: "id", args.time_column: "time",
        args.lat_column: "lat", args.lon_column: "lon",
    }
    beach = read_input_table(
        args.beach_parquet,
        columns=list(col_map.keys()),
    ).rename(columns=col_map)
    beach = ensure_required_columns(beach, ["id", "time", "lat", "lon"])

    logger.info("Loaded: %d rows, %.2f GB",
                len(beach), beach.memory_usage(deep=True).sum() / 1e9)

    beach["lon"] = normalize_lon(beach["lon"].to_numpy())
    trajectory_sorted = beach.sort_values(["id", "time"]).reset_index(drop=True)

    last_points_input = (
        trajectory_sorted
        .drop_duplicates(subset="id", keep="last")[["id", "lat", "lon"]]
        .reset_index(drop=True)
    )
    trajectory_start = (
        trajectory_sorted
        .drop_duplicates(subset="id", keep="first")[["id", "lat", "lon"]]
        .rename(columns={"lat": "start_lat", "lon": "start_lon"})
        .set_index("id").loc[last_points_input["id"]].reset_index()
    )

    ids = last_points_input["id"].to_numpy()
    beach_lat = last_points_input["lat"].to_numpy(dtype=np.float32).reshape(-1, 1)
    beach_lon = last_points_input["lon"].to_numpy(dtype=np.float32).reshape(-1, 1)
    start_lat = trajectory_start["start_lat"].to_numpy(dtype=np.float32).reshape(-1, 1)
    start_lon = trajectory_start["start_lon"].to_numpy(dtype=np.float32).reshape(-1, 1)

    logger.info("Number of trajectories: %d", len(ids))

    # --- Ocean path data ---
    ocean_path_data = build_ocean_path_data(args.coast_grid_resolution)
    coast_detour_features = build_coast_detour_features(
        last_points_input, ocean_path_data)

    # --- Trajectory history features ---
    history_features_by_sc = {}
    for sc in sorted(set(trajectory_sample_counts)):
        history_features_by_sc[sc] = build_trajectory_history_features(
            trajectory_sorted=trajectory_sorted,
            ordered_ids=ids,
            trajectory_sample_count=sc,
        )

    del beach

    # --- Run all experiments ---
    futures = []
    if use_dask:
        logger.info("Scattering reusable arrays to Dask workers")
        ids_f = client.scatter(ids, broadcast=True)
        blat_f = client.scatter(beach_lat, broadcast=True)
        blon_f = client.scatter(beach_lon, broadcast=True)
        slat_f = client.scatter(start_lat, broadcast=True)
        slon_f = client.scatter(start_lon, broadcast=True)
        cdf_f = client.scatter(coast_detour_features, broadcast=True)
        hf_futures = {
            sc: client.scatter(mat, broadcast=True)
            for sc, mat in history_features_by_sc.items()
        }

        logger.info("Submitting %d experiments", len(all_param_records))
        for params in all_param_records:
            fut = client.submit(
                run_single_experiment,
                params, ids_f, blat_f, blon_f, slat_f, slon_f,
                hf_futures[params["trajectory_sample_count"]],
                cdf_f, args.start_weight, args.coast_detour_weight,
                pure=False,
            )
            futures.append(fut)

    summaries = []
    labels_by_run = {}
    probs_by_run = {}
    completed = 0

    if use_dask:
        iterator = ((f.result(), len(futures)) for f in as_completed(futures))
    else:
        iterator = (
            (run_single_experiment(
                p, ids, beach_lat, beach_lon, start_lat, start_lon,
                history_features_by_sc[p["trajectory_sample_count"]],
                coast_detour_features, args.start_weight, args.coast_detour_weight,
            ), len(all_param_records))
            for p in all_param_records
        )

    for result, total in iterator:
        s = result["summary"]
        rid = s["run_id"]
        summaries.append(s)
        labels_by_run[rid] = result["labels"]
        probs_by_run[rid] = result["probabilities"]
        completed += 1
        if completed % 25 == 0 or completed == total:
            logger.info("Completed %d/%d; latest=%s watersheds=%d",
                        completed, total, rid, s["n_watersheds"])

    all_summaries = pd.DataFrame(summaries).sort_values([
        "trajectory_history_weight", "trajectory_sample_count",
        "min_cluster_size", "cluster_selection_epsilon",
    ]).reset_index(drop=True)
    all_summaries.to_csv(outdir / "run_summary_all.csv", index=False)

    # --- Per-regime consensus ---
    logger.info("Building per-regime consensus clustering")

    regime_results = {}
    for regime_name, run_ids in regime_run_id_map.items():
        # Only keep run_ids that actually completed (in case of dedup)
        valid_ids = [r for r in run_ids if r in labels_by_run]
        logger.info("Regime %s: %d runs", regime_name, len(valid_ids))

        regime_dir, cluster_df, regime_summary, watershed_summary = run_regime_consensus(
            regime_name=regime_name,
            regime_run_ids=valid_ids,
            all_summaries=all_summaries,
            all_labels=labels_by_run,
            all_probabilities=probs_by_run,
            ids=ids,
            outdir=outdir,
            noise_grouped_threshold=args.noise_grouped_threshold,
            min_consensus_cluster_size=args.min_consensus_cluster_size,
            min_consensus_existence_probability=args.min_consensus_existence_probability,
            jaccard_threshold=args.watershed_existence_jaccard_threshold,
        )
        regime_results[regime_name] = {
            "dir": regime_dir,
            "cluster_df": cluster_df,
            "run_summary": regime_summary,
            "watershed_summary": watershed_summary,
        }

    # --- Figures ---
    if not args.no_figures:
        trajectory_sorted_for_figs = None
        if ccrs is not None:
            logger.info("Reloading trajectory data for figures")
            trajectory_sorted_for_figs = read_input_table(
                args.beach_parquet, columns=list(col_map.keys()),
            ).rename(columns=col_map)
            trajectory_sorted_for_figs = ensure_required_columns(
                trajectory_sorted_for_figs, ["id", "time", "lat", "lon"])
            trajectory_sorted_for_figs["lon"] = normalize_lon(
                trajectory_sorted_for_figs["lon"].to_numpy())
            trajectory_sorted_for_figs = trajectory_sorted_for_figs.sort_values(
                ["id", "time"]).reset_index(drop=True)

        for regime_name, res in regime_results.items():
            if res["cluster_df"] is None:
                    continue
            make_figures(
                run_summary=res["run_summary"],
                cluster_df=res["cluster_df"],
                watershed_summary=res["watershed_summary"],
                outdir=res["dir"],
                last_points_input=last_points_input,
                trajectory_sorted=trajectory_sorted_for_figs,
                world=ocean_path_data["world"],
            )

    # --- Cleanup ---
    logger.info("Total time elapsed: %.3f minutes", (time.perf_counter() - t0) / 60)

    if args.scalene_profile and scalene_profiler is not None:
        scalene_profiler.stop()

    if client is not None:
        client.close()
    if cluster is not None:
        cluster.close()


if __name__ == "__main__":
    main()