


function [A,Avec] = CalculateAngle(Afun,xn,yn,Us,Vs,lat,lon,ids,t2b,beached)
% DESCRIPTION
% for each drifter, finds the angle* between the wind and land where it beached. 
% *angle can also be dot product, magnitude, or whatever other parameter
% you'd like. The parameter you are calculating is specified by function "Afun"
% 
% INPUTS
% Afun: anon func for finding whatever angle equivalent. Inputs must always
%       be @(xn,yn,Us,Vs)
% xn,yn: land normal vectors from 'Data/land_normal.mat'
% Us,Vs: wind vectors from 'Data/wind_coarse.mat'
% lat,lon,ids,t2b,beached: data from 'Data/undrogued_beach.mat'
%   Afun = @(xn,yn,Us,Vs)
% note: the names from the files are different :P
% 
% OUTPUT
% A: matrix of all of the angles
% Avec: vector of all of the angles per beached drifter
% 
% Afun examples
% Afun = @(xn,yn,Us,Vs) sign(xn.*Us).*abs(Us);

%% Calculate things

% find angle-like parameter from input function
A = Afun(xn,yn,Us,Vs);

% matrix where each input is # drifters beached in each cell
[~,~,nB_out,IDcellout] = gridit(lon(beached),lat(beached),-180:180,-90:90,'ZData',ids(beached));
nB = [nB_out(:,181:end),nB_out(:,1:180)];
IDcell = [IDcellout(:,181:end),IDcellout(:,1:180)];

% indices of beachedmat where there is beaching
beachindx = find(nB~=0);

% tabulate all the angles where a drifter beached
Avec = [];
n=1;
for ib = beachindx'

    % find the row and column from the index
    [ir,ic] = ind2sub(size(nB),ib);

    % if the beached cell is not nan
    if ~isnan(A(ir,ic))
        Avec_app = A(ir,ic)*ones(nB(ib),1);
    
    % if it is nan, do some things
    else

        % rows and columns of nearby cells
        irs = ir-1:ir+1;
        ics = ic-1:ic+1;
        
        % wrap the columns around the world
        if ics(1) == 0;ics(1) = size(nB,2);end
        if ics(3)>size(nB,2);ics(3) = 1;end
        
        % matrix around the beached drifters
        As = A(irs,ics);
        
        % where are the nans
        nanss = ~isnan(As);
        
        % middle of the ocean
        if sum(nanss(:))==0
            Avec_app = ones(nB(ib),1)*NaN;
        
        % just one nearby cell
        elseif sum(nanss(:))==1
            % this_angle=anglemat(matnan);
            Avec_app = ones(nB(ib),1) * As(nanss);
        
        % if that didn't work, place using trajectory
        else
            
            % try to find the correct cell via trajectory
            % fprintf('Checking Trajectory...\n%1.0f of %1.0f cells\n\n',n,numel(beachindx))
            Avec_app = placewithtrajs(As,IDcell{ir,ic},lat,lon,t2b,ids);
            
            % if some didn't have cells, do averaging
            if any(isnan(Avec_app))

                % check if it's in the middle of land
                sumc = sum(nanss,2);
                sumr = sum(nanss,1);

                % if it's not sandwiched between land, take the angle as
                % the average of those around it
                if sumc(2)~=2 && sumr(2)~=2
                    avg_angle = mean(As(nanss));
                    Avec_app(isnan(Avec_app)) = avg_angle;
                end

            end %any(isnan(Avec_app))
        end %sum(nanss(:))==0
    end % ~isnan(A(ir,ic))

    % assign the angle
    Avec = [Avec; Avec_app]; 
    n=n+1;
end

%% try to use previous trajectories to place

function anglesout = placewithtrajs(As,IDs,lat,lon,time2beach,ids)
% go through each id in a cell and try to find the nearest cell where there
% is land data

% initialize
anglesout = nan(numel(IDs),1);

% function to find boundaries of point
boundfunc = @(pt) [floor(pt),ceil(pt)];

% loop through each id in this cell
for ii = 1:numel(IDs)

    % find the trajectory for just this data
    thisIDindx = ids == IDs(ii);
    lons = lon(thisIDindx);
    lats = lat(thisIDindx);
    tims = time2beach(thisIDindx);

    % find the last time
    tzero = find(tims==0);
    latzero = lats(tzero);
    lonzero = lons(tzero);

    % define the center cell
    loncenter = boundfunc(lonzero);
    latcenter = boundfunc(latzero);

    % find what the time to beach in center cell is
    lonsin = lons>=loncenter(1) & lons<=loncenter(2);
    latsin = lats>=latcenter(1) & lats<=latcenter(2);
    ptsin = lonsin & latsin;
    timsin = tims(ptsin);

    % find cell right before beaching
    timbef = max(timsin)+1;
    befindx = tims == timbef;
    lonbef = boundfunc(lons(befindx));
    latbef = boundfunc(lats(befindx));

    % check that there is a time before
    if any(befindx)

        % find how it shifted
        eastwest = loncenter-lonbef;
        norsouth = latcenter-latbef;

        % check the column values
        if eastwest(1)==eastwest(2)
            ic2 = eastwest(1);
            if ic2==0;ic2=360;end
            if ic2==361;ic2=1;end
        else;keyboard %idk what happens here?
        end
    
        % check the row values
        if norsouth(1)==norsouth(2)
            ir2 = norsouth(1);
        else;keyboard %idk what happens here?
        end
    
        % assign
        if abs(ir2)>1 || abs(ic2)>1
            anglesout(ii) = NaN;
        else
            anglesout(ii) = As(ir2+2,ic2+2);
        end

    else
        anglesout(ii) = NaN;
    end
end %ii = 1:numel(IDs)
end %function

end
