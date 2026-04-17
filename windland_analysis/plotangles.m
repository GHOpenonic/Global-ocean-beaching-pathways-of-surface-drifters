
%% GENERATE DATA
clc;clear
%% Load data

% load the normal data
load('Data/land_normal.mat')

% load the wind data
load('Data/wind_coarse.mat')

% load parts of the undrogued data
load('Data/undrogued_beach_wind.mat','beached','lat','lon','ids','time2beach')

%% Angle
% calculate the angle between the wind and land
Afun = @(x,y,u,v) acos((x.*u + y.*v)./ sqrt(u.^2+v.^2));
[A,Avec] = CalculateAngle(Afun,xnorm,ynorm,Usmall,Vsmall,lat,lon,ids,time2beach,beached);

% bootstrap the data
[angle.x,angle.sp,angle.y,angle.ma,angle.mc,angle.sa,angle.sc] = myboot(A,Avec);

% some drifters are not within a cell that has land normal data, this is
% the number of drifters which are used in analysis. 
ncaught = sum(~isnan(Avec));

%% magnitude

% calculate the dot product between the wind and the land
Afun = @(xn,yn,Us,Vs) xn.*Us+yn.*Vs;
[A,Avec] = CalculateAngle(Afun,xnorm,ynorm,Usmall,Vsmall,lat,lon,ids,time2beach,beached);

% bootstrap data
[dot.x,dot.sp,dot.y,dot.ma,dot.mc,dot.sa,dot.sc] = myboot(A,Avec);

%% save data
save('plotdata.mat','angle','dot','ncaught')

%% create figure
figure(3)

% plot angle resutls
subplot(121);hold on

% plots specs
errorspecs = {'LineStyle','none','LineWidth',.5,'Color','k'};
pltspecs = {'.-','LineWidth',1,'MarkerSize',20,'Color','m'}; %#D37FA8

% right plot: original y value so I dont have to propagate error
yyaxis('right')
% errorbar(xplt_angle,y_angle,sp_angle,errorspecs{:})
fill([xplt_angle,flip(xplt_angle)],[y_angle+sp_angle(:,1);flip(y_angle-sp_angle(:,2))],'m','FaceAlpha',.2,'LineStyle','none')
% plot(xplt_angle,y_angle+sp_angle(:,1),'m--')
% plot(xplt_angle,y_angle-sp_angle(:,2),'m--')

% left plot: actual y ticks
dx = (xplt_angle(2)-xplt_angle(1));
yyaxis('left')
yplt = y_angle/sum(y_angle)/dx;
plot(xplt_angle,yplt,pltspecs{:})

% set yticks so everything syncs
ylims = [0,.9];
yyaxis('left')
set(gca,'YColor','k')
ylim(ylims)
yyaxis('right')
set(gca,'YTick',[],'YColor','k')
ylim(ylims*max(y_angle)/max(yplt))

% xticks 
set(gca,'XTick',0:pi/4:pi,'XTickLabel',{'0','\frac{\pi}{4}','\frac{\pi}{2}','\frac{3\pi}{4}','\pi'})
hellpilabs(0:pi/4:pi,[1,0,0]);

% labels
yyaxis('left')
ylabel('normalized histogram')
xlabel('Coastal wind angle $\theta_w$')

%plot magnitude results
subplot(122);hold on

% separate data
xpos = xplt_dot(xplt_dot>0);
ypos = y_dot(xplt_dot>0);
spos = sp_dot(xplt_dot>0,:);
xneg = abs(xplt_dot(xplt_dot<0));
yneg = y_dot(xplt_dot<0);
sneg = sp_dot(xplt_dot<0,:);

% Label axes
xlabel('$|V_w\cdot \hat{n}|$ (m/s)');

% colors
awcol = '#183E61';
tocol = '#BA6462';
awerrorcol = '#5B5D79';
toerrorcol = '#E99973';

% plot specs 
pltspecs = {'.-','LineWidth',1,'MarkerSize',20,'Color'};
errorspecs = {'LineStyle','none','LineWidth',1,'Color'};

% right plot: original y value so I dont have to propagate error
yyaxis('right')
plot(xpos,ypos,pltspecs{:},tocol)
plot(xneg,yneg,pltspecs{:},awcol)
fill([xpos,flip(xpos)],[ypos+spos(:,1);flip(ypos-spos(:,2))],'r','FaceAlpha',.2,'LineStyle','none','Marker','none')
fill([xneg,flip(xneg)],[yneg+sneg(:,1);flip(yneg-sneg(:,2))],'b','FaceAlpha',.2,'LineStyle','none','Marker','none')


% left plot: correct y values
yyaxis('left')
dx = (xplt_dot(2)-xplt_dot(1));
ypos2 = ypos/(sum(ypos)+sum(yneg))/dx;
yneg2 = yneg/(sum(ypos)+sum(yneg))/dx;
plot(xpos,ypos2,pltspecs{:},tocol)
plot(xneg,yneg2,pltspecs{:},awcol)

% set yticks so everything syncs
ylims = [0,.5];
yyaxis('left')
set(gca,'YColor','k')
ylim(ylims)
yyaxis('right')
set(gca,'YTick',[],'YColor','k')
ylim(ylims*max([ypos;yneg])/max([ypos2;yneg2]))

% legend
leg=legend('offshore', 'onshore','Interpreter','latex','Location','northwest');
title(leg,'Wind direction')


%% function
function count = myhist(xvec,A)

count = zeros(length(xvec)-1,1);

for i=1:length(xvec)-1
    if i==1
        count(i) = sum(A>=xvec(i) & A<=xvec(i+1));
    end
    count(i) = sum(A>xvec(i) & A<=xvec(i+1));
end

end

%% function
function [xplt,sp,y,ma,mc,sa,sc] = myboot(A,Avec)
% data
AA = A(~isnan(A));
AC = Avec(~isnan(Avec));

AA(abs(AA)>4)=[];
AC(abs(AC)>4)=[];

% Parameters
nboot = 100; 
nbins = 20; 

% edges of the data
edges = linspace(min(AA),max(AA),nbins+1);

% bootstrap angledata
Aboot = zeros(nboot,nbins);AbootC = Aboot;
for i = 1:nboot

    % for the angle data
    Abootsample = AA(randi(length(AA), length(AA), 1)); 
    Aboot(i,:) = myhist(edges,Abootsample);

    % for the angle counts
    AbootCsample = AC(randi(length(AC), length(AC), 1)); 
    AbootC(i,:) = myhist(edges,AbootCsample);
    
end
xplt = (edges(1:end-1)+edges(2:end))/2;

% Calculate mean and standard deviation for each bin
ma = mean(Aboot)';
for i=1:size(Aboot,2)
    sortaboot = sort(Aboot(:,i));
    sa(i,:) = [sortaboot(round(nboot*.05)),sortaboot(round(nboot*.95))]-ma(i);
end
% sa = std(Aboot);

mc = mean(AbootC)';
% sc = std(AbootC);
for i=1:size(AbootC,2)
    sortaboot = sort(AbootC(:,i));
    sc(i,:) = [sortaboot(round(nboot*.05)),sortaboot(round(nboot*.95))]-mc(i);
end

% find the covariance
sac = sc*0;
for i=1:nbins
    covout = cov(Aboot(:,i),AbootC(:,i));
    sac(i,:) = covout(1,1)*[1,1];
end

% propagate uncertaintiy
sp = mc./ma .* sqrt(abs( (sa./ma).^2 + (sc./mc).^2 - 2*sac./ma./mc ));

% normalized
y = mc./ma;

end
