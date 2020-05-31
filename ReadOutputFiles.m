input=dlmread('input');
dx=input(3);
tend=input(4);
rnum1=input(7);
rint=input(8);
D50=input(10);

xgrid=dlmread('xgrid');
ygrid=dlmread('ygrid');

% Read Output Files From C Code
topo=dlmread('./topo');
vel=dlmread('./vel');
velx=dlmread('./u');
vely=dlmread('./v');
depth=dlmread('depth');
solid=dlmread('solid');
stage=dlmread('stage');
maxvel=dlmread('maxvel');
maxdepth=dlmread('maxdepth');
rain1=dlmread('rain1');

topo=topo';
vel=vel';
velx=velx';
vely=vely';
depth=depth';
solid=solid';
maxvel=maxvel';
maxdepth=maxdepth';
stage=stage';

modeltime=stage(1,:);           % Time
modeleddepth=stage(16,:);       % Time series of depth at lowest elevation

modeleddepth1=stage(12,:);       % Time series of max flow depth along model edge 1 
modeleddepth2=stage(13,:);       % Time series of max flow depth along model edge 2 
modeleddepth3=stage(14,:);       % Time series of max flow depth along model edge 3 
modeleddepth4=stage(15,:);       % Time series of max flow depth along model edge 4 

% Plot Some Results
figure(1)
subplot(2,1,1)
set(gca,'FontName','Arial','FontSize',14)
plot((0:rint:rint*(rnum1-1))/60,rain1*3600*1000,'-','LineWidth',2,'color','black')
hold on
ylim([0 max(1.1*rain1*3600*1000)])
xlim([0 max(modeltime)/60])
xlabel('Time (min)')
ylabel('Rainfall Intensity (mm/hr)')
set(gca,'FontName','Arial','FontSize',14)
subplot(2,1,2)
set(gca,'FontName','Arial','FontSize',14)
plot(stage(1,:)/60,modeleddepth4,'-','LineWidth',2,'color','black')
xlim([0 max(modeltime)/60])
xlabel('Time (min)')
ylabel('Flow Depth (m)')
set(gca,'FontName','Arial','FontSize',14)

% If points are outside the computational domain, set their values to NaN
maxvel(solid==1)=nan;
maxdepth(solid==1)=nan;

% define red-blue color map
rbcolor=dlmread('redbluecolormap.txt');

figure(2)
subplot(2,1,1)
set(gca,'FontName','Arial','FontSize',14)
surf(xgrid,ygrid,maxvel,maxvel)
colormap(rbcolor)
shading interp
axis tight
xlabel('Easting (m)')
ylabel('Northing (m)')
title('Velocity (m/s)')
view([0 90])
title('Maximum Velocity [m/s]')
set(gca,'DataAspectRatio',[1 1 1])
grid off
box on
colorbar
set(gca,'FontName','Arial','FontSize',14)
grid off
subplot(2,1,2)
set(gca,'FontName','Arial','FontSize',14)
surf(xgrid,ygrid,maxdepth,maxdepth)
colormap(rbcolor)
shading interp
axis tight
xlabel('Easting (m)')
ylabel('Northing (m)')
title('Depth (m)')
view([0 90])
title('Maximum Depth [m]')
set(gca,'DataAspectRatio',[1 1 1])
grid off
box on
colorbar
grid off
set(gca,'FontName','Arial','FontSize',14)
