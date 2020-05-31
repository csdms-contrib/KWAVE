% File names
fname_topo=['topofilteredfilled_reavis2.txt'];
fname_channelmask=['ChannelMask2000m2_reavis2.txt'];

% Grid Parameters
dx=1;                                           % Grid Spacing [m]

% Model Parameters
tend=2400;                                      % Total Simulation Time [s]
ks=8/(3600*1000);                               % Saturated Hydraulic Conductivity [m/s]
hf=0.01;                                        % Suction Head [m]
theta_init=0.1;                                 % Initial volumetric soil moisture content [-]
theta_sat=0.4;                                  % Volumetric soil moisture at saturation [-]
hydraulicroughness=0.1;                         % Manning's n [s/m^1/3]
h0=0.003;                                       % Empirical parameter in depth-dependent manning equation
epsilon=-0.33;                                  % Empirical parameter in depth-dependent manning equation
d84=0.01;                                       % D84 if using VPE

% Parameters for Interception Model (Rutter, 1972)
fractionveg=0;                                      % Fraction vegetaion cover [-]
pi=0.35;                                            % Free throughfall Coefficient (Kim et al., 2013)
Si=1.51/1000;                                       % Canopy capacity [m] (Thompson et al., 2011)
Ki=5.62e-7;                                         % Canopy Drainage Rate Coefficient [m]/[s] (Thompson et al., 2011)
gi=5.04*1000;                                       % Canopy Drainage Rate Exponent 1/[m] (Thompson et al., 2011)

% Initial Topography, Grid, Mask for Internal Boundary Conditions
[xgrid,ygrid,topo]=arcgridread(fname_topo);                     % topography (m)
topo=topo(1:end,5:end);                                         % Clip the model domain here if needed 
xgrid=xgrid(1:end,5:end);                                       % Clip the model domain here if needed 
ygrid=ygrid(1:end,5:end);                                       % Clip the model domain here if needed 
solid=zeros(size(topo));                                        % Mask for boundary conditions, 1 indicates areas outside of the model domain, 0 indicates valid grid nodes where flow will be computed           
solid(isnan(topo))=1;
topo(isnan(topo))=0;                                            % if there is no data set elevation to 0
[~,~,channelmask]=arcgridread(fname_channelmask);               % Mask identifying bedrock channels within the model domain, 1 indicates a channel location
channelmask=channelmask(1:end,5:end);
itopo=topo;
solidin=solid;
vegcover=fractionveg*ones(size(topo));

% Initial Flow Depth
[nx,ny]=size(topo);
depth=1e-5*ones(nx,ny);
idepth=depth;

% More Parameters and Initial Conditions for Green-Ampt
theta0in=theta_init*ones(nx,ny);               % Initial Volumetric Water Content [-]
thetasin=theta_sat*ones(nx,ny);              % Saturated Soil Volumetric Water Content [-]
vinfin=1e-4*ones(nx,ny);                % Depth of Water Infiltrated Prior to Start of Simulation [m]
    
% Input Rainfall Intensity From Idealized Storm
idealizedstorm_peakI15=40;                      % Peak 15-minute average rainfall intensity for input rainfall time series
rain1=dlmread('15min_50mm.txt');                % Read idealized rainstorm from file
rain1=idealizedstorm_peakI15/50*rain1;          % Scale rainfall intensity to get desired peak value
rint=60;                                        % Time interval between rainfall entries (s)
rnum1=length(rain1);                            % Length of rainfall intensity vector is needed as model input
rain2=rain1;

% Create Hydraulic Roughness Coefficient Grid
manning=hydraulicroughness*ones(nx,ny);

% Create Ks and Hf grids
ksin=ks*ones(nx,ny);
hfin=hf*ones(nx,ny);

% Reshape and Write Input Files
topoin=reshape(itopo',1,nx*ny);
ksin=reshape(ksin',1,nx*ny);
theta0in=reshape(theta0in',1,nx*ny);
thetasin=reshape(thetasin',1,nx*ny);
hfin=reshape(hfin',1,nx*ny);
vinfin=reshape(vinfin',1,nx*ny);
depthin=reshape(idepth',1,nx*ny);
channelin=reshape(channelmask',1,nx*ny);
manningin=reshape(manning',1,nx*ny);
vegcoverin=reshape(vegcover',1,nx*ny);

dlmwrite('topoin',topoin,'delimiter','\t','precision','%.3f');
dlmwrite('depthin',depthin,'delimiter','\t');
dlmwrite('solidin',solidin,'delimiter','\t');
dlmwrite('rain1',rain1,'delimiter','\t');
dlmwrite('rain2',rain2,'delimiter','\t');
dlmwrite('ksin',ksin,'delimiter','\t');
dlmwrite('theta0in',theta0in,'delimiter','\t');
dlmwrite('thetasin',thetasin,'delimiter','\t');
dlmwrite('hfin',hfin,'delimiter','\t');
dlmwrite('manningin',manningin,'delimiter','\t');
dlmwrite('vinfin',vinfin,'delimiter','\t');
dlmwrite('vegcoverin',vegcoverin,'delimiter','\t');
dlmwrite('channelin',channelin,'delimiter','\t');
dlmwrite('xgrid',xgrid,'delimiter','\t','precision','%.2f');
dlmwrite('ygrid',ygrid,'delimiter','\t','precision','%.2f');

dlmwrite('input',[nx ny dx tend epsilon h0 rnum1 rint d84 pi Si Ki gi],'delimiter','\t');

figure(1)
plot((0:rint:rint*(rnum1-1))/60,rain1*3600*1000,'color','black','LineWidth',2)
xlim([0 rint*(rnum1-1)/60])
set(gca,'FontSize',14)
ylabel('Rainfall Intensity')
xlabel('Time (min)')

% Commands to compile the C code using gcc
% gcc -o kinematic kinematic.c -lm
%./kinematic

