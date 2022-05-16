%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D vapor chamber evporation wick modeling
% resistance-network-based approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data and global variables

clear; close all; clc;
clear global

global DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all Aavg_all TFarea_all K p_k p_a;

load '10D50H_gradedData.mat' DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all ...
    Aavg_all TFarea_all K p_k p_a;

% DPH key
% 1: 10 micron diameter, 15 micron pitch, 50 micron height
% 2: 10 micron diameter, 20 micron pitch, 50 micron height
% 3: 10 micron diameter, 25 micron pitch, 50 micron height
% 4: 10 micron diameter, 30 micron pitch, 50 micron height
% 5: 10 micron diameter, 35 micron pitch, 50 micron height
% 6: 10 micron diameter, 40 micron pitch, 50 micron height
% 7: 10 micron diameter, 45 micron pitch, 50 micron height
% 8: 10 micron diameter, 50 micron pitch, 50 micron height

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% given parameters and properties 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% water properties
rho = 997; % kg/m^3, density of water at 100C;
hfg = 2264705; % J/kg latent heat of water

% wick dimensions, assume square for simplicity
Lx = 2e-3; % m, length of wick
Ly = 2e-3; % m, width of wick

% heating parameters 
q = 1e6; % W/m^2, uniform heat flux
Q = q*Lx*Ly; % W, total heat input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hybrid segment configuration
DPH_vec = [8 7 6 5 4 3 2 1];
seg_vec = [0.25 0.25 0.5 0 0 0 0 0];
%seg_vec = [1.0 0.0 0.0 0 0 0 0 0];
assert(sum(seg_vec)==1, 'Segment length proportions incorrect!')

% non-zero hybrid segment configuration
ind = find(seg_vec > 0);
DPH_vec = DPH_vec(ind);
seg_vec = seg_vec(ind);

% segment micropillar geometry
d_vec = DPH_Key(DPH_vec, 1)'; % [m], diameter
p_vec = DPH_Key(DPH_vec, 2)'; % [m], pitch
h_vec = DPH_Key(DPH_vec, 3)'; % [m], height

% segmented region
seg_len = Lx.*seg_vec; % [m], length
n_seg = round(seg_len./p_vec); % num node

% segment with largest pitch
p_max = max(p_vec); 

% num nodes in wick with largest pitch segment
n_max = round(Lx./p_max); 

% scale all segments wrt max pitch segment
scale = p_vec/p_max;
n_scale = round(n_seg.*scale); % effective num nodes

% scaling factor grid
n_cum = 0;
scale_grid = zeros(n_max, n_max); 
DPH_grid = zeros(n_max,n_max);
for i = 1:length(n_scale)
    n = n_scale(i);
    scale_grid(n_cum+1:n_cum+n,:) = 1/scale(i);
    DPH_grid(n_cum+1:n_cum+n,:) = DPH_vec(i);
    n_cum = n_cum + n;
end
scale_grid = triu(scale_grid)' + triu(scale_grid,1);
DPH_grid = triu(DPH_grid)' + triu(DPH_grid,1);
    
% wick nodes
nx = n_max; % num nodes in x-direction
ny = n_max; % num nodes in y-direction

% heat per node
qn = Q/(n*n); % [W/node]

Pinit = 0; % [Pa]
P = zeros(n,n); % [Pa]
for j = 1:ny
    for i = 1:nx
        P(j,i) = 0.01*(Pinit-1-2*(i-1)*(j-1));
    end
end
P_flattened = reshape(P,[ny*nx,1]);

% guess 
x0 = P_flattened;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving non-linear equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

fprintf("Solving system of non-linear equations:\n");
options = optimset('Display','iter','TolFun',1e-12,'TolX',1e-12, ...
                   'MaxIter',1e6,'MaxFunEvals',1e6,'UseParallel',false);
f = @(x)MeniscusDist(x,nx,ny,qn,scale_grid,DPH_grid);
[x,fval,exitflag,output] = fsolve(f,x0,options);

toc

Poptim_flattened = x(1:(ny*nx));
Poptim = reshape(Poptim_flattened,[ny,nx]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xcoord = linspace(0,Lx,nx)*10^3;
ycoord = linspace(0,Ly,ny)*10^3; 
[X,Y] = meshgrid(xcoord, ycoord);

% pressure
figure(1)
[~, c] = contour(X, Y, Poptim, 'ShowText','on');
xlabel('x[mm]')
ylabel('y[mm]')
title('P_{liq}-P_{vap} [Pa]')
c.LineWidth = 3;

figure(2)
surf(X, Y, Poptim)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('P_{liq}-P_{vap} [Pa]');
view(45,10)
zlim([min(Poptim,[],'all'),max(Poptim,[],'all')])