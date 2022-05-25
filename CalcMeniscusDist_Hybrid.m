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
DPH_vec = [8 7 6 5 4 3 2 1]; % order of micropillars in hybrid configuration


% Lvec should eventually be used as an input to this function during the
% genetic algorithm optimization, with the function output being wick
% thermal resistance and dryout heat flux
Lvec = [0.25 0.25 0.5 0 0 0 0 0]; % example distribution of hybrid wick, quarter with geomtery 8, quarter geometry 7, half geometry 6
LRegion = Lx.*Lvec; %absolute length per region, assuming square wick

%Find number of unit cells per region for each micropillar geometry
for i = 1:length(DPH_vec)
    d(i) = DPH_Key(DPH_vec(i),1);
    p(i) = DPH_Key(DPH_vec(i),2); %vector of p for each patch based on x(i)
    h(i) = DPH_Key(DPH_vec(i),3);
    n(i) = round(LRegion(i)/p(i)); %num unit cells per region, x
    LRegion(i) = n(i)*p(i);%scale L_region to whole number of pitch multiples
end

%Don't solve for super short segments (n < 3)
nz_DPH  = DPH_vec(find(n>3)); %nz stands for non zero, find all micropillar geometries that are used in this configuration
nz_Indices = find(n>3); %indices corresponding to nonzero segments
nz_n = n(find(n>3)); %number of unit cells in each segment, don't count ones with unit cells < 3
nz_p = p(nz_Indices); %pitches corresponding to non zero segments
nz_d = d(nz_Indices); %diameters corresponding to non zero segments
nz_h = h(nz_Indices); %heights corresponding to non zero segments

nz_LRegion = LRegion(nz_Indices); %length of each non zero region
numRegions = length(nz_DPH); % total number of non zero regions
L = sum(nz_LRegion); %total length, rescaled by non zero segments (technically should have different Lx and Ly, but we're assuming square

%create matrix of micropillar geometry properties for entire wick area
totalp(1:nz_n(1)) = nz_p(1);
totald(1:nz_n(1)) = nz_d(1);
totalh(1:nz_n(1)) = nz_h(1);
totalDPH(1:nz_n(1)) = nz_DPH(1);
totalx(1:nz_n(1)) = linspace(0,nz_LRegion(1),nz_n(1));
nsum = nz_n(1);
nx = sum(nz_n);
ny = nx;
ntot = ny*nx;

if length(nz_Indices)> 1
    for i = 2:length(nz_Indices) %create vector of all pitches
        totalp(nsum+1:nsum+nz_n(i)) = nz_p(i);
        totald(nsum+1:nsum+nz_n(i)) = nz_d(i);
        totalh(nsum+1:nsum+nz_n(i)) = nz_h(i);
        totalDPH(nsum+1:nsum+nz_n(i)) = nz_DPH(i);
        totalx(nsum+1:nsum+nz_n(i)) = linspace(totalx(nsum)+nz_p(i),totalx(nsum)+nz_LRegion(i),nz_n(i));
        nsum = nsum + nz_n(i);
    end
end

%create matrices of DPH_num and p,d,h assuming square wick geometry
for i = 1:ny
    for j = 1:nx
    pmat(i,j) = totalp(i);
    dmat(i,j) = totald(i);
    hmat(i,j) = totalh(i);
    DPHmat(i,j) = totalDPH(i);
    end
end
pmat = triu(pmat)' + triu(pmat,1);
dmat = triu(dmat)' + triu(dmat,1);
hmat = triu(hmat)' + triu(hmat,1);
DPHmat = triu(DPHmat)' + triu(DPHmat,1);

% number for (diameter, pitch, height) set in database
% DPH_num = 8; 

% % set of micropillar paramters corresponding to DPH_num
% d = DPH_Key(DPH_num,1); % diameter 
% p = DPH_Key(DPH_num,2); % pitch
% h = DPH_Key(DPH_num,3); % height
% fprintf("Micropillar geometry: diameter = %.6f m, height = %.6f m, " + ...
%     "pitch = %.6f m\n", d, h, p);

% % number of cells
% nx = floor(L/p); % #nodes in x-direction
% ny = floor(W/p); % #nodes in y-direction
% fprintf("nx = %d, & ny = %d \n",nx,ny);
% ntot = nx*ny; % total #nodes

% heat per node
qn = Q/ntot; % W/cells, no evaporation at border nodes
fprintf("Heat per node = %.3f W\n",qn);

Pinit = 0; % Pa

% initialize pressure
P = zeros(ny,nx);
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
f = @(x)MeniscusDist_Hybrid(x,nx,ny, qn,DPHmat);
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

% figure(1)
% P = diag(Poptim);
% d = linspace(0, sqrt(2)*L, length(P));
% plot(d, P,'LineWidth',2)
% xlabel('Length along diagonal')
% ylabel('P_{liq}-P_{vap} [Pa]')

DPHmat_flattened = reshape(DPHmat, [ny*nx,1]);
for i = 1:length(DPHmat_flattened)
    A(i) = polyval(p_a(:,DPHmat_flattened(i)),Poptim_flattened(i));
    R(i) = polyval(p_k(:,DPHmat_flattened(i)),Poptim_flattened(i));
end
A = reshape(A,[ny,nx]);
R = reshape(R,[ny,nx]);
% 
% 
% % cell resistance
% figure(3)
% [~, c] = contour(X, Y, R, 'ShowText','on');
% xlabel('x[mm]')
% ylabel('y[mm]')
% title('Cell Resistance [Pa.s/m]')
% c.LineWidth = 3;
% 
% figure(4)
% surf(X, Y, R);
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('Cell Resistance [Pa.s/m]');
% view(45,10)
% zlim([min(R,[],'all'),max(R,[],'all')])
% 
% % cell area
% figure(5)
% [~, c] = contour(X, Y, A, 'ShowText','on');
% xlabel('x[mm]')
% ylabel('y[mm]')
% title('Area [m^2]')
% c.LineWidth = 3;
% 
% figure(6)
% surf(X, Y, A);
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('Area [m^2]');
% view(45,10)
% zlim([min(A,[],'all'),max(A,[],'all')])
% 
[Ux,Mx,Uy,My,Umag,M] = P2U(rho,Poptim,R,A);
% 
% xcoord = linspace(0,L,nx-1)*10^3; 
% ycoord = linspace(0,W,ny)*10^3;
% [X,Y] = meshgrid(xcoord, ycoord);
% 
% % unit cell velocity in x-direction
% figure(7)
% [~, c] = contour(X, Y, Ux, 'ShowText','on');
% xlabel('x[mm]')
% ylabel('y[mm]')
% title('U_x [m/s]');
% c.LineWidth = 3;
% 
% figure(8)
% surf(X, Y, Ux);
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('U_x [m/s]');
% view(45,10)
% zlim([min(Ux,[],'all'),max(Ux,[],'all')])
% 
% % mass flow rate in x-direction
% figure(9)
% [~, c] = contour(X, Y, Mx, 'ShowText','on');
% xlabel('x[mm]')
% ylabel('y[mm]')
% title('M_x [kg/s]')
% c.LineWidth = 3;
% 
% figure(10)
% surf(X, Y, Mx);
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('M_x [kg/s]');
% view(45,10)
% zlim([min(Mx,[],'all'),max(Mx,[],'all')])
% 
% xcoord = linspace(0,L,nx)*10^3; 
% ycoord = linspace(0,W,ny-1)*10^3;
% [X,Y] = meshgrid(xcoord, ycoord);
% 
% % unit cell velocity in y-direction
% figure(11)
% [~, c] = contour(X, Y, Uy, 'ShowText','on');
% xlabel('x[mm]')
% ylabel('y[mm]')
% title('U_y [m/s]');
% c.LineWidth = 3;
% 
% figure(12)
% surf(X, Y, Uy);
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('U_y [m/s]');
% view(45,10)
% zlim([min(Uy,[],'all'),max(Uy,[],'all')])
% 
% % mass flow rate in y-direction
% figure(13)
% [~, c] = contour(X, Y, My, 'ShowText','on');
% xlabel('x[mm]')
% ylabel('y[mm]')
% title('M_y [kg/s]');
% c.LineWidth = 3;
% 
% figure(14)
% surf(X, Y, My);
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('M_y [kg/s]');
% view(45,10)
% zlim([min(My,[],'all'),max(My,[],'all')])
% 
xcoord = linspace(0,Lx,nx-1)*10^3; 
ycoord = linspace(0,Ly,ny-1)*10^3;
[X,Y] = meshgrid(xcoord, ycoord);
% 
% % magnitude of unit cell velocity
% figure(15)
% [~, c] = contour(X, Y, Umag, 'ShowText','on');
% xlabel('x[mm]')
% ylabel('y[mm]')
% title('U_{mag} [m/s]')
% c.LineWidth = 3;
% 
figure(3)
surf(X, Y, Umag);
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('U_{mag} [m/s]')
view(45,10)
zlim([min(Umag,[],'all'),max(Umag,[],'all')])
% 
% 
% % unit cell mass flow rate
% figure(17)
% [~, c] = contour(X, Y, M, 'ShowText','on');
% xlabel('x[mm]')
% ylabel('y[mm]')
% title('M [kg/s]')
% c.LineWidth = 3;
% 
% figure(18)
% surf(X, Y, M);
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('M [kg/s]')
% view(45,10)
% zlim([min(M,[],'all'),max(M,[],'all')])