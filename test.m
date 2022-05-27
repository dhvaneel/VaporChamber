clc; clear; close all;

global DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all Aavg_all TFarea_all ...
       K p_k p_a Pcap_vec_15 teff2um_all TF_area2um_10D;


load '10D50H_gradedData.mat' DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all ...
      Aavg_all TFarea_all K p_k p_a Pcap_vec_15 teff2um_all TF_area2um_10D;

% Lx = 2e-3;
% Ly = 2e-3;
% DPH_vec = [8, 7, 6, 5, 4, 3, 2, 1];
% seg_vec = [1.0, 0, 0, 0, 0, 0, 0, 0];
% q = 1e6;
% plot = true;
% verbose = true;
% % Poptim = solver(Lx,Ly,DPH_vec,seg_vec,q,plot,verbose);
% q1 = dryout(Lx,Ly,DPH_vec,seg_vec);
% 
Lx = 0.5e-3;
Ly = 0.5e-3;
DPH_vec = [8, 7, 6, 5, 4, 3, 2, 1];
seg_vec = [0    0.3004    0.3553    0.3443   -0.0000   -0.0000   -0.0000   -0.0000];
%seg_vec = [1    0    0    0         0         0         0         0];
q2 = dryout(Lx,Ly,DPH_vec,seg_vec);
P = solver(Lx,Ly,DPH_vec,seg_vec,q2,false,true);
r = resistance(Lx,Ly,DPH_vec,seg_vec,P);

% Lx = 2e-3;
% Ly = 2e-3;
% DPH_vec = [8, 7, 6, 5, 4, 3, 2, 1];
% seg_vec = [0.25, 0.25, 0.5, 0, 0, 0, 0, 0];
% q3 = dryout(Lx,Ly,DPH_vec,seg_vec);
% P = solver(Lx,Ly,DPH_vec,seg_vec,q3,true,true);
