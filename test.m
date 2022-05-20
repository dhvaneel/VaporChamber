clc; clear; close all;

global DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all Aavg_all TFarea_all K p_k p_a Pcap_vec_15;

load '10D50H_gradedData.mat' DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all ...
        Aavg_all TFarea_all K p_k p_a Pcap_vec_15;

% Lx = 2e-3;
% Ly = 2e-3;
% DPH_vec = [8, 7, 6, 5, 4, 3, 2, 1];
% seg_vec = [1.0, 0, 0, 0, 0, 0, 0, 0];
% q1 = dryout(Lx,Ly,DPH_vec,seg_vec);
% 
% Lx = 2e-3;
% Ly = 2e-3;
% DPH_vec = [8, 7, 6, 5, 4, 3, 2, 1];
% seg_vec = [0, 0, 1.0, 0, 0, 0, 0, 0];
% q2 = dryout(Lx,Ly,DPH_vec,seg_vec);

Lx = 2e-3;
Ly = 2e-3;
DPH_vec = [8, 7, 6, 5, 4, 3, 2, 1];
seg_vec = [0.25, 0.25, 0.5, 0, 0, 0, 0, 0];
q3 = dryout(Lx,Ly,DPH_vec,seg_vec);