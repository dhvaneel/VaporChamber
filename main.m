%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main script for multiobjective GA optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

global DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all Aavg_all TFarea_all ...
       K p_k p_a Pcap_vec_15 teff2um_all TF_area2um_10D;


load '10D50H_gradedData.mat' DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all ...
      Aavg_all TFarea_all K p_k p_a Pcap_vec_15 teff2um_all TF_area2um_10D;

num = 8;
lb = zeros(num,1);
ub = ones(num,1);

Aeq = ones(1,num);
beq = 1;

pop_size = 10; % 200
gen_lim = 20; % 200
num_gen = 30; % 300

tic

opts = optimoptions('gamultiobj','Populationsize',pop_size, ...
       'Generations',num_gen,'MutationFcn',{@mutationadaptfeasible},...
       'CrossoverFcn',@crossoverintermediate,'StallGenLimit',gen_lim, ...
       'PlotFcn',{@gaplotpareto,@gaplotscorediversity},'Display','iter');

%rng default % for reproducibility
[x,fval,exitflag,output,population,scores] = gamultiobj(@fitness,num, ...
                                             [],[],Aeq,beq,lb,ub,[],opts);
toc