function [f] = fitness(seg_vec)

    global DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all Aavg_all TFarea_all ...
       K p_k p_a Pcap_vec_15 teff2um_all TF_area2um_10D;
    
    seg_vec = round(abs(seg_vec),4);
    Lx = 5e-4;
    Ly = 5e-4;
    DPH_vec = [8, 7, 6, 5, 4, 3, 2, 1];

%     disp(seg_vec)
%     try
%         f(1) = -dryout(Lx,Ly,DPH_vec,seg_vec);
%         P = solver(Lx,Ly,DPH_vec,seg_vec,-f(1),false,false);
%         f(2) = resistance(Lx,Ly,DPH_vec,seg_vec,P); 
%     catch
%         f(1) = -1e10;
%         f(2) = 1e10;
%     end
    f(1) = -dryout(Lx,Ly,DPH_vec,seg_vec);
    P = solver(Lx,Ly,DPH_vec,seg_vec,-f(1),false,false);
    f(2) = resistance(Lx,Ly,DPH_vec,seg_vec,P);
end
