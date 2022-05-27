function [q_dryout] = dryout(Lx,Ly,DPH_vec,seg_vec)

    global DPH_Key Pgrad_Key Prel_Key Uavg_all CA_all Aavg_all ...
           TFarea_all K p_k p_a Pcap_vec_15;
    
    % solver settings
    plot = false;
    verbose = false;
    
    % non-zero hybrid segment configuration
    p_vec = DPH_Key(DPH_vec, 2)'; % [m], pitch
    p_max = max(p_vec);
    scale = p_vec./p_max;

    seg_len = Lx.*seg_vec; % [m], length
    n_seg = floor(seg_len./p_vec); % num node
    n_scale = floor(n_seg.*scale); % effective num nodes

    ind = find(n_scale > 1);
    DPH_vec = DPH_vec(ind);
    seg_vec = seg_vec(ind);

    % segment micropillar geometry
    d_vec = DPH_Key(DPH_vec, 1)'; % [m], diameter
    p_vec = DPH_Key(DPH_vec, 2)'; % [m], pitch
    h_vec = DPH_Key(DPH_vec, 3)'; % [m], height

    % segmented region
    seg_len = Lx.*seg_vec; % [m], length
    n_seg = floor(seg_len./p_vec); % num node
    seg_len = n_seg.*p_vec; % updated
    
    % segment with largest pitch
    p_max = max(p_vec); 
    
    % scale all segments wrt max pitch segment
    scale = p_vec/p_max;
    n_scale = floor(n_seg.*scale); % effective num nodes
    n_max = sum(n_scale);
    
    q_dryout = 1e10; % random initialization
    n_cum = 0;
    for i = 1:length(n_scale)
        n = n_scale(i);
        DPH_num = DPH_vec(i);
    
        % max capillary pressure
        P_dryout = Pcap_vec_15(DPH_num);
    
        % newton raphson method
        eps = 1e-2;
        q = 1e6;
        h = 1e3;
        f1 = 1e6;
        for iter = 1:10
            Poptim_1 = solver(Lx,Ly,DPH_vec,seg_vec,q,plot,verbose);
            Pseg_1 = abs(triu(Poptim_1)); % assuming Lx = Ly symmetry
            Pseg_1 = Pseg_1(n_cum+1:n_cum+n,:);
            P1 = -max(max(Pseg_1));
            f1 = (P1 - P_dryout)/P_dryout;
   
            Poptim_2 = solver(Lx,Ly,DPH_vec,seg_vec,q+h,plot,verbose);
            Pseg_2 = abs(triu(Poptim_2)); % assuming Lx = Ly symmetry
            Pseg_2 = Pseg_2(n_cum+1:n_cum+n,:);
            P2 = -max(max(Pseg_2));
            f2 = (P2 - P_dryout)/P_dryout;
    
            df = (f2 - f1)/h;
            step = f1 / df;
            q = q - step;
            if abs(f1) < eps
                break
            end
        end
        if q < q_dryout
            q_dryout = q;
        end
        n_cum = n_cum + n;
    end
    
    %disp(q_dryout);

end