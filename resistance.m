function [r] = resistance(Lx,Ly,DPH_vec,seg_vec,P)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    % @ Lx = length of wick in x-direction [m]
    % @ Ly = length of wick in y-direction [m]
    % @ DPH_vec = DPH_num of each segment (8 segments)
    % @ seg_vec = length proportions of each segment (8 segments)
    % @ P = pressure profile [Pa]
    % 
    % Output
    % @ r = resistance grid [K/W]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    global DPH_Key;

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
    p_max = min(p_vec); 
    
    % num nodes in wick with largest pitch segment
    n_max = round(Lx./p_max); 
    
    % scale all segments wrt max pitch segment
    scale = p_vec/p_max;
    n_scale = round(n_seg.*scale); % effective num nodes
    
    % scaling factor grid
    n_cum = 0;
    scale_grid = zeros(n_max, n_max); 
    DPH_grid = zeros(n_max,n_max);
    p_grid = zeros(n_max,n_max);
    h_grid = zeros(n_max,n_max);
    d_grid = zeros(n_max,n_max);
    for i = 1:length(n_scale)
        n = n_scale(i);
        scale_grid(n_cum+1:n_cum+n,:) = 1/scale(i);
        DPH_grid(n_cum+1:n_cum+n,:) = DPH_vec(i);
        p_grid(n_cum+1:n_cum+n,:) = p_vec(i);
        h_grid(n_cum+1:n_cum+n,:) = h_vec(i);
        d_grid(n_cum+1:n_cum+n,:) = d_vec(i);
        n_cum = n_cum + n;
    end
    scale_grid = triu(scale_grid)' + triu(scale_grid,1);
    DPH_grid = triu(DPH_grid)' + triu(DPH_grid,1);
    p_grid = triu(p_grid)' + triu(p_grid,1);
    h_grid = triu(h_grid)' + triu(h_grid,1);
    d_grid = triu(d_grid)' + triu(d_grid,1);
        
    % wick nodes
    nx = n_max; % num nodes in x-direction
    ny = n_max; % num nodes in y-direction
    
    % resistance constants
    Tvap_avg = 100; % C
    alpha = 0.1;
    kSi = 149; % W/m-K
    kWater = 0.6; % W/-K

    teff = t_eff(nx,ny,P,DPH_grid);
    tfarea = tf_area(nx,ny,P,DPH_grid);
    
    for j = 1:ny
        for i = 1:nx
            % micropillar conduction resistance
            rpin(j,i) = h_grid(j,i)/(kSi*0.25*pi*d_grid(j,i)^2);
            % bulk water resistance
            rwater(j,i) = h_grid(j,i)/(kWater*(p_grid(j,i)^2-0.25*pi* ...
                          d_grid(j,i)^2));
            % thin film resistance
            rtf(j,i) = teff(j,i)/(kWater*tfarea(j,i));
            % liquid-vapor interface resistance
            rwaterint(j,i) = rinterface(p_grid(j,i)^2-0.25*pi* ...
                             d_grid(j,i)^2,Tvap_avg,alpha);
            % thin film interface resistance
            rint(j,i) = rinterface(tfarea(j,i),Tvap_avg, alpha);
        end
    end
    
    % equivalent resistance grid
    r = ((rpin+rtf+rint).^-1+(rwater+rwaterint).^-1).^-1;
end