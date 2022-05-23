function [tfarea] = tf_area(nx, ny, P, DPH_grid)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    % @ nx = number of nodes in x-direction
    % @ ny = number of nodes in y-direction
    % @ P = pressure profile [Pa]
    % @ DPH_grid = DPH num of each node

    % Output
    % @ tfarea = thin film are
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    global Prel_Key TF_area2um_10D
    
    tfarea = zeros(ny,nx);
    for j = 1:ny
        for i = 1:nx
            P_vec = Prel_Key(:,DPH_grid(j,i));
            tfarea_vec = TF_area2um_10D(:,DPH_grid(j,i));
            tfarea(j,i) = interp1(P_vec(2:end),tfarea_vec(2:end),P(j,i) ...
                   ,'spline');
        end
    end
end