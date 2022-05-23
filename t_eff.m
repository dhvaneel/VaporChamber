function [teff] = t_eff(nx, ny, P, DPH_grid)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    % @ nx = number of nodes in x-direction
    % @ ny = number of nodes in y-direction
    % @ P = pressure profile [Pa]
    % @ DPH_grid = DPH num of each node

    % Output
    % @ teff = thin film thickness [m]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global Prel_Key teff2um_all

    teff = zeros(ny,nx);
    for j = 1:ny
        for i = 1:nx
            P_vec = Prel_Key(:,DPH_grid(j,i));
            teff_vec = teff2um_all(:,DPH_grid(j,i));
            teff(j,i) = interp1(P_vec(2:end),teff_vec,P(j,i),'spline');
        end
    end

end