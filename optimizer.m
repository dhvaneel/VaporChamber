function [eqn] = optimizer(x, nx, ny, qn, scale_grid, DPH_grid)
    
    global p_k p_a;
    
    % nx/ny = number of unit cells
    % qn = heat per node
    % scale = factor to correct for changed area and cell resistance
    % DPH_num = number corresponding to look up table of specific diameter,
    % pitch, and height micropillar combination
    
    rho = 997; % kg/m3, density of water at 100C;
    hfg = 2264705; % J/kg, latent heat of water
    
    % initialize unknowns [size = ny*nx]
    
    P_flattened = x(1:ny*nx);
    P = reshape(P_flattened,[ny,nx]);
    A = zeros(ny,nx);
    R = zeros(ny,nx);
    for j = 1:ny
        for i = 1:nx
            A(j,i) = polyval(p_a(:,DPH_grid(j,i)),P(j,i))*scale_grid(j,i);
            R(j,i) = polyval(p_k(:,DPH_grid(j,i)),P(j,i))*scale_grid(j,i);
        end
    end

    % equations for unkowns [size = ny*nx]
    eqn = zeros(ny*nx, 1);
    
    % boundary conditions
    Pinit = 0; 
    
    % top row: P = 0 [size = nx]
    eqn(1:nx) = P(1,:) - Pinit;
    
    % left column [size = ny-1]: 
    eqn(nx+1:nx+ny-1) = P(2:ny,1) - Pinit;
    
    % bottom row [size = nx-2] 
    eqn(nx+ny:2*nx+ny-3) = ((A(end,1:nx-2).*R(end,2:nx-1))./(R(end,1:nx-2).*A(end,2:nx-1))).*(P(end,1:nx-2)-P(end,2:nx-1)) ...
                         + ((A(end-1,2:nx-1).*R(end,2:nx-1))./(R(end-1,2:nx-1).*A(end,2:nx-1))).*(P(end-1,2:nx-1)-P(end,2:nx-1)) ...
                         - (P(end,2:nx-1)-P(end,3:nx)) - (qn*R(end,2:nx-1))./(rho*hfg*A(end,2:nx-1));
    
    % right column [size = ny-2] 
    eqn(2*nx+ny-2:2*nx+2*ny-5) = (A(2:ny-1,end-1).*R(2:ny-1,end))./(R(2:ny-1,end-1).*A(2:ny-1,end)).*(P(2:ny-1,end-1)-P(2:ny-1,end)) ...
                               + (A(1:ny-2,end).*R(2:ny-1,end))./(R(1:ny-2,end).*A(2:ny-1,end)).*(P(1:ny-2,end)-P(2:ny-1,end)) ...
                               - (P(2:ny-1,end)-P(3:ny,end)) - (qn*R(2:ny-1,end))./(rho*hfg*A(2:ny-1,end));
    
    % interior nodes [size = (ny-2)*(nx-2)] 
    unflatten_eqn = (A(2:ny-1,1:nx-2).*R(2:ny-1,2:nx-1))./(R(2:ny-1,1:nx-2).*A(2:ny-1,2:nx-1)).*(P(2:ny-1,1:nx-2)-P(2:ny-1,2:nx-1)) ...
                  + (A(1:ny-2,2:nx-1).*R(2:ny-1,2:nx-1))./(R(1:ny-2,2:nx-1).*A(2:ny-1,2:nx-1)).*(P(1:ny-2,2:nx-1)-P(2:ny-1,2:nx-1)) ...
                  - (2*P(2:ny-1,2:nx-1)-P(2:ny-1,3:nx)-P(3:ny,2:ny-1)) - (qn*R(2:ny-1,2:nx-1))./(rho*hfg*A(2:ny-1,2:nx-1));

    eqn(2*ny+2*nx-4:nx*ny-1) = reshape(unflatten_eqn,[(ny-2)*(nx-2),1]);

    % center [size = 1] 
    eqn(ny*nx) = (A(ny,nx-1)*R(ny,nx))/(R(ny,nx-1)*A(ny,nx))*(P(ny,nx-1)-P(ny,nx)) ...
               + (A(ny-1,nx)*R(ny,nx))/(R(ny-1,nx)*A(ny,nx))*(P(ny-1,nx)-P(ny,nx)) ...
               - (qn*R(ny,nx))/(rho*hfg*A(ny,nx));
end 