function [rint] = rinterface(area,Tvap_avg,alpha)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs
    % @ area = area of interface [m^2]
    % @ Tvap_avg = average temperature of vaporization [C]
    % @ alpha = ?

    % Output
    % @ rint = interface thermal resistance [K/W]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rho_v = TempDep_rho(Tvap_avg+273,1); % [kg/m^3], density of vapor
    rho_l = TempDep_rho(Tvap_avg+273,0); % [kg/m^3], density of liquid
    vfg = (rho_v.^-1-rho_l.^-1); % [m^3/kg]
    hfg = TempDep_hfg(Tvap_avg+273); % [J/kg], latent heat of water
    R = 8.314; % [J/mol.K], universal gas constant
    M = 0.01801; % [kg/mol], molecular mass
    Pvap = ClausiusClapeyron(Tvap_avg+273); % [Pa], vapor pressure
    rint = (1./(area))*(2-alpha/(2*alpha)).*((Tvap_avg+273).*vfg./ ...
           (hfg.^2)).*sqrt(2*pi*R*(Tvap_avg+273)/M).*(1-Pvap.*vfg./ ...
           (2*hfg)).^-1;
end