function [Pv] = ClausiusClapeyron(Tsat)
% load 'SteamProp.mat' Thfg_Key hfg_Key
% hfg = interp1(Thfg_Key,hfg_Key,Tsat)*10^3; %J/kg
hfg = TempDep_hfg(Tsat);
Rg = 8.3144598; %J/K/mol
Mv = 0.01801528; %kg/mol
Patm = 101325; %atmospheric pressure, pascals
Tatm = 100+273; %Tsat at Patm
Pv = Patm.*exp(-(hfg*Mv/Rg).*(Tsat.^-1-Tatm^-1));