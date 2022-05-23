function hfg = TempDep_hfg(Tv)
load 'SteamProp.mat'  Thfg_Key hfg_Key
hfg = interp1(Thfg_Key,hfg_Key,Tv)*10^3; %J/kg