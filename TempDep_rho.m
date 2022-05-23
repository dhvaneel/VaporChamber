function rho = TempDep_rho(Tv, vap)
load 'SteamProp.mat'  SteamProp LiqProp Tvap_Key Tliq_Key
if vap %interpolate vapor density
rho = interp1(Tvap_Key,SteamProp(:,2),Tv); 
else %interpolate liquid density
   rho = interp1(Tliq_Key,LiqProp(:,1),Tv);
end