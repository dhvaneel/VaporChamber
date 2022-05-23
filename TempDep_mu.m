function mu = TempDep_mu(Tv, vap)
load 'SteamProp.mat'  SteamProp LiqProp Tvap_Key Tliq_Key
if vap %interpolate vapor density
mu = interp1(Tvap_Key,SteamProp(:,5),Tv); 
else %interpolate liquid density
    mu = interp1(Tliq_Key,LiqProp(:,3),Tv);
end