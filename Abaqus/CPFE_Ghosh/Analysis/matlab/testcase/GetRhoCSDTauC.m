function [rhoCSD,TauC,xi]=GetRhoCSDTauC(theta,H,xi0,thetaC,A,kB,TauCC,G,rho0,b)

xi=xi0.*exp(A./(theta-thetaC));
rhoCSD=rho0.*exp(-H./(kB.*theta));
TauC(13:18,1)=TauCC;
TauC(1:12,1)=xi.*G.*b.*sqrt(rhoCSD(1:12,1));

end