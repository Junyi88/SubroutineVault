function [v0,TauPass,TauCut]=GetV0TauP(theta,rhoP,rhoF,rhoM,c1,c2,c3,c4,kB,G,b)
    
    v0=c1.*(theta.^c2)./sqrt(rhoP.*rhoF);
    TauPass=c3.*G.*b.*sqrt(rhoP+rhoM);
    TauCut=c4.*kB.*theta.*sqrt(rhoF)./(b.^2);
end