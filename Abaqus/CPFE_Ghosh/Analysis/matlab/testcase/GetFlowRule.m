function [GammaDot,v]=GetFlowRule(rhoM,v0,tau,theta,...
    tauPass,tauCut,tauC,Q,kB,P,b)

v=zeros(18,1);
M=exp(-Q./(kB.*theta));

for n1=1:18
   if (abs(tau(n1))-tauPass(n1))>tauC(n1)
       if abs(tau(n1))>tauPass(n1)
           v(n1)=M.*v0(n1).*...
               sinh(((abs(tau(n1))-tauPass(n1))./tauCut(n1)).^P).*...
               sign(tau(n1));
       end
   end
    
end

GammaDot=rhoM.*b.*v;

end