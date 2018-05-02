function [Tau,TauPE,TauSE,TauCB]=GetTausAllen(Stress,FCCSlips,CubicSlips)

Tau=zeros(18,1);
TauPE=zeros(12,1);
TauSE=zeros(12,1);
TauCB=zeros(12,1);

for n1=1:12
   Tau(n1)=FCCSlips(n1).n.'*Stress*FCCSlips(n1).s; 
   TauPE(n1)=FCCSlips(n1).npe.'*Stress*FCCSlips(n1).spe; 
   TauSE(n1)=FCCSlips(n1).nse.'*Stress*FCCSlips(n1).sse; 
   TauCB(n1)=FCCSlips(n1).ncb.'*Stress*FCCSlips(n1).scb; 
   
end

for n1=1:6
   Tau(n1+12)=CubicSlips(n1).n.'*Stress*CubicSlips(n1).s; 
end


end