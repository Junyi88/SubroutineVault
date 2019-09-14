function [Tau,TauPE,TauSE,TauCB]=GetTaus(Stress,FCCSlips,CubicSlips)

Tau=zeros(18,1);
TauPE=zeros(12,1);
TauSE=zeros(12,1);
TauCB=zeros(12,1);

for n1=1:12
   Tau(n1)=FCCSlips(n1).np.'*Stress*FCCSlips(n1).sb; 
   TauPE(n1)=FCCSlips(n1).np.'*Stress*FCCSlips(n1).se; 
   TauSE(n1)=FCCSlips(n1).ns.'*Stress*FCCSlips(n1).se; 
   TauCB(n1)=FCCSlips(n1).nc.'*Stress*FCCSlips(n1).sb; 
   
end

for n1=1:6
   Tau(n1+12)=CubicSlips(n1).np.'*Stress*CubicSlips(n1).sb; 
end


end