clear;

p=0.5;
TauPass=300;
Tau=linspace(-500,500,1000);
TauCut=100;
TauEff=zeros(size(Tau));
TauEff2=zeros(size(Tau));
for n1=1:length(TauEff)
   TauEff(n1)=abs(Tau(n1))-TauPass;
   if TauEff(n1)<=0
       TauEff(n1)=0;
   end
   TauEff2(n1)=abs(Tau(n1))-TauPass;
   
end

rho1=sinh((TauEff./TauCut).^p).*sign(Tau)./TauEff;
rho2=sinh((TauEff./TauCut).^p).*sign(Tau)./TauEff2;

figure(1);
clf;
hold on;
plot(Tau,rho1,'r-');
plot(Tau,rho2,'b-');

figure(2);
clf;
hold on;
plot(Tau,sinh((TauEff./TauCut).^p),'r-');
plot(Tau,1./TauEff2,'b-');


figure(3);
clf;
hold on;
plot(Tau,(TauEff2./TauCut),'r-');