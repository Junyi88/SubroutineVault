clear;

x=linspace(30,100,5001);
Lx=length(x);

y0=zeros(size(x));
dy0=zeros(size(x));
y1=zeros(size(x));

A=7.2;
Tcut=1.4;
Tpass=35;
P=0.5;

for nx=1:Lx
    
    if (x(nx) >= Tpass)
        Teff=x(nx)-Tpass;
        
        y0(nx)=A.*sinh((Teff./Tcut).^P);
        y1(nx)=A.*sinh((Teff./Tcut)).^P;
        dy0(nx)=A.*cosh((Teff./Tcut).^P).*P.*((Teff./Tcut).^(P-1))./Tcut;
    end
    
end

ft1=fit(x.',y0.','smoothingspline');
dy0f=differentiate(ft1,x.');

figure(1);
clf;
hold on;
plot(x,y0,'r-');
plot(x,y1,'b-');
ylim([0 3000])

figure(2);
clf;
hold on;
plot(x,dy0f,'r-');
plot(x,dy0,'b--');

PStrain=zeros(100,1);
Stress=zeros(100,1);
Stress(1)=1.7.*TM;
Gdot=zeros(100,1);
TotalStrain=1.8;
YM=20;
StressTrial=YM.*TotalStrain;
dt=0.1;
for n1=2:length(st)
    F=StressTrial-Stress(n1-1)-TM.*(PSTrain(n1-1));   
    if (Stress(n1-1) >= Tpass)
        
        
        Teff=Stress(n1-1)-Tpass;
        
        Gdot(n1)=A.*sinh((Teff./Tcut).^P);
        
        
        y1(nx)=A.*sinh((Teff./Tcut)).^P;
        dy0(nx)=A.*cosh((Teff./Tcut).^P).*P.*((Teff./Tcut).^(P-1))./Tcut;
    end
    
end