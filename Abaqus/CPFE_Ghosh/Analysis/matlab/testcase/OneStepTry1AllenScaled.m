clear;
theta=273.15+25;
% Stress=[100 100 100;...
%     100 200 400;...
%     100 400 300];

Stress=[0 0 0;...
    0 0 0;...
    0 0 800e6].*(1e-6);
ttx=0;
R=[cosd(ttx) -sind(ttx) 0; sind(ttx) cosd(ttx) 0; 0 0 1];
Stress=R.'*Stress.*R;

dStress=[1e-2 0 0;0 0 0; 0 0 0].*(1e-6);
rhoSSD=(7.5e9-(8.36e6).*theta).*(1e-12);

%% Setup Elastic

C11=1e9.*(325-0.096.*theta).*(1e-6);
C12=1e9.*(209-0.035.*theta).*(1e-6);
C44=1e9.*(144-0.057.*theta).*(1e-6);

CVoigt=zeros(6,6);

CVoigt(1,1)=C11;
CVoigt(2,2)=C11;
CVoigt(3,3)=C11;

CVoigt(1,2)=C12;
CVoigt(2,1)=C12;
CVoigt(1,3)=C12;
CVoigt(3,1)=C12;
CVoigt(2,3)=C12;
CVoigt(3,2)=C12;

CVoigt(4,4)=C44;
CVoigt(5,5)=C44;
CVoigt(6,6)=C44;

L=zeros(3,3,3,3);

ToVoigtUMAT=[1 4 5;
    4 2 6;
    5 6 3];

for ni=1:3
    for nj=1:3
        for nk=1:3
            for nl=1:3
                L(ni,nj,nk,nl)=CVoigt(ToVoigtUMAT(ni,nj),ToVoigtUMAT(nk,nl));
            end
        end
    end
end

%%

load SlipSystemsAllen;
[Tau,TauPE,TauSE,TauCB]=GetTausAllen(Stress,FCCSlips,CubicSlips);

%%
c10=25;
kB=(1.38064852e-23).*(1e12);
G=C44;
b=(2.49e-10).*(1e6); % meters 

[rhoP,rhoF,rhoM]=CalculateRhoPFMAllen(rhoSSD,c10,kB,theta,G,b,FCCSlips,CubicSlips);

%%
c1=1.7e16.*(1e-6);
c2=-3.77;
c3=4.0;
c4=100;

[v0,TauPass,TauCut]=GetV0TauP(theta,rhoP,rhoF,rhoM,c1,c2,c3,c4,kB,G,b);

%%
Ga111=0.083;
Ga010=0.3;
h=0.3;
k1=0.5;
k2=0.2;
Gx=142.2e9.*(1e-6);
[H,tPE,tSE,tCB,Inx]=GetHForABS(TauPE,TauSE,TauCB,b,Ga111,Ga010,Gx,h,k1,k2);


%%
xi0=2.1;
thetaC=1400;
A=325;
TauCC=(330e6).*(1e-6);
rho0=(5e15).*(1e-12);
Gxx=Gx;
[rhoCSD,TauC,xi]=GetRhoCSDTauC(theta,H,xi0,thetaC,A,kB,TauCC,Gxx,rho0,b);
%TauC=TauC.*0-1000; %HACK

%%
Q=(1.1e-20).*(1e12);
P=0.5;
[GammaDot,v]=GetFlowRule(rhoM,v0,Tau,theta,...
    TauPass,TauCut,TauC,Q,kB,P,b);

%% GET PROPERTIES

Pr=zeros(30,1);
Temp=300;

Pr(1)=c10.*kB.*Temp./(G.*(b.^3));

Pr(2)=c3.*G.*b;
Pr(3)=c4.*kB.*Temp./(b.^2);
Pr(4)=c1.*(Temp.^c2);

Pr(5)=b./Ga111;
Pr(6)=G.*(b.^3)./(4.*pi);
Pr(7)=G.*(b.^2)./(2.*pi.*Ga111);
Pr(8)=xi0.*exp(A./(Temp-thetaC));
Pr(9)=TauCC;
Pr(10)=G.*b.*b.*b./(4.*pi);
Pr(11)=h;
Pr(12)=k1;
Pr(13)=k2;
Pr(14)=(1./sqrt(3))-Ga010./Ga111;
B=G.*b.*b./(2.*pi.*Ga111);
Pr(15)=b./B; 
Pr(16)=rho0;
Pr(17)=kB.*Temp;

Pr(18)=exp(-Q./(kB.*Temp)); 
Pr(19)=P;
Pr(20)=b;

c5=1e-3;
c6=1e-4;
c7=10;
c8=10;
c9=0.3;
c10=25;

Pr(21)=c5./b; 
Poisson=C12./(2.*(C12+C44));
Pr(22)=c6.*(sqrt(3).*G)./(16.*(1-Poisson));
Pr(23)=c7;
Pr(24)=0.0; 
Pr(25)=c9;
Pr(26)=1.0;

Pr(27)=C11; 
Pr(28)=C12;
Pr(29)=C44;
Pr(30)=rhoSSD;

FileName='PROPSSCALED.csv';
fstr='%.15e';

dlmwrite(FileName,Pr,'precision',fstr);

save AllenResult1Scaled;

%% Pr2
Pr2=zeros(30+49,1);
Pr2(1)=1;
Pr2(2)=0;
Pr2(3)=0;
Pr2(4)=0;
Pr2(5)=1;
Pr2(6)=0;
Pr2(7)=0;
Pr2(8)=0;
Pr2(9)=1;

Pr2(10:27)=rhoSSD;

Pr2(28:30)=C11;
Pr2(31:33)=C44;

Pr2([34 35 39])=C12;

Pr2(50:79)=Pr;
FileName='PROPS2SCALED.csv';
fstr='%.15e';

dlmwrite(FileName,Pr2,'precision',fstr);