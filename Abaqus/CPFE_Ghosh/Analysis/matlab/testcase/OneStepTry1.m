clear;
theta=273.15+700;
% Stress=[100 100 100;...
%     100 200 400;...
%     100 400 300];

Stress=[100e6 0 0;...
    0 0 0;...
    0 0 0];
ttx=90;
R=[cosd(ttx) -sind(ttx) 0; sind(ttx) cosd(ttx) 0; 0 0 1];
Stress=R.'*Stress.*R;

dStress=[1e-2 0 0;0 0 0; 0 0 0];
rhoSSD=7.5e9-8.36e6.*theta;

%% Setup Elastic

C11=1e9.*(325-0.096.*theta);
C12=1e9.*(209-0.035.*theta);
C44=1e9.*(144-0.057.*theta);

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

load SlipSystems;
[Tau,TauPE,TauSE,TauCB]=GetTaus(Stress,FCCSlips,CubicSlips);

%%
c10=25;
kB=1.38064852e-23;
G=C44;
b=2.49e-10; % meters 

[rhoP,rhoF,rhoM]=CalculateRhoPFM(rhoSSD,c10,kB,theta,G,b,FCCSlips,CubicSlips);

%%
c1=1.7e16;
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
Gx=142.2e9;
[H,tPE,tSE,tCB,Inx]=GetHForABS(TauPE,TauSE,TauCB,b,Ga111,Ga010,Gx,h,k1,k2);


%%
xi0=2.1;
thetaC=1400;
A=325;
TauCC=330e6;
rho0=5e15;
Gxx=Gx;
[rhoCSD,TauC,xi]=GetRhoCSDTauC(theta,H,xi0,thetaC,A,kB,TauCC,Gxx,rho0);
%TauC=TauC.*0-1000; %HACK

%%
Q=1.1e-20;
P=0.5;
[GammaDot,v]=GetFlowRule(rhoM,v0,Tau,theta,...
    TauPass,TauCut,TauC,Q,kB,P,b);
