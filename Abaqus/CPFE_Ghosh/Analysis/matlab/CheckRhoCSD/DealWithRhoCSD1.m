clear;

load SlipSystemsAllen.mat FCCSlips CubicSlips;

MagStress=3000;

Para.h=0.3;
Para.G010=0.083e-3;
Para.G111=0.3e-3;
Para.b=2.49e-7;
Para.G=142.2e3;
Para.rho0=5e9;
Para.kb=1.38064852e-20;
Para.Temp=273.15;
Para.k1=0.5;
Para.k2=0.6;
Para.tauCo=Para.b/Para.G111;
Para.B=Para.G*(Para.b.^2)/(2.*pi.*Para.G111);
Para.CH=Para.G.*(Para.b.^3)./(4.*pi);

Stress=[MagStress 0 0;0 0 0;0 0 0];

%%
[MyOut]=CalculateRhoCSDSingle(Stress,Para,FCCSlips);
thetaX=18*5*pi/180;
Rx=[cos(thetaX) -sin(thetaX) 0;...
    sin(thetaX) cos(thetaX) 0;...
    0 0 1];

Stresses=linspace(-30000,30000,121).';
LStresses=length(Stresses);

TauPE=zeros(LStresses,12);
TauSE=zeros(LStresses,12);
TauCB=zeros(LStresses,12);

tPE=zeros(LStresses,12);
tSE=zeros(LStresses,12);
tCB=zeros(LStresses,12);

A1=zeros(LStresses,12);
A2=zeros(LStresses,12);
A3=zeros(LStresses,12);
A4=zeros(LStresses,12);
H=zeros(LStresses,12);
A5=zeros(LStresses,12);
A6=zeros(LStresses,12);

%% 


for nStress=1:LStresses
    Stress=Rx*[Stresses(nStress) 0 0;0 0 0;0 0 0]*Rx.';
    [MyOut]=CalculateRhoCSDSingle(Stress,Para,FCCSlips);
    
    
    TauPE(nStress,:)=MyOut.TauPE.';
    TauSE(nStress,:)=MyOut.TauSE.';
    TauCB(nStress,:)=MyOut.TauCB.';

    tPE(nStress,:)=MyOut.tPE.';
    tSE(nStress,:)=MyOut.tSE.';
    tCB(nStress,:)=MyOut.tCB.';

    A1(nStress,:)=MyOut.A1.';
    A2(nStress,:)=MyOut.A2.';
    A3(nStress,:)=MyOut.A3.';
    A4(nStress,:)=MyOut.A4.';
    H(nStress,:)=MyOut.H.';
    A5(nStress,:)=MyOut.A5.';
    A6(nStress,:)=MyOut.A6.';
    
end

%%
save Temporary1;

%%
MyLineSpec(1).s=['r-'];
MyLineSpec(2).s=['g-'];
MyLineSpec(3).s=['b-'];
MyLineSpec(4).s=['k-'];
MyLineSpec(5).s=['c-'];
MyLineSpec(6).s=['m-'];

MyLineSpec(7).s=['r--'];
MyLineSpec(8).s=['g--'];
MyLineSpec(9).s=['b--'];
MyLineSpec(10).s=['k--'];
MyLineSpec(11).s=['c--'];
MyLineSpec(12).s=['m--'];

%%
figure(1);
clf;
hold on;
for n1=1:12
   plot(Stresses,A1(:,n1),MyLineSpec(n1).s);
end

%%
figure(2);
clf;
hold on;
for n1=1:12
   plot(Stresses,A3(:,n1),MyLineSpec(n1).s);
end

%%
figure(3);
clf;
hold on;
for n1=1:12
   plot(Stresses,A4(:,n1),MyLineSpec(n1).s);
end

%%
figure(4);
clf;
hold on;
for n1=1:12
   plot(Stresses,A5(:,n1),MyLineSpec(n1).s);
end

%%
figure(5);
clf;
hold on;
for n1=1:12
   plot(Stresses,Para.rho0.*A6(:,n1),MyLineSpec(n1).s);
end