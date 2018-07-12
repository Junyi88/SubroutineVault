function [MyOut]=CalculateRhoCSDSingle(Stress,Para,FCCSlips)


NSlips=length(FCCSlips);

MyOut.TauPE=zeros(12,1);
MyOut.TauSE=zeros(12,1);
MyOut.TauCB=zeros(12,1);

MyOut.tPE=zeros(12,1);
MyOut.tSE=zeros(12,1);
MyOut.tCB=zeros(12,1);

MyOut.A1=zeros(12,1);
MyOut.A2=zeros(12,1);
MyOut.A3=zeros(12,1);
MyOut.A4=zeros(12,1);
MyOut.H=zeros(12,1);
MyOut.A5=zeros(12,1);
MyOut.A6=zeros(12,1);

for nS=1:NSlips
   
    MyOut.TauPE(nS)=(FCCSlips(nS).npe.')*(Stress*FCCSlips(nS).npe);
    MyOut.TauSE(nS)=(FCCSlips(nS).nse.')*(Stress*FCCSlips(nS).nse);
    MyOut.TauCB(nS)=(FCCSlips(nS).ncb.')*(Stress*FCCSlips(nS).ncb);
    
    MyOut.tPE(nS)=MyOut.TauPE(nS)*Para.tauCo;
    MyOut.tSE(nS)=MyOut.TauSE(nS)*Para.tauCo;
    MyOut.tCB(nS)=MyOut.TauCB(nS)*Para.tauCo;
    
    MyOut.A1(nS)=Para.k1*(MyOut.tPE(nS)-Para.k2*MyOut.tSE(nS));
    MyOut.A2(nS)=(1/sqrt(3))-Para.G010/Para.G111;
    MyOut.A3(nS)=MyOut.A2(nS)+abs(MyOut.tCB(nS));
    MyOut.A4(nS)=Para.h+MyOut.A2(nS)+sqrt(MyOut.A3(nS)*Para.b/Para.B);
    MyOut.H(nS)=Para.CH*MyOut.A4(nS);
    MyOut.A5(nS)=MyOut.H(nS)/(Para.kb*Para.Temp);
    MyOut.A6(nS)=exp(-MyOut.A5(nS));
    
end
end