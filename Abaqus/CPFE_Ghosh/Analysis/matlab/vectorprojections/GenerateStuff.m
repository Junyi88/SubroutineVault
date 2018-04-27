clear;

load FCCSlips;
% Secondary(n1).test=[];
XSlips=Slips;

for n1=1:12
   Secondary(n1).Edge=cross(XSlips(n1).n,XSlips(n1).s);
end

for n1=1:12
    [~,nI]=max(abs(Secondary(n1).Edge));
    Sg1=sign((Secondary(n1).Edge(nI)));
    Secondary(n1).C=[0;0;0];
    Secondary(n1).C(nI)=Sg1;
end

for n1=1:12
    [~,nI]=max(abs(Secondary(n1).Edge));
    Sg1=sign((Secondary(n1).Edge(nI)));
    Secondary(n1).SeP=XSlips(n1).n;
    Secondary(n1).SeS=XSlips(n1).n;
    Secondary(n1).SeS(nI)=-Secondary(n1).SeS(nI);
end

%%
theta=10;
R=[cosd(theta) -sind(theta) 0;sind(theta) cosd(theta) 0;0 0 1];

Stress=[1 0 0;0 0 0; 0 0 0];
StreeRot=R.'*(Stress*R);

tauPB=zeros(12,1);
tauPE=zeros(12,1);

tauCB=zeros(12,1);
tauSEP=zeros(12,1);
tauSES=zeros(12,1);


for alpha=1:12
    tauPB(alpha)=XSlips(alpha).s.'*(StreeRot*XSlips(alpha).n);
    tauPE(alpha)=Secondary(alpha).Edge.'*(StreeRot*XSlips(alpha).n);
    tauCB(alpha)=XSlips(alpha).s.'*(StreeRot*Secondary(alpha).C);
    
    tauSEP(alpha)=Secondary(alpha).Edge.'*(StreeRot*Secondary(alpha).SeP);
    tauSES(alpha)=Secondary(alpha).Edge.'*(StreeRot*Secondary(alpha).SeS);

end
Tau=[tauPB tauPE tauCB tauSEP tauSES];
Taus=[sum(tauPB) sum(tauPE) sum(tauCB) sum(tauSEP) sum(tauSES)];