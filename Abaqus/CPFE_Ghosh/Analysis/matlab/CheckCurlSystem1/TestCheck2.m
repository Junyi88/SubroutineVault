clear;

load SlipSystemsAllen.mat;

Stress=[-300.0 0 0;0 0 0;0 0 0];
TauC=115.0./sqrt(6);
Gdref=1e-3;
TauCut=300.0;

MyLoc=[0 1;0 1;0 1];
%%
Tau=zeros(18,1);

for n1=1:12
   Tau(n1)=FCCSlips(n1).n.'*(Stress*FCCSlips(n1).s) ;
end
for n1=1:6
   Tau(n1+12)=CubicSlips(n1).n.'*(Stress*CubicSlips(n1).s) ;
end

%%
Gdot=zeros(18,1);

for n1=1:18
   TauEff=abs(Tau(n1))-TauC;
   if TauEff>0.0
       Gdot(n1)=Gdref*sign(Tau(n1))*sinh(TauEff./TauCut);
   end
end

%% Now Start With Points

NatCoords=[...
    -1 1 1;...
    -1 -1 1;...
    -1 1 -1;...
    -1 -1 -1;...
    1 1 1;...
    1 -1 1;...
    1 1 -1;...
    1 -1 -1;...
    ];

GaussNatCoords=NatCoords./(sqrt(3));

FullPos=[
    MyLoc(1,1),MyLoc(2,1+1),MyLoc(3,1+1);...
    MyLoc(1,1),MyLoc(2,1),MyLoc(3,1+1);...
    MyLoc(1,1),MyLoc(2,1+1),MyLoc(3,1);...
    MyLoc(1,1),MyLoc(2,1),MyLoc(3,1);...
    MyLoc(1,1+1),MyLoc(2,1+1),MyLoc(3,1+1);...
    MyLoc(1,1+1),MyLoc(2,1),MyLoc(3,1+1);...
    MyLoc(1,1+1),MyLoc(2,1+1),MyLoc(3,1);...
    MyLoc(1,1+1),MyLoc(2,1),MyLoc(3,1);...
    ];

FullPosShape=zeros(8,3);

for n1=1:3
    for n2=1:8
        for n3=1:8
            FullPosShape(n2,n1)=FullPosShape(n2,n1)+...
                FullPos(n3,n1)*0.125*...
                (1+NatCoords(n2,1)*NatCoords(n3,1))*...
                (1+NatCoords(n2,2)*NatCoords(n3,2))*...
                (1+NatCoords(n2,3)*NatCoords(n3,3));
        end
    end
end

GaussPosShape=zeros(8,3);

for n1=1:3
    for n2=1:8
        for n3=1:8
            GaussPosShape(n2,n1)=GaussPosShape(n2,n1)+...
                FullPos(n3,n1)*0.125*...
                (1+GaussNatCoords(n2,1)*NatCoords(n3,1))*...
                (1+GaussNatCoords(n2,2)*NatCoords(n3,2))*...
                (1+GaussNatCoords(n2,3)*NatCoords(n3,3));
        end
    end
end


figure(1);
clf;
hold on;
plot3(FullPos(:,1),FullPos(:,2),FullPos(:,3),'rx');
plot3(GaussPosShape(:,1),GaussPosShape(:,2),GaussPosShape(:,3),'bs');

%%


for n1=1:12
   AllSlips(n1).n=FCCSlips(n1).n;
   AllSlips(n1).s=FCCSlips(n1).s;
   AllSlips(n1).t=FCCSlips(n1).t;
end
for n1=1:6
   AllSlips(n1+12).n=CubicSlips(n1).n;
   AllSlips(n1+12).s=CubicSlips(n1).s;
   AllSlips(n1+12).t=CubicSlips(n1).t;
end
%%
myseed=23;
rng(myseed);
%[CurlOut]=GetGaussianCurl(fc,xnat8,gauss,gausscoords);
Co=rand(3,8);

Co=zeros(3,8);
Co(:,1)=1;

%%
RhoS=zeros(8,18);
RhoEN=zeros(8,18);
RhoET=zeros(8,18);

RhoTS=zeros(8,18);
RhoTEN=zeros(8,18);
RhoTET=zeros(8,18);

%%


for nSl=1:18
    fc=zeros(3,8);
    fcT=zeros(3,8);
    for nG=1:8
        x=GaussPosShape(nG,1);
        y=GaussPosShape(nG,2);
        z=GaussPosShape(nG,3);

    [U,Fp,FpT,FpN,FpTN,dFpN,dFpTN]=...
        GetAnalVals(x,y,z,AllSlips(nSl).n,Co);
    
    ASol(nG,nSl).Fp=Fp;
    ASol(nG,nSl).FpT=FpT;
    ASol(nG,nSl).FpN=FpN;
    ASol(nG,nSl).FpTN=FpTN;
    ASol(nG,nSl).dFpN=dFpN;
    ASol(nG,nSl).dFpTN=dFpTN;
    
    fc(:,nG)=FpN;
    fcT(:,nG)=FpTN;
    end
    [CurlOut]=GetGaussianCurl(fc,NatCoords,NatCoords,GaussPosShape.');
    [CurlOutT]=GetGaussianCurl(fcT,NatCoords,NatCoords,GaussPosShape.');
    Xsol(nSl).CurlOut=CurlOut;
    Xsol(nSl).CurlOutT=CurlOutT;
    
    for nG=1:8
    RhoS(nG,nSl)=fc(:,nG).'*AllSlips(nSl).s;
    RhoEN(nG,nSl)=fc(:,nG).'*AllSlips(nSl).n;
    RhoET(nG,nSl)=fc(:,nG).'*AllSlips(nSl).t;

    RhoTS(nG,nSl)=fcT(:,nG).'*AllSlips(nSl).s;
    RhoTEN(nG,nSl)=fcT(:,nG).'*AllSlips(nSl).n;
    RhoTET(nG,nSl)=fcT(:,nG).'*AllSlips(nSl).t;
    end
end


RhoXS=[RhoS;RhoTS];
RhoXEN=[RhoEN;RhoTEN];
RhoXET=[RhoET;RhoTET];

RS=zeros(8,1);
REN=zeros(8,1);
RET=zeros(8,1);

RTS=zeros(8,1);
RTEN=zeros(8,1);
RTET=zeros(8,1);

for nG=1:8
   for nS=1:18
       RS(nG)=RS(nG)+RhoS(nG,nS).^2;
       REN(nG)=REN(nG)+RhoEN(nG,nS).^2;
       RET(nG)=RET(nG)+RhoET(nG,nS).^2;
       
       RTS(nG)=RTS(nG)+RhoTS(nG,nS).^2;
       RTEN(nG)=RTEN(nG)+RhoTEN(nG,nS).^2;
       RTET(nG)=RTET(nG)+RhoTET(nG,nS).^2;
   end
end

RA=RS+REN+RET;
RTA=RTS+RTEN+RTET;
