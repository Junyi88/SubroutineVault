clear;

load SlipSystemsAllen.mat;

Stress=[-300.0 0 0;0 0 0;0 0 0];
TauC=115.0./sqrt(6);
Gdref=1e-3;
TauCut=300.0;

MyLoc=[-1 1;-1 1;-1 1];
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
myseed=23;
rng(myseed);
%[CurlOut]=GetGaussianCurl(fc,xnat8,gauss,gausscoords);
Co=rand(3,8);

Co=zeros(3,8);
Co(:,1)=1;
for n1=1:8
x=GaussPosShape(n1,1);
y=GaussPosShape(n1,2);
z=GaussPosShape(n1,3);

[An(n1).U,An(n1).Fp,An(n1).FpT,An(n1).FpN,...
    An(n1).FpTN,An(n1).dFpN,An(n1).dFpTN]=GetAnalVals(x,y,z,[1;1;1],Co);

end

save Test1;
fc=zeros(3,8);
for n1=1:8
    fc(:,n1)=An(n1).FpN;
end

CurlAn11=zeros(8,3);
for n1=1:8
    CurlAn11(n1,:)=An(n1).dFpN;
end

[CurlOut11]=GetGaussianCurl(fc,NatCoords,NatCoords,GaussPosShape.');

