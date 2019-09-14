clear;

c44=27;
c12=40;
c11=2.*c44+c12;

CVoigt0=[c11 c12 c12 0 0 0;...
         c12 c11 c12 0 0 0;...
         c12 c12 c11 0 0 0;...
         0 0 0 c44 0 0;...
         0 0 0 0 c44 0;...
         0 0 0 0 0 c44;...
         ];
     
%%
ThetaX=11;
ThetaY=-98;
ThetaZ=30;

Rx=[1 0 0;0 cosd(ThetaX) -sind(ThetaX);0 sind(ThetaX) cosd(ThetaX)];
Ry=[cosd(ThetaY) 0 sind(ThetaY);0 1 0; -sind(ThetaY) 0 cosd(ThetaY)];
Rz=[cosd(ThetaZ) -sind(ThetaZ) 0; sind(ThetaZ) cosd(ThetaZ) 0;0 0 1];

Q=Rz*Ry*Rz;

%%

[CFull0]=ConvertVoigt2Full(CVoigt0);
[CFull1]=RotateFull(Q,CFull0);
[CVoigtList1,PosList1]=ConvertFullToVoigtList(CFull1);

[CFull0]=ConvertVoigt2Full(CVoigt0);
[CVoigtList0,PosList0]=ConvertFullToVoigtList(CFull0);
for ni=1:6
    for nj=1:6
[CType0(ni,nj).C]=ExtractNumberVoigt(ni,nj,CVoigtList0,PosList0);
    end
end
%%
for ni=1:6
    for nj=1:6
[CType(ni,nj).C]=ExtractNumberVoigt(ni,nj,CVoigtList1,PosList1);
    end
end

C1=zeros(6,6);

for ni=1:6
    for nj=1:6
C1(ni,nj)=CType(ni,nj).C(1);
    end
end

Strain1=[1;0;0;0;0;0];
Stress0=CVoigt0*Strain1;
Stress1=C1*Strain1;