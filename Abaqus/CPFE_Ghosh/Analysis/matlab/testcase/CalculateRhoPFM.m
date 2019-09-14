function [rhoP,rhoF,rhoM]=CalculateRhoPFM(rhoSSD,c10,kB,theta,G,b,FCCSlips,CubicSlips)

rhoP=zeros(18,1);
rhoF=zeros(18,1);
rhoM=zeros(18,1);

for n1=1:12 
    for n2=1:12
        [~,rF,rP]=...
            CosineProject(FCCSlips(n1).np,FCCSlips(n2).tpb);
        rhoF(n1)=rhoF(n1)+rF;
        rhoP(n1)=rhoP(n1)+rP;
    end
    for n2=1:6
        [~,rF,rP]=...
            CosineProject(FCCSlips(n1).np,CubicSlips(n2).tpb);
        rhoF(n1)=rhoF(n1)+rF;
        rhoP(n1)=rhoP(n1)+rP;
    end
end

%%
for n1=1:6 
    for n2=1:12
        [~,rF,rP]=...
            CosineProject(CubicSlips(n1).np,FCCSlips(n2).tpb);
        rhoF(n1+12)=rhoF(n1+12)+rF;
        rhoP(n1+12)=rhoP(n1+12)+rP;
    end
    for n2=1:6
        [~,rF,rP]=...
            CosineProject(CubicSlips(n1).np,CubicSlips(n2).tpb);
        rhoF(n1+12)=rhoF(n1+12)+rF;
        rhoP(n1+12)=rhoP(n1+12)+rP;
    end
end

%% 
rhoF=rhoSSD.*rhoF;
rhoP=rhoSSD.*rhoP;
M=c10*kB*theta./(G.*(b.^3));

rhoM=M.*sqrt(rhoP.*rhoF);

end