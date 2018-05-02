function [StrainP]=GetStrainPlastic(GammaDot,Stress,FCCSlips,CubicSlips)

for n1=1:18
            M1(n1).M=zeros(3,3);
            M2(n1).M=zeros(3,3);
            M3(n1).M=zeros(3,3);
            StrainP(n1).M=zeros(3,3);
            StrainP(n1).Ep=zeros(3,3);
end

%%
for n1=1:12
    for nk=1:3
        for nl=1:3
            M1(n1).M=L(:,:,nk,nl)*FCCSlips(n1).mu;
            M2(n1).M=FCCSlips(n1).ohm*Stress.';
            M3(n1).M=Stress.'*FCCSlips(n1).ohm.';
        end
    end
    
end

for n1=1:6
    for nk=1:3
        for nl=1:3
            M1(n1).M=L(:,:,nk,nl)*CubicSlips(n1).mu;
            M2(n1).M=CubicSlips(n1).ohm*Stress.';
            M3(n1).M=Stress.'*CubicSlips(n1).ohm.';
        end
    end
end

%%
for n1=1:18
    StrainP(n1).M=M1(n1).M+M2(n1).M+M3(n1).M;
    StrainP(n1).Ep=StrainP(n1).M.*GammaDot(n1);
end


end