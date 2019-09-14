function [STRESS]=CalcStress(C,Strain)

    STR=[Strain(1) Strain(2) Strain(3);
        Strain(4) Strain(5) Strain(6);
        Strain(7) Strain(8) Strain(9)];
    S=zeros(3,3);

   
    for ni=1:3
    for nj=1:3
    for nk=1:3
    for nl=1:3
      S(ni,nj)=S(ni,nj)+C(ni,nj,nk,nl).*STR(nk,nl);
    end
    end
    end
    end

    STRESS=[S(1,1);S(1,2);S(1,3); ...
        S(2,1);S(2,2);S(2,3); ...
        S(3,1);S(3,2);S(3,3)];
end