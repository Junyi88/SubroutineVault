function [TotalStrain]=GetTotalStrains(dStress,EPlastic,L)

Map33to9=[1 2 3;4 5 6;7 8 9];
Map9to33=[1,1;1 2; 1 3;2 1; 2 2;2 3;3 1;3 2;3 3];

AMap=zeros(9,9);



%% Fill L

for n1=1:9
    for n2=1:9
       ni=Map9to33(n1,1);
       nj=Map9to33(n1,2);
       nk=Map9to33(n2,1);
       nl=Map9to33(n2,2);
       
       AMap(n1,n2)=AMap(n1,n2)+L(ni,nj,nk,nl);
       
       if nk==nl
         AMap(n1,n2)=AMap(n1,n2)-Stress(ni,nk);
       end
    end
end

St=zeros(9,1);
Ep=zeros(9,1);
for n1=1:3
    for n2=1:3
        Ep(Map33to9(n1,n2))=EPlastic(n1,n2);
        St(Map33to9(n1,n2))=dStress(n1,n2);
    end
end

LHS=St+Ep;
RHS=AMap\LHS;

TotalStrain=zeros(3,3);
for n1=1:9
    TotalStrain(Map9to33(n1,1),Map9to33(n1,2))=[RHS(n1)];
end

end
