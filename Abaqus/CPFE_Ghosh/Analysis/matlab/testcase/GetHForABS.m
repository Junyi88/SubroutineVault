function [H,tPE,tSE,tCB,Inx]=GetHForABS(TauPE,TauSE,TauCB,b,Ga111,Ga010,G,h,k1,k2)

M=b./Ga111;
tPE=TauPE.*M;
tSE=TauSE.*M;
tCB=TauCB.*M;

B=G.*b.*b./(2.*pi.*Ga111);
CH=G.*b.*b.*b./(4.*pi);

M=(((1./sqrt(3))-(Ga010./Ga111)+abs(tCB)).*b./B);
for n1=1:12
    if M(n1)<0
       M(n1)=0; 
    else
        M(n1)=sqrt(M(n1));
    end
    M(n1)=sqrt(M(n1));
end
Inx=(h+k1.*(tPE-k2.*tSE)+M);
H=CH.*Inx;
end