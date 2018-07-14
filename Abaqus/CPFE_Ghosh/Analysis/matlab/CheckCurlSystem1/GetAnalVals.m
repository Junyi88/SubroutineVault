function [U,Fp,FpT,FpN,FpTN,dFpN,dFpTN]=GetAnalVals(x,y,z,N,Co)

%%
U=zeros(3,1);
for n1=1:3
U(n1)=Co(n1,1)*x*y*z+...
    Co(n1,2)*y*z+...
    Co(n1,3)*x*z+...
    Co(n1,4)*x*y+...
    Co(n1,5)*x+...
    Co(n1,6)*y+...
    Co(n1,7)*z+...
    Co(n1,8);
end

%%
Fp=zeros(3,3);

for n1=1:3
Fp(n1,1)=Co(n1,1)*y*z+Co(n1,3)*z+Co(n1,4)*y+Co(n1,5);
Fp(n1,2)=Co(n1,1)*x*z+Co(n1,2)*z+Co(n1,4)*x+Co(n1,6);
Fp(n1,3)=Co(n1,1)*x*y+Co(n1,2)*y+Co(n1,3)*x+Co(n1,7);
end

Fp=Fp+eye(3);

%%
FpT=Fp.';

FpN=Fp*N;
FpTN=FpT*N;

%%
dFpN=zeros(3,1);
dFpN(1)=(Co(3,1)*z+Co(3,1)*x+Co(3,2)+Co(3,4))*N(3)-...
    (Co(2,1)*y+Co(2,1)*x+Co(2,2)+Co(2,3))*N(2);
dFpN(2)=(Co(1,1)*y+Co(1,1)*x+Co(1,2)+Co(1,3))*N(1)-...
    (Co(3,1)*z+Co(3,1)*y+Co(3,3)+Co(3,4))*N(3);
dFpN(3)=(Co(2,1)*z+Co(2,1)*y+Co(2,3)+Co(2,4))*N(2)-...
    (Co(1,1)*x+Co(1,1)*y+Co(1,2)+Co(1,4))*N(1);

%%
dFpTN=zeros(3,1);
A1=Co(1,1)+Co(2,1)+Co(3,1);
A2=Co(1,2)+Co(2,2)+Co(3,2);
B1=Co(1,1)+Co(2,1)+Co(3,1);
B2=Co(1,2)+Co(2,2)+Co(3,2);
dFpTN(1)=(A1*x+A2)*N(3)-...
    (B1*x+B2)*N(2);


A1=Co(1,1)+Co(2,1)+Co(3,1);
A2=Co(1,3)+Co(2,3)+Co(3,3);
B1=Co(1,1)+Co(2,1)+Co(3,1);
B2=Co(1,3)+Co(2,3)+Co(3,3);
dFpTN(2)=(A1*y+A2)*N(1)-...
    (B1*y+B2)*N(3);

A1=Co(1,1)+Co(2,1)+Co(3,1);
A2=Co(1,4)+Co(2,4)+Co(3,4);
B1=Co(1,1)+Co(2,1)+Co(3,1);
B2=Co(1,4)+Co(2,4)+Co(3,4);
dFpTN(3)=(A1*z+A2)*N(2)-...
    (B1*z+B2)*N(1);

end