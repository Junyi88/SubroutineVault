clear;

load abaqusdata;


F=[2 -2 0; 1 1 0; 0 0 1];
% F=Fa;
C=(F.')*F;

[A,Cs]=eig(C);
Us=sqrt(Cs);
UsI=zeros(3,3);
for n1=1:3
UsI(n1,n1)=1./Us(n1,n1);
end

U1=(A.')*(Us*A);
U2=(A)*(Us*A.'); % Correct

UI1=(A.')*(UsI*A);
UI2=(A)*(UsI*A.'); % Correct

R=F*(UI2);
eul = rotm2eul(R)
