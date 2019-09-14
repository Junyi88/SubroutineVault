function [CurlOut,dfdx]=GetGaussianCurl(fc,xnat8,gauss,gausscoords)

xnmat=zeros(8,8);
fnode=zeros(3,8);
dfdx=zeros(3,8);
CurlOut=zeros(8,3);


for n1=1:8
yp=gauss(n1,1);
yq=gauss(n1,2);
yr=gauss(n1,3);
[xn,dndloc]=kshapes8(yp,yq,yr,xnat8);
xnmat(n1,:)=xn;
xj=gausscoords*dndloc;
xjinv=inv(xj);
dndx=dndloc*xjinv;

for n2=1:8
    fnode(1,n2)=fc(1,n2);
    fnode(2,n2)=fc(2,n2);
    fnode(3,n2)=fc(3,n2);
end

fmat=fnode*dndx;

CurlOut(n1,1)=fmat(3,2)-fmat(2,3);
CurlOut(n1,2)=fmat(1,3)-fmat(3,1);
CurlOut(n1,3)=fmat(2,1)-fmat(1,2);
% dfdx(:,1)=fmat;
end





end

function [xn,dndloc]=kshapes8(yp,yq,yr,xnat8)

xn=zeros(8,1);
dndloc=zeros(8,3);

gn=zeros(8,1);
gnr=zeros(8,1);
gns=zeros(8,1);
gnt=zeros(8,1);

for n1=1:8
    gr=0.5*(1+yp*xnat8(n1,1));
    gs=0.5*(1+yq*xnat8(n1,2));
    gt=0.5*(1+yr*xnat8(n1,3));
    
    dgr=0.5*(xnat8(n1,1));
    dgs=0.5*(xnat8(n1,2));
    dgt=0.5*(xnat8(n1,3));
    
    gn(n1)=gr*gs*gt;
    xn(n1)=gn(n1);
    
    gnr(n1)=dgr*gs*gt;
    gns(n1)=gr*dgs*gt;
    gnt(n1)=gr*gs*dgt;
end

dndloc(:,1)=gnr;
dndloc(:,2)=gns;
dndloc(:,3)=gnt;

end