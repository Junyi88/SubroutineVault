function [FEVals]=GetFEVals(X,Y,Z)

FEVals.NatCoords=[-1 1 1;...
    -1 -1 1;...
    -1 1 -1;...
    -1 -1 -1;...
    1 1 1;...
    1 -1 1;...
    1 1 -1;...
    1 -1 -1];

FEVals.NatCoordsQ=FEVals.NatCoords./sqrt(3);

%%
FEVals.Coords=[X.';Y.';Z.'];

FEVals.GaussCoords=zeros(3,8);
for nG2=1:8
for nG=1:8
    for nI=1:3
        FEVals.GaussCoords(nI,nG2)=FEVals.GaussCoords(nI,nG2)+...
            FEVals.Coords(nI,nG).*...
            (1+FEVals.NatCoords(nG,1).*FEVals.NatCoordsQ(nG2,nI)).*...
            (1+FEVals.NatCoords(nG,2).*FEVals.NatCoordsQ(nG2,nI)).*...
            (1+FEVals.NatCoords(nG,3).*FEVals.NatCoordsQ(nG2,nI));
    end
end
end
FEVals.GaussCoords=FEVals.GaussCoords./8;

%%
for n1=1:8
      yp = FEVals.NatCoordsQ(n1,1);
      yq = FEVals.NatCoordsQ(n1,2);
      yr = FEVals.NatCoordsQ(n1,3);
    
   [xn,dndloc]=kshapes8(yp,yq,yr,FEVals.NatCoords);
   FEVals.Abq(n1).xn=xn; 
   FEVals.Abq(n1).dndloc=dndloc; 
end

%%
for n1=1:8
   FEVals.Abq(n1).xj=FEVals.GaussCoords*FEVals.Abq(n1).dndloc;
   FEVals.Abq(n1).xjI=inv(FEVals.Abq(n1).xj);
   
   FEVals.Abq(n1).dndx=FEVals.Abq(n1).dndloc*FEVals.Abq(n1).xjI;
end





end


%%

function [xn,dndloc]=kshapes8(yp,yq,yr,xnat)

xn=zeros(8,1);
dndloc=zeros(8,3);
gn=zeros(8,1);
gnr=zeros(8,1);
gns=zeros(8,1);
gnt=zeros(8,1);

for ni=1:8
    ri=xnat(ni,1);
    si=xnat(ni,2);
    ti=xnat(ni,3);
    
    if(ri == 1. || ri==-1.) 
        gr = 0.5*(1.+ri*yp);
        dgr = 0.5*ri;
    else
        gr = (1.-yp*yp);
        dgr = -2.0*yp;
    end 

    if(si == 1. || si==-1.)
        gs = 0.5*(1.+si*yq);
        dgs = 0.5*si;
    else
        gs = (1.-yq*yq);
        dgs = -2.0*yq;
    end
    
    if(ti == 1. || ti==-1.)
        gt = 0.5*(1.+ti*yr);
        dgt = 0.5*ti;
    else
        gt = (1.-yr*yr);
        dgt = -2.0*yr;
    end
    
    gn(ni) = gr*gs*gt;
    
    gnr(ni) = dgr*gs*gt;
    gns(ni) = gr*dgs*gt;
    gnt(ni) = gr*gs*dgt;
    
    
end
   
      xn(1) = gn(1);
      xn(2) = gn(2);
      xn(3) = gn(3);
      xn(4) = gn(4);
      xn(5) = gn(5);
      xn(6) = gn(6);
      xn(7) = gn(7);
      xn(8) = gn(8);
      

      dndloc(1,1) = gnr(1);
      dndloc(2,1) = gnr(2);
      dndloc(3,1) = gnr(3);
      dndloc(4,1) = gnr(4);
      dndloc(5,1) = gnr(5);
      dndloc(6,1) = gnr(6);
      dndloc(7,1) = gnr(7);
      dndloc(8,1) = gnr(8);
      
      dndloc(1,2) = gns(1);
      dndloc(2,2) = gns(2);
      dndloc(3,2) = gns(3);
      dndloc(4,2) = gns(4);
      dndloc(5,2) = gns(5);
      dndloc(6,2) = gns(6);
      dndloc(7,2) = gns(7);
      dndloc(8,2) = gns(8);
      
      dndloc(1,3) = gnt(1);
      dndloc(2,3) = gnt(2);
      dndloc(3,3) = gnt(3);
      dndloc(4,3) = gnt(4);
      dndloc(5,3) = gnt(5);
      dndloc(6,3) = gnt(6);
      dndloc(7,3) = gnt(7);
      dndloc(8,3) = gnt(8);

end