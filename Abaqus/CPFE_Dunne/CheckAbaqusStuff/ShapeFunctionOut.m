function [N,Np,Nq,Nr,...
    dNdp,dNdq,dNdr,...
    dNpdp,dNqdq,dNrdr...
    ]=ShapeFunctionOut(p,q,r,NatNode)

 L=size(NatNode,1);
 
 Np=zeros(L,1);
 Nq=zeros(L,1);
 Nr=zeros(L,1);
 
 dNpdp=zeros(L,1);
 dNqdq=zeros(L,1);
 dNrdr=zeros(L,1);


 for n=1:L
    Np(n)=0.5*(1+NatNode(n,1)*p);
    Nq(n)=0.5*(1+NatNode(n,2)*q);
    Nr(n)=0.5*(1+NatNode(n,3)*r);
    
    dNpdp(n)=0.5*(NatNode(n,1));
    dNqdq(n)=0.5*(NatNode(n,2));
    dNrdr(n)=0.5*(NatNode(n,3));
 end
 
 N=Np.*Nq.*Nr;
 dNdp=dNpdp.*Nq.*Nr;
 dNdq=Np.*dNqdq.*Nr;
 dNdr=Np.*dNrdr.*Nr;
end