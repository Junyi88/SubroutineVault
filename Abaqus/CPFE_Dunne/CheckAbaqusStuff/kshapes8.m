function [xn,dndloc]=kshapes8(yp,yq,yr,xnat8)  


gn=zeros(8,1);

% gp=zeros(8,1);
% gq=zeros(8,1);
% gr=zeros(8,1);
% 
% dgp=zeros(8,1);
% dgq=zeros(8,1);
% dgr=zeros(8,1);

gnp=zeros(8,1);
gnq=zeros(8,1);
gnr=zeros(8,1);

dndloc=zeros(8,3);

%%
for n1=1:8
    gp=0.5*(1+yp*xnat8(n1,1));
    gq=0.5*(1+yq*xnat8(n1,2));
    gr=0.5*(1+yr*xnat8(n1,3));
    
    dgp=0.5*(xnat8(n1,1));
    dgq=0.5*(xnat8(n1,2));
    dgr=0.5*(xnat8(n1,3));
    
    gn(n1)=gp*gq*gr;
    gnp(n1)=dgp*gq*gr;
    gnq(n1)=gp*dgq*gr;
    gnr(n1)=gp*gq*dgr;
    
    dndloc(n1,1)=gnp(n1);
    dndloc(n1,2)=gnq(n1);
    dndloc(n1,3)=gnr(n1);
end
xn=gn;

end