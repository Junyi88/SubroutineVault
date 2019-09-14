clear;

%% Properties
Prop_C=[296378,198565,127005];
C66=zeros(6,6);

for n1=1:3
   C66(n1,n1)=Prop_C(1); 
end

C66(1,2)=Prop_C(2);
C66(2,1)=Prop_C(2);
C66(1,3)=Prop_C(2);
C66(3,1)=Prop_C(2);
C66(2,3)=Prop_C(2);
C66(3,2)=Prop_C(2);

for n1=4:6
   C66(n1,n1)=Prop_C(3); 
end

% --
h0=200;
Taus=67.1;
Tau0=47.1;
g0=47.1;
Nex=10;
adot=0.01;

%%
S_FCC0=zeros(3,12);

N_FCC0=[...
    1,1,1;...
    1,1,1;...
    1,1,1;...
    1,-1,-1;...
    1,-1,-1;...
    1,-1,-1;...    
    -1,1,-1;...
    -1,1,-1;...
    -1,1,-1;...
    -1,-1,1;...
    -1,-1,1;...
    -1,-1,1];

S_FCC0=[...
    0,1,-1;...
    -1,0,1;...
    1,-1,0;...
    0,-1,1;...
    -1,0,-1;...
    1,1,0;...
    0,1,1;...
    1,0,-1;...
    -1,-1,0;...
    0,-1,-1;...
    1,0,1;...
    -1,1,0];

for n1=1:12
   N_FCC0(n1,:)=N_FCC0(n1,:)./sqrt(dot(N_FCC0(n1,:),N_FCC0(n1,:)));
   S_FCC0(n1,:)=S_FCC0(n1,:)./sqrt(dot(S_FCC0(n1,:),S_FCC0(n1,:)));
end

%%
dt=0.1;
dEp=[0.001;0;0;0;0;0];
Stress0=zeros(6,1);

StressTr=C66*dEp;

%%
S_FCC=S_FCC0;

%% Calculate CSNNS

CSNNS=zeros(6,12);

for n1=1:12
    for n2=1:6
      for n3=1:6
        
          if (n3<=3)
              CSNNS(n2,n1)=CSNNS(n2,n1)+C66(n2,n3).*...
                  (S_FCC(n1,n3).*N_FCC0(n1,n3)).*2;
          elseif n3==4
              CSNNS(n2,n1)=CSNNS(n2,n1)+2.*C66(n2,n3).*...
                  (S_FCC(n1,1).*N_FCC0(n1,2)+S_FCC(n1,2).*N_FCC0(n1,1));
          elseif n3==5
              CSNNS(n2,n1)=CSNNS(n2,n1)+2.*C66(n2,n3).*...
                  (S_FCC(n1,1).*N_FCC0(n1,3)+S_FCC(n1,3).*N_FCC0(n1,1));
          elseif n3==6
              CSNNS(n2,n1)=CSNNS(n2,n1)+2.*C66(n2,n3).*...
                  (S_FCC(n1,2).*N_FCC0(n1,3)+S_FCC(n1,3).*N_FCC0(n1,2));
              
          end
          
          
      end
    end
end

CSNNS=CSNNS./2;
    
%% Calculate TAu
Tau=zeros(12,1);
dGamma=zeros(12,1);
dG=zeros(12,1);
Res=zeros(30,1);
dRes=zeros(30,30);
g=ones(12,1).*g0;
Gamma=zeros(12,1);
%%
for n1=1:12
    for n2=1:6
        if (n2<=3)
            Tau(n1)=Tau(n1)+StressTr(n2).*(S_FCC(n1,n2).*N_FCC0(n1,n2));
        elseif n2==4
            Tau(n1)=Tau(n1)+StressTr(n2).*(...
                S_FCC(n1,1).*N_FCC0(n1,2)+S_FCC(n1,2).*N_FCC0(n1,1));
        elseif n2==5
            Tau(n1)=Tau(n1)+StressTr(n2).*(...
                S_FCC(n1,1).*N_FCC0(n1,3)+S_FCC(n1,3).*N_FCC0(n1,1));
        elseif n2==6
            Tau(n1)=Tau(n1)+StressTr(n2).*(...
                S_FCC(n1,2).*N_FCC0(n1,3)+S_FCC(n1,3).*N_FCC0(n1,2));
            
        end
   
    end
end

%%
Res(1:6)=StressTr-Stress0-C66*dEp+C66*(CSNNS*dGamma);

for n1=1:12
Res(6+n1)=dGamma(n1)-dt.*adot.*((abs(Tau(n1)./(g(n1)+dG(n1)))).^(Nex)).*...
    sign(Tau(n1));

Res(18+n1)=dG(n1)-h0.*(...
    sech(abs(h0.*(Gamma(n1)+dGamma(n1))./(Taus-Tau0))).^2).*abs(dGamma(n1));
end    

%% 
dRdx=eye(30,30);

dRdx(1:6,7:18)=CSNNS;

for n1=1:12
    for n2=1:6
        if (n2<=3)
             dRdx(6+n1,n2)=-dt.*adot.*Nex.*...
                 (abs(Tau(n1)./...
                 (g(n1)+dG(n1))).^(Nex-1)).*...
                 S_FCC(n1,n2).*N_FCC0(n1,n2);
        elseif n2==4
            dRdx(6+n1,n2)=-dt.*adot.*Nex.*...
                 (abs(Tau(n1)./...
                 (g(n1)+dG(n1))).^(Nex-1)).*...
                 S_FCC(n1,1).*N_FCC0(n1,2);
             
        elseif n2==5
            dRdx(6+n1,n2)=-dt.*adot.*Nex.*...
                 (abs(Tau(n1)./...
                 (g(n1)+dG(n1))).^(Nex-1)).*...
                 S_FCC(n1,1).*N_FCC0(n1,3);
             
        elseif n2==6
            dRdx(6+n1,n2)=-dt.*adot.*Nex.*...
                 (abs(Tau(n1)./...
                 (g(n1)+dG(n1))).^(Nex-1)).*...
                 S_FCC(n1,2).*N_FCC0(n1,3);
        end
        
        %----------------

        dRdx(6+n1,n1+18)=-dt.*adot.*Nex.*...
            (abs(Tau(n1)./...
            (g(n1)+dG(n1))).^(Nex)).*...
            sign(Tau(n1))./(g(n1)+dG(n1));

    end
end



