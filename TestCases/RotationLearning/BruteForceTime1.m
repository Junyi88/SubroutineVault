clear;

YM=70e3;
nu=0.3;
CVoigt=zeros(6,6);
CVoigt(1,1)=(1-nu);
CVoigt(2,2)=(1-nu);
CVoigt(3,3)=(1-nu);
CVoigt(4,4)=(1-2.*nu);
CVoigt(5,5)=(1-2.*nu);
CVoigt(6,6)=(1-2.*nu);

CVoigt(1,2)=nu;
CVoigt(1,3)=nu;
CVoigt(2,3)=nu;
CVoigt(2,1)=nu;
CVoigt(3,1)=nu;
CVoigt(3,2)=nu;

CVoigt=CVoigt.*YM./((1+nu).*(1-2.*nu));
[CFull]=ConvertCVoigt2Full(CVoigt);

t01=linspace(0,1,51);
t12=linspace(0,1,51);

t=[t01,t12+t01(end)];
a=0.1;
omega=pi./2;

x01=zeros(3,length(t01));
x12=zeros(3,length(t01));

x01(1,:)=1+a.*t01;
x12(1,:)=(1+a).*cos(omega.*t12);
x12(2,:)=-(1+a).*sin(omega.*t12);

x=[x01,x12];

%%
figure(1);
clf;
hold on;
plot(t,x(1,:),'r-');
plot(t,x(2,:),'b-');
xlabel('Time');
ylabel('x');

%% 

F01=zeros(9,length(t01));
F12=zeros(9,length(t01));

F01(5,:)=1;
F01(9,:)=1;
F01(1,:)=(1+a.*t01);

F12(1,:)=(1+a).*cos(omega.*t12);
F12(2,:)=-sin(omega.*t12);
F12(4,:)=(1+a).*sin(omega.*t12);
F12(5,:)=cos(omega.*t12);
F12(9,:)=1;

F=[F01,F12];
%%
figure(2);
clf;
hold on;
plot(t,F(1,:),'rx-');
plot(t,F(2,:),'go-');
plot(t,F(4,:),'b-');
plot(t,F(5,:),'k-');

xlabel('Time');
ylabel('F');
%%
STRAIN=zeros(9,length(t));

for n1=1:length(t)
   f=F(:,n1); 
   [STR]=CalcStrain(f);
   STRAIN(:,n1)=STR;
end

figure(3);
clf;
hold on;
plot(t,STRAIN(1,:),'r>-');
plot(t,STRAIN(2,:),'gs-');
plot(t,STRAIN(3,:),'b>-');

plot(t,STRAIN(4,:),'rx-');
plot(t,STRAIN(5,:),'gx-');
plot(t,STRAIN(6,:),'bx-');

plot(t,STRAIN(7,:),'ro-');
plot(t,STRAIN(8,:),'go-');
plot(t,STRAIN(9,:),'bo-');

xlabel('Time');
ylabel('Strain');

%%
STRESS=zeros(9,length(t));

for n1=1:length(t)
   [STRESS(:,n1)]=CalcStress(CFull,STRAIN(:,n1));
end

figure(4);
clf;
hold on;
plot(t,STRESS(1,:),'r-');
plot(t,STRESS(2,:),'gx-');
plot(t,STRESS(3,:),'b-');
plot(t,STRESS(5,:),'k-');
plot(t,STRESS(6,:),'c-');
plot(t,STRESS(9,:),'mp-');

xlabel('Time');
ylabel('Stress');

%%

Ro01=zeros(9,length(t01));
Ro12=zeros(9,length(t12));

Ro01(1,:)=1;
Ro01(5,:)=1;
Ro01(9,:)=1;


Ro12(1,:)=cos(omega.*t12);
Ro12(2,:)=-sin(omega.*t12);
Ro12(4,:)=sin(omega.*t12);
Ro12(5,:)=cos(omega.*t12);

Ro12(9,:)=1;

Ro=[Ro01,Ro12];

STRESSR=zeros(9,length(t));

for n1=1:length(t)
   [STRESSR(:,n1)]=RotateStress(Ro(:,n1),STRESS(:,n1));
end
figure(5);
clf;
hold on;
plot(t,STRESSR(1,:),'r-');
plot(t,STRESSR(2,:),'gx-');
plot(t,STRESSR(3,:),'b-');
plot(t,STRESSR(5,:),'k-');
plot(t,STRESSR(6,:),'c-');
plot(t,STRESSR(9,:),'mp-');

xlabel('Time');
ylabel('Stress');