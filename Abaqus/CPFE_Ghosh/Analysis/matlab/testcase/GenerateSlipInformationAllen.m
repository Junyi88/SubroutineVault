clear;

%%
FCCSlips(1).n=[1;1;1];
FCCSlips(2).n=[1;1;1];
FCCSlips(3).n=[1;1;1];

FCCSlips(6).n=[1;-1;-1];
FCCSlips(4).n=[1;-1;-1];
FCCSlips(5).n=[1;-1;-1];

FCCSlips(7).n=-[1;-1;1];
FCCSlips(9).n=-[1;-1;1];
FCCSlips(8).n=-[1;-1;1];

FCCSlips(10).n=-[1;1;-1];
FCCSlips(11).n=-[1;1;-1];
FCCSlips(12).n=-[1;1;-1];

%%
FCCSlips(1).s=-[0;-1;1];
FCCSlips(2).s=-[1;0;-1];
FCCSlips(3).s=-[-1;1;0];

FCCSlips(6).s=[0;-1;1];
FCCSlips(4).s=[-1;0;-1];
FCCSlips(5).s=[1;1;0];

FCCSlips(7).s=[0;1;1];
FCCSlips(9).s=[1;0;-1];
FCCSlips(8).s=-[1;1;0];

FCCSlips(10).s=-[0;1;1];
FCCSlips(11).s=[1;0;1];
FCCSlips(12).s=[-1;1;0];

%%
for n1=1:12
FCCSlips(n1).npe=FCCSlips(n1).n;
end

FCCSlips(1).spe=[-2;1;1];
FCCSlips(2).spe=[1;-2;1];
FCCSlips(3).spe=[1;1;-2];

FCCSlips(6).spe=[-2;-1;-1];
FCCSlips(4).spe=[1;2;-1];
FCCSlips(5).spe=[1;-1;2];

FCCSlips(7).spe=[2;1;-1];
FCCSlips(9).spe=[-1;-2;-1];
FCCSlips(8).spe=[-1;1;2];

FCCSlips(10).spe=[2;-1;1];
FCCSlips(11).spe=[-1;2;1];
FCCSlips(12).spe=[-1;-1;-2];

%%
FCCSlips(1).nse=[1;-1;-1];
FCCSlips(2).nse=[-1;1;-1];
FCCSlips(3).nse=[-1;-1;1];

FCCSlips(6).nse=[1;1;1];
FCCSlips(4).nse=[-1;-1;1];
FCCSlips(5).nse=[-1;1;-1];

FCCSlips(7).nse=[-1;-1;1];
FCCSlips(9).nse=[1;1;1];
FCCSlips(8).nse=[1;-1;-1];

FCCSlips(10).nse=[-1;1;-1];
FCCSlips(11).nse=[1;-1;-1];
FCCSlips(12).nse=[1;1;1];

%
FCCSlips(1).sse=[-2;-1;-1];
FCCSlips(2).sse=[-1;-2;-1];
FCCSlips(3).sse=[-1;-1;-2];

FCCSlips(6).sse=[-2;1;1];
FCCSlips(4).sse=[-1;2;1];
FCCSlips(5).sse=[-1;1;2];

FCCSlips(7).sse=[2;-1;1];
FCCSlips(9).sse=[1;-2;1];
FCCSlips(8).sse=[1;-1;2];

FCCSlips(10).sse=[2;1;-1];
FCCSlips(11).sse=[1;2;-1];
FCCSlips(12).sse=[1;1;-2];

%%
FCCSlips(1).ncb=[1;0;0];
FCCSlips(2).ncb=[0;1;0];
FCCSlips(3).ncb=[0;0;1];

FCCSlips(6).ncb=[1;0;0];
FCCSlips(4).ncb=[0;-1;0];
FCCSlips(5).ncb=[0;0;-1];

FCCSlips(7).ncb=[-1;0;0];
FCCSlips(8).ncb=[0;1;0];
FCCSlips(9).ncb=[0;0;-1];

FCCSlips(10).ncb=[-1;0;0];
FCCSlips(11).ncb=[0;-1;0];
FCCSlips(12).ncb=[0;0;1];

%
FCCSlips(1).scb=[0;1;-1];
FCCSlips(2).scb=[-1;0;1];
FCCSlips(3).scb=[1;-1;0];

FCCSlips(6).scb=[0;-1;1];
FCCSlips(4).scb=[-1;0;-1];
FCCSlips(5).scb=[1;1;0];

FCCSlips(7).scb=[0;1;1];
FCCSlips(9).scb=[1;0;-1];
FCCSlips(8).scb=[-1;-1;0];

FCCSlips(10).scb=[0;-1;-1];
FCCSlips(11).scb=[1;0;1];
FCCSlips(12).scb=[-1;1;0];

%% Normallise 
for n1=1:12
[FCCSlips(n1).n,~]=NormalisedVector(FCCSlips(n1).n);
[FCCSlips(n1).s,~]=NormalisedVector(FCCSlips(n1).s);

[FCCSlips(n1).npe,~]=NormalisedVector(FCCSlips(n1).npe);
[FCCSlips(n1).spe,~]=NormalisedVector(FCCSlips(n1).spe);

[FCCSlips(n1).nse,~]=NormalisedVector(FCCSlips(n1).nse);
[FCCSlips(n1).sse,~]=NormalisedVector(FCCSlips(n1).sse);

[FCCSlips(n1).ncb,~]=NormalisedVector(FCCSlips(n1).ncb);
[FCCSlips(n1).scb,~]=NormalisedVector(FCCSlips(n1).scb);

FCCSlips(n1).t=cross(FCCSlips(n1).s,FCCSlips(n1).n);
FCCSlips(n1).mu=0.5*(FCCSlips(n1).s*FCCSlips(n1).n.'+FCCSlips(n1).n*FCCSlips(n1).s.');
FCCSlips(n1).ohm=0.5*(FCCSlips(n1).s*FCCSlips(n1).n.'-FCCSlips(n1).n*FCCSlips(n1).s.');

end

%% ----------------------------
%%
CubicSlips(1).n=[1;0;0];
CubicSlips(2).n=[1;0;0];

CubicSlips(3).n=[0;1;0];
CubicSlips(4).n=[0;1;0];

CubicSlips(5).n=[0;0;1];
CubicSlips(6).n=[0;0;1];


%%
CubicSlips(1).s=[0;1;1];
CubicSlips(2).s=[0;1;-1];

CubicSlips(3).s=[1;0;1];
CubicSlips(4).s=[1;0;-1];

CubicSlips(5).s=[1;1;0];
CubicSlips(6).s=[1;-1;0];

%% Normallise 
for n1=1:6
[CubicSlips(n1).n,~]=NormalisedVector(CubicSlips(n1).n);
[CubicSlips(n1).s,~]=NormalisedVector(CubicSlips(n1).s);

CubicSlips(n1).t=cross(CubicSlips(n1).s,CubicSlips(n1).n);
CubicSlips(n1).mu=0.5*(CubicSlips(n1).s*CubicSlips(n1).n.'+CubicSlips(n1).n*CubicSlips(n1).s.');
CubicSlips(n1).ohm=0.5*(CubicSlips(n1).s*CubicSlips(n1).n.'-CubicSlips(n1).n*CubicSlips(n1).s.');
end


%% WriteOutput
FileName='FCC_S.csv';
Out1=zeros(12,3);
for n1=1:12
   Out1(n1,:)=FCCSlips(n1).s.'; 
end
dlmwrite(FileName,Out1,'precision','%12.12e');

save SlipSystemsAllen;