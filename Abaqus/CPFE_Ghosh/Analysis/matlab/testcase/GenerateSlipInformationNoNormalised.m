clear;

%%
FCCSlips(1).np=[1;1;1];
FCCSlips(2).np=[1;1;1];
FCCSlips(3).np=[1;1;1];

FCCSlips(4).np=[-1;1;1];
FCCSlips(5).np=[-1;1;1];
FCCSlips(6).np=[-1;1;1];

FCCSlips(7).np=[1;-1;1];
FCCSlips(8).np=[1;-1;1];
FCCSlips(9).np=[1;-1;1];

FCCSlips(10).np=[1;1;-1];
FCCSlips(11).np=[1;1;-1];
FCCSlips(12).np=[1;1;-1];

%%
FCCSlips(1).sb=[0;-1;1];
FCCSlips(2).sb=[1;0;-1];
FCCSlips(3).sb=[-1;1;0];

FCCSlips(4).sb=[1;0;1];
FCCSlips(5).sb=[1;1;0];
FCCSlips(6).sb=[0;-1;1];

FCCSlips(7).sb=[0;1;1];
FCCSlips(8).sb=[1;1;0];
FCCSlips(9).sb=[1;0;-1];

FCCSlips(10).sb=[0;1;1];
FCCSlips(11).sb=[1;0;1];
FCCSlips(12).sb=[-1;1;0];

%%
FCCSlips(1).ns=[-1;1;1];
FCCSlips(2).ns=[1;-1;1];
FCCSlips(3).ns=[1;1;-1];

FCCSlips(4).ns=[-1;-1;1];
FCCSlips(5).ns=[-1;1;-1];
FCCSlips(6).ns=[1;1;1];

FCCSlips(7).ns=[-1;-1;1];
FCCSlips(8).ns=[1;-1;-1];
FCCSlips(9).ns=[1;1;1];

FCCSlips(10).ns=[-1;1;-1];
FCCSlips(11).ns=[1;-1;-1];
FCCSlips(12).ns=[1;1;1];

%%
FCCSlips(1).nc=[1;0;0];
FCCSlips(2).nc=[0;1;0];
FCCSlips(3).nc=[0;0;1];

FCCSlips(4).nc=[0;1;0];
FCCSlips(5).nc=[0;0;-1];
FCCSlips(6).nc=[1;0;0];

FCCSlips(7).nc=[-1;0;0];
FCCSlips(8).nc=[0;0;1];
FCCSlips(9).nc=[0;1;0];

FCCSlips(10).nc=[1;0;0];
FCCSlips(11).nc=[0;-1;0];
FCCSlips(12).nc=[0;0;1];

%%
FCCSlips(1).se=[2;-1;-1];
FCCSlips(2).se=[-1;2;-1];
FCCSlips(3).se=[-1;1;2];

FCCSlips(4).se=[1;2;-1];
FCCSlips(5).se=[-1;1;-2];
FCCSlips(6).se=[2;1;1];

FCCSlips(7).se=[-2;-1;1];
FCCSlips(8).se=[-1;1;2];
FCCSlips(9).se=[1;2;1];

FCCSlips(10).se=[2;-1;1];
FCCSlips(11).se=[1;2;-1];
FCCSlips(12).se=[1;1;2];

%% Normallise 
for n1=1:12
% [FCCSlips(n1).np,~]=NormalisedVector(FCCSlips(n1).np);
% [FCCSlips(n1).ns,~]=NormalisedVector(FCCSlips(n1).ns);
% [FCCSlips(n1).sb,~]=NormalisedVector(FCCSlips(n1).sb);
% [FCCSlips(n1).se,~]=NormalisedVector(FCCSlips(n1).se);
% [FCCSlips(n1).nc,~]=NormalisedVector(FCCSlips(n1).nc);

FCCSlips(n1).tpb=cross(FCCSlips(n1).sb,FCCSlips(n1).np);
FCCSlips(n1).mu=0.5*(FCCSlips(n1).sb*FCCSlips(n1).np.'+FCCSlips(n1).np*FCCSlips(n1).sb.');
FCCSlips(n1).ohm=0.5*(FCCSlips(n1).sb*FCCSlips(n1).np.'-FCCSlips(n1).np*FCCSlips(n1).sb.');

FCCSlips(n1).muPE=0.5*(FCCSlips(n1).se*FCCSlips(n1).np.'+FCCSlips(n1).np*FCCSlips(n1).se.');
FCCSlips(n1).muSE=0.5*(FCCSlips(n1).se*FCCSlips(n1).ns.'+FCCSlips(n1).ns*FCCSlips(n1).se.');
FCCSlips(n1).muCB=0.5*(FCCSlips(n1).sb*FCCSlips(n1).nc.'+FCCSlips(n1).nc*FCCSlips(n1).sb.');
end

%% ----------------------------
%%
CubicSlips(1).np=[1;0;0];
CubicSlips(2).np=[1;0;0];

CubicSlips(3).np=[0;1;0];
CubicSlips(4).np=[0;1;0];

CubicSlips(5).np=[0;0;1];
CubicSlips(6).np=[0;0;1];


%%
CubicSlips(1).sb=[0;1;1];
CubicSlips(2).sb=[0;1;-1];

CubicSlips(3).sb=[1;0;1];
CubicSlips(4).sb=[1;0;-1];

CubicSlips(5).sb=[1;1;0];
CubicSlips(6).sb=[1;-1;0];

%% Normallise 
for n1=1:6
% [CubicSlips(n1).np,~]=NormalisedVector(CubicSlips(n1).np);
% [CubicSlips(n1).sb,~]=NormalisedVector(CubicSlips(n1).sb);

CubicSlips(n1).tpb=cross(CubicSlips(n1).sb,CubicSlips(n1).np);
CubicSlips(n1).mu=0.5*(CubicSlips(n1).sb*CubicSlips(n1).np.'+CubicSlips(n1).np*CubicSlips(n1).sb.');
CubicSlips(n1).ohm=0.5*(CubicSlips(n1).sb*CubicSlips(n1).np.'-CubicSlips(n1).np*CubicSlips(n1).sb.');
end

save SlipSystemsNoNormalised;