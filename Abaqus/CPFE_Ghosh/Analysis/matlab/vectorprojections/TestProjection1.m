clear;

V1=[1;1;1];
[v1,m1]=NormalisedVector(V1);

theta=50;

x=cosd(theta);
y=sind(theta);
v=[x;y;0];
X=[1;0;0];

[p,c,s]=CosineProject(X,v);

%% 