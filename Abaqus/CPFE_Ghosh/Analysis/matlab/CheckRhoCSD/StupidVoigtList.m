clear;

%%
Map_V2L=zeros(6,6);

Map_V2L(1,1)=1;
Map_V2L(2,2)=2;
Map_V2L(3,3)=3;
Map_V2L(4,4)=4;
Map_V2L(5,5)=5;
Map_V2L(6,6)=6;

Map_V2L(1,2)=7;
Map_V2L(2,1)=7;

Map_V2L(1,3)=8;
Map_V2L(3,1)=8;

Map_V2L(1,4)=9;
Map_V2L(4,1)=9;

Map_V2L(1,5)=10;
Map_V2L(5,1)=10;

Map_V2L(1,6)=11;
Map_V2L(6,1)=11;

%
Map_V2L(2,3)=12;
Map_V2L(3,2)=12;

Map_V2L(2,4)=13;
Map_V2L(4,2)=13;

Map_V2L(2,5)=14;
Map_V2L(5,2)=14;

Map_V2L(2,6)=15;
Map_V2L(6,2)=15;

%
Map_V2L(3,4)=16;
Map_V2L(4,3)=16;

Map_V2L(3,5)=17;
Map_V2L(5,3)=17;

Map_V2L(3,6)=18;
Map_V2L(6,3)=18;

%
Map_V2L(4,5)=19;
Map_V2L(5,4)=19;

Map_V2L(4,6)=20;
Map_V2L(6,4)=20;

Map_V2L(5,6)=21;
Map_V2L(6,5)=21;

%% 
Map_L2V=[...
    1 1; ...
    2 2; ...
    3 3; ...
    4 4; ...
    5 5; ...
    6 6; ...
    1 2; ...
    1 3; ...
    1 4; ...
    1 5; ...
    1 6; ...
    2 3; ...
    2 4; ...
    2 5; ...
    2 6; ...
    3 4; ...
    3 5; ...
    3 6; ...
    4 5; ...
    4 6; ...
    5 6 ...
    ];


%%
Map_V2F=[...
    1 1;...
    2 2;...
    3 3;...
    1 2;...
    1 3;...
    2 3];

%%

Map_L2F=zeros(21,4);

for n1=1:21
   ni=Map_L2V(n1,1);
   nj=Map_L2V(n1,2);
   
   Map_L2F(n1,1:2)=Map_V2F(ni,:);
   Map_L2F(n1,3:4)=Map_V2F(nj,:);
end

%%
Map_F2V=[...
    1 4 5;...
    4 2 6;...
    5 6 3];

Map_F2L=zeros(3,3,3,3);

for ni=1:3
for nj=1:3
for nk=1:3    
for nl=1:3
    
    NI=Map_F2V(ni,nj);
    NJ=Map_F2V(nk,nl);
    Map_F2L(ni,nj,nk,nl)=Map_V2L(NI,NJ);
    
end    
end    
end    
end

Map_F2Lsim=zeros(9,9);


for ni=1:3
for nj=1:3
for nk=1:3    
for nl=1:3
    
    NI=(ni-1)*3+nj;
    NJ=(nk-1)*3+nl;
    
    Map_F2Lsim(NI,NJ)=Map_F2L(ni,nj,nk,nl);
end    
end    
end    
end