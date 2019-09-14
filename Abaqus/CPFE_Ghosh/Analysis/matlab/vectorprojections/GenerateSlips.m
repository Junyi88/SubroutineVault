clear;

%%
Slips(1).n=[1;1;1];
Slips(2).n=[1;1;1];
Slips(3).n=[1;1;1];

Slips(4).n=[-1;1;1];
Slips(5).n=[-1;1;1];
Slips(6).n=[-1;1;1];

Slips(7).n=[1;-1;1];
Slips(8).n=[1;-1;1];
Slips(9).n=[1;-1;1];

Slips(10).n=[1;1;-1];
Slips(11).n=[1;1;-1];
Slips(12).n=[1;1;-1];

%%
Slips(1).s=[0;-1;1];
Slips(2).s=[1;0;-1];
Slips(3).s=[-1;1;0];

Slips(4).s=[1;0;1];
Slips(5).s=[1;1;0];
Slips(6).s=[0;-1;1];

Slips(7).s=[0;1;1];
Slips(8).s=[1;1;0];
Slips(9).s=[1;0;-1];

Slips(10).s=[0;1;1];
Slips(11).s=[1;0;1];
Slips(12).s=[-1;1;0];

%%

for n1=1:12
   NSlips(n1).s=NormalisedVector(Slips(n1).s); 
   NSlips(n1).n=NormalisedVector(Slips(n1).n); 
   NSlips(n1).t=cross(Slips(n1).n,Slips(n1).s); 
end

save FCCSlips;
