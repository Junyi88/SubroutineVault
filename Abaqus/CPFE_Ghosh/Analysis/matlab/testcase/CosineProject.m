function [p,c,s]=CosineProject(a,b)

[abar,maga]=NormalisedVector(a);
[bbar,magb]=NormalisedVector(b);

p=dot(b,abar);
c2=dot(p,p);
c=sqrt(c2);
s=sqrt(1-c2);

end