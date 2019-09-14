function [vbar,magv]=NormalisedVector(v)

    magv=sqrt(dot(v,v));
    vbar=v./magv;

end