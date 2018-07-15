function [curlFpN,curlFpTN]=GetSquaredAnal(x,y,z,N,A0,B0,C0)

    curlFpN=zeros(3,1);
    curlFpTN=zeros(3,1);

    
    %%
    curlFpN(1)=2*C0*(2*x*y*z*z+2*x*x*y*z+x*x*z*z)*N(3)-...
        2*B0*(2*x*y*y*z+2*x*x*y*z+x*x*y*y)*N(2);
    curlFpN(2)=2*A0*(2*x*y*y*z+2*x*x*y*z+x*x*y*y)*N(1)-...
        2*C0*(2*x*y*z*z+2*x*y*y*z+z*z*y*y)*N(3);
    curlFpN(3)=2*B0*(2*x*y*z*z+2*x*y*y*z+x*x*z*z)*N(2)-...
        2*A0*(2*x*y*z*z+2*x*y*y*z+z*z*x*x)*N(1);
    
    %%
    curlFpTN(1)=4*(A0+B0+C0)*()
end