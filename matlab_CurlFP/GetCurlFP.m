function [FEIn,curlFP]=GetCurlFP(kFP,FEVals)

%%

for n1=1:8
    FEIn(n1).fnode1=zeros(3,8);
    FEIn(n1).fnode2=zeros(3,8);
    FEIn(n1).fnode3=zeros(3,8);
    
    %%
    for n2=1:8
        FEIn(n1).fnode1(1,n2)=kFP(n2).FP(1,1);
        FEIn(n1).fnode1(2,n2)=kFP(n2).FP(2,1);
        FEIn(n1).fnode1(3,n2)=kFP(n2).FP(3,1);
        
        FEIn(n1).fnode2(1,n2)=kFP(n2).FP(1,2);
        FEIn(n1).fnode2(2,n2)=kFP(n2).FP(2,2);
        FEIn(n1).fnode2(3,n2)=kFP(n2).FP(3,2);
        
        FEIn(n1).fnode3(1,n2)=kFP(n2).FP(1,3);
        FEIn(n1).fnode3(2,n2)=kFP(n2).FP(2,3);
        FEIn(n1).fnode3(3,n2)=kFP(n2).FP(3,3);
    end
    
    FEIn(n1).fmat1=FEIn(n1).fnode1*FEVals.Abq(n1).dndx;
    FEIn(n1).fmat2=FEIn(n1).fnode2*FEVals.Abq(n1).dndx;
    FEIn(n1).fmat3=FEIn(n1).fnode3*FEVals.Abq(n1).dndx;
    
    curlFP(n1).cFP=zeros(3,3);
    
    curlFP(n1).cFP(1,1)=FEIn(n1).fmat1(3,2)-FEIn(n1).fmat1(2,3);
    curlFP(n1).cFP(2,1)=FEIn(n1).fmat1(1,3)-FEIn(n1).fmat1(3,1);
    curlFP(n1).cFP(3,1)=FEIn(n1).fmat1(2,1)-FEIn(n1).fmat1(1,2);
    
    curlFP(n1).cFP(1,2)=FEIn(n1).fmat2(3,2)-FEIn(n1).fmat2(2,3);
    curlFP(n1).cFP(2,2)=FEIn(n1).fmat2(1,3)-FEIn(n1).fmat2(3,1);
    curlFP(n1).cFP(3,2)=FEIn(n1).fmat2(2,1)-FEIn(n1).fmat2(1,2);
    
    curlFP(n1).cFP(1,3)=FEIn(n1).fmat3(3,2)-FEIn(n1).fmat3(2,3);
    curlFP(n1).cFP(2,3)=FEIn(n1).fmat3(1,3)-FEIn(n1).fmat3(3,1);
    curlFP(n1).cFP(3,3)=FEIn(n1).fmat3(2,1)-FEIn(n1).fmat3(1,2);
end

end