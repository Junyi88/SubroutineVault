clear;

C=zeros(6,6);
for ni=1:6
    for nj=1:6
        C(ni,nj)=ni.*1+nj;
    end
end

%%
L=zeros(3,3,3,3);

Map1=[1 4 5;4 2 6;5 6 3];
for ni=1:6
    for nj=1:6
        CL(ni,nj).c=[];
    end
end

CLN=[];
for ni=1:3
    for nj=1:3
       for nk=1:3
           for nl=1:3
               
               na=Map1(ni,nj);
               nb=Map1(nk,nl);
               
               P=1000.*ni+100.*nj+10.*nk+nl;
               CL(na,nb).c=[CL(na,nb).c; ni,nj,nk,nl];
               CLN=[CLN;ni,nj,nk,nl,na,nb];
           end
       end 
    end
end


%% Need functions

% ConvertVoigtToFull X
% RotateFull
% ConvertFullToVoigtList
% ExtractNumberVoigt