clear;

L=zeros(3,3,3,3);

ToVoigtUMAT=[1 4 5;
    4 2 6;
    5 6 3];

FromVoigtUMAT=[1 1;2 2;3 3;1 2;1 3;2 3];

for ni=1:3
    for nj=1:3
        for nk=1:3
            for nl=1:3
                L(ni,nj,nk,nl)=ToVoigtUMAT(ni,nj).*10+ToVoigtUMAT(nk,nl);
            end
        end
    end
end

%%
