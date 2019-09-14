clear;


L=zeros(3,3,3,3);

ToVoigtUMAT=[1 4 5;
    4 2 6;
    5 6 3];

for ni=1:3
    for nj=1:3
        for nk=1:3
            for nl=1:3
                L(ni,nj,nk,nl)=10.*ToVoigtUMAT(ni,nj)+ToVoigtUMAT(nk,nl);
            end
        end
    end
end

L11(:,:)=L(1,1,:,:);
L12(:,:)=L(1,2,:,:);