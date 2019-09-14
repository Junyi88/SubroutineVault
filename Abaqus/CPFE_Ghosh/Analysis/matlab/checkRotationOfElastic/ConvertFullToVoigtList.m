function [CVoigtList,PosList]=ConvertFullToVoigtList(CFull)

Map1=[1 4 5;4 2 6;5 6 3];
CVoigtList=zeros(81,1);
PosList=zeros(81,6);
count=1;
    for ni=1:3
        for nj=1:3
           for nk=1:3
               for nl=1:3

                   na=Map1(ni,nj);
                   nb=Map1(nk,nl);

                   PosList(count,:)=[ni,nj,nk,nl,na,nb];
                   CVoigtList(count,1)=CFull(ni,nj,nk,nl);
                   count=count+1;
               end
           end 
        end
    end
end