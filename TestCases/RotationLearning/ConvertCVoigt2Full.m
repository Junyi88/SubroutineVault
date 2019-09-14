function [CFull]=ConvertCVoigt2Full(CVoigt)
    
    CFull=zeros(3,3,3,3);
    Map=[1 4 5;
        4 2 6;
        5 6 3];
    for ni=1:3
    for nj=1:3
    for nk=1:3
    for nl=1:3
      na=Map(ni,nj);
      nb=Map(nk,nl);
      CFull(ni,nj,nk,nl)=CVoigt(na,nb);
    end
    end
    end
    end

end