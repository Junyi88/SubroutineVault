function [COut]=ExtractNumberVoigt(PI,PJ,CVoigtList,PosList)

COut=[];

for n1=1:81
    if PosList(n1,5)==PI
    if PosList(n1,6)==PJ
        COut=[COut;CVoigtList(n1)];
    end
    end
    
end

end