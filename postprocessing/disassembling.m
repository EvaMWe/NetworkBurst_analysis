% this is a function to disassemble the result Cell that contains all
% featuers to several features, to have n > p in the end

function result_fragments = disassembling(resultCell)
fragMat = [1 6;7 20;21 49; 50 91;92 120; 121 132];
nbFrag = length(fragMat);
nbTimeBins =size(resultCell,1);
Variables = resultCell{2,2};
wells = resultCell{1,2};

result_fragments = cell(nbTimeBins,nbFrag,2);
result_fragments{1,2,2} = wells;
for f = 1:nbFrag
    fr1 = fragMat(f,1);
    fr2 = fragMat(f,2);    
    for b = 1: nbTimeBins
        cellT = resultCell{b,1};
        cellT = cellT(fr1:fr2,:);
        result_fragments{b,f,1} = cellT;
    end
    result_fragments{f,1,2} = Variables(fr1:fr2);    
end
end