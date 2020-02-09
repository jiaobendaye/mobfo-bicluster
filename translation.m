function submatI = translation(nestI,data)
    %input:
    %   nestI        (n+m)   0/1 bits
    %   data        gene expression data n by m
    %output:
    %   submatI     (submat)      
        n = size(data,1);
        rInd = find(nestI(1:n));
        cInd = find(nestI(n+1:size(nestI,2)));
        subR = data(rInd,:);
        submatI = subR(:,cInd);
    end
    