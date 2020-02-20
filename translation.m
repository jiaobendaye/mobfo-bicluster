function submatI = translation(nestI)
    %input:
    %   nestI        (data_n+data_m)   0/1 bits
    %   data        gene expression data n by m
    %output:
    %   submatI     (submat)      
        global data data_n
        rInd = find(nestI(1:data_n));
        cInd = find(nestI(data_n+1:size(nestI,2)));
        submatI = data(rInd,cInd);
    end
    