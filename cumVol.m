function vol = cumVol(nest)
    %input:
    %   nest        nPop*(data_n+data_m)   0/1 bits     
    %output:
    %          
    %   vol        nPop*3  (bic,gene and sample Volume)
        global data_n
        nPop = size(nest,1);
        vol = zeros(nPop,3);
        for i = 1:nPop
            nestI = nest(i,:);
            r = nestI(1:data_n);
            c = nestI(data_n+1:size(nestI,2));
            countR = sum(r);
            countC = sum(c);
            bicVol = countC*countR;
    
            vol(i,1) = bicVol;
            vol(i,2) = countR;
            vol(i,3) = countC;
    
        end
end
