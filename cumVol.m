function vol = cumVol(nest,n)
    %input:
    %   nest        nPop*(n+m)   0/1 bits     
    %output:
    %          
    %   vol        nPop*3  (bic,gene and sample Volume)
        nPop = size(nest,1);
        vol = zeros(nPop,3);
        for i = 1:nPop
            nestI = nest(i,:);
            r = nestI(1:n);
            c = nestI(n+1:size(nestI,2));
            countR = sum(r);
            countC = sum(c);
            bicVol = countC*countR;
    
            vol(i,1) = bicVol;
            vol(i,2) = countR;
            vol(i,3) = countC;
    
        end
end
