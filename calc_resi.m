function [resi]=calc_resi(nest,data)
    %input:
    %   nest        nPop*(n+m)   0/1 bits        
    %
    %output:
    %   resi   nPop*1
        nPop  = size(nest,1);
        resi = 500*ones(nPop,1);    
        for i=1:nPop
            submatI = translation(nest(i,:),data);
            resi(i) = residu(submatI);
        end
        
end