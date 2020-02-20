function [vari]=calc_var(nest)
    %input:
    %   nest        nPop*(n+m)   0/1 bits        
    %
    %output:
    %   vari   nPop*1
        nPop  = size(nest,1);
        vari = zeros(nPop,1);    
        for i = 1:nPop
            bicI = translation(nest(i,:));
            v = var(bicI,0,2);
            vari(i) = mean(v);
        end       
end