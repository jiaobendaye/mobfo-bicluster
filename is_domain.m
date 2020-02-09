function is_ =is_domain(i, j)
    %if Sj domain Si, return 1, else return 0
    %input:
    %      i    [msrs, -gv,  -cv]
    %      j    [msrs, -gv,  -cv]
    %   
    if sum(j <= i) == size(i,2) && sum(j < i) > 0 
        is_ = 1;
    else
        is_ = 0;
    end

end