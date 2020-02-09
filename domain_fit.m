function [fitness]=domain_fit(nests, scores)
    % input:
    %     nests       nPop*(n+m) 0/1bits
    %     scores      nPop*3  [msr, -gv, -cv]     
    % return:
    %     fitness       nPop*1   domained_cnt :  smaller is better
    nPop = size(nests,1);
    D = ones(nPop, nPop);

    for row=1:nPop
        for col=1:nPop
            if row == col 
                D(row, col) = 0;
                continue
            end
            D(row,col) = is_domain(scores(row,:), scores(col,:));
        end
    end
    fitness = sum(D,2);
end
    
