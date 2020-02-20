function scores = calc_scores(nests)
    % input:
    %     nests       nPop*(n+m) 0/1 bits
    %     data        n*m         datasets
    %output:
    %     scores      nPop*3 [msr, -gv, -cv]

    % [msr,-gv,-cv]
    MSRs = calc_resi(nests);
    vols = cumVol(nests);
    scores = [MSRs,-vols(:,2:3)];
end