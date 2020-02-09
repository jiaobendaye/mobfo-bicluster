function scores = calc_scores(nests, data)
    % input:
    %     nests       nPop*(n+m) 0/1 bits
    %     data        n*m         datasets
    %output:
    %     scores      nPop*3 [msr, -gv, -cv]
    n = size(data,1);
    % [msr,-gv,-cv]
    MSRs = calc_resi(nests, data);
    vols = cumVol(nests, n);
    scores = [MSRs,-vols(:,2:3)];
end