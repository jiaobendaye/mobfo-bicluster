function result = calc_bench(pop)
%pop    nPop*(n+m)   0/1 
%data   n*m
%result nPop*5        residue,bicV,geneV,condV,var

    resi = calc_resi(pop);
    Vol = cumVol(pop);
    vari = calc_var(pop);

    result = [resi Vol vari];
end