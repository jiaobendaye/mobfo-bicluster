function result = calc_bench(pop,data)
%pop    nPop*(n+m)   0/1 
%data   n*m
%result nPop*5        residue,bicV,geneV,condV,var
    n = size(data,1);

    resi = calc_resi(pop,data);
    Vol = cumVol(pop,n);
    vari = calc_var(pop,data);

    result = [resi Vol vari];
end