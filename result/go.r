library('clusterProfiler')
library('topGO')
library('Rgraphviz')
library('pathview')
library('org.Hs.eg.db')

data.names = c('BCLL', 'RAT', 'YC', 'PBC')
algs = c('CS', 'FA', 'CSFA', 'PSO', 'QPSO')
root = './'
p_th = 0.01
q_th = 0.2

go_csv <- function(csv.name, out.dir) {
    df <- read.csv(csv.name, stringsAsFactors = FALSE)
    bic.cnt = dim(df)[1]
    for (i in 1:bic.cnt) {
        i_format = sprintf("%04d", i)
        out.csv = paste(out.dir, i_format,'.csv', sep='')
        print(out.csv)
        if (file.exists(out.csv)) {next}

        genes = unlist(strsplit(df[i,'Genes'], split='+', fixed=T))
        genes_entrez_id = mapIds(x=org.Hs.eg.db, column = 'ENTREZID',
                                    keys=as.character(genes), keytype='SYMBOL')

        enrich_go <- enrichGO(gene = genes_entrez_id,
                            OrgDb = 'org.Hs.eg.db',
                            ont = 'ALL', pvalueCutoff = p_th, 
                            qvalueCutoff = q_th) 
        res.df <- enrich_go@result
        selected <- res.df[res.df$p.adjust<p_th,]
        write.csv(selected, out.csv)
    }
}
for (dataN in data.names) {
    for (alg in algs) {
        csv.name = paste(root, dataN, '/', alg, '.csv', sep='')
        out.dir = paste(root, dataN, '/', alg, '_go/', sep='')
        go_csv(csv.name, out.dir)
    }
break
}