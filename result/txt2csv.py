import os
import pandas as pd
import numpy as np
import glob

datas = ['BCLL', 'PBC', 'RAT', 'YC']
columns = ['MSR', 'GV', 'CV', 'Var', 'Genes', 'Conditions']
ident_dict = {}
for data in datas:
    dataset_csv = glob.glob(os.path.join('/home/jiaobendaye/Documents/lab/fanzhenhao/biye/experiment/data', '*' + data + '_norm.csv'))
    df = pd.read_csv(dataset_csv[0])
    if 'Gene Symbol' in df.columns:
        ident_dict[data] = df['Gene Symbol'].to_numpy()
    else:
        ident_dict[data] = df['IDENTIFIER'].to_numpy()

def txt2csv(txt:str):
    if not os.path.exists(txt): return 
    csvfile = txt.replace('.txt', '.csv')
    items = []
    # read txt
    with open(txt, 'r') as fid:
        lines = fid.readlines()
        bic_cnt = len(lines) // 6
        #msr, gv, cv, var, genes, condis
        for _ in range(bic_cnt):
            msr = lines.pop(0).strip()
            gv = lines.pop(0).strip()
            cv = lines.pop(0).strip()
            var = lines.pop(0).strip()
            genes = ident_dict[data][list(map(int,lines.pop(0).strip().split()))]
            genes_str = '+'.join(genes) 
            condis = lines.pop(0).strip()
            items.append([msr, gv, cv, var, genes_str, condis])
    fid.close()
    #write csv
    df = pd.DataFrame(items)
    df.columns = columns
    df.to_csv(csvfile, index=None)

if __name__ == "__main__":
    root = './'
    for data in datas:
        txt = root + data + '.txt'
        print(txt)
        txt2csv(txt)