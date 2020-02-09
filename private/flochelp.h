#ifndef FLOCHELP_H_INCLUDED
#define FLOCHELP_H_INCLUDED

void echange(double *t, int *tab, int i, int j);

void tri(double *t, int *tab, int G, int D);

int count_row_col(int z, int nrowColData, double *bicRowCol);

void sum(int z, int nrowData, int ncolData, double *Data, double *bicRow, double *bicCol, double *sumBic, double *sumRow, double *sumCol);

double residue(int z, int nrowData, int ncolData, double *Data, double *bicRow, double *bicCol, double *sumBic, double *sumRow, double *sumCol);

void bestgain(int k, double r, int nrowData, int ncolData, double *Data, double *bicRow, double *bicCol, double *bicRow2, double *bicCol2, double *sumBic, double *sumRow, double *sumCol, double *sumBic2, double *sumRow2, double *sumCol2, double *vecBestGain, int *vecBestBic, double *inv2R, double *vecResvolBic, int N, int M, int *vecBlocGene, int *vecBlocSample);

void order(double *inv2R,int nrowData,int ncolData,double *vecBestGain, int *vecOrder);

void action(int k, int nrowData, int ncolData, double *Data, int *vecOrder, int *vecBestBic, double *bicRow, double *bicCol, double *bicRow2, double *bicCol2, double r, int *valbreak, double *vecResvolBic, double *sumBic, double *sumRow, double *sumCol, double *sumBic2, double *sumRow2, double *sumCol2, int N, int M, int zz, int *aa, int *vecBlocGene, int *vecBlocSample);

#endif
