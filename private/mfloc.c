#include <stdlib.h>
#include <string.h>
#include <stdio.h>
/* #include "matrix.h" */
#include "flochelp.h"
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  double *data = mxGetPr(prhs[0]);
  int nrowData = (int)mxGetScalar(prhs[1]);
  int ncolData = (int)mxGetScalar(prhs[2]);
  double *bicRow  = mxGetPr(prhs[3]);
  double *bicCol  = mxGetPr(prhs[4]);
  double *vecResvolBic = mxGetPr(prhs[5]);
  double *r    = mxGetPr(prhs[6]);
  int k        = (int) mxGetScalar(prhs[7]);
  int N        = (int) mxGetScalar(prhs[8]);
  int M        = (int) mxGetScalar(prhs[9]);
  int t        = (int) mxGetScalar(prhs[10]);
  int *vecBlocGene = (int *) mxGetPr(prhs[11]);
  int *vecBlocSample = (int *) mxGetPr(prhs[12]);

  int zz, i, valbreak, j;
  int total = nrowData + ncolData;
  double inv2R;

  double *bicRow2 = (double *) mxCalloc(k*nrowData, sizeof(double));
  double *bicCol2 = (double *) mxCalloc(k*ncolData, sizeof(double));
  int *vecOrder = (int *) mxCalloc(total, sizeof(int));
  double *vecBestGain = (double *) mxCalloc(total, sizeof(double));
  int *vecBestBic = (int *) mxCalloc(total, sizeof(int));
  double *sumRow = (double *) mxCalloc(k*nrowData, sizeof(double));
  double *sumCol = (double *) mxCalloc(k*ncolData, sizeof(double));
  double *sumBic = (double *) mxCalloc(k, sizeof(double));
  double *sumRow2 = (double *) mxCalloc(k*nrowData, sizeof(double));
  double *sumCol2 = (double *) mxCalloc(k*ncolData, sizeof(double));
  double *sumBic2 = (double *) mxCalloc(k, sizeof(double));

  double invk = 1/(double)k;

  memcpy(bicRow2, bicRow, k*nrowData*sizeof(double));
  memcpy(bicCol2, bicCol, k*ncolData*sizeof(double));

  for (zz=0; zz<k; zz++)
    {
      vecResvolBic[zz*4+2] = count_row_col(zz, nrowData, bicRow);
      vecResvolBic[zz*4+3] = count_row_col(zz, ncolData, bicCol);
      sum(zz, nrowData, ncolData, data, bicRow, bicCol, sumBic, sumRow, sumCol);
      sum(zz, nrowData, ncolData, data, bicRow, bicCol, sumBic2, sumRow2, sumCol2);
      vecResvolBic[zz*4] = residue(zz, nrowData, ncolData, data, bicRow, bicCol, sumBic, sumRow, sumCol);
    }
  j = 0;

  for (zz=0; zz<t; zz++)
    {
      valbreak = 0;

      bestgain(k, *r, nrowData, ncolData, data, bicRow, bicCol,bicRow2, bicCol2, sumBic, sumRow, sumCol, sumBic2, sumRow2, sumCol2, vecBestGain, vecBestBic, &inv2R, vecResvolBic, N, M, vecBlocGene, vecBlocSample);

      for (i=0; i<total; i++) vecOrder[i] = i;

      tri(vecBestGain, vecOrder, 0, total-1);
      order(&inv2R, nrowData, ncolData, vecBestGain, vecOrder);

      action(k, nrowData, ncolData, data, vecOrder, vecBestBic, bicRow, bicCol, bicRow2, bicCol2, *r, &valbreak, vecResvolBic, sumBic, sumRow, sumCol, sumBic2, sumRow2, sumCol2, N, M, zz, &j, vecBlocGene, vecBlocSample);

      if (valbreak==0)
        {
          break;
        }
    }
}
