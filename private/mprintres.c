#include <stdlib.h>
#include <stdio.h>
/* #include "matrix.h" */
#include "flochelp.h"
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
  if (nrhs!=5){
    mexErrMsgTxt("Wrong number of inputs");
  }
  int nrowData = mxGetScalar(prhs[0]);
  int ncolData = mxGetScalar(prhs[1]);
  double * data = mxGetPr(prhs[2]);
  double *bicRow =  mxGetPr(prhs[3]); 
  double *bicCol =  mxGetPr(prhs[4]);
  
  /* Result to be returned */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  double *res = mxGetPr(plhs[0]);
  
  double *sumRow = (double *) mxCalloc(nrowData, sizeof(double));
  double *sumCol = (double *) mxCalloc(ncolData, sizeof(double));
  double *sumBic = (double *) mxCalloc(1, sizeof(double));
    
  sum(0, nrowData, ncolData, data, bicRow, bicCol, sumBic, sumRow, sumCol);
  *res = residue(0, nrowData, ncolData, data, bicRow, bicCol, sumBic, sumRow, sumCol);
}
