#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include "flochelp.h"
#include "mex.h"

/**************/
/* A version of Marsaglia-MultiCarry */

static unsigned int I1=1234, I2=5678;

void set_seed(unsigned int i1, unsigned int i2)
{
  I1 = i1; I2 = i2;
}

void get_seed(unsigned int *i1, unsigned int *i2)
{
  *i1 = I1; *i2 = I2;
}


double unif_rand(void)
{
  I1= 36969*(I1 & 0177777) + (I1>>16);
  I2= 18000*(I2 & 0177777) + (I2>>16);
  return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
}

/*****************/

void echange (double *t, int *tab, int i, int j)   
{
  int tampon1 = tab[i];
  double tampon2 = t[i];
  tab[i] = tab[j];
  tab[j] = tampon1;
  t[i] = t[j];
  t[j] = tampon2;
}

void tri(double *t, int *tab, int G, int D)    
{
  int g,d;
  double val;
    
  if(D <= G)
    return;
  val = t[D];
  g = G - 1;
  d = D;
  do
    {
	    while ((t[++g]>val));
	    while((t[--d]<val) && (d>G));
	    
	    if (g < d) echange(t, tab, g, d);
    } while (g < d);
  echange(t, tab, g, D);
  tri(t, tab, G, g-1);
  tri(t, tab, g+1, D);
}

/**************************************************************************/


int count_row_col(int z, int nrowColData, double *bicRowCol)
{
  int N = 0, i;
  for (i = 0 ; i < nrowColData ; i++) N += (int)bicRowCol[nrowColData * z + i];
  return(N);
}

/***************************************************************************/

void sum(int z, int nrowData, int ncolData, double *Data, double *bicRow, double *bicCol, double *sumBic, double *sumRow, double *sumCol)
{
  int i, j, first = 0;
    
  sumBic[z] = 0;                                                  
    
  for (i = 0 ; i < nrowData ; i++)
    {
	    if ((int)bicRow[nrowData * z + i])                             
        {
          sumRow[nrowData * z + i] = 0;   
          for (j = 0 ; j < ncolData ; j++)
            {
              if ((int)bicCol[ncolData * z + j])
                {
                  if (first == 0) sumCol[ncolData * z + j] = 0;  
                  sumRow[z * nrowData + i] += Data[ncolData * i + j];
                  sumCol[z * ncolData + j] += Data[ncolData * i + j];
                }
            }
          sumBic[z] += sumRow[z * nrowData + i];
          first = 1;
        }
    }
}

/********************************************************************/

double residue(int z, int nrowData, int ncolData, double *Data, double *bicRow, double *bicCol, double *sumBic, double *sumRow, double *sumCol)
{
  int i, j;
  int nrowBic = count_row_col(z, nrowData, bicRow);
  int ncolBic = count_row_col(z, ncolData, bicCol);
  double res = 0, result;
    
  double const invrow = 1/(double)nrowBic;
  double const invcol = 1/(double)ncolBic;
  double const invvol = 1/(double)(nrowBic * ncolBic);
  double const meanbic = sumBic[z] * invvol;
    
  for (i = 0 ; i < nrowData ; i++)
    {
	    if ((int)bicRow[nrowData * z + i])
        { 
          for (j = 0 ; j < ncolData ; j++)
            {
              if ((int)bicCol[ncolData * z + j])
                {
                  result = Data[ncolData * i + j]  - sumRow[z * nrowData + i] * invcol - sumCol[z * ncolData + j] * invrow + meanbic;   
                  res += result * result;
                }
            }
        }
    }
  res = res * invvol;
  return(res);
}

/*************************************************************************/

void bestgain(int k, double r, int nrowData, int ncolData, double *Data, double *bicRow, double *bicCol, double *bicRow2, double *bicCol2, double *sumBic, double *sumRow, double *sumCol, double *sumBic2, double *sumRow2, double *sumCol2, double *vecBestGain, int *vecBestBic, double *inv2R, double *vecResvolBic, int N, int M, int *vecBlocGene, int *vecBlocSample)
{
  int z, i, j, u;
  int nrowBic, ncolBic, volBic, volBic2;
  double invvol, invr2onres;
  double valGain, res1, res2;
  double gainMax = -DBL_MAX, gainMin = DBL_MAX;
    
  for (z = 0 ; z < k ; z++)              
    {
	    nrowBic = vecResvolBic[z * 4 + 2];
	    ncolBic = vecResvolBic[z * 4 + 3];
	    volBic = nrowBic * ncolBic;
	    invvol = 1 / (double)volBic;
	    
	    res1 = vecResvolBic[z * 4];
	    invr2onres = 1 / (r * r / res1);
	    
	    for (i = 0 ; i < nrowData ; i++) 
        {
          if (!vecBlocGene[nrowData * z + i]){
            bicRow2[nrowData * z + i] = 1 - bicRow2[nrowData * z + i];
            u = 2 * (int)bicRow2[nrowData * z + i] - 1; 
            if ((nrowBic + u) >= N)
              { 
                sumBic2[z] = sumBic[z]; 
                sumRow2[nrowData * z + i] = 0;
				
                for (j = 0 ; j < ncolData ; j++)
                  {
                    if ((int)bicCol[ncolData * z + j] == 1)
                      {
                        sumRow2[nrowData * z + i] += Data[ncolData * i + j];
                        sumCol2[ncolData * z + j] = sumCol[ncolData * z +j] + u * Data[ncolData * i +j];
                      } 
                  }
                sumBic2[z] += u * sumRow2[nrowData * z + i];
                volBic2 = (nrowBic + u) * ncolBic;
                res2 = residue(z, nrowData, ncolData, Data, bicRow2, bicCol, sumBic2, sumRow2, sumCol2);
                valGain = (res1 - res2) * invr2onres + (volBic2 - volBic) * invvol;
              }
            else valGain = -DBL_MAX;
          }
		    
          else valGain = -DBL_MAX;

          if (z == 0)
            {
              vecBestGain[i] = valGain;
              vecBestBic[i] = 0;
            }
          else if (valGain >= vecBestGain[i])
            {
              vecBestGain[i] = valGain;
              vecBestBic[i] = z;
            }
		    
          if (vecBestGain[i] > gainMax) gainMax = vecBestGain[i];
          if (vecBestGain[i] < gainMin) gainMin = vecBestGain[i];
          bicRow2[nrowData * z + i] = bicRow[nrowData * z + i];
        }
	    
	    for (j = 0 ; j < ncolData ; j++)
        {    
          if ((int)bicCol[ncolData * z + j] == 1) sumCol2[ncolData * z + j] = sumCol[ncolData * z + j]; 
        }
	    sum(z, nrowData, ncolData, Data, bicRow, bicCol, sumBic2, sumRow2, sumCol2);
	    
	    for (j = 0; j < ncolData; j++)
        {
          if (!vecBlocSample[ncolData *z + j]){
            bicCol2[ncolData * z + j] = 1 - bicCol2[ncolData * z + j];
            u = 2 * (int)bicCol2[ncolData * z + j] - 1;
			
            if ((ncolBic + u) >= M)
              {
                sumBic2[z] = sumBic[z];
                sumCol2[ncolData * z + j] = 0;
                for (i = 0; i < nrowData; i++)
                  {
                    if ((int)bicRow[nrowData * z + i] == 1)
                      {
                        sumCol2[ncolData * z + j] += Data[ncolData * i + j];
                        sumRow2[nrowData * z + i] = sumRow[nrowData * z + i] + u * Data[ncolData * i + j];
                      }
                  }
                sumBic2[z] += u * sumCol2[ncolData * z + j];
                volBic2 = nrowBic * (ncolBic + u);
                res2 = residue(z, nrowData, ncolData, Data, bicRow, bicCol2, sumBic2, sumRow2, sumCol2);
                valGain = (res1 - res2) * invr2onres + (volBic2 - volBic) * invvol;
              }
            else valGain = -DBL_MAX;
          }
          else valGain= -DBL_MAX;
		    
          if (z == 0)
            {
              vecBestGain[j + nrowData] = valGain;
              vecBestBic[j + nrowData] = 0;
            }
          else if (valGain >= vecBestGain[j + nrowData])
            {
              vecBestGain[j + nrowData] = valGain;
              vecBestBic[j + nrowData] = z;
            }
          if (vecBestGain[j + nrowData] > gainMax) gainMax = vecBestGain[i];
          if (vecBestGain[j + nrowData] < gainMin) gainMin = vecBestGain[i];
          bicCol2[ncolData * z +j] = bicCol[ncolData * z +j]; 
		    
        }

	    for (i = 0 ; i < nrowData ; i++)
        {
          if ((int)bicRow[nrowData * z + i] == 1) sumRow2[nrowData * z + i] = sumRow[nrowData * z + i];
        }
	    
	    vecResvolBic[z * 4] = res1;
	    vecResvolBic[z * 4 + 1] = volBic;
	    vecResvolBic[z * 4 + 2] = nrowBic;
	    vecResvolBic[z * 4 + 3] = ncolBic;
    }
  *inv2R = 1/(gainMax - gainMin);
}

/*****************************************************************************/

void order(double *inv2R,int nrowData,int ncolData,double *vecBestGain, int *vecOrder)
{
  int i, rand1, rand2, rand1b, rand2b; 
  double proba, pij;
    
  for (i = 0 ; i < 2 * (nrowData + ncolData) ; i++) 
    {
	    rand1 = (int)((float)unif_rand() / RAND_MAX *(nrowData + ncolData -1)); 
	    rand2 = (int)((float)unif_rand() / RAND_MAX *(nrowData + ncolData -1));
	    rand1b = vecOrder[rand1];
	    rand2b = vecOrder[rand2];
	    
	    pij = 0.5 + (vecBestGain[rand2b] - vecBestGain[rand1b]) * *inv2R;
	    
	    proba = (double)((float)unif_rand() / RAND_MAX );
	    
	    if (pij >= proba)  
        {
          vecOrder[rand1] = rand2b;
          vecOrder[rand2] = rand1b;
        }
    } 
}

/*****************************************************************************/

void action(int k, int nrowData, int ncolData, double *Data, int *vecOrder, int *vecBestBic, double *bicRow, double *bicCol, double *bicRow2, double *bicCol2, double r, int *valbreak, double *vecResvolBic, double *sumBic, double *sumRow, double *sumCol, double *sumBic2, double *sumRow2, double *sumCol2, int N, int M, int zz, int *aa, int *vecBlocGene, int *vecBlocSample)
{
  int x, b, z, i, j, u, nrowBic, ncolBic, nrowBic2, ncolBic2, volBic, volBic2, amelio, contrainte;
  double res1, res2;
    
  for (x = 0 ; x < (nrowData + ncolData) ; x++)  
    {
	    b = vecOrder[x]; 
	    z = vecBestBic[b]; 
	    amelio = contrainte = 0;
	    
	    res1 = vecResvolBic[z * 4];
	    nrowBic = vecResvolBic[z * 4 + 2];
	    ncolBic = vecResvolBic[z * 4 + 3];
	    volBic = vecResvolBic[z * 4 + 1];
	    
	    sumBic2[z] = sumBic[z];
	    
	    if (b < nrowData) 
        {
          i = b;
          bicRow2[nrowData * z + i] = 1 - bicRow[nrowData * z + i];
          u = 2 * (int)bicRow2[nrowData * z + i] - 1;
          nrowBic2 = nrowBic + u;
          ncolBic2 = ncolBic;
		    
          if(nrowBic2 >= N && !vecBlocGene[nrowData * z + i]) 
            {  
              contrainte = 1;
              sumRow2[nrowData * z + i] = 0;
              for (j = 0 ; j < ncolData ; j++)
                {
                  if ((int)bicCol[ncolData * z + j] == 1)
                    {
                      sumRow2[nrowData * z + i] += Data[ncolData * i + j];
                      sumCol2[ncolData * z + j] += u * Data[ncolData * i + j];
                    }
                }
              sumBic2[z] += u * sumRow2[nrowData * z + i];
            }
        }
	    else  
        {
          j = b - nrowData;
          bicCol2[ncolData * z + j] = 1 - bicCol[ncolData * z + j];
          u = 2 * (int)bicCol2[ncolData * z + j] - 1;
          ncolBic2 = ncolBic + u;
          nrowBic2 = nrowBic;
		    
          if (ncolBic2 >= M && !vecBlocSample[ncolData * z + j])
            {
              contrainte = 1;
              sumCol2[ncolData * z + j] = 0;
              for (i = 0; i < nrowData; i++)
                {
                  if ((int)bicRow[nrowData * z + i] == 1)
                    {
                      sumCol2[ncolData * z + j] += Data[ncolData * i + j];
                      sumRow2[nrowData * z + i] += u * Data[ncolData * i + j];
                    }
                }
              sumBic2[z] += u * sumCol2[ncolData * z + j];
            }
        }
	    
	    if (contrainte)
        {
          res2 = residue(z, nrowData, ncolData, Data, bicRow2, bicCol2, sumBic2, sumRow2, sumCol2); 
          volBic2 = nrowBic2 * ncolBic2;
		    
          if ((res2 < res1) || ((res2 < r) && (volBic2 > volBic))) 
            {
              *aa += 1;
			    
              amelio = 1;           
              vecResvolBic[z * 4] = res2;
              vecResvolBic[z * 4 + 1] = volBic2;
              vecResvolBic[z * 4 + 2] = nrowBic2;
              vecResvolBic[z * 4 + 3] = ncolBic2;
              *valbreak = 1;
              sumBic[z] = sumBic2[z];
			    
              if (b < nrowData)  
                {
                  bicRow[nrowData * z + i] = bicRow2[nrowData * z + i];
                  sumRow[nrowData * z + i] = sumRow2[nrowData * z + i];
                  for (j = 0 ; j < ncolData ; j++)
                    {
                      if ((int)bicCol[ncolData * z + j] == 1)
                        {
                          sumCol[ncolData * z + j] = sumCol2[ncolData * z + j];
                        }
                    }
                }
              else 
                {
                  bicCol[ncolData * z + j] = bicCol2[ncolData * z +j];
                  sumCol[ncolData * z + j] = sumCol2[ncolData * z + j];
                  for(i = 0 ; i < nrowData ; i++)
                    {
                      if ((int)bicRow[nrowData * z + i] == 1)
                        {
                          sumRow[nrowData * z + i] = sumRow2[nrowData * z + i];
                        }
                    }
                }
            }
        }	     
	    if(!amelio)
        {
          if (b < nrowData)
            {
              bicRow2[nrowData * z + i] = bicRow[nrowData * z + i];
              for (j=0; j<ncolData; j++)
                {
                  if ((int)bicCol[ncolData * z + j] == 1)
                    {
                      sumCol2[ncolData * z + j] = sumCol[ncolData * z + j];
                    }
                }
            }
          else
            {
              bicCol2[ncolData * z + j] = bicCol[ncolData * z +j];
              for(i = 0 ; i < nrowData ; i++)
                {
                  if ((int)bicRow[nrowData * z + i] == 1)
                    {
                      sumRow2[nrowData * z + i] = sumRow[nrowData * z + i];
                    }
                }
            }
        }
    }
}
