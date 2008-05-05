//Changes
//
//	12/18/07	add imputation() for Hctrl/Hcase with nsam cases/ctronls and nsnp SNPs with Roberts et al 2007 approach


/* filename: powerstudy.h
 * content: header file for the program "powerstudy", to study power of association testings.
 * Note: adapted from as-sim.h, which is for simulating sample using Durrant Alogrithm
 */
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* global variables */
double *scalepos; //positions scaled by recombintion between them
float DAF; //Disease Allele Freq.
int pstart; //start position in Durrant Algorithm
//extern float AF;
extern int nimpute, ncorrect;

struct hapinfo{
  int nseg; //number of segregating sites
  int nhap; //number of haplotypes
  float maf[100]; //distri. of MAFs on a grid of values, by = 2.5%
  float pdiff[100]; //distri. of pairwise differences, by = 3 bp
  float dprime[100]; //distr. of D' decay, by = 5 kb
  float rsquare[100]; //distr. of r^2 decay, by = 5 kb
  float corr[20][20]; //correlation


};


/* function prototypes */

void recomb_sim(char **simarr, int nsam, int nsnp, char **haparr, int nhap);
//void shuffle_sim(char **simarr, int nsam, int nsnp, char **haparr, int nhap);
void shuffleN_sim(char **simarr, int nsam, int nsnp, char **haparr, int nhap, int wsize);

int findDA(double *maf, int nmarkers, double daf);
void simgeno(char geno[], double daf);
//void simhap(char **simarr, int index, int snpindex, char phase, double **freq, int types, int centers);
void simhapN(char **simarr, int index, int snpindex, char phase, double **freq, int types, int centers, int wsize);





/* function definitions */


/* function: imputation()
 *	Icase and Ictrl store the SNP info (0s and 1s)
 *	Icase_snpind = 0 for missing and 1 for present; same for Ictrl_snpind
 *	Hcase, Hctrl, case_snpind, and ctrl_snpind store imputated SNP calls and missing/present info
 * note: 
 *      You can focus on the code of Case;
 *      Check the original bioinformatics paper for the hashing trick (not implemented here);
 */
float imputation(char **Icase, char **Ictrl, int **Icase_snpind, int **Ictrl_snpind, int nsam, int nsnp, char **Hcase, char **Hctrl, int **case_snpind, int **ctrl_snpind)
{
  //a base call is missing for the j-th SNP in the i-th chromosome when xxxx_snpind[i][j] = 0
  int i, j, k, ii, jj,  **imatrix();
  
  int **MVA; //mismatch vector array
  int n = nsam; //first we consider cases/controls separately
  int *z, *o, *q, *w; //temporary mismatch vectors
  int ncol, min, start, len, low, high, L;
  //int  nimpute, ncorrect;
  int tmpstart, tmplen;
  char base;
  char **sample, **result;
  int **snpind;
  float accuracy;
  
  ncol = n*(n-1)/2;	
  MVA = imatrix(nsnp+1, ncol);
  MALLOC(z, sizeof(int)*n, int *);
  MALLOC(o, sizeof(int)*n, int *);
  MALLOC(q, sizeof(int)*n, int *);
  MALLOC(w, sizeof(int)*ncol, int *);
  
  /////////  (I) Cases   ////////////
  
  //initialization
  sample = Hcase;
  snpind = case_snpind;
  result = Icase;
  
  //generate Mismatch Vectors
  for(k = 0; k < ncol; ++k)
    MVA[0][k] = 0;
  
  for(j = 0; j < nsnp; ++j)
    {
      for(i = 0; i < n; ++i)
	{
	  if (snpind[i][j])
	    {
	      if (sample[i][j] == '0')
		{
		  z[i] = 0;
		  o[i] = 2;
		} else
		  {
		    z[i] = 2;
		    o[i] = 0;
		  }
	    }
	  else
	    {
	      //missing 
	      z[i] = 1;
	      o[i] = 1;
	    }
	  
	  q[i] = 1;
	}
      
      start = 0;
      len = 0;
      for(i = 0; i < n-1; ++i)
	{
	  start += len;
	  len = n-i-1;
	  
	  if (snpind[i][j])
	    {
	      if (sample[i][j] == '0')
		for(k = 0; k < len; ++k)
		  MVA[j+1][start+k] = z[i+1+k];
	      else
		for(k = 0; k < len; ++k)
		  MVA[j+1][start+k] = o[i+1+k];
	    }
	  else
	    {
	      for(k = 0; k < len; ++k)
		MVA[j+1][start+k] = q[i+1+k];
	    }
	  
	  
	}
    }	
  
  
  //construct Mismatch vector array
  for(j = 0; j < nsnp; ++j)
    {
      for(k = 0; k < ncol; ++k)
	MVA[j+1][k] = MVA[j+1][k] + MVA[j][k];
    }
  
  L = 50;
  nimpute = ncorrect = 0;
  for(j = 0; j < nsnp; ++j)
    {
      //compute W(j, L)
      low = j - L;
      if (low < 0)
	low = 0;
      high = j + L+1;
      if (high > nsnp)
	high = nsnp;
      
      for(k = 0; k < ncol; ++k)
	w[k] = MVA[high][k] - MVA[low][k];
      
      start = 0;
      len = 0;
      for(i = 0; i < n; ++i)
	{
	  start += len;
	  len = n-i-1;
	  if (snpind[i][j] == 0)  //impute
	    {
	      //find the "closest" chromosome with known allele calls
	      min = 2*L;
	      base = '2';
	      for(k = 0; k < len; ++k)
		if (w[start+k] < min && snpind[i+k+1][j])
		  {
		    min = w[start+k];
		    base = sample[i+k+1][j];
		  }		
	      
	      //
	      tmpstart = 0; tmplen = 0;
	      for(ii = 0; ii < i; ++ii)
		{
		  tmpstart += tmplen;
		  tmplen = n - ii - 1;
		  if (w[tmpstart+i-ii-1] < min && snpind[ii][j])
		    {
		      min = w[tmpstart+i-ii-1];
		      base = sample[ii][j];
		    }
		}
	      
	      if (base == '0' || base == '1')
		{
		  result[i][j] = base;
		  nimpute++;
		  Icase_snpind[i][j] = 1;
		  
		  if (base == sample[i][j])
		    ncorrect++;
		}
	      else
		Icase_snpind[i][j] = 0;
	    }
	  else
	    {
	      result[i][j] = sample[i][j];
	      Icase_snpind[i][j] = 1;
	    }
	}
    }
  
  //obtain accuracy estimator
  accuracy = ncorrect*1.0/nimpute;
  //printf("total = %d\t imputed = %d\t corrct = %d\n", nsnp*n, nimpute, ncorrect);
  
  
  
  /////////  (II) Controls   ////////////
  
  //initialization
  sample = Hctrl;
  snpind = ctrl_snpind;
  result = Ictrl;
  
  //generate Mismatch Vectors
  for(k = 0; k < ncol; ++k)
    MVA[0][k] = 0;
  
  for(j = 0; j < nsnp; ++j)
    {
      for(i = 0; i < n; ++i)
	{
	  if (snpind[i][j])
	    {
	      if (sample[i][j] == '0')
		{
		  z[i] = 0;
		  o[i] = 2;
		} else
		  {
		    z[i] = 2;
		    o[i] = 0;
		  }
	    }
	  else
	    {
	      //missing 
	      z[i] = 1;
	      o[i] = 1;
	    }
	  
	  q[i] = 1;
	}
      
      start = 0;
      len = 0;
      for(i = 0; i < n-1; ++i)
	{
	  start += len;
	  len = n-i-1;
	  
	  if (snpind[i][j])
	    {
	      if (sample[i][j] == '0')
		for(k = 0; k < len; ++k)
		  MVA[j+1][start+k] = z[i+1+k];
	      else
		for(k = 0; k < len; ++k)
		  MVA[j+1][start+k] = o[i+1+k];
	    }
	  else
	    {
	      for(k = 0; k < len; ++k)
		MVA[j+1][start+k] = q[i+1+k];
	    }
	  
	  
	}
    }	
  
  
  //construct Mismatch vector array
  for(j = 0; j < nsnp; ++j)
    {
      for(k = 0; k < ncol; ++k)
	MVA[j+1][k] = MVA[j+1][k] + MVA[j][k];
    }
  
  
  L = 50;
  //nimpute = ncorrect = 0;
  for(j = 0; j < nsnp; ++j)
    {
      //compute W(j, L)
      low = j - L;
      if (low < 0)
	low = 0;
      high = j + L+1;
      if (high > nsnp)
	high = nsnp;
      
      for(k = 0; k < ncol; ++k)
	w[k] = MVA[high][k] - MVA[low][k];
      
      start = 0;
      len = 0;
      for(i = 0; i < n; ++i)
	{
	  start += len;
	  len = n-i-1;
	  if (snpind[i][j] == 0)  //impute
	    {
	      //find the "closest" chromosome with known allele calls
	      min = 2*L;
	      base = '2';
	      for(k = 0; k < len; ++k)
		if (w[start+k] < min && snpind[i+k+1][j])
		  {
		    min = w[start+k];
		    base = sample[i+k+1][j];
		  }		
	      
	      //
	      tmpstart = 0; tmplen = 0;
	      for(ii = 0; ii < i; ++ii)
		{
		  tmpstart += tmplen;
		  tmplen = n - ii - 1;
		  if (w[tmpstart+i-ii-1] < min && snpind[ii][j])
		    {
		      min = w[tmpstart+i-ii-1];
		      base = sample[ii][j];
		    }
		}
	      
	      if (base == '0' || base == '1')
		{
		  result[i][j] = base;
		  nimpute++;
		  Ictrl_snpind[i][j] = 1;
		  
		  if (base == sample[i][j])
		    ncorrect++;
		}
	      else
		Ictrl_snpind[i][j] = 0;
	    }
	  else
	    {
	      result[i][j] = sample[i][j];
	      Ictrl_snpind[i][j] = 1;
	    }
	}
    }
  
  //obtain accuracy estimator
  accuracy = ncorrect*1.0/nimpute;
  //printf("total = %d\t imputed = %d\t corrct = %d\n", nsnp*n, nimpute, ncorrect);
  
  
  free(z);
  free(o);
  free(q);
  
  return (accuracy);
}

