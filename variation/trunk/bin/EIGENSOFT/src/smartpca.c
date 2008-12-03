#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "eigsubs.h"  
#include "egsubs.h"  

#define WVERSION   "7521" 
/** 
Simple eigenvector analysis
Options to look at groups (simple ANOVA)
Weights allowed for individuals 
missing mode
dotpops added
recompiled with new twtail. Output form at changed
Cleaned up twestxx
fancynorm mode (divide by sqrt(p*(1-p)) 
poplistname supported.  Eigenanalysis just on individuals in population 
But all individuals figure in eigenvector output
New way of computing effective marker size  (twl2mode)
popdifference implemented
ldregression  ldlimit (genetic distance in Morgans)
nostatslim added  
dotpop has new format if many groups
uses new I/O
Supports packmode
Alkes style outlier removal added
Only half XTX computed
xdata (huge array) removed

fst calculation added 
popsizelimit added
divergence added (not useful?) 

SNPs discarded if no data.
Phylipfile now supported

Preparations for parallelization made
Various fixups for EIGENSTRAT and altnormstyle

output capability added (like convertf)

bug fixed (a last iteration needed for outlier removal)
bug fixed (numindivs unlimited)
output files fixed up (NULL OK)

Many Alkes style options added
Support for outliername added (outlier info)
familyname added (ped files)

bugfix: jackrat dies (outlier removes all of population)
bugfix: See ~ap92/REICH/FALL06/EIGENSOFTCODE/bugfix  bugs 1, 3 fixed

nrows, ncols output added
nrows, ncols set each outlier iteration
indivs with no data removed

writesnpeig added 

bugfix: popsize of 1 no anova done
minallelecnt added
chrom:  added
latest greatest handling of chromosome number added.
bad bugfix: numvalidgtypes 

checksizemode added.  Bug fix to version 7520 (> 2B genotypes fixed)

*/


#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL ;
char *twxtabname = NULL ;
char *trashdir = "/var/tmp" ;
int verbose = NO ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;

int numsnps, numindivs ; 
int numeigs = 10 ;  /// default
int markerscore = NO ;
int seed = 0 ;
int chisqmode = NO ;  // approx p-value better to use F-stat
int missingmode = NO ;
int plotmode = NO ;
int fancynorm = YES ;
int dotpopsmode = YES ;
int noxdata = YES ;  /* default as pop structure dubious if Males and females */
int pcorrmode = NO ;
int pcpopsonly = YES ;
int nostatslim = 10 ;
int znval = -1 ;
int popsizelimit = -1 ;
int altnormstyle = YES ;  // affects subtle details in normalization formula
int minallelecnt = 1 ;

int packout = -1 ;
extern enum outputmodetype outputmode  ;
extern int checksizemode ;
extern int packmode ;
int ogmode = NO ;

int numoutliter = 5, numoutleigs = 10 ;  
double outlthresh = 6.0 ; 
OUTLINFO  **outinfo ;
char *outinfoname = NULL ;  


double plo = .001 ;
double phi = .999 ;
double pvhit = .001 ;
double pvjack = 1.0e-6 ;
double *chitot ;
int    *xpopsize ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *poplistname = NULL ;
char *outliername = NULL ;
char *phylipname = NULL ;
char *snpeigname = NULL ;

char *indoutfilename = NULL ;
char *snpoutfilename = NULL ;
char  *genooutfilename = NULL ;
char  *omode = "packedancestrymap" ;

double r2thresh = -1.0 ;  
double r2genlim = 0.01 ; // Morgans 
double r2physlim = 5.0e6 ;     
int killr2 = NO ;  

int xchrom = -1 ;
// list of outliers

int ldregress = 0 ;
double ldlimit = 9999.0 ;  /* default is infinity */
/* we only consider markers as in possible LD if gdis <= ldlimit */

char *outputname = NULL ;
char *outputvname = NULL ;
char *weightname = NULL ;
FILE *ofile, *ovfile ;

double twestxx(double *lam, int m, double *pzn,  double *pzvar)  ;
double twnorm(double lam, double m, double n) ;

void readcommands(int argc, char **argv) ;
int loadindx(Indiv **xindlist, int *xindex, Indiv **indivmarkers, int numindivs) ;
int loadsnpx(SNP **xsnplist, SNP **snpmarkers, int numsnps, Indiv **indivmarkers) ;
void loadxdataind(double *xrow, SNP **snplist, int ind,  int ncols) ;           
void fixxrow(double *xrow, double *xmean, double *xfancy, int len)  ;
void dofancy(double *cc, int n, double *fancy) ;
int vadjust(double *rr, int n, double *pmean, double *fancy)  ;
void getcol(double *cc, double *xdata, int col, int nrows, int ncols)  ;
void getcolx(double *xcol, SNP *cupt, int *xindex, 
  int nrows, int col, double *xmean, double *xfancy)  ;          
void getcolxz(double *xcol, SNP *cupt, int *xindex, 
  int nrows, int col, double *xmean, double *xfancy, int *n0, int *n1)  ;          
void putcol(double *cc, double *xdata, int col, int nrows, int ncols)  ;
double dottest(char *sss, double *vec, char **eglist, int numeg, int *xtypes, int len) ;
double yll(double x1, double x2, double xlen) ;
void calcmean(double *wmean, double *vec, int len, int *xtypes, int numeg)  ;
double anova1(double *vec, int len, int *xtypes, int numeg) ;
double anova(double *vec, int len, int *xtypes, int numeg) ;
void publishit(char *sss, int df, double chi) ;

void setmiss(SNP **snpm, int numsnps)  ;
void setfvecs(double *fvecs, double *evecs, int nrows, int numeigs) ;
void dotpops(double *X, char **eglist, int numeg, int *xtypes, int nrows)  ;
void printxcorr(double *X, int nrows, Indiv **indxx)  ;      

void ldreg(double *ldmat, double *ldmat2, 
  double *vv, double *vv2, double *ldvv, 
  double *ldvv2, int rsize, int n)  ;

void clearld(double *ldmat, double *ldvv, int rsize, int n, int nclear)  ;

void
setidmat(double *mat, int n, double scal)  ; 

void fixrho(double *a, int n) ;
void printdiag(double *a, int n) ;

int
ridoutlier(double *evecs, int n, int neigs, 
 double thresh, int *badlist) ;

void addoutersym(double *X, double *v, int n)  ;
void symit(double *X, int n)  ;

double dofst(double *fstans, double *fstsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg) ;

double fstcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) ;

double fst(SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int type1, int type2, double *psd) ;

void
 setplimit(Indiv **indivmarkers, int numindivs, 
 char **eglist, int numeg, int plimit)  ;

void jackrat(double *xmean, double *xsd, double *top, double *bot,  int len)  ;
void domult(double  *tvecs, double  *tblock, int numrow, int len)  ;
void writesnpeigs(char *snpeigname, SNP **xsnplist, double *ffvecs, int numeigs, int ncols)  ;

int main(int argc, char **argv)
{

  char sss[MAXSTR] ;
  int **snppos ;
  int *snpindx ;
  char **snpnamelist, **indnamelist ;
  char **eglist ;
  int  lsnplist, lindlist, numeg ;
  int i, j, k, k1, k2; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;
  double y1, y2, y ;
  FILE *twxtestfp;

  int ch1, ch2 ;
  int fmnum , lmnum ;
  int num, n0, n1, nkill ;

  int nindiv = 0, e, f, lag=1  ;
  double xc[9], xd[4], xc2[9] ;
  double ychi,  tail, tw ;
  int nignore, numrisks = 1 ;
  double  *xrow, *xpt ; 
  SNP **xsnplist  ;
  Indiv **xindlist ;
  int *xindex, *xtypes ;
  int nrows, ncols, m ;
  double *XTX, *cc, *evecs, *ww, weight, *qvec, *qcoord ; 
  double *lambda ;
  double *tvecs ;
  double zn, zvar ;
  double *fvecs, *fxvecs, *fxscal ;
  double *ffvecs ;
  int weightmode = NO ;
  double chisq, ynrows ;
  int *numhits, t ;  
  double *xmean, *xfancy ;
  double *ldmat, *ldmat2, *ldvv, *ldvv2, *vv2  ;
  double *fstans, *fstsd ; 
  int chrom,  numclear ;
  double gdis ;
  int outliter, numoutiter, *badlist, nbad ;
  int a, b, kmax, kmin ;
  FILE *outlfile, *phylipfile  ;
  double **eigmoment, **eigindmoment ;
  double *eigkurt, *eigindkurt ;
  double *snpsc ; 
  

  int xblock, blocksize=10000 ;   
  double *tblock ;  

  OUTLINFO *outpt ;

  readcommands(argc, argv) ;
  printf("## smartpca version: %s\n", WVERSION) ;
  packmode = YES ;
  setomode(&outputmode, omode) ;

  if (parname == NULL) return 0 ;

  if (fancynorm) printf("norm used\n\n") ;
  else printf("no norm used\n\n") ;

  nostatslim = MAX(nostatslim, 3) ;

  outlfile = ofile = stdout; 

  if (outputname != NULL)  openit(outputname, &ofile, "w") ;
  if (outliername != NULL) openit(outliername, &outlfile, "w") ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;
  k = getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;


  if (poplistname != NULL) 
  { 
    ZALLOC(eglist, numindivs, char *) ; 
    numeg = loadlist(eglist, poplistname) ;
    seteglist(indivmarkers, numindivs, poplistname);
  }
  else
  {
    setstatus(indivmarkers, numindivs, NULL) ;
    ZALLOC(eglist, MAXPOPS, char *) ;
    numeg = makeeglist(eglist, MAXPOPS, indivmarkers, numindivs) ;
  }
  for (i=0; i<numeg; i++) 
  {  
    /* printf("%3d %s\n",i, eglist[i]) ; */
  }

  nindiv=0 ;
  for (i=0; i<numindivs; i++) 
  {
    indx = indivmarkers[i] ;
    if(indx -> affstatus == YES) ++nindiv  ;
  }

  /** X not supported. */
  for (i=0; i<numsnps; i++)  
  {  
    cupt = snpmarkers[i] ; 
    chrom = cupt -> chrom ;
    if ((noxdata) && (chrom == 23)) cupt-> ignore = YES ;
    if (chrom == 0) cupt -> ignore = YES ;
    if (chrom > 23) cupt -> ignore = YES ;
  }
  for (i=0; i<numsnps; i++)  
  {
    cupt = snpmarkers[i] ; 
    if ((xchrom>0) && (cupt -> chrom != xchrom)) cupt -> ignore = YES ;
    if (cupt -> ignore) continue ;
    if (numvalidgtx(indivmarkers, cupt, YES) <= 1) 
    { 
      printf("nodata: %20s\n", cupt -> ID) ;
      cupt -> ignore = YES ;
    }
  }

  if (killr2) {
   nkill = killhir2(snpmarkers, numsnps, numindivs, r2physlim, r2genlim, r2thresh) ;
   if (nkill>0) printf("killhir2.  number of snps killed: %d\n", nkill) ;
  }

  for (i=0; i<numindivs; ++i)  { 
   if (numvalidgtind(snpmarkers, numsnps, i) ==0) { 
    indx = indivmarkers[i] ;
    indx -> ignore = YES ; 
   }
  }

  numsnps = rmsnps(snpmarkers, numsnps) ;  //  rid ignorable snps
  if (missingmode) 
  {
    setmiss(snpmarkers, numsnps) ;
    fancynorm = NO ;
  }

  if  (weightname != NULL)   
  {  
    weightmode = YES ;
    getweights(weightname, snpmarkers, numsnps) ;
  }
  if (ldregress>0) 
  {  
    ZALLOC(ldvv,  ldregress*numindivs, double) ;
    ZALLOC(ldvv2,  ldregress*numindivs, double) ;
    ZALLOC(vv2,  numindivs, double) ;
    ZALLOC(ldmat,  ldregress*ldregress, double) ;
    ZALLOC(ldmat2,  ldregress*ldregress, double) ;
    setidmat(ldmat, ldregress, 1.0e-6) ;
  }

  ZALLOC(xindex, numindivs, int) ;
  ZALLOC(xindlist, numindivs, Indiv *) ;
  ZALLOC(xsnplist, numsnps, SNP *) ;

  if (popsizelimit > 0) 
  {  
    setplimit(indivmarkers, numindivs, eglist, numeg, popsizelimit) ; 
  }

  nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
  ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;
  printf("number of samples used: %d number of snps used: %d\n", nrows, ncols) ;

  /* printf("## nrows: %d  ncols  %d\n", nrows, ncols) ; */
  ZALLOC(xmean, ncols, double) ;
  ZALLOC(xfancy, ncols, double) ;
  ZALLOC(XTX, nrows*nrows, double) ;
  ZALLOC(evecs, nrows*nrows, double) ;
  ZALLOC(tvecs, nrows*nrows, double) ;
  ZALLOC(lambda, nrows, double) ;
  ZALLOC(cc, nrows, double) ;
  ZALLOC(ww, nrows, double) ;
  ZALLOC(badlist, nrows, int) ;

  blocksize = MIN(blocksize, ncols) ; 
  ZALLOC(tblock, nrows*blocksize, double) ;

  // xfancy is multiplier for column xmean is mean to take off
  // badlist is list of rows to delete (outlier removal) 

  numoutiter = 1 ;  

  if (numoutliter>=1) 
  {
    numoutiter = numoutliter+1 ;
    ZALLOC(outinfo, nrows,  OUTLINFO *) ;  
    for (k=0; k<nrows; k++) 
    {  
      ZALLOC(outinfo[k], 1, OUTLINFO) ;
    }
    /* fprintf(outlfile, "##%18s %4s %6s %9s\n", "ID", "iter","eigvec", "score") ; */
  }

  for (outliter = 1; outliter <= numoutiter ; ++outliter)  
  { 
    if (outliter>1) {
     ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;
    }
    vzero(XTX, nrows*nrows) ;
    vzero(tblock, nrows*blocksize) ;
    xblock = 0 ; 

    vzero(xmean, ncols) ;
    vclear(xfancy, 1.0, ncols) ;

    for (i=0; i<ncols; i++) 
    { 
      cupt = xsnplist[i] ;
      chrom = cupt -> chrom ;
      getcolxz(cc, cupt, xindex, nrows, i, xmean, xfancy, &n0, &n1) ;
      t = MIN(n0, n1) ; 

      if (t <= minallelecnt)  {  
       cupt -> ignore = YES ;
       vzero(cc, nrows) ; 
      }

      if (weightmode) 
      {
        weight = xsnplist[i] -> weight ;
        vst(cc, cc, xsnplist[i] -> weight, nrows) ;
      }
      if (ldregress>0) 
      {  
        numclear = 0 ;
        for (k=1; k<= ldregress; ++k)  
        {  
          j = i-k ;  
          if (j<0) 
          { 
            numclear = ldregress-k+1 ; 
            break ;
          }
          cupt2 = xsnplist[j] ;  
          if (cupt2 -> chrom != chrom) gdis = ldlimit + 1.0 ; 
          else gdis = cupt -> genpos - cupt2 -> genpos ;
          if (gdis>=ldlimit) 
          {   
            numclear = ldregress-k+1 ; 
            break ;
          }
        }
        if (numclear>0) clearld(ldmat, ldvv, ldregress, nrows, numclear) ; 
        ldreg(ldmat, ldmat2, cc, vv2, ldvv, ldvv2, ldregress, nrows) ;
        copyarr(ldmat2, ldmat, ldregress*ldregress) ;
        copyarr(vv2, cc, nrows) ;
        copyarr(ldvv2, ldvv, ldregress*nrows) ;
      }
      copyarr(cc, tblock+xblock*nrows, nrows) ;
      ++xblock ; 

/** this is the key code to parallelize */
      if (xblock==blocksize) 
      {  
        domult(tvecs, tblock, xblock, nrows) ;
        vvp(XTX, XTX, tvecs, nrows*nrows) ;
        xblock = 0 ;
        vzero(tblock, nrows*blocksize) ;
      }
    }

    if (xblock>0) 
    { 
     domult(tvecs, tblock, xblock, nrows) ;
     vvp(XTX, XTX, tvecs, nrows*nrows) ;
    }
    symit(XTX, nrows) ;

    /**
    a = 0; b=0 ;
    printf("zz1 %12.6f ", XTX[a*nrows+b]) ;
    a = nrows-1; b=nrows-1 ;
    printf(" %12.6f %15.9g\n", XTX[a*nrows+b], asum(XTX, nrows*nrows)) ;
    */

    if (verbose) 
    {
      printdiag(XTX, nrows) ;
    }

    y = trace(XTX, nrows) / (double) (nrows-1) ;
    if (isnan(y)) fatalx("bad XTX matrix\n") ;
    /* printf("trace:  %9.3f\n", y) ; */
    if (y<=0.0) fatalx("XTX has zero trace (no data??)\n") ;
    vst(XTX, XTX, 1.0/y, nrows * nrows) ;
    eigvecs(XTX, lambda, evecs, nrows) ;
// eigenvalues are in decreasing order 

    if (outliter > numoutliter) break ;  
    // last pass skips outliers 
    numoutleigs = MIN(numoutleigs, nrows-1) ;
    nbad = ridoutlier(evecs, nrows, numoutleigs, outlthresh, badlist) ;
    if (nbad == 0) break ; 
    for (i=0; i<nbad; i++) 
    {  
      j = badlist[i] ;
      indx = xindlist[j] ;
      outpt = outinfo[j] ;
      fprintf(outlfile, "REMOVED outlier %s iter %d evec %d sigmage %.3f\n", indx -> ID, outliter, outpt -> vecno, outpt -> score) ;
      indx -> ignore = YES ;
    }
    nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
    printf("number of samples after outlier removal: %d\n", nrows) ;
  }

  if (outliername != NULL) fclose(outlfile) ;

  m = numgtz(lambda, nrows)  ;
  /* printf("matrix rank: %d\n", m) ; */
  if (m==0) fatalx("no data\n") ;

  /* Now, print Tracy-Widom stats, if twtable is valid */
   if (settwxtable(twxtabname)<0) 
  {
    printf("\n## To get Tracy-Widom statistics: recompile smartpca with");
    printf(" TWTAB correctly specified in Makefile, or\n");
    printf("   just run twstats (see README file in POPGEN directory)\n");
  }
  else
  {
    /* *** START of code to print Tracy-Widom statistics */
    printf("\n## Tracy-Widom statistics: rows: %d  cols: %d\n", nrows, ncols);
    y = -1.0 ;
    printf("%4s  %12s", "#N", "eigenvalue") ;
    printf("%12s",  "difference") ;
    printf(" %9s %12s", "twstat", "p-value") ;
    printf(" %9s", "effect. n") ;
    printf("\n") ;

    ynrows = (double) nrows ;  

    for (i=0; i<m; ++i) 
    { 
      zn = znval ;
      if (zn>0) zn = MAX(ynrows, zn) ;
      tail = dotwcalc(lambda+i, m-i, &tw, &zn, &zvar, nostatslim) ;
      printf("%4d  %12.6f", i+1, lambda[i]) ;
      if (i==0) printf( "%12s", "NA") ;
      else printf("%12.6f", lambda[i]-lambda[i-1]) ;
      if (tail>=0.0) printf( " %9.3f %12.6g", tw, tail) ;
      else printf( " %9s %12s", "NA", "NA") ;
      if (zn>0.0) 
      {
        printf( " %9.3f", zn) ;
      }
      else 
      {
        printf( " %9s", "NA") ;
      } 
      printf( "\n") ;
    }
    /* END of code to print Tracy-Widom statistics */
  }

  numeigs = MIN(numeigs, nrows) ;
  numeigs = MIN(numeigs, ncols) ;

  /* fprintf(ofile, "##genotypes: %s\n", genotypename) ; */
  /* fprintf(ofile, "##numrows(indivs):: %d\n", nrows) ; */
  /* fprintf(ofile, "##numcols(snps):: %d\n", ncols) ; */
  /* fprintf(ofile, "##numeigs:: %d\n", numeigs) ; */
  fprintf(ofile, "%20s ", "#eigvals:") ;
  for (j=0; j<numeigs; j++) { 
	  fprintf(ofile, "%9.3f ", lambda[j]) ;
  }
  fprintf(ofile, "\n") ;

  if (outputvname != NULL)  {  
    openit(outputvname, &ovfile, "w") ;
    for (j=0; j<nrows; j++) { 
      fprintf(ovfile, "%12.6f\n", lambda[j]) ;
    }
    fclose(ovfile) ;
  }


  ZALLOC(xtypes, nrows, int) ;

  ZALLOC(fvecs, nrows*numeigs, double) ;
  ZALLOC(fxvecs, nrows*numeigs, double) ;
  ZALLOC(fxscal, numeigs, double) ;

  ZALLOC(ffvecs, ncols*numeigs, double) ;
  ZALLOC(xrow, ncols, double) ;
  setfvecs(fvecs, evecs, nrows, numeigs) ;

   for (i=0; i<ncols; i++)  { 
    cupt = xsnplist[i] ;
    getcolx(cc, cupt, xindex, nrows, i, NULL, NULL) ;
    for (j=0; j<numeigs; j++)  { 
     for (k=0; k<nrows; k++)  {
      ffvecs[j*ncols+i] += fvecs[j*nrows+k]*cc[k] ;
     }
    }
   }

  eigmoment    = initarray_2Ddouble(numeigs, 5, 0.0) ;
  eigindmoment = initarray_2Ddouble(numeigs, 5, 0.0) ;

  ZALLOC(eigkurt, numeigs, double) ;
  ZALLOC(eigindkurt, numeigs, double) ;

  for (j=0; j<numeigs; ++j) {  
   eigkurt[j] = kurtosis(ffvecs+j*ncols, ncols) ;
   eigindkurt[j] = kurtosis(fvecs+j*nrows, nrows) ;
  }

  for (i=0; i<nrows; i++) {

   indx = xindlist[i] ;
   k = indxindex(eglist, numeg, indx -> egroup) ;
   xtypes[i] = k ;

   loadxdataind(xrow, xsnplist, xindex[i],  ncols) ;            
   fixxrow(xrow, xmean, xfancy, ncols) ;

   for (j=0; j<numeigs; j++) { 

     xpt = ffvecs+j*ncols ;
     y = fxvecs[j*nrows+i] = vdot(xrow, xpt, ncols) ;
     fxscal[j] += y*y ;
   
   }
  }

  for (j=0; j<numeigs; j++) { 
    y = fxscal[j] ;                                            
//  fxscal[j] = 10.0/sqrt(y) ; // eigenvectors have norm 10 (perhaps eccentric)
    fxscal[j] = 1.0/sqrt(y) ;  // standard
  }

  for (i=0; i < numindivs ; i++)  { 
          indx = indivmarkers[i] ;
          if (indx -> ignore) continue ;
          loadxdataind(xrow, xsnplist, i,  ncols) ;            
          fixxrow(xrow, xmean, xfancy, ncols) ;
	  fprintf(ofile, "%20s ", indx -> ID) ;
	  for (j=0; j<numeigs; j++) { 
           y = fxscal[j]*vdot(xrow, ffvecs+j*ncols, ncols) ;


           fprintf(ofile, "%10.4f  ", y) ;
	  }
	  fprintf(ofile, "%15s\n", indx -> egroup) ;
  }

 printf("%12s %4s %9s %9s\n", "kurtosis", "", "snps", "indivs") ;

  for (j=0; j<numeigs; ++j)  {  
   y1 = eigkurt[j] ;              
   y2 = eigindkurt[j] ;              
   printf("%12s %4d %9.3f %9.3f\n", "eigenvector", j+1, y1, y2) ;
  }


// output files
  settersemode(YES) ;

  outfiles(snpoutfilename, indoutfilename, genooutfilename, 
   snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode) ;

  ZALLOC(xpopsize, numeg, int) ;
  for (i=0; i<nrows; i++) { 
    k = xtypes[i] ;
    ++xpopsize[k] ;
  }

  for (i=0; i<numeg; i++) 
  {  
    printf("population: %3d %s %4d\n",i, eglist[i], xpopsize[i]) ; 
  }


  if (numeg==1) dotpopsmode = NO ;

  if (dotpopsmode == NO)  {
    writesnpeigs(snpeigname, xsnplist, ffvecs, numeigs, ncols) ;
    printxcorr(XTX, nrows, xindlist) ;
    printf("##end of smartpca run\n") ;
    return 0 ;
  }

  ZALLOC(chitot, numeg*numeg, double) ;

  dotpops(XTX, eglist, numeg, xtypes, nrows) ;
  ZALLOC(fstans, numeg*numeg, double) ;
  ZALLOC(fstsd , numeg*numeg, double) ;

  dofst(fstans, fstsd, xsnplist, xindex, xtypes, nrows, ncols, numeg) ;

  if (phylipname == NULL) {
   printf("## Fst statistics between populations:         fst       std error\n");
   for (k1=0; k1<numeg; ++k1) {  
    for (k2=k1+1; k2<numeg; ++k2) {  
     printf(" %20s %20s %9.3f %10.4f\n", eglist[k1], eglist[k2], 
      fstans[k1*numeg+k2], fstsd[k1*numeg+k2]) ;
    }
   }
   printf("\n");
  }
  
  else {
   openit(phylipname, &phylipfile, "w") ;
   fprintf(phylipfile, "%6d\n",numeg) ; 
   sss[10] = CNULL ;
    for (k1=0; k1<numeg; ++k1) { 
     strncpy(sss, eglist[k1], 10) ;  
     fprintf(phylipfile, "%10s", sss) ;
     for (k2=0; k2<numeg; ++k2) {  
      y1 = fstans[k1*numeg+k2] ; 
      y2 = fstans[k2*numeg+k1] ; 
      fprintf(phylipfile, "%6.3f", (0.5*(y1+y2))) ;
     }
     fprintf(phylipfile, "\n") ;
    }
   fclose(phylipfile) ;
  }

  


  for (j=0; j< numeigs; j++)  {  
   sprintf(sss, "eigenvector %d", j+1) ;
   y=dottest(sss, evecs+j*nrows, eglist, numeg, xtypes, nrows) ;
  }

  printf("\n## Statistical significance of differences beween populations:\n");
  printf("                                pop1                  pop2      chisq          p-value   |pop1|   |pop2|\n");
  for (k1=0; k1<numeg; ++k1) {  
   for (k2=k1+1; k2<numeg; ++k2) {  
    ychi = chitot[k1*numeg+k2] ;
    tail = rtlchsq(numeigs, ychi) ;
    printf("popdifference:  %20s  %20s  %12.3f  %12.6g", eglist[k1], eglist[k2], ychi, tail) ;
    printf ("   %5d", xpopsize[k1]) ;
    printf ("   %5d", xpopsize[k2]) ;
    printf("\n") ;
   }
  }
  printf("\n");
  for (i=0; i<ncols; i++)  {  
   if (markerscore == NO) break;
   cupt = xsnplist[i] ;
   getcolx(cc, cupt, xindex, nrows, i, NULL, NULL) ;
   sprintf(sss, "%s raw", cupt -> ID) ;
   dottest(sss, cc, eglist, numeg, xtypes, nrows) ;
   for (j=0; j< numeigs; j++)  {  
    sprintf(sss, "%s subtract sing vec %d", cupt ->ID, j+1) ;
    y = vdot(cc, evecs+j*nrows, nrows) ;
    vst(ww, evecs+j*nrows, y, nrows) ;
    vvm(cc, cc, ww, nrows) ;
    dottest(sss, cc, eglist, numeg, xtypes, nrows) ;
   }
  }

  printxcorr(XTX, nrows, xindlist)  ; 


  writesnpeigs(snpeigname, xsnplist, ffvecs, numeigs, ncols) ;
  printf("##end of smartpca run\n") ;
  return 0 ;
}

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   if (parname==NULL) { 
    fprintf(stderr, "no parameters\n") ;
    return ;
   }

   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "snpeigname:", &snpeigname) ;
   getstring(ph, "snpweightoutname:", &snpeigname) ; /* changed 09/18/07 */
   getstring(ph, "output:", &outputname) ;
   getstring(ph, "outputvecs:", &outputname) ;
   getstring(ph, "evecoutname:", &outputname) ; /* changed 11/02/06 */
   getstring(ph, "outputvals:", &outputvname) ;
   getstring(ph, "evaloutname:", &outputvname) ; /* changed 11/02/06 */
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "outliername:", &outliername) ;
   getstring(ph, "outlieroutname:", &outliername) ; /* changed 11/02/06 */
   getstring(ph, "phylipname:", &phylipname) ;
   getstring(ph, "phylipoutname:", &phylipname) ; /* changed 11/02/06 */
   getstring(ph, "weightname:", &weightname) ;
   getint(ph, "numeigs:", &numeigs) ;
   getint(ph, "numoutevec:", &numeigs) ; /* changed 11/02/06 */
   getint(ph, "markerscore:", &markerscore) ; 
   getint(ph, "chisqmode:", &chisqmode) ; 
   getint(ph, "missingmode:", &missingmode) ; 
   getint(ph, "fancynorm:", &fancynorm) ; 
   getint(ph, "usenorm:", &fancynorm) ;  /* changed 11/02/06 */
   getint(ph, "dotpopsmode:", &dotpopsmode) ; 
   getint(ph, "pcorrmode:", &pcorrmode) ;  /* print correlations */
   getint(ph, "pcpopsonly:", &pcpopsonly) ;  /* but only within population */
   getint(ph, "altnormstyle:", &altnormstyle) ;  
   getint(ph, "hashcheck:", &hashcheck) ;
   getint(ph, "popgenmode:", &altnormstyle) ;  
   getint(ph, "noxdata:", &noxdata) ; 
   getint(ph, "nostatslim:", &nostatslim) ; 
   getint(ph, "popsizelimit:", &popsizelimit) ; 
   getint(ph, "minallelecnt:", &minallelecnt) ; 
   getint(ph, "chrom:", &xchrom) ;  
   getint(ph, "checksizemode:", &checksizemode) ;  

   getint(ph, "ldregress:", &ldregress) ;
   getint(ph, "nsnpldregress:", &ldregress) ; /* changed 11/02/06 */
   getdbl(ph, "ldlimit:", &ldlimit) ;  /* in morgans */
   getdbl(ph, "maxdistldregress:", &ldlimit) ;  /* in morgans */ /* changed 11/02/06 */
   getint(ph, "minleneig:", &nostatslim) ;
   getint(ph, "malexhet:", &malexhet) ;
   getint(ph, "nomalexhet:", &malexhet) ; /* changed 11/02/06 */
   getint(ph, "familynames:", &familynames) ;

   getint(ph, "numoutliter:", &numoutliter) ;
   getint(ph, "numoutlieriter:", &numoutliter) ; /* changed 11/02/06 */
   getint(ph, "numoutleigs", &numoutleigs) ;
   getint(ph, "numoutlierevec:", &numoutleigs) ; /* changed 11/02/06 */
   getdbl(ph, "outlthresh:", &outlthresh) ;
   getdbl(ph, "outliersigmathresh:", &outlthresh) ; /* changed 11/02/06 */

   getstring(ph, "indoutfilename:", &indoutfilename) ;
   getstring(ph, "indivoutname:", &indoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "snpoutfilename:", &snpoutfilename) ;
   getstring(ph, "snpoutname:", &snpoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "genooutfilename:", &genooutfilename) ;
   getstring(ph, "genotypeoutname:", &genooutfilename) ; /* changed 11/02/06 */
   getstring(ph, "outputformat:", &omode) ;
   getstring(ph, "outputmode:", &omode) ;
   getint(ph, "outputgroup:", &ogmode) ;
   getint(ph, "packout:", &packout) ; /* now obsolete 11/02/06 */
   getstring(ph, "twxtabname:", &twxtabname) ;

   getdbl(ph, "r2thresh:", &r2thresh) ;
   getdbl(ph, "r2genlim:", &r2genlim) ;
   getdbl(ph, "r2physlim:", &r2physlim) ;
   getint(ph, "killr2:",  &killr2) ;


   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);

}
int
loadindx(Indiv **xindlist, int *xindex, Indiv **indivmarkers, int numindivs) 
{
  int i, n=0 ;  
  Indiv *indx ;
  for (i=0; i<numindivs; i++)  {
   indx = indivmarkers[i] ;
   if (indx -> ignore) continue ;
   if (indx -> affstatus == NO) continue ;
   xindex[n] = i ;
   xindlist[n] = indx ;
   ++n ;
  }
  return n ;
}
int
loadsnpx(SNP **xsnplist, SNP **snpmarkers, int numsnps, Indiv **indivmarkers) 
// snps with <= 1 valid genotype forced to ignore
{
  int i, n=0, t ;  
  SNP *cupt ;
  for (i=0; i<numsnps; i++)  {
   cupt = snpmarkers[i] ;
   t = numvalidgtypes(cupt) ;
// printf("zzz %s %d\n", cupt -> ID, t) ;
   if (t<=1) cupt -> ignore = YES ;
   if (cupt -> ignore) continue ;
   xsnplist[n] = cupt ;
   ++n ;
  }
  return n ;
}

void
getrawcol(int *rawcol, SNP *cupt, int *xindex, int nrows) 
{
  int t, j ;
  for (j=0; j< nrows; j++) {  
   t = xindex[j] ;
   rawcol[j] = getgtypes(cupt, t) ;          
   if (verbose) printf("www %d %d %d\n", j, t, rawcol[j]) ;
  }
}

void
getcolxz(double *xcol, SNP *cupt, int *xindex, int nrows, int col,  
 double *xmean, double *xfancy, int *n0, int *n1)             
// side effect set xmean xfancy and count variant and reference alleles
{
 Indiv *indx  ;
 int  j,  n, g, t ;
 double y, pmean, p, yfancy ;
 int *rawcol ;
 int c0, c1 ;

  c0 = c1 = 0 ;
  ZALLOC(rawcol, nrows, int) ;
  n = cupt -> ngtypes ;
  if (n<nrows) fatalx("bad snp: %s %d\n", cupt -> ID, n) ;
  getrawcol(rawcol, cupt, xindex, nrows) ;
  for (j=0; j<nrows; ++j) { 
   g = rawcol[j] ;  
   if (g<0) continue ;  
   c0 += g   ;
   c1 += 2-g ;
  }
  floatit(xcol, rawcol, nrows) ;

  vadjust(xcol, nrows, &pmean, &yfancy) ;
  vst(xcol, xcol, yfancy, nrows) ;
  if (xmean != NULL) {
   xmean[col] = pmean*yfancy ; 
   xfancy[col] = yfancy ;
  }
  free(rawcol) ;
  *n0 = c0 ; 
  *n1 = c1 ;
}

void
getcolx(double *xcol, SNP *cupt, int *xindex, int nrows, int col,  
 double *xmean, double *xfancy)             
// side effect set xmean xfancy
{
 Indiv *indx  ;
 int  j,  n, g, t ;
 double y, pmean, p, yfancy ;
 int *rawcol ;

  if (xmean != NULL) {
   xmean[col] = xfancy[col] = 0.0 ;
  }
  
  if (cupt -> ignore) {  
   vzero(xcol, nrows) ;
   return ;
  }

  ZALLOC(rawcol, nrows, int) ;
  n = cupt -> ngtypes ;
  if (n<nrows) fatalx("bad snp: %s %d\n", cupt -> ID, n) ;
  getrawcol(rawcol, cupt, xindex, nrows) ;
  floatit(xcol, rawcol, nrows) ;

  vadjust(xcol, nrows, &pmean, &yfancy) ;
  vst(xcol, xcol, yfancy, nrows) ;
  if (xmean != NULL) {
   xmean[col] = pmean*yfancy ; 
   xfancy[col] = yfancy ;
  }
  free(rawcol) ;
}
void
loadxdataind(double *xrow, SNP **snplist, int ind,  int ncols)             
{
 SNP *cupt ;
 Indiv *indx ;
 int i, j, k, n,  g ;

 for (i=0; i<ncols; i++) {  
   cupt = snplist[i] ;
   g = getgtypes(cupt, ind) ; 
   xrow[i] = (double) g ;
 }
}
void fixxrow(double *xrow, double *xmean, double *xfancy, int len) 
{
    int i ;

    vvt(xrow, xrow, xfancy, len) ;
    for (i=0; i<len; i++) { 
     if (xrow[i] < -0.1) xrow[i] = 0.0 ; 
     else xrow[i] -= xmean[i] ;
    }
}

void dofancy(double *cc, int n, double *fancy) 
{
  int i, t, nmiss=0 ;
  int top, bot ;  
  double p, yvar, y ;

  top = bot = 0 ;
  for (i=0; i<n; i++) {  
   t = nnint(cc[i]) ;
   if (t<0) { 
    ++nmiss ;
    continue ;
   }
   top += t ;
   bot += 2 ;
  }
  if (bot==0) return ;
  if (top == 0) return ;
  if (top == bot) return ;
  p = (double) top / ((double) bot) ;
  yvar = p*(1.0-p) ;
  y = 1.0/sqrt(yvar) ;
  vst(cc, cc, y, n )  ;                    
  *fancy = y ;
}

int vadjust(double *cc, int n, double *pmean, double *fancy) 
/* take off mean  force missing to zero */
/* set up fancy norming  */
{
 double p, ynum, ysum, y, ymean, yfancy = 1.0 ;
 int i, nmiss=0 ;

 ynum = ysum = 0.0 ;
 for (i=0; i<n; i++) {  
  y = cc[i] ;
  if (y < 0.0) { 
   ++nmiss ;
   continue ;
  }
  ++ynum ;
  ysum += y ;
 }
 if (ynum==0.0) { 
  fatalx("(vadjust) snp has no data\n") ;
 }
 ymean = ysum/ynum ;
 for (i=0; i<n; i++) {  
  y = cc[i] ;
  if (y < 0.0) cc[i] = 0.0 ;
  else cc[i] -= ymean ;
 }
 if (pmean != NULL) *pmean = ymean ;
  if (fancynorm) {  
   p = 0.5*ymean ;   // autosomes
   if (altnormstyle == NO) p = (ysum+1.0)/(2.0*ynum+2.0) ;
   y  = p * (1.0-p) ;  
   if (y>0.0) yfancy = 1.0/sqrt(y) ;
  }
 if (fancy != NULL) *fancy = yfancy ;
 return nmiss ;
}

double
dottest(char *sss, double *vec, char **eglist, int numeg, int *xtypes, int len) 
// vec will always have mean 0 
// perhaps should rewrite to put xa1 etc in arrays
{
   double *w1 ; 
   int *xt ;
   int i, k1, k2, k, j, n, x1, x2 ;
   double y1, y2, ylike, yl0, yl1, yl2 ;
   double ychi ;
   double *wmean ;
   int imax, imin ;
   static int ncall = 0 ;

   char ss1[MAXSTR] ;
   char ss2[MAXSTR] ;
   char sshit[4] ;
   double tail, ans, ftail, ftailx, ansx ; 

   ZALLOC(wmean, numeg, double) ;
   ZALLOC(w1, len, double) ;
   ZALLOC(xt, len, int) ;
   strcpy(ss1, "") ;

   calcmean(wmean, vec, len, xtypes, numeg) ;
   vlmaxmin(wmean, numeg, &imax, &imin) ;  
    if (chisqmode) {
     ylike = anova1(vec, len, xtypes, numeg) ;
     ans = 2.0*ylike ;
    }
    else {
     ans = ftail = anova(vec, len, xtypes, numeg) ;
    }
    ++ncall ;

    
    if (numeg>2) {
     sprintf(ss2, "%s %s ", sss, "overall") ;
     publishit(ss2, numeg-1, ans) ;
     printf(" %20s minv: %9.3f %20s maxv: %9.3f\n", 
     eglist[imin], wmean[imin], eglist[imax], wmean[imax]) ;
    }


    for (k1 = 0; k1<numeg; ++k1) { 
     for (k2 = k1+1; k2<numeg; ++k2) { 
      n = 0 ;
      x1 = x2 = 0 ; 
      for (i=0; i<len ; i++)   {  
        k = xtypes[i] ;
        if (k == k1) {  
         w1[n] = vec[i] ; 
         xt[n] = 0 ; 
         ++n ;
         ++x1 ;
        }
        if (k == k2) {  
         w1[n] = vec[i] ; 
         xt[n] = 1 ; 
         ++n ;
         ++x2 ;
        }
      }

     if (x1 <= 1) continue ;
     if (x2 <= 1) continue ;

     ylike = anova1(w1, n, xt, 2) ;
     ychi  = 2.0*ylike ;
     chitot[k1*numeg + k2]  += ychi ;
     if (chisqmode) {
      ansx = ychi ;
     }
     else {
      ansx = ftailx = anova(w1, n, xt, 2) ;
     }

      sprintf(ss2,"%s %s %s ", sss, eglist[k1], eglist[k2]) ;
      publishit(ss2, 1, ansx) ;
      
     }
    }
    free(w1) ;
    free(xt) ;
    free(wmean) ;
    return ans ;
}
double yll(double x1, double x2, double xlen) 
{
  double m1, m2, var ;

  if (xlen < 1.5 ) return 0.0 ;
  m1 = x1/xlen ;
  m2 = x2/xlen ;  
  var = m2-m1*m1 ;
  if (var <= 0.0) fatalx("bad yll\n") ;
  return -0.5*xlen*log(var) ;
}
void
calcmean(double *wmean, double *vec, int len, int *xtypes, int numeg) 
{

   int i, k ;
   double y1 ;
   double *w0, *popsize ;

   ZALLOC(w0, len, double) ;
   ZALLOC(popsize, numeg, double) ;

   y1 = asum(vec, len)/ (double) len ;  // mean
   vsp(w0, vec, -y1, len) ;

    for (i=0; i<len; i++)  { 
     k = xtypes[i] ;
     ++popsize[k] ;
     wmean[k] += w0[i] ;
    }



    vsp(popsize, popsize, 1.0e-12, numeg) ;
    vvd(wmean, wmean, popsize, numeg) ;

    free(w0) ;
    free(popsize) ;



}

double anova(double *vec, int len, int *xtypes, int numeg)
// anova 1 but f statistic
{
   int i, k ; 
   double y1, y2, ylike, top, bot, ftail ;  
   double *w0, *w1, *popsize, *wmean ;

   static int ncall2  = 0 ;

   if (numeg >= len) fatalx("bad anova\n") ;
   ZALLOC(w0, len, double) ;
   ZALLOC(w1, len, double) ;
   ZALLOC(wmean, numeg, double) ;
   ZALLOC(popsize, numeg, double) ;

   y1 = asum(vec, len)/ (double) len ;  // mean
   vsp(w0, vec, -y1, len) ;

    for (i=0; i<len; i++)  { 
     k = xtypes[i] ;
     ++popsize[k] ;
     wmean[k] += w0[i] ;
    }

/* debug */
    if (numeg == 2)  {  
     ++ncall2 ;
     for (i=0; i<len; ++i) {  
      if (ncall2<0) break ;
      k = xtypes[i] ;
//    printf("yy %4d %4d %12.6f %12.6f\n", i, k, vec[i], w0[i]) ;
     }
    }

    vsp(popsize, popsize, 1.0e-12, numeg) ;
    vvd(wmean, wmean, popsize, numeg) ;

    vvt(w1, wmean, wmean, numeg) ;
    top = vdot(w1, popsize, numeg) ;
    
    for (i=0; i<len ; i++)   {  
     k = xtypes[i] ;
     w1[i] = w0[i] - wmean[k] ;
    }
    bot = asum2(w1, len) / (double) (len-numeg) ;
    bot *= (double) (numeg-1) ;
    ftail = rtlf(numeg-1, len-numeg, top/bot) ;

    free(w0) ; 
    free(w1) ; 
    free(popsize) ;
    free(wmean) ;

    return ftail ;

}
double anova1(double *vec, int len, int *xtypes, int numeg)
{
   int i, k ; 
   double y1, y2, ylike ;  
   double *w0, *w1, *popsize, *wmean ;

   ZALLOC(w0, len, double) ;
   ZALLOC(w1, len, double) ;
   ZALLOC(wmean, numeg, double) ;
   ZALLOC(popsize, numeg, double) ;

   y1 = asum(vec, len)/ (double) len ;  // mean
   vsp(w0, vec, -y1, len) ;

    for (i=0; i<len; i++)  { 
     k = xtypes[i] ;
     ++popsize[k] ;
     wmean[k] += w0[i] ;
    }

    vsp(popsize, popsize, 1.0e-12, numeg) ;
    vvd(wmean, wmean, popsize, numeg) ;
    
    for (i=0; i<len ; i++)   {  
     k = xtypes[i] ;
     w1[i] = w0[i] - wmean[k] ;
    }

    y1 = asum2(w0, len) / (double) len ;
    y2 = asum2(w1, len) / (double) len ;
    ylike = 0.5*((double) len)*log(y1/y2) ;

    free(w0) ; 
    free(w1) ; 
    free(popsize) ;
    free(wmean) ;

    return ylike ;

}
void publishit(char *sss, int df, double chi) 
{
      double tail ;
      char sshit[4] ;
      char ss2[MAXSTR] ;
      int i, n ;
      char cblank, cunder ;
      static int ncall = 0 ;

      ++ncall ;
      cblank = ' ' ;
      cunder = '_' ;
      n = strlen(sss) ;

      strcpy(ss2, sss) ;  
      for (i=0; i< n; ++i)  { 
       if (ss2[i] == cblank)  ss2[i] = cunder ;
      }

      if (chisqmode) {
       if (ncall==1) printf("## Anova statistics for population differences along each eigenvector:\n");
       if (ncall==1) printf("%40s %6s %9s %12s\n", "", "dof", "chisq", "p-value") ;
       printf("%40s %6d %9.3f",ss2, ss2, df, chi) ;
       tail = rtlchsq(df, chi) ;  
       printf(" %12.6g", tail) ;
      }
      else {  
       if (ncall==1) printf("## Anova statistics for population differences along each eigenvector:\n");
       if (ncall==1) printf("%40s %12s\n", "",  "p-value") ;
       printf("%40s ", ss2) ;
       tail = chi ;  
       printf(" %12.6g", tail) ;
      }
      strcpy(sshit, "") ;  
      if (tail < pvhit) strcpy(sshit, "***") ;
      if (tail < pvjack) strcpy(sshit, "+++") ;
      printf(" %s", sshit) ;
      printf("\n") ;
}

void
setmiss(SNP **snpm, int numsnps) 
{
   SNP *cupt ;
   int i, j, t, n, tot ;

   for (i=0; i<numsnps; i++)  {  
    cupt = snpm[i] ;
    n = cupt -> ngtypes ;
    if (n <= 0) continue ;
    tot = 0 ;
    for (j=0; j<n; j++) {  
     if (getgtypes (cupt, j) >= 0) { 
      t = 1 ;
     }
     else { 
      t = 0 ;
     }
     putgtypes(cupt, j, t) ;
     tot +=t ;
    }
    if (verbose)
     printf("Valids: %s %d\n", cupt -> ID, tot) ;
   }
}

void
setfvecs(double *fvecs, double *evecs, int nrows, int numeigs) 
{

   double *w ;
   double xmax, xmin ;
   int i, j ;

   ZALLOC(w, nrows, double) ;

// NJP: clean this up sometime
   for (j=0; j<numeigs; j++) {  
    copyarr(evecs+j*nrows, w, nrows) ;
    copyarr(w, fvecs+j*nrows, nrows) ;
   }
   free(w) ;
}
void countpops(int ***counts, SNP **xsnplist, int *xindex, int *xtypes, int nrows, int ncols) 
// countpops is int [ncols][npops][2]  already zero filled
{
     int col, i, g1, g2, g, k1 ; 
     SNP *cupt ;
     int *rawcol ;

     ZALLOC(rawcol, nrows, int) ;
     for (col = 0; col < ncols; ++col) {  
      cupt = xsnplist[col] ;
      getrawcol(rawcol, cupt, xindex, nrows) ;
      for (i=0; i<nrows; i++) { 
       g = rawcol[i] ;
       if (g<0) continue ;
       k1 = xtypes[i] ;
       if (k1<0) continue ;
       ++counts[col][k1][g] ;
       ++counts[col][k1][2-g] ;
      }
     }
     free(rawcol) ;
}
double fstx(int ***counts, int ncols, int ntypes, int popa, int popb) 
// NP style fst estimator unsymmetric in A, B 
{
   double top, bot, xtop, xbot, pa, pb, ha, hb  ;
   int col, a0, a1, b0, b1, aa, bb ;

    
   xtop = xbot = 0.0 ;
   for (col=0; col<ncols; col++) {  

    a0 = counts[col][popa][0] ;
    a1 = counts[col][popa][1] ;

    b0 = counts[col][popb][0] ;
    b1 = counts[col][popb][1] ;

    aa = a0 + a1 ;
    bb = b0 + b1 ;

    if (aa==0) continue ;  
    if (bb==0) continue ;  
    if (a0==0) continue ;  
    if (a1==0) continue ;  // polymorphic in a 
      
    --a0 ; 
    --a1 ;  
    aa -= 2 ;
    if (aa==0) continue ;
/* ascertained in A */
    pa = (double) a0 / (double) aa ;
    pb = (double) b0 / (double) bb ;
    ha = (double) (a0*a1) / (double) (aa*(aa-1)) ;
    hb = (double) (b0*b1) / (double) (bb*(bb-1)) ;
    top = (pa-pb)*(pa-pb) - (ha/(double) aa) - (hb/(double) bb) ;
    bot = top + ha + hb ;  
    xtop += top ;
    xbot += bot ;
   }
   if (xbot < 1.0e-6) {  
    fprintf(stderr, "(fst) insufficient data for pops %d %d\n", popa, popb) ;
    return 0.0 ;
   }
   return xtop/xbot ;
}
void
dotpops(double *X, char **eglist, int numeg, int *xtypes, int nrows) 
{
      double *pp, *npp, val, yy ;
      int *popsize ;
      int i, j, k1, k2 ;


     ZALLOC(pp, numeg * numeg, double) ;
     ZALLOC(npp, numeg * numeg, double) ;
     popsize = xpopsize; 

     ivzero(popsize, numeg) ;

     for (i=0; i<nrows; i++) { 
      k1 = xtypes[i] ;
      ++popsize[k1] ;
      for (j=i+1; j<nrows; j++) { 
       k2 = xtypes[j] ;
       val = X[i*nrows+i] + X[j*nrows+j] - 2.0*X[i*nrows+j] ;
       pp[k1*numeg+k2] += val ;
       pp[k2*numeg+k1] += val ;
       ++npp[k1*numeg+k2]  ;
       ++npp[k2*numeg+k1]  ;
      }
     }
     vsp(npp, npp, 1.0e-8, numeg*numeg) ;
     vvd(pp, pp, npp, numeg*numeg) ;
// and normalize so that mean on diagonal is 1 
     yy = trace(pp, numeg) / (double) numeg ;
     vst(pp, pp, 1.0/yy, numeg*numeg) ;
     printf("\n## Average divergence between populations:");
     if (numeg<=10) {
      printf("\n") ;
      printf("%10s", "") ;
      for (k1=0; k1<numeg; ++k1) {  
       printf(" %10s", eglist[k1]) ;
      }
      printf("  %10s", "popsize") ;
      printf("\n") ;
      for (k2=0; k2<numeg; ++k2) {  
       printf("%10s", eglist[k2]) ;
       for (k1=0; k1<numeg; ++k1) {  
        val = pp[k1*numeg+k2] ;
        printf(" %10.3f", val) ;
       }
       printf("  %10d", popsize[k2]) ;
       printf("\n") ;
      }
     }
     else {   // numeg >= 10 
      printf("\n") ;
      for (k2=0; k2<numeg; ++k2) {  
       for (k1=k2; k1<numeg; ++k1) {  
        printf("dotp: %10s", eglist[k2]) ;
        printf(" %10s", eglist[k1]) ;
        val = pp[k1*numeg+k2] ;
        printf(" %10.3f", val) ;
        printf("    %10d", popsize[k2]) ;
        printf(" %10d", popsize[k1]) ;
        printf("\n") ;
       }
      }
     }
    printf("\n") ;
    printf("\n") ;
     

    free(pp) ;
    free(npp) ;

}
void printxcorr(double *X, int nrows, Indiv **indxx) 
{
   int k1, k2, t ; 
   double y1, y2, yy, rho ;
   Indiv *ind1, *ind2 ;

   if (pcorrmode == NO) return ;
   for (k1=0; k1<nrows; ++k1) {  
    for (k2=k1+1; k2<nrows; ++k2) {  

     ind1 = indxx[k1] ;
     ind2 = indxx[k2] ;

     t = strcmp(ind1 -> egroup, ind2 -> egroup) ;
     if (pcpopsonly && (t != 0))  continue ;
     

     y1 = X[k1*nrows+k1] ;
     y2 = X[k2*nrows+k2] ;
     yy = X[k1*nrows+k2] ;

     rho = yy/sqrt(y1*y2+1.0e-20) ;
     printf("corr: %20s %20s %20s %20s %9.3f\n", 
      ind1 -> ID, ind2 -> ID, ind1 -> egroup, ind2 -> egroup, rho) ;

    }
   }
}
void clearld(double *ldmat, double *ldvv, int rsize, int n, int nclear)    
{
  int i, j, lo, hi ;

  if (nclear>rsize) fatalx("bad nclear\n") ;

  lo = rsize-nclear ;  
  hi = rsize-1 ;

  for (i=lo; i<=hi; ++i)  { 
   vzero(ldvv+i*n, n) ;
    for (j=0; j<=hi; ++j)  { 
      ldmat[j*rsize+i] = ldmat[i*rsize+j] = 0.0 ;
    }
// force matrix non-singular 
    ldmat[i*rsize+i] = 1.0e-8 ;  
  }
}

void
ldreg(double *ldmat, double *ldmat2, double *vv, double *vv2, double *ldvv, 
  double *ldvv2, int rsize, int n) 
/** ldmat2 is inner product matrix for last rsize columns on exit */
{
   int i, j, k1, k2 ; 
   double *rr, *ans, *tt ;
   double y ;

   ZALLOC(rr, rsize, double) ;
   ZALLOC(ans, rsize, double) ;
   ZALLOC(tt, n, double) ;

   if (rsize>1) 
    copyarr(ldvv, ldvv2+n, n*(rsize-1)) ;
   for (i=0; i<rsize-1 ; i++)  { 
    for (j=0; j<rsize-1 ; j++)  { 
     k1 = i*rsize+j  ;
     k2 = (i+1)*rsize+j+1 ;
     ldmat2[k2] = ldmat[k1] ;
    }
   }
   copyarr(vv, ldvv2, n) ;
   i = 0 ;
   for (j=0; j<rsize ; j++)  { 
     y = rr[j] = vdot(vv, ldvv+j*n, n) ;
     y =  vdot(vv, ldvv2+j*n, n) ;
     if (j==0) y += 1.0e-6 ;
     ldmat2[i*rsize+j] = ldmat2[j*rsize+i] = y ;
   }
  solvit(ldmat, rr, rsize, ans) ; /* solve normal equations */
  copyarr(vv, vv2, n) ;
  for (i=0; i<rsize; i++) {  
   vst(tt, ldvv+i*n, -ans[i], n) ;
   vvp(vv2, vv2, tt, n) ;
  }
  free(rr) ;
  free(ans) ;
  free(tt) ;
}
void
setidmat(double *mat, int n, double scal) 
{
  int i ;
  
  vzero(mat, n*n) ; 
  for (i=0; i<n; i++) { 
   mat[i*n+i] = scal ;
  }
}
void fixrho (double *a, int n) 
// turn a into correlation matrix
{
   double *d, *tt, y ;
   int i ;  

   ZALLOC(d, n, double) ;
   ZALLOC(tt, n*n, double) ;

   getdiag(d, a, n) ;

   vsqrt(d, d, n) ;
   addouter(tt, d, n) ;

   vvd(a, a, tt, n*n) ;


  free(d) ;
  free(tt) ;


}
void printdiag(double *a, int n) 
{
   double *d, *tt, y ;
   int i ;  

   ZALLOC(d, n, double) ;
   getdiag(d, a, n) ;
    y = asum(d,n) / (double) n ; 
    for (i=0; i<n; i++) {  
     printf("diag: %9.3f\n", d[i]/y) ;
    }


   free(d) ;

}
int
ridoutlier(double *evecs, int n, int neigs, 
 double thresh, int *badlist) 
{
/* badlist contains list of outliers */
 double *ww, y1 , y2 ; 
 int *vbad ;
 int i, j  ;
 int nbad = 0 ; 
 OUTLINFO *outpt ;

 ZALLOC(ww, n, double) ;
 ZALLOC(vbad, n, int) ;
 for (j=0; j<n; j++)  {  
  outpt = outinfo[j] ;
  outpt -> vecno = -1 ;
 }
 for (i=0; i<neigs; ++i) {  
  copyarr(evecs+i*n, ww, n) ;
  y1 = asum(ww, n) / (double) n ;
  vsp(ww, ww, -y1, n) ;
  y2 = asum2(ww, n) / (double)  n ;
  y2 = sqrt(y2) ;
  vst(ww, ww, 1.0/y2, n) ;

  for (j=0; j<n; j++) {  
   if (fabs(ww[j])>thresh) { 
    vbad[j] = 1 ;
    outpt = outinfo[j] ;
    if (outpt -> vecno < 0) {  
     outpt -> vecno = i ;
     outpt -> score = ww[j] ; 
    }
   }
  }
 }
 for (j=0; j<n; j++) {  
  if (vbad[j] == 1) { 
   badlist[nbad] = j ;
   ++nbad ;
  }
 }
 free(ww) ; 
 free(vbad) ;
 return nbad ;

}
/** this is the code to parallelize */
void
domult(double  *tvecs, double  *tblock, int numrow, int len) 
{
  int i ;
  double ycheck ;
  vzero(tvecs, len*len) ;
  for (i=0; i<numrow; i++) {  
    ycheck = asum(tblock+i*len, len) ;  
    if (fabs(ycheck)>.00001) fatalx("bad ycheck\n") ;
    addoutersym(tvecs, tblock+i*len, len) ;
  }
}

void
addoutersym(double *X, double *v, int n)  
{
  int i, j ;

  for (i=0; i<n; i++) { 
    for (j=i; j<n; j++) {
      X[i*n+j] += v[i]*v[j] ;
    }
  }
}

void symit(double *X, int n) 
{
  int i, j ;

  for (i=0; i<n; i++) { 
    for (j=i+1; j<n; j++) {
      X[j*n+i] = X[i*n+j]  ;
    }
  }
}

double dofst(double *fstans, double *fstsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg) 

{

   int t1, t2 ;
   double y, sd ; 

   vclear(fstans, 0.0, numeg*numeg) ;
   vclear(fstsd, 1.0, numeg*numeg) ;

   for (t1=0; t1<numeg; ++t1) {  

     fstans[t1*numeg+t1] = 0.0 ;

      for (t2=t1+1; t2<numeg; ++t2) {  

       if (xpopsize[t1] <= 0) continue ;
       if (xpopsize[t2] <= 0) continue ;

       y =  fst(xsnplist, xindex, xtypes, nrows, ncols, t1, t2, &sd) ;

       fstans[t1*numeg+t2] = fstans[t2*numeg+t1] =  y ;
       fstsd[t1*numeg+t2] = fstsd[t2*numeg+t1] =  sd ;

      }
   }


}

double fst(SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int type1, int type2, double *psd) 
{
  double en, ed, tn, td ; 
  int col ;
  SNP *cupt ;
  double *xn, *xd, xdsum = 0.0 ;  
  double xmean, xsd ;
  int num=0 ;

  en = ed = 0.0 ;

  ZALLOC(xn, ncols, double) ;
  ZALLOC(xd, ncols, double) ;

  for (col=0; col<ncols; ++col) {  
   cupt = xsnplist[col] ;
   fstcol(&tn, &td, cupt, xindex, xtypes, nrows, type1, type2) ;  
   if (td<0.0) continue ;
   xn[num] = tn ; 
   xd[num] = td ;
   ++num ;
   xdsum += td ;
   if (td<0.0) fatalx("bad td\n") ;
   if (tn>td) fatalx("bad tn\n") ;
  }

  if ((num<2) || (xdsum <=0.0))  {
   xmean = 0.0 ; 
   xsd = 100.0 ; 
  }

  else {
   jackrat(&xmean, &xsd, xn, xd,  num) ;
  }

  *psd = xsd ;

  free(xn) ; 
  free(xd) ; 

  return (xmean) ;
}

void jackrat(double *xmean, double *xsd, double *top, double *bot, int len) 
{
  double yn, yt, yb, xt, xb,  yh, y, ylen ; 
  double *hh, *ttau, *tjack, tmean, ttmean, ymean  ; 
  int i ;


  yt = asum(top, len) ;
  yb = asum(bot, len) ; 
  
  if (yb<=0.0) {  
   fatalx("(jackrat) no data\n") ;
  }

  if (len<=1) fatalx("(jackrat) length must be 2 or more\n") ;

  ZALLOC(tjack, len, double) ;  
  ZALLOC(ttau, len, double) ;  

  tmean = yt/yb ;  
  for (i=0; i<len; i++) { 
   xt = yt-top[i] ;
   xb = yb-bot[i] ;
   if (xb<=0.0) {
    tjack[i] = tmean ; 
    continue ;
   }
   else  {
    ttmean = xt/xb ;
    tjack[i] = ttmean ;
   }
  }
  ymean = asum(tjack, len) / (double) len ;
/* ymean is now jackknifed estimator */
  ylen = (double) len ;
  for (i=0; i<len; i++) {  
   ttau[i] = ylen*ymean-(ylen-1.0)*tjack[i] ;
/* pseudovalue */
  }
  vsp(ttau, ttau, -ymean, len) ; 
  xt = asum2(ttau, len) ; 
  y = xt  / (ylen*(ylen-1.0)) ; /* var */
  if (y<=0.0) y = 0.0 ; 
  else y = sqrt(y) ; 
  
  *xmean = ymean ; 
  *xsd = y ; 

  free(tjack) ;
  free(ttau) ;

}    


double fstcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) 
{
   int c1[2], c2[2], *cc ;
   int *rawcol ;
   int k, g, i ; 
   double ya, yb, yaa, ybb, p1, p2, en, ed ;
   double z, zz, h1, h2, yt ;


   ZALLOC(rawcol, nrows, int) ;

   getrawcol(rawcol, cupt, xindex, nrows)  ;

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;
    cc = NULL ;
    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (cc == NULL) continue ;
    g = rawcol[i] ;
    if (g<0) continue ;  
    cc[0] += g ; 
    cc[1] += 2-g ;
   }
/**
   printf("zz\n") ;
   printimat(c1, 1, 2) ;
   printimat(c2, 1, 2) ;
*/

   ya = c1[0] ;
   yb = c1[1] ;
   yaa = c2[0] ;
   ybb = c2[1] ;
   z = ya + yb ;
   zz = yaa+ybb ;
   if ((z<0.1) || (zz<0.1)) { 
    *estn = 0.0 ;
    *estd = -1.0 ;
    free(rawcol) ;
    return 0.0;
   }

   yt = ya+yb ;
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;

   yt = yaa+ybb ;
   p2 = yaa/yt ;         
   h2 = yaa*ybb/(yt*(yt-1.0)) ;

   en = (p1-p2)*(p1-p2) ;  
   en -= h1/z ; 
   en -= h2/zz ; 
 
   ed = en ; 
   ed += h1 ; 
   ed += h2 ;

   *estn = en ; 
   *estd = ed ; 
   

   free(rawcol) ;
   return z + zz ;

}

void
 setplimit(Indiv **indivmarkers, int numindivs, 
 char **eglist, int numeg, int plimit)  
{
     int *indnums ;
     int *psize ;
     int i, k, kk ;
     Indiv *indx ; 

     ZALLOC(indnums, numindivs, int) ;
     ZALLOC(psize, numeg, int) ;


     idperm(indnums, numindivs) ;
     ranperm(indnums, numindivs) ;

     for (i=0; i<numindivs; i++)  {  
      k = indnums[i] ; 
      indx = indivmarkers[k] ;
      if (indx -> ignore) continue ;
      kk = indxindex(eglist, numeg, indx -> egroup) ;
      if (kk<0) continue ; 
      ++psize[kk] ; 
      if (psize[kk] > plimit) indx -> ignore = YES ;
     }



     free(psize) ;
     free(indnums) ;

}
void
writesnpeigs(char *snpeigname, SNP **xsnplist, double *ffvecs, int numeigs, int ncols) 
{
// this is called at end and ffvecs overwritten
  double *xpt, y, yscal, *snpsc ;
  int i, j, k, kmax, kmin ;
  SNP * cupt  ;
  FILE *fff ;
  
  for (j=0; j<numeigs; ++j) {  
   xpt = ffvecs+j*ncols ;  
   y = asum2(xpt, ncols) ;  
   yscal = (double) ncols / y ;
   yscal = sqrt(yscal) ;
   vst(xpt, xpt, yscal, ncols) ;
  }


  ZALLOC(snpsc, ncols, double) ;
  for (j=0; j<numeigs; ++j) {  
   for (i=0; i<ncols; ++i) {  
    y = ffvecs[j*ncols+i] ;
    snpsc[i] = fabs(y) ; 
   }
   for (k=0; k<=10; ++k) { 
    vlmaxmin(snpsc, ncols, &kmax, &kmin) ;
    cupt = xsnplist[kmax] ;
    printf("eigbestsnp %4d %20s %2d %12d %9.3f\n", j+1, cupt -> ID, cupt -> chrom, nnint(cupt -> physpos), snpsc[kmax]) ;
    snpsc[kmax] = -1.0 ;
   }
  }
  free(snpsc) ;


  if (snpeigname == NULL) return ;
  openit (snpeigname, &fff, "w") ;

  for (i=0; i<ncols; ++i) {  
   cupt = xsnplist[i] ;
   if (cupt -> ignore) continue ;

   fprintf(fff, "%20s", cupt -> ID) ;
   fprintf(fff,  " %2d", cupt -> chrom) ;
   fprintf(fff,  " %12d", nnint(cupt -> physpos)) ;

   for (j=0; j<numeigs; ++j) {  
    fprintf(fff, " %9.3f", ffvecs[j*ncols+i]) ;  
   }
   fprintf(fff, "\n") ;
  }

  fclose(fff) ;

}
