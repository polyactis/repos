 
/////MCMC for Mixture of Polya Trees -Meijuan Li July 2008 --C++ used for the MTP paper///////
/* standard libraries to include */
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

#pragma hdrstop
#define FREE_ARG char*

#include "ranker.h"


static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b)>= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
   (dmaxarg1) : (dmaxarg2))
static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
   (dminarg1) : (dminarg2))
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
   (maxarg1) : (maxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
   (iminarg1) : (iminarg2))

double minFcn(double *bb, double sigma, double tauzu, double *zu,double *prob);

void dfunc(double *p,double *df);
double gasdev(long *idum);
double ran1(long *idum);
int **imatrix(long nrl,long nrh,long ncl,long nch);
void free_dmatrix(double **m, long nrl,long nrh,long ncl,long nch);
void free_matrix(float **m, long nrl,long nrh,long ncl,long nch);
int *ivector(long nl,long nh);
void free_dvector(double *v,long nl,long nh);
void free_ivector(int *v,long nl, long nh);
void choldc(double **a, int n, double p[]);
void nrerror(char error_text[]);
double mean(double *data, int n);
double var(double *data, int n);

void print1d(double *oned, int sizee);
void print2d(double **twod, int row, int col);
void print2di(int **twod, int row, int col);

double *getVector(double **a, int row_b, int row_e, int col_b, int col_e);
double *getVector_row(double **a, int row_b, int row_e, int col_b, int col_e);
double **getMatrix(double **a, int row_b, int row_e, int col_b, int col_e);
void sort_vector(double *a,int size);

double **dmatrix(long nrl,long nrh,long ncl,long nch);
void copy_dmatrix(double **orig_mat, int rb, int re, int cb, int ce, double **copy_mat );
double gamdev(double ia, long *idum);

double normal_pdf ( double x, double a, double b );
double normal_dist(double x);
double betadev(double a,double b, long *idum);
double gammln(double xx);
double ran1(long *idum);
double *dvector(long nl,long nh);
void get_initValu(long setSeed,int numpar, double *initValu_Mu, double **initValu_V,double *initValu);

 //////////////////////////////////////////////////////////// glob variable
int noPar=11;
int nGlob=100;
int M=4;
double mfactor=20, cfactor=1;
double a=0.01, b=0.01;
double *XX1, *XX2, *XX3, *XX4, *XX5, *XX6,*XX7, *XX8, *XX9, *y1Glob, **kin;
int M2=(int)pow(2,M);
int M2_1=(int)pow(2,(M-1));

int main()
{

  double **Xdat, **initValall, **covMat,**snpAdat, **bigyGlob;
  int N_gene=500;
  int nr, i,  j, jj, gibbs_N;
  Xdat=dmatrix(1,nGlob,1,8);
  initValall=dmatrix(1,N_gene,1,noPar);
  covMat=dmatrix(1,(N_gene*noPar),1,noPar);
  snpAdat=dmatrix(1,nGlob,1,N_gene);
  kin=dmatrix(1,nGlob,1,nGlob);
  bigyGlob=dmatrix(1,nGlob,1,N_gene);

  //This holds the response data for all genes
  ifstream filer0("Ydat.txt");
   if(!filer0){
      cerr << "Failure to open filer0\n";
      exit(EXIT_FAILURE);
   }

 filer0.seekg(0);
  for(i=1;i<=nGlob;i++){
	for(j=1;j<=N_gene;j++){
     filer0 >> bigyGlob[i][j];
    }
  }

  //This holds the candatae gene genotypes
 ifstream filer1("snpdat.txt");
   if(!filer1){
      cerr << "Failure to open filer1\n";
      exit(EXIT_FAILURE);
   }   
   
  filer1.seekg(0);
   for(i=1;i<=nGlob;i++){
      for(j=1;j<=N_gene;j++){
         filer1 >> snpAdat[i][j];
      }
   }


   //This holds the population structure
   ifstream filer2("PCA.txt");
   if(!filer2){
      cerr << "Failure to open filer2\n";
      exit(EXIT_FAILURE);
   }

  filer2.seekg(0);
   for(i=1;i<=nGlob;i++){
      for(j=1;j<=8;j++){
         filer2 >> Xdat[i][j];
      }
   }
  
  XX1=dvector(1,nGlob);
  XX2=dvector(1,nGlob);
  XX3=dvector(1,nGlob);
  XX4=dvector(1,nGlob);
  XX5=dvector(1,nGlob);
  XX6=dvector(1,nGlob);
  XX7=dvector(1,nGlob);
  XX8=dvector(1,nGlob);
  XX1=getVector(Xdat, 1, nGlob, 1,1);
  XX2=getVector(Xdat, 1, nGlob, 2,2);
  XX3=getVector(Xdat, 1, nGlob, 3,3);
  XX4=getVector(Xdat, 1, nGlob, 4,4);
  XX5=getVector(Xdat, 1, nGlob, 5,5);
  XX6=getVector(Xdat, 1, nGlob, 6,6);
  XX7=getVector(Xdat, 1, nGlob, 7,7);
  XX8=getVector(Xdat, 1, nGlob, 8,8);

  // cout << "X8 "<<nGlob<<endl;
  //print1d(XX8, nGlob);

  //This holds initial values for the all genes*/
  ifstream filer3("PCApar.txt");
   if(!filer3){
      cerr << "Failure to open filer3\n";
      exit(EXIT_FAILURE);
   }   

 filer3.seekg(0);
   for(i=1;i<=N_gene;i++){
      for(j=1;j<=noPar;j++){
         filer3 >> initValall[i][j];
      }
   }
 
   //This holds covaraince matrix for parameters/
  ifstream filer4("PCAparcov.txt");
   if(!filer4){
      cerr << "Failure to open filer4\n";
      exit(EXIT_FAILURE);
   }  


  filer4.seekg(0);
   for(i=1;i<=(N_gene*noPar);i++){
      for(j=1;j<=noPar;j++){
         filer4 >> covMat[i][j];
      }
   }

  
   //This hold P_K^c
   ifstream filer5("ZPCAdat.txt");
   if(!filer5){
      cerr << "Failure to open filer5\n";
      exit(EXIT_FAILURE);
   }

  filer5.seekg(0);
   for(i=1;i<=nGlob;i++){
      for(j=1;j<=nGlob;j++){
         filer5 >> kin[i][j];
      }
   }
  
  //This is file for quantile and standard deviation of candidate gene regression coefficients
  ofstream filewgene("gnbetaPCA.txt");
  if(!filewgene){
     cerr << "Failure to open filewgene\n";
     exit(EXIT_FAILURE);
  }	 
 
  //This is a file for accept rate
  ofstream filewrate("acceptRPCA.txt");
  if(!filewrate){
     cerr << "Failure to open filewrate\n";
     exit(EXIT_FAILURE);
  }	


  //////////////////////////// initialize some parameter values within LOOP/////////////////////////
  long setSeed=-8623;
  int   GI=5000, thinfac1=5, thinfac2=2, tree=2,chain, chain_N=2, GG=(GI*chain_N);
  int   burnin=6000;
  ////////////////////////////////////////// start gene loop//////////////////////
  for (nr=1;nr<=N_gene;nr++){
  double **covb=dmatrix(1,noPar,1,noPar);
  XX9=dvector(1,nGlob);
  y1Glob=dvector(1,nGlob);
  int **tally;
  int	  tallyb, thin, tallyu, r, k, j;
  double *initVal=dvector(1, noPar); 
  double *initValu=dvector(1, noPar); 
  double  so, sn,llnew, llold;
  double  llposterior;
  double  **zo=dmatrix(1,M,1,M2);
  double  **zd=dmatrix(1,M,1,M2);
  double *po=dvector(1,M2);
  double *pn=dvector(1,M2);
  double *betait1=dvector(1,GG);
  double *bsigmao=dvector(1,noPar);
  double *bsigman=dvector(1,noPar);
  double *bo=dvector(1,(noPar-1));
  double *bn=dvector(1,(noPar-1));
  double *zuo=dvector(1,nGlob);
  double *zun=dvector(1,nGlob);
  double *pp=dvector(1,noPar);
  double ranbeta, sum, mu,tauzu;
  double *temp=dvector(1, nGlob);
  double *mm=dvector(1,noPar);	
  double scalingF=0.02;
  double scalingR=0.002;

 tally=imatrix(1,M,1,M2_1);
  for( i=1;i<=M;i++){
     for(j=1;j<=M2_1;j++){
        tally[i][j]= 0;
      }
  }	 
 tallyb=0;
 tallyu=0;

	XX9=getVector(snpAdat, 1, nGlob, nr,nr);
    cout << "XX9 ="<< nGlob << endl;
    print1d(XX9, nGlob);


	y1Glob=getVector(bigyGlob, 1, nGlob, nr,nr);

    cout << "y1Glob ="<< nGlob << endl;
    print1d(y1Glob, nGlob);


    covb=getMatrix(covMat, (nr-1)*noPar+1, nr*noPar, 1,noPar);
	initVal=getVector_row(initValall, nr, nr, 1,noPar);
	cout<<"initVal= "<< noPar<<endl;
	print1d(initVal, noPar);

  /******** chold decomposition for covariance matrix ****/
   for(i=1;i<=noPar;i++){
    for(j=1;j<=noPar;j++){
      covb[i][j]=scalingF*covb[i][j];
    }
  }
  choldc(covb,noPar,pp);

  ////////////////////////////////// Start GI loop //////////////////////////////
  for(chain=1; chain<=chain_N; chain++){
	  double rho=0,tauzu=0.01;
	  for(i=1;i<=M;i++){
		for(j=1;j<=M2;j++){
			zo[i][j]= 0.5;
					}
				}
	 for (i=1; i<=nGlob;i++){
	   zuo[i]=0;
		}
	 for(i=1;i<=M2;i++) {
		 po[i]= (double)1/M2; }
	 get_initValu(setSeed,noPar, initVal, covb,initValu);
	 cout<<"initial value for chain "<<chain<<endl;
  	 print1d(initValu, noPar);

	 for(i=1;i<=(noPar-1);i++) { bo[i]=initValu[i];}
	  so=initValu[noPar];


	  llposterior=minFcn(bo, so, tauzu, zuo,po);
	  cout << "llposterior at initil value --\n";
	  cout << llposterior << endl;

     for(gibbs_N=1;gibbs_N<=(GI+burnin);gibbs_N++){
     for(thin=1;thin<=thinfac1;thin++){
		for (i=tree; i<=M; i++)
        {
            for (j=1; j<=(int)pow(2,(i-1)); j++)
            {
              /*copy zo matrix to zd matrx: zd=zo*/
              copy_dmatrix(zo, 1, M, 1, M2, zd);		 
			  ranbeta= betadev(mfactor*zo[i][2*j-1],mfactor*zo[i][2*j],&setSeed);
              zd[i][2*j-1]=ranbeta;
              zd[i][2*j]=1-ranbeta;

                for (k=2; k<=M; k++)
                {
                  for (int s=1; s<=(int)pow(2,k); s++)
                  {
                    zd[k][s]=zd[k-1][(int)(ceil((double)(1.000*s/2)))]*zd[k][s];
                  }
                }
                for ( k=1; k<=M2; k++)
                {
                  pn[k]=zd[M][k];
                }
          
		
			   llnew=minFcn(bo, so, tauzu, zuo,pn);
               llold=minFcn(bo, so, tauzu, zuo,po);
              
			rho=llnew-llold+(mfactor*ranbeta-cfactor*((int)pow(i,2)))*log(zo[i][2*j-1])
				+(mfactor*(1-ranbeta)-cfactor*((int)pow(i,2)))*log(zo[i][2*j])
				-(mfactor*zo[i][2*j-1]-cfactor*((int)pow(i,2)))*log(ranbeta)
				-(mfactor*zo[i][2*j]-cfactor*((int)pow(i,2)))*log(1-ranbeta)
				-gammln(mfactor*ranbeta)-gammln(mfactor*(1-ranbeta))
				+gammln(mfactor*zo[i][j*2-1])+gammln(mfactor*zo[i][j*2]);
				 
            //sample random_number(u)
               if (ran1(&setSeed)<exp(rho)){
                   for ( int k=1; k<=M2; k++) { po[k]=pn[k]; }          
                   zo[i][2*j-1]=ranbeta;
                   zo[i][2*j]=1-ranbeta;
                   tally[i][j]=tally[i][j]+1;
               }  //end if

             } // end of for j=1<=pow
         } // end of  for tree =M

      //sample beta, sigma, kinship effect

     for (i=1; i<=thinfac2; i++)
	  {
		///################################sample zu
	   for(r=1;r<=nGlob;r++) temp[r]=gasdev(&setSeed);
	   for(r=1;r<=nGlob;r++) {
		   zun[r]=zuo[r]+scalingR*temp[r]/sqrt(tauzu);
	   }		

	   /*	 mu=mean(zun, nGlob);
		 for(r=1;r<=nGlob;r++) {
		   zun[r] -=mu;
		 }*/	
	
		 llnew=minFcn(bo, so,tauzu,zun, po);
		 llold=minFcn(bo, so,tauzu,zuo,po);
			
		 if (ran1(&setSeed)<exp(llnew-llold)){
			 for(r=1;r<=nGlob;r++) {
					 zuo[r]=zun[r];
						}
					 tallyu +=1;
			}  //end if
	
 	  	 //#############################sample tauzu
		   sum=0;
		   for (i=1;i<=nGlob; i++)
				{
				sum+=0.5*DSQR(zuo[i]);
					}
		  tauzu=gamdev((a+0.5*nGlob), &setSeed);
		  tauzu=tauzu/(b+sum);

		//####################### sample beta and sigma

          for (j=1; j<=(noPar-1); j++) {bsigmao[j]= bo[j];} 
          bsigmao[noPar]= so;
          for(j=1;j<=noPar;j++) mm[j]=gasdev(&setSeed);
          for(j=1;j<=noPar;j++){
            sum=0.0;
            for(k=1;k<j;k++) sum+=covb[j][k]*mm[k];
				bsigman[j]=bsigmao[j]+sum+pp[j]*mm[j];
			}
   
	      if (bsigman[noPar]>0) {       
	      for (j=1; j<=(noPar-1); j++) {bn[j]=bsigman[j];}
		  sn=bsigman[noPar];  
   	
          llnew=minFcn(bn, sn, tauzu,zuo, po);
          llold=minFcn(bo, so, tauzu,zuo, po);

            //sample random_number(u)
            if (ran1(&setSeed)<exp(llnew-llold)){
                 so=sn;
   			   for (jj=1; jj<=(noPar-1); jj++) { bo[jj]=bn[jj];}  
               tallyb +=1;
          }  //endof if 

		  } // end if bsigman[noPar]>
      } //end of thinfac2 loop 

    } // end of thin loop

 
    ///////////////////////////////// save MCMC for sigma beta, and prob
  	if(gibbs_N>burnin){
	  betait1[((chain-1)*GI+(gibbs_N-burnin))]=bo[(noPar-1)];
      cout << " gene "<<nr<<"  "<<bo[noPar-1]<<endl;
		}  // end if gibbs_N>burin
   }//end gibbs_N

}  //end chain_N
  
 ////////////////////// end GI loop for(gibbs_N=1;gibbs_N<=GI;gibbs_N++) //MCMC loop ends
  
  cout<<"# accpet for zu= "<<tallyu<<endl;
  cout<<"accept rate for zu= "<<(double)(1.0*tallyu/(chain_N*(GI+burnin)*thinfac1*thinfac2))<<endl;

  cout<<"# accept for beta and sigma: "<<tallyb<<endl;
  cout<<" \n accept rate beta and sigma= "<<(double)(1.0000*tallyb/(chain_N*(GI+burnin)*thinfac1*thinfac2));
  cout<<endl;
  
  for (i=tree; i<=M; i++){
    for (j=1; j<=(int)pow(2,(i-1)); j++)
    {
      cout<<"\nlevel= "<<i<<" place= "<<j<<" accept rate= "<<(double)(1.0000*tally[i][j]/(chain_N*(GI+burnin)*thinfac1))<<" ";
    }
  }
  cout<<endl;

  ////////////////// get qualtile for beta of gene

 filewgene<<nr<<" "<<quantile(betait1, GG, 0.025)<<" "<<quantile(betait1, GG, 0.5)<<" "<<quantile(betait1, GG, 0.975)<<" "<<sqrt(var(betait1,GG))<<" "<<quantile(betait1, GG, 0.005)<<" "<<quantile(betait1, GG, 0.995)<<" ";
 filewgene<<endl; 

  ////////////////////////// 
  free_dvector(zuo,1,nGlob);
  free_dvector(zun,1,nGlob);
  free_dmatrix(zo,1,M,1,M2);
  free_dmatrix(zd,1,M,1,M2);
  free_dvector(po,1,M2);
  free_dvector(pn,1,M2);
  free_dvector(bsigmao,1,noPar);
  free_dvector(bsigman,1,noPar);
  free_dvector(bo,1,(noPar-1));
  free_dvector(bn,1,(noPar-1));
  free_dvector(pp,1,noPar);
  free_dvector(initVal,1, noPar);
  free_dmatrix(covb,1,noPar,1,noPar);
  free_dvector(betait1,1,GG);
  free_dvector(XX9,1,nGlob);
  free_dvector(y1Glob,1,nGlob);
  free_dvector(temp,1, nGlob);
  free_dvector(mm,1,noPar);
                       } /////////////////////////////////////// end of loop for N_gen
  free_dmatrix(kin,1,nGlob,1,nGlob);
  free_dvector(XX1,1,nGlob);
  free_dvector(XX2,1,nGlob);
  free_dvector(XX3,1,nGlob);
  free_dvector(XX4,1,nGlob);
  free_dvector(XX5,1,nGlob);
  free_dvector(XX6,1,nGlob);
  free_dvector(XX7,1,nGlob);
  free_dvector(XX8,1,nGlob);
  free_dmatrix(Xdat, 1,nGlob,1,8);
  free_dmatrix(initValall,1,N_gene,1,noPar);
  free_dmatrix(snpAdat,1,nGlob,1,N_gene);
  free_dmatrix(covMat,1,(N_gene*noPar),1,noPar);
  free_dmatrix(bigyGlob,1,nGlob,1,N_gene);

} ///////////////////////////////////////////////////////////////////////////////////////////end of main loop 
 

/////////////////////////////////////////////// minFunction //////////////////////////////////////////

double minFcn(double *bb, double sigma, double tauzu, double *zu,double *prob)
{
  int i,j,r,k;
  int Ni;
  double res=0.0;
  double mui, cd; 
  double PNi;

  /******** log(posterior), sigma is on original-scale, zu ~Gamma(a,b)***********/  
  for(i=1;i<=nGlob;i++)
  {   
	cd=0;
	for(j=1; j<=nGlob; j++){
		   cd+=kin[i][j]*zu[j];
			}
    mui=bb[1]+XX1[i]*bb[2]+XX2[i]*bb[3]+XX3[i]*bb[4]+XX4[i]*bb[5]+XX5[i]*bb[6]+XX6[i]*bb[7]+XX7[i]*bb[8]+XX8[i]*bb[9]+XX9[i]*bb[10]+cd;
	Ni=(int)(M2*normal_dist((y1Glob[i]-mui)/sqrt(sigma))+1);
	if(Ni>M2) Ni=M2;
	PNi=prob[Ni]; 
	res+=log(M2)+log(PNi)+log(normal_pdf(y1Glob[i],mui, sqrt(sigma)))+log(normal_pdf(zu[i],0, 1/sqrt(tauzu)));
  }
  res+=(a-1)*log(tauzu)-b*tauzu;
  res+=-log(sigma);
  return res;
  }


/* this function draws a uniform [0,1] deviate */
#define IAran1 16807
#define IMran1 2147483647
#define AMran1 (1.0/IMran1)
#define IQran1 127773
#define IRran1 2836
#define NDIVran1 (1+(IMran1-1)/32)
#define EPSran1 1.2e-7
#define RNMXran1 (1.0-EPSran1)
double ran1(long *idum)
{
   const int NTABran1=32;
   int j;
   long k;
   static long iy=0;
   static long iv[NTABran1];
   double temp;
   if (*idum <= 0 || !iy) {
      if (-(*idum)<1) *idum=1;
      else *idum = -(*idum);
      for(j=NTABran1+7;j>=0;j--) {
         k=(*idum)/IQran1;
         *idum=IAran1*(*idum-k*IQran1)-IRran1*k;
         if (*idum<0) *idum += IMran1;
         if (j<NTABran1) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQran1;
   *idum=IAran1*(*idum-k*IQran1)-IRran1*k;
   if (*idum<0) *idum += IMran1;
   j=iy/NDIVran1;
   iy=iv[j];
   iv[j] = *idum;
   if ((temp=AMran1*iy)>RNMXran1) return RNMXran1;
   else return temp;
}

double gasdev(long *idum)
{
   /* draws a standard normal random variable */
   double ran1(long *idum);
   static int iset=0;
   static double gset;
   double fac,rsq,v1,v2;
   if(iset==0){
      do{
         v1=2.0*ran1(idum)-1.0;
         v2=2.0*ran1(idum)-1.0;
         rsq=v1*v1+v2*v2;
      } while (rsq>=1.0 || rsq==0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   } else{
      iset=0;
      return gset;
      }
}

/*
The following code is for beta distribution random generator. The link is
http://www.biostat.umn.edu/~cavanr/geneNormRepHier.txt.
*/

double gamdev(double ia, long *idum)
{
   /* gets a gamma deviate */
   double am,e,s,v1,v2,x,y,p,q;
   if(ia<=0.0) {
     cout << "you input parameter is "<<ia;
     nrerror("error in gamdev: parameter non-positive");
   }
   if(ia==1.0) x=-log(ran1(idum));
   else{
      if(ia<1.0){
         p=exp(1.0)/(ia+exp(1.0));
         do{
            if(ran1(idum)<p){
               x=pow(ran1(idum),1.0/ia);
               q=exp(-x);
            }
            else{
               x=1.0-log(ran1(idum));
               q=pow(x,ia-1.0);
            }
         } while(ran1(idum)>=q);
      }
      else{
         do{
            do{
               do{
                  v1=ran1(idum);
                  v2=2.0*ran1(idum)-1.0;
               } while (v1*v1+v2*v2 > 1.0);
               y=v2/v1;
               am=ia-1.0;
               s=sqrt(2.0*am+1.0);
               x=s*y+am;
            } while(x <= 0);
            e=(1.0+y*y)*exp(am*log(x/am)-s*y);
         } while(ran1(idum) > e);
      }
   }
   return x;
}


double betadev(double a,double b, long *idum)
{
   /* gets a beta deviate */
   double xa,xb;
  // cout << " in betadev " << a << endl;
   //cout << " in betadev " << b << endl;
   xa=gamdev(a,idum);
   xb=gamdev(b,idum);
   return xa/(xa+xb);
}


/* the log(gamma) is*/
double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
      24.01409824083091,-1.231739572450155,
      0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}


double normal_dist(double x)
{

  int pdf;
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}


void choldc(double **a, int n, double p[])
{
/* this finds the cholesky decomp. returns result in lower half
   of a and the diagonal in p */
   void nrerror(char error_text[]);
   int i,j,k;
   double sum;


   // print2d(a, noPar, noPar);

   for(i=1;i<=n;i++){

      for(j=i;j<=n;j++){
         for(sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
         if (i==j){
             //cout<< "sum is "<< sum <<endl;
             if (sum<=0.0){
                 nrerror("choldc failed");
              }
             p[i]=sqrt(sum);
         } else a[j][i]=sum/p[i];
      }
   }
}


 /* routines for dealing with memory and creating vectors and matrices */

double *dvector(long nl,long nh)
/* allocate a vector with index range v[nl...nh] */
{
   void nrerror(char error_text[]);
   double *v;
   v=(double *)malloc((size_t) ((nh-nl+2)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+1;
}



double *getVector_row(double **a, int row_b, int row_e, int col_b, int col_e)
{
  double *retn;
    retn=dvector(row_e-row_b+1, col_e-col_b+1);
      for(int i=1;i<=col_e-col_b+1;i++){
          retn[i]=a[row_b][col_b-1+i];
      }
  return retn;
}


double *getVector(double **a, int row_b, int row_e, int col_b, int col_e)
{
  double *retn;
    retn=dvector(col_e-col_b+1,row_e-row_b+1);
      for(int i=1;i<=row_e-row_b+1;i++){
          retn[i]=a[row_b-1+i][col_b];
      }
  return retn;
}

void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-1));
}

int *ivector(long nl,long nh)
/* allocate a vector with index range v[nl...nh] */
{
   void nrerror(char error_text[]);
   int *v;
   v=(int *)malloc((size_t) ((nh-nl+2)*sizeof(int)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+1;
}

void free_ivector(int *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-1));
}


int **imatrix(long nrl,long nrh,long ncl,long nch)
/* allocate a int matrix with subscript range m[nrl...nrh][ncl...nch] */
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   int **m;
   void nrerror(char error_text[]);
   /* allocate pointers to rows */
   m=(int **) malloc((size_t)((nrow+1)*sizeof(int*)));
   if (!m) nrerror("allocation failure 1 in dmatrix()");
   m += 1;
   m -= nrl;
   /* allocate rows and set pointers to them */
   m[nrl]=(int *) malloc((size_t)((nrow*ncol+1)*sizeof(int)));
   if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
   m[nrl] += 1;
   m[nrl] -= ncl;
   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
   /* return pointer to array of pointers to rows */
   return m;
}


double **getMatrix(double **a, int row_b, int row_e, int col_b, int col_e)
{
  double **retn;
  retn=dmatrix(1, row_e-row_b+1, 1, col_e-col_b+1);

  for(int i=1;i<=row_e-row_b+1;i++){
    for(int j=1;j<=col_e-col_b+1;j++){
        retn[i][j]=a[row_b-1+i][col_b-1+j];
    }
  }
  return retn;
}


void free_dmatrix(double **m, long nrl,long nrh,long ncl,long nch)
/* free a double matrix allocated by dmatrix() */
{
   free((FREE_ARG) (m[nrl]+ncl-1));
   free((FREE_ARG) (m+nrl-1));
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
   cerr << "Run-time error...\n";
   cerr << error_text;
   cerr << "...now exiting to system...\n";
   exit(EXIT_FAILURE);
//   fprintf(stderr,"Run-time error...\n");
//   fprintf(stderr,"%s\n",error_text);
//   fprintf(stderr,"...now exiting to system...\n");
//   exit(1);
}


double **dmatrix(long nrl,long nrh,long ncl,long nch)
/* allocate a double matrix with subscript range m[nrl...nrh][ncl...nch] */
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;
   void nrerror(char error_text[]);
   /* allocate pointers to rows */
   m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
   if (!m) nrerror("allocation failure 1dd in dmatrix()");
   m += 1;
   m -= nrl;
   /* allocate rows and set pointers to them */
   m[nrl]=(double *) malloc((size_t)((nrow*ncol+1)*sizeof(double)));
   if (!m[nrl]) nrerror("allocation failure 2dd- in dmatrix()");
   m[nrl] += 1;
   m[nrl] -= ncl;
   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
   return m;
}

void copy_dmatrix(double **orig_mat, int rb, int re, int cb, int ce, double **copy_mat)
{
  int i, j;
  for(i=rb;i<=re;i++){
     for(j=cb;j<=ce;j++){
        copy_mat[i][j]= orig_mat[i][j];
      }
  }
}

void sort_vector(double *a,int size)
{
		double hold;
		for ( int pass=1; pass < size; pass++)
			for (int j=1; j<size; j++)
				if(a[j] > a[j+1]){
					hold=a[j];
					a[j]=a[j+1];
					a[j+1]=hold;
				}
}


void print1d(double *oned, int sizee)
{
  for(int i=1;i<=sizee;i++){
    cout << oned[i] << " ";
  }
  cout << endl;
}


void print2d(double **twod, int row, int col)
{
  for(int i=1;i<=row;i++){
    for(int j=1;j<=col;j++){
      cout << twod[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}


void print2di(int **twod, int row, int col)
{
  for(int i=1;i<=row;i++){
    for(int j=1;j<=col;j++){
      cout << twod[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}
double normal_pdf ( double x, double a, double b )

{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;
  y = ( x - a ) / b;
  pdf = exp ( -0.5 * y * y )  / ( b * sqrt ( 2.0 * pi ) );
  return pdf;
}

  double mean(double *data, int n)
{
   int j;
   double ave=0.0;
   for(j=1;j<=n;j++) ave += data[j];
   ave /= n;
   return ave;
}

double var(double *data, int n)
{
   double ave,dev,variance;
   double s=0.0;
   int j;
   ave=mean(data,n);
   for(j=1;j<=n;j++){
      dev = data[j]-ave;
      s += dev*dev;
   }
   variance=s/(n-1.0);
   return variance;
}

void get_initValu(long setSeed,int numpar, double *initValu_Mu, double **initValu_V,double *initValu)
{
   int j,k;
   double sum,*pp,*mm;
   pp=dvector(1,numpar);
   mm=dvector(1,numpar);
   choldc(initValu_V,numpar,pp);

   for(j=1;j<=numpar;j++) mm[j]=gasdev(&setSeed);
     for(j=1;j<=numpar;j++){
        sum=0.0;
        for(k=1;k<j;k++) sum += initValu_V[j][k]*mm[k];
        initValu[j] = initValu_Mu[j] + sum+pp[j]*mm[j];
     }
  
   free_dvector(mm,1,numpar);
   free_dvector(pp,1,numpar);
}


