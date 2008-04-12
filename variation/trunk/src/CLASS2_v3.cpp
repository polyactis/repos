/*
Sung Kim
version 1.0
last edit 6/2003
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <iostream.h>
#include <time.h>

int numind;
/*
routines :

  pnormal
	input : double Z
	output : double p-value

  ipnormal
	input : double p-value
	output : double Z

  pchisq
	input : double x, int df
	output : double p-value

*/

/************************************************************************
	
	pnormal( Z ) = p-value :

	  compute p_value s.t. integral [-oo, Z] of normal pdf

				                 / Z
                        1       |       -x*x/2
          p-value = ----------- |      e       dx
                    sqrt(2 pi)  |
                               /-oo

	ipnormal( p-value ) = Z :

								 / Z
                        1       |       -x*x/2
          p-value = ----------- |      e       dx
                    sqrt(2 pi)  |
                               /-oo

  

	include math.h
	include dist.h


************************************************************************/

double pnormal( double );
double ipnormal( double );
double gammp( double, double );
void gser( double, double );
void gcf( double, double );
double gammln( double );
double pchisq( double, int );
// variable for pvalue
static double gamser, gammcf, gln;

void gser(  double a, double x )
{
	int n, ITMAX=100;
	double sum, del, ap, eps = 3.0e-7;
	gln = gammln( a );
	if( x <= 0.0 )
	{
		if( x < 0.0 )
		{
			printf( "x less than 0 in routine gser" );
			exit( 1 );
		}
		gamser = 0.0;
		return;
	}
	else
	{
		ap = a;
		del = sum = 1.0 / a;

		for( n = 1; n <= ITMAX; n++ )
		{
			++ap;
			del *= x / ap;
			sum += del;
			if( fabs( del ) < fabs( sum ) * eps)
			{
				gamser = sum * exp( -x + a * log( x ) - gln );
				return;
			}
		}
		printf( "Error:a too large, ITMAX too small in routine gser" );
		exit( 1 );
	}
	return;
//	returns the incomplete gamma function P(a,x) evaluated by its series reprentation
//  as gamser, also returns ln(gamma(a) as gln;
//  the incomplete gamma function P(a,x) is efined:
//		P(a,x) = 1/gamma(a)* integral(0,x) exp(-t)t^(a-1) dt (a>0);
//	it has the limiting values P(a,0)=0 and P(a,infinity)=1;
//  it relates many cumulative probability functions, such as chi-square;
}

double gammp( double a, double x )
{

	// x = chisquare, a = df/2
	if( x < 0.0 || a <= 0.0 )
	{
		printf( "error: invalid arguments in routine gammp" );
		exit( 1 );
	}
	if( x < ( a + 1.0 ) )
	{
		gser( a, x );
		return( gamser );
	}
	else
	{
		gcf( a, x );
		return( 1.0 - gammcf );
	}
//	returns the incomplete gamma function P(a,x);
}



void gcf(  double a, double x )
{
	int i, ITMAX=100;
	double an, b, c, d, del, h, eps = 3.0e-7, fpmin = 1.0e-30;
	gln = gammln( a );
	b = x + 1.0 - a;
	c = 1.0 / fpmin;
	d = 1.0 / b;
	h = d;
	for( i = 1; i < ITMAX; i++ )
	{
		
		an = -i * ( i-a );
		b += 2.0;
		d = an * d + b;
		if( fabs( d ) < fpmin )
			d = fpmin;
		c = b + an / c;
		if( fabs( c ) < fpmin )
			c = fpmin;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if( fabs( del-1.0 ) < eps )
			break;
	}
	if( i > ITMAX )
	{
		printf( "Error:a too large, ITMAX too small in routine gser" );
		exit( 1 );
	}
	gammcf = exp( -x + a * log( x ) - gln ) * h;
//	returns the incomplete gamma function P(a,x) evaluated by its continued fraction
//	reprentation as gammcf, also returns ln(gamma(a) as gln;
//  the definition of incomplete function P(a,x) can see routine gser();
}


double gammln( double xx )
{
	double x, y, tmp, ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,
		-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	
	y = x = xx;
	tmp = x + 5.5;
	tmp -= ( x + 0.5 ) * log( tmp );
	ser = 1.000000000190015;
	for( j = 0; j < 6; j++ )
		ser += cof[ j ] / ++y;
	return( -tmp + log( 2.5066282746310005 * ser / x ) );
//	Returns the logarithmic value of gamma function of xx for xx>0;
//  the gamma function is defined by the intergral:
// 		gamma(z)=inter(0,infinity) t^(z-1)epx(-t)dt;
//  n!=gamma(n+1);
}


// return the culmulative function of standard normal distribution; 
//=========distribution of standard normal distribution
// return the probability of P(N(0,1)<x);
double pnormal( double x )
{
	if(x>=0)
		return 0.5+0.5*gammp(0.5,x*x/2.0);
	else 
		return 0.5-0.5*gammp(0.5,x*x/2.0);
}

//return the inverse of standard normal distribution
double ipnormal( double p )
{
	if(p>1.0||p<0)
	{
		printf("The P-value must be greater than 0 and less than 1!\n");
		return(0);
	}
	if(1.0-p<=1e-10)
		return(5);
	else if(p<=1e-10)
		return(-5);
	else if(fabs(p-0.5)<=1e-10)
		return(0);
	double a,b,c;
	if(p>0.5)
	{
		a=5;
		b=0;
	}
	else
	{
		a=0;
		b=-5;
	}
	while(a-b>1e-5)
	{
		c=b+(a-b)/2;
		if(pnormal(c)>p)
			a=c;
		else
			b=c;
	}
	return(c);
}

/************************************************************************
	
	computes pvalue for entered chi-square statistic 
	and degree of freedom

************************************************************************/

double pchisq( double x, int df )
{

	// df = degrees of freedom
	// x = chi square value
	return( gammp( df/2.0, x/2.0 ) );
//	return the culmulative function of chi-square distribuition
//  with df freedom;
}

class NJ
{
public:
	int * accession;
	int trace;
	double edge_c1;
	double edge_c2;
	double sumdist;

	double * distance;
	double * delta;

	int num_acc;
	int num_child;
	
	NJ ** child; 

	NJ( void );
	~NJ( );	

};

// constructor
NJ::NJ( void )
{

	int i;
	accession=new int [numind];
	for ( i=0; i<numind; i++ ){accession[i] = -99;}
	trace=-99;
	edge_c1=0;
	edge_c2=0;

	num_child=0;
	num_acc=1;

	distance = new double [numind];
	delta = new double [numind];
		
	child = new NJ * [2];
	child[0]=NULL;
	child[1]=NULL;

}

// destructor
NJ::~NJ(  )
{
	
}

/*
	distance matrix based on haplotype sharing around a snp
*/
void DM_share( int ** hap, double ** distance, int snp, int numind, int numsnp, int * snppos  )
{
	int i, k;

	int countp;
	double score;
	double norm;
	double max=0;
	int r1, r2, l1, l2;

	for ( i=0; i<numind; i++ )
	{
		for ( k=i; k<numind; k++ )
		{
			if ( i!=k )
			{
				if ( hap[i][snp]!=hap[k][snp] )
				{
					score=0;
					r1=0;
					r2=0;
					l1=0;
					l2=0;
				}
				else 
				{
					countp=snp;
					score=0;								
					// searches to right
					do 
					{
						score++;
						if ( score>0 )
						{
							countp++; // indicator for position of snp
							// when indicator exceeds total number of snps within fragment
							if ( countp>=numsnp ) 
							{
								break;
							}
						}
					}
					// continue as long as match or match to missing data
					while( hap[i][countp]==hap[k][countp] || ( hap[i][countp]==0 || hap[k][countp]==0 ) );
					
					if ( countp>=(numsnp-1) ){r2=snppos[numsnp-1];r1=snppos[numsnp-1]+1;}
					else {r2=snppos[countp];r1=snppos[countp+1];}
					score--;
					countp=snp;
					// searches to left
					do 
					{
						score++;
						if ( score>0 )
						{
							countp--; // indicator for position of snp
							// when indicator exceeds minimum number of snps within fragment
							if ( countp<0 )
							{
								break;
							}
						}
					}
					// continue as long as match or match to missing data
					while( hap[i][countp]==hap[k][countp] || ( hap[i][countp]==0 || hap[k][countp]==0 ) );

					if ( countp<0){l2=snppos[0];l1=snppos[0]-1;}
					else {l2=snppos[countp];l1=snppos[countp+1];}

				}				
				
				r1=(r1+r2)/2;
				l1=(l1+l2)/2;

				score=(r1-l1);
			//	if ( score>2000000 )score=2000000;

				if ( score<0 )
				{
					//printf( "score<0\t frag: %d snp: %d acc: %d %d\n", frag, snp, i, k );
					score=fabs(score);
				}

				distance[i][k]=score;
				distance[k][i]=score;
				if ( score>=max )
				{
					max=score;
				}
			}
		}
	}
//*
	for ( i=0; i<numind; i++ )
	{
		for ( k=0; k<numind; k++ )
		{
			if ( distance[i][k]!=0 )
			{
				distance[i][k]/=max;
				distance[i][k]=1.0-distance[i][k];
			}
			else if ( i!=k && distance[i][k]== 0 )
			{
				//distance[i][k]=max;
				distance[i][k]/=max;
				distance[i][k]=1.0-distance[i][k];
			}
			norm+=distance[i][k]*distance[i][k];
		}
	}

}



/************************************************************************
	function to quicksort -  partition
************************************************************************/
int Partition2( double ** t, int p, int r )
{
	int a, b;
	double x;

	double y1, y2, y3;
	
	x = t[p][1];
	
	a = p;

	for ( b = p+1; b <= r; b++ )
	{
		// < small to big
		if ( t[b][1] <= x )
		{
			
			a++;
			y1 = t[a][1];
			t[a][1] = t[b][1];
			t[b][1] = y1;
			
			y2 = t[a][0];
			t[a][0] = t[b][0];
			t[b][0] = y2;

			y3 = t[a][2];
			t[a][2] = t[b][2];
			t[b][2] = y3;
		}
	}
	
	y1 = t[p][1];
	t[p][1] = t[a][1];
	t[a][1] = y1;

	y2 = t[p][0];
	t[p][0] = t[a][0];
	t[a][0] = y2;

	y3 = t[p][2];
	t[p][2] = t[a][2];
	t[a][2] = y3;
	
	return a;
	
}

double ** QuickSort2( double ** t, int p, int r )
{
	int q;

	if ( p < r )
	{
		q = Partition2( t, p, r );
		
		QuickSort2( t, p, q-1 );
		QuickSort2( t, q+1, r );
	}

	return t;
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//* to be used for quantitative traits
double kw_quant ( int ** cladeslist, int * cladesusage, int * numperclade, int * cladeacc,
					  int numclades, double * pheno, int df, int clade, int numind,
					  double ** score, double * sumranks, double tie )
{
	int a,i,m,k;
	double index;
	double test_stat;
	double x, z;
	for ( i=0; i!=50; i++ ){sumranks[i]=0;}
	int * temp = new int [numind];
	for ( i=0; i<numind; i++ )
	{
		temp[i]=cladeacc[i];
	}
	for ( i=0; i<numperclade[clade]; i++ )
	{
		temp[cladeslist[clade][i]]=df;
	}

	//FILE * fp1;fp1=fopen( "cladeindicator.txt", "w");for ( i=0; i<numind; i++ ){fprintf( fp1, "%d %d %d\n", i, cladeacc[i], temp[i] );}fclose(fp1);


	// Using the rank scores computes the following kruskal-wallis statistic
	//				         flag5	    					             
	//			    12	  |	  ---       sumranks[j]*sumranks[j]   |      
	//		T =     ----------|	   \	    -----------------------   |  -   3(N+1)
	//			  N(N+1)  |        /            red_table[j][2]       |      
	//			          |	  ---				      |      
	//				          j=1							        
	
	for ( i=0; i<=df; i++ )
	{
		for ( m=0; m<numind; m++ )
		{
			if ( temp[(int)score[m][0]]==i )
			{
				sumranks[i]+=score[m][2];
			}
		}
	}
	
	index=0;
	for ( i=0; i<=df; i++ )
	{
		sumranks[i]*=sumranks[i];
		x=0;
		for ( m=0; m<numind; m++ )
		{
			if ( temp[(int)score[m][0]]==i )
			{ 
				x++;
			}
		}
		sumranks[i]/=x;
		index+=sumranks[i];
	}

	test_stat=(12.0/(double)(numind*(numind+1)));
	test_stat*=index;
	test_stat-=(double)3.0*(double)(numind+1);
	tie/=(pow(numind, 3)-numind );
	test_stat/=(1.0-tie);

	delete [] temp;
	return test_stat;
}
//*/
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/* to be used for binary traits
double chisq ( int ** cladeslist, int * cladesusage, int * numperclade, int * cladeacc,
					  int numclades, double * pheno, int df, int clade, int numind,
					  double ** score, double * sumranks, double tie )
{
	int i,j;

	double test_stat;
	double x;
	for ( i=0; i!=50; i++ ){sumranks[i]=0;}
	int * temp = new int [numind];
	for ( i=0; i<numind; i++ )
	{
		temp[i]=cladeacc[i];
	}
	for ( i=0; i<numperclade[clade]; i++ )
	{
		temp[cladeslist[clade][i]]=df;
	}

	//FILE * fp1;fp1=fopen( "cladeindicator.txt", "w");for ( i=0; i<numind; i++ ){fprintf( fp1, "%d %d %d\n", i, cladeacc[i], temp[i] );}fclose(fp1);


	double ** exp_table = new double * [df+2];
	double ** obs_table = new double * [df+2];
	for ( i=0; i<=df+1; i++ )
	{
		exp_table[i] = new double [3];
		obs_table[i] = new double [3];
		for ( j=0; j<3; j++ )
		{
			exp_table[i][j]=0;
			obs_table[i][j]=0;
		}
	}

	for ( i=0; i<numind; i++ )
	{
 		obs_table[temp[i]][(int)pheno[i]]+=1.0;
		obs_table[temp[i]][2]+=1.0;
		obs_table[df+1][(int)pheno[i]]+=1.0;
		obs_table[df+1][2]+=1.0;
//			printf( "%d ", l );
	}
//	  printf( "\n---\n" );

//	printf( "\nclades: %d\n", df+1 );
//	for ( j=0; j<3; j++ ){for ( i=0; i<=df+1; i++ ){printf( "%.0f ", obs_table[i][j] );}printf( "\n" );}
//	printf( "-----------------------------------------\n-----------------------------------------\n" );

	exp_table[df+1][0] = obs_table[df+1][0];
	exp_table[df+1][1] = obs_table[df+1][1];
	for ( i=0; i<=df+1; i++ ){exp_table[i][2]=obs_table[i][2];}
	for ( i=0; i<df+1; i++ )
	{
		for ( j=0; j<2; j++ )
		{
			exp_table[i][j]=(exp_table[i][2]*exp_table[df+1][j])/exp_table[df+1][2];
		}	
	}

	test_stat=0;
	for ( i=0; i<df+1; i++ )
	{
		for ( j=0; j<2; j++ )
		{
			x=pow((obs_table[i][j]-exp_table[i][j]), 2 );
			x/=exp_table[i][j];
			test_stat+=x;
		}
	}

	for ( i=0; i<df+2; i++ )
	{
		delete [] exp_table[i];
		delete [] obs_table[i];
	}
	delete [] exp_table;
	delete [] obs_table;
	delete [] temp;

	return test_stat;
}
//*/
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
void claderemove( int * cladesusage, int ** cladeslist, int * numperclades, int numclades, int df, int * cladeacc )
{
	int i,j,k,l;
	for ( i=0; i<numclades; i++ )
	{
		if ( cladesusage[i]==df )
		{
			for ( j=0; j<numclades; j++ )
			{
				if ( i!=j )
				{
					for ( k=0; k<numperclades[i]; k++ )
					{
						for ( l=0; l<numperclades[j]; l++ )
						{
							if ( cladeslist[i][k]==cladeslist[j][l] )
							{
								cladesusage[j]=-df;								
								break;
							}
						}
						if ( cladesusage[j]==-df ){break;}
					}
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
int cladecount( int numclades, int * cladesusage )
{
	int i,j=0;
	for ( i=0; i<numclades; i++ )
	{
		if( cladesusage[i]==0 )
		{
			j++;
		}
	}
	return j;
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
double test ( int ** hapallele, int snp, double * pheno, int numind, 
			 int ** cladeslist, int * cladesusage, int * numperclade, int numclades, int * cladeacc,
			 double ** score, double * sumranks, double tie, double ** stat )
{
	double pvalue, maxstat, temppv;
	double kwstat;
	double rawkw;
	double kw;

	int clade;
	int i;
	int df=1;
	int sigclade;
	pvalue=1;
	do 
	{
		sigclade=-9;
		maxstat=1;
		for ( clade=0; clade<numclades; clade++ )
		{
			if ( cladesusage[clade]==0 )
			{
				// test ith clade with previously found clades and its complement set
				kwstat = kw_quant ( cladeslist, cladesusage, numperclade, cladeacc, numclades,pheno, df, clade, numind, score, sumranks, tie );
				//kwstat = chisq ( cladeslist, cladesusage, numperclade, cladeacc, numclades,pheno, df, clade, numind, score, sumranks, tie );
				temppv = 1-pchisq(kwstat, df );

				if ( temppv<maxstat )
				{
					rawkw=kwstat;
					maxstat=temppv;
					sigclade=clade;
				}
			}
		}
		if ( maxstat<pvalue )
		{
			pvalue=maxstat;
			kw=rawkw;
		}
		else
		{
			break;
		}

		if ( sigclade!=-9 )
		{
			cladesusage[sigclade]=df;
			for ( i=0; i<numperclade[sigclade]; i++ )
			{
				cladeacc[cladeslist[sigclade][i]]=df;
			}
			claderemove( cladesusage, cladeslist, numperclade, numclades, df, cladeacc );
			df++;
		}
		i=i;
		
	}while( cladecount( numclades, cladesusage )!=0 );

	for ( i=0; i<numind; i++ )
	{
		if ( cladeacc[i]==0 )
		{
			cladeacc[i]=df;
		}
		hapallele[i][snp]=cladeacc[i];		
	}
	stat[snp][0]=kw;
	stat[snp][1]=df-1;

	return pvalue;

}


int main ( int argc, char* argv[] )
{

	FILE *fp1, *fp2, *fp3;

	char c = 'A', d='A';

	char * temp = new char [100];
	char * file = new char [100];

	int a, b, i, j, k, l, m;
	double x=0, y=0, z=0;

	int num_frag=1; // total number of fragment under study
	//int numind;
	int numsnp;

	file = argv[1];
	//file = "chrom4_FRI2.txt";
	fp1 = fopen( file, "r");

	fscanf( fp1, "%d", &numind );
	fscanf( fp1, "%d", &numsnp );
	//numind/=2;
	int ** hap = new int * [numind];
	int ** hapallele = new int * [numind];
	for ( i=0; i<numind; i++ )
	{
		hap[i]=new int [numsnp];
		hapallele[i]=new int [numsnp];
		for ( j=0; j<numsnp; j++ )
		{
			hap[i][j]=0;
			hapallele[i][j]=-99;
		}
	}
	int * snppos=new int [numsnp];
	double * pvalue=new double [numsnp];for ( j=0; j<numsnp; j++ ){pvalue[j]=-99;}
	double ** stat=new double * [numsnp];for ( j=0; j<numsnp; j++ ){stat[j]=new double [2];}
	double ** distance = new double * [numind];
	for ( i=0; i<numind; i++ ){distance[i]=new double [numind];for ( j=0; j<numind; j++ ){distance[i][j]=0;}}

	int ** cladeslist = new int * [numind];for ( i=0; i<numind; i++ ){cladeslist[i]=new int [numind];for ( j=0; j<numind; j++ ){cladeslist[i][j]=0;}}
	int * cladesusage = new int [numind];for ( i=0;i<numind;i++ ){cladesusage[i]=0;}
	int * numperclade = new int [numind];for ( i=0;i<numind;i++ ){numperclade[i]=0;}
	int * cladeacc = new int [numind];for ( i=0;i<numind;i++ ){cladeacc[i]=0;}
	int numclades=0;

	double * pheno = new double [numind];
	for ( i=0; i<numind; i++ )
	{
//		for ( k=0; k<2; k++)
//		{
			for ( j=0; j<numsnp; j++ )
			{
				fscanf( fp1, "%d",&l );
				hap[i][j]=l;
//				hap[i][j]+=l;
			}
//		}
		for ( j=0; j<numsnp; j++ )
		{
			if ( hap[i][j]>100 )
			{
				hap[i][j]=0;
			}
//			else
//			{
//				hap[i][j]++; // 0 : missing data
//			}
		}

	}
	for( i=0; i<numsnp; i++ )
	{
		fscanf( fp1, "%d", &snppos[i] );
	}
	for ( i=0; i<numind; i++ )
	{
		fscanf( fp1, "%lf", &pheno[i] );
	//	fscanf( fp1, "%lf", &x );
	}
	fclose( fp1 );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double ** score = new double * [numind];for ( i=0; i<numind; i++ ){score[i]=new double [3];}
	double * sumranks = new double [50];for ( i=0; i!=50; i++ ){sumranks[i]=0;}
	double tie=0;
//*
	for ( i=0; i<numind; i++ )
	{
		score[i][0]=(double)i;
		score[i][1]=(double)pheno[i];
		score[i][2]=(double)0;	
	}

	// sorts the array according to phenotype
	score = QuickSort2( score, 0, numind-1 );			
	// given sorted values, determines rank score.
	// if phenotypes are equal, average of rank scores are computed
	for ( i=0; i<numind; i++ )
	{
		if ( i==(numind-1) )
		{
			score[i][2]=(double)i+1;
			break;
		}
		if ( score[i][1]!=score[i+1][1] && i<numind )
		{
			score[i][2]=(double)i+1;
		}
		else
		{
			z=(double)i+1;
			a=i;
			m=1;
			while( score[i][1]==score[i+1][1] && i<=(numind-2) )
			{
				m++;
				i++;
				z+=(double)i+1;
				if ( i==(numind-1) )
				{
					break;
				}
			}
			z/=(double)m;
			if ( m>1 )
			{
				tie+=(pow(m,3)-m);
			}

			for ( k=0; k!=m; k++ )
			{
				score[a+k][2]=z;
			}
		}
	}

	sprintf( temp, "%s_phenranks.txt", file );	
	fp1=fopen( temp, "w" );
	for ( i=0; i<numind; i++ )
	{
		fprintf( fp1, "%.0f\t%f\t%f\n", score[i][0], score[i][1], score[i][2] );
	}
	fclose( fp1 );
//*/

	int numacc=i;
	int level=0;
	int index=0;
	double ave_dist1, ave_dist2;
	double R1, R2;
	double min_dist;
	int child1, child2;
	int num_taxa;
	int snp;

	for ( snp=0; snp<numsnp; snp++ )
	//for ( snp=6; snp<7; snp++ )
	{

		for ( i=0; i<numind; i++ ){for ( j=0; j<numind; j++ ){cladeslist[i][j]=0;}}
		for ( i=0;i<numind;i++ ){cladesusage[i]=0;}
		for ( i=0;i<numind;i++ ){numperclade[i]=0;}
		for ( i=0;i<numind;i++ ){cladeacc[i]=0;}
		int numclades=0;
		
		NJ ** data2; // entire dataset
		data2 = new NJ * [2];
		for ( i=0; i<2; i++ )
		{
			data2[i]=new NJ [numind+1-i];
		}
		
		for ( i=0; i<numind; i++ )
			for ( j=0; j<numind; j++ )			
				distance[i][j]=0;
					
		// computes a distance metric based on scoring scheme:
		// match : +1
		// mismatch : -2
		// aligned to missing data: 0
		// computes a distance metric based on haplotype sharing
		DM_share( hap, distance, snp, numind, numsnp, snppos );

		//sprintf( temp, "%s_clades_%d.txt", file, snppos[snp] );				
		//fp2=fopen( temp, "w" );
//*
		sprintf( temp, "%s_dist_%d.txt", file, snppos[snp] );
		fp3=fopen( temp, "w" );

		for ( i=0; i<numind; i++ )
		{
			fprintf( fp3, "%d ", i );
			for ( j=0; j<numind; j++ )
			{
				fprintf( fp3, "%.3f ", distance[i][j] );
			}
			fprintf( fp3, "\n" );
		}
		fclose( fp3 );
//*/
		// transfers information from distance matrix to class data2
		num_taxa=numind;
		for ( i=0; i<num_taxa; i++ )
		{
			for (j=0; j<num_taxa ; j++ )
			{
				if ( j==0 )
				{
					data2[0][i].accession[0]=i;
					data2[0][i].trace=i;
				}
				data2[0][i].distance[j]=distance[i][j];
				data2[0][i].child[0]=NULL;
				data2[0][i].child[1]=NULL;
				data2[0][i].delta[j]=0;
				data2[0][i].num_acc=1;
				data2[0][i].num_child=0;	
				data2[0][i].edge_c1=0;
				data2[0][i].edge_c2=0;

			}
		}

		level=0;
		// initialization step for neighbor joining
		// compute corrected distance for all pairs i and j
		// the corrected distance is the pairwise distance between i and j subtracted by 
		// Ri and Rj, where Ri and R j are the averaged distance to all other accessions for 
		// i and j, respectively
		min_dist=1e8;
		for ( i=0; i<num_taxa; i++ )
		{
			for ( j=(i+1); j<num_taxa; j++ )
			{
				ave_dist1=0;
				ave_dist2=0;
				// averaged distance of i to all other non j nodes
				for( a=0; a<num_taxa; a++ )
				{
					if ( a!=i && a!=j )
					{
						ave_dist1+=data2[0][i].distance[a];
					}
				}
				ave_dist1/=(num_taxa-2);
				// averaged distance of j to all other non i nodes
				for( a=0; a<num_taxa; a++ )
				{
					if ( a!=i && a!=j )
					{
						ave_dist2+=data2[0][j].distance[a];
					}
				}
				ave_dist2/=(num_taxa-2);

				// computes the corrected distance
				data2[0][i].delta[j]=data2[0][i].distance[j]-ave_dist1-ave_dist2;
				
				//finds the pair that gives the minimum corrected distance
				if ( data2[0][i].delta[j] <= min_dist )
				{
					min_dist = data2[0][i].delta[j];
					child1=i;
					child2=j;
					R1=ave_dist1;
					R2=ave_dist2;
				}
			}
		}

		// initial output of each accession
		// repeat until number of accession is 2
		while ( num_taxa>2 )
		{	
			// reference index for first position of the pair that gave the minimum corrected distance
			index=child1;
			if ( child1>child2 )
			{
				index=child2;
				j=child1;
			}
			
			// data2 merge for next level, i.e. takes all accessions and merges the
			// minimum corrected distance pair
			for ( i=0; i<num_taxa; i++ )
			{
				// when a mergeable cell is found, 
				if ( child1==i || child2==i )
				{
					if ( i==index )
					{
						// first child, i.e. mergeable cell
						// transfer address to parent cell
						data2[1][index].child[0] = &data2[0][i];
						data2[1][index].num_child=2;
						// computes edge of parent cell to child 1
						data2[1][index].edge_c1=(data2[0][i].distance[j]+R1-R2)/2;
						// keeps track of an accession as it moves along the tree
						data2[1][i].trace=i;
						// stores information of accessions being merged
						
						for ( a=0; a<data2[0][index].num_acc; a++ )
						{
							data2[1][index].accession[a]=data2[0][i].accession[a];
						}
						data2[1][index].num_acc=data2[0][i].num_acc;
					}
					else
					{
						// found second child, i.e. mergeable cell
						data2[1][index].child[1] = &data2[0][i];
						data2[1][index].edge_c2=( data2[0][i].distance[i]-data2[1][index].edge_c1 );
					
						// stores information of accessions being merged
						for ( b=0; b<data2[0][i].num_acc; b++ )
						{
							data2[1][index].accession[data2[0][index].num_acc+b]=data2[0][i].accession[b];

						}
						data2[1][index].num_acc+=data2[0][i].num_acc;

						if ( data2[1][index].num_acc>=5 && data2[1][index].num_acc<numind-5 )
						{							
							numperclade[numclades]=data2[1][index].num_acc;
							for ( a=0; a<data2[1][index].num_acc; a++ )
							{
								//fprintf( fp2, "%d ", data2[1][index].accession[a] );
								cladeslist[numclades][a]=data2[1][index].accession[a];
							}
							//fprintf( fp2, "\n" );
							numclades++;
						}

						// shift adjustment to account for removed later child, just passing address of previous node
						for ( b=i; b<num_taxa; b++ )
						{
							data2[1][b].child[0] = &data2[0][b+1];
							data2[1][b].num_child=1;
							data2[1][b].trace=b;
							for ( a=0; a<data2[0][b+1].num_acc; a++ )
							{
								data2[1][b].accession[a]=data2[0][b+1].accession[a];
							}
							data2[1][b].num_acc=data2[0][b+1].num_acc;
						}						
						i=num_taxa;	
						
					}				
				}
				// when no mergeable child is found, just pass address of previous node
				else
				{
					data2[1][i].child[0] = &data2[0][i];
					data2[1][i].num_child=1;
					data2[1][i].trace=i;
					for ( a=0; a<data2[0][i].num_acc; a++ )
					{
						data2[1][i].accession[a]=data2[0][i].accession[a];
					}
					data2[1][i].num_acc=data2[0][i].num_acc;
				}
			}


			num_taxa--;
			// generate new distance matrix 
			for ( i=0; i<num_taxa; i++ )
			{
				for (j=i; j<num_taxa ; j++ )
				{
					if ( i==j )
					{
						data2[1][i].distance[j]=0;
					}
					// when either accession i or j matches minimum pair k and l that formed index,
					// then must recalculate distance between i or j to index adjusting for the merge
					// in other words, to compute the distance between i=index and j
					// add the distance(j, k) and distance(j, l)
					// subtract with distance(k, l) and divide by 2
					else if ( i==index && j!=index )
					{
						data2[1][i].distance[j] = 
							( 
							data2[1][i].child[0]->distance[data2[1][j].child[0]->trace] +
							data2[1][i].child[1]->distance[data2[1][j].child[0]->trace] -
							data2[1][i].child[0]->distance[data2[1][i].child[1]->trace] 
							)/2.0;
					}
					else if ( i!=index && j==index )
					{
						data2[1][i].distance[j] = 
							( 
							data2[1][j].child[0]->distance[data2[1][i].child[0]->trace] +
							data2[1][j].child[1]->distance[data2[1][i].child[0]->trace] -
							data2[1][j].child[0]->distance[data2[1][j].child[1]->trace] 
							)/2.0;
					}
					// when neither of accesion i or j are from the minimum pair,
					// then the distance remains unchanged between i and j
					else if ( i!=index && j!=index )
					{
						data2[1][i].distance[j] = 
						data2[1][i].child[0]->distance[data2[1][j].child[0]->trace];						
					}
					data2[1][j].distance[i]=data2[1][i].distance[j];
				}
			}

			for ( i=0; i<num_taxa; i++ )
			{
				data2[0][i]=data2[1][i];
			}
		
			//finds the minimum pair for new distances computed after merging
			min_dist=1e8;
			//finds the minimum pair for new distances computed after merging
			for ( i=0; i<num_taxa; i++ ){data2[level][i].sumdist=0;for ( j=0; j<num_taxa; j++ ){data2[level][i].sumdist+=data2[level][i].distance[j];}}
			for ( i=0; i<num_taxa; i++ )
			{
				for ( j=(i+1); j<num_taxa; j++ )
				{
					// averaged distance of i to all other non j nodes
					ave_dist1=data2[level][i].sumdist;
					ave_dist1-=data2[level][i].distance[i];
					ave_dist1-=data2[level][i].distance[j];
					ave_dist1/=(num_taxa-2);

					// averaged distance of j to all other non i nodes
					ave_dist2=data2[level][j].sumdist;
					ave_dist2-=data2[level][j].distance[i];
					ave_dist2-=data2[level][j].distance[j];
					ave_dist2/=(num_taxa-2);

					data2[level][i].delta[j]=data2[level][i].distance[j]-ave_dist1-ave_dist2;

					if ( data2[level][i].delta[j] <= min_dist )
					{
						min_dist = data2[level][i].delta[j];
						child1=i;
						child2=j;
						R1=ave_dist1;
						R2=ave_dist2;
					}
				}
			}
			
		}
		// fclose( fp2 );

		pvalue[snp] = test( hapallele, snp, pheno, numind, 
			cladeslist, cladesusage, numperclade, numclades, cladeacc, 
			score, sumranks, tie, stat );
	
		printf( "%d %f\n", snp, -log(pvalue[snp]) );			
	
		sprintf( temp, "%s_results.txt", file);
		fp3=fopen( temp, "a" );
		fprintf( fp3, "%d %.20f %f %.10f %.0f\n", snp, pvalue[snp], -log(pvalue[snp]), stat[snp][0], stat[snp][1] );
		fclose( fp3 );

		for ( i=0; i<2; i++ )
		{
			delete [] data2[i];
		}
		delete []data2;

	}

	sprintf( temp, "%s_results.txt", file);
	fp3=fopen( temp, "w" );
	for ( i=0; i<numsnp; i++ )
	{
		fprintf( fp3, "%d %.20f %f %.10f %.0f\n", i, pvalue[i], -log(pvalue[i]), stat[i][0], stat[i][1] );
	}
	fclose( fp3 );

	sprintf( temp, "%s_hap_results.txt", file );
	fp3=fopen( temp, "w" );
	for ( i=0; i<numind; i++ )
	{
		fprintf( fp3, "%d ", i );
		for ( j=0; j<numsnp; j++ )
		{
			fprintf( fp3, "%d ", hapallele[i][j] );
		}
		fprintf( fp3, "\n" );

	}
	fclose( fp3 );

	sprintf( temp, "%s_haps.txt", file);
	fp3=fopen( temp, "w" );
	for ( i=0; i<numind; i++ )
	{
		fprintf( fp3, "%d ", i );
		for ( j=0; j<numsnp; j++ )
		{
			fprintf( fp3, "%d ", hap[i][j] );
		}
		fprintf( fp3, "%.0f ", pheno[i] );
		fprintf( fp3, "\n" );

	}
	fclose( fp3 );

	

 	return 0;
}

