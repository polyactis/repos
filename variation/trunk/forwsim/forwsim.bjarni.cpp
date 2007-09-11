#include"mtrand.h"
#include<iostream.h>
#include<math.h>
#include<stdio.h>
#include<fstream.h>
#include<time.h>


#define      SSIZE            100    // Number of samples to output from the final population.
#define    POPSIZE            20000    // Total number of chromsomes in the diploid population.  (Has to be divisible by two)
#define  SEQLENGTH          50000    // Length of the DNA sequence in basepairs.Choose a large enough value >= DELETE*mu*POPSIZE so that you dont run out of locations under the infinite-sites model.
#define         re            0.001    // Per-generation per-sequence rate of recombination.
#define         mu            0.001    // Per-generation per-sequence rate of mutation.
//#define          s            0.00    // Probability of self-fertilization.
#define         ns            0.95    // Probability of self-fertilization.
#define        GEN           10000    // Number of generations.
#define    SELFGEN            5000    // Number of generations before the selfing rate changes.
#define     DELETE             500    // Intervals after which fixed mutations are removed.  (Should be =GEN to get an exact algorithm?)
#define       ZERO               0    // Constant. Do not change.

double s = 0.00;    // Probability of self-fertilization.


unsigned long init[4] = {0x123, 0x234, 0x345, 0x456},length = 4;
MTRand drand;//double in [0, 1) generator, already init
MTRand_int32 irand(init, length);// 32-bit int generator
MTRand_int32 seed (init, length);// 32-bit int generator
MTRand mt;
int i,j,k,l,m,n,o,r,trp,trk,start,end,ind;
int n2count[POPSIZE],n2r[POPSIZE],n3count[POPSIZE],n3r[POPSIZE];
int n4count[POPSIZE],n4r[POPSIZE],n5count[POPSIZE],n5r[POPSIZE];
int n6count[POPSIZE],n6r[POPSIZE],n7count[POPSIZE],n7r[POPSIZE];
int n8count[POPSIZE],n8r[POPSIZE],n9count[POPSIZE],n9r[POPSIZE];
int current[POPSIZE],size[POPSIZE],siz[POPSIZE];
int temp1[POPSIZE],temp2[POPSIZE],recount[POPSIZE];
int indic[POPSIZE],flag1,flag2,indicator[SEQLENGTH],multhit[SEQLENGTH];
double p,rec,x,y;

static inline double  minimum(double i,double j){return(i>j) ? j:i;}
static inline double  maximum(double i,double j){return(i>j) ? i:j;}
static inline int poisson(float mean){int count=0;double q,bound;bound = exp(-mean);
	for(q=1;q>=bound;q*=drand()){++count;}
	return(count-1);
}

class next{
public:
	double np[POPSIZE];int n[POPSIZE];
};next *gen1,*gen2,*gen3,*gen4,*gen5,*gen6,*gen7,*gen8,*gen9,*gtemp;

class member{
public:
	int q[25000 + int(GEN*mu)];//Choose a suitable value depending on the expected number of mutations per individual
};member *temp;

struct node{
	member *chr;
};node odd[POPSIZE],even[POPSIZE];




main(){
	rec = 1 - exp(-re);for(i=0;i<length;i++){init[i]=time(NULL);}
	mt.seed(init,length);
	gen1 = new next;gen2 = new next;gen3 = new next;gen4 = new next;gen5 = new next;gen6 = new next;gen7 = new next;gen8 = new next;
	gen9 = new next;gtemp = new next;temp = new member;

	for(i=0;i<POPSIZE;i++){odd[i].chr = new member;even[i].chr = new member;}
	for(i=0;i<POPSIZE;i++){size[i]=0;siz[i]=0;}
	for(i=0;i<SEQLENGTH;i++){multhit[i]=0;}

	m=0;
	while(m<POPSIZE){  //Initializing without recombination!!
		gen1->np[m] = drand();gen1->np[m+1] = drand();x=drand();
		if(x>=s){gen1->n[m]=int(POPSIZE*drand());gen1->n[m+1]=int(POPSIZE*drand());}  //Non-selfing
		else{y=drand();if(y<=0.50){gen1->n[m]=int(POPSIZE*drand());gen1->n[m+1]= gen1->n[m];}  //Selfing
			else{gen1->n[m]=int(POPSIZE*drand());if((gen1->n[m])%2==0){gen1->n[m+1]=(gen1->n[m])+1;}
				if((gen1->n[m])%2==1){gen1->n[m+1]=(gen1->n[m])-1;}

			}//Selfing

		}

		gen2->np[m] = drand();gen2->np[m+1] = drand();x=drand();
		if(x>=s){gen2->n[m]=int(POPSIZE*drand());gen2->n[m+1]=int(POPSIZE*drand());}
		else{y=drand();if(y<=0.50){gen2->n[m]=int(POPSIZE*drand());gen2->n[m+1]= gen2->n[m];}
			else{gen2->n[m]=int(POPSIZE*drand());if((gen2->n[m])%2==0){gen2->n[m+1]=(gen2->n[m])+1;}
				if((gen2->n[m])%2==1){gen2->n[m+1]=(gen2->n[m])-1;}

			}

		}

		gen3->np[m] = drand();gen3->np[m+1] = drand();x=drand();
		if(x>=s){gen3->n[m]=int(POPSIZE*drand());gen3->n[m+1]=int(POPSIZE*drand());}
		else{y=drand();if(y<=0.50){gen3->n[m]=int(POPSIZE*drand());gen3->n[m+1]= gen3->n[m];}
			else{gen3->n[m]=int(POPSIZE*drand());if((gen3->n[m])%2==0){gen3->n[m+1]=(gen3->n[m])+1;}
				if((gen3->n[m])%2==1){gen3->n[m+1]=(gen3->n[m])-1;}

			}

		}

		gen4->np[m] = drand();gen4->np[m+1] = drand();x=drand();
		if(x>=s){gen4->n[m]=int(POPSIZE*drand());gen4->n[m+1]=int(POPSIZE*drand());}
		else{y=drand();if(y<=0.50){gen4->n[m]=int(POPSIZE*drand());gen4->n[m+1]= gen4->n[m];}
			else{gen4->n[m]=int(POPSIZE*drand());if((gen4->n[m])%2==0){gen4->n[m+1]=(gen4->n[m])+1;}
				if((gen4->n[m])%2==1){gen4->n[m+1]=(gen4->n[m])-1;}

			}

		}

		gen5->np[m] = drand();gen5->np[m+1] = drand();x=drand();
		if(x>=s){gen5->n[m]=int(POPSIZE*drand());gen5->n[m+1]=int(POPSIZE*drand());}
		else{y=drand();if(y<=0.50){gen5->n[m]=int(POPSIZE*drand());gen5->n[m+1]= gen5->n[m];}
			else{gen5->n[m]=int(POPSIZE*drand());if((gen5->n[m])%2==0){gen5->n[m+1]=(gen5->n[m])+1;}
				if((gen5->n[m])%2==1){gen5->n[m+1]=(gen5->n[m])-1;}

			}

		}

		gen6->np[m] = drand();gen6->np[m+1] = drand();x=drand();
		if(x>=s){gen6->n[m]=int(POPSIZE*drand());gen6->n[m+1]=int(POPSIZE*drand());}
		else{y=drand();if(y<=0.50){gen6->n[m]=int(POPSIZE*drand());gen6->n[m+1]= gen6->n[m];}
			else{gen6->n[m]=int(POPSIZE*drand());if((gen6->n[m])%2==0){gen6->n[m+1]=(gen6->n[m])+1;}
				if((gen6->n[m])%2==1){gen6->n[m+1]=(gen6->n[m])-1;}

			}

		}

		gen7->np[m] = drand();gen7->np[m+1] = drand();x=drand();
		if(x>=s){gen7->n[m]=int(POPSIZE*drand());gen7->n[m+1]=int(POPSIZE*drand());}
		else{y=drand();if(y<=0.50){gen7->n[m]=int(POPSIZE*drand());gen7->n[m+1]= gen7->n[m];}
			else{gen7->n[m]=int(POPSIZE*drand());if((gen7->n[m])%2==0){gen7->n[m+1]=(gen7->n[m])+1;}
				if((gen7->n[m])%2==1){gen7->n[m+1]=(gen7->n[m])-1;}

			}

		}

		gen8->np[m] = drand();gen8->np[m+1] = drand();x=drand();
		if(x>=s){gen8->n[m]=int(POPSIZE*drand());gen8->n[m+1]=int(POPSIZE*drand());}
		else{y=drand();if(y<=0.50){gen8->n[m]=int(POPSIZE*drand());gen8->n[m+1]= gen8->n[m];}
			else{gen8->n[m]=int(POPSIZE*drand());if((gen8->n[m])%2==0){gen8->n[m+1]=(gen8->n[m])+1;}
				if((gen8->n[m])%2==1){gen8->n[m+1]=(gen8->n[m])-1;}

			}

		}

		gen9->np[m] = drand();gen9->np[m+1] = drand();x=drand();
		if(x>=s){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]=int(POPSIZE*drand());}
		else{y=drand();if(y<=0.50){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]= gen9->n[m];}
			else{gen9->n[m]=int(POPSIZE*drand());if((gen9->n[m])%2==0){gen9->n[m+1]=(gen9->n[m])+1;}
				if((gen9->n[m])%2==1){gen9->n[m+1]=(gen9->n[m])-1;}

			}

		}

		n2count[m]=0;n2r[m]=0;n2count[m+1]=0;n2r[m+1]=0;
		n3count[m]=0;n3r[m]=0;n3count[m+1]=0;n3r[m+1]=0;
		n4count[m]=0;n4r[m]=0;n4count[m+1]=0;n4r[m+1]=0;
		n5count[m]=0;n5r[m]=0;n5count[m+1]=0;n5r[m+1]=0;
		n6count[m]=0;n6r[m]=0;n6count[m+1]=0;n6r[m+1]=0;
		n7count[m]=0;n7r[m]=0;n7count[m+1]=0;n7r[m+1]=0;
		n8count[m]=0;n8r[m]=0;n8count[m+1]=0;n8r[m+1]=0;
		n9count[m]=0;n9r[m]=0;n9count[m+1]=0;n9r[m+1]=0;
		m=m+2;
	}


	//Start main for-loop
	for(i=1;i<=GEN;i++){

		//Change of selfing.
		if (i==SELFGEN){
			s = 0.95;
		}

		//Create odd generation if statement  (odd-even categorisation needed for pointer swithcing?)
		if(i%2 == 1){

			//What to remember at current generation looking at future ones!
			for(m=0;m<POPSIZE;m++){ //(For every m!)
				temp1[m] = gen2->n[m];n2count[temp1[m]] = 1;
				if(gen2->np[m]<rec){                          //Recombination   (It follows a geometric distribution (a discrete counterpart of the exp. dist.).)
					if(gen2->n[m]%2==0){n2r[gen2->n[m]+1] = 1;}
					if(gen2->n[m]%2==1){n2r[gen2->n[m]-1] = 1;}


				}
				current[m] = m;size[m]=0;indic[m]=0;recount[m]=0; //current keeps count of which genes are segregating in the current population?
			}

			for(m=0;m<POPSIZE;m++){
				temp2[m] = temp1[gen3->n[m]];n3count[temp2[m]] = 1;
				if(gen3->np[m]<rec){
					if(gen3->n[m]%2==0){n3r[temp1[gen3->n[m]+1]] = 1;}
					if(gen3->n[m]%2==1){n3r[temp1[gen3->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp1[m] = temp2[gen4->n[m]];n4count[temp1[m]] = 1;
				if(gen4->np[m]<rec){
					if(gen4->n[m]%2==0){n4r[temp2[gen4->n[m]+1]] = 1;}
					if(gen4->n[m]%2==1){n4r[temp2[gen4->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp2[m] = temp1[gen5->n[m]];n5count[temp2[m]] = 1;
				if(gen5->np[m]<rec){
					if(gen5->n[m]%2==0){n5r[temp1[gen5->n[m]+1]] = 1;}
					if(gen5->n[m]%2==1){n5r[temp1[gen5->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp1[m] = temp2[gen6->n[m]];n6count[temp1[m]] = 1;
				if(gen6->np[m]<rec){
					if(gen6->n[m]%2==0){n6r[temp2[gen6->n[m]+1]] = 1;}
					if(gen6->n[m]%2==1){n6r[temp2[gen6->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp2[m] = temp1[gen7->n[m]];n7count[temp2[m]] = 1;
				if(gen7->np[m]<rec){
					if(gen7->n[m]%2==0){n7r[temp1[gen7->n[m]+1]] = 1;}
					if(gen7->n[m]%2==1){n7r[temp1[gen7->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp1[m] = temp2[gen8->n[m]];n8count[temp1[m]] = 1;
				if(gen8->np[m]<rec){
					if(gen8->n[m]%2==0){n8r[temp2[gen8->n[m]+1]] = 1;}
					if(gen8->n[m]%2==1){n8r[temp2[gen8->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp2[m] = temp1[gen9->n[m]];n9count[temp2[m]] = 1;
				if(gen9->np[m]<rec){
					if(gen9->n[m]%2==0){n9r[temp1[gen9->n[m]+1]] = 1;}
					if(gen9->n[m]%2==1){n9r[temp1[gen9->n[m]-1]] = 1;}


				}
			}


			// While loop starts
			m=0;
			while(m<POPSIZE){
				flag1=0;
				while(ZERO==0){
					if(n2count[m]==0){flag1=1;break;}
					if(n3r[m]==0){
						if(n3count[m]==0){flag1=1;break;}
						if(n4r[m]==0){
							if(n4count[m]==0){flag1=1;break;}
							if(n5r[m]==0){
								if(n5count[m]==0){flag1=1;break;}
								if(n6r[m]==0){
									if(n6count[m]==0){flag1=1;break;}
									if(n7r[m]==0){
										if(n7count[m]==0){flag1=1;break;}
										if(n8r[m]==0){
											if(n8count[m]==0){flag1=1;break;}
											if(n9r[m]==0){
												if(n9count[m]==0){flag1=1;break;}


											}
										}
									}
								}
							}
						}
					}
					break;
				}

				flag2=0;
				while(ZERO==0){
					if(n2count[m+1]==0){flag2=1;break;}
					if(n3r[m+1]==0){
						if(n3count[m+1]==0){flag2=1;break;}
						if(n4r[m+1]==0){
							if(n4count[m+1]==0){flag2=1;break;}
							if(n5r[m+1]==0){
								if(n5count[m+1]==0){flag2=1;break;}
								if(n6r[m+1]==0){
									if(n6count[m+1]==0){flag2=1;break;}
									if(n7r[m+1]==0){
										if(n7count[m+1]==0){flag2=1;break;}
										if(n8r[m+1]==0){
											if(n8count[m+1]==0){flag2=1;break;}
											if(n9r[m+1]==0){
												if(n9count[m+1]==0){flag2=1;break;}


											}
										}
									}
								}
							}
						}
					}
					break;
				}

				if(flag1==1 && flag2==1){indic[m]=1;indic[m+1]=1;} // The chromasomes are not required in future generations.
				else{
					if(flag1==1 && n2r[m]==0){indic[m]=1;}
					if(flag2==1 && n2r[m+1]==0){indic[m+1]=1;}


				}

				//m case
				if(indic[m]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){  //WHAT IS DELETE??
					p = gen1->np[m];j = gen1->n[m];  // j denotes the index of the chromosome at pos m in 1 gen.
					if(p < rec){if(j%2==0){k = j+1;}
						if(j%2==1){k = j-1;}
						++recount[(j+k)/2];
					}
					else{
						if(current[j]!=j){memcpy(&odd[m].chr->q[0],&odd[current[j]-POPSIZE].chr->q[0],4*siz[j]);size[m]=siz[j];}//2. Else copy stuff.
						if(current[j]==j){temp = odd[m].chr;odd[m].chr = even[j].chr;even[j].chr = temp;current[j]=m + POPSIZE;size[m] = siz[j];}  //1. Move pointers if possible


					}
				}

				//m+1 case
				if(indic[m+1]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
					p = gen1->np[m+1];j = gen1->n[m+1];
					if(p < rec){if(j%2==0){k = j+1;}
						if(j%2==1){k = j-1;}
						++recount[(j+k)/2];
					}
					else{
						if(current[j]!=j){memcpy(&odd[m+1].chr->q[0],&odd[current[j]-POPSIZE].chr->q[0],4*siz[j]);size[m+1]=siz[j];}
						if(current[j]==j){temp = odd[m+1].chr;odd[m+1].chr = even[j].chr;even[j].chr = temp;current[j]=m+1+POPSIZE;size[m+1] = siz[j];}


					}
				}
				m=m+2;
			}
			// While loop ends


			// for-loop starts
			for(m=0;m<POPSIZE;m++){

				p = gen1->np[m];j = gen1->n[m];// Again j denotes the index of the chromosome at pos m in 1 gen.  And p is the associated random val.

				if(i%DELETE > 0  && i%DELETE < DELETE-10 && i < GEN-10){if(indic[m]==1){continue;}

				}

				if(p<rec){
					r = int((SEQLENGTH-1)*drand());
					if(j%2==0){k = j+1;}
					if(j%2==1){k = j-1;}
					--recount[(j+k)/2];

					if(current[j]==j && siz[j] > 0){
						if(r>=even[j].chr->q[siz[j]-1]){
							if(recount[(j+k)/2]==0){temp=odd[m].chr;odd[m].chr=even[j].chr;even[j].chr=temp;size[m]=siz[j];}
							else{memcpy(&odd[m].chr->q[0],&even[j].chr->q[0],4*siz[j]);size[m]=siz[j];}


						}

						else{
							start=0;
							end=siz[j]-1;
							trk=end;
							ind=0;
							while(trk > 0 && trk <=siz[j]-1){
								if(r >=even[j].chr->q[trk-1] && r <even[j].chr->q[trk]){ind=1;break;}
								if(r >=even[j].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(r <even[j].chr->q[trk]){end=trk;trk=(start+end)/2;}


							}
							if(ind==1){
								if(recount[(j+k)/2]==0){temp=odd[m].chr;odd[m].chr=even[j].chr;even[j].chr=temp;size[m]=trk;}
								else{memcpy(&odd[m].chr->q[0],&even[j].chr->q[0],4*trk);size[m]=trk;}


							}
						}
					}

					if(current[j]!=j && siz[j] > 0){
						trp=current[j]-POPSIZE;
						if(r>=odd[trp].chr->q[siz[j]-1]){memcpy(&odd[m].chr->q[0],&odd[trp].chr->q[0],4*siz[j]);size[m]=siz[j];}
						else{
							start=0;end=siz[j]-1;trk=end;ind=0;
							while(trk > 0 && trk<=siz[j]-1){
								if(r >=odd[trp].chr->q[trk-1] && r <odd[trp].chr->q[trk]){ind=1;break;}
								if(r >=odd[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(r <odd[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}


							}
							if(ind==1){memcpy(&odd[m].chr->q[0],&odd[trp].chr->q[0],4*trk);size[m]=trk;}


						}
					}

					if(current[k]==k && siz[k] > 0){
						if(r<even[k].chr->q[0]){memcpy(&odd[m].chr->q[size[m]],&even[k].chr->q[0],4*siz[k]);size[m]=size[m]+siz[k];}
						else{
							start=0;end=siz[k]-1;trk=end;ind=0;
							while(trk > 0 && trk <=siz[k]-1){
								if(r >=even[k].chr->q[trk-1] && r <even[k].chr->q[trk]){ind=1;break;}
								if(r >=even[k].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(r <even[k].chr->q[trk]){end=trk;trk=(start+end)/2;}


							}
							if(ind==1){memcpy(&odd[m].chr->q[size[m]],&even[k].chr->q[trk],4*(siz[k]-trk));size[m]=size[m]+siz[k]-trk;}


						}
					}

					if(current[k]!=k && siz[k] > 0){trp=current[k]-POPSIZE;
						if(r<odd[trp].chr->q[0]){memcpy(&odd[m].chr->q[size[m]],&odd[trp].chr->q[0],4*siz[k]);size[m]=size[m]+siz[k];}
						else{
							start=0;end=siz[k]-1;trk=end;ind=0;
							while(trk > 0 && trk<=siz[k]-1){
								if(r >=odd[trp].chr->q[trk-1] && r <odd[trp].chr->q[trk]){ind=1;break;}
								if(r >=odd[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(r <odd[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}


							}
							if(ind==1){memcpy(&odd[m].chr->q[size[m]],&odd[trp].chr->q[trk],4*(siz[k]-trk));size[m]=size[m]+siz[k]-trk;}


						}
					}
				}
			}
			// for-loop ends




			for(m=0;m<POPSIZE;m++){

				if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m]==1){continue;}

				}

				//INSERT NEW MUTATIONS
				k = poisson(mu);
				for(l=0;l<k;l++){
					trp = int(SEQLENGTH*drand());while(multhit[trp]==1){trp = int(SEQLENGTH*drand());}
					multhit[trp]=1;
					if((size[m] > 0 && trp>=odd[m].chr->q[size[m]-1]) || size[m]==0){odd[m].chr->q[size[m]]=trp;}
					else{
						if(trp<=odd[m].chr->q[0]){memmove(&odd[m].chr->q[1],&odd[m].chr->q[0],4*size[m]);odd[m].chr->q[0]=trp;}
						else{
							start=0;end=size[m]-1;trk=end;ind=0;
							while(trk > 0 && trk<=size[m]-1){
								if(trp > odd[m].chr->q[trk-1] && trp <=odd[m].chr->q[trk]){ind=1;break;}
								if(trp > odd[m].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(trp <= odd[m].chr->q[trk]){end=trk;trk=(start+end)/2;continue;}


							}
							if(ind==1){memmove(&odd[m].chr->q[trk+1],&odd[m].chr->q[trk],4*(size[m]-trk));odd[m].chr->q[trk]=trp;}


						}
					}++size[m];
				}
			}


			gtemp = gen1;gen1 = gen2;gen2 = gen3;gen3 = gen4;gen4 = gen5;gen5 = gen6;gen6 = gen7;gen7 = gen8;gen8=gen9;gen9 = gtemp;
			m=0;
			while(m<POPSIZE){
				gen9->np[m] = drand();gen9->np[m+1] = drand();x=drand();
				if(x>=s){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]=int(POPSIZE*drand());}  //Selfing
				else{y=drand();if(y<=0.50){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]= gen9->n[m];}
					else{gen9->n[m]=int(POPSIZE*drand());if((gen9->n[m])%2==0){gen9->n[m+1]=(gen9->n[m])+1;}
						else{gen9->n[m+1]=(gen9->n[m])-1;}


					}
				}
				n2count[m]=0;n2r[m]=0;n2count[m+1]=0;n2r[m+1]=0;
				n3count[m]=0;n3r[m]=0;n3count[m+1]=0;n3r[m+1]=0;
				n4count[m]=0;n4r[m]=0;n4count[m+1]=0;n4r[m+1]=0;
				n5count[m]=0;n5r[m]=0;n5count[m+1]=0;n5r[m+1]=0;
				n6count[m]=0;n6r[m]=0;n6count[m+1]=0;n6r[m+1]=0;
				n7count[m]=0;n7r[m]=0;n7count[m+1]=0;n7r[m+1]=0;
				n8count[m]=0;n8r[m]=0;n8count[m+1]=0;n8r[m+1]=0;
				n9count[m]=0;n9r[m]=0;n9count[m+1]=0;n9r[m+1]=0;
				m=m+2;
			}
		}
		//End odd number if statement.




		if(i%2 == 0){
			for(m=0;m<POPSIZE;m++){
				temp1[m] = gen2->n[m];n2count[temp1[m]] = 1;
				if(gen2->np[m]<rec){
					if(gen2->n[m]%2==0){n2r[gen2->n[m]+1] = 1;}
					if(gen2->n[m]%2==1){n2r[gen2->n[m]-1] = 1;}


				}current[m] = m;siz[m]=0;indic[m]=0;recount[m]=0;
			}

			for(m=0;m<POPSIZE;m++){
				temp2[m] = temp1[gen3->n[m]];n3count[temp2[m]] = 1;
				if(gen3->np[m]<rec){
					if(gen3->n[m]%2==0){n3r[temp1[gen3->n[m]+1]] = 1;}
					if(gen3->n[m]%2==1){n3r[temp1[gen3->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp1[m] = temp2[gen4->n[m]];n4count[temp1[m]] = 1;
				if(gen4->np[m]<rec){
					if(gen4->n[m]%2==0){n4r[temp2[gen4->n[m]+1]] = 1;}
					if(gen4->n[m]%2==1){n4r[temp2[gen4->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp2[m] = temp1[gen5->n[m]];n5count[temp2[m]] = 1;
				if(gen5->np[m]<rec){
					if(gen5->n[m]%2==0){n5r[temp1[gen5->n[m]+1]] = 1;}
					if(gen5->n[m]%2==1){n5r[temp1[gen5->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp1[m] = temp2[gen6->n[m]];n6count[temp1[m]] = 1;
				if(gen6->np[m]<rec){
					if(gen6->n[m]%2==0){n6r[temp2[gen6->n[m]+1]] = 1;}
					if(gen6->n[m]%2==1){n6r[temp2[gen6->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp2[m] = temp1[gen7->n[m]];n7count[temp2[m]] = 1;
				if(gen7->np[m]<rec){
					if(gen7->n[m]%2==0){n7r[temp1[gen7->n[m]+1]] = 1;}
					if(gen7->n[m]%2==1){n7r[temp1[gen7->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp1[m] = temp2[gen8->n[m]];n8count[temp1[m]] = 1;
				if(gen8->np[m]<rec){
					if(gen8->n[m]%2==0){n8r[temp2[gen8->n[m]+1]] = 1;}
					if(gen8->n[m]%2==1){n8r[temp2[gen8->n[m]-1]] = 1;}


				}
			}

			for(m=0;m<POPSIZE;m++){
				temp2[m] = temp1[gen9->n[m]];n9count[temp2[m]] = 1;
				if(gen9->np[m]<rec){
					if(gen9->n[m]%2==0){n9r[temp1[gen9->n[m]+1]] = 1;}
					if(gen9->n[m]%2==1){n9r[temp1[gen9->n[m]-1]] = 1;}


				}
			}



			m=0;
			while(m<POPSIZE){
				flag1=0;
				while(ZERO==0){
					if(n2count[m]==0){flag1=1;break;}
					if(n3r[m]==0){
						if(n3count[m]==0){flag1=1;break;}
						if(n4r[m]==0){
							if(n4count[m]==0){flag1=1;break;}
							if(n5r[m]==0){
								if(n5count[m]==0){flag1=1;break;}
								if(n6r[m]==0){
									if(n6count[m]==0){flag1=1;break;}
									if(n7r[m]==0){
										if(n7count[m]==0){flag1=1;break;}
										if(n8r[m]==0){
											if(n8count[m]==0){flag1=1;break;}
											if(n9r[m]==0){
												if(n9count[m]==0){flag1=1;break;}


											}
										}
									}
								}
							}
						}
					}
					break;
				}


				flag2=0;
				while(ZERO==0){
					if(n2count[m+1]==0){flag2=1;break;}
					if(n3r[m+1]==0){
						if(n3count[m+1]==0){flag2=1;break;}
						if(n4r[m+1]==0){
							if(n4count[m+1]==0){flag2=1;break;}
							if(n5r[m+1]==0){
								if(n5count[m+1]==0){flag2=1;break;}
								if(n6r[m+1]==0){
									if(n6count[m+1]==0){flag2=1;break;}
									if(n7r[m+1]==0){
										if(n7count[m+1]==0){flag2=1;break;}
										if(n8r[m+1]==0){
											if(n8count[m+1]==0){flag2=1;break;}
											if(n9r[m+1]==0){
												if(n9count[m+1]==0){flag2=1;break;}


											}
										}
									}
								}
							}
						}
					}
					break;
				}

				if(flag1==1 && flag2==1){indic[m]=1;indic[m+1]=1;}
				else{
					if(flag1==1 && n2r[m]==0){indic[m]=1;}
					if(flag2==1 && n2r[m+1]==0){indic[m+1]=1;}


				}


				if(indic[m]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
					p = gen1->np[m];j = gen1->n[m];if(p < rec){if(j%2==0){k = j+1;}
						if(j%2==1){k = j-1;}
						++recount[(j+k)/2];
					}
					else{
						if(current[j]!=j){memcpy(&even[m].chr->q[0],&even[current[j]-POPSIZE].chr->q[0],4*size[j]);siz[m]=size[j];}
						if(current[j]==j){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+POPSIZE;siz[m] = size[j];}


					}
				}

				if(indic[m+1]==0 || (i%DELETE>=DELETE-10 || i%DELETE==0 || i>=GEN-10)){
					p = gen1->np[m+1];j = gen1->n[m+1];if(p < rec){if(j%2==0){k = j+1;}
						if(j%2==1){k = j-1;}
						++recount[(j+k)/2];
					}
					else{
						if(current[j]!=j){memcpy(&even[m+1].chr->q[0],&even[current[j]-POPSIZE].chr->q[0],4*size[j]);siz[m+1]=size[j];}
						if(current[j]==j){temp = even[m+1].chr;even[m+1].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+1+POPSIZE;siz[m+1] = size[j];}


					}
				}
				m=m+2;
			}



			for(m=0;m<POPSIZE;m++){

				p = gen1->np[m];j = gen1->n[m];

				if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m]==1){continue;}

				}

				if(p<rec){
					r = int((SEQLENGTH-1)*drand());if(j%2==0){k = j+1;}
					if(j%2==1){k = j-1;}
					--recount[(j+k)/2];
					if(current[j]==j && size[j]>0){
						if(r>=odd[j].chr->q[size[j]-1]){
							if(recount[(j+k)/2]==0){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j] = m+POPSIZE;siz[m] = size[j];}
							else{memcpy(&even[m].chr->q[0],&odd[j].chr->q[0],4*size[j]);siz[m]=size[j];}


						}
						else{
							start=0;end=size[j]-1;trk=end;ind=0;
							while(trk > 0 && trk <=size[j]-1){
								if(r >=odd[j].chr->q[trk-1] && r <odd[j].chr->q[trk]){ind=1;break;}
								if(r >=odd[j].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(r <odd[j].chr->q[trk]){end=trk;trk=(start+end)/2;}


							}
							if(ind==1){
								if(recount[(j+k)/2]==0){temp = even[m].chr;even[m].chr = odd[j].chr;odd[j].chr = temp;current[j]=m+POPSIZE;siz[m] = trk;}
								else{memcpy(&even[m].chr->q[0],&odd[j].chr->q[0],4*trk);siz[m]=trk;}


							}
						}
					}


					if(current[j]!=j && size[j]>0){
						trp=current[j]-POPSIZE;
						if(r>=even[trp].chr->q[size[j]-1]){memcpy(&even[m].chr->q[0],&even[trp].chr->q[0],4*size[j]);siz[m]=size[j];}
						else{
							start=0;end=size[j]-1;trk=end;ind=0;
							while(trk > 0 && trk<=size[j]-1){
								if(r >=even[trp].chr->q[trk-1] && r <even[trp].chr->q[trk]){ind=1;break;}
								if(r >=even[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(r <even[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}


							}
							if(ind==1){memcpy(&even[m].chr->q[0],&even[trp].chr->q[0],4*trk);siz[m]=trk;}


						}
					}

					if(current[k]==k && size[k]>0){
						if(r<odd[k].chr->q[0]){memcpy(&even[m].chr->q[siz[m]],&odd[k].chr->q[0],4*size[k]);siz[m]=siz[m]+size[k];}
						else{
							start=0;end=size[k]-1;trk=end;ind=0;
							while(trk > 0 && trk <=size[k]-1){
								if(r >=odd[k].chr->q[trk-1] && r <odd[k].chr->q[trk]){ind=1;break;}
								if(r >=odd[k].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(r <odd[k].chr->q[trk]){end=trk;trk=(start+end)/2;}


							}
							if(ind==1){memcpy(&even[m].chr->q[siz[m]],&odd[k].chr->q[trk],4*(size[k]-trk));siz[m]=siz[m]+size[k]-trk;}


						}
					}

					if(current[k]!=k && size[k]>0){trp=current[k]-POPSIZE;
						if(r<even[trp].chr->q[0]){memcpy(&even[m].chr->q[siz[m]],&even[trp].chr->q[0],4*size[k]);siz[m]=siz[m]+size[k];}
						else{
							start=0;end=size[k]-1;trk=end;ind=0;
							while(trk > 0 && trk<=size[k]-1){
								if(r >=even[trp].chr->q[trk-1] && r <even[trp].chr->q[trk]){ind=1;break;}
								if(r >=even[trp].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(r <even[trp].chr->q[trk]){end=trk;trk=(start+end)/2;}


							}
							if(ind==1){memcpy(&even[m].chr->q[siz[m]],&even[trp].chr->q[trk],4*(size[k]-trk));siz[m]=siz[m]+size[k]-trk;}


						}
					}
				}
			}




			for(m=0;m<POPSIZE;m++){

				if(i%DELETE < DELETE-10 && i%DELETE > 0 && i < GEN-10){if(indic[m]==1){continue;}

				}

				//INSERT NEW MUTATIONS
				k = poisson(mu);
				for(l=0;l<k;l++){
					trp = int(SEQLENGTH*drand());while(multhit[trp]==1){trp = int(SEQLENGTH*drand());}
					multhit[trp]=1;
					if((siz[m]>0 && trp>=even[m].chr->q[siz[m]-1]) || siz[m]==0){even[m].chr->q[siz[m]]=trp;}
					else{
						if(trp<=even[m].chr->q[0]){memmove(&even[m].chr->q[1],&even[m].chr->q[0],4*siz[m]);even[m].chr->q[0]=trp;}
						else{
							start=0;end=siz[m]-1;trk=end;ind=0;
							while(trk > 0 && trk<=siz[m]-1){
								if(trp > even[m].chr->q[trk-1] && trp <=even[m].chr->q[trk]){ind=1;break;}
								if(trp > even[m].chr->q[trk]){start=trk;trk=(start+end)/2;if(trk==start){++trk;}
									continue;
								}
								if(trp <=even[m].chr->q[trk]){end=trk;trk=(start+end)/2;continue;}


							}
							if(ind==1){memmove(&even[m].chr->q[trk+1],&even[m].chr->q[trk],4*(siz[m]-trk));even[m].chr->q[trk]=trp;}


						}
					}
					++siz[m];
				}
			}


			gtemp = gen1;gen1 = gen2;gen2 = gen3;gen3 = gen4;gen4 = gen5;gen5 = gen6;gen6 = gen7;gen7 = gen8;gen8=gen9;gen9 = gtemp;
			m=0;
			while(m<POPSIZE){
				gen9->np[m] = drand();gen9->np[m+1] = drand();x=drand();
				if(x>=s){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]=int(POPSIZE*drand());}
				else{y=drand();if(y<=0.50){gen9->n[m]=int(POPSIZE*drand());gen9->n[m+1]= gen9->n[m];}
					else{gen9->n[m]=int(POPSIZE*drand());if((gen9->n[m])%2==0){gen9->n[m+1]=(gen9->n[m])+1;}
						else{gen9->n[m+1]=(gen9->n[m])-1;}


					}
				}

				n2count[m]=0;n2r[m]=0;n2count[m+1]=0;n2r[m+1]=0;
				n3count[m]=0;n3r[m]=0;n3count[m+1]=0;n3r[m+1]=0;
				n4count[m]=0;n4r[m]=0;n4count[m+1]=0;n4r[m+1]=0;
				n5count[m]=0;n5r[m]=0;n5count[m+1]=0;n5r[m+1]=0;
				n6count[m]=0;n6r[m]=0;n6count[m+1]=0;n6r[m+1]=0;
				n7count[m]=0;n7r[m]=0;n7count[m+1]=0;n7r[m+1]=0;
				n8count[m]=0;n8r[m]=0;n8count[m+1]=0;n8r[m+1]=0;
				n9count[m]=0;n9r[m]=0;n9count[m+1]=0;n9r[m+1]=0;
				m=m+2;
			}
		}


		//REMOVE FIXED MUTATIONS
		if(i%DELETE==0 || i==GEN){
			for(l=0;l<SEQLENGTH;l++){indicator[l]=0;}
			for(m=0;m<POPSIZE;m++){
				if(siz[m] > 2500 + int(GEN*mu)){cout<<"Out of Memory. Increase array size.\n";return(0);}
				for(l=0;l<siz[m];l++){trk=even[m].chr->q[l];++indicator[trk];}


			}

			for(l=0;l<SEQLENGTH;l++){if(indicator[l]==0 || indicator[l]==POPSIZE){multhit[l]=0;}

			}

			int count=0;
			for(m=0;m<POPSIZE;m++){
				k=0;trp=siz[m];l=0;
				while(l< trp){
					trk = even[m].chr->q[l];
					if(indicator[trk]!=POPSIZE){even[m].chr->q[k] = trk;++k;++l;continue;}
					if(indicator[trk]==POPSIZE){++l;continue;}


				}siz[m]=k;count = count + siz[m];
			}
			cout<<i<<" "<<count<<"  "<< s <<"\n";
		}
	}
	//End main for-loop!

	ofstream output1("finalpopulation.txt",ios::out);
	for(m=0;m<SSIZE;m++){
		for(l=0;l<siz[m];l++){output1<<even[m].chr->q[l]<<" ";}
		output1<<"\n";
	}
	output1.close();







}














