#include <cmath>
#include "gtype_cel_to_pq.h"
#include "log_average.h"

void log_average(float (*intensity)[N_MATCH_TYPES][MAX_PQ], int (*count)[N_MATCH_TYPES], double *A_il, double *B_il)
{
	double vA_il = 0.0, vB_il = 0.0;
	/*
	   for(int a=0; a < N_ALLELES; a++) {
	   for(int m=0; m < N_MATCH_TYPES; m++) {
	*/
	/* 03/16/08 yh: this block assumes presence of MM, which is not available in 250k.
	for(int pq=0; pq < count[0][0]; pq++)	//03/16/08 yh: definition of intensity: float intensity[N_ALLELES][N_MATCH_TYPES][MAX_PQ];
	{
		double M_il = (log(intensity[0][1][pq]) + log(intensity[1][1][pq]))/2.0;	//03/16/08 yh: average of log(MM of allele A&B) as background. 
		if ( intensity[0][0][pq] > intensity[0][1][pq] )	//03/16/08 yh: if PM>MM for allele A. this pair is good. add it.
		{
			vA_il += (log(intensity[0][0][pq]) -  M_il);
		} // else add nothing
		if ( intensity[1][0][pq] > intensity[1][1][pq] )	//03/16/08 yh: ditto for allele B
		{
			vB_il += (log(intensity[1][0][pq]) - M_il);
		} // else add nothing
	}
	*/
	vA_il += log(intensity[0][0][0]) + log(intensity[0][0][1]);	//03/16/08 yh: sum of log(antisense)+ log(sense) of allele A 
	vB_il += log(intensity[1][0][0]) + log(intensity[1][0][1]);	//03/16/08 yh: ditto for allele B
	*A_il = vA_il / count[0][0];
	*B_il = vB_il / count[0][0];
	return;
}
