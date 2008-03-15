#include <cmath>
#include "gtype_cel_to_pq.h"
#include "log_average.h"

void log_average(float (*intensity)[N_MATCH_TYPES][MAX_PQ], int (*count)[N_MATCH_TYPES], double *A_il, double *B_il) {
  double vA_il = 0.0, vB_il = 0.0;
  /* 
     for(int a=0; a < N_ALLELES; a++) {
     for(int m=0; m < N_MATCH_TYPES; m++) {
  */
  for(int pq=0; pq < count[0][0]; pq++) {
    double M_il = (log(intensity[0][1][pq]) + log(intensity[1][1][pq]))/2.0;
    if ( intensity[0][0][pq] > intensity[0][1][pq] ) {
      vA_il += (log(intensity[0][0][pq]) -  M_il);
    } /* else add nothing */
    if ( intensity[1][0][pq] > intensity[1][1][pq] ) {
      vB_il += (log(intensity[1][0][pq]) - M_il);
    } /* else add nothing */
  }
  
  *A_il = vA_il / count[0][0];
  *B_il = vB_il / count[0][0]; 
  return;
}
