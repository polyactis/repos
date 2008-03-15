
/*
void log_average(float intensity[N_ALLELES][N_MATCH_TYPES][MAX_PQ], int count[N_ALLELES][N_MATCH_TYPES], double *A_il, double *B_il);

void log_average(float ***intensity, int **count, double *A_il, double *B_il);
*/

void log_average(float (*intensity)[N_MATCH_TYPES][MAX_PQ], int (*count)[N_MATCH_TYPES], double *A_il, double *B_il);
