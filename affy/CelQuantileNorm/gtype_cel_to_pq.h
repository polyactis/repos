/*
  Extracted from gtype_cel_to_pq.cpp, because some of the defines 
  are needed in log_average.c. 
  
  Values slightly modified. 

  Original value of MAX_CEL       10000
  Original value of MAX_PQ        10000

  Reason of change: they are slightly wasteful.
  
  MAX_CEL:
  The maximum number of files one can read currently depends on max-number 
  of characters one can fit on a command line, which is 65536, I think.
  Only 6.5 characters per file name is a bit unrealistic. Setting it too large
  is dangerous - better the program crashes during the initial stage
  of enumerating the cel files, than later trying to open non-existent 
  file (after processing all the previous).

  MAX_PQ:
  The code requires that nCell < MAX_PQ/nGroups , where nCell is typically
  24 or 40, and nGroups is what? Anyway, the warning 
  "too many groups in probeset - recompile with larger MAX_PQ"
  is retained has been updated to be more verbose.
*/

#define MAX_CEL       2000      // Max number of cel files which can be processed in one call
#define MAX_PQ        200       // Max number of probe quartets in a probeset
#define MAX_LINE_LEN  4096      // Max length of line in subset file
#define N_ALLELES     2         // Because we deal with bi-allelic SNPs
#define N_MATCH_TYPES 2         // PM and MM

#define MAX_SNPS      270000    // one of NSP/STY pair in the 500k chip
// MAX_CEL * MAX_SNPS *sizeof(double) * 2 is the max allocation.

#define CEL_INTENSITY_ARRAY_SIZE 6553600

// should be defined in linux/limits.h, just in case.
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif
