/**
 * Copyright (c) 2005-2006 Hin-Tak Leung, modifications
 *
 * **** Be very careful that while snps and samples are both within 32-bit,
 * their product can go over. ****
 */

////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
/////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include "hash.h"
#include "CELFileData.h"
#include "CDFFileData.h"

#include "stable_compare.h"

using namespace std;

#include "gtype_cel_to_pq.h"
#include "log_average.h"

enum MATCH_TYPES {MATCH_PM, MATCH_MM};

struct opt {
  int help;
  int debug;
  int subset;
  int no_header;
  int n_random;
  std::string cdfFileName;
  std::string outDirName;
  std::string subsetFileName;
  std::string outFileSuffix;
  bool with_stdv;
  bool log_average;
  int n_cel;
  char *celFile[MAX_CEL];
  FILE *intfile;
  FILE *single;
  FILE *split;
};

#include "intensity_rank.h"

int process_args(int argc, char *argv[], struct opt *o);
void usage(char *argv[]);
std::string getOutFileName(std::string celName, struct opt *o);
std::string getStem(std::string celName);
int getMatchType (char x, char y);
int base2int (char b);

int main(int argc, char *argv[]) {
  /* process options and arguments */
  struct opt o;
  double *refintensity = NULL;
  int *reforder = NULL;
  t_intensity_rank *ranks = NULL;
  if(process_args(argc,argv,&o)==EXIT_FAILURE) {
    cerr << argv[0] << ": problem processing command-line args\n";
    exit(EXIT_FAILURE);
  }

  if (o.intfile) {
    refintensity = (double *)calloc(CEL_INTENSITY_ARRAY_SIZE, sizeof(double));
    size_t readsize = fread(refintensity, sizeof(double), CEL_INTENSITY_ARRAY_SIZE, o.intfile);
    cout << "read " << readsize << " of " << CEL_INTENSITY_ARRAY_SIZE << endl;
    if (readsize != CEL_INTENSITY_ARRAY_SIZE) {
      cout << "read error, aborting " << endl;
      exit(EXIT_FAILURE);
    }
    ranks = (t_intensity_rank *)calloc(CEL_INTENSITY_ARRAY_SIZE, sizeof(t_intensity_rank));
    reforder = (int *)calloc(CEL_INTENSITY_ARRAY_SIZE, sizeof(int));
  }

  // If subset file was specified, read it and install probeset names in hash
  struct nlist *subsetHash[HASH_SIZE];
  double *Alog = NULL, *Blog = NULL;
  /* snp_idx is the order of in which a specific snp appear,
     which can be different from lineNo if
     there is bad entries */
  int snp_idx = 0;
  std::vector<std::string> snpname(MAX_SNPS);
  if(o.subset) {
    char line[MAX_LINE_LEN];
    char snpID[MAX_LINE_LEN];
    hashtab_init(subsetHash,HASH_SIZE);
    ifstream subsetFile(o.subsetFileName.c_str());
    if(!subsetFile.is_open()) {
      cerr << argv[0] << ": problem opening file " << o.subsetFileName << "\n";
      exit(EXIT_FAILURE);
    }
    int lineNo=0;
    while(!subsetFile.eof()) {
      subsetFile.getline(line,MAX_LINE_LEN-1);
      lineNo++;
      if(strlen(line) == 0)
        continue;
      if(1 != sscanf(line,"%s",snpID)) {
        cerr << argv[0] << ": skipping bad entry in line " << lineNo << " of file " << o.subsetFileName << "\n";
      } else if(NULL==hashtab_lookup(snpID,subsetHash)) {
	snpname[snp_idx] = snpID;
        if(NULL==hashtab_install(snpID,snp_idx++,subsetHash)) {
          cerr << argv[0] << ": problem inserting value " << snpID << " into hash table\n";
          exit(EXIT_FAILURE);
        } else if(o.debug) {
          cout << "inserted value " << snpID << " into hash table\n";
        }
      }
    }
    subsetFile.close();
    snpname.resize(snp_idx);
    if ((o.single) || (o.split)){
      Alog = (double *) malloc((size_t) sizeof(double) * (size_t) snp_idx * (size_t) o.n_cel); /* the product can go over 32-bit */
      Blog = (double *) malloc((size_t) sizeof(double) * (size_t) snp_idx * (size_t) o.n_cel); /* the product can go over 32-bit */
      if ((!Alog) || (!Blog)) {
	cerr << "Can't allocate enough memory\n";
	exit(EXIT_FAILURE);
      }
    }
  }

  // read CDF file
  affxcdf::CCDFFileData cdf;
  cdf.SetFileName(o.cdfFileName.c_str());
  if(o.debug)
    cout << "reading cdf file " << o.cdfFileName << "\n";
  if(!(cdf.Exists())) {
    cerr << "ERROR: CDF file " << o.cdfFileName << " not found\n";
    exit(EXIT_FAILURE);
  }
  cdf.Read();
  affxcdf::CCDFFileHeader cdfHeader = cdf.GetHeader();

  int nPsets = cdfHeader.GetNumProbeSets();
  // Iterate over all probesets and if its a genotyping probeset (and is in the
  // specified subset) put it in the vector of elegible probesets.
  std::vector<int> gtype_psets(nPsets);
  int gtype_psets_index=0;
  for(int probeset_index=0; probeset_index < nPsets; probeset_index++) {
    affxcdf::CCDFProbeSetInformation probeSet;
    cdf.GetProbeSetInformation(probeset_index, probeSet);
    // check if it is a genotyping probeset
    if(probeSet.GetProbeSetType() != affxcdf::GenotypingProbeSetType)
      continue;
    int count[N_ALLELES][N_MATCH_TYPES];
    for(int a=0; a < N_ALLELES; a++)
      for(int m=0; m < N_MATCH_TYPES; m++)
        count[a][m] = 0;
    std::string SNP_ID = cdf.GetProbeSetName(probeset_index);
    char SNP_ID_cstr[MAX_LINE_LEN];
    strcpy(SNP_ID_cstr,SNP_ID.c_str());
    if((o.subset) && (NULL==hashtab_lookup(SNP_ID_cstr,subsetHash)))
      continue;
    affxcdf::CCDFProbeGroupInformation probesetBlock;
    int nGroups = probeSet.GetNumGroups();
    if(o.debug)
      cout << "On probeset " << SNP_ID << " index " << probeset_index << " (" << nGroups <<" groups)\n";
    // Loop over all each block in the probeset to count the number of probes for each allele and of each type (PM/MM)
    for(int group_index=0; group_index < nGroups; group_index++) {
      probeSet.GetGroupInformation(group_index, probesetBlock);
      int nCell = probesetBlock.GetNumCells();
      if(nCell >= (MAX_PQ/nGroups)) {
        cerr << argv[0] << ": " << nCell << " probes in "<< nGroups <<
	  " groups in probeset - too many (recompile with larger MAX_PQ than currently " << MAX_PQ << "\n";
        exit(EXIT_FAILURE);
      }
      if(0 != (nCell % N_MATCH_TYPES)) {
        cerr << argv[0] << ": problem reading probeset " << SNP_ID << " group " << group_index << ": number of probes expected to be a multiple of 2\n";
        continue;
      }
      int allele = group_index % N_ALLELES;
      affxcdf::CCDFProbeInformation feature;
      for(int probe_index=0; probe_index < nCell; probe_index++) {
        probesetBlock.GetCell(probe_index, feature);
        int matchType = getMatchType(feature.GetPBase(),feature.GetTBase());
        count[allele][matchType] += 1;
      }
    }
    bool is_in_quartet_structure = true;
    for(int a=0; a < N_ALLELES; a++)
      for(int m=0; m < N_MATCH_TYPES; m++)
        if(count[a][m] != count[0][0])
          is_in_quartet_structure = false;
    if(is_in_quartet_structure)
      gtype_psets[gtype_psets_index++] = probeset_index;
  }

  if ((o.subset) && (gtype_psets_index != snp_idx)) {
    cerr << "Horror - snps in the subset file = " << snp_idx <<
      ", only " << gtype_psets_index << " found in the cdf file." << endl;
    exit(EXIT_FAILURE);
  }
  gtype_psets.resize(gtype_psets_index);

  // Valid probesets have been identified, now pick the ones to work on and put them in the vector probeset_indices
  std::vector <int> probeset_indices;
  int n_available = gtype_psets.size();
  if(o.debug)
    cout << "Found " << n_available << " eligible psets\n";
  if(o.n_random > 0) {
    if(o.n_random > n_available) {
      cerr << "ERROR: request was for " << o.n_random << " random probesets but only " << n_available << " available\n";
      exit(EXIT_FAILURE);
    }
    cerr << "ERROR: -random option not currently implemented\n";
    exit(EXIT_FAILURE);
    //probeset_indices.resize(o.n_random);
    //for(int i=0; i < o.n_random; i++)
    //  probeset_indices[i] = gtype_psets[i];
  } else {
    probeset_indices.resize(n_available);
    std::vector<int>::iterator pset_it;
    int i;
    for(i=0, pset_it=gtype_psets.begin(); pset_it != gtype_psets.end(); pset_it++, i++)
      probeset_indices[i] = *pset_it;
  }

  if(o.debug) {
    cout << "Working on " << probeset_indices.size() << " probesets\n";
    for(unsigned int i=0; i<probeset_indices.size(); i++)
      cout << "  " << probeset_indices[i] << "\n";
  }

  // read each CEL file and open a file for output intensities
  for(int cel_index=0; cel_index < o.n_cel; cel_index++) {
    /* locate and read CEL file */
    affxcel::CCELFileData cel;
    std::string celFileName;

    celFileName = o.celFile[cel_index];
    cel.SetFileName(celFileName.c_str());
    if(!cel.Exists()) {
      cerr << "ERROR: CEL file " << cel.GetFileName() << " not found\n";
      exit(EXIT_FAILURE);
    }
    cout << "Processing " << cel.GetFileName() << endl;
    cel.Read();

    if (o.intfile) {
      for (int icel=0; icel< CEL_INTENSITY_ARRAY_SIZE; icel++) {
	ranks[icel].intensity = cel.GetIntensity(icel);
	ranks[icel].rank = icel;
      }
      qsort(ranks, CEL_INTENSITY_ARRAY_SIZE, sizeof(t_intensity_rank), stable_compare_intensity_rank);
      for (int icel=0; icel< CEL_INTENSITY_ARRAY_SIZE; icel++) {
	reforder[ranks[icel].rank] = icel;
      }
    }
    // Get name of output file and open it
    ofstream outFile;
    if (!(o.single) && !(o.split)) {
      std::string outFileName = getOutFileName(celFileName,&o);
      if(o.debug)
	cout << "opening file " << outFileName << "\t(" << (1+cel_index) << " of " << o.n_cel << ")\n";
      outFile.open(outFileName.c_str());
    }
    if(cel.GetChipType() != cdf.GetChipType()) {
      cerr << "ERROR: mismatch between CHP file type and CDF file type for " << celFileName << "\n";
      exit(EXIT_FAILURE);
    }

    // Loop over all probeset and write out intensities
    float intensity[N_ALLELES][N_MATCH_TYPES][MAX_PQ];
    float stdv[N_ALLELES][N_MATCH_TYPES][MAX_PQ];
    int count[N_ALLELES][N_MATCH_TYPES];
    for(int i=0; i < (int) probeset_indices.size(); i++) {
      if(o.debug && i==10)
        break;
      // We assume that we are iterating over only genotyping probesets as that was checked earlier
      int probeset_index = probeset_indices[i];
      affxcdf::CCDFProbeSetInformation probeSet;
      cdf.GetProbeSetInformation(probeset_index, probeSet);
      for(int a=0; a < N_ALLELES; a++)
        for(int m=0; m < N_MATCH_TYPES; m++)
          count[a][m] = 0;
      std::string SNP_ID = cdf.GetProbeSetName(probeset_index);
      char SNP_ID_cstr[MAX_LINE_LEN];
      strcpy(SNP_ID_cstr,SNP_ID.c_str());
      affxcdf::CCDFProbeGroupInformation probesetBlock;
      int nGroups = probeSet.GetNumGroups();
      if(o.debug)
        cout << "On probeset " << SNP_ID << " (" << nGroups <<" groups)\n";
      for(int group_index=0; group_index < nGroups; group_index++) {
        probeSet.GetGroupInformation(group_index, probesetBlock);
        int nCell = probesetBlock.GetNumCells();
        int allele = group_index % N_ALLELES;
        affxcdf::CCDFProbeInformation feature;
        for(int probe_index=0; probe_index < nCell; probe_index++) {
          probesetBlock.GetCell(probe_index, feature);
          int matchType = getMatchType(feature.GetPBase(),feature.GetTBase());
          int x = feature.GetX();
          int y = feature.GetY();

	  if (o.intfile) {
	    intensity[allele][matchType][count[allele][matchType]] = refintensity[reforder[cel.XYToIndex(x,y)]];
	  }
	  else {
	    intensity[allele][matchType][count[allele][matchType]] = cel.GetIntensity(x,y);
	  }
	
          stdv[allele][matchType][count[allele][matchType]] = cel.GetStdv(x,y);
          count[allele][matchType] += 1;
        }
      }

      if (o.single || o.split) {
	double A_il = 0.0, B_il = 0.0;
	unsigned long idx = (hashtab_lookup(SNP_ID_cstr,subsetHash))->value;	
	log_average(intensity, count, &A_il, &B_il);
	Alog[idx * (size_t) o.n_cel + cel_index] = A_il; /* the product can go over 32-bit */
	Blog[idx * (size_t) o.n_cel + cel_index] = B_il; /* the product can go over 32-bit */
      } else {
	outFile << SNP_ID;
	if(o.log_average) {
	  double A_il = 0.0, B_il = 0.0;
	  log_average(intensity, count, &A_il, &B_il);
	  outFile << "\t" << A_il << "\t" << B_il ;
	} else {
	  for(int pq=0; pq < count[0][0]; pq++) {
	    for(int a=0; a < N_ALLELES; a++) {
	      for(int m=0; m < N_MATCH_TYPES; m++) {
		outFile << "\t" << intensity[a][m][pq];
		if (o.with_stdv) {
		  outFile << "\t" << stdv[a][m][pq];
		}
	      }
	    }
	  }
	}
	outFile << "\n";
      }
    }

    cel.Close();
  }

  // cleanup
  cdf.Close();
  if(o.subset) {
    if (o.single) {
      for (int i_cel = 0; i_cel < o.n_cel; i_cel ++) {
	std::string stem = getStem(o.celFile[i_cel]);
	fprintf(o.single, "\t%s%s\t%s%s", stem.c_str(), "_A", stem.c_str(), "_B");
      }
      fprintf(o.single, "\n");
      for (int i_snp = 0; i_snp < snp_idx; i_snp ++) {
	fprintf(o.single, "%s", snpname[i_snp].c_str());
	for (int i_cel = 0; i_cel < o.n_cel; i_cel ++) {
	  fprintf(o.single, "\t%f\t%f",
		  Alog[i_snp * (size_t) o.n_cel + i_cel],  /* the product can go over 32-bit */
		  Blog[i_snp * (size_t) o.n_cel + i_cel]); /* the product can go over 32-bit */
	}
	fprintf(o.single, "\n");
      }
      fclose(o.single);
    }
    
    if (o.split) {
      char stem[PATH_MAX];
      char ext[PATH_MAX];
      char currentfilename[PATH_MAX];
      FILE *current_fp = NULL;
      int count_snp = 0;
      int idx_snp = 0;
      int scanned = 0;
      
      fscanf(o.split, "%s", stem);
      
      scanned = fscanf(o.split, "%d %s", &count_snp, ext);
      sprintf(currentfilename, "%s_%s.txt", stem, ext);
      current_fp = fopen(currentfilename, "wb");
      
      for (int i_cel = 0; i_cel < o.n_cel; i_cel ++) {
	std::string stem = getStem(o.celFile[i_cel]);
	fprintf(current_fp, "\t%s%s\t%s%s", stem.c_str(), "_A", stem.c_str(), "_B");
      }
      fprintf(current_fp, "\n");
      
      for (int i_snp = 0; i_snp < snp_idx; i_snp ++) {
	fprintf(current_fp, "%s", snpname[i_snp].c_str());
	for (int i_cel = 0; i_cel < o.n_cel; i_cel ++) {
	  fprintf(current_fp, "\t%f\t%f",
		  Alog[i_snp * (size_t) o.n_cel + i_cel],  /* the product can go over 32-bit */
		  Blog[i_snp * (size_t) o.n_cel + i_cel]); /* the product can go over 32-bit */
	}
	fprintf(current_fp, "\n");
	idx_snp++;
	if(idx_snp == count_snp) {
	  idx_snp = 0;
	  fclose(current_fp);
	  scanned = fscanf(o.split, "%d %s", &count_snp, ext);
	  if (scanned == 2) {
	    sprintf(currentfilename, "%s_%s.txt", stem, ext);
	    current_fp = fopen(currentfilename, "wb");
	    
	    for (int i_cel = 0; i_cel < o.n_cel; i_cel ++) {
	      std::string stem = getStem(o.celFile[i_cel]);
	      fprintf(current_fp, "\t%s%s\t%s%s", stem.c_str(), "_A", stem.c_str(), "_B");
	    }
	    fprintf(current_fp, "\n");
	  } else {
	    if (scanned == EOF) {
	      cerr << "EOF of split file reached" << endl;
	    } else if (scanned == 0 && feof(o.split)) {
	      cerr << "EOF of split file reached under mscrt/wine reached" << endl;	      
	    } else {
	      cerr << "Split file seems to be broken" << endl;
	    }
	  }
	}
      }      
      assert(idx_snp==0); // make sure things add up
      fclose(o.split);
    }
    
    hashtab_free(subsetHash);
    
    if (o.single || o.split) {
      free(Alog);
      free(Blog);
    }
  }

  if (o.intfile) {
    // house keeping, probably not necessary since we are exiting soon...
    fclose(o.intfile);
    free(refintensity);
    free(ranks);
    free(reforder);
  }

  exit(EXIT_SUCCESS);
}


/* Print usage info to STDOUT */
void usage(char *argv[]) {
  cout << "\n";
  cout << argv[0] << "\n";
  cout << "\n";
  cout << "  Extracts probe quartet intensities from genotyping CEL files.\n";
  cout << "\n";
  cout << "DESCRIPTION:\n";
  cout << "\n";
  cout << "  Probe intensities are written to a tab-delimited text file for each CEL file\n";
  cout << "  processed.  The text file with have the extension .raw and by default will be\n";
  cout << "  located in the same directory as the CEL file.  The format is tab-delimted text\n";
  cout << "  with one row per SNP.  The first field is the probeset identifier, followed by\n";
  cout << "  quartets of probe intensities.  Each quartet is in the form\n";
  cout << "  (PM_A,MM_A,PM_B,MM_B).  Only genotyping probesets in quartet structure will be\n";
  cout << "  reported.\n";
  cout << "\n";
  cout << "USAGE:\n  [options] -cdf cdfFile file1.cel ... fileN.cel\n";
  cout << "\n";
  cout << "OPTIONS:\n";
  cout << "  -h                Print usage information\n";
  cout << "  -cdf file.cdf     CDF file\n";
  cout << "  -outdir           output files will be written to outdir\n";
  cout << "  -single <file>    output files will be written to one single file (must be used together with -log-average and -subset)\n";
  cout << "  -split <file>     files will be written splitted as directed (must be used together with -log-average and -subset)\n";
  cout << "  -subset  <file>   text file with list of probeset IDs to extract\n";
  cout << "  -no_header        Supress header line in output\n";
  cout << "  -debug            Run in debug mode (extra debug info written to STDOUT)\n";
  cout << "  -n_random n       (CURRENTLY NOT IMPLEMENTED) Randomly sample n probesets from the CDF (intersected with subset if specified)\n";
  cout << "  -with-stdv        with stdv values next to the intensity values\n";
  cout << "  -log-average      output log average values instead of all the intensity of all the probe quartets\n";
  cout << "  -refintensity <file>  reference intensity file generated by cel_qnorm_pass1\n";
  cout << "\n";
}

/* Processes command-line options */
int process_args(int argc, char *argv[], struct opt *o) {
  int i;
  char *thisarg;

  /* set defaults */
  o->help = 0;
  o->debug = 0;
  o->subset = 0;
  o->no_header = 0;
  o->n_random = 0;
  o->cdfFileName = "";
  o->outDirName = "";
  o->subsetFileName = "";
  o->outFileSuffix = ".raw";
  o->n_cel = 0;
  o->intfile = NULL;
  o->single = NULL;
  o->split = NULL;
  o->with_stdv = false;
  o->log_average = false;

  /* read args */
  if(argc == 1)
    o->help = 1;
  for (i=1; i<argc; i++) {
    thisarg = argv[i];
    if(thisarg[0] == '-') {
      if (!(strcmp(thisarg+1,"h") && strcmp(thisarg+1,"help"))) {
        o->help = 1;
      } else if (!strcmp(thisarg+1,"debug")) {
        o->debug = 1;
      } else if (!strcmp(thisarg+1,"cdf")) {
        i++;
        if (i>=argc) {
          cerr << "must provide a filename with -cdf option\n";
          return(EXIT_FAILURE);
        } else {
          o->cdfFileName = argv[i];
        }
      } else if (!strcmp(thisarg+1,"subset")) {
        i++;
        if (i>=argc) {
          cerr << "must provide a subset file with -subset option\n";
          return(EXIT_FAILURE);
        } else {
          o->subsetFileName = argv[i];
          o->subset = 1;
        }
      } else if (!strcmp(thisarg+1,"n_random")) {
        i++;
        if (i>=argc) {
          cerr << "must provide an integer value with -n_random option\n";
          return(EXIT_FAILURE);
        } else {
          o->n_random = atoi(argv[i]);
        }
      } else if (!strcmp(thisarg+1,"outdir")) {
        i++;
        if (i>=argc) {
          cerr << "must provide an output directory with -outdir option\n";
          return(EXIT_FAILURE);
        } else {
          o->outDirName = argv[i];
        }
      } else if (!strcmp(thisarg+1,"with-stdv")) {
	o->with_stdv = true;
      } else if (!strcmp(thisarg+1,"log-average")) {
	o->log_average = true;
	o->outFileSuffix = ".logavg";
      } else if (!strcmp(thisarg+1,"refintensity")) {
        i++;
        if (i>=argc) {
          cerr << "must provide a reference file with -refintensity option\n";
          return(EXIT_FAILURE);
        } else {
          o->intfile = fopen(argv[i], "rb");
	  if (!(o->intfile)) {
	    cerr << "intfile does not exist" << endl;
	    return(EXIT_FAILURE);
	  }
        }
      } else if (!strcmp(thisarg+1,"single")) {
        i++;
        if (i>=argc) {
          cerr << "must provide an output file with -single option\n";
          return(EXIT_FAILURE);
        } else {
          o->single = fopen(argv[i], "wb");
        }
      } else if (!strcmp(thisarg+1,"split")) {
        i++;
        if (i>=argc) {
          cerr << "must provide an output file with -split option\n";
          return(EXIT_FAILURE);
        } else {
          o->split = fopen(argv[i], "rb");
	  if (!(o->split)) {
	    cerr << "split does not exist" << endl;
	    return(EXIT_FAILURE);
	  }
        }
      } else if (!strcmp(thisarg+1,"no_header")) {
        o->no_header = 1;
      } else {
        cerr << "ERROR: " << argv[0] << ": unrecognized option: " << thisarg << "\n";
        return(EXIT_FAILURE);
      }
    } else {
      (o->celFile)[o->n_cel] = thisarg;
      o->n_cel++;
    }
  }

  if(o->help) {
    usage(argv);
    exit(EXIT_SUCCESS);
  }

  if(o->single || o->split) {
    if ((!(o->log_average)) || (!o->subset)){
      cerr << argv[0] << ": must specify -log-average and -subset option with -single/-split (for help, use -h option)\n";
      exit(EXIT_FAILURE);
    }
  }

  if( (o->n_cel > 0) && (!std::strcmp(o->cdfFileName.c_str(),"")) ) {
    cerr << argv[0] << ": must specify cdf file with -cdf option (for help, use -h option)\n";
    exit(EXIT_FAILURE);
  }

  return(EXIT_SUCCESS);
}

/* Translates call codes in genotyping CHP file to AA/AB/BB/NN */
// Looks for .CEL or .cel suffix in celName and replaces it by the suffix as specified in the opt struct
// If expected suffix not found, just appends .txt
std::string getOutFileName(std::string celName, struct opt *o) {
  int n;
  std::string outName;

  int len = celName.length();
  if(strcmp(o->outDirName.c_str(),"")) {
    if((n = celName.rfind('/')) > 0)
      celName = celName.substr(n);
    celName = o->outDirName + celName;
  }

  if(((n = celName.rfind(".CEL")) > 0) && (n == (len-4))) {
    outName = celName.substr(0,n) + o->outFileSuffix;
  } else if(((n = celName.rfind(".cel")) > 0) && (n == (len-4))) {
    outName = celName.substr(0,n) + o->outFileSuffix;
  } else {
    outName = celName + o->outFileSuffix;
  }

  return(outName);
}

// remove directory path and .CEL or .cel suffix,and returning only 
// the sample ID.
std::string getStem(std::string celName) {
  int n;
  std::string outName;

  if((n = celName.rfind('/')) > 0)
    celName = celName.substr(n + 1); // skip the "/" itself

  if( ((n = celName.rfind("_NSP")) > 0) || 
      ((n = celName.rfind("_STY")) > 0) ||
      ((n = celName.rfind("_Nsp")) > 0) || 
      ((n = celName.rfind("_Sty")) > 0) ||
      ((n = celName.rfind("_nsp")) > 0) || 
      ((n = celName.rfind("_sty")) > 0) ) {
    outName = celName.substr(0,n);
  } else if ( ((n = celName.rfind(".cel")) > 0) || 
	      ((n = celName.rfind(".CEL")) > 0) ) {
    outName = celName.substr(0,n);
  } else {
    outName = celName;
  }
  return(outName);
}

int getMatchType (char x, char y) {
  int matchType = ((base2int(x) + base2int(y))==3) ? MATCH_PM : MATCH_MM;
  return(matchType);
}

int base2int (char b) {

  switch(b) {

  case 'a':
  case 'A':
    return(0);
    break;

  case 'c':
  case 'C':
    return(1);
    break;

  case 'g':
  case 'G':
    return(2);
    break;

  case 't':
  case 'T':
    return(3);
    break;

  default:
    cerr << "Base base " << b << "\n";
    exit(EXIT_FAILURE);
  }
}
