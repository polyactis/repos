// BiCluster.h: interface for the BiCluster class.
//
#include "datastructure.h"

typedef struct HAPLOTYPE
{
    int repID;
    int hCn;
    double bFreq;
    int code;
} HAPLOTYPE;

typedef struct PARAMETER
{
    vector<HAPLOTYPE> haps;
    vector<int> hapId;

    vector<HAPLOTYPE> hapsU;
    vector<int> hapIdU;
    bool Dsp;

    vector<vector<int> > nlCn;
    vector<vector<int> > nlSCn;
    vector<int> Dloci;
    vector<int> DSloci;
    vector<int> lociStatus;
    vector<bool> indState;
    int Dcn;
    double logP;
    double logpart[3];

    bool operator=(PARAMETER &a)
    {
        Dcn = a.Dcn;
        Dloci = a.Dloci;
        DSloci = a.DSloci;
        Dsp = a.Dsp;
        hapId = a.hapId;
        hapIdU = a.hapIdU;
        haps = a.haps;
        hapsU = a.hapsU;
        indState = a.indState;
        logP = a.logP;
        nlCn = a.nlCn;
        nlSCn = a.nlSCn;
        return true;
    }

} PARAMETER;

typedef struct OUTPUT
{
    vector<int> markers;
    double pvalue;
} OUTPUT;

typedef struct MARKERSET
{
    vector<int> markers;
    int code;
    int count;
} MARKERSET;

class BiCluster
{

public:
    BiCluster(int acn);
    virtual ~BiCluster();

public:
    double mPrior1, mPrior2;

    const gsl_rng_type *T;
    gsl_rng *gammar;

    vector<vector<int> > lociCluster;
    int clusterCn;
    int minDistance; //minimum distance to allow for interactions
    int mMCMC, mBurnin, mThin;
    vector<int> mOrder, mOrderMap; //order SNPs according to there marginal effects
    double mRate; //pdf(x)=x^rate, rate = 0 means uniform distribution
    vector<bool> UMark; //to replace indState for controls

    vector<vector<double> > proBase;
    vector<double> logGamma; //log factorial indeed
    vector<double> logGammaH;
    vector<double> logGammaS;
    vector<vector<double> > baseGamma;
    vector<vector<double> > alleleFreq;
    vector<vector<int> > MCcount, MCcountU;
    vector<vector<int> > SCcount, SCcountU;
    vector<double> singleDLP, singleLP;
    int alleleCn;
    bool allVectors;	//whether consider control data as genotype vectors,
    //it is needed for permutation tests
    //when using this model, we need an additional boolean variable to indicate
    //that whether this cluster of interactions are specific to the disease sample or not
    bool verbose;

private:
    inline bool _invalidInt(int i, int j, int L)
    {
        return ((int)lociCluster.size() >= L && lociCluster[i][0] == lociCluster[j][0] && abs(lociCluster[i][1] - lociCluster[j][1]) < minDistance);
    };

public:
    double testIndependence_twogroups_LR(vector<vector<char> > const &dataU, vector<int> const &loci1, vector<int> const &loci2, vector<int> const &taken,
                                         double &mle, double &df);
    bool testIndependence_full_LR(vector<vector<char> > const &dataU, vector<int> const &loci, vector<int> const &taken, double P, int radius);
    double LLR(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> const &dloci);
    double chiSquare(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> &dloci);
    void hierarchicalPower(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, double pvalue, int method,
                           vector<vector<vector<int> > > &markers, vector<OUTPUT> &output);

    double conditionalLLR(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> const &Dloci, vector<bool> const &mask);

    void computeConstant1(int Nd, int Nu, double alpha, int intsz, vector<double> &results);

    bool updateLoci_S_Multiple(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                               PARAMETER &para, double p0, double p1, double T);

    double getHapLP1(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, PARAMETER &para, bool model);

    double getHapLP(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, PARAMETER &para);
    double getHapLP_sub(vector<vector<char> > const &data, vector<bool> const &indState, vector<int> const &Dloci,
                        vector<int> const &hapCn, vector<int> const &hapId, int datatype);
    void ChiSquare(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> const &Dloci,
                   vector<int> &countC, vector<int> &countU, double &chi, double &pvalue);
    void getAssociations(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU,
                         vector<int> const &Ssamples, vector<int> const &Isamples, vector<vector<int> > const &positions,
                         double Pvalue, int validationmethod, vector<OUTPUT> &output);
    int runChains(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, double P, vector<OUTPUT> &output,
                  vector<int> &totalcounts, vector<double> &group1, vector<double> &group2, int indRuns, vector<vector<int> > const &positions = vector<vector<int> >(), bool multiTry = false, bool singleOnly = false, int validationmethod = 0,
                  int startTries = 0, int tryLen = 10);
    double computeConditionalLnr(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU,
                                 vector<int> const &Dloci, vector<bool> const &mask);
    double computeLnr(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> const &Dloci);
    double computeLnr(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, PARAMETER &para, vector<int> const &Dloci);

    void initializeGamma(int NC, int NU, int L, int alleleCn, int locilimit);
    void updateDsp(PARAMETER &para);
    double getHapLP_step1(PARAMETER &para);
    double getHapLP_step2(PARAMETER &para);
    void combineHapCount(vector<HAPLOTYPE> const &hapsC, vector<HAPLOTYPE> const &hapsU, vector<int> const &hapIdC, vector<int> const &hapIdU,
                         int lcn, vector<int> &hapcount, vector<int> &hapid);

    bool updateLoci_DS_Multiple(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                                PARAMETER &para, double T, int locilimit);
    bool updateLoci_Multiple(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                             PARAMETER &para, double p0, double p1, double p2, double T, int locilimit);
    bool updateLoci_DS(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                       PARAMETER &para, double T, int locilimit);
    bool updateLoci_S(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                      PARAMETER &para, double p0, double p1, double T);
    double getSLP_MC(vector<vector<char> > const &data, vector<bool> const &indState, int loci);
    double getCompleteLP_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                            PARAMETER &para);
    void computeAlleleMC(vector<vector<char> > const &data, vector<vector<char> > const &dataU);
    int parallelTempering_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU, vector<vector<int> > const &positions,
                             vector<int> &Ssamples, vector<int> &Isamples, vector<double> &g1Sz, vector<double> &g2Sz, bool multiTry, bool first, bool singleOnly = false, int run = 1,
                             vector<vector<int> > const &configures = vector<vector<int> >(), vector<double> const &lps = vector<double>());
    void randomStarts(vector<vector<char> > const &data, vector<vector<char> > const &dataU, vector<vector<int> > const &positions,
                      bool singleOnly, int runs, int tryLen, vector<vector<int> > &configures, vector<double> &lps);

    bool updateLoci_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                       PARAMETER &para, double p0, double p1, double p2, double T, int locilimit);
    double getULP_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                     vector<bool> const &indState, int loci, int prevloci);
    double getDLP_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                     vector<bool> const &indState, int loci);

    void getHapRep(vector<vector<char> > const &data, vector<int> const &Dloci, vector<bool> const &indState,
                   vector<HAPLOTYPE> &haps, vector<int> &hapid);
    void incHapRep(vector<vector<char> > const &data, int loci, vector<bool> const &indState,
                   vector<HAPLOTYPE> &haps, vector<int> &hapId);
    void decHapRep(vector<vector<char> > const &data, int loci, vector<int> const &Dloci, vector<bool> const &indState,
                   vector<HAPLOTYPE> &haps, vector<int> &hapId);
    void computeAlleleFreq(vector<vector<char> > const &data);
    void getHapRep(vector<vector<char> > const &data, vector<int> const &Dloci, vector<bool> const &indState,
                   vector<int> &haprep, vector<int> &hapid, vector<int> &hapcn);
    void incHapRep(vector<vector<char> > const &data, int loci, vector<bool> const &indState,
                   vector<int> &hapRep, vector<int> &hapId, vector<int> &hapCn);
    void decHapRep(vector<vector<char> > const &data, int loci, vector<int> const &Dloci, vector<bool> const &indState,
                   vector<int> &hapRep, vector<int> &hapId, vector<int> &hapCn);


    void logicR_counts(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU,
                       vector<int> const &dloci, vector<int> &countC, vector<int> &countU);
    void logicR_estPara(vector<int> const &countC, vector<int> const &countU, vector<double> &coeff, vector<double> &lps);
    double logicR_logLikelihood(vector<int> const &countC, vector<int> const &countU, vector<double> const &lps);
    void associationScore_test(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU,
                               int degree, char *filename);//vector<double> &statistics);
};
