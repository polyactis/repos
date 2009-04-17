// BiCluster.cpp: implementation of the BiCluster class.
//
#include "BiCluster.h"

double logprop[3];
double Alphan;
double delta;
double gsldelta;

char *LnPfile = "lnp.txt";
char *Posteriorfile = "posterior.txt";

bool useConstMarginPrior = true;
int mNC = 0, mNU = 0;

int minlimit;
int TrueLimit;
int mixture = 1; //0:	use independence marker model
//1:	use a mixture of independence model and 1st order Markov chain
//2:	use 1st order Markov chain
double autoRestart = 10;

typedef struct MySortTypeb {
    double value;
    int indn;
    int indices[5];
} MySortTypeb;

bool operator<(const MySortTypeb &a, const MySortTypeb &b)
{
    return a.value < b.value;
}

BiCluster::BiCluster(int acn)
{
    mPrior1 = mPrior2 = 0.001;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    gammar = gsl_rng_alloc(T);
    gsl_rng_set(gammar, (unsigned int)time(NULL));

    verbose = true;

    alleleCn = acn;
    Alphan = 0.5 * (double)acn;
    delta = Alphan / (double)alleleCn;
    gsldelta = gsl_sf_lngamma(delta);

    allVectors = true;//
    minDistance = 1000000; //1Mb

    mBurnin = 100000;//5000;
    mMCMC = 900000;//15000;
    mThin = 1000;
}

BiCluster::~BiCluster()
{
    gsl_rng_free(gammar);
}

void BiCluster::getHapRep(vector<vector<char> > const &data, vector<int> const &Dloci, vector<bool> const &indState,
                          vector<int> &haprep, vector<int> &hapid, vector<int> &hapcn)
{
    int i, j, k;
    int N = (int)data.size();
    int l = (int)Dloci.size();

    haprep.clear();
    hapid.clear();
    hapcn.clear();
    hapid.resize(N, -1);
    if (l == 0) return;
    for (i = 0; i < N; i++)
        if ((int)indState[i])
        {
            for (k = 0; k < (int)haprep.size(); k++)
            {
                for (j = 0; j < l; j++)
                    if (data[i][Dloci[j]] != data[haprep[k]][Dloci[j]])
                        break;
                if (j >= l)
                    break;
            }
            if (k >= (int)haprep.size())
            {
                haprep.push_back(i);
                hapcn.push_back(1);
            }
            else
                hapcn[k]++;
            hapid[i] = k;
        }
}


void BiCluster::getHapRep(vector<vector<char> > const &data, vector<int> const &Dloci, vector<bool> const &indState,
                          vector<HAPLOTYPE> &haps, vector<int> &hapid)
{
    int i, j, k;
    int N = (int)data.size();
    int l = (int)Dloci.size();

    haps.clear();
    hapid.clear();
    hapid.resize(N, -1);
    for (i = 0; i < N; i++)
        if (indState[i])
        {
            for (k = 0; k < (int)haps.size(); k++)
            {
                for (j = 0; j < l; j++)
                    if (data[i][Dloci[j]] != data[haps[k].repID][Dloci[j]])
                        break;
                if (j >= l)
                    break;
            }
            if (k >= (int)haps.size())
            {
                HAPLOTYPE nh;
                nh.repID = i;
                nh.hCn = 1;
                nh.bFreq = 1.;
                nh.code = 0;
                for (j = 0; j < (int)Dloci.size(); j++)
                {	//tmpremove nh.bFreq *= alleleFreq[Dloci[j]][(int)data[i][Dloci[j]]];
                    nh.code = nh.code * alleleCn + (int)data[i][Dloci[j]];
                }
                haps.push_back(nh);
            }
            else
                haps[k].hCn++;
            hapid[i] = k;
        }
}

void BiCluster::incHapRep(vector<vector<char> > const &data, int loci, vector<bool> const &indState,
                          vector<int> &hapRep, vector<int> &hapId, vector<int> &hapCn)
{
    int sz = (int)hapRep.size();
    int N = (int)data.size();
    int hsz = alleleCn * sz;
    if (hsz == 0) hsz = alleleCn;
    vector<int> newrep(hsz, -1);
    vector<int> newid(N);
    vector<int> newcn(hsz, 0);
    vector<int> pointer(hsz, -1);
    int i, j = 0, k;
    for (i = 0; i < N; i++)
        if (indState[i])
        {
            if (sz == 0) k = (int)data[i][loci];
            else
            {
                k = hapId[i];
                k += sz * (int)data[i][loci];
            }
            if (pointer[k] < 0)
            {
                newrep[j] = i;
                pointer[k] = j;
                j++;
            }
            newcn[pointer[k]]++;
            newid[i] = pointer[k];
        }

    for (i = alleleCn * sz - 1; i > 0; i--)
        if (newrep[i] >= 0)
            break;
    newrep.resize(i + 1);
    newcn.resize(i + 1);

    hapRep = newrep;
    hapCn = newcn;
    hapId = newid;
}

void BiCluster::incHapRep(vector<vector<char> > const &data, int loci, vector<bool> const &indState,
                          vector<HAPLOTYPE> &haps, vector<int> &hapId)
{
    int sz = (int)haps.size();

    int N = (int)data.size();
    int hsz = alleleCn * sz;
    if (hsz == 0) hsz = alleleCn;
    vector<HAPLOTYPE> newhaps(hsz);
    vector<int> newid(N);
    vector<int> pointer(hsz, -1);
    int i, j = 0, k, t;

    for (i = 0; i < N; i++)
        if (indState[i])
        {
            if (sz == 0) t = (int)data[i][loci];
            else
            {
                k = hapId[i];
                t = (int)data[i][loci] * sz + k;
            }
            if (pointer[t] < 0)
            {
                newhaps[j].repID = i;
                newhaps[j].hCn = 0;
                newhaps[j].code = (int)data[i][loci];
                if (sz > 0)
                    newhaps[j].code += haps[k].code * alleleCn;
                //tmpremove
                /*if((int)haps.size() > 0)
                	newhaps[j].bFreq = haps[k].bFreq * alleleFreq[loci][data[i][loci]];
                else newhaps[j].bFreq = alleleFreq[loci][data[i][loci]];
                */
                pointer[t] = j;
                j++;
            }
            if (pointer[t] < 0 || pointer[t] >= (int)newhaps.size())
                i=i;
            newhaps[pointer[t]].hCn++;
            newid[i] = pointer[t];
        }
    newhaps.resize(j);

    haps = newhaps;
    hapId = newid;
}

void BiCluster::decHapRep(vector<vector<char> > const &data, int loci, vector<int> const &Dloci, vector<bool> const &indState,
                          vector<int> &hapRep, vector<int> &hapId, vector<int> &hapCn)
{
    int i, j, k;
    if ((int)Dloci.size() == 1)
    {
        hapRep.clear();
        hapCn.clear();
        for (i = 0; i < (int)hapId.size(); i++)
            hapId[i] = -1;
        return;
    }

    vector<int> idchange((int)hapRep.size(), -1);
    for (i = 0; i < (int)hapRep.size() - 1; i++)
    {
        int x = hapRep[i];
        if (idchange[i] >= 0) continue;
        for (j = i + 1; j < (int)hapRep.size(); j++)
        {
            int y = hapRep[j];
            for (k = 0; k < (int)Dloci.size(); k++)
                if (Dloci[k] != loci && data[x][Dloci[k]] != data[y][Dloci[k]])
                    break;
            if (k >= (int)Dloci.size()) //two are identical
            {
                idchange[j] = i;
                hapCn[i] += hapCn[j];
                //	break;
            }
        }
    }
    vector<int> newRep;
    vector<int> newCn;
    vector<int> newid((int)hapRep.size(), -1);
    for (i = 0; i < (int)hapRep.size(); i++)
    {
        if (idchange[i] < 0)
        {
            newRep.push_back(hapRep[i]);
            newCn.push_back(hapCn[i]);
            newid[i] = (int)newRep.size() - 1;
        }
    }
    for (i = 0; i < (int)hapId.size(); i++)
        if (hapId[i] >= 0)
        {
            if (idchange[hapId[i]] >= 0)
                hapId[i] = newid[idchange[hapId[i]]];
            else hapId[i] = newid[hapId[i]];
        }
    hapRep = newRep;
    hapCn = newCn;
}

void BiCluster::decHapRep(vector<vector<char> > const &data, int loci, vector<int> const &Dloci, vector<bool> const &indState,
                          vector<HAPLOTYPE> &haps, vector<int> &hapId)
{
    int i, j, k;
    if ((int)Dloci.size() == 1)
    {
        haps.clear();
        for (i = 0; i < (int)hapId.size(); i++)
            hapId[i] = -1;
        return;
    }

    vector<int> idchange((int)haps.size(), -1);
    //tmpremove
    //for(i = 0; i < (int)haps.size(); i++)
    //	haps[i].bFreq /= alleleFreq[loci][(int)data[haps[i].repID][loci]];
    for (i = 0; i < (int)Dloci.size(); i++)
        if (Dloci[i] == loci)
            break;
    k = i;
    int mod = (int)(pow((double)alleleCn, (double)Dloci.size() - k - 1.) + 0.5);
    int sz = (int)(pow((double)alleleCn, (double)Dloci.size()) + 0.5);
    vector<int> ch(sz, -1);
    for (i = 0; i < (int)haps.size(); i++)
    {
        k = haps[i].code - (int)data[haps[i].repID][loci] * mod;
        if (ch[k] >= 0)
        {
            idchange[i] = ch[k];
            haps[ch[k]].hCn += haps[i].hCn;
        }
        else ch[k] = i;
    }

    vector<HAPLOTYPE> newhaps;
    vector<int> newid((int)haps.size(), -1);
    j = 0;
    for (i = 0; i < (int)haps.size(); i++)
    {
        if (idchange[i] < 0)
        {
            newhaps.push_back(haps[i]);
            newid[i] = (int)newhaps.size() - 1;
            newhaps[j].code = newhaps[j].code % mod + (newhaps[j].code - newhaps[j].code % (mod * alleleCn)) / alleleCn;
            j++;
        }
    }
    for (i = 0; i < (int)hapId.size(); i++)
        if (hapId[i] >= 0)
        {
            if (idchange[hapId[i]] >= 0)
                hapId[i] = newid[idchange[hapId[i]]];
            else hapId[i] = newid[hapId[i]];
        }
    haps = newhaps;
}

void BiCluster::computeAlleleFreq(vector<vector<char> > const &data)
{
    int N = (int)data.size();
    int L = (int)data[0].size();
    alleleFreq.clear();
    alleleFreq.resize(L, vector<double>(alleleCn, 0));
    int i, j;
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < N; j++)
            alleleFreq[i][data[j][i]]++;
        for (j = 0; j < alleleCn; j++)
            alleleFreq[i][j] /= (double)N;
    }
}

//compute non-disease loci likelihood based on Markov Chain
double BiCluster::getULP_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                            vector<bool> const &indState, int loci, int prevloci)
{
    int i, j;
    int N = (int)data.size();
    int NU = (int)dataU.size();
    double rt = 0;

    if (loci >= (int)data[0].size())
        return 0;

    double rt1 = 0;
    if (mixture < 2)
    {
        vector<int> cn(alleleCn, 0);
        vector<int> cnU(alleleCn, 0);
        cn = SCcount[loci];
        cnU = SCcountU[loci];

        for (i = 0; i < alleleCn; i++)
        {
            rt1 += logGammaS[cn[i] + cnU[i]]; //
            //if(cn[i] + cnU[i] > 0)
            rt1 -= gsldelta;
        }
        rt1 -= logGamma[N + NU]; //changed
        rt1 += logGamma[0];
    }
    if (mixture > 0)
    {
        if (loci == 0 || prevloci < 0)
        {
            vector<int> cn(alleleCn, 0);
            vector<int> cnU(alleleCn, 0);
            cn = MCcount[loci];
            cnU = MCcountU[loci];

            for (i = 0; i < alleleCn; i++)
            {
                rt += logGammaS[cn[i] + cnU[i]]; //
                //if(cn[i] + cnU[i] > 0)
                rt -= gsldelta;
            }
            rt -= logGamma[N + NU]; //changed
            rt += logGamma[0];
        }
        else
        {
            vector<int> cn(alleleCn * alleleCn, 0);
            vector<int> cnU(alleleCn * alleleCn, 0);
            cn = MCcount[loci];
            cnU = MCcountU[loci];

            vector<int> sumn(alleleCn, 0);
            vector<int> sumnd(alleleCn, 0);//
            for (i = 0; i < (int)cn.size(); i += alleleCn)
                for (j = i; j < i + alleleCn; j++)
                {
                    sumn[i / alleleCn] += cnU[j];
                    sumnd[i / alleleCn] += cn[j]; //
                }
            for (i = 0; i < (int)cn.size(); i++)
            {
                rt += logGammaS[cn[i] + cnU[i]];//
                //if(cn[i] + cnU[i] > 0)
                rt -= gsldelta; //updated
            }
            for (i = 0; i < alleleCn; i++)//
            {
                rt -= logGamma[sumn[i] + sumnd[i]];//changed
                rt += logGamma[0];
            }
        }
    }
    if (mixture == 0) return rt1;
    else if (mixture == 1)
    {
        double mmax = rt;
        if (mmax < rt1) mmax = rt1;
        rt = log((exp(rt - mmax) + exp(rt1 - mmax)) / 2.) + mmax;
        return rt;
    }
    else return rt;
}

//compute disease loci likelihood based on Markov Chain
double BiCluster::getDLP_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                            vector<bool> const &indState, int loci)
{
    int i, j;
    //int N = (int)data.size();
    int NU = (int)dataU.size();
    double rt = 0;

    if (loci >= (int)data[0].size())
        return 0;

    double rt1 = 0;
    if (mixture < 2)
    {
        vector<int> cn(alleleCn, 0);
        vector<int> cnU(alleleCn, 0);
        int k = 0;

        cnU = SCcountU[loci];

        for (i = 0; i < alleleCn; i++)
        {
            rt1 += logGammaS[cn[i] + cnU[i]]; //
            //if(cn[i] + cnU[i] > 0)
            rt1 -= gsldelta;
        }
        rt1 -= logGamma[k + NU]; //changed
        rt1 += logGamma[0];
    }
    if (mixture > 0)
    {
        if (loci == 0)
        {
            vector<int> cn(alleleCn, 0);
            vector<int> cnU(alleleCn, 0);
            int k = 0;
            cnU = MCcountU[loci];

            for (i = 0; i < alleleCn; i++)
            {
                rt += logGammaS[cn[i] + cnU[i]]; //
                //if(cn[i] + cnU[i] > 0)
                rt -= gsldelta;
            }
            rt -= logGamma[k + NU]; //changed
            rt += logGamma[0];
        }
        else
        {
            vector<int> cn(alleleCn * alleleCn, 0);
            vector<int> cnU(alleleCn * alleleCn, 0);

            cnU = MCcountU[loci];

            vector<int> sumn(alleleCn, 0);
            vector<int> sumnd(alleleCn, 0);//
            for (i = 0; i < (int)cn.size(); i += alleleCn)
                for (j = i; j < i + alleleCn; j++)
                {
                    sumn[i / alleleCn] += cnU[j];
                    sumnd[i / alleleCn] += cn[j]; //
                }
            for (i = 0; i < (int)cn.size(); i++)
            {
                rt += logGammaS[cn[i] + cnU[i]];//
                //if(cn[i] + cnU[i] > 0)
                rt -= gsldelta;//updated
            }
            for (i = 0; i < alleleCn; i++)//
            {
                rt -= logGamma[sumn[i] + sumnd[i]];//changed
                rt += logGamma[0];
            }
        }
    }
    if (mixture == 0) return rt1;
    else if (mixture == 1)
    {
        double mmax = rt;
        if (mmax < rt1) mmax = rt1;
        rt = log((exp(rt - mmax) + exp(rt1 - mmax)) / 2.) + mmax;
        return rt;
    }
    else return rt;
}

//update loci using three moves
bool BiCluster::updateLoci_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                              PARAMETER &para, double p0, double p1, double p2, double T, int locilimit)
{	//int N = (int)data.size();
    int L = (int)data[0].size();
    int l;
    PARAMETER newpara;
    newpara = para;

    double logratio = 0;
    int movetype;
    int ad, rm, idrm;
    double lp = 0, p;

    ad = rm = -1;
    if ((int)para.Dloci.size() == 0)
        p = p0 + p1 + gsl_ran_flat(gammar, 0, p2);
    else if ((int)para.Dloci.size() < minlimit)
        p = p0 + p1 + gsl_ran_flat(gammar, 0, p2);
    else if ((int)para.Dloci.size() > locilimit)
        p = gsl_ran_flat(gammar, 0, p0);
    else if (minlimit == locilimit)
        p = p0 + gsl_ran_flat(gammar, 0, p1);
    else if ((int)para.Dloci.size() == minlimit)
        p = p0 + gsl_ran_flat(gammar, 0, p1 + p2);
    else if ((int)para.Dloci.size() == locilimit)
        p = gsl_ran_flat(gammar, 0, p0 + p1);
    else p = gsl_ran_flat(gammar, 0, 1);

    if (L <= (int)para.Dloci.size() + (int)para.DSloci.size() || clusterCn <= (int)para.Dloci.size())
        p = 0;

    if (p < p0)
    {
        idrm = (int)gsl_ran_flat(gammar, 0, (double)para.Dloci.size());
        rm = para.Dloci[idrm];
        movetype = 0;
        if ((int)para.Dloci.size() == 1)
            logratio = log(1. / p0);
        else if ((int)para.Dloci.size() == minlimit + 1)
            logratio = log(p2 / (p1 + p2) / p0);
        else if ((int)para.Dloci.size() == locilimit)
            logratio = log(p2 * (p0 + p1) / p0);
        else logratio = log(p2 / p0);
        logratio += log((double)para.Dloci.size() / ((double)L - (double)para.Dloci.size() - (double)para.DSloci.size() + 1.));
    }
    else if (p < p0 + p1)
    {
        idrm = (int)gsl_ran_flat(gammar, 0, (double)para.Dloci.size());
        rm = para.Dloci[idrm];
        bool a_inner;
        int ccc = 0;
        do {
            ad = (int)gsl_ran_flat(gammar, 0, L);
            a_inner = true;
            if (para.lociStatus[ad] == 0)
            {
                int t;
                for (t = 0; t < (int)para.Dloci.size(); t++)
                {
                    if (para.Dloci[t] == rm) continue;
                    //if(lociCluster[ad] == lociCluster[para.Dloci[t]])
                    if (_invalidInt(ad, para.Dloci[t], L))
                        break;
                }
                if (t >= (int)para.Dloci.size())
                    a_inner = false;
            }
            ccc ++;
        } while (a_inner && ccc < 200);
        if (a_inner) return false;
        movetype = 1;
    }
    else
    {
        bool a_inner;
        int ccc = 0;
        do {
            ad = (int)gsl_ran_flat(gammar, 0, L);
            a_inner = true;
            if (para.lociStatus[ad] == 0)
            {
                int t;
                for (t = 0; t < (int)para.Dloci.size(); t++)
                {	//if(lociCluster[ad] == lociCluster[para.Dloci[t]])
                    if (_invalidInt(ad, para.Dloci[t], L))
                        break;
                }
                if (t >= (int)para.Dloci.size())
                    a_inner = false;
            }
            ccc++;
        } while (a_inner && ccc < 200);
        if (a_inner) return false;
        movetype = 2;
        if ((int)para.Dloci.size() == 0)
            logratio = log(p0);
        else if ((int)para.Dloci.size() == minlimit)
            logratio = log(p0 * (p1 + p2) / p2);
        else if ((int)para.Dloci.size() == locilimit - 1)
            logratio = log(p0 / (p0 + p1) / p2);
        else logratio = log(p0 / p2);
        logratio += log(((double)L - (double)para.Dloci.size() - (double)para.DSloci.size()) / ((double)para.Dloci.size() + 1.));
    }

    if (movetype < 2) //genetic -> phenocopy
    {
        decHapRep(data, rm, newpara.Dloci, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
            decHapRep(dataU, rm, newpara.Dloci, UMark, newpara.hapsU, newpara.hapIdU);

        if (!allVectors) lp -= getDLP_MC(data, dataU, para.indState, rm);
        lp += getULP_MC(data, dataU, para.indState, rm, rm - 1);

        l = idrm;
        newpara.Dloci.erase(newpara.Dloci.begin() + l, newpara.Dloci.begin() + l + 1);
        newpara.nlCn.erase(newpara.nlCn.begin() + l, newpara.nlCn.begin() + l + 1);
        //	newpara.lociStatus[rm] = 0;
    }

    if (movetype > 0) //phenocopy -> genetic
    {
        incHapRep(data, ad, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
            incHapRep(dataU, ad, UMark, newpara.hapsU, newpara.hapIdU);

        lp -= getULP_MC(data, dataU, para.indState, ad, ad - 1);
        if (!allVectors) lp += getDLP_MC(data, dataU, para.indState, ad);

        newpara.Dloci.push_back(ad);
        vector<int> tc(alleleCn, 0);

        newpara.nlCn.push_back(tc);
        //	newpara.lociStatus[ad] = 2;
    }

    //update log likelihood
    double lp2 = lp + getHapLP_step1(newpara) - getHapLP_step1(para) + getHapLP_step2(newpara) - getHapLP_step2(para);

    if (movetype == 0)
    {
        double dp = logprop[0] - logprop[2];
//		lp1 += dp;
        lp2 += dp;
    }
    else if (movetype == 2)
    {
        double dp = logprop[2] - logprop[0];
//		lp1 += dp;
        lp2 += dp;
    }
    double un = gsl_ran_flat(gammar, 0, 1);
    if (lp2 / T + logratio >= 0 || exp((lp2 / T + logratio)) >= un || (movetype == 0 && (int)para.Dloci.size() > TrueLimit))
    {
        newpara.logP += lp2;
        newpara.Dsp = true;
        para = newpara;
        if (movetype < 2) para.lociStatus[rm] = 0;
        if (movetype > 0) para.lociStatus[ad] = 2;
        return true;
    }
    return false;
}

//main loop
int BiCluster::parallelTempering_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU, vector<vector<int> > const &positions,
                                    vector<int> &Ssamples, vector<int> &Isamples, vector<double> &g1Sz, vector<double> &g2Sz, bool multiTry, bool first, bool singleOnly, int run,
                                    vector<vector<int> > const &configures, vector<double> const &lps)
{
    int NU = (int)dataU.size();
    int N = (int)data.size();
    int L = (int)data[0].size();
    mNC = N;
    mNU = NU;

    int meansz1, meansz2;
    meansz1 = max(1, (int)((double)L * mPrior1));
    meansz2 = max(1, (int)((double)L * mPrior2));
    if (singleOnly) meansz2 = 0;

    logprop[0] = log((double)(L - meansz1 - meansz2) / (double)L + 0.000000001);
    logprop[1] = log((double)meansz1 / (double)L + 0.000000001);
    logprop[2] = log((double)meansz2 / (double)L + 0.000000001);

    int locilimit = (int)(log((double)N) / log((double)alleleCn)) - 1;
    TrueLimit = (int)(log((double)N / 10.) / log((double)alleleCn));
    if (locilimit <= 0) locilimit = 1;
    if (TrueLimit <= 0) TrueLimit = 1;

    int chain = 1;
    vector<PARAMETER> para(chain);
    vector<int> chainmap(chain);
    vector<double> T(chain);
    int i, j, k, l, m;

    //define positions if not specified, and obtain the number of groups
    if ((int)positions.size() == 0)
    {
        lociCluster.resize(L, vector<int>(2));
        for (i = 0; i < L; i++)
        {
            lociCluster[i][0] = i;
            lociCluster[i][1] = 0;
        }
        clusterCn = L;
    }
    else
    {
        lociCluster = positions;
        m = positions[0][0];
        for (i = 1; i < (int)positions.size(); i++)
            if (m < positions[i][0]) m = positions[i][0];
        vector<bool> map(m + 1, false);
        for (i = 0; i < (int)positions.size(); i++)
            map[positions[i][0]] = true;
        clusterCn = 0;
        for (i = 0; i < (int)map.size(); i++)
            if (map[i]) clusterCn++;
    }
    if (locilimit > clusterCn) locilimit = clusterCn;

    //initialize useful numbers
    double un, p;

    if (first)
    {
        initializeGamma(N, NU, L, alleleCn, locilimit);

        computeAlleleFreq(data);
        computeAlleleMC(data, dataU);

        para[0].indState.clear();
        para[0].indState.resize(N, true);

        singleDLP.clear();
        singleDLP.resize(L);
        singleLP.clear();
        singleLP.resize(L);
        for (i = 0; i < L; i++)
            singleDLP[i] = getDLP_MC(data, dataU, para[0].indState, i);
    }

    if (singleOnly) locilimit = TrueLimit = 0;

    vector<double> cumsum = lps;
    if ((int)cumsum.size() > 0)
    {
        double maxlp = cumsum[0];
        for (j = 1; j < (int)cumsum.size(); j++)
            if (maxlp < cumsum[j]) maxlp = cumsum[j];
        cumsum[0] = exp(cumsum[0] - maxlp);
        for (j = 1; j < (int)cumsum.size(); j++)
            cumsum[j] = cumsum[j - 1] + exp(cumsum[j] - maxlp);
    }
    //initialize random starting points for MCMC chains
    for (i = 0; i < chain; i++)
    {
        chainmap[i] = i;
        T[i] = (double)(i * 0.5 + 1);

        para[i].lociStatus.resize(L, 0);
        para[i].indState.resize(N, false);
        para[i].Dcn = 0;
        if ((int)configures.size() > 0)
        {
            double un = gsl_ran_flat(gammar, 0, *(cumsum.end() - 1));
            for (j = 0; j < (int)cumsum.size() - 1; j++)
                if (un < cumsum[j]) break;
            for (k = 0; k < L; k++)
            {
                if (configures[j][k] == 1)
                {
                    para[i].lociStatus[k] = 1;
                    para[i].DSloci.push_back(k);
                    vector<int> tc(alleleCn, 0);
                    para[i].nlSCn.push_back(tc);
                }
                else if (configures[j][k] == 2)
                {
                    para[i].lociStatus[k] = 2;
                    para[i].Dloci.push_back(k);
                }
            }
        }
        else
        {
            for (j = 0; j < L; j++)
            {
                double uu = gsl_ran_flat(gammar, 0, 1);
                if (uu < (double)meansz1 / L)
                {
                    para[i].lociStatus[j] = 1;
                    para[i].DSloci.push_back(j);
                    vector<int> tc(alleleCn, 0);
                    para[i].nlSCn.push_back(tc);
                }
            }
            for (j = 0; j < locilimit; j++)
            {
                bool a_inner;
                do{
                    k = (int)gsl_ran_flat(gammar, 0, L);

                    a_inner = true;
                    if (para[i].lociStatus[k] == 0)
                    {
                        for (m = 0; m < (int)para[i].Dloci.size(); m++)
                            if (_invalidInt(k, para[i].Dloci[m], L))
                                break;
                        if (m >= (int)para[i].Dloci.size())
                            a_inner = false;
                    }
                } while (a_inner);
                para[i].lociStatus[k] = 2;
                para[i].Dloci.push_back(k);
            }
        }

        p = 1.;
        for (j = 0; j < N; j++)
        {
            un = gsl_ran_flat(gammar, 0, 1.);
            if (un < p)
            {
                para[i].indState[j] = true;
                para[i].Dcn++;
            }
        }
        para[i].nlCn.resize((int)para[i].Dloci.size(), vector<int>(alleleCn, 0));
        getHapRep(data, para[i].Dloci, para[i].indState, para[i].haps, para[i].hapId);
        if (allVectors)
            getHapRep(dataU, para[i].Dloci, UMark, para[i].hapsU, para[i].hapIdU);
        para[i].Dsp = true;

        for (j = 0; j < (int)para[i].Dloci.size(); j++)
        {
            for (k = 0; k < N; k++)
                if (!para[i].indState[k])
                    para[i].nlCn[j][(int)data[k][para[i].Dloci[j]]]++;
        }

        getCompleteLP_MC(data, dataU, para[i]);
    }

    vector<int> swc(chain, 0), swn(chain, 0);
    vector<int> chainact(chain, 0), chainatt(chain, 0);

    //gibbs-sampling
    int burnin = mBurnin;
    int mcmc = mMCMC;
    int thin = mThin;
    int lociCn = 1;
    Ssamples.clear();
    Ssamples.resize(L, 0);
    Isamples.clear();
    Isamples.resize(L, 0);
    g1Sz.clear();
    g2Sz.clear();

    //vector<int> dcnsz;
    int actn[3] = {0}, attn[3] = {0}, actn1, attn1;
    int actb, attb;
    actb = actn1 = 0;
    attb = attn1 = 1;
    vector<vector<int> > dszfn(chain, vector<int>(locilimit + 1, 0));
    vector<vector<int> > dszfnS(chain, vector<int>(L + 1, 0));

    if (verbose && run == 1)
    {
        FILE *f = fopen(LnPfile, "w");
        fclose(f);
    }
    vector<double> tlogp(1000);
    int samplesize = 0;
    double initialLp = para[chainmap[0]].logP;

    printf("Run%d:\n", run);
    FILE *ff = fopen(LnPfile, "a");
    fprintf(ff, "0\t%f\n", para[chainmap[0]].logP);
    fclose(ff);
    for (i = 0; i < burnin + mcmc; i++)
    {
        T[0] = 1.;
        minlimit = 0;
        if (chain == 1 && i < burnin / 1.5 && (int)configures.size() == 0) //if run a single chain, annealing is used
            T[0] += (double)log((double)(N + NU)) * (double)(burnin / 1.5 - i) * 1.5 / (double)burnin;

        if (i < burnin && !singleOnly) //lower bound of interaction size is applied in burnin
            minlimit = (int)((double)(locilimit - 1) * (1. - exp(-(double)(burnin - i) / burnin)) / (1. - exp(-1.))) + 1;

        if ((i % 10000) == 0) //debugging only, output progress
        {
            printf("%d:\tlnP= %f\tl0= %d, l1= %d, l2=%d (", i, para[chainmap[0]].logP, L - (int)para[chainmap[0]].DSloci.size() - (int)para[chainmap[0]].Dloci.size(), (int)para[chainmap[0]].DSloci.size(), (int)para[chainmap[0]].Dloci.size());
            for (j = 0; j < (int)para[chainmap[0]].Dloci.size(); j++)
                printf("%d,", para[chainmap[0]].Dloci[j]);
            printf(")\n");
        }
        if (verbose && (i % 1000) == 999)
        {
            tlogp[(i / 1000) % 1000] = para[chainmap[0]].logP;
            if ((i / 1000) % 1000 == 999)
            {
                ff = fopen(LnPfile, "a");
                for (j = 0; j < 1000; j++)
                    //fprintf(ff, "%d\t%f\n", i-999000+j * 1000, tlogp[j]);
                    fprintf(ff, "%d\t%f\n", i + 1 - ((i / 1000) % 1000) * 1000+j * 1000, tlogp[j]);
                fclose(ff);
            }
        }

        //ran parallele chains with different temporatures
        for (k = 0; k < chain; k++)
        {
            m = chainmap[k];

            //sample loci
            for (j = 0; j < lociCn; j++)
            {
                double sp0 = 0.6;	//sp: frequency of different moves
                double sp1 = 0.3;
                double p0 = 1. / 10.; //p: frequency of different moves for updating interaction group
                double p1 = 8. / 10.;
                double p2 = 1. - p0 - p1;
                if (i < burnin)
                {
                    p0 = 1. / 3. * (double)(i + 1) / burnin;
                    p1 = 1. - p0 * 2.;
                    p2 = 1. - p0 - p1;
                }

                double un = gsl_ran_flat(gammar, 0, 1.);
                bool rt;
                int v = 0;
                if (un < sp0 && !singleOnly)
                {
                    if (multiTry && para[m].Dsp)
                        rt = updateLoci_Multiple(data, dataU, para[m], p0, p1, p2, T[k], locilimit);
                    else rt = updateLoci_MC(data, dataU, para[m], p0, p1, p2, T[k], locilimit);
                    v = 0;
                }
                else if (un < sp0 + sp1 || singleOnly)
                {
                    rt = updateLoci_S(data, dataU, para[m], 0.33333, 0.33333, T[k]);
                    v = 1;
                }
                else
                {
                    if (multiTry && para[m].Dsp)
                        rt = updateLoci_DS_Multiple(data, dataU, para[m], T[k], locilimit);
                    else rt = updateLoci_DS(data, dataU, para[m], T[k], locilimit);
                    v = 2;
                }
                //	if(i >= burnin && allVectors && !singleOnly) updateDsp(para[m]); ///!!!!

                //monitor the MCMC process
                if (i >= burnin)
                {
                    dszfn[k][(int)para[m].Dloci.size()]++;
                    if ((int)para[m].Dloci.size() == 1)
                        dszfnS[k][(int)para[m].DSloci.size() + 1]++;
                    else dszfnS[k][(int)para[m].DSloci.size()]++;
                }
                if (k == 0)
                {
                    if (i < burnin)
                    {
                        if (rt) actb++;
                        attb++;
                    }
                    else
                    {
                        if (rt) actn[v]++;
                        attn[v]++;
                    }
                }
                chainatt[k]++;
                if (rt) chainact[k]++;
                if (k == 0 && i >= burnin && ((i - burnin) % thin) == 0 ) //take one sample every thin iterations
                {
                    samplesize++;
                    if (!singleOnly)
                    {
                        if (para[m].Dsp)
                        {
                            if ((int)para[m].Dloci.size() > 1)
                            {
                                double a = getHapLP_step1(para[m]) + getHapLP_step2(para[m]);
                                para[m].Dsp = false;
                                double b = getHapLP_step1(para[m]) + getHapLP_step2(para[m]);
                                para[m].Dsp = true;
                                double mm = max(a, b);
                                double x = exp(a - mm);
                                double u = gsl_ran_flat(gammar, 0, exp(a - mm) + exp(b - mm));
                                if (u < x) //take the sample only when dsp = true /////////////////////////////////may use better implementation to speed up
                                {
                                    for (int t = 0; t < (int)para[m].Dloci.size(); t++)
                                        Isamples[para[m].Dloci[t]]++;
                                }
                            }
                            else if ((int)para[m].Dloci.size() == 1)
                                Ssamples[para[m].Dloci[0]]++;
                        }
                    }
                    if ((int)para[m].DSloci.size() > 0)
                    {
                        for (int t = 0; t < (int)para[m].DSloci.size(); t++)
                            Ssamples[para[m].DSloci[t]]++;
                    }
                }
            }
        }

        //exchange chains
        if (chain > 1)
        {
            for (k = 0; k < chain / 2; k++)
            {
                l = (int)gsl_ran_flat(gammar, 0, (double)chain - 1.);
                int c1 = chainmap[l];
                int c2 = chainmap[l + 1];
                double rc = (para[c1].logP - para[c2].logP) * (1. / T[l + 1] - 1. / T[l]);

                vector<int> set1, set2, kset1, kset2;
                double lp1, lp2;
                if (true) //ONLY exchange the interaction set (which is equivalent to ONLY exchange the single marker set)
                {
                    lp1 = para[c1].logP;
                    lp2 = para[c2].logP;
                    {
                        int t;
                        for (t = 0; t < (int)para[c1].DSloci.size(); t++)
                        {
                            if (para[c2].lociStatus[para[c1].DSloci[t]] == 0)
                            {
                                set1.push_back(t);
                                lp1 -= singleLP[para[c1].DSloci[t]];
                                lp2 += singleLP[para[c1].DSloci[t]];
                            }
                            else kset1.push_back(para[c1].DSloci[t]);
                        }
                        for (t = 0; t < (int)para[c2].DSloci.size(); t++)
                        {
                            if (para[c1].lociStatus[para[c2].DSloci[t]] == 0)
                            {
                                set2.push_back(t);
                                lp2 -= singleLP[para[c2].DSloci[t]];
                                lp1 += singleLP[para[c2].DSloci[t]];
                            }
                            else kset2.push_back(para[c2].DSloci[t]);
                        }
                        int l1 = (int)para[c1].DSloci.size();
                        int l0 = L - l1 - (int)para[c1].Dloci.size();
                        int d = (int)set2.size() - (int)set1.size();
                        lp1 += (double)d * (logprop[1] -  logprop[0]);

                        l1 = (int)para[c2].DSloci.size();
                        l0 = L - l1 - (int)para[c2].Dloci.size();
                        lp2 += - (double)d * (logprop[1] -  logprop[0]);
                    }
                    rc = (lp1 / T[l + 1] + lp2 / T[l]) - (para[c1].logP / T[l] + para[c2].logP / T[l + 1]);
                }

                un = gsl_ran_flat(gammar, 0, 1.);
                swn[l]++;
                if (rc > 0 || exp(rc) >= un)
                {
                    swc[l]++;
                    m = chainmap[l];
                    chainmap[l] = chainmap[l + 1];
                    chainmap[l + 1] = m;
                    if (true)
                    {
                        para[c1].logP = lp1;
                        para[c2].logP = lp2;
                        int t, m;
                        vector<int> tc(alleleCn, 0);
                        for (t = (int)set1.size() - 1; t >= 0; t--)
                        {
                            m = para[c1].DSloci[set1[t]];
                            para[c1].lociStatus[m] = 0;
                            para[c2].lociStatus[m] = 1;
                            set1[t] = m;
                        }
                        for (t = (int)set2.size() - 1; t >= 0; t--)
                        {
                            m = para[c2].DSloci[set2[t]];
                            para[c1].lociStatus[m] = 1;
                            para[c2].lociStatus[m] = 0;
                            set2[t] = m;
                        }
                        para[c1].DSloci = kset1;
                        para[c1].DSloci.insert(para[c1].DSloci.end(), set2.begin(), set2.end());
                        para[c2].DSloci = kset2;
                        para[c2].DSloci.insert(para[c2].DSloci.end(), set1.begin(), set1.end());
                        para[c1].nlSCn.resize((int)para[c1].DSloci.size());
                        para[c2].nlSCn.resize((int)para[c2].DSloci.size());
                    }
                }
            }
        }
        if (i == burnin && initialLp < para[chainmap[0]].logP) initialLp = para[chainmap[0]].logP;
        if (i > burnin && para[chainmap[0]].logP > initialLp + autoRestart)
        {
            printf("Found a better mode: %f > %f + %f.\nRestart MCMC sampling\n", para[chainmap[0]].logP, initialLp, autoRestart);
            i = burnin;
            initialLp = para[chainmap[0]].logP;
        }
    }
    ff = fopen(LnPfile, "a");
    for (j = 0; j < (i / 1000) % 1000; j++)
        fprintf(ff, "%d\t%f\n", i - ((i / 1000) % 1000) * 1000 + (j + 1) * 1000, tlogp[j]);
    fclose(ff);

    //output summary results
    if (run == 1) ff = fopen(Posteriorfile, "w");
    else ff = fopen(Posteriorfile, "a");
    fprintf(ff, "Run%d:\n", run);
    printf("\n[Summary]\n");
    if (chain > 1)
    {
        for (i = 0; i < chain; i++)
        {
            fprintf(ff, "chain %d: swap frequency= %f\tacceptance= %f\n", i, (double)swc[i] / (double)(swn[i] + 0.001), (double)chainact[i] / (double)chainatt[i]);
            printf("chain %d: swap frequency= %f\tacceptance= %f\n", i, (double)swc[i] / (double)(swn[i] + 0.001), (double)chainact[i] / (double)chainatt[i]);
        }
    }
    fprintf(ff, "Estimate of Posterior Sizes:\n");
    printf("Estimate of Posterior Sizes:\n");
    fprintf(ff, "Marginal Association Size:\n");
    printf("Marginal Association Size:\n");
    for (j = (int)dszfnS[0].size() - 1; j >= 0; j--)
        if (dszfnS[0][j] > 0)
            break;
    for (k = 0; k <= j + 1; k++)
    {
        fprintf(ff, "%d: %1.3f, ", k, (double)dszfnS[0][k] / mcmc);
        printf("%d: %1.3f, ", k, (double)dszfnS[0][k] / mcmc);
        g1Sz.push_back((double)dszfnS[0][k] / mcmc);
    }

    fprintf(ff, "\nInteraction Association Size:\n");
    printf("\nInteraction Association Size:\n");
    for (j = (int)dszfn[0].size() - 1; j >= 0; j--)
        if (dszfn[0][j] > 0)
            break;
    for (k = 0; k <= j; k++)
    {
        if (k == 0) continue;
        if (k < 2)
        {
            fprintf(ff, "0: %1.3f, ", (double)(dszfn[0][0] + dszfn[0][1]) / mcmc);
            printf("0: %1.3f, ", (double)(dszfn[0][0] + dszfn[0][1]) / mcmc);
            g2Sz.push_back((double)(dszfn[0][0] + dszfn[0][1]) / mcmc);
        }
        else
        {
            fprintf(ff, "%d: %1.3f, ", k, (double)dszfn[0][k] / mcmc);
            printf("%d: %1.3f, ", k, (double)dszfn[0][k] / mcmc);
            g2Sz.push_back((double)dszfn[0][k] / mcmc);
        }
    }
    fprintf(ff, "\nMetropolis-Hastings Acceptance Rates:\n");
    printf("\nMetropolis-Hastings Acceptance Rates:\n");

    fprintf(ff, "burnin: %1.3f,\t", (double)actb/(attb + 1));
    printf("burnin= %1.3f,\t", (double)actb/(attb + 1));
    for (i = 0; i < 3; i++)
    {
        fprintf(ff, "move%d= %1.3f,\t", i, (double)actn[i]/(attn[i] + 1));
        printf("move%d= %1.3f,\t", i, (double)actn[i]/(attn[i] + 1));
    }
    fprintf(ff, "\n---------------------------------------------------\n");
    printf("\n---------------------------------------------------\n");
    fclose(ff);

    return samplesize;
}

void BiCluster::randomStarts(vector<vector<char> > const &data, vector<vector<char> > const &dataU, vector<vector<int> > const &positions,
                             bool singleOnly, int runs, int tryLen, vector<vector<int> > &configures, vector<double> &lps)
{
    int NU = (int)dataU.size();
    int N = (int)data.size();
    int L = (int)data[0].size();
    mNC = N;
    mNU = NU;

    int meansz1, meansz2;
    meansz1 = max(1, (int)((double)L * mPrior1));
    meansz2 = max(1, (int)((double)L * mPrior2));
    if (singleOnly) meansz2 = 0;

    logprop[0] = log((double)(L - meansz1 - meansz2) / (double)L + 0.000000001);
    logprop[1] = log((double)meansz1 / (double)L + 0.000000001);
    logprop[2] = log((double)meansz2 / (double)L + 0.000000001);

    int locilimit = (int)(log((double)N) / log((double)alleleCn)) - 1;
    TrueLimit = (int)(log((double)N / 10.) / log((double)alleleCn));
    if (locilimit <= 0) locilimit = 1;
    if (TrueLimit <= 0) TrueLimit = 1;

    int chain = 1;
    vector<PARAMETER> para(chain);
    vector<int> chainmap(chain);
    vector<double> T(chain);
    int i, j, k, l, m;

    //define positions if not specified, and obtain the number of groups
    if ((int)positions.size() == 0)
    {
        lociCluster.resize(L, vector<int>(2));
        for (i = 0; i < L; i++)
        {
            lociCluster[i][0] = i;
            lociCluster[i][1] = 0;
        }
        clusterCn = L;
    }
    else
    {
        lociCluster = positions;
        m = positions[0][0];
        for (i = 1; i < (int)positions.size(); i++)
            if (m < positions[i][0]) m = positions[i][0];
        vector<bool> map(m + 1, false);
        for (i = 0; i < (int)positions.size(); i++)
            map[positions[i][0]] = true;
        clusterCn = 0;
        for (i = 0; i < (int)map.size(); i++)
            if (map[i]) clusterCn++;
    }
    if (locilimit > clusterCn) locilimit = clusterCn;

    //initialize useful numbers
    double un, p;

    {
        initializeGamma(N, NU, L, alleleCn, locilimit);

        computeAlleleFreq(data);
        computeAlleleMC(data, dataU);

        singleDLP.clear();
        singleDLP.resize(L);
        singleLP.clear();
        singleLP.resize(L);
        for (i = 0; i < L; i++)
            singleDLP[i] = getDLP_MC(data, dataU, para[0].indState, i);
    }

    if (singleOnly) locilimit = TrueLimit = 0;

    configures.clear();
    configures.resize(runs);
    lps.clear();
    lps.resize(runs, -100000000);
    printf("Search for good starting points ...\n");
    for (int r = 0; r < runs; r++)
    {	//initialize random starting points for MCMC chains
        for (i = 0; i < chain; i++)
        {
            chainmap[i] = i;
            T[i] = (double)(i * 0.5 + 1);

            PARAMETER pp;
            para[i] = pp;
            para[i].lociStatus.clear();
            para[i].lociStatus.resize(L, 0);
            para[i].indState.resize(N, false);
            para[i].Dcn = 0;
            for (j = 0; j < L; j++)
            {
                double uu = gsl_ran_flat(gammar, 0, 1);
                if (uu < (double)meansz1 / L)
                {
                    para[i].lociStatus[j] = 1;
                    para[i].DSloci.push_back(j);
                    vector<int> tc(alleleCn, 0);
                    para[i].nlSCn.push_back(tc);
                }
            }
            for (j = 0; j < locilimit; j++)
            {
                bool a_inner;
                do{
                    k = (int)gsl_ran_flat(gammar, 0, L);

                    a_inner = true;
                    if (para[i].lociStatus[k] == 0)
                    {
                        for (m = 0; m < (int)para[i].Dloci.size(); m++)
                            if (_invalidInt(k, para[i].Dloci[m], L))
                                break;
                        if (m >= (int)para[i].Dloci.size())
                            a_inner = false;
                    }
                } while (a_inner);
                para[i].lociStatus[k] = 2;
                para[i].Dloci.push_back(k);
            }

            p = 1.;
            for (j = 0; j < N; j++)
            {
                un = gsl_ran_flat(gammar, 0, 1.);
                if (un < p)
                {
                    para[i].indState[j] = true;
                    para[i].Dcn++;
                }
            }
            para[i].nlCn.resize((int)para[i].Dloci.size(), vector<int>(alleleCn, 0));
            getHapRep(data, para[i].Dloci, para[i].indState, para[i].haps, para[i].hapId);
            if (allVectors)
                getHapRep(dataU, para[i].Dloci, UMark, para[i].hapsU, para[i].hapIdU);
            para[i].Dsp = true;

            for (j = 0; j < (int)para[i].Dloci.size(); j++)
            {
                for (k = 0; k < N; k++)
                    if (!para[i].indState[k])
                        para[i].nlCn[j][(int)data[k][para[i].Dloci[j]]]++;
            }

            getCompleteLP_MC(data, dataU, para[i]);
        }

        vector<int> swc(chain, 0), swn(chain, 0);
        vector<int> chainact(chain, 0), chainatt(chain, 0);

        //gibbs-sampling
        int burnin = tryLen;
        int mcmc = 0;
        int thin = mThin;
        int lociCn = 1;

        printf("%d: ", r);
        //vector<int> dcnsz;
        int last = -1;
        for (i = 0; i < burnin + mcmc; i++)
        {
            T[0] = 1.;
            minlimit = 0;
            if (chain == 1 && i < burnin / 1.5) //if run a single chain, annealing is used
                T[0] += (double)log((double)(N + NU)) * (double)(burnin / 1.5 - i) * 1.5 / (double)burnin;

            if (i < burnin && !singleOnly) //lower bound of interaction size is applied in burnin
                minlimit = (int)((double)(locilimit - 1) * (1. - exp(-(double)(burnin - i) / burnin)) / (1. - exp(-1.))) + 1;

            //ran parallele chains with different temporatures
            for (k = 0; k < chain; k++)
            {
                m = chainmap[k];
                //sample loci
                for (j = 0; j < lociCn; j++)
                {
                    double sp0 = 0.6;	//sp: frequency of different moves
                    double sp1 = 0.3;
                    double p0 = 1. / 10.; //p: frequency of different moves for updating interaction group
                    double p1 = 8. / 10.;
                    double p2 = 1. - p0 - p1;
                    if (i < burnin)
                    {
                        p0 = 1. / 3. * (double)(i + 1) / burnin;
                        p1 = 1. - p0 * 2.;
                        p2 = 1. - p0 - p1;
                    }

                    double un = gsl_ran_flat(gammar, 0, 1.);
                    bool rt;
                    int v = 0;
                    if (un < sp0 && !singleOnly)
                    {
                        rt = updateLoci_MC(data, dataU, para[m], p0, p1, p2, T[k], locilimit);
                        v = 0;
                    }
                    else if (un < sp0 + sp1 || singleOnly)
                    {
                        rt = updateLoci_S(data, dataU, para[m], 0.33333, 0.33333, T[k]);
                        v = 1;
                    }
                    else
                    {
                        rt = updateLoci_DS(data, dataU, para[m], T[k], locilimit);
                        v = 2;
                    }
                }
                if (k == 0 && (i == 0 || lps[r] < para[m].logP))
                {
                    configures[r] = para[m].lociStatus;
                    lps[r] = para[m].logP;
                    last = i;
                }
            }
            if (last >= 0 && last > 4 * L && i - last > 4 * L) break; //if stuck in a mode for too long

            //exchange chains
            if (chain > 1)
            {
                for (k = 0; k < chain / 2; k++)
                {
                    l = (int)gsl_ran_flat(gammar, 0, (double)chain - 1.);
                    int c1 = chainmap[l];
                    int c2 = chainmap[l + 1];
                    double rc = (para[c1].logP - para[c2].logP) * (1. / T[l + 1] - 1. / T[l]);

                    vector<int> set1, set2, kset1, kset2;
                    double lp1, lp2;
                    if (true) //ONLY exchange the interaction set (which is equivalent to ONLY exchange the single marker set)
                    {
                        lp1 = para[c1].logP;
                        lp2 = para[c2].logP;
                        {
                            int t;
                            for (t = 0; t < (int)para[c1].DSloci.size(); t++)
                            {
                                if (para[c2].lociStatus[para[c1].DSloci[t]] == 0)
                                {
                                    set1.push_back(t);
                                    lp1 -= singleLP[para[c1].DSloci[t]];
                                    lp2 += singleLP[para[c1].DSloci[t]];
                                }
                                else kset1.push_back(para[c1].DSloci[t]);
                            }
                            for (t = 0; t < (int)para[c2].DSloci.size(); t++)
                            {
                                if (para[c1].lociStatus[para[c2].DSloci[t]] == 0)
                                {
                                    set2.push_back(t);
                                    lp2 -= singleLP[para[c2].DSloci[t]];
                                    lp1 += singleLP[para[c2].DSloci[t]];
                                }
                                else kset2.push_back(para[c2].DSloci[t]);
                            }
                            int l1 = (int)para[c1].DSloci.size();
                            int l0 = L - l1 - (int)para[c1].Dloci.size();
                            int d = (int)set2.size() - (int)set1.size();
                            lp1 += (double)d * (logprop[1] -  logprop[0]);

                            l1 = (int)para[c2].DSloci.size();
                            l0 = L - l1 - (int)para[c2].Dloci.size();
                            lp2 += - (double)d * (logprop[1] -  logprop[0]);
                        }
                        rc = (lp1 / T[l + 1] + lp2 / T[l]) - (para[c1].logP / T[l] + para[c2].logP / T[l + 1]);
                    }

                    un = gsl_ran_flat(gammar, 0, 1.);
                    swn[l]++;
                    if (rc > 0 || exp(rc) >= un)
                    {
                        swc[l]++;
                        m = chainmap[l];
                        chainmap[l] = chainmap[l + 1];
                        chainmap[l + 1] = m;
                        if (true)
                        {
                            para[c1].logP = lp1;
                            para[c2].logP = lp2;
                            int t, m;
                            vector<int> tc(alleleCn, 0);
                            for (t = (int)set1.size() - 1; t >= 0; t--)
                            {
                                m = para[c1].DSloci[set1[t]];
                                para[c1].lociStatus[m] = 0;
                                para[c2].lociStatus[m] = 1;
                                set1[t] = m;
                            }
                            for (t = (int)set2.size() - 1; t >= 0; t--)
                            {
                                m = para[c2].DSloci[set2[t]];
                                para[c1].lociStatus[m] = 1;
                                para[c2].lociStatus[m] = 0;
                                set2[t] = m;
                            }
                            para[c1].DSloci = kset1;
                            para[c1].DSloci.insert(para[c1].DSloci.end(), set2.begin(), set2.end());
                            para[c2].DSloci = kset2;
                            para[c2].DSloci.insert(para[c2].DSloci.end(), set1.begin(), set1.end());
                            para[c1].nlSCn.resize((int)para[c1].DSloci.size());
                            para[c2].nlSCn.resize((int)para[c2].DSloci.size());
                        }
                    }
                }
            }
        }
        m = chainmap[0];
        printf("lnP= %f\tl0= %d, l1= %d, l2=%d (", para[m].logP, L - (int)para[m].DSloci.size() - (int)para[m].Dloci.size(), (int)para[m].DSloci.size(), (int)para[m].Dloci.size());
        for (j = 0; j < (int)para[m].Dloci.size(); j++)
            printf("%d,", para[m].Dloci[j]);
        printf("), #itern = %d\n", i);
    }
    printf("done\n\n");
}

void BiCluster::computeAlleleMC(vector<vector<char> > const &data, vector<vector<char> > const &dataU)
{
    int N = (int)data.size();
    int NU = (int)dataU.size();
    int L = (int)data[0].size();
    int i, j, k;

    SCcount.clear();
    SCcount.resize(L, vector<int>(alleleCn, 0));
    SCcountU.clear();
    SCcountU.resize(L, vector<int>(alleleCn, 0));
    MCcount.clear();
    MCcount.resize(L, vector<int>(alleleCn * alleleCn, 0));
    MCcountU.clear();
    MCcountU.resize(L, vector<int>(alleleCn * alleleCn, 0));

    for (j = 0; j < N; j++)
    {
        k = (int)data[j][0];
        MCcount[0][k]++;
    }
    for (j = 0; j < NU; j++)
    {
        k = (int)dataU[j][0];
        MCcountU[0][k]++;
    }
    SCcount[0] = MCcount[0];
    SCcountU[0] = MCcountU[0];
    SCcount[0].resize(alleleCn);
    SCcountU[0].resize(alleleCn);
    for (i = 1; i < L; i++)
    {
        for (j = 0; j < N; j++)
        {
            k = (int)data[j][i - 1] * alleleCn + (int)data[j][i];
            MCcount[i][k]++;
            SCcount[i][(int)data[j][i]]++;
        }
        for (j = 0; j < NU; j++)
        {
            k = (int)dataU[j][i - 1] * alleleCn + (int)dataU[j][i];
            MCcountU[i][k]++;
            SCcountU[i][(int)dataU[j][i]]++;
        }
    }
}

double BiCluster::getCompleteLP_MC(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                                   PARAMETER &para)
{
    int j;
    int L = (int)data[0].size();
    para.logP = 0;

    double p1;
    para.logpart[0] = para.logpart[1] = para.logpart[2] = 0;
    for (j = 0; j < L; j++)
    {
        p1 = 0;
        if (para.lociStatus[j] > 0)
        {
            if (!allVectors || para.lociStatus[j] < 2)
            {
                p1 = getDLP_MC(data, dataU, para.indState, j);
                if (para.lociStatus[j] == 1)
                    para.logpart[1] += p1;
                else para.logpart[2] += p1;
            }
        }
        else
        {
            p1 = getULP_MC(data, dataU, para.indState, j, j - 1);
            para.logpart[0] += p1;
        }
        para.logP += p1;
    }

    //double p = delta;
    double olp = para.logP;

    para.logP += getHapLP_step1(para);
    para.logP += getHapLP_step2(para);

    para.logpart[2] += para.logP - olp;

    for (j = 0; j < (int)para.DSloci.size(); j++)
    {
        p1 = getSLP_MC(data, para.indState, para.DSloci[j]);
        para.logP += p1;
        para.logpart[1] += p1;
    }

    int l1, l2;
    l1 = (int)para.DSloci.size();
    l2 = (int)para.Dloci.size();
    para.logP += (double)l1 * logprop[1] + (double)l2 * logprop[2] + (double)(L - l1 - l2) * logprop[0];

    return para.logP;
}

double BiCluster::getSLP_MC(vector<vector<char> > const &data, vector<bool> const &indState, int loci)
{
    int i, j, k;
    int N = (int)data.size();
    double rt = 0;

    if (loci >= (int)data[0].size())
        return 0;

    double rt1 = 0;
    if (mixture < 2)
    {
        vector<int> cn(alleleCn, 0);
        cn = SCcount[loci];
        k = N;

        for (i = 0; i < alleleCn; i++)
        {
            rt1 += logGammaS[cn[i]]; //
            //if(cn[i] > 0)
            rt1 -= gsldelta;
        }
        rt1 -= logGamma[k]; //changed
        rt1 += logGamma[0];
    }
    if (mixture > 0)
    {
        if (loci == 0)
        {
            vector<int> cn(alleleCn, 0);
            cn = MCcount[loci];
            k = N;

            for (i = 0; i < alleleCn; i++)
            {
                rt += logGammaS[cn[i]]; //
                //if(cn[i] > 0)
                rt -= gsldelta;
            }
            rt -= logGamma[k]; //changed
            rt += logGamma[0];
        }
        else
        {
            vector<int> cn(alleleCn * alleleCn, 0);
            /*	for(i = 0; i < N; i++)
            		if(indState[i])
            		{	j = (int)data[i][loci - 1] * alleleCn + (int)data[i][loci];
            			cn[j]++;
            		}
            		*/
            cn = MCcount[loci];

            vector<int> sumnd(alleleCn, 0);//
            for (i = 0; i < (int)cn.size(); i += alleleCn)
            {
                for (j = i; j < i + alleleCn; j++)
                    sumnd[i / alleleCn] += cn[j]; //
            }
            for (i = 0; i < (int)cn.size(); i++)
            {
                rt += logGammaS[cn[i]];//
                //if(cn[i] > 0)
                rt -= gsldelta;//updated
            }
            for (i = 0; i < alleleCn; i++)//
            {
                rt -= logGamma[sumnd[i]];//changed
                rt += logGamma[0];
            }
        }
    }
    if (mixture == 0) return rt1;
    else if (mixture == 1)
    {
        double mmax = rt;
        if (mmax < rt1) mmax = rt1;
        rt = log((exp(rt - mmax) + exp(rt1 - mmax)) / 2.) + mmax;
        return rt;
    }
    else return rt;
}

bool BiCluster::updateLoci_S(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                             PARAMETER &para, double p0, double p1, double T)
{	//int N = (int)data.size();
    int L = (int)data[0].size();
    double p2 = 1. - p0 - p1;
//	PARAMETER newpara;
//	newpara = para;

    double logratio = 0;
    int movetype;
    int rm, idrm, ad;
    double lp = 0, p;
    int bound = L - (int)para.Dloci.size();

    ad = rm = -1;
    if ((int)para.DSloci.size() == 0)
        p = 1.;
    else if ((int)para.DSloci.size() >= bound)
        p = gsl_ran_flat(gammar, 0, p0 + p1);
    else p = gsl_ran_flat(gammar, 0, 1);
    if (p < p0)
    {
        idrm = (int)gsl_ran_flat(gammar, 0, (double)para.DSloci.size());
        rm = para.DSloci[idrm];
        movetype = 0;
        if ((int)para.DSloci.size() == 1)
            logratio = log(1. / p0);
        else if ((int)para.DSloci.size() >= bound)
            logratio = log(p2);
        else logratio = log(p2 / p0);
        logratio += log((double)para.DSloci.size() / ((double)L - (double)para.Dloci.size() - (double)para.DSloci.size() + 1.));
    }
    else if (p < p0 + p1)
    {
        idrm = (int)gsl_ran_flat(gammar, 0, (double)para.DSloci.size());
        rm = para.DSloci[idrm];
        bool a_inner;
        do {
            ad = (int)gsl_ran_flat(gammar, 0, L);
            a_inner = true;
            if (para.lociStatus[ad] == 0)
                a_inner = false;
        } while (a_inner);
        movetype = 1;
    }
    else
    {
        bool a_inner;
        do {
            ad = (int)gsl_ran_flat(gammar, 0, L);
            a_inner = true;
            if (para.lociStatus[ad] == 0)
                a_inner = false;
        } while (a_inner);
        movetype = 2;
        if ((int)para.DSloci.size() == 0)
            logratio = log(p0);
        else if ((int)para.DSloci.size() == bound - 1)
            logratio = log(1. / p2);
        else logratio = log(p0 / p2);
        logratio += log(((double)L - (double)para.Dloci.size() - (double)para.DSloci.size()) / ((double)para.DSloci.size() + 1.));
    }

    if (movetype < 2) //genetic -> phenocopy
    {
        lp -= getSLP_MC(data, para.indState, rm);
        lp -= getDLP_MC(data, dataU, para.indState, rm);
        lp += getULP_MC(data, dataU, para.indState, rm, rm - 1);
    }
    if (movetype > 0)//phenocopy -> genetic
    {
        lp -= getULP_MC(data, dataU, para.indState, ad, ad - 1);
        lp += getDLP_MC(data, dataU, para.indState, ad);
        lp += getSLP_MC(data, para.indState, ad);
    }
    if (movetype == 0)
    {
        lp += /*log((double)para.DSloci.size() / ((double)L - (double)para.DSloci.size() - (double)para.Dloci.size() + 1)) + */
            logprop[0] - logprop[1];
    }
    else if (movetype == 2)
    {
        lp += /*log(((double)L - (double)para.DSloci.size() - (double)para.Dloci.size()) / ((double)para.DSloci.size() + 1.)) + */
            logprop[1] - logprop[0];
    }

    //newpara.logP += lp;
    double un = gsl_ran_flat(gammar, 0, 1);
    if (lp / T + logratio >= 0 || exp((lp / T + logratio)) >= un)
    {

        //para = newpara;
        para.logP += lp;
        if (movetype < 2)
        {
            para.DSloci.erase(para.DSloci.begin() + idrm, para.DSloci.begin() + idrm + 1);
            para.nlSCn.erase(para.nlSCn.begin() + idrm, para.nlSCn.begin() + idrm + 1);
            para.lociStatus[rm] = 0;
        }
        if (movetype > 0)
        {
            para.DSloci.push_back(ad);
            vector<int> tc(alleleCn, 0);
            para.nlSCn.push_back(tc);
            para.lociStatus[ad] = 1;
        }
        return true;
    }
    else return false;
}

bool BiCluster::updateLoci_DS(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                              PARAMETER &para, double T, int locilimit)
{
    int N = (int)data.size();
    int L = (int)data[0].size();
    int l;
    PARAMETER newpara;
    newpara = para;

    int l1, l2, idl1, idl2;
    double lp = 0;

    if ((int)para.DSloci.size() == 0 || (int)para.Dloci.size() == 0)
        return false;

    int ccc = 0;
    bool a_inner;
    do{
        idl1 = (int)gsl_ran_flat(gammar, 0, (double)para.DSloci.size());
        l1 = para.DSloci[idl1];
        idl2 = (int)gsl_ran_flat(gammar, 0, (double)para.Dloci.size());
        l2 = para.Dloci[idl2];
        a_inner = true;
        for (l = 0; l < (int)para.Dloci.size(); l++)
        {
            if (para.Dloci[l] == l2) continue;
            //if(lociCluster[l1] == lociCluster[para.Dloci[l]])
            if (_invalidInt(l1, para.Dloci[l], L))
                break;
        }
        if (l >= (int)para.Dloci.size())
            a_inner = false;
        ccc++;
    } while (a_inner && ccc < 1000);
    if (a_inner && ccc >= 1000) return false;

    {
        decHapRep(data, l2, newpara.Dloci, newpara.indState, newpara.haps, newpara.hapId);
        incHapRep(data, l1, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
        {
            decHapRep(dataU, l2, newpara.Dloci, UMark, newpara.hapsU, newpara.hapIdU);
            incHapRep(dataU, l1, UMark, newpara.hapsU, newpara.hapIdU);
        }

        lp -= getDLP_MC(data, dataU, para.indState, l1);
        lp -= getSLP_MC(data, para.indState, l1);
        if (!allVectors) lp += getDLP_MC(data, dataU, para.indState, l1);

        if (!allVectors) lp -= getDLP_MC(data, dataU, para.indState, l2);
        lp += getDLP_MC(data, dataU, para.indState, l2);
        lp += getSLP_MC(data, para.indState, l2);

        //	for(l = 0; l < (int)newpara.Dloci.size(); l++)
        //		if(newpara.Dloci[l] == l2)
        //			break;
        l = idl2;
        newpara.Dloci.erase(newpara.Dloci.begin() + l, newpara.Dloci.begin() + l + 1);
        newpara.nlCn.erase(newpara.nlCn.begin() + l, newpara.nlCn.begin() + l + 1);
        newpara.DSloci.push_back(l2);
        vector<int> tc(alleleCn, 0);
        /*for(l = 0; l < N; l++)
        	if(!newpara.indState[l])
        		tc[(int)data[l][l2]]++;
        		*/ //single cluster
        newpara.nlSCn.push_back(tc);
        //	newpara.lociStatus[l2] = 1;

        //for(l = 0; l < (int)newpara.DSloci.size(); l++)
        //	if(newpara.DSloci[l] == l1)
        //		break;
        l = idl1;
        newpara.DSloci.erase(newpara.DSloci.begin() + l, newpara.DSloci.begin() + l + 1);
        newpara.nlSCn.erase(newpara.nlSCn.begin() + l, newpara.nlSCn.begin() + l + 1);
        newpara.Dloci.push_back(l1);
        tc.clear();
        tc.resize(alleleCn, 0);
        for (l = 0; l < N; l++)
            if (!newpara.indState[l])
                tc[(int)data[l][l2]]++;
        newpara.nlCn.push_back(tc);
        //	newpara.lociStatus[l1] = 2;
    }
//-------------------------------------
    lp += getHapLP_step1(newpara) - getHapLP_step1(para);

    newpara.logP += lp;
    double un = gsl_ran_flat(gammar, 0, 1);
    if (lp >= 0 || exp(lp / T) >= un)
    {
        para = newpara;
        para.lociStatus[l2] = 1;
        para.lociStatus[l1] = 2;
        return true;
    }
    else return false;
}

bool BiCluster::updateLoci_Multiple(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                                    PARAMETER &para, double p0, double p1, double p2, double T, int locilimit)
{	//int N = (int)data.size();
    int L = (int)data[0].size();
    int l;
    PARAMETER newpara;

    double logratio = 0;
    int movetype;
    int ad, rm;
    double lp = 0, p;

    ad = rm = -1;
    if ((int)para.Dloci.size() == 0)
        p = p0 + p1 + gsl_ran_flat(gammar, 0, p2);
    else if ((int)para.Dloci.size() < minlimit)
        p = p0 + p1 + gsl_ran_flat(gammar, 0, p2);
    else if ((int)para.Dloci.size() > locilimit)
        p = gsl_ran_flat(gammar, 0, p0);
    else if (minlimit == locilimit)
        p = p0 + gsl_ran_flat(gammar, 0, p1);
    else if ((int)para.Dloci.size() == minlimit)
        p = p0 + gsl_ran_flat(gammar, 0, p1 + p2);
    else if ((int)para.Dloci.size() == locilimit)
        p = gsl_ran_flat(gammar, 0, p0 + p1);
    else p = gsl_ran_flat(gammar, 0, 1);
    if (p < p0)
    {
        rm = (int)gsl_ran_flat(gammar, 0, (double)para.Dloci.size());
        rm = para.Dloci[rm];
        movetype = 0;
        if ((int)para.Dloci.size() == 1)
            logratio = log(1. / p0);
        else if ((int)para.Dloci.size() == minlimit + 1)
            logratio = log(p2 / (p1 + p2) / p0);
        else if ((int)para.Dloci.size() == locilimit)
            logratio = log(p2 * (p0 + p1) / p0);
        else logratio = log(p2 / p0);
        logratio += log(1. / ((double)L - (double)para.Dloci.size() - (double)para.DSloci.size() + 1.));
    }
    else if (p < p0 + p1)
    {	//rm = (int)gsl_ran_flat(gammar, 0, (double)para.Dloci.size());
        //rm = para.Dloci[rm];
        bool a_inner;
        int cc = 0;
        do {
            ad = (int)gsl_ran_flat(gammar, 0, L);
            a_inner = true;
            if (para.lociStatus[ad] == 0) //consider ad in other groups, otherwise no need for multitry
            {
                int t;
                for (t = 0; t < (int)para.Dloci.size(); t++)
                {	//if(para.Dloci[t] == rm) continue;
                    //if(lociCluster[ad] == lociCluster[para.Dloci[t]])
                    if (_invalidInt(ad, para.Dloci[t], L))
                        break;
                }
                if (t >= (int)para.Dloci.size())
                    a_inner = false;
            }
            cc++;
        } while (a_inner && cc < 1000);
        if (cc >= 1000) return false; //no exchangable markers
        movetype = 1;
    }
    else
    {
        bool a_inner;
        int cc = 0;
        do {
            ad = (int)gsl_ran_flat(gammar, 0, L);
            a_inner = true;
            if (para.lociStatus[ad] == 0)
            {
                int t;
                for (t = 0; t < (int)para.Dloci.size(); t++)
                {	//if(lociCluster[ad] == lociCluster[para.Dloci[t]])
                    if (_invalidInt(ad, para.Dloci[t], L))
                        break;
                }
                if (t >= (int)para.Dloci.size())
                    a_inner = false;
            }
            cc++;
        } while (a_inner && cc < 1000);
        if (cc >= 1000) return false; //no exchangable markers
        movetype = 2;
        if ((int)para.Dloci.size() == 0)
            logratio = log(p0);
        else if ((int)para.Dloci.size() == minlimit)
            logratio = log(p0 * (p1 + p2) / p2);
        else if ((int)para.Dloci.size() == locilimit - 1)
            logratio = log(p0 / (p0 + p1) / p2);
        else logratio = log(p0 / p2);
        logratio += log(((double)L - (double)para.Dloci.size() - (double)para.DSloci.size()));
    }

    //int i;
    if (movetype == 1)
    {
        vector<double> lps((int)para.Dloci.size(), 0);
        double lpmax;
        for (l = 0; l < (int)para.Dloci.size(); l++)
        {
            newpara = para;
            decHapRep(data, newpara.Dloci[l], newpara.Dloci, newpara.indState, newpara.haps, newpara.hapId);
            if (allVectors)
                decHapRep(dataU, newpara.Dloci[l], newpara.Dloci, UMark, newpara.hapsU, newpara.hapIdU);

            //for(i = 0; i < (int)para.haps.size(); i++)
            //	lps[l] -= logGammaH[para.haps[i].hCn];
            if (!allVectors) lps[l] -= getDLP_MC(data, dataU, para.indState, para.Dloci[l]);
            lps[l] += getULP_MC(data, dataU, newpara.indState, para.Dloci[l], para.Dloci[l] - 1);

            incHapRep(data, ad, newpara.indState, newpara.haps, newpara.hapId);
            if (allVectors)
                incHapRep(dataU, ad, UMark, newpara.hapsU, newpara.hapIdU);

            //for(i = 0; i < (int)newpara.haps.size(); i++)
            //	lps[l] += logGammaH[newpara.haps[i].hCn];
            lps[l] -= getULP_MC(data, dataU, newpara.indState, ad, ad - 1);
            if (!allVectors) lps[l] += getDLP_MC(data, dataU, newpara.indState, ad);

            newpara.Dloci.erase(newpara.Dloci.begin() + l, newpara.Dloci.begin() + l + 1);
            newpara.Dloci.push_back(ad);

            lps[l] += getHapLP_step1(newpara) - getHapLP_step1(para);

            //lps[l] -= ((double)newpara.haps.size()) * gsldelta;
            //lps[l] += ((double)para.haps.size()) * gsldelta;
            lps[l] /= T;
            if (l == 0 || lpmax < lps[l]) lpmax = lps[l];
        }

        vector<double> p = lps;
        vector<double> psum = p;
        for (l = 0; l < (int)p.size(); l++)
        {
            if (p[l] - lpmax > -20.)
                p[l] = exp(p[l] - lpmax);
            else p[l] = 0;
            if (l == 0) psum[l] = p[l];
            else psum[l] = p[l] + psum[l - 1];
        }
        double sum = *(psum.end() - 1);
        double un = gsl_ran_flat(gammar, 0, sum);
        for (l = 0; l < (int)psum.size() - 1; l++)
            if (psum[l] > un)
                break;
        rm = para.Dloci[l];
        lp = lps[l] * T;
        double a0 = exp(-lpmax);
        if (sum - p[l] + a0 > 0 && p[l] > 0)
        {
            if (a0 * sum > 0)
                logratio += log(a0 * sum / p[l] / (sum - p[l] + a0)); //debug
            else logratio += -1000000;
        }

//logratio += mRate * (log(mOrderMap[l]) - log(mOrderMap[ad]));

        newpara = para;
        decHapRep(data, rm, newpara.Dloci, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
            decHapRep(dataU, rm, newpara.Dloci, UMark, newpara.hapsU, newpara.hapIdU);
        newpara.Dloci.erase(newpara.Dloci.begin() + l, newpara.Dloci.begin() + l + 1);
        newpara.nlCn.erase(newpara.nlCn.begin() + l, newpara.nlCn.begin() + l + 1);
        //	newpara.lociStatus[rm] = 0;

        incHapRep(data, ad, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
            incHapRep(dataU, ad, UMark, newpara.hapsU, newpara.hapIdU);
        newpara.Dloci.push_back(ad);
        vector<int> tc(alleleCn, 0);
        newpara.nlCn.push_back(tc);
        //	newpara.lociStatus[ad] = 2;
    }
    else if (movetype == 0) //genetic -> phenocopy
    {
        vector<double> lps((int)para.Dloci.size(), 0);
        double lpmax;
        for (l = 0; l < (int)para.Dloci.size(); l++)
        {
            newpara = para;
            decHapRep(data, newpara.Dloci[l], newpara.Dloci, newpara.indState, newpara.haps, newpara.hapId);
            if (allVectors)
                decHapRep(dataU, newpara.Dloci[l], newpara.Dloci, UMark, newpara.hapsU, newpara.hapIdU);

            if (!allVectors) lps[l] -= getDLP_MC(data, dataU, para.indState, para.Dloci[l]);
            lps[l] += getULP_MC(data, dataU, newpara.indState, para.Dloci[l], para.Dloci[l] - 1);

            newpara.Dloci.erase(newpara.Dloci.begin() + l, newpara.Dloci.begin() + l + 1);
            lps[l] += getHapLP_step1(newpara) - getHapLP_step1(para);

            lps[l] /= T;
            if (l == 0 || lpmax < lps[l]) lpmax = lps[l];
        }

        vector<double> p = lps;
        vector<double> psum = p;
        for (l = 0; l < (int)p.size(); l++)
        {
            if (p[l] - lpmax > -20.)
                p[l] = exp(p[l] - lpmax);
            else p[l] = 0;
            if (l == 0) psum[l] = p[l];
            else psum[l] = p[l] + psum[l - 1];
        }
        double sum = *(psum.end() - 1);
        double un = gsl_ran_flat(gammar, 0, sum);
        for (l = 0; l < (int)psum.size() - 1; l++)
            if (psum[l] > un)
                break;
        rm = para.Dloci[l];
        lp = lps[l] * T;
        if (sum > 0)
        {
            if (p[l] > 0)
                logratio += log(sum / p[l]);
            else logratio += 1000000;
        }

//logratio += log(mRate + 1.) + mRate * log(mOrderMap[l]) - (mRate + 1.) * log(L);

        newpara = para;
        decHapRep(data, rm, newpara.Dloci, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
            decHapRep(dataU, rm, newpara.Dloci, UMark, newpara.hapsU, newpara.hapIdU);
        newpara.Dloci.erase(newpara.Dloci.begin() + l, newpara.Dloci.begin() + l + 1);
        newpara.nlCn.erase(newpara.nlCn.begin() + l, newpara.nlCn.begin() + l + 1);
        //	newpara.lociStatus[rm] = 0;
    }
    else //phenocopy -> genetic
    {
        newpara = para;
        incHapRep(data, ad, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
            incHapRep(dataU, ad, UMark, newpara.hapsU, newpara.hapIdU);

        newpara.Dloci.push_back(ad);
        vector<int> tc(alleleCn, 0);
        newpara.nlCn.push_back(tc);
        //	newpara.lociStatus[ad] = 2;
        lp -= getULP_MC(data, dataU, para.indState, ad, ad - 1);
        if (!allVectors) lp += getDLP_MC(data, dataU, para.indState, ad);
        lp += getHapLP_step1(newpara) - getHapLP_step1(para); //why originally I don't have the gsldelta part?
        //lp += getHapLP_step2(newpara) - getHapLP_step2(para);

        vector<double> lps((int)newpara.Dloci.size(), 0);
        double lpmax;
        for (l = 0; l < (int)newpara.Dloci.size(); l++)
        {
            if (l < (int)newpara.Dloci.size() - 1)
            {
                PARAMETER backpara = newpara;
                decHapRep(data, backpara.Dloci[l], backpara.Dloci, backpara.indState, backpara.haps, backpara.hapId);
                if (allVectors)
                    decHapRep(dataU, backpara.Dloci[l], backpara.Dloci, UMark, backpara.hapsU, backpara.hapIdU);

                if (!allVectors) lps[l] -= getDLP_MC(data, dataU, newpara.indState, newpara.Dloci[l]);
                lps[l] += getULP_MC(data, dataU, backpara.indState, newpara.Dloci[l], newpara.Dloci[l] - 1);

                backpara.Dloci.erase(backpara.Dloci.begin() + l, backpara.Dloci.begin() + l + 1);
                lps[l] += getHapLP_step1(backpara) - getHapLP_step1(newpara);
            }
            else lps[l] = -lp;
            lps[l] /= T;
            if (l == 0 || lpmax < lps[l]) lpmax = lps[l];
        }

        double sum = 0;
        for (l = 0; l < (int)lps.size(); l++)
        {
            if (lps[l] - lpmax > -20.)
                lps[l] = exp(lps[l] - lpmax);
            else lps[l] = 0;
            sum += lps[l];
        }
        if (sum > 0)
        {
            if (*(lps.end() - 1) > 0)
                logratio += log(*(lps.end() - 1) / sum);
            else logratio += -1000000;
        }
//logratio -= log(mRate + 1.) + mRate * log(mOrderMap[ad]) - (mRate + 1.) * log(L);
    }
//-------------------------------------
    if (movetype != 1)
    {
        lp += getHapLP_step2(newpara) - getHapLP_step2(para);
        if (movetype == 0)
        {
            lp += /*log((double)para.Dloci.size() / ((double)L - (double)para.DSloci.size() - (double)para.Dloci.size() + 1)) + */
                logprop[0] - logprop[2];
        }
        else if (movetype == 2)
        {
            lp += /*log(((double)L - (double)para.DSloci.size() - (double)para.Dloci.size()) / (double)newpara.Dloci.size()) + */
                logprop[2] - logprop[0];
        }
    }
    newpara.logP += lp;

    double un = gsl_ran_flat(gammar, 0, 1);
    if (lp / T + logratio >= 0 || exp((lp / T + logratio)) >= un || (movetype == 0 && (int)para.Dloci.size() > TrueLimit))
    {
        para = newpara;
        if (movetype > 0) para.lociStatus[ad] = 2;
        if (movetype < 2) para.lociStatus[rm] = 0;
        return true;
    }
    else return false;
}

bool BiCluster::updateLoci_DS_Multiple(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                                       PARAMETER &para, double T, int locilimit)
{	//int N = (int)data.size();
    int L = (int)data[0].size();
    int l;
    //int i;
    PARAMETER newpara;
    newpara = para;

    int l1, l2;
    double lp = 0;
    double logratio = 0;

    if ((int)para.DSloci.size() == 0 || (int)para.Dloci.size() == 0)
        return false;

    bool a_inner = true;
    int cc = 0;
    do {
        l1 = (int)gsl_ran_flat(gammar, 0, (double)para.DSloci.size());
        l1 = para.DSloci[l1];
        for (l = 0; l < (int)para.Dloci.size(); l++)
        {	//if(lociCluster[l1] == lociCluster[para.Dloci[l]])
            if (_invalidInt(l1, para.Dloci[l], L))
                break;
        }
        if (l >= (int)para.Dloci.size())
            a_inner = false;
        cc++;
    } while (a_inner && cc < 1000);
    if (cc >= 1000) return false;

    vector<double> lps((int)para.Dloci.size(), 0), p, psum;
    {
        double lpmax;
        for (l = 0; l < (int)para.Dloci.size(); l++)
        {
            newpara = para;
            decHapRep(data, newpara.Dloci[l], newpara.Dloci, newpara.indState, newpara.haps, newpara.hapId);
            if (allVectors)
                decHapRep(dataU, newpara.Dloci[l], newpara.Dloci, UMark, newpara.hapsU, newpara.hapIdU);

            //for(i = 0; i < (int)para.haps.size(); i++)
            //	lps[l] -= logGammaH[para.haps[i].hCn];
            if (!allVectors)
                lps[l] -= getDLP_MC(data, dataU, para.indState, para.Dloci[l]);

            incHapRep(data, l1, newpara.indState, newpara.haps, newpara.hapId);
            if (allVectors)
                incHapRep(dataU, l1, UMark, newpara.hapsU, newpara.hapIdU);

            //for(i = 0; i < (int)newpara.haps.size(); i++)
            //	lps[l] += logGammaH[newpara.haps[i].hCn];
            if (!allVectors)
                lps[l] += getDLP_MC(data, dataU, newpara.indState, l1);

            lps[l] -= getDLP_MC(data, dataU, para.indState, l1);
            lps[l] -= getSLP_MC(data, para.indState, l1);
            lps[l] += getDLP_MC(data, dataU, newpara.indState, para.Dloci[l]);
            lps[l] += getSLP_MC(data, newpara.indState, para.Dloci[l]);

            newpara.Dloci.erase(newpara.Dloci.begin() + l, newpara.Dloci.begin() + l + 1);
            newpara.Dloci.push_back(l1);

            lps[l] += getHapLP_step1(newpara) - getHapLP_step1(para);
            //lps[l] -= ((double)newpara.haps.size()) * gsldelta;
            //lps[l] += ((double)para.haps.size()) * gsldelta;
            lps[l] /= T;
            if (l == 0 || lpmax < lps[l]) lpmax = lps[l];
        }
        p = lps;
        psum = p;
        for (l = 0; l < (int)p.size(); l++)
        {
            if (p[l] - lpmax > -20.)
                p[l] = exp(p[l] - lpmax);
            else p[l] = 0;
            if (l == 0) psum[l] = p[l];
            else psum[l] = p[l] + psum[l - 1];
        }
        double sum = *(psum.end() - 1);
        double un = gsl_ran_flat(gammar, 0, sum);
        for (l = 0; l < (int)psum.size() - 1; l++)
            if (psum[l] > un)
                break;
        l2 = para.Dloci[l];
        lp = lps[l] * T;
        double a0 = exp(-lpmax);
        if (p[l] > 0 && sum - p[l] + a0 > 0)
        {
            if (a0 * sum > 0)
                logratio += log(a0 * sum / p[l] / (sum - p[l] + a0)); //debug
            else logratio += -1000000;
        }
    }
//-------------------------------------

    double un = gsl_ran_flat(gammar, 0, 1);

    if (lp / T + logratio >= 0 || exp(lp / T + logratio) >= un)
    {
        newpara = para;
        decHapRep(data, l2, newpara.Dloci, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
            decHapRep(dataU, l2, newpara.Dloci, UMark, newpara.hapsU, newpara.hapIdU);

        newpara.Dloci.erase(newpara.Dloci.begin() + l, newpara.Dloci.begin() + l + 1);
        newpara.nlCn.erase(newpara.nlCn.begin() + l, newpara.nlCn.begin() + l + 1);

        incHapRep(data, l1, newpara.indState, newpara.haps, newpara.hapId);
        if (allVectors)
            incHapRep(dataU, l1, UMark, newpara.hapsU, newpara.hapIdU);

        for (l = 0; l < (int)newpara.DSloci.size(); l++)
            if (newpara.DSloci[l] == l1)
                break;
        newpara.DSloci.erase(newpara.DSloci.begin() + l, newpara.DSloci.begin() + l + 1);
        newpara.nlSCn.erase(newpara.nlSCn.begin() + l, newpara.nlSCn.begin() + l + 1);

        newpara.Dloci.push_back(l1);
        newpara.DSloci.push_back(l2);
        vector<int> tc(alleleCn, 0);
        newpara.nlSCn.push_back(tc);
        newpara.nlCn.push_back(tc);

        para = newpara;

        para.logP += lp;
        para.lociStatus[l1] = 2;
        para.lociStatus[l2] = 1;

        return true;
    }
    else return false;
}

//compute total counts of haplotypes in Case + Control data
void BiCluster::combineHapCount(vector<HAPLOTYPE> const &hapsC, vector<HAPLOTYPE> const &hapsU,
                                vector<int> const &hapIdC, vector<int> const &hapIdU,
                                int lcn, vector<int> &hapcount, vector<int> &hapid)
{
    int i, j;
    int sz = (int)(pow((double)alleleCn, (double)lcn) + 0.5);
    vector<int> ch(sz, -1);

    j = 0;
    for (i = 0; i < (int)hapsC.size(); i++)
        ch[hapsC[i].code] = j++;
    for (i = 0; i < (int)hapsU.size(); i++)
        if (ch[hapsU[i].code] < 0)
            ch[hapsU[i].code] = j++;

    hapcount.clear();
    hapcount.resize(j, 0);
    for (i = 0; i < (int)hapsC.size(); i++)
        hapcount[ch[hapsC[i].code]] += hapsC[i].hCn;
    for (i = 0; i < (int)hapsU.size(); i++)
        hapcount[ch[hapsU[i].code]] += hapsU[i].hCn;

    hapid.clear();
    hapid.resize((int)hapIdC.size() + (int)hapIdU.size());
    int NC = (int)hapIdC.size();
    for (i = 0; i < (int)hapIdC.size(); i++)
        hapid[i] = ch[hapsC[hapIdC[i]].code];
    for (i = 0; i < (int)hapIdU.size(); i++)
        hapid[i + NC] = ch[hapsU[hapIdU[i]].code];
}

double BiCluster::getHapLP_step1(PARAMETER &para)
{
    double lp = 0;

    if ((int)para.Dloci.size() > 0)
    {
        int l;
        if (!allVectors)
        {
            if (useConstMarginPrior)
            {
                double cd = delta / pow((double)alleleCn, (double)para.Dloci.size() - 1.);
                double gcd = gsl_sf_lngamma(cd);
                for (l = 0; l < (int)para.haps.size(); l++)
                    lp += gsl_sf_lngamma(para.haps[l].hCn + cd) - gcd;
            }
            else
            {
                for (l = 0; l < (int)para.haps.size(); l++)
                    lp += logGammaH[para.haps[l].hCn];
                lp -= ((double)para.haps.size()) * gsldelta;
            }
        }
        else
        {
            if (para.Dsp)
            {
                if (useConstMarginPrior)
                {
                    double cd = delta / pow((double)alleleCn, (double)para.Dloci.size() - 1.);
                    double gcd = gsl_sf_lngamma(cd);
                    for (l = 0; l < (int)para.haps.size(); l++)
                        lp += gsl_sf_lngamma(para.haps[l].hCn + cd) - gcd;
                }
                else
                {
                    for (l = 0; l < (int)para.haps.size(); l++)
                        lp += logGammaH[para.haps[l].hCn];
                    lp -= ((double)para.haps.size()) * gsldelta;
                }
                //	for(l = 0; l < (int)para.hapsU.size(); l++)//
                //		lp += logGammaH[para.hapsU[l].hCn];//
                //	lp -= ((double)para.hapsU.size()) * gsldelta;//

                for (l = 0; l < (int)para.Dloci.size(); l++)
                    lp += singleDLP[para.Dloci[l]];
            }
            else
            {
                vector<int> thc, thid;
                combineHapCount(para.haps, para.hapsU, para.hapId, para.hapIdU, (int)para.Dloci.size(), thc, thid);
                if (useConstMarginPrior)
                {
                    double cd = delta / pow((double)alleleCn, (double)para.Dloci.size() - 1.);
                    double gcd = gsl_sf_lngamma(cd);
                    for (l = 0; l < (int)thc.size(); l++)
                        lp += gsl_sf_lngamma(thc[l] + cd) - gcd;
                }
                else
                {
                    for (l = 0; l < (int)thc.size(); l++)
                        lp += logGammaH[thc[l]];
                    lp -= ((double)thc.size()) * gsldelta;
                }
            }
        }
    }

    return lp;
}

//when using markov models, the interaction also needs to consider markov properties, otherwise likelihoods will be downtilled
//however, even if we consider Markov properties, it's still not working very well, still two modes in lnr, because the model itself is not accurate.
double BiCluster::getHapLP(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, PARAMETER &para)
{
    double lp = 0;
    int l;

    if ((int)para.Dloci.size() > 0)
    {
        vector<int> hapCn((int)para.haps.size());
        for (l = 0; l < (int)hapCn.size(); l++)
            hapCn[l] = para.haps[l].hCn;

        if (!allVectors)
            lp = getHapLP_sub(dataC, para.indState, para.Dloci, hapCn, para.hapId, 1);
        else
        {
            if (para.Dsp)
            {
                lp = getHapLP_sub(dataC, para.indState, para.Dloci, hapCn, para.hapId, 1);

                //	vector<int> hapCnU((int)para.hapsU.size());
                //	for(l = 0; l < (int)hapCnU.size(); l++)
                //		hapCnU[l] = para.hapsU[l].hCn;
                //	lp += getHapLP_sub(dataU, UMark, para.Dloci, hapCnU, para.hapIdU, 2);
                /////////debugging
                for (l = 0; l < (int)para.Dloci.size(); l++)
                    lp += getDLP_MC(dataC, dataU, para.indState, para.Dloci[l]);
                //lp += singleDLP[para.Dloci[l]];
            }
            else
            {
                vector<int> thc, thid;
                combineHapCount(para.haps, para.hapsU, para.hapId, para.hapIdU, (int)para.Dloci.size(), thc, thid);
                vector<vector<char> > data = dataC;
                data.insert(data.end(), dataU.begin(), dataU.end());
                vector<bool> tindState((int)data.size(), true);
                lp = getHapLP_sub(data, tindState, para.Dloci, thc, thid, 3);
            }
        }
    }

    return lp;
}

double BiCluster::getHapLP1(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, PARAMETER &para, bool model)
{
    double lp = 0;
    int l;

    if ((int)para.Dloci.size() > 0)
    {
        vector<int> hapCn((int)para.haps.size());
        for (l = 0; l < (int)hapCn.size(); l++)
            hapCn[l] = para.haps[l].hCn;

        {
            if (para.Dsp)
            {	//if(model)
                lp = getHapLP_sub(dataC, para.indState, para.Dloci, hapCn, para.hapId, 1);

                if (model)
                {
                    vector<int> hapCnU((int)para.hapsU.size());
                    for (l = 0; l < (int)hapCnU.size(); l++)
                        hapCnU[l] = para.hapsU[l].hCn;
                    lp += getHapLP_sub(dataU, UMark, para.Dloci, hapCnU, para.hapIdU, 2);
                }
                else
                {
                    for (l = 0; l < (int)para.Dloci.size(); l++)
                        lp += getDLP_MC(dataC, dataU, para.indState, para.Dloci[l]);
                }
            }
            else
            {
                if (model)
                {
                    vector<int> thc, thid;
                    combineHapCount(para.haps, para.hapsU, para.hapId, para.hapIdU, (int)para.Dloci.size(), thc, thid);
                    //	vector<vector<char> > data = dataC;
                    //	data.insert(data.end(), dataU.begin(), dataU.end());
                    //	vector<bool> tindState((int)data.size(), true);
                    //	lp = getHapLP_sub(data, tindState, para.Dloci, thc, thid, 3);
                    int lsz = (int)para.Dloci.size();
                    lp = 0;
                    if (useConstMarginPrior)
                    {
                        double cd = delta / pow((double)alleleCn, (double)para.Dloci.size() - 1.);
                        double gcd = gsl_sf_lngamma(cd);
                        for (int l = 0; l < (int)thc.size(); l++)
                            lp += gsl_sf_lngamma(thc[l] + cd) - gcd;
                        lp -= gsl_sf_lngamma((double)dataC.size() + (double)dataU.size() + Alphan) - gsl_sf_lngamma(Alphan);
                    }
                    else
                    {
                        for (int l = 0; l < (int)thc.size(); l++)
                            lp += logGammaH[thc[l]];
                        lp -= (double)thc.size() * gsldelta;
                        lp -= baseGamma[lsz][3];
                        lp += baseGamma[lsz][0];
                    }
                }
                else
                {
                    for (l = 0; l < (int)para.Dloci.size(); l++)
                        lp += getULP_MC(dataC, dataU, para.indState, para.Dloci[l], -1);
                }
            }
        }
    }

    return lp;
}

double BiCluster::getHapLP_sub(vector<vector<char> > const &data, vector<bool> const &indState, vector<int> const &Dloci,
                               vector<int> const &hapCn, vector<int> const &hapId, int datatype)
{
    int i, j, l, lsz = (int)Dloci.size();
    double rts, rtm;
    //int alpha = (int)(pow((double)alleleCn, (double)lsz) * delta + 0.5);
    int hN = (int)hapCn.size();

    vector<int> tloci = Dloci, nloci;
    sort(tloci.begin(), tloci.end());
    for (l = 0; l < lsz; l++)
    {
        if ((l == 0 && tloci[l] > 0) || (l > 0 && tloci[l] > tloci[l - 1] + 1))
            nloci.push_back(tloci[l] - 1);
    }
    if (mixture < 2 || (int)nloci.size() == 0)
    {
        rts = 0;
        if (useConstMarginPrior)
        {
            double cd = delta / pow((double)alleleCn, (double)Dloci.size() - 1.);
            double gcd = gsl_sf_lngamma(cd);
            for (l = 0; l < hN; l++)
                rts += gsl_sf_lngamma(hapCn[l] + cd) - gcd;
            rts -= gsl_sf_lngamma((double)data.size() + Alphan) - gsl_sf_lngamma(Alphan);
        }
        else
        {
            for (l = 0; l < hN; l++)
                rts += logGammaH[hapCn[l]];
            rts -= (double)hN * gsldelta;
            rts -= baseGamma[lsz][datatype];
            rts += baseGamma[lsz][0];
        }
    }
    if (mixture > 0 && (int)nloci.size() > 0)
    {
        vector<HAPLOTYPE> nhaps;
        vector<int> nhapId;
        getHapRep(data, nloci, indState, nhaps, nhapId);

        int N = (int)data.size();
        int nN = (int)nhaps.size();
        vector<vector<int> > newCn(nN, vector<int>(hN, 0));
        vector<vector<int> > pointer(nN, vector<int>(hN, -1));
        vector<int> cn(nN, 0);
        for (i = 0; i < N; i++)
            if (indState[i])
            {
                l = hapId[i];
                j = nhapId[i];
                if (pointer[j][l] < 0)
                {
                    newCn[j][cn[j]] = 0;
                    pointer[j][l] = cn[j];
                    cn[j]++;
                }
                newCn[j][pointer[j][l]]++;
            }
        for (i = 0; i < nN; i++)
            newCn[i].resize(cn[i]);

        rtm = 0;
        for (i = 0; i < nN; i++)
        {
            int sum = 0;
            if (useConstMarginPrior)
            {
                double cd = delta / pow((double)alleleCn, (double)Dloci.size() - 1.);
                double gcd = gsl_sf_lngamma(cd);
                for (j = 0; j < (int)newCn[i].size(); j++)
                {
                    rtm += gsl_sf_lngamma(newCn[i][j] + cd) - gcd;
                    sum += newCn[i][j];
                }
                rtm -= gsl_sf_lngamma(sum + Alphan) - gsl_sf_lngamma(Alphan);
            }
            else
            {
                for (j = 0; j < (int)newCn[i].size(); j++)
                {
                    rtm += logGammaH[newCn[i][j]];
                    sum += newCn[i][j];
                }
                rtm -= (double)newCn[i].size() * logGammaH[0];
                rtm -= logGamma[sum];
                rtm += logGamma[0];
            }
        }
    }

    if (mixture == 0 || (int)nloci.size() == 0) return rts;
    else if (mixture == 1)
    {
        double maxrt = rts;
        if (maxrt < rtm) maxrt = rtm;
        rts = exp(rtm - maxrt) + exp(rts - maxrt);
        return(log(rts / 2.) + maxrt);
    }
    else
    {
        return rtm;
    }
}

//note: prior of indicators are not accounted in here, it's seperately calculated and added to the total likelihood
double BiCluster::getHapLP_step2(PARAMETER &para)
{
    double lp = 0;

    if ((int)para.Dloci.size() > 0)
    {
        if (!allVectors)
        {
            if (useConstMarginPrior)
                lp -= gsl_sf_lngamma(mNC + Alphan) - gsl_sf_lngamma(Alphan);
            else
            {
                lp -= baseGamma[(int)para.Dloci.size()][1];
                lp += baseGamma[(int)para.Dloci.size()][0];
            }
        }
        else
        {
            if (para.Dsp)
            {
                if (useConstMarginPrior)
                    lp -= gsl_sf_lngamma(mNC + Alphan) - gsl_sf_lngamma(Alphan);
                else
                {
                    lp -= baseGamma[(int)para.Dloci.size()][1];
                    lp += baseGamma[(int)para.Dloci.size()][0];
                }
                //	lp -= baseGamma[(int)para.Dloci.size()][2];//
                //	lp += baseGamma[(int)para.Dloci.size()][0];//
            }
            else
            {
                if (useConstMarginPrior)
                    lp -= gsl_sf_lngamma(mNC + mNU + Alphan) - gsl_sf_lngamma(Alphan);
                else
                {
                    lp -= baseGamma[(int)para.Dloci.size()][3];
                    lp += baseGamma[(int)para.Dloci.size()][0];
                }
            }
        }
    }

    return lp;
}

void BiCluster::updateDsp(PARAMETER &para)
{
    bool odsp = para.Dsp;
    para.Dsp = true;
    double a = getHapLP_step1(para) + getHapLP_step2(para);
    para.Dsp = false;
    double b = getHapLP_step1(para) + getHapLP_step2(para);
    double oa = a, ob = b;

    if (a > b + 30.) para.Dsp = true;
    else if (a < b - 30.) para.Dsp = false;
    else
    {
        a = exp(a - b) / (exp(a - b) + 1.);
        double un = gsl_ran_flat(gammar, 0, 1.);
        if (un < a) para.Dsp = true;
        else para.Dsp = false;
    }
    if (para.Dsp && !odsp)
        para.logP += oa - ob;
    else if (!para.Dsp && odsp)
        para.logP += ob - oa;
}

void BiCluster::initializeGamma(int NC, int NU, int L, int alleleCn, int locilimit)
{
    int i;
    double p;

    delta = Alphan / (double)alleleCn;
    gsldelta = gsl_sf_lngamma(delta);

    p = gsl_sf_lngamma(Alphan);
    logGamma.resize(NC + L + NU + (int)pow((double)alleleCn, (double)locilimit) + alleleCn - 1 + (int)(Alphan + 1.));
    for (i = 0; i < (int)logGamma.size(); i++)
    {
        logGamma[i] = p;
        p += log((double)i + Alphan);
    }
    logGammaS.resize(NC + L + NU + (int)pow((double)alleleCn, (double)locilimit) + alleleCn - 1 + (int)(Alphan + 1.));
    p = gsldelta;
    for (i = 0; i < (int)logGammaS.size(); i++)
    {
        logGammaS[i] = p;
        p += log((double)(i + delta));
    }
    logGammaH = logGammaS;

    baseGamma.resize(locilimit + 1, vector<double>(4, 0));
    int k = alleleCn;
    for (i = 1; i < (int)baseGamma.size(); i++)
    {
        baseGamma[i][0] = gsl_sf_lngamma(delta * (double)k);
        baseGamma[i][1] = gsl_sf_lngamma((double)NC + delta * (double)k);
        baseGamma[i][2] = gsl_sf_lngamma((double)NU + delta * (double)k);
        baseGamma[i][3] = gsl_sf_lngamma((double)NC + (double)NU + delta * (double)k);
        k *= alleleCn;
    }

    UMark.clear();
    UMark.resize(NU, true);
}

//before using this function, some initialization is required! such as logGamma[], SCcounts/MCcounts, and singleDLP
double BiCluster::computeLnr(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> const &Dloci)
{
    int NC = (int)dataC.size();
    int NU = (int)dataU.size();
//	int L = (int)dataC[0].size();

    int om = mixture;
    mixture = 0;

    PARAMETER para;
    para.indState.resize(NC, true);
    if ((int)UMark.size() < NU)
    {
        UMark.clear();
        UMark.resize(NU, true);
    }
    para.Dcn = NC;
//	para.lociStatus.resize(L, 0);
    para.Dsp = true;

    double lp = 0;
    PARAMETER newpara = para;
    int j, k;

    if ((int)Dloci.size() == 1)
    {
        k = Dloci[0];
        lp -= getULP_MC(dataC, dataU, para.indState, k, k - 1);
        lp += getDLP_MC(dataC, dataU, para.indState, k);
        lp += getSLP_MC(dataC, para.indState, k);
    }
    else
    {
        newpara.Dloci = Dloci;
        for (j = 0; j < (int)Dloci.size(); j++)
        {
            k = Dloci[j];
            incHapRep(dataC, k, newpara.indState, newpara.haps, newpara.hapId);
            if (allVectors)
                incHapRep(dataU, k, UMark, newpara.hapsU, newpara.hapIdU);
            lp -= getULP_MC(dataC, dataU, newpara.indState, k, k - 1);
            if (!allVectors)
                lp += getDLP_MC(dataC, dataU, newpara.indState, k);
        }
        newpara.Dsp = true;
        double a = getHapLP1(dataC, dataU, newpara, false);
        double b = getHapLP1(dataC, dataU, newpara, true);
        newpara.Dsp = false;
        double c = getHapLP1(dataC, dataU, newpara, false);
        double d = getHapLP1(dataC, dataU, newpara, true);
        double max1 = a>b?a:b;
        double max2 = c>d?c:d;
        lp = log((exp(a-max1)+exp(b-max1))/2.) + max1;
        lp -= log((exp(c-max2)+exp(d-max2))/2.) + max2;
    }

    mixture = om;

    return lp;
}

//call this from outside
int BiCluster::runChains(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, double P, vector<OUTPUT> &output,
                         vector<int> &totalcounts, vector<double> &group1, vector<double> &group2, int indRuns, vector<vector<int> > const &positions, bool multiTry, bool singleOnly, int validationmethod,
                         int startTries, int tryLen)
{
    int L = (int)dataC[0].size();
    int i, j, k;
    vector<vector<int> > Ssamples, Isamples;

    totalcounts.clear();
    totalcounts.resize(2 * L, 0);
    group1.clear();
    group2.clear();
    int samplesize = 0;

    vector<vector<int> > configures;
    vector<double> lps;
    if (startTries > 0)
        randomStarts(dataC, dataU, positions, singleOnly, startTries, tryLen, configures, lps);
    for (i = 0; i < indRuns; i++)
    {
        vector<int> ts, ti;
        vector<double> g1, g2;
        samplesize += parallelTempering_MC(dataC, dataU, positions, ts, ti, g1, g2, multiTry, (i == 0), singleOnly, i + 1, configures, lps);

        for (k = 0; k < L; k++)
        {
            totalcounts[k] += ts[k];
            totalcounts[k + L] += ti[k];
        }
        for (k = 0; k < (int)g1.size(); k++)
        {
            if (k < (int)group1.size())
                group1[k] += g1[k];
            else group1.push_back(g1[k]);
        }
        for (k = 0; k < (int)g2.size(); k++)
        {
            if (k < (int)group2.size())
                group2[k] += g2[k];
            else group2.push_back(g2[k]);
        }
    }
    vector<int> sc = totalcounts;
    vector<int> ic = totalcounts;
    sc.resize(L);
    ic.erase(ic.begin(), ic.begin() + L);
    getAssociations(dataC, dataU, sc, ic, positions, P, validationmethod, output);

    for (i = 0; i < (int)group1.size(); i++)
        group1[i] /= (double)indRuns;
    for (i = 0; i < (int)group2.size(); i++)
        group2[i] /= (double)indRuns;

    return samplesize;
}

//used by outside calls where only the postior counts of single markers (interactions assigned to single markers) are available
void BiCluster::getAssociations(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU,
                                vector<int> const &Ssamples, vector<int> const &Isamples, vector<vector<int> > const &positions,
                                double Pvalue, int validationmethod, vector<OUTPUT> &output)
{
    int i, j, k;
    int L = (int)dataC[0].size();

    vector<vector<vector<int> > > markers(3);
    output.clear();
    vector<OUTPUT> tmpoutput;
    //collect significant outputs
    for (i = 1; i <= (int)markers.size(); i++)
    {
        vector<MySortTypeb> counts(L);
        for (j = 0; j < L; j++)
        {
            if (i == 1) counts[j].value = (double)Ssamples[j];
            else counts[j].value = (double)Isamples[j];
            counts[j].indices[0] = j;
        }
        sort(counts.begin(), counts.end());
        if (i == 1)
        {
            for (j = (int)counts.size() - 1; j > 0; j--)
            {
                markers[i - 1].push_back(vector<int>(1, counts[j].indices[0]));
                if (counts[j].value < 1 || counts[j] < counts[(int)counts.size() * 4 / 5]) break; //only check top ones
            }
        }
        if (i > 1)
        {	//first, select a candidate set of interacting markers
            vector<int> selected;
            for (j = (int)counts.size() - 1; j >= (int)counts.size() * 4 / 5; j--)
            {
                if (counts[j].value < 1) break; ////////
                selected.push_back(counts[j].indices[0]);
            }
            sort(selected.begin(), selected.end());
            reverse(selected.begin(), selected.end());
            int sz = (int)selected.size();
            //second, based on these set of markers, we obtain all combinations of them
            if (sz >= i)
            {
                vector<int> pt(i);
                for (j = 0; j < i; j++)
                    pt[j] = j;
                bool a_inner = true;
                do{
                    bool valid = true;
                    if ((int)positions.size() == L)
                    {
                        for (k = 0; k < i - 1; k++)
                        {
                            int t;
                            for (t = k + 1; t < i; t++)
                            {
                                if (_invalidInt(selected[pt[k]], selected[pt[t]], L))
                                    break;
                            }
                            if (t < i) break;
                        }
                        if (k < i - 1) valid = false;
                    }
                    if (valid)
                    {
                        vector<int> a(i);
                        for (k = 0; k < i; k++)
                            a[k] = selected[pt[k]];
                        markers[i - 1].push_back(a);
                    }

                    k = i - 1;
                    pt[k]++;
                    while (pt[k] > sz - i + k)
                    {
                        k--;
                        if (k >= 0)
                        {
                            pt[k]++;
                            for (j = k + 1; j < i; j++)
                                pt[j] = pt[k] + j - k;
                        }
                        else
                        {
                            a_inner = false;
                            break;
                        }
                    }
                } while (a_inner);
            }
        }
    }
    hierarchicalPower(dataC, dataU, Pvalue, validationmethod, markers, output);
}

bool BiCluster::updateLoci_S_Multiple(vector<vector<char> > const &data, vector<vector<char> > const &dataU,
                                      PARAMETER &para, double p0, double p1, double T)
{	//int N = (int)data.size();
    int L = (int)data[0].size();
    int l;
    PARAMETER newpara;
    int bound = L - (int)para.Dloci.size();
    double p2 = 1. - p0 - p1;

    double logratio = 0;
    int movetype;
    int ad, rm;
    double lp = 0, p;

    ad = rm = -1;
    if ((int)para.DSloci.size() == 0)
        p = 1.;
    else if ((int)para.DSloci.size() >= bound)
        p = gsl_ran_flat(gammar, 0, p0 + p1);
    else p = gsl_ran_flat(gammar, 0, 1);
    if (p < p0)
    {
        movetype = 0;
        if ((int)para.DSloci.size() == 1)
            logratio = log(1. / p0);
        else if ((int)para.DSloci.size() >= bound)
            logratio = log(p2);
        else logratio = log(p2 / p0);
        logratio += - log((double)L - (double)para.Dloci.size() - (double)para.DSloci.size() + 1.);
    }
    else if (p < p0 + p1)
    {
        bool a_inner;
        do {
            ad = (int)gsl_ran_flat(gammar, 0, L);
            a_inner = true;
            if (para.lociStatus[ad] == 0)
                a_inner = false;
        } while (a_inner);
        movetype = 1;
    }
    else
    {
        bool a_inner;
        do {
            ad = (int)gsl_ran_flat(gammar, 0, L);
            a_inner = true;
            if (para.lociStatus[ad] == 0)
                a_inner = false;
        } while (a_inner);
        movetype = 2;
        if ((int)para.DSloci.size() == 0)
            logratio = log(p0);
        else if ((int)para.DSloci.size() == bound - 1)
            logratio = log(1. / p2);
        else logratio = log(p0 / p2);
        logratio += log((double)L - (double)para.Dloci.size() - (double)para.DSloci.size());
    }

    //int i;
    if (movetype == 1)
    {
        vector<double> lps((int)para.DSloci.size(), 0);
        double lpmax;
        for (l = 0; l < (int)para.DSloci.size(); l++)
        {
            newpara = para;

            lps[l] -= getDLP_MC(data, dataU, para.indState, para.DSloci[l]);
            lps[l] -= getSLP_MC(data, para.indState, para.DSloci[l]);
            lps[l] += getULP_MC(data, dataU, newpara.indState, para.DSloci[l], para.DSloci[l] - 1);

            lps[l] -= getULP_MC(data, dataU, newpara.indState, ad, ad - 1);
            lps[l] += getDLP_MC(data, dataU, newpara.indState, ad);
            lps[l] += getSLP_MC(data, newpara.indState, ad);

            lps[l] /= T;
            if (l == 0 || lpmax < lps[l]) lpmax = lps[l];
        }

        vector<double> p = lps;
        vector<double> psum = p;
        for (l = 0; l < (int)p.size(); l++)
        {
            if (p[l] - lpmax > -20.)
                p[l] = exp(p[l] - lpmax);
            else p[l] = 0;
            if (l == 0) psum[l] = p[l];
            else psum[l] = p[l] + psum[l - 1];
        }
        double sum = *(psum.end() - 1);
        double un = gsl_ran_flat(gammar, 0, sum);
        for (l = 0; l < (int)psum.size() - 1; l++)
            if (psum[l] > un)
                break;
        rm = para.DSloci[l];
        lp = lps[l] * T;
        double a0 = exp(-lpmax);
        if (p[l] > 0 && sum - p[l] + a0 > 0)
        {
            if (a0 * sum > 0)
                logratio += log(a0 * sum / p[l] / (sum - p[l] + a0)); //debug
            else logratio += -1000000;
        }
        newpara = para;
        newpara.DSloci.erase(newpara.DSloci.begin() + l, newpara.DSloci.begin() + l + 1);
        newpara.nlSCn.erase(newpara.nlSCn.begin() + l, newpara.nlSCn.begin() + l + 1);

        newpara.DSloci.push_back(ad);
        vector<int> tc(alleleCn, 0);
        newpara.nlSCn.push_back(tc);
    }
    else if (movetype == 0) //genetic -> phenocopy
    {
        vector<double> lps((int)para.DSloci.size(), 0);
        double lpmax;
        for (l = 0; l < (int)para.DSloci.size(); l++)
        {
            newpara = para;

            lps[l] -= getDLP_MC(data, dataU, para.indState, para.DSloci[l]);
            lps[l] -= getSLP_MC(data, para.indState, para.DSloci[l]);
            lps[l] += getULP_MC(data, dataU, newpara.indState, para.DSloci[l], para.DSloci[l] - 1);

            lps[l] /= T;
            if (l == 0 || lpmax < lps[l]) lpmax = lps[l];
        }

        vector<double> p = lps;
        vector<double> psum = p;
        for (l = 0; l < (int)p.size(); l++)
        {
            if (p[l] - lpmax > -20.)
                p[l] = exp(p[l] - lpmax);
            else p[l] = 0;
            if (l == 0) psum[l] = p[l];
            else psum[l] = p[l] + psum[l - 1];
        }
        double sum = *(psum.end() - 1);
        double un = gsl_ran_flat(gammar, 0, sum);
        for (l = 0; l < (int)psum.size() - 1; l++)
            if (psum[l] > un)
                break;
        rm = para.DSloci[l];
        lp = lps[l] * T;
        if (sum > 0)
        {
            if (p[l] > 0)
                logratio += log(sum / p[l]);
            else logratio += 1000000;
        }

        newpara = para;
        newpara.DSloci.erase(newpara.DSloci.begin() + l, newpara.DSloci.begin() + l + 1);
        newpara.nlSCn.erase(newpara.nlSCn.begin() + l, newpara.nlSCn.begin() + l + 1);
    }
    else //phenocopy -> genetic
    {
        newpara = para;

        newpara.DSloci.push_back(ad);
        vector<int> tc(alleleCn, 0);
        newpara.nlSCn.push_back(tc);
        lp -= getULP_MC(data, dataU, para.indState, ad, ad - 1);
        lp += getDLP_MC(data, dataU, para.indState, ad);
        lp += getSLP_MC(data, para.indState, ad);

        vector<double> lps((int)newpara.DSloci.size(), 0);
        double lpmax;
        for (l = 0; l < (int)newpara.DSloci.size(); l++)
        {
            if (l < (int)newpara.DSloci.size() - 1)
            {
                PARAMETER backpara = newpara;

                lps[l] -= getDLP_MC(data, dataU, newpara.indState, newpara.DSloci[l]);
                lps[l] -= getSLP_MC(data, newpara.indState, newpara.DSloci[l]);
                lps[l] += getULP_MC(data, dataU, backpara.indState, newpara.DSloci[l], newpara.DSloci[l] - 1);
            }
            else lps[l] = -lp;
            lps[l] /= T;
            if (l == 0 || lpmax < lps[l]) lpmax = lps[l];
        }

        double sum = 0;
        for (l = 0; l < (int)lps.size(); l++)
        {
            if (lps[l] - lpmax > -20.)
                lps[l] = exp(lps[l] - lpmax);
            else lps[l] = 0;
            sum += lps[l];
        }
        if (sum > 0)
        {
            if (*(lps.end() - 1) > 0)
                logratio += log(*(lps.end() - 1) / sum);
            else logratio += -1000000;
        }
    }
//-------------------------------------
    if (movetype != 1)
    {
        if (movetype == 0)
            lp += logprop[0] - logprop[1];
        else if (movetype == 2)
            lp += logprop[1] - logprop[0];
    }
    newpara.logP += lp;

    double un = gsl_ran_flat(gammar, 0, 1);

    if (lp / T + logratio >= 0 || exp((lp / T + logratio)) >= un)
    {
        para = newpara;
        if (movetype > 0) para.lociStatus[ad] = 1;
        if (movetype < 2) para.lociStatus[rm] = 0;
        return true;
    }
    else return false;
}

void BiCluster::logicR_counts(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU,
                              vector<int> const &dloci, vector<int> &countC, vector<int> &countU)
{
    int i, j, k;
    int sz = (int)pow((double)alleleCn, (double)dloci.size());
    countC.clear();
    countU.clear();
    countC.resize(sz + 1, 0);
    countU.resize(sz + 1, 0);

    for (i = 0; i < (int)dataC.size(); i++)
    {
        k = 0;
        for (j = 0; j < (int)dloci.size(); j++)
            k = k * alleleCn + (int)dataC[i][dloci[j]];
        countC[k]++;
    }
    for (i = 0; i < (int)dataU.size(); i++)
    {
        k = 0;
        for (j = 0; j < (int)dloci.size(); j++)
            k = k * alleleCn + (int)dataU[i][dloci[j]];
        countU[k]++;
    }
    countC[sz] = (int)dataC.size();
    countU[sz] = (int)dataU.size();
}

void BiCluster::logicR_estPara(vector<int> const &countC, vector<int> const &countU,
                               vector<double> &coeff, vector<double> &lps)
{
    int sz = (int)countC.size();
    coeff.resize(sz);
    int i = 0;
    for (i = 0; i < sz; i++)
    {
        if (countC[i] == 0 && countU[i] > 0) coeff[i] = -1000;
        else if (countC[i] > 0 && countU[i] == 0) coeff[i] = 1000;
        else if (countC[i] == 0 && countU[i] == 0) coeff[i] = 0;
        else coeff[i] = log((double)countC[i] / (double)countU[i]);
    }
    lps.resize(2 * sz);
    for (i = 0; i < sz; i++)
    {
        lps[i] = - log(1. + exp(-coeff[i]));
        if (coeff[i] <= -1000) lps[i] = -1000.;
        lps[i + sz] = - log(1. + exp(coeff[i]));
        if (coeff[i] >= 1000) lps[i + sz] = -1000;
    }
}

double BiCluster::logicR_logLikelihood(vector<int> const &countC, vector<int> const &countU, vector<double> const &lps)
{
    double LLR = 0;
    int i;
    int sz = (int)countC.size();

    for (i = 0; i < sz - 1; i++)
        LLR += (double)countC[i] * lps[i];
    for (i = 0; i < sz - 1; i++)
        LLR += (double)countU[i] * lps[i + sz];
    return LLR;
}

void BiCluster::associationScore_test(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU,
                                      int degree, char *filename)
{
    int L = (int)dataC[0].size();
    int NC = (int)dataC.size();
    int NU = (int)dataU.size();
    vector<int> dloci(degree);
    int i;

    int sz = (int)pow((double)alleleCn, (double)degree);
    bool a_inner;

//	statistics.clear();
    for (i = 0; i < degree; i++)
        dloci[i] = i;
    a_inner = true;
    FILE *f = fopen(filename, "w");
    fclose(f);
    vector<double> statistics(1000000);
    int ct = 0;
    do {
        double llr = computeLnr(dataC, dataU, dloci);
        //statistics.push_back(llr);
        statistics[ct++] = llr;
        if (ct >= 1000000)
        {
            f = fopen(filename, "a");
            for (int k = 0; k < ct; k++)
                fprintf(f, "%f\n", statistics[k]);
            fclose(f);
            ct = 0;
        }

        i = degree;
        do {
            i--;
            dloci[i]++;
        } while (i > 0 && dloci[i] > L - degree + i);
        if (i == 0 && dloci[i] > L - degree)
            a_inner = false;
        else
        {
            for (i = i + 1; i < degree; i++)
                dloci[i] = dloci[i - 1] + 1;
        }
    } while (a_inner);

    {
        f = fopen(filename, "a");
        for (int k = 0; k < ct; k++)
            fprintf(f, "%f\n", statistics[k]);
        fclose(f);
    }
}

double BiCluster::computeLnr(vector<vector<char> > const &dataC,
                             vector<vector<char> > const &dataU,
                             PARAMETER &para, vector<int> const &Dloci)
{
    int om = mixture;
    mixture = 0;

    int NC = (int)dataC.size();
    int NU = (int)dataU.size();

    double lp = 0;
    int k;

    if ((int)Dloci.size() == 1)
    {
        k = Dloci[0];
        lp -= getULP_MC(dataC, dataU, para.indState, k, k - 1);
        lp += getDLP_MC(dataC, dataU, para.indState, k);
        lp += getSLP_MC(dataC, para.indState, k);
    }
    else
    {
        para.Dloci = Dloci;

        para.Dsp = true;
        double a = getHapLP1(dataC, dataU, para, false);
        double b = getHapLP1(dataC, dataU, para, true);
        para.Dsp = false;
        double c = getHapLP1(dataC, dataU, para, false);
        double d = getHapLP1(dataC, dataU, para, true);

        double max1 = a>b?a:b;
        double max2 = c>d?c:d;
        lp = log((exp(a-max1)+exp(b-max1))/2.) + max1;
        lp -= log((exp(c-max2)+exp(d-max2))/2.) + max2;
    }
    mixture = om;

    return lp;
}

//conditional on the masked loci, which are significant ones
double BiCluster::computeConditionalLnr(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU,
                                        vector<int> const &Dloci, vector<bool> const &mask)
{
    int om = mixture;
    mixture = 0;

    int NC = (int)dataC.size();
    int NU = (int)dataU.size();

    double lp = 0;
    int j, k;

    PARAMETER newpara;
    newpara.indState.resize(NC, true);
    if ((int)UMark.size() < NU)
    {
        UMark.clear();
        UMark.resize(NU, true);
    }
    newpara.Dcn = NC;

    vector<int> tloci, rloci = Dloci;
    for (int i = 0; i < (int)mask.size(); i++)
        if (mask[i] == true)
        {
            tloci.push_back(Dloci[i]);
            for (int j = 0; j < (int)rloci.size(); j++)
                if (rloci[j] == Dloci[i])
                {
                    rloci.erase(rloci.begin() + j);
                    break;
                }
        }

    if ((int)Dloci.size() == 1 && (int)tloci.size() == 0)
    {
        k = Dloci[0];
        lp -= getULP_MC(dataC, dataU, newpara.indState, k, k - 1);
        lp += getDLP_MC(dataC, dataU, newpara.indState, k);
        lp += getSLP_MC(dataC, newpara.indState, k);
    }
    else if ((int)rloci.size() > 0)
    {
        newpara.Dloci = tloci;
        for (j = 0; j < (int)tloci.size(); j++)
        {
            k = tloci[j];
            incHapRep(dataC, k, newpara.indState, newpara.haps, newpara.hapId);
            if (allVectors)
                incHapRep(dataU, k, UMark, newpara.hapsU, newpara.hapIdU);
        }
        newpara.Dsp = true;
        double a = -getHapLP1(dataC, dataU, newpara, false);
        double b = -getHapLP1(dataC, dataU, newpara, true);
        newpara.Dsp = false;
        double c = -getHapLP1(dataC, dataU, newpara, true);
        double d = -getHapLP1(dataC, dataU, newpara, false);

        newpara.Dloci.insert(newpara.Dloci.end(), rloci.begin(), rloci.end());
        for (j = 0; j < (int)rloci.size(); j++)
        {
            k = rloci[j];
            incHapRep(dataC, k, newpara.indState, newpara.haps, newpara.hapId);
            if (allVectors)
                incHapRep(dataU, k, UMark, newpara.hapsU, newpara.hapIdU);
        }
        newpara.Dsp = true;
        a += getHapLP1(dataC, dataU, newpara, false);
        b += getHapLP1(dataC, dataU, newpara, true);
        newpara.Dsp = false;
        c += getHapLP1(dataC, dataU, newpara, true);
        d += getHapLP1(dataC, dataU, newpara, false);

        double max1 = a>b?a:b;
        double max2 = c>d?c:d;
        lp = log((exp(a-max1)+exp(b-max1))/2.) + max1;
        lp -= log((exp(c-max2)+exp(d-max2))/2.) + max2;
    }
    mixture = om;

    return lp;
}

void BiCluster::computeConstant1(int Nd, int Nu, double alpha, int intsz, vector<double> &results)
{
    int i;
    results.resize(intsz * 2);
    for (i = 1; i <= intsz; i++)
    {
        int d = (int)pow((double)alleleCn, i);
        double k = (double)(d - 1) / 2.;
        results[i - 1] = k * (log((double)Nd) + log((double)Nu) - log((double)Nd + Nu) - log(2. * 3.141593));
        results[i - 1] -= gsl_sf_lngamma(alpha * (double)d) - (double)d * gsl_sf_lngamma(alpha);
        results[intsz + i - 1] = results[i - 1]; //for LD markers
        results[i - 1] -= log((double)Nu / (double)(Nd + Nu)) * (double)(d - (alleleCn - 1) * i - 1) / 2.;//////for non-LD markers
    }
}

double BiCluster::conditionalLLR(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> const &Dloci, vector<bool> const &mask)
{
    vector<double> coeff, lps;
    vector<int> countC, countU;
    int NC = (int)dataC.size(), NU = (int)dataU.size();
    int i, j;

    vector<int> tloci, rloci;
    for (i = 0; i < (int)mask.size(); i++)
        if (mask[i]) tloci.push_back(Dloci[i]);
    for (i = 0; i < (int)Dloci.size(); i++)
    {
        for (j = 0; j < (int)tloci.size(); j++)
            if (tloci[j] == Dloci[i])
                break;
        if (j >= (int)tloci.size())
            rloci.push_back(Dloci[i]);
    }
    if (rloci.size() == 0) return 0;
    int tsz = (int)pow((double)alleleCn, (double)tloci.size());
    int rsz = (int)pow((double)alleleCn, (double)rloci.size());
    double lnrl = 0;
    vector<int> dloci = tloci;
    dloci.insert(dloci.end(), rloci.begin(), rloci.end());
    logicR_counts(dataC, dataU, dloci, countC, countU);
    for (i = 0; i < tsz; i++)
    {
        vector<int> cc(rsz + 1,0), cu(rsz + 1,0);
        for (j = 0; j < rsz; j++)
        {
            cc[j] = countC[i * rsz + j];
            cu[j] = countU[i * rsz + j];
            cc[rsz] += cc[j];
            cu[rsz] += cu[j];
        }
        logicR_estPara(cc, cu, coeff, lps);
        double ll1 = logicR_logLikelihood(cc, cu, lps);
        double ll2 = (double)cc[rsz] * lps[rsz] + (double)cu[rsz] * lps[rsz + rsz + 1];
        lnrl += 2 * (ll1 - ll2);
    }
    return lnrl;
}

//method: 0: B-stat, 1: chi-square, 2: logistic reg
void BiCluster::hierarchicalPower(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, double pvalue, int method,
                                  vector<vector<vector<int> > > &markers, vector<OUTPUT> &output)
{
    int N = (int)dataC.size();
    int NU = (int)dataU.size();
    int L = (int)dataC[0].size();
    vector<double> offset;
    int i, j, k;

    int maxdegree = (int)markers.size();
    computeConstant1((int)dataC.size(), (int)dataU.size(), Alphan / (double)alleleCn, maxdegree, offset);
    initializeGamma(N, NU, L, alleleCn, 3);
    computeAlleleMC(dataC, dataU);
    singleDLP.clear();
    singleDLP.resize(L);
    for (i = 0; i < L; i++)
        singleDLP[i] = getDLP_MC(dataC, dataU, vector<bool>(), i);

    vector<int> found;
    output.clear();

    for (i = 1; i <= maxdegree; i++)
    {
        vector<int> loci(i);

        for (int t = 0; t < (int)markers[i - 1].size(); t++)
        {
            loci = markers[i - 1][t];
            double score;
            if (method == 0)	score = computeLnr(dataC, dataU, loci);
            else if (method == 1) score = chiSquare(dataC, dataU, loci);
            else score = LLR(dataC, dataU, loci);

            int df = (int)pow((double)alleleCn, (double)i) - 1;
            vector<int> taken;
            bool independence = true;
            if (method == 0)
            {
                independence = testIndependence_full_LR(dataU, loci, vector<int>(), 0.01, 1);
                if (i == 1 || independence)
                    score += offset[i - 1]; ///////////////////////////////////////////////////
                else score += offset[maxdegree + i - 1];
            }
            for (j = 0; j < (int)found.size(); j++)
            {
                for (k = 0; k < i; k++)
                    if (loci[k] == found[j])
                        break;
                if (k < i) taken.push_back(loci[k]);
            }
            if (taken.size() > 0)
            {
                if (method == 1) score -= chiSquare(dataC, dataU, taken);
                else
                {
                    vector<bool> mask((int)loci.size(), false);
                    for (j = 0; j < (int)taken.size(); j++)
                    {
                        for (k = 0; k < (int)loci.size(); k++)
                            if (loci[k] == taken[j])
                            {
                                mask[k] = true;
                                break;
                            }
                    }
                    if (method == 0)
                    {
                        score = computeConditionalLnr(dataC, dataU, loci, mask);
                        k = (int)taken.size();
                        if (i == 1 || independence)
                            score += offset[(int)loci.size() - 1] - offset[(int)taken.size() - 1]; /////////////////////////////////////////////////////
                        else score += offset[maxdegree + (int)loci.size() - 1] - offset[maxdegree + (int)taken.size() - 1];
                        //	score += log((double)dataU.size() / (double)(dataC.size() + dataU.size())) * (pow(3., (double)k) - 2. * k - 1.) / 2.; /////////////////////////////
                    }
                    else if (method == 2)
                        score = conditionalLLR(dataC, dataU, loci, mask);
                }
                df -= (int)pow((double)alleleCn, (double)taken.size())-1;
            }

            if (method == 0) score = score * 2;
            double p = 1 - gsl_cdf_chisq_P(score, df);
            if (df <= 0) p = 1;
            for (k = 0; k < i; k++)
                p *= (double)(L - k) / (k + 1);
            if (p > 1) p = 1;

            FILE *f = fopen(Posteriorfile, "a");
            vector<int> untt, tt = taken;
            for (k = 0; k < (int)loci.size(); k++)
            {
                int t = 0;
                for (t = 0; t < (int)taken.size(); t++)
                    if (taken[t] == loci[k]) break;
                if (t >= (int)taken.size()) untt.push_back(loci[k]);
            }
            sort(tt.begin(), tt.end());
            sort(untt.begin(), untt.end());
            for (k = 0; k < (int)tt.size(); k++)
            {
                if (k == 0) fprintf(f, "(");
                fprintf(f, "%d", tt[k]);
                if (k < (int)tt.size() - 1) fprintf(f, ", ");
                else fprintf(f, "), ");
            }
            for (k = 0; k < (int)untt.size(); k++)
            {
                fprintf(f, "%d", untt[k]);
                if (k < (int)untt.size() - 1) fprintf(f, ", ");
                else fprintf(f, " :\t");
            }
            fprintf(f, " B-stat = %f\t P = %f\n", score, p);
            fclose(f);

            if (p <= pvalue)
            {
                OUTPUT rt;
                rt.markers = loci;
                rt.pvalue = p;
                output.push_back(rt);
                int x, y;
                for (x = 0; x < (int)loci.size(); x++)
                {
                    for (y = 0; y < (int)found.size(); y++)
                        if (found[y] == loci[x]) break;
                    if (y >= (int)found.size())
                        found.push_back(loci[x]);
                }
            }
        }
    }
}

double BiCluster::chiSquare(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> &dloci)
{
    int i, j, k;
    int NC = (int)dataC.size();
    int NU = (int)dataU.size();

    int hsz = (int)pow((double)alleleCn, (double)dloci.size());
    double stat = 0;
    vector<int> cn(hsz, 0), un(hsz, 0);
    for (k = 0; k < NC; k++)
    {
        j = 0;
        for (i = 0; i < (int)dloci.size(); i++)
            j = j * alleleCn + dataC[k][dloci[i]];
        cn[j]++;
    }
    for (k = 0; k < NU; k++)
    {
        j = 0;
        for (i = 0; i < (int)dloci.size(); i++)
            j = j * alleleCn + dataU[k][dloci[i]];
        un[j]++;
    }

    double w = (double)NC / (double)(NC + NU);
    for (k = 0; k < alleleCn; k++)
        if (cn[k] + un[k] > 0) stat += pow(((double)cn[k] - w * (double)(cn[k] + un[k])), 2.) / (w * (double)(cn[k] + un[k]));
    for (k = 0; k < alleleCn; k++)
        if (cn[k] + un[k] > 0) stat += pow(((double)un[k] - (1. - w) * (double)(cn[k] + un[k])), 2.) / ((1. - w) * (double)(cn[k] + un[k]));

    return stat;
}

double BiCluster::LLR(vector<vector<char> > const &dataC, vector<vector<char> > const &dataU, vector<int> const &dloci)
{
    vector<double> coeff, lps;
    vector<int> countC, countU;
    int NC = (int)dataC.size(), NU = (int)dataU.size();

    int sz = (int)pow((double)alleleCn, (double)dloci.size());
    logicR_counts(dataC, dataU, dloci, countC, countU);
    logicR_estPara(countC, countU, coeff, lps);
    double ll1 = logicR_logLikelihood(countC, countU, lps);
    double ll2 = (double)NC * lps[sz] + (double)NU * lps[sz + sz + 1];

    return 2 * (ll1 - ll2);
}

bool BiCluster::testIndependence_full_LR(vector<vector<char> > const &dataU, vector<int> const &loci, vector<int> const &taken, double P, int radius)
{
    int i, j, k, l, m;
    if ((int)loci.size() <= (int)taken.size() + 1) return true;
    vector<int> rloci;
    int L = (int)dataU[0].size();

    for (i = 0; i < (int)loci.size(); i++)
    {
        for (j = 0; j < (int)taken.size(); j++)
            if (loci[i] == taken[j])
                break;
        if (j >= (int)taken.size()) rloci.push_back(loci[i]);
    }
    double rsz = (double)rloci.size();
    P /= (rsz * (rsz - 1) / 2.);

    bool ind = true;
    vector<int> loci1(1), loci2(1);
    for (j = 0; j < (int)rloci.size() - 1; j++)
        for (k = j + 1; k < (int)rloci.size(); k++)
        {
            double tmle = 0, tdf = 0;
            for (l = -radius; l <= radius; l++)
                for (m = -radius; m <= radius; m++)
                {
                    loci1[0] = rloci[j] + l;
                    loci2[0] = rloci[k] + m;
                    if (loci1[0] < 0 || loci1[0] >= L || loci2[0] < 0 || loci2[0] >= L) continue;
                    if ((rloci[j] < rloci[k] && loci1[0] < loci2[0]) || (rloci[j] > rloci[k] && loci1[0] > loci2[0]))
                    {
                        double mle, df;
                        testIndependence_twogroups_LR(dataU, loci1, loci2, taken, mle, df);
                        tmle += mle;
                        tdf += df;
                    }
                }
            if (1. - gsl_cdf_chisq_P(tmle, tdf) < P) return false;
        }

    return true;
}

double BiCluster::testIndependence_twogroups_LR(vector<vector<char> > const &dataU, vector<int> const &loci1, vector<int> const &loci2, vector<int> const &taken,
        double &mle, double &df)
{
    int i, j, t, kt, k1, k2;
    int szt = (int)taken.size();
    int lsz1 = (int)loci1.size();
    int lsz2 = (int)loci2.size();
    if (lsz1 <= 0 || lsz2 <= 0) return true;

    int hszt = (int)pow((double)alleleCn, (double)szt);
    int hsz1 = (int)pow((double)alleleCn, (double)lsz1);
    int hsz2 = (int)pow((double)alleleCn, (double)lsz2);
    vector<vector<vector<double> > > counts(hszt);
    vector<double> n(hszt, 0);
    for (t = 0; t < hszt; t++)
        counts[t].resize(hsz1, vector<double>(hsz2, 0));

    for (i = 0; i < (int)dataU.size(); i++)
    {
        kt = k1 = k2 = 0;
        for (j = 0; j < szt + lsz1 + lsz2; j++)
        {
            if (j < szt) kt = kt * alleleCn + (int)dataU[i][taken[j]];
            else if (j < szt + lsz1) k1 = k1 * alleleCn + (int)dataU[i][loci1[j - szt]];
            else k2 = k2 * alleleCn + (int)dataU[i][loci2[j - szt - lsz1]];
        }
        counts[kt][k1][k2]++;
        n[kt]++;
    }

    vector<vector<double> > rsum(hszt, vector<double>(hsz1, 0)), csum(hszt, vector<double>(hsz2, 0));
    for (t = 0; t < hszt; t++)
        for (i = 0; i < hsz1; i++)
            for (j = 0; j < hsz2; j++)
            {
                rsum[t][i] += counts[t][i][j];
                csum[t][j] += counts[t][i][j];
            }
    mle = 0;
    df = (double)hszt * (hsz1 - 1) * (hsz2 - 1);
    for (t = 0; t < hszt; t++)
        for (i = 0; i < hsz1; i++)
            for (j = 0; j < hsz2; j++)
            {
                double mean = (double)rsum[t][i] * csum[t][j] / n[t];
                if (counts[t][i][j] > 0)
                    mle += (double)counts[t][i][j] * (log((double)counts[t][i][j]) - log(mean + 0.00001));
            }
    mle = 2. * mle;
    double pv = 1. - gsl_cdf_chisq_P(mle, df);

    return pv;
}
