// main.cpp : Defines the entry point for the console application.
//

#include "datastructure.h"
#include "BiCluster.h"

void outputResult(const char *filename, vector<double> const &para, vector<OUTPUT> const &output, vector<int> const &totalcounts, vector<double> const &group1, vector<double> const &group2, vector<string> const &rs, vector<vector<int> > const &positions, int samplesize);
int loadCaseControl(const char *filename, vector<vector<char> > &dataC, vector<vector<char> > &dataU, vector<string> &rs, vector<vector<int> > &positions);
void loadParameters(const char *filename);

bool SNPid = true;
bool SNPpos = true;
int burnin = 100000;
int mcmc = 100000;
int thin = 1000;
int indRuns = 3;
bool multiTry = false;
bool singleOnly = false;
double prior1 = 0.01;
double prior2 = 0.01;
int mindistance = 1000000;
char outputfile[200];
char inputfile[200];
int startTries = 20;
int tryLen = 20;
double pThreshold = 0.1;
extern double autoRestart;

typedef struct MySortTypeP {
    int chr;
    int pos;
    int index;
} MySortTypeP;

bool operator<(const MySortTypeP &a, const MySortTypeP &b)
{
    if (a.chr != b.chr) return a.chr < b.chr;
    else return a.pos < b.pos;
}

int main(int argc, char* argv[])
{
    sprintf(inputfile, "%s", "data.txt");
    sprintf(outputfile, "%s", "result.txt");
    loadParameters("parameters.txt");
    if (argc > 1)
    {
        int i;
        sprintf(inputfile, "%s", argv[1]);
        sprintf(outputfile, "%s", argv[2]);
    }

    vector<vector<char> > dataC, dataU;
    vector<string> rs;
    vector<vector<int> > positions;
    printf("Loading data...\n");
    int alleleCn = loadCaseControl(inputfile, dataC, dataU, rs, positions);
    int NC = (int)dataC.size();
    int NU = (int)dataU.size();
    int L = 0;
    if (NC > 0) L = (int)dataC[0].size();
    if (burnin <= 0) burnin = L * 10;
    if (mcmc <= 0) mcmc = L * max(100, L);
    if (thin <= 0) thin = L;
    if (startTries <= 0) burnin = L * 100;
    if (tryLen <= 0) tryLen = L * 20;

    printf("___________________________________________________________\n");
    printf("Input file: \"%s\"\n", inputfile);
    printf("Maximum number of distinct alleles per marker: %d\n", alleleCn);
    printf("\n[PARAMETERS]\n");
    printf("# of Cases = %d\n# of Controls = %d\n# of Markers = %d\n", NC, NU, L);
    printf("# of chains = %d\nBurnin in each chain = %d\nIteration after burnin = %d\nSample at every %d updates\n", indRuns, burnin, mcmc, thin);
    printf("Priors: p1 = %f, p2 = %f\n", prior1, prior2);
    printf("Threhold for Bonferroni corrected p-values: %f\n", pThreshold);
    if (singleOnly) printf("Test for marginal associations only.\n");
    else printf("Test for both marginal and epistatic associations.\n");
    printf("___________________________________________________________\n");
    if (L <= 0) return 0;

    BiCluster bcluster(alleleCn);
    vector<OUTPUT> output;
    vector<int> totalcounts;
    vector<double> group1, group2;

    bcluster.mBurnin = burnin;
    bcluster.mMCMC = mcmc;
    bcluster.mThin = thin;
    bcluster.minDistance = mindistance;
    bcluster.mPrior1 = prior1;
    bcluster.mPrior2 = prior2;

    time_t st, ed;
    time(&st);

    int samplesize = bcluster.runChains(dataC, dataU, pThreshold, output, totalcounts, group1, group2, indRuns, positions, multiTry, singleOnly, 0, startTries, tryLen);
    vector<double> para;
    para.push_back((int)dataC.size());
    para.push_back((int)dataU.size());
    para.push_back((int)dataC[0].size());
    para.push_back(indRuns);
    para.push_back(bcluster.mBurnin);
    para.push_back(bcluster.mMCMC);
    para.push_back(bcluster.mThin);
    para.push_back(bcluster.mPrior1);
    para.push_back(bcluster.mPrior2);
    outputResult(outputfile, para, output, totalcounts, group1, group2, rs, positions, samplesize);

    time(&ed);
    printf("\nDone!\nTotal Time = %d sec.\nResults are in the file \"", (int)ed - st);
    printf("%s\".\n", outputfile);

    return 0;
}

void outputResult(const char *filename, vector<double> const &para, vector<OUTPUT> const &output, vector<int> const &totalcounts, vector<double> const &group1, vector<double> const &group2,
                  vector<string> const &rs, vector<vector<int> > const &positions, int samplesize)
{
    int i, j, k;

    FILE *f = fopen(filename, "w");
    fprintf(f, "[INPUT_FILE]\n\"%s\"\n", inputfile);

    fprintf(f, "\n[PARAMETERS]\n");
    fprintf(f, "# of Cases = %d\n# of Controls = %d\n# of Markers = %d\n", (int)para[0], (int)para[1], (int)para[2]);
    fprintf(f, "# of chains = %d\nBurnin in each chain = %d\nIteration after burnin = %d\nSample at every %d updates\n", (int)para[3], (int)para[4], (int)para[5], (int)para[6]);
    fprintf(f, "Priors: p1 = %f, p2 = %f\n", para[7], para[8]);
    fprintf(f, "Threhold for Bonferroni corrected p-values: %f\n", pThreshold);
    if (singleOnly) printf("Test for marginal associations only.\n");
    else printf("Test for both marginal and epistatic associations.\n");

    fprintf(f, "\n[POSTERIOR_SIZES]\n");
    fprintf(f, "# of Marginal Associations:\n");
    for (k = 0; k < (int)group1.size(); k++)
        fprintf(f, "%d: %1.3f, ", k, group1[k]);

    fprintf(f, "\nInteraction Sizes:\n");
    for (k = 0; k < (int)group2.size(); k++)
    {
        if (k == 0) fprintf(f, "0: %1.3f, ", group2[k]);
        else fprintf(f, "%d: %1.3f, ", k + 1, group2[k]);
    }

    fprintf(f, "\n\n[DETECTED_ASSOCIATIONS]\nMarker_Group\tPvalue\n");
    for (i = 0; i < (int)output.size(); i++)
    {
        fprintf(f, "( ");
        for (j = 0; j < (int)output[i].markers.size(); j++)
            fprintf(f, "%d ", output[i].markers[j]);
        fprintf(f, " )\t%f\n", output[i].pvalue);
    }
    fprintf(f, "\n[POSTERIOR_DISTRIBUTION]\nMarginal + Interaction = Total\n");
    int L = (int)totalcounts.size() / 2;
    for (i = 0; i < L; i++)
    {
        fprintf(f, "%d:", i);
        if ((int)rs.size() > 0)
            fprintf(f, "%s\t", rs[i].c_str());
        if ((int)positions.size() > 0)
        {
            if (positions[i][0] < 23)
                fprintf(f, "Chr%d:%d\t", positions[i][0], positions[i][1]);
            else if (positions[i][0] == 23)
                fprintf(f, "ChrX:%d\t", positions[i][1]);
            else fprintf(f, "ChrY:%d\t", positions[i][1]);
        }
        fprintf(f, "%f\t+\t%f\t=\t%f\n", (double)totalcounts[i] / samplesize, (double)totalcounts[i + L] / samplesize, (double)(totalcounts[i] + totalcounts[i + L]) / samplesize);
    }
    fclose(f);
}

//disease status are represented by integers, take mean, and larger values will be taken as cases
int loadCaseControl(const char *filename, vector<vector<char> > &dataC, vector<vector<char> > &dataU, vector<string> &rs, vector<vector<int> > &positions)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL)
    {
        printf("Cannot open the data file \"%s\"\n", filename);
        return 0;
    }
    dataC.clear();
    dataU.clear();
    rs.clear();
    positions.clear();

    int i, j;
    char tmp[100000];
    vector<int> status;
    int alleleCn = 2;

    double mean = 0;
    fgets(tmp, 100000, f);
    for (i = 0; i < (int)strlen(tmp); i++)
        if ((tmp[i] >= 48 && tmp[i] < 58) || tmp[i] == '-') //integers
        {
            j = atoi(&tmp[i]);
            status.push_back(j);
            mean += (double)j;
        }
    mean /= (double)status.size();

    vector<vector<char> > tmpgenotype;
    while (fgets(tmp, 100000, f) != NULL)
    {
        if ((int)strlen(tmp) < 4) continue;
        i = 0;
        if (SNPid)
        {
            char str[100];
            j = 0;
            for (i = i; i < (int)strlen(tmp); i++)
            {
                str[j++] = tmp[i];
                if (tmp[i] == ' ' || tmp[i] == '\t')
                    break;
            }
            str[j] = 0;
            rs.push_back(str);
        }
        if (SNPpos)
        {
            vector<int> p(2);
            for (i = i; i < (int)strlen(tmp); i++)
                if (tmp[i] == 'C' || tmp[i] == 'c')
                    break;
            if (tmp[i + 3] == 'X' || tmp[i + 3] == 'x') p[0] = 23;
            else if (tmp[i + 3] == 'Y' || tmp[i + 3] == 'y') p[0] = 24;
            else p[0] = atoi(&tmp[i + 3]);
            for (i = i + 3; i < (int)strlen(tmp); i++)
                if (tmp[i] == ' ' || tmp[i] == '\t')
                    break;
            for (i = i; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            p[1] = atoi(&tmp[i]);
            for (i = i; i < (int)strlen(tmp); i++)
                if (tmp[i] == ' ' || tmp[i] == '\t')
                    break;
            positions.push_back(p);
        }
        vector<char> row;
        do
        {
            for (i = i; i < (int)strlen(tmp); i++)
                if ((tmp[i] >= 48 && tmp[i] < 58) || tmp[i] == '-')
                    break;
            if (i < (int)strlen(tmp))
            {
                j = atoi(&tmp[i]);
                row.push_back(j);
                if (j >= 0 && j >= alleleCn) alleleCn = j + 1;
            }
            i++;
        } while (i < (int)strlen(tmp));

        //makeup missing alleles
        vector<int> cn(alleleCn, 0);
        for (i = 0; i < (int)row.size(); i++)
            if (row[i] >= 0) cn[row[i]]++;
        for (i = 1; i < (int)cn.size(); i++)
            cn[i] += cn[i - 1];
        for (i = 0; i < (int)row.size(); i++)
            if (row[i] < 0)
            {
                int k = (int)((double)rand() / RAND_MAX * (double)cn[alleleCn - 1]);
                for (j = 0; j < alleleCn - 1; j++)
                    if (k < cn[j])
                        break;
                row[i] = j;
            }
        tmpgenotype.push_back(row);
    }
    fclose(f);

    int L = (int)tmpgenotype.size();
    //if positions are provided, sort SNPs
    vector<MySortTypeP> posmap;
    if ((int)positions.size() == (int)tmpgenotype.size())
    {
        posmap.resize(L);
        for (i = 0; i < L; i++)
        {
            posmap[i].chr = positions[i][0];
            posmap[i].pos = positions[i][1];
            posmap[i].index = i;
        }
        sort(posmap.begin(), posmap.end());
        for (i = 0; i < L; i++)
        {
            positions[i][0] = posmap[i].chr;
            positions[i][1] = posmap[i].pos;
        }
        if ((int)rs.size() > 0)
        {
            vector<string> rsbak = rs;
            for (i = 0; i < L; i++)
                rs[i] = rsbak[posmap[i].index];
        }
    }

    for (i = 0; i < (int)status.size(); i++)
    {
        vector<char> row(L);
        for (j = 0; j < (int)tmpgenotype.size(); j++)
        {
            if ((int)posmap.size() > 0)
                row[j] = tmpgenotype[posmap[j].index][i];
            else row[j] = tmpgenotype[j][i];
        }
        if ((double)status[i] <= mean)
            dataU.push_back(row);
        else dataC.push_back(row);
    }
    return alleleCn;
}

void loadParameters(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL)
    {
        printf("Cannot open the parameter file \"%s\"\n", filename);
        return;
    }
    int i;
    char tmp[1000];
    while (fgets(tmp, 1000, f) != NULL)
    {
        char str[1000];
        sprintf(str, "%s", tmp);
        str[11] = 0;
        if (strcmp(str, "INC_SNP_POS") == 0)
        {
            for (i = 11; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            if (tmp[i] == 48) SNPpos = false;
            else SNPpos = true;
        }
        else if (strcmp(str, "SINGLE_ONLY") == 0)
        {
            for (i = 11; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            if (tmp[i] == 48) singleOnly = false;
            else singleOnly = true;
        }
        else if (strcmp(str, "MINDISTANCE") == 0)
        {
            for (i = 11; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            mindistance = atoi(&tmp[i]);
        }
        else if (strcmp(str, "INITIALTRYS") == 0)
        {
            for (i = 11; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            startTries = atoi(&tmp[i]);
        }
        else if (strcmp(str, "AUTORESTART") == 0)
        {
            for (i = 11; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            autoRestart = (double)atoi(&tmp[i]);
        }
        else if (strcmp(str, "P_THRESHOLD") == 0)
        {
            for (i = 11; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            pThreshold = atof(&tmp[i]);
        }

        str[10] = 0;
        if (strcmp(str, "INC_SNP_ID") == 0)
        {
            for (i = 10; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            if (tmp[i] == 48) SNPid = false;
            else SNPid = true;
        }
        else if (strcmp(str, "TRY_LENGTH") == 0)
        {
            for (i = 10; i < (int)strlen(tmp); i++)
                if (tmp[i] >= 48 && tmp[i] < 58)
                    break;
            tryLen = atoi(&tmp[i]);
        }

        str[7] = 0;
        if (strcmp(str, "INPFILE") == 0)
        {
            for (i = 7; i < (int)strlen(tmp); i++)
                if (tmp[i] != ' ' && tmp[i] != '\t')
                    break;
            int j;
            for (j = i + 1; j < (int)strlen(tmp); j++)
                if (tmp[j] == ' ' || tmp[j] == '\t' || tmp[j] == '"')
                    break;
            tmp[j] = 0;
            if (tmp[i] == '"') sprintf(inputfile, "%s", &tmp[i + 1]);
            else sprintf(inputfile, "%s", &tmp[i]);
        }
        else if (strcmp(str, "OUTFILE") == 0)
        {
            for (i = 7; i < (int)strlen(tmp); i++)
                if (tmp[i] != ' ' && tmp[i] != '\t')
                    break;
            int j;
            for (j = i + 1; j < (int)strlen(tmp); j++)
                if (tmp[j] == ' ' || tmp[j] == '\t' || tmp[j] == '"')
                    break;
            tmp[j] = 0;
            if (tmp[i] == '"') sprintf(outputfile, "%s", &tmp[i + 1]);
            else sprintf(outputfile, "%s", &tmp[i]);
        }
        str[6] = 0;
        if (strcmp(str, "PRIOR1") == 0)
        {
            for (i = 6; i < (int)strlen(tmp); i++)
                if (tmp[i] != ' ' && tmp[i] != '\t')
                    break;
            prior1 = atof(&tmp[i]);
        }
        else if (strcmp(str, "PRIOR2") == 0)
        {
            for (i = 6; i < (int)strlen(tmp); i++)
                if (tmp[i] != ' ' && tmp[i] != '\t')
                    break;
            prior2 = atof(&tmp[i]);
        }
        else if (strcmp(str, "BURNIN") == 0)
        {
            for (i = 6; i < (int)strlen(tmp); i++)
                if (tmp[i] != ' ' && tmp[i] != '\t')
                    break;
            burnin = atoi(&tmp[i]);
        }

        str[5] = 0;
        if (strcmp(str, "CHAIN") == 0)
        {
            for (i = 5; i < (int)strlen(tmp); i++)
                if (tmp[i] != ' ' && tmp[i] != '\t')
                    break;
            indRuns = atoi(&tmp[i]);
        }
        str[4] = 0;
        if (strcmp(str, "THIN") == 0)
        {
            for (i = 4; i < (int)strlen(tmp); i++)
                if (tmp[i] != ' ' && tmp[i] != '\t')
                    break;
            thin = atoi(&tmp[i]);
        }
        else if (strcmp(str, "MCMC") == 0)
        {
            for (i = 4; i < (int)strlen(tmp); i++)
                if (tmp[i] != ' ' && tmp[i] != '\t')
                    break;
            mcmc = atoi(&tmp[i]);
        }
    }

    fclose(f);
}
