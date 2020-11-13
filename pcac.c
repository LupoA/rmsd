#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpfr.h>
#include <math.h>
//#include "corefcts.h"
#include <pcac.h>
#include "utils.h"

int tmx = 0;
int nms = 1;
int prec = 128;

mpfr_t *corr, *cov;
mpfr_t *mcorr, *smcorr;
mpfr_t *AXcorr, *AXcov;
mpfr_t *AXmcorr, *AXsmcorr;
mpfr_t *mass, *bmass, *deltam;
int nfit = 0;
int getmass = 0;
char path[128];
char filePP[128];
char fileAX[128];

FILE *fopen_path(char *name)
{
    char fn[128];
    FILE *fp;

    sprintf(fn, "%s/%s", path, name);
    fp = fopen(fn, "w");
    if(fp == NULL)
    {
        error("Failed to open file for writing");
    }

    return fp;
}

void read_corr_PP()
{
    double res, dres, norm;
    char *token;
    char *buf;
    int offset, th;
    FILE *fp;
    th = tmx/2;

    fp = fopen(filePP, "r");
    
    if(fp == NULL)
    {
        error("Unable to open P-P file for reading");
    }
  
    buf = malloc(2048*sizeof(char));

    offset = 0;

    while(fgets(buf, 2048, fp))
    {
  
        // token = strtok(buf, " ");
        token = strtok(buf, "//");
        if(strlen(buf) == 0)
        {
          
            break;
        }

        while(token)
        {
   
            mpfr_set_d(mcorr[offset], atof(token), ROUNDING);
           // printf("%1.8e", atof(token) );

            token = strtok(NULL, "//");
            offset++;
        }
    }

    if(offset != tmx*nms)
    {
        error("Error reading P-P correlator data, number of elements is incorrect");
    }

    fclose(fp);
    fp = fopen_path("PPcorrelator.txt");
    //sym_corr();
    //bootstrap_symm_sample();
    bootstrap_sample(mcorr, corr, cov);
    //full_sample();
    norm = mpfr_get_d(corr[0], ROUNDING);
    
    
    for(int i = 0; i < tmx; i++)
    {
        res = mpfr_get_d(corr[i], ROUNDING);
        dres = mpfr_get_d(cov[i], ROUNDING);
        dres = norm*sqrt(dres);
        fprintf(fp, "%3d %1.8e %1.8e\n", i, res, dres);
    }

    fclose(fp);
    /*
    fs = fopen_path("PP_symm_corr.txt");
    sym_corr();
    bootstrap_symm_sample();
    
    for(int i = 0; i < th+1; i++)
    {
        res = mpfr_get_d(corr[i], ROUNDING);
        dres = mpfr_get_d(cov[i], ROUNDING);
        dres = norm*sqrt(dres);
        fprintf(fs, "%3d %1.8e %1.8e\n", i, res, dres);
    }
    
    fclose(fs);*/
    free(buf);
}

void read_corr_AXIAL()
{
    double res, dres, norm;
    char *token;
    char *buf;
    int offset, th;
    FILE *fp;
    th = tmx/2;

    fp = fopen(fileAX, "r");
    
    if(fp == NULL)
    {
        error("Unable to open AX file for reading");
    }
  
    buf = malloc(2048*sizeof(char));

    offset = 0;

    while(fgets(buf, 2048, fp))
    {/Users/s2002838/Desktop/analyse/16x16x16x32/F2/Mesons
  
        // token = strtok(buf, " ");
        token = strtok(buf, "//");
        if(strlen(buf) == 0)
        {
          
            break;
        }

        while(token)
        {
   
            mpfr_set_d(AXmcorr[offset], atof(token), ROUNDING);
           // printf("%1.8e", atof(token) );

            token = strtok(NULL, "//");
            offset++;
        }
    }

    if(offset != tmx*nms)
    {
        error("Error reading Axial correlator data, number of elements is incorrect");
    }

    fclose(fp);
    fp = fopen_path("Axial_correlator.txt");
    //sym_corr();
    //bootstrap_symm_sample();
    bootstrap_sample(AXmcorr, AXcorr, AXcov);
    
    //full_sample(AXcorr, AXmcorr, AXcov);
    norm = mpfr_get_d(AXcorr[0], ROUNDING);
    
    
    for(int i = 0; i < tmx; i++)
    {
        res = mpfr_get_d(AXcorr[i], ROUNDING);
        dres = mpfr_get_d(AXcov[i], ROUNDING);
        dres = norm*sqrt(dres);
        fprintf(fp, "%3d %1.8e %1.8e\n", i, res, dres);
    }

    fclose(fp);
    /*
    fs = fopen_path("Axial_symm_corr.txt");
    sym_corr();
    bootstrap_symm_sample();
    
    for(int i = 0; i < th+1; i++)
    {
        res = mpfr_get_d(AXcorr[i], ROUNDING);
        dres = mpfr_get_d(AXcov[i], ROUNDING);
        dres = norm*sqrt(dres);
        fprintf(fs, "%3d %1.8e %1.8e\n", i, res, dres);
    }
    
    fclose(fs); */
    free(buf);
}

void full_sample(mpfr_t *thematrixcorr, mpfr_t *thecorr, mpfr_t *thecov)
{
    mpfr_t tmp;
    mpfr_init(tmp);

    if(nms <= 1)
    {
        return;
    }

    for(int i = 0; i < tmx; i++)
    {
        mpfr_set_zero(thecorr[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_add(thecorr[i], thecorr[i], thematrixcorr[j*tmx+i], ROUNDING);
        }
        mpfr_div_d(thecorr[i], thecorr[i], (double)nms, ROUNDING);

        mpfr_set_zero(thecov[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_sub(tmp, thematrixcorr[j*tmx+i], thecorr[i], ROUNDING);
            mpfr_mul(tmp, tmp, tmp, ROUNDING);
            mpfr_add(thecov[i], thecov[i], tmp, ROUNDING);
        }
        mpfr_div_d(thecov[i], thecov[i], (double)(nms-1), ROUNDING);
        mpfr_mul(tmp, thecorr[0], thecorr[0], ROUNDING);
        mpfr_div(thecov[i], thecov[i], tmp, ROUNDING);
    }

    mpfr_clear(tmp);
}


void bootstrap_sample(mpfr_t *thematrixcorr, mpfr_t *thecorr, mpfr_t *thecov)
{
    mpfr_t tmp;
    int *rd;

    if(nms <= 1)
    {
        return;
    }

    mpfr_init(tmp);
    rd = malloc(sizeof(int)*nms);
    randi(rd, nms, 0, nms);

    for(int i = 0; i < tmx; i++)
    {
        mpfr_set_zero(thecorr[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_add(thecorr[i], thecorr[i], thematrixcorr[rd[j]*tmx+i], ROUNDING);
        }
        mpfr_div_d(thecorr[i], thecorr[i], (double)nms, ROUNDING);

        mpfr_set_zero(thecov[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_sub(tmp, thematrixcorr[rd[j]*tmx+i], corr[i], ROUNDING);
            mpfr_mul(tmp, tmp, tmp, ROUNDING);
            mpfr_add(thecov[i], thecov[i], tmp, ROUNDING);
        }
        mpfr_div_d(thecov[i], thecov[i], (double)(nms-1), ROUNDING);
        mpfr_mul(tmp, thecorr[0], thecorr[0], ROUNDING);
        mpfr_div(thecov[i], thecov[i], tmp, ROUNDING);
    }

    mpfr_clear(tmp);
    free(rd);
}


void evaluatePCAC()
{
    FILE *fpcac;
    mpfr_t tmp;
    double *mpcac, dtmp, *pcac, *errpcac;
    mpfr_init(tmp);
    int th=tmx/2;
    
    mpcac = malloc(nms*(th-1)*sizeof(double));

     if(mpcac == 0)
     {
         error("Failed to allocate auxiliary array");
     }
    
    pcac = malloc((th-1)*sizeof(double));
    if(mpcac == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    
    errpcac = malloc((th-1)*sizeof(double));
    if(errpcac == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    fpcac = fopen_path("maxx.txt");
    
    for (int j = 0; j < nms; j++)
    {
        bootstrap_sample(AXcorr, AXmcorr, AXcov);
        bootstrap_sample(corr, mcorr, cov);
        for (int i = 1; i < th; i++)
        {
            mpfr_sub(tmp, AXcorr[i+1], AXcorr[i-1], ROUNDING);
            mpfr_div(tmp, tmp, corr[i], ROUNDING);
            mpfr_div_d(tmp, tmp, (double)(4), ROUNDING);
            mpcac[j*(th-1)+i-1] = mpfr_get_d(tmp, ROUNDING);
        }
    }
    
    for (int i = 0; i < th-1; i++)
    {
        pcac[i]=0;
    
        for (int j = 0; j < nms; j++)
        {
            pcac[i] = pcac[i] + mpcac[j*(th-1)+i-1];
        }
        pcac[i] = pcac[i]/nms;
    }
    
    for (int i = 0; i < th-1; i++)//compute error on the mass
     {
         errpcac[i]=0;
         for (int j = 0; j < nms; j++)
         {
             dtmp=pcac[i]-mpcac[j*(th-1)+i];
             //mpfr_sub(tmp, pcac[i], mpcac[j*(th-1)+i], ROUNDING);
             //mpfr_mul(tmp, tmp, tmp, ROUNDING);
             dtmp=dtmp*dtmp;
             errpcac[i]=errpcac[i]+dtmp;
             //mpfr_add(errpcac[i], errpcac[i], tmp, ROUNDING);
         }
         errpcac[i] /= nms;
         //mpfr_div_d(errpcac[i], errpcac[i], (double)(nms), ROUNDING);
         errpcac[i]=sqrt(errpcac[i]);
        // mpfr_sqrt(errpcac[i], errpcac[i], ROUNDING);
     }
    
    for (int i = 0; i < th-1; i++)
    {
        fprintf(fpcac, "%2d %3.8e %3.8e\n", i+1, pcac[i], errpcac[i]);
    }
    
    fclose(fpcac);

    mpfr_clear(tmp);
}



void prepare_path()
{
    FILE *fp;

    sprintf(path, "./%ld", time(NULL));
    if(mkdir(path, 0755))
    {
        error("Failed to create directory");
    }

    fp = fopen_path("params.txt");
    fprintf(fp, "T      = %d\n", tmx);
    fprintf(fp, "nms    = %d\n", nms);
    fprintf(fp, "filePP   = %s\n", filePP);
    fprintf(fp, "fileAX   = %s\n", fileAX);
    fprintf(fp, "prec   = %d\n", prec);
    fclose(fp);
}


void allocate_data()
{
    int count;


    count = 2*(tmx + tmx + tmx*nms);
    corr = malloc(count*sizeof(mpfr_t)); // tmx
    cov = corr+tmx; //tmx
    mcorr = cov+tmx;    // tmx*nms

    AXcorr=mcorr+tmx*nms; //tmx
    AXcov=AXcorr+tmx; //tmx
    AXmcorr=AXcov+tmx;
    

    if(corr == NULL)
    {
        error("Failed to allocate MPFR variables");
    }

    for(int i = 0; i < count; i++)
    {
        mpfr_init(corr[i]);
    }
    
    
}


int main(int argc, char *argv[])
{
    memset(filePP, 0, sizeof(filePP));
    memset(fileAX, 0, sizeof(fileAX));
    find_int(argc, argv, "-T", &tmx);
    find_int(argc, argv, "-nms", &nms);
    find_str(argc, argv, "-filePP", filePP);
    find_str(argc, argv, "-fileAX", fileAX);
    find_int(argc, argv, "-prec", &prec);

    

    if(tmx <= 0 || nms <= 0 || strlen(filePP) == 0 || strlen(fileAX) == 0)
    {
        printf("Usage: %s -T <int> -nms <int> -file <string> [ options ]\n", argv[0]);
        printf("Options:\n");
        printf("  -T       <int>     temporal extent of correlator\n");
        printf("  -nms     <int>     number of measurements\n");
        printf("  -filePP    <string>  path for reading the PP correlator data\n");
        printf("  -fileAX    <string>  path for reading the AX correlator data\n");
        printf("  -prec    <int>     working precision (default %d decimal places)\n", prec);
        printf("\n");
        return 0;
    }

    srand48(1337);
    mpfr_set_default_prec(3.322*prec);

    //de = (ef-ei)/nsteps;
    //set_params(sigma, lambda, kid, lambdap);
    //set_time_parms(1, tmx/2, tmx);

    prepare_path();
    allocate_data();
    read_corr_PP();
    read_corr_AXIAL();
    evaluatePCAC();
    return 0;
}

