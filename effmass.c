#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpfr.h>
#include <math.h>
//#include "corefcts.h"
#include <effmass.h>
#include "utils.h"

int tmx = 0;
int nms = 1;
int prec = 128;

mpfr_t *corr, *cov;
mpfr_t *mcorr, *smcorr;
mpfr_t *mass, *bmass, *deltam;
int nfit = 0;
int getmass = 0;
char path[128];
char file[128];

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

void read_corr()
{
    double res, dres, norm;
    char *token;
    char *buf;
    int offset, th;
    FILE *fp, *fs;
    th = tmx/2;

    fp = fopen(file, "r");
    
    if(fp == NULL)
    {
        error("Unable to open file for reading");
    }
  
    buf = malloc(4*2048*sizeof(char));

    offset = 0;

    while(fgets(buf, 4*2048, fp))
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
        error("Error reading correlator data, number of elements is incorrect");
    }

    fclose(fp);
    fp = fopen_path("correlator.txt");
    //sym_corr();
    //bootstrap_symm_sample();
    //bootstrap_sample();
    full_sample();
    norm = mpfr_get_d(corr[0], ROUNDING);
    
    
    for(int i = 0; i < tmx; i++)
    {
        res = mpfr_get_d(corr[i], ROUNDING);
        dres = mpfr_get_d(cov[i], ROUNDING);
        dres = norm*sqrt(dres);
        fprintf(fp, "%3d %1.8e %1.8e\n", i, res, dres);
    }

    fclose(fp);
    
    fs = fopen_path("symm_corr.txt");
    sym_corr();
    bootstrap_symm_sample();
    
    for(int i = 0; i < th+1; i++)
    {
        res = mpfr_get_d(corr[i], ROUNDING);
        dres = mpfr_get_d(cov[i], ROUNDING);
        dres = norm*sqrt(dres);
        fprintf(fs, "%3d %1.8e %1.8e\n", i, res, dres);
    }
    
    fclose(fs);
    free(buf);
}


void full_sample()
{
    mpfr_t tmp;
    mpfr_init(tmp);

    if(nms <= 1)
    {
        return;
    }

    for(int i = 0; i < tmx; i++)
    {
        mpfr_set_zero(corr[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_add(corr[i], corr[i], mcorr[j*tmx+i], ROUNDING);
        }
        mpfr_div_d(corr[i], corr[i], (double)nms, ROUNDING);

        mpfr_set_zero(cov[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_sub(tmp, mcorr[j*tmx+i], corr[i], ROUNDING);
            mpfr_mul(tmp, tmp, tmp, ROUNDING);
            mpfr_add(cov[i], cov[i], tmp, ROUNDING);
        }
        mpfr_div_d(cov[i], cov[i], (double)(nms-1), ROUNDING);
        mpfr_mul(tmp, corr[0], corr[0], ROUNDING);
        mpfr_div(cov[i], cov[i], tmp, ROUNDING);
    }

    mpfr_clear(tmp);
}

void bootstrap_sample()
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
        mpfr_set_zero(corr[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_add(corr[i], corr[i], mcorr[rd[j]*tmx+i], ROUNDING);
        }
        mpfr_div_d(corr[i], corr[i], (double)nms, ROUNDING);

        mpfr_set_zero(cov[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_sub(tmp, mcorr[rd[j]*tmx+i], corr[i], ROUNDING);
            mpfr_mul(tmp, tmp, tmp, ROUNDING);
            mpfr_add(cov[i], cov[i], tmp, ROUNDING);
        }
        mpfr_div_d(cov[i], cov[i], (double)(nms-1), ROUNDING);
        mpfr_mul(tmp, corr[0], corr[0], ROUNDING);
        mpfr_div(cov[i], cov[i], tmp, ROUNDING);
    }

    mpfr_clear(tmp);
    free(rd);
}


void sym_corr()
{

    int th = tmx/2;

    
    for (int j=0; j<nms; j++)
    {
        mpfr_set(smcorr[j*(th+1)], mcorr[j*tmx], ROUNDING); // c(0) --> c(0)
        mpfr_set(smcorr[j*(th+1)+th], mcorr[j*tmx+th], ROUNDING); // c(T/2) --> c(T/2)
        for (int i=1; i<th; i++)
        {
            mpfr_add(smcorr[j*(th+1)+i], mcorr[j*tmx+i], mcorr[j*tmx+(tmx-i)], ROUNDING);   // c(1)=c(1)+c(31) ... c(T/2-1) = c(T/2-1)+c(T/2+1)
            mpfr_div_d(smcorr[j*(th+1)+i], smcorr[j*(th+1)+i], (double)(2.), ROUNDING);
        }
    }
    

    
}

void bootstrap_symm_sample()
{
    mpfr_t tmp;
    int *rd, th;
    th = tmx/2;

    if(nms <= 1)
    {
        return;
    }
    
    mpfr_init(tmp);
    rd = malloc(sizeof(int)*nms);
    randi(rd, nms, 0, nms);

    for(int i = 0; i < th+1; i++)
    {
        mpfr_set_zero(corr[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_add(corr[i], corr[i], smcorr[rd[j]*(th+1)+i], ROUNDING);
        }
        mpfr_div_d(corr[i], corr[i], (double)nms, ROUNDING);

        mpfr_set_zero(cov[i], 1);
        for(int j = 0; j < nms; j++)
        {
            mpfr_sub(tmp, smcorr[rd[j]*(th+1)+i], corr[i], ROUNDING);
            mpfr_mul(tmp, tmp, tmp, ROUNDING);
            mpfr_add(cov[i], cov[i], tmp, ROUNDING);
        }
        mpfr_div_d(cov[i], cov[i], (double)(nms-1), ROUNDING);
        mpfr_mul(tmp, corr[0], corr[0], ROUNDING);
        mpfr_div(cov[i], cov[i], tmp, ROUNDING);
    }
    
    

    mpfr_clear(tmp);
    free(rd);
}

void bootmass()
{
    FILE *fp;
    mpfr_t tmp;
    
    int th = tmx/2;
    mpfr_init(tmp);
    
    fp = fopen_path("mass.txt");
    
    for (int j = 0; j < nms; j++)  // number of bootstrap samples = nms
    {
        bootstrap_symm_sample(); // fills corr[] with a bootstrapped symmetrised sample
        //full_sample();
        for (int i = 1; i < th; i++)
        {
            mpfr_add(tmp, corr[i+1], corr[i-1], ROUNDING);

            mpfr_div(tmp, tmp, corr[i], ROUNDING);
            
            mpfr_div_d(tmp, tmp, (double)2, ROUNDING);
            
            mpfr_acosh(mass[j*(th-1)+i-1], tmp, ROUNDING); // I am shifting the masses so they start from m(0) instead of m(1)
        }
    }
    
    for (int i = 0; i < th-1; i++) //compute the mass
    {
        mpfr_set_zero(tmp, 1);
        for (int j = 0; j < nms; j++)
        {
            mpfr_add(tmp, tmp, mass[j*(th-1)+i], ROUNDING);
        }
        mpfr_div_d(bmass[i], tmp, (double)(nms), ROUNDING); // mass computed as the average of the nms bootstrap samples
    }
    
     for (int i = 0; i < th-1; i++)//compute error on the mass
     {
         mpfr_set_zero(deltam[i], 1);
         for (int j = 0; j < nms; j++)
        {
             mpfr_sub(tmp, bmass[i], mass[j*(th-1)+i], ROUNDING);
             mpfr_mul(tmp, tmp, tmp, ROUNDING);
             mpfr_add(deltam[i], deltam[i], tmp, ROUNDING);
        }
         
         mpfr_div_d(deltam[i], deltam[i], (double)(nms), ROUNDING);
         mpfr_sqrt(deltam[i], deltam[i], ROUNDING);
     }
    // now print the mass. printing i+1 because I shifted the pointer to start from zero
    for (int i = 0; i < th-1; i++)
    {
        fprintf(fp, "%2d %3.8e %3.8e\n", i+1, mpfr_get_d(bmass[i], ROUNDING), mpfr_get_d(deltam[i], ROUNDING));
    }
    
    fclose(fp);
    mpfr_clear(tmp);
}


void fitmass()
{
    char *end;
    char buf[LINE_MAX];
    mpfr_t mfit, deltamfit, tmp;
    
    mpfr_inits(mfit, deltamfit, tmp, NULL);
    int th = tmx/2;
    
    printf("Enter number of points for linear fit, starting from the last\n");

    
    do
    {
         if (!fgets(buf, sizeof buf, stdin))
            break;

         buf[strlen(buf) - 1] = 0;

         nfit = strtol(buf, &end, 10);
    } while (end != buf + strlen(buf));
    
    
    mpfr_set_zero(mfit, 1);

    for(int i = th-1-nfit; i < th-1; i++)
    {
        mpfr_add(mfit, mfit, bmass[i], ROUNDING);
    }
    mpfr_div_d(mfit, mfit, (double)(nfit), ROUNDING);
        
    mpfr_set_zero(deltamfit, 1);
    for(int i = th-1-nfit; i < th-1; i++)
    {
        mpfr_sub(tmp, mfit, bmass[i], ROUNDING);
        mpfr_mul(tmp, tmp, tmp, ROUNDING);
        mpfr_add(deltamfit, deltamfit, tmp, ROUNDING);
    }
    mpfr_div_d(deltamfit, deltamfit, (double)(nfit), ROUNDING);
    mpfr_sqrt(deltamfit, deltamfit, ROUNDING);
        
    
    printf("m = %3.8e ± %3.8e\n", mpfr_get_d(mfit, ROUNDING), mpfr_get_d(deltamfit, ROUNDING));
    
    
    mpfr_clears(mfit, deltamfit, tmp, NULL);
    
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
    fprintf(fp, "file   = %s\n", file);
    fprintf(fp, "prec   = %d\n", prec);
    fclose(fp);
}


void allocate_data()
{
    int count;
    int th = tmx/2;

    //count = 4*tmx+2*tmx*nms;
    count = tmx + tmx + tmx*nms + nms*th-nms + tmx-2 + nms*th+nms; //+ nsteps;
    corr = malloc(count*sizeof(mpfr_t)); // tmx
    cov = corr+tmx; //tmx
    mcorr = cov+tmx;    // tmx*nms
    
    mass = mcorr+tmx*nms;   //nms*(T/2-1) = nms*T/2-nms
    bmass = mass+th*nms-nms;      // T/2-1
    deltam = bmass+th-1;         // T/2-1
    smcorr = deltam+tmx/2-1;    // (T/2+1)nms = nms*T/2 + nms
    //deltazero = smcorr + nms*th+nms;    // nsteps

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
    memset(file, 0, sizeof(file));
    find_int(argc, argv, "-T", &tmx);
    find_int(argc, argv, "-nms", &nms);
    find_str(argc, argv, "-file", file);
    find_int(argc, argv, "-prec", &prec);

    

    if(tmx <= 0 || nms <= 0 || strlen(file) == 0)
    {
        printf("Usage: %s -T <int> -nms <int> -file <string> [ options ]\n", argv[0]);
        printf("Options:\n");
        printf("  -T       <int>     temporal extent of correlator\n");
        printf("  -nms     <int>     number of measurements\n");
        printf("  -file    <string>  path for reading the correlator data\n");
 
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
    read_corr();
    bootmass();
    fitmass();
    return 0;
}

