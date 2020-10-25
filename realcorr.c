#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpfr.h>
#include <math.h>
#include "corefcts.h"
#include "utils.h"
// imp
double sigma = 0.1;
double lambda = 0.05;
double lambdap = 0.05;
double lambdastar = 0.05;
int nsteps = 100;
double ei = 0.0;
double ef = 1.0;
int tmx = 0;
int nms = 1;
int prec = 128;
double emin = 0.0;
double escan = 0.0;
int kid = -1;
int iter = 0;
double *bootrho0, *bootdeltarho0, *bootrho1, *bootdeltarho1;
double *rho1, *rho0, *drho1, *drho0;
int nfit = 0;
int getmass = 0;
int nboot = 40;

mpfr_t *corr, *cov;
mpfr_t *expected;
mpfr_t *mcorr, *smcorr, *deltazero;
mpfr_t *mass, *bmass, *deltam;


char path[128];
char file[128];
double de;

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




void get_spectral_density()
{
	double val, ag, mag, bg, mbg;
	double *res, rho, stat, sys;
	mpfr_t estar, ec, tmp0, tmp1, e0, difff;
	FILE *sfp, *rfp;

	mpfr_inits(estar, ec, tmp0, tmp1, e0, difff, NULL);
	mpfr_set_d(e0, emin, ROUNDING);
	res = malloc(nms*sizeof(double));

	if(res == 0)
	{
		error("Failed to allocate auxiliary array");
	}

	rfp = fopen_path("results.txt");
	sfp = fopen_path("systematic.txt");

	for(int n = 0; n < nsteps; n++)
	{
		mpfr_set_d(estar, ei+n*de, ROUNDING);
		rho = 0;
		mag = 0;
		mbg = 0;
		sys = 0;
        mpfr_set_zero(deltazero[n], 1);

		printf("\rStarting step %d/%d", n+1, nsteps);
		fflush(stdout);

		for(int k = 0; k < nms; k++)
		{
			//bootstrap_sample();
            bootstrap_symm_sample();
			rm_method_cosh(e0, estar, cov);

			transform(e0, estar, corr, cov, &val, &ag, &bg);
			res[k] = val;
            rho0[n+nsteps*k] = val;
			rho += val;
			mag += ag;
			mbg += bg;

			mpfr_set(ec, estar, ROUNDING);
			if(kid == 3)
			{
				mpfr_add_d(ec, ec, sigma, ROUNDING);
			}

			deltabar(tmp1, ec);
			delta(tmp0, estar, ec);
			
            mpfr_sub(tmp0, tmp1, tmp0, ROUNDING);   //deltabar-delta
            mpfr_add(deltazero[n], deltazero[n], tmp0, ROUNDING); // NB the sign of this;
		//	mpfr_div(tmp0, tmp0, tmp1, ROUNDING); //1 -delta\bardelta
            delta(tmp1, estar, ec);
            mpfr_div(tmp0, tmp0, tmp1, ROUNDING);      // 1 - deltabar\delta

			val = mpfr_get_d(tmp0, ROUNDING);
			val = fabs(val * res[k]);
            drho0[n+nsteps*k] = val;
			sys += val;
            
		}

		rho /= nms;
        mpfr_div_d(deltazero[n], deltazero[n], (double)(nms), ROUNDING);
		mag /= nms;
		mbg /= nms;
		sys *= 0.6827;
		sys /= nms;
		stat = 0;

		if(nms > 1)
		{
			for(int k = 0; k < nms; k++)
			{
				stat += (res[k]-rho)*(res[k]-rho);
			}
			stat /= (nms-1);
			stat = sqrt(stat);
		}

		fprintf(rfp, "%1.8e %1.8e %1.8e %1.8e %1.8e %1.8e\n",
					ei+n*de, rho, stat, sys, sqrt(stat*stat+sys*sys), mag+mbg);

		//bootstrap_sample();
        bootstrap_symm_sample();
		rm_method_cosh(e0, estar, cov);

		for(int k = 0; k < nsteps; k++)
		{
            // remember that estar = ei+n*de
			mpfr_set_d(ec, ei+k*de, ROUNDING);
			deltabar(tmp0, ec);
			delta(tmp1, estar, ec);
			fprintf(sfp, "%1.8e %1.8e %1.8e %1.8e\n",
						ei+n*de, ei+k*de, mpfr_get_d(tmp1, ROUNDING),
						mpfr_get_d(tmp0, ROUNDING));

		}

		fflush(sfp);
		fflush(rfp);
	}

	fclose(rfp);
	fclose(sfp);

	printf("\n");
	mpfr_clears(estar, ec, tmp0, tmp1, e0, NULL);
	free(res);
}


void get_spectral_density_improv()
{
    double val, ag, mag, bg, mbg, diff;
    double *res, rho, stat, sys, *resbis, rhobis, valbis;
    mpfr_t estar, ec, tmp0, tmp1, e0;
    FILE *sfp, *rfp;

    mpfr_inits(estar, ec, tmp0, tmp1, e0, NULL);
    mpfr_set_d(e0, emin, ROUNDING);
    res = malloc(nms*sizeof(double));

    if(res == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    resbis = malloc(nms*sizeof(double));
    
    if(resbis == 0)
    {
        error("Failed to allocate auxiliary array");
    }


    rfp = fopen_path("results_lvl1.txt");
    sfp = fopen_path("systematic_lvl1.txt");

    for(int n = 0; n < nsteps; n++)
    {
        mpfr_set_d(estar, ei+n*de, ROUNDING);
        rho = 0;
        mag = 0;
        mbg = 0;
        sys = 0;
        diff = 0;
        rhobis = 0;

        printf("\rStarting step %d/%d", n+1, nsteps);
        fflush(stdout);

        for(int k = 0; k < nms; k++)
        {
           // bootstrap_sample();
            bootstrap_symm_sample();
            
            rm_method_cosh_lvl1(e0, estar, cov);
            transform_lvl1(e0, estar, corr, cov, &val, &ag, &bg);
            res[k] = val;
            transform_only_step1(corr, &valbis);
            resbis[k] = valbis;
            
            rho += val;
            rho1[n+nsteps*k] = val;
            rhobis += valbis;
            
            mag += ag;
            mbg += bg;

            mpfr_set(ec, estar, ROUNDING);
            if(kid == 3)
            {
                mpfr_add_d(ec, ec, sigma, ROUNDING);
            }

 
            
            deltabar_lvl1(tmp1, ec);
            deltabar(tmp0, ec);
            mpfr_add(tmp0, tmp0, tmp1, ROUNDING);
            delta(tmp1, estar, ec);
           // mpfr_sub(tmp1, tmp1, tmp0, ROUNDING);        // these gives bar quantities
           // mpfr_div(tmp1, tmp1, tmp0, ROUNDING);     // in the denominator
            mpfr_sub(tmp0, tmp1, tmp0, ROUNDING);
            mpfr_div(tmp0, tmp0, tmp1, ROUNDING);
          //  diff = mpfr_get_d(tmp1, ROUNDING);
            diff = mpfr_get_d(tmp0, ROUNDING);
            diff = fabs(diff);
            diff = diff*(val+valbis);
            
            drho1[n+nsteps*k] = diff;
            
            sys += diff;
            
            
        }

        rho /= nms;
        rhobis /= nms;
        mag /= nms;
        mbg /= nms;
        sys *= 0.6827;
        sys /= nms;

        stat = 0;

        if(nms > 1)
        {
            for(int k = 0; k < nms; k++)
            {
                stat += (res[k]-rho)*(res[k]-rho);
            }
            stat /= (nms-1);
            stat = sqrt(stat);
        }

        fprintf(rfp, "%1.8e %1.8e %1.8e %1.8e %1.8e %1.8e\n",
                    ei+n*de, rho, stat, sys, sqrt(stat*stat+sys*sys), mag+mbg);

        //bootstrap_sample();
        bootstrap_symm_sample();
        rm_method_cosh(e0, estar, cov);
        
        for(int k = 0; k < nsteps; k++)
        {
            mpfr_set_d(ec, ei+k*de, ROUNDING);
            deltabar_lvl1(tmp0, ec);
            fprintf(sfp, "%1.8e %1.8e %1.8e %1.8e\n",
                        ei+n*de, ei+k*de, mpfr_get_d(deltazero[k], ROUNDING),
                        mpfr_get_d(tmp0, ROUNDING));
        }

        fflush(sfp);
        fflush(rfp);
    }

    fclose(rfp);
    fclose(sfp);

    printf("\n");
    mpfr_clears(estar, ec, tmp0, tmp1, e0, NULL);
    free(res);
}


void scan_lambda()
{
    double val, ag, bg, dl, lstar, wg, ll;
	mpfr_t estar, e0;
	FILE *fp;

	fp = fopen_path("lambda.txt");

	mpfr_inits(e0, estar, NULL);
	mpfr_set_d(e0, emin, ROUNDING);
	mpfr_set_d(estar, escan, ROUNDING);
    
    ll = 0;
    
    dl = 0.001;

	
        wg = 0;

	for(double l = dl; l < 0.5+dl; l += dl)
	{
		set_params(sigma, l, kid, ll);
		rm_method_cosh(e0, estar, cov);
		transform(e0, estar, corr, cov, &val, &ag, &bg);
       	        if (ag+bg > wg) {
               		 wg = ag+bg;
               	         lstar = l;
        	}
		fprintf(fp, "%1.8e %1.8e %1.8e %1.8e\n", l, ag, bg, ag+bg);
        
        
	}
	
        printf("Suggested choice for lambda (lvl 0): %1.8e\n", lstar);
        lambdastar = lstar;
	
	mpfr_clears(e0, estar, NULL);
	fclose(fp);
}


void scan_lambda_lvl1()
{
    double val, ag, bg, dl, lstar, wg, ll;
    mpfr_t estar, e0;
    FILE *fp;

    fp = fopen_path("lambda_lvl_1.txt");

    mpfr_inits(e0, estar, NULL);
    
    mpfr_set_d(e0, emin, ROUNDING);
    mpfr_set_d(estar, escan, ROUNDING);
    
    lstar = 0.31415926;
    ll = 0;
    dl = 0.001;

    
        wg = 0;

    for(double l = dl; l < 1+dl; l += dl)
    {
        set_params(sigma, lambdastar, kid, l);
        rm_method_cosh_lvl1(e0, estar, cov);
        transform_lvl1(e0, estar, corr, cov, &val, &ag, &bg);
            if (ag+bg > wg)
            {
                wg = ag+bg;
                lstar = l;
            }
        fprintf(fp, "%1.8e %1.8e %1.8e %1.8e\n", l, ag, bg, ag+bg);
        
    }
    
        printf("Suggested choice for lambda (lvl 1): %1.8e\n", lstar);
    
    mpfr_clears(e0, estar, NULL);
    fclose(fp);
}


void peak_bootstrap()
{
    double *peak, psum, *epeak, esum, dsum, tmp;
    int *rd;
    int counter, flag;
    FILE *fp;

    fp = fopen_path("PeakLocation.txt");
    
    peak = malloc(nms*sizeof(double));

    if(peak == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    epeak = malloc(nms*sizeof(double));

    if(epeak == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    dsum = 0;
    counter = 0;
    psum = 0;
    esum = 0;
    flag = 0;
    
    tmp = 0.61 - ei;
    tmp = tmp/de;
    int cut = tmp;

    

    
    rd = malloc(sizeof(int)*nms);
    randi(rd, nms, 0, nms);
    
    
    

    
    
    for (int k = 0; k < nms; k++)
    {

            peak[k] = 0;
            for (int n = 0; n < cut; n++)
            {
                if (peak[k] < rho0[n+nsteps*k])
               // if (peak[k] < rho0[n+nsteps*k])
                {
                    peak[k] = rho0[n+nsteps*k];
                 //   peak[k] = rho0[n+nsteps*k];
                    epeak[k] = ei+n*de;


                }
                
            }

       
    }
    
    for (int k = 0; k < nms; k++)
    {
        psum = psum + peak[rd[k]];
        esum = esum + epeak[rd[k]];
    }
    
    psum /= nms;
    esum /= nms;
    
    for (int k = 0; k < nms; k++) {
        tmp = esum - epeak[k];
        tmp = tmp*tmp;
        dsum += tmp;
    }
    
    dsum /= nms;
    dsum = sqrt(dsum);
    
    
    printf("Peak at %3.6f ± %3.6f\n", esum, dsum);
    fprintf(fp, "Peak at %3.6f ± %3.6f\n", esum, dsum);

    
    free(rd);
    free(peak);
    free(epeak);
    fclose(fp);
}

void peak_bootstrap_iter()
{
      double *peak, psum, *epeak, esum, dsum, tmp;
       int *rd;
       int counter, flag;
    FILE *fp;

    fp = fopen_path("iterPeakLocation.txt");


       
       peak = malloc(nms*sizeof(double));

       if(peak == 0)
       {
           error("Failed to allocate auxiliary array");
       }
       
       epeak = malloc(nms*sizeof(double));

       if(epeak == 0)
       {
           error("Failed to allocate auxiliary array");
       }
       
       dsum = 0;
       counter = 0;
       psum = 0;
       esum = 0;
       flag = 0;
       
       tmp = 0.61 - ei;
       tmp = tmp/de;
       int cut = tmp;
      // cut = cut+1;
       

       
       rd = malloc(sizeof(int)*nms);
       randi(rd, nms, 0, nms);
       
       
       

       
       
       for (int k = 0; k < nms; k++)
       {

               peak[k] = 0;
               for (int n = 0; n < cut; n++)
               {
                   if (peak[k] < rho1[n+nsteps*k])
                  // if (peak[k] < rho1[n+nsteps*k])
                   {
                       peak[k] = rho1[n+nsteps*k];
                       epeak[k] = ei+n*de;


                   }
                   
               }

          
       }
       
       for (int k = 0; k < nms; k++)
       {
           psum = psum + peak[rd[k]];
           esum = esum + epeak[rd[k]];
       }
       
       psum /= nms;
       esum /= nms;
       
       for (int k = 0; k < nms; k++) {
           tmp = esum - epeak[k];
           tmp = tmp*tmp;
           dsum += tmp;
       }
       
       dsum /= nms;
       dsum = sqrt(dsum);
       
       
       printf("Peak at %3.6f ± %3.6f\n", esum, dsum);
       fprintf(fp, "Peak at %3.6f ± %3.6f\n", esum, dsum);


       fclose(fp);
       free(rd);
       free(peak);
       free(epeak);
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
	fprintf(fp, "kernel = %d\n", kid);
	fprintf(fp, "file   = %s\n", file);
	fprintf(fp, "prec   = %d\n", prec);
	fprintf(fp, "lambda = %lg\n", lambda);
    fprintf(fp, "lambdap = %lg\n", lambdap);
	fprintf(fp, "sigma  = %lg\n", sigma);
	fprintf(fp, "ei     = %lg\n", ei);
	fprintf(fp, "ef     = %lg\n", ef);
	fprintf(fp, "nsteps = %d\n", nsteps);
	fprintf(fp, "emin   = %lg\n", emin);
	fprintf(fp, "scan   = %lg\n", escan);
	fclose(fp);
}

void allocate_data()
{
	int count;
    int th = tmx/2;

	//count = 4*tmx+2*tmx*nms;
    count = tmx + tmx + tmx*nms + nms*th-nms + tmx-2 + nms*th+nms + nsteps;
	corr = malloc(count*sizeof(mpfr_t)); // tmx
	cov = corr+tmx; //tmx
	mcorr = cov+tmx;    // tmx*nms
    
    mass = mcorr+tmx*nms;   //nms*(T/2-1) = nms*T/2-nms
    bmass = mass+th*nms-nms;      // T/2-1
    deltam = bmass+th-1;         // T/2-1
    smcorr = deltam+tmx/2-1;    // (T/2+1)nms = nms*T/2 + nms
    deltazero = smcorr + nms*th+nms;    // nsteps

	if(corr == NULL)
	{
		error("Failed to allocate MPFR variables");
	}

	for(int i = 0; i < count; i++)
	{
		mpfr_init(corr[i]);
	}
    
 

    
    
    bootrho0 = malloc(nsteps*sizeof(double));

    if(bootrho0 == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    bootdeltarho0 = malloc(nsteps*sizeof(double));

    if(bootdeltarho0 == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    bootrho1 = malloc(nsteps*sizeof(double));

    if(bootrho1 == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    bootdeltarho1 = malloc(nsteps*sizeof(double));

    if(bootdeltarho1 == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    
    rho0 = malloc(nms*nsteps*sizeof(double));
    if(rho0 == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    rho1 = malloc(nms*nsteps*sizeof(double));
    if(rho1 == 0)
    {
        error("Failed to allocate auxiliary array");
    }
    drho0 = malloc(nms*nsteps*sizeof(double));
    if(drho0 == 0)
     {
         error("Failed to allocate auxiliary array");
     }
    drho1 = malloc(nms*nsteps*sizeof(double));
    if(drho1 == 0)
     {
         error("Failed to allocate auxiliary array");
     }

    
    
}

int main(int argc, char *argv[])
{
	memset(file, 0, sizeof(file));
	find_int(argc, argv, "-T", &tmx);
	find_int(argc, argv, "-nms", &nms);
	find_int(argc, argv, "-kernel", &kid);
	find_str(argc, argv, "-file", file);

	if(tmx <= 0 || nms <= 0 || kid == -1 || strlen(file) == 0)
	{
		printf("Usage: %s -T <int> -nms <int> -kernel <int> -file <string> [ options ]\n", argv[0]);
		printf("Options:\n");
		printf("  -T       <int>     temporal extent of correlator\n");
		printf("  -nms     <int>     number of measurements\n");
		printf("  -file    <string>  path for reading the correlator data\n");
		printf("  -lambda  <float>   trade-off parameter (default %1.6f)\n", lambda);
        printf("  -lambdap  <float>   level 1 trade-off parameter (default %1.6f)\n", lambdap);
		printf("  -sigma   <float>   width for smearing function (default %1.6f)\n", sigma);
		printf("  -prec    <int>     working precision (default %d decimal places)\n", prec);
		printf("  -ei      <float>   start of energy range for spectral density (default %1.6f)\n", ei);
		printf("  -ef      <float>   end of energy range for spectral density (default %1.6f)\n", ef);
		printf("  -nsteps  <int>     number of steps used in the energy range (default %d)\n", nsteps);
		printf("  -emin    <float>   lower limit for integration over energy (default %1.6f)\n", emin);
		printf("  -scan    <float>   perform a scan in lambda at the specified energy (not default)\n");
        printf("  -iter  <int>       enter 1 for an extra iteration (not default)\n");
        printf("  -getmass              only compute the effective mass on the bootstrapped sample (not default)\n");
		printf("\n");
		printf("Kernels:\n");
		printf("  0        Gaussian\n");
		printf("  1        Sinc\n");
		printf("  2        i-epsilon (real part)\n");
		printf("  3        i-epsilon (imaginary part)\n");
		return 0;
	}

	find_dbl(argc, argv, "-sigma", &sigma);
	find_dbl(argc, argv, "-lambda", &lambda);
    find_dbl(argc, argv, "-lambdap", &lambdap);
	find_dbl(argc, argv, "-ei", &ei);
	find_dbl(argc, argv, "-ef", &ef);
	find_int(argc, argv, "-nsteps", &nsteps);
	find_int(argc, argv, "-prec", &prec);
	find_dbl(argc, argv, "-emin", &emin);
	find_dbl(argc, argv, "-scan", &escan);
    find_int(argc, argv, "-iter", &iter);
    find_opt(argc, argv, "-getmass", &getmass);

	srand48(1337);
	mpfr_set_default_prec(3.322*prec);

	de = (ef-ei)/nsteps;
	set_params(sigma, lambda, kid, lambdap);
	set_time_parms(1, tmx/2, tmx);

	prepare_path();
	allocate_data();
	read_corr();
    
    if (getmass)
    {
        bootmass();
        fitmass();
        exit(1);
    }
    

    if (iter)
    {
        if(escan)
        {
            scan_lambda();
            scan_lambda_lvl1();
        }
        else
        {
            get_spectral_density();
            peak_bootstrap();
            get_spectral_density_improv();
            peak_bootstrap_iter();
            
        }
    }
    else
    {
        if(escan)
        {
            scan_lambda();
        }
        else
        {
            get_spectral_density();
            peak_bootstrap();
        }
        
    }
    
	

	return 0;
}
