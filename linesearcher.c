#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_interp.h>
#define max(a,b) (((a)>(b))?(a):(b))


// This code is a variation of ARES code. The output of this code is a list of lines found in a given interval for a spectra.
// As in ARES it is better if the spectra has a preliminary normalization. And the user should be carefull with spectral gaps, or irregulations 
// such as cosmic rays.
// 
// The input file 'input.search' is also used. Altough the parameters readlinedat, space, miniline are not used by the code.
//
// 
//
//
//

// ARES v1.1 Automatic Routine for Equivalent widths of Spectra
// method and implemention :Sergio G. Sousa (CAUP & CAAUL) sousasag@astro.up.pt
// supervised by Nuno C. Santos and Mario J. P. F. G. Monteiro
//
//Compilations instructions:
//
// gcc -o runsearch linesearcher.c -L/usr/local/cfitsio/lib/ -I/usr/local/cfitsio/include/ -lcfitsio -lgsl -lgslcblas -lm
// this works if cfitsio is in the directory: /usr/local/cfitsio



/*compilar com: gcc -o runsearch linesearcher.c -L/usr/local/cfitsio/lib/ -I/usr/local/cfitsio/include/ -lcfitsio -lgsl -lgslcblas -lm */

// É necessario ter a gnu science library instalada, gsl
// Para graficos é necessario estar instalado o plotutils, os plots sao gerados atraves da linha de comandos usando o comando graph.
// Este programa nao foi criado para mostrar plots...


long file_lines(char *);
void arraysubcp(double *, double *, long, long);
void plotxy(double *, double *, long, double, double);
void plotxyover(double *, double *, long, double *, double *, long, double, double);
void plotxyover2(double *, double *, long, double *, double *, long, double, double);
void plotxyover3(double *, double *, long, double *, double *, long, double *, double *, long, double, double);
void poly_fitn(double *, double *, double *, long, long, double *);
void continuum_det5 (double *, double *, double *, long, double *, double, double, int);
void deriv(double *, double *, double *, long);
void smooth(double *, long, int, double *);
void zeroscenterfind(double *, double *, double *, double *, double *, long, double *, long, double *, long, long *, long *, double);
double maxele_vec(double *, long);
void fitngauss(double *, double *, double *, long,  double *, int, int *);

int main()
{

/*Declaraçao de variaveis*/

    FILE * pFile;
    FILE * pFile2;
    FILE * pFile3;
    char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
    int status = 0, status2;   /* CFITSIO status value MUST be initialized to zero! */
    int bitpix, naxis, ii, anynul;
    long naxes[2] = {1,1}, fpixel[2] = {1,1},npoints;
    double *pixels, *xpixels, caga, cdelta1, crval1;
    double *atesttot;
    char format[20], hdformat[20];
    int single = 0, hdupos, nkeys;
/* As minhas variaveis*/
    char filetest[200], fileleitura[200], fileout[200];
    double lambdai, lambdaf, smoothder, tree, distlinha, miniline,snratio;
    double space;
    int plots_flag, plots_stop,nchar;
    FILE * fopt;
    char str[200],str2[200];


/* FIM DE Declaracao de variaveis*/


/* leitura das opcoes , mine.opt file read*/

	fopt = fopen("input.search","rt");

	system("clear");
	printf("Input Parameters:\n\n");

	fgets (str , 200 , fopt);
	char *pch;
	pch = strtok (str,"' ");
	pch = strtok (NULL, "' ");
	strcpy (filetest,pch);
	printf("specfits: %s\n",filetest);


//	fgets (str , 200 , fopt);
//	pch = strtok (str,"' ");
//	pch = strtok (NULL, "' ");
//	strcpy (fileleitura,pch);
//	printf("readlinedat: %s\n",fileleitura);


	fgets (str , 200 , fopt);
	pch = strtok (str,"' ");
	pch = strtok (NULL, "' ");
	strcpy (fileout,pch);
	printf("fileout: %s\n",fileout);


	fgets (str , 200 , fopt);
	pch = strtok (str,"=");
	pch = strtok (NULL, " ");
	lambdai=atof(pch);
	printf("lambdai: %6.1f\n",lambdai);


	fgets (str , 200 , fopt);
	pch = strtok (str,"=");
	pch = strtok (NULL, " ");
	lambdaf=atof(pch);
	printf("lambdaf: %6.1f\n",lambdaf);


	fgets (str , 200 , fopt);
	pch = strtok (str,"=");
	pch = strtok (NULL, " ");
	smoothder=atof(pch);
	printf("smoothder: %3.1f\n",smoothder);


//	fgets (str , 200 , fopt);
//	pch = strtok (str,"=");
//	pch = strtok (NULL, " ");
//	space=atol(pch);
//	printf("space: %.2f \n",(float) space);


	fgets (str , 200 , fopt);
	pch = strtok (str,"=");
	pch = strtok (NULL, " ");
	tree=atof(pch);
	printf("tree: %5.3f\n",tree);


//	fgets (str , 200 , fopt);
//	pch = strtok (str,"=");
//	pch = strtok (NULL, " ");
//	snratio=atof(pch);
//	printf("snratio: %8.1f\n",snratio);


	fgets (str , 200 , fopt);
	pch = strtok (str,"=");
	pch = strtok (NULL, " ");
	distlinha=atof(pch);
	printf("lineresol: %4.3f\n",distlinha);


//	fgets (str , 200 , fopt);
//	pch = strtok (str,"=");
//	pch = strtok (NULL, " ");
//	miniline=atof(pch);
//	printf("miniline: %3.1f\n",miniline);


//	fgets (str , 200 , fopt);
//	pch = strtok (str,"=");
//	pch = strtok (NULL, " ");
//	plots_flag=atoi(pch);
//	printf("plots_flag: %i\n",plots_flag);

	fclose (fopt);

	plots_stop=0;
	if (plots_flag == 2)
	{
		plots_flag=1;
		plots_stop=1;
	}


    if (!fits_open_file(&fptr, filetest, READONLY, &status))
    {
        if (!fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) )
        {
          if (naxis > 2 || naxis == 0)
             printf("Error: only 1D or 2D images are supported\n");
          else
          {
            /* get memory for 1 row */
            pixels = (double *) malloc(naxes[0] * sizeof(double));
	    xpixels = (double *) malloc(naxes[0] * sizeof(double));
            if (pixels == NULL) {
                printf("Memory allocation error\n");
                return(1);
            }

            if (bitpix > 0) {  /* set the default output format string */
               strcpy(hdformat, " %7d");
               strcpy(format,   " %7.0f");
            } else {
               strcpy(hdformat, " %15d");
               strcpy(format,   " %15.5f");
            }

            printf("\n");          /* print column header */
            status=fits_read_key(fptr, TDOUBLE, "CDELT1", &caga, card, &status);
	    cdelta1=caga;
            status=fits_read_key(fptr, TDOUBLE, "CRVAL1", &caga, card, &status);
	    crval1=caga;
	    npoints=naxes[0];

            for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--)
            {
               if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL,
                    pixels, NULL, &status) )  /* read row of pixels */
                  break;  /* jump out of loop on error */
	       for (ii = 0; ii < naxes[0]; ii++)
               {
		xpixels[ii]= (double)(ii*cdelta1+crval1);
	       }
            }
          }
        }
        fits_close_file(fptr, &status);

	double restest[2];
	restest[0]=1./cdelta1;
	restest[1]=-crval1/cdelta1;
	double testA=5000.;
	double test=restest[0]*testA+restest[1];
	long testi=(long)test;

	/* leitura do ficheiro das riscas a medir laboratory.dat*/
	
	long nl;
//	nl=file_lines(fileleitura);
//	double lab[3][nl-2];
//	char labstr[nl-2][200];
//	pFile = fopen (fileleitura,"rt");
//	fgets (str , 200, pFile);
//	fgets (str , 200, pFile);
//	int il;
//	for (il=0; il < nl-3; il++)
//	{
//		fgets (str , 200, pFile);
//		strcpy (labstr[il],str);
//		pch = strtok (str," ");
//		char *stopstring;
//		lab[0][il]=strtod(pch,&stopstring);
//	}
//	fclose (pFile);

	

// Abrir o ficheiro para escrita dos resultados

	pFile2 = fopen (fileout,"wt");
	pFile3 = fopen ("logsearch.txt","wt");

	int i=0,fgh;


//	espaços de 8AA (4 para cada lado) para encontrar riscas nos 6AA (3 para cada lado)interiores
	space=4.;
	double spacel=space-1.;
	/*truque para tornar rapidos o acesso aos pontos*/
//	long spacepix=(long)500*space;
	long spacepix=(long)space/cdelta1;

	double xfim=xpixels[npoints-1]-1, xini=xpixels[0]+1;
	long ninter=(xfim-xini)/(2*spacel)+1;
	double linhasvec[ninter];
	
	for (i=0; i<ninter; i++)
		linhasvec[i]=xini+3.+i*(2.*spacel);
	
		
//	printf("%f  %f\n", xini,xfim);
//	printf("%f %f %f\n", linhasvec[0]-3.,linhasvec[0],linhasvec[0]+3.);
//	printf("%f %f %f\n", linhasvec[1]-3.,linhasvec[1],linhasvec[1]+3.);
//	printf("%f %f %f\n", linhasvec[ninter-1]-3.,linhasvec[ninter-1],linhasvec[ninter-1]+3.);
//		int pausav;
//		printf("\n ENTER to continue\n");
//		scanf("%i", &pausav);
	
	

//    for (fgh=0;fgh<nl-3;fgh++)
      for (fgh=0;fgh<ninter;fgh++)
    {

//	float linha=lab[0][fgh];
	float linha=linhasvec[fgh];

      if (lambdai < xpixels[0]) lambdai=xpixels[0];
      if (lambdaf > xpixels[npoints-1]) lambdaf=xpixels[npoints-1];

if ( (linha > lambdai+2*space) && (linha < lambdaf-2*space) )

       {
	printf("\n\nline nº %i searching for line %.2f in the interval [%.2f,%.2f]  \n\n",fgh+1, linha, lambdai, lambdaf);
	fprintf(pFile3,"\n\nline nº %i searching for line %.2f in the interval [%.2f,%.2f]  \n\n",fgh+1, linha, lambdai, lambdaf);
	long nctest=(long) (restest[0]*linha+restest[1]);
	nctest++;
	int nx1test, nx2test;
	if (nctest < spacepix) {nx1test = 0;} else {nx1test=nctest-spacepix;}
	if (nctest > npoints-spacepix) {nx2test = npoints-1;} else {nx2test=nctest+spacepix;}
	double retarr[nx2test-nx1test];
	arraysubcp(retarr, pixels,nx1test,nx2test );
	int t;
	double xltest[nx2test-nx1test], atest[nx2test-nx1test];
	arraysubcp(xltest, xpixels,nx1test,nx2test );
	arraysubcp(atest ,  pixels,nx1test,nx2test );

	double x2[2],y2[2];
	x2[0]=x2[1]=linha;
	y2[0]=0.;
	y2[1]=.2;

// encontrar o continuum

	long n1cont= (long) (nctest - space/cdelta1);
	long n2cont= (long) (nctest + space/cdelta1);
	long nx=n2cont-n1cont, nxt=nx2test-nx1test;
	double x[nx],y[nx], atesttotnorm[nxt], ynorm[nx],ynorm2[nxt];
	arraysubcp(x, xpixels,n1cont,n2cont );
	arraysubcp(y,  pixels,n1cont,n2cont );
	double res[4];
	
	continuum_det5(x,y,ynorm,nx,res,tree,snratio,plots_flag);



	for (i=0; i<nxt; i++)
		atesttotnorm[i]=atest[i]/(res[0]+res[1]*xltest[i]+res[2]*xltest[i]*xltest[i]+res[3]*xltest[i]*xltest[i]*xltest[i]);
	for (i=0; i<nx; i++)
		y[i]=y[i]/(res[0]+res[1]*x[i]+res[2]*x[i]*x[i]+res[3]*x[i]*x[i]*x[i]);




	double x3[4],y3[4];
	x3[0]=x3[1]=linha;
	x3[2]=linha-space;
	x3[3]=linha+space;

	y3[0]=0.;
	y3[1]=1.;
	y3[2]=1.;
	y3[3]=1.;

//	calculo das derivadas

	double dy[nx],d2y[nx],d3y[nx];
	deriv(x,y,dy,nx);
	deriv(x,dy,d2y,nx);
	deriv(x,d2y,d3y,nx);

	int xind1=0,xind2=nx-1,hjk;
	float klo=0.1;
//	for (hjk=0; hjk < nx; hjk++){
//		if ( (y[hjk] > tree) && (x[hjk]-(linha-klo) > x[xind1] - (linha-klo)) && (x[hjk] - (linha-klo) < 0) )
//			xind1=hjk;
//		if ( (y[hjk] > tree) && (x[hjk]-(linha+klo) < x[xind2] - (linha+klo)) && (x[hjk] - (linha+klo) > 0) )
//			xind2=hjk;
//	}

	xind1=0;
	xind2=nx-1;


	int nlin=xind2-xind1;
	double xlin[nlin], iylin[nlin], ylin[nlin], dylin[nlin], ddylin[nlin], tmp[nlin];
	double ylincaga[nx], dylincaga[nx], ddylincaga[nx], tmpcaga[nx];
	printf("\nPASSOU AQUI 222\n");
		
	arraysubcp(xlin, x,xind1,xind2 );
	arraysubcp(iylin, y,xind1,xind2 );

	printf("\nPASSOU AQUI 333\n");

//	deriv(xlin,iylin,ylin,nlin);
//		deriv(xlin,iylin,tmp,nlin);
		deriv(x,y,tmpcaga,nx);
		smooth(tmpcaga, nx, (int)smoothder, ylincaga);
		arraysubcp(ylin, ylincaga,xind1,xind2 );

		deriv(x,ylincaga,tmpcaga,nx);
		smooth(tmpcaga, nx, (int)smoothder, dylincaga);
		arraysubcp(dylin, dylincaga,xind1,xind2 );

		deriv(x,dylincaga,tmpcaga,nx);
		smooth(tmpcaga, nx, (int)smoothder, ddylincaga);
		arraysubcp(ddylin, ddylincaga,xind1,xind2 );
	

//	deriv(xlin,ylin,tmp,nlin);
//	smooth(tmp, nlin, (int)smoothder, dylin);
//	deriv(xlin,dylin,tmp,nlin);
//	smooth(tmp, nlin, (int)smoothder, ddylin);


//	if (plots_flag == 1)
//	{
//		plotxy(xlin,ylin,nlin,linha-1.,linha+1.);
//		plotxy(xlin,dylin,nlin,linha-1.,linha+1.);
//		plotxy(xlin,ddylin,nlin,linha-1.,linha+1.);
//	}	

//	procura das riscas que ha a volta da risca que queremos

	double cont[nlin], zeros[nlin], tutezeros[nlin][2];
	long ncont=nlin, nzeros=nlin, ntutezeros=nlin, ncenter=nlin, center[nlin];

	zeroscenterfind(xlin, ylin, iylin, dylin,ddylin,nlin, zeros, nzeros, cont, ncont, center, &ncenter , tree);

//	calculo, interpolacao da posicao das riscas no espectro

	if (center[0] != -1 & ncenter != 0)
		{
		double xlinhas[ncenter], ylinhas[ncenter];
		int i1, i2;
		for (i=0; i<ncenter; i++)
			{
			i1=center[i];
			xlinhas[i]= ( -ddylin[center[i]-1] + ( ddylin[center[i]] - ddylin[center[i]-1] )/( xlin[center[i]]-xlin[center[i]-1] ) * xlin[center[i]] )  / ( ( ddylin[center[i]] - ddylin[center[i]-1] )/(xlin[center[i]]-xlin[center[i]-1]) );
			ylinhas[i]= ( iylin[center[i]] - iylin[center[i]-1] )/( xlin[center[i]] -xlin[center[i]-1] ) * xlinhas[i] + iylin[center[i]-1] - ( iylin[center[i]]- iylin[center[i]-1] )/( xlin[center[i]] -xlin[center[i]-1]) * xlin[center[i]] ;
			}

		printf("\n LINES FOUND TO FIT \n");
		for (i=0; i<ncenter; i++)	printf("%.2f ", xlinhas[i]);
		printf("\n");

		fprintf(pFile3,"\n LINES FOUND TO FIT \n");		
		for (i=0; i<ncenter; i++)	fprintf(pFile3,"%.2f ", xlinhas[i]);
		fprintf(pFile3,"\n");

		double xvec2[ncenter], yvec2[ncenter];
		int nvec2,j;
		xvec2[0]=xlinhas[0];
		yvec2[0]=ylinhas[0];
		j=0;
		for(i=1;i<ncenter;i++)
		{
			if (fabs(xvec2[j]-xlinhas[i]) < distlinha )
			{
			xvec2[j]=(xvec2[j]+xlinhas[i])/2.;
			yvec2[j]=(yvec2[j]+ylinhas[i])/2.;
			}
			else
			{
			j++;
			xvec2[j]=xlinhas[i];
			yvec2[j]=ylinhas[i];
			}
		}
		nvec2=j+1;

		printf("\n RESAMPLING \n");
		for (i=0; i<nvec2; i++)	printf("%.2f ", xvec2[i]);
		printf("\n");

		fprintf(pFile3,"\n RESAMPLING \n");
		for (i=0; i<nvec2; i++)     fprintf(pFile3,"%.2f ", xvec2[i]);
		fprintf(pFile3,"\n");

		ncenter=nvec2;
		int para=3*ncenter;
		int npara=0;
		double acoef[para];
		for (i=0;i<ncenter;i++)
		{
			acoef[3*npara]=yvec2[i]-1.;
			acoef[3*npara+1]=400.;
			acoef[3*npara+2]=xvec2[i];
			npara++;
		}

		double xfit[nlin], yfit[nlin], sigma[nlin];
		for (i=0;i<nlin;i++)
			{
			xfit[i]=xlin[i];
			yfit[i]=iylin[i]-1.0;
			sigma[i]=1.0;
			}
		printf("\n GUESS COEFS :\n");
		for (i=0;i<para;i+=3)
			printf("acoef[%2i]:  %.5f acoef[%2i]:  %9.5f acoef[%2i]:  %7.2f \n", i, acoef[i]+1., i+1, acoef[i+1], i+2, acoef[i+2]);
		
		fprintf(pFile3,"\n GUESS COEFS :\n");
		for (i=0;i<para;i+=3)
			fprintf(pFile3,"acoef[%2i]:  %.5f acoef[%2i]:  %9.5f acoef[%2i]:  %7.2f \n", i, acoef[i]+1., i+1, acoef[i+1], i+2, acoef[i+2]);


//		fitngauss(xfit,yfit,sigma,nlin,acoef,para,&status2);

		printf("\n FITTED COEFS :\n");
		for (i=0;i<para;i+=3)
			printf("acoef[%2i]:  %.5f acoef[%2i]:  %9.5f acoef[%2i]:  %7.2f \n", i, acoef[i]+1., i+1, acoef[i+1], i+2, acoef[i+2]);

		fprintf(pFile3,"\n FITTED COEFS :\n");
		for (i=0;i<para;i+=3)
			fprintf(pFile3,"acoef[%2i]:  %.5f acoef[%2i]:  %9.5f acoef[%2i]:  %7.2f \n", i, acoef[i]+1., i+1, acoef[i+1], i+2, acoef[i+2]);

		double yfit2[nxt];
		for (i=0;i<nxt;i++)
			{
			yfit2[i]=1.0;
			for (j=0;j<ncenter;j++)
				yfit2[i]+=acoef[j*3]* exp (- acoef[j*3+1] * (xltest[i]-acoef[j*3+2]) * (xltest[i]-acoef[j*3+2]) );
			}

		double medida=0;
		int nmed=0, hj, hjl;

		for (hj=0; hj<ncenter;hj++)
			{
			if ( fabs(linha-acoef[3*hj+2]) < distlinha )
				{
				medida+=acoef[3*hj]*sqrt(3.1415927/acoef[3*hj+1]);
				nmed++;
				hjl=hj;
				}
			}

		medida=medida*(-1000.);
		printf("\n---------------------------\nline result: %.5f \n", linha);
		printf("ew (mA)  :  %.5f \n", medida);
		printf("nfit : %i \n", ncenter);

		fprintf(pFile3,"\n---------------------------\nline result: %.5f \n", linha);
		fprintf(pFile3,"ew (mA)  :  %.5f \n", medida);
		fprintf(pFile3,"nfit : %i \n", ncenter);



		if (nmed == 1)
		{
			printf("line depth : %.5f \n", -acoef[3*hjl]);
			// FWHM para a gaussiana defenida: F(X)=Aexp(-Lambda(x-c)^2) => FWHM=2*sqrt(ln(2)/lambda)
			printf("FWHM : %.5f \n-------------------------\n", 2.*sqrt(log(2)/acoef[3*hjl+1]));
			fprintf(pFile3,"line depth : %.5f \n", -acoef[3*hjl]);
			// FWHM para a gaussiana defenida: F(X)=Aexp(-Lambda(x-c)^2) => FWHM=2*sqrt(ln(2)/lambda)
			fprintf(pFile3,"FWHM : %.5f \n-------------------------\n", 2.*sqrt(log(2)/acoef[3*hjl+1]));
		}

		printf("int 2 status: %i", status2); 
		fprintf(pFile3,"int 2 status: %i", status2); 




//		if (medida > miniline && medida < 500. && status2 == 0)
		if (para > 0)
		{
//				lambda nfit depth FWHM EW coef1 coef2 coef3
		printf("AQUI ENTROU 123xpto");
			for (hjl=0;hjl<para;hjl+=3)
				{
				i=hjl;
//				printf("acoef[%2i]:  %.5f acoef[%2i]:  %9.5f acoef[%2i]:  %7.2f \n", i, acoef[i]+1., i+1, acoef[i+1], i+2, acoef[i+2]);
//				printf("%.5f  %.5f %.5f \n", linha-3, acoef[hjl+2],linha+3);			
				if ( acoef[hjl+2] > linha-spacel && acoef[hjl+2] < linha+spacel)
//					printf("%.5f ENTROU\n", acoef[hjl+2]);
					fprintf(pFile2,"%7.2f  %10.5f\n",acoef[hjl+2], acoef[hjl]+1);
					printf("%7.2f  %10.5f\n",acoef[hjl+2], acoef[hjl]+1);

				}
		}


		if (plots_flag == 1)
		{	
			double xcvec[ncenter], ycvec[ncenter];
			for (i=0;i<ncenter;i++)
			{
				xcvec[i]=acoef[i*3+2];
				ycvec[i]=acoef[i*3]+1.;
			}	

			plotxyover2(xltest,atesttotnorm,nxt,xltest,yfit2,nxt,linha-space,linha+space);
			int pausav;
			printf ("\n\nTo Close the plots, click on it.\n 1-continue to show plots, 0-stop plots\n Make your choise:");
			scanf("%i", &pausav);
			plots_flag=pausav;
		}
		}
		else
		{
		printf("\n line not found\n");
		fprintf(pFile3,"\n line not found\n");
//		int pausav;
//		printf("\n ENTER to continue\n");
//		scanf("%i", &pausav);
		}

       }
    }
	fclose (pFile2);
	fclose (pFile3);
	free(pixels);
    } 

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}


long file_lines(char * file) {
    	FILE * pFile;
	long n = 0;
	char str[200];
	pFile = fopen (file,"rt");
	if (pFile==NULL) perror ("Error opening file");
	else
	{
	while (!feof(pFile)) {
		fgets (str , 200, pFile);
		n++;
		}
    	fclose (pFile);
  	}
return n;
}

void arraysubcp(double retarr[], double arr[],long a, long b){
	long i;
	for (i=0; i < b-a; i++)
	{
		retarr[i]=arr[a+i];
	}
}

void plotxy(double xvec[], double yvec[], long np, double xi, double xf){
	FILE * pFile2;
	long t;
	char str[200];
	pFile2 = fopen ("tmp","wt");
	for(t=0;t<np;t++) fprintf(pFile2," %15.6f  %15.6f \n",xvec[t],yvec[t]);
	fclose (pFile2);
	char *buffer;
	int decimal, sign;
	buffer = fcvt (xi, 0, &decimal, &sign);
	strcpy (str,"graph -T X -x ");
	strcat (str,buffer);
	strcat (str," ");
	buffer = fcvt (xf, 0, &decimal, &sign);
	strcat (str,buffer);
	strcat (str," < tmp");
	system(str);
	system("rm tmp");	
}

void plotxyover(double xvec[], double yvec[],long np, double xvec2[], double yvec2[],long np2,double xi, double xf){
	FILE * pFile;
	FILE * pFile2;
	long t;
	char str[200];
	pFile2 = fopen ("tmp","wt");
	for(t=0;t<np;t++) fprintf(pFile2," %15.6f  %15.6f \n",xvec[t],yvec[t]);
	fclose (pFile2);
	pFile = fopen ("tmp2","wt");
	for(t=0;t<np2;t++) fprintf(pFile," %15.6f  %15.6f \n",xvec2[t],yvec2[t]);
	fclose (pFile);
	char *buffer;
	int decimal, sign;
	buffer = fcvt (xi, 0, &decimal, &sign);
	strcpy (str,"graph -T X -x ");
	strcat (str,buffer);
	strcat (str," ");
	buffer = fcvt (xf, 0, &decimal, &sign);
	strcat (str,buffer);
//	strcat (str," tmp -S 2 -C -m 42 tmp2");
	strcat (str," tmp -S 1 -C -m 43 tmp2");
	system(str);
	system("rm tmp tmp2");	
}

void plotxyover2(double xvec[], double yvec[],long np, double xvec2[], double yvec2[],long np2,double xi, double xf){
	FILE * pFile;
	FILE * pFile2;
	FILE * pFile3;
	long t;
	char str[200];
	pFile2 = fopen ("tmp","wt");
	for(t=0;t<np;t++) fprintf(pFile2," %15.6f  %15.6f \n",xvec[t],yvec[t]);
	fclose (pFile2);
	pFile = fopen ("tmp2","wt");
	for(t=0;t<np2;t++) fprintf(pFile," %15.6f  %15.6f \n",xvec2[t],yvec2[t]);
	fclose (pFile);
	pFile3 = fopen ("tmp3","wt");
	fprintf(pFile3," %15.6f  %15.6f \n",xi,1.);
	fprintf(pFile3," %15.6f  %15.6f \n",xf,1.);
	fclose (pFile3);
	char *buffer;
	int decimal, sign;
	buffer = fcvt (xi, 0, &decimal, &sign);
	strcpy (str,"graph -T X -x ");
	strcat (str,buffer);
	strcat (str," ");
	buffer = fcvt (xf, 0, &decimal, &sign);
	strcat (str,buffer);
//	strcat (str," tmp tmp3 -S 2 -C -m 42 tmp2 ");
	strcat (str," tmp tmp3 -S 1 -C -m 43 tmp2 ");
	system(str);
	system("rm tmp tmp2 tmp3");	
}

void plotxyover3(double xvec[], double yvec[],long np, double xvec2[], double yvec2[], long np2, double xvec3[], double yvec3[], long np3, double xi, double xf){
	FILE * pFile;
	FILE * pFile2;
	FILE * pFile3;
	long t;
	char str[200];
	pFile2 = fopen ("tmp","wt");
	for(t=0;t<np;t++) fprintf(pFile2," %15.6f  %15.6f \n",xvec[t],yvec[t]);
	fclose (pFile2);

	pFile = fopen ("tmp2","wt");
	for(t=0;t<np2;t++) fprintf(pFile," %15.6f  %15.6f \n",xvec2[t],yvec2[t]);
	fclose (pFile);
	
	pFile3 = fopen ("tmp3","wt");
	for(t=0;t<np3;t++) fprintf(pFile," %15.6f  %15.6f \n",xvec3[t],yvec3[t]);
	fclose (pFile3);
	char *buffer;
	int decimal, sign;
	buffer = fcvt (xi, 0, &decimal, &sign);
	strcpy (str,"graph -T X -x ");
	strcat (str,buffer);
	strcat (str," ");
	buffer = fcvt (xf, 0, &decimal, &sign);
	strcat (str,buffer);
//	strcat (str," tmp -S 2 -C -m 42 tmp2 -S 4 -C -m 0 tmp3");
	strcat (str," tmp -S 1 -C -m 43 tmp2 -S 1 -C -m 0 tmp3");
	printf("%s\n",str);
	system(str);
	system("rm tmp tmp2 tmp3");	

}

void poly_fitn(double xvec[], double yvec[], double err[], long n, long ord, double coefs[])
{
  int i, j, k;
  double xi, yi, ei, chisq,xi2;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;
  ord++;
  X = gsl_matrix_alloc (n, ord);
  y = gsl_vector_alloc (n);
  w = gsl_vector_alloc (n);
  c = gsl_vector_alloc (ord);
  cov = gsl_matrix_alloc (ord, ord);
  for (i = 0; i < n; i++)
    {
      xi=xvec[i];
      yi=yvec[i];
      ei=err[i];
      for (j = 0; j < ord; j++)
	{
	xi2=1.0;
	for (k=0; k<j; k++) xi2*=xi;
	gsl_matrix_set (X, i, j, xi2);
	}
      gsl_vector_set (y, i, yi);
      gsl_vector_set (w, i, 1.0/(ei*ei));
    }
  {
    gsl_multifit_linear_workspace * work 
      = gsl_multifit_linear_alloc (n, ord);
    gsl_multifit_wlinear (X, w, y, c, cov,
                          &chisq, work);
    gsl_multifit_linear_free (work);
  }

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  {

    for (j = 0; j < ord; j++)
	{
	coefs[j]=C(j);
	}
  }
}


void deriv(double x[], double y[], double dy[], long n)
{
int i;
gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_interp *interp = gsl_interp_alloc (gsl_interp_cspline, n);
//gsl_interp *interp = gsl_interp_alloc (gsl_interp_akima, n);
gsl_interp_init (interp, x, y, n);
for (i=0; i<n; i++){
	dy[i]=gsl_interp_eval_deriv (interp, x, y,x[i],acc);
	}
gsl_interp_free (interp);
gsl_interp_accel_free (acc);
}


void continuum_det5 (double x[], double y[], double ynorm[], long nxele, double res[], double tree, double snratio, int plots_flag){
int order,i,j;
order=2;
double err[nxele], coefs[order+1];
long nvec;

for(i=0;i<nxele;i++)
	err[i]=1.;
poly_fitn(x,y,err,nxele,order,coefs);

double xi=1.;
for(i=0;i<nxele;i++)
	{
	ynorm[i]=0.;
	xi=1.;
	for (j=0;j<order+1;j++)
		{
		ynorm[i]+=coefs[j]*xi;
		xi*=x[i];
		}
	}
double vecx[nxele],vecy[nxele];
int jk;
for (jk=0; jk<5; jk++){
	nvec=0;
	for (i=0; i<nxele-1; i++){
// testes foram feitos com 0.01, nao deve causar problemas por maior. Puz 0.1
		if (y[i] > ynorm[i]*tree && fabs(y[i]-y[i+1]) < 0.1*y[i]){
			vecx[nvec]=x[i];
			vecy[nvec]=y[i];
			nvec++;
		}
	}
	poly_fitn(vecx,vecy,err,nvec,order,coefs);

	for(i=0;i<nxele;i++)
		{
		ynorm[i]=0.;
		xi=1.;
		for (j=0;j<order+1;j++)
			{
			ynorm[i]+=coefs[j]*xi;
			xi*=x[i];
			}
		}
}

	for (i=0; i < order+1; i++)
	{	
		res[i]=coefs[i];
	}

res[3]=0.;
	if (plots_flag == 1)
		plotxyover3(x,y,nxele,x,ynorm,nxele,vecx,vecy,nvec,x[0],x[nxele-1]);


}


void smooth(double vec[], long n, int w, double svec[])
{
int i,j;
double soma;

if (w%2 != 1)
	w++;
for (i=0; (i < (w-1)/2);i++)
	svec[i]=vec[i];
for (i=(w-1)/2;i<n-((w-1)/2);i++)
{
	soma=0.;
	for (j=i-((w-1)/2); j<=i+((w-1)/2);j++)
		soma+=vec[j];
	svec[i]=soma/w;
}
for (i=n-((w-1)/2); i<n;i++)
	svec[i]=vec[i];
}


void zeroscenterfind(double x[], double y[], double iy[], double dy[], double ddy[], long n, double zeros[], long nzeros, double cont[], long ncont, long center[], long *ncenter, double tree)
{
double zerostot[n], contot[n], tutezerostot[n][2], maxdy;
long ntot=0, nctot=0, ctot=0, i, centertot[n];
int signal=0, signalc=0, signal_ant, signalc_ant;
if (y[0] == abs(y[0]))
	signal=1;
if (ddy[0] == abs(ddy[0]))
	signalc=1;
signal_ant=signal;
signalc_ant=signalc;
maxdy=maxele_vec(dy,n);

//printf("\n n : %i  max dy:  %f10.3  -- 0.05*maxdy: %f10.3 \n", n, maxdy, 0.05*maxdy);

for (i=0; i<n; i++)
{

//	printf("\n  dy[i]:  %10.3f  -- 0.05*maxdy: %10.3f", dy[i], 0.05*maxdy);
	
	signalc=0;
	if ( (float) ddy[i] == fabs( (float) ddy[i]) )	signalc=1;
	if ( (signalc != signalc_ant) && (dy[i] > 0.01*maxdy) && (iy[i] < 0.98) && (ddy[i] < -0.1) )
		{
//		printf("\n signal: %i, x: %6.3f, iy: %6.3f, y: %6.3f, dy: %6.3f, ddy: %6.3f", signalc, x[i], iy[i],y[i],dy[i],ddy[i]);
		centertot[ctot]=i;
		ctot++;
		}
	signal=0;
	if ( (float) y[i] == (float) fabs(y[i]))		signal=1;
	if ( signal != signal_ant)
		{
		tutezerostot[ntot+nctot][0]=i;
		if (iy[i] < 0.98)
			{
			zerostot[ntot]=i;
			if (dy[i] <= 0)		tutezerostot[ntot+nctot][1]=0;
			else			tutezerostot[ntot+nctot][1]=0.5;
			ntot++;
			}
		else
			{
			contot[nctot]=i;
			tutezerostot[ntot+nctot][1]=1.;
			nctot++;
			}
		}
	signal_ant=signal;
	signalc_ant=signalc;
		
}


if (ntot != 0)
	{
	nzeros=ntot;
	arraysubcp(zeros, zerostot,0,ntot);
	}
else
	{
	nzeros=0;
	}

if (ctot != 0)
	{
	*ncenter=ctot;
	for (i=0;i<ctot;i++) 	center[i]=centertot[i];
	}
else
	{
	center[0]=-1;
	*ncenter=0;
	}
}


double maxele_vec(double vec[], long nvec)
{
long i;
double maxi=vec[1+nvec/20];
//double maxi=vec[0];
//printf("\n maxdy : nvec: %i , 1+nvec/10: %i nvec-nvec/10: %i \n", nvec, 1+nvec/10, nvec-nvec/10);
for (i=1+nvec/20; i<nvec-nvec/20; i++)
//for (i=1; i<nvec; i++)
	maxi = max(maxi,vec[i]);
return maxi;
}






