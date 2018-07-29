
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <process.h>
#include <unistd.h>
#include <time.h>
//#include "cplex.h"



/* Change any of these parameters to match your needs */
/* This is a bi-objective MIP example */

//#define n 80          /* # decision variables = nc + nb + ni */
//#define m 80          /* # constraints        */
//#define nc 40		  /* # continuous variables */
//#define nb 0          /* # binary variables   */
//#define ni 40          /* # integer variables  */

//#define seedno 33      /* seed number, change for each problem */


#define MODLUS   2147483647      /*Required parameters for Marse&Roberts RNG*/
#define MULT1   24112
#define MULT2   26143


FILE *testp1;
FILE *testp2;


float rand_(int seed);   



static long zrng[] =
{         1,
 1973272912, 281629770,  20006270,1280689831,2096730329,1933576050,
  913566091, 246780520,1363774876, 604901985,1511192140,1259851944,
  824064364, 150493284, 242708531,  75253171,1964472944,1202299975,
  233217322,1911216000, 726370533, 403498145, 993232223,1103205531,
  762430696,1922803170,1385516923,  76271663, 413682397, 726466604,
  336157058,1432650381,1120463904, 595778810, 877722890,1046574445,
   68911991,2088367019, 748545416, 622401386,2122378830, 640690903,
 1774806513,2132545692,2079249579,  78130110, 852776735,1187867272,
 1351423507,1645973084,1997049139, 922510944,2045512870, 898585771,
  243649545,1004818771, 773686062, 403188473, 372279877,1901633463,
  498067494,2087759558, 493157915, 597104727,1530940798,1814496276,
  536444882,1663153658, 855503735,  67784357,1432404475, 619691088,
  119025595, 880802310, 176192644,1116780070, 277854671,1366580350,
 1142483975,2026948561,1053920743, 786262391,1792203830,1494667770,
 1923011392,1433700034,1244184613,1147297105, 539712780,1545929719,
  190641742,1645390429, 264907697, 620389253,1502074852, 927711160,
  364849192,2049576050, 638580085, 547070247 };




void main(int argc, char **argv)
{
	int seedno = atoi(argv[1]);
	int o = atoi(argv[2]);          // number objectives
	int n = atoi(argv[3]);          // number variables
	int m = atoi(argv[4]);          /* # constraints        */
	int nc = atoi(argv[5]);		    /* # continuous variables */
	int nb = atoi(argv[6]);         /* # binary variables   */
	int ni = atoi(argv[7]);         // number integer variables
	unsigned char *filename[10] = {"testp1.lp", "testp2.lp", "testp3.lp", "testp4.lp", "testp5.lp", "testp6.lp", "testp7.lp", "testp8.lp", "testp9.lp", "testp10.lp"};
	
	srand(seedno);
	
	printf("seedno: %d\n",seedno);
/*	printf("%c\n", filename[5]);*/
	
/*	exit(0);*/
	
	int  i, j, ccount=0, upper=0, lower=0;
	double obj_coef1[n], obj_coef2[n], coef[m][n], rhs[m];


	//CPXLPptr   lpclone=NULL;
	
/*	exit(0);*/
	
	/*** rand_omly assign technological coefficients ***/

    upper = 2100;
    lower = -100;
	for(j=0;j<m;j++)
		for(i=0;i<n;i++)
			coef[j][i] = (double)((rand() % (upper - lower + 1)) + lower)/100.; //-1 + rand_(rand()+6)*(21);

/*    exit(0);*/
	
	/*** rand_omly assign RHS    ***/

    upper = 10000;
    lower = 5000;
	for(j=0;j<m;j++)
		rhs[j] = (double)((rand() % (upper - lower + 1)) + lower)/100.; //50 + rand_(rand()+7)*100;
    
    for(int k = 0; k < o; k++)
    {
        ccount = 0;
		if (( testp1 = fopen(filename[k],"w+"))==NULL)
			  {
			  printf("could not open file: %s\n",filename[k]);
			  exit(1);
			  }

	

		/*** rand_omly assign objective coefficients  ***/

        upper = 2000;
        lower = -1000;
		for(i=0; i<nc;i++)   /* continuous var. */
		{
			obj_coef1[i] = (double)((rand() % (upper - lower + 1)) + lower)/100.; //-10 + rand_(rand() + 8*k)*20;
/*			obj_coef2[i] = -10 + rand_(rand()+1)*20;*/

		}

        upper = 40000;
        lower = -20000;
		for(i=nc; i<nc+nb; i++)  /* binary var. */
		{
			obj_coef1[i] = (double)((rand() % (upper - lower + 1)) + lower)/100.; //-200 + rand_(rand()+2 + 8*k)*400;
/*			obj_coef2[i] = -200 + rand_(rand()+3)*400;*/

		}

        upper = 10000;
        lower = -5000;
		for(i=nc+nb; i<n; i++)  /* integer var. */
		{
			obj_coef1[i] = (double)((rand() % (upper - lower + 1)) + lower)/100.; //-50 + rand_(rand()+4 + 8*k)*100;
/*			obj_coef2[i] = -50 + rand_(rand()+5)*100;*/

		}


		/*** write the model1  ***/


	    /********************************************************/
	    /* Objective      
	    /********************************************************/

	    fprintf(testp1,"MAXIMIZE\n");
	    fprintf(testp1,"OBJ: ");
	    for(i=0;i<n;i++)
		    fprintf(testp1," %+.2fx%d",obj_coef1[i],i);
	
       

	    /*********************************************************/
	    /* Constraints
	    /*********************************************************/


	    fprintf(testp1,"\nSUBJECT TO");

	    for(j=0;j<m;j++)
	    {
	
		    fprintf(testp1,"\nC%d: ",ccount);
			    ccount++;
		
		    for(i=0;i<n;i++)
			    fprintf(testp1," %+.2fx%d",coef[j][i],i); 


		    fprintf(testp1," + s%d = %.2f",j, rhs[j]);  /* add slack */	
		
	    }

	    /*********************************************************/
	    /* Bounds
	    /*********************************************************/

	    fprintf(testp1,"\nBinary");
	    for(i=nc; i< nc+nb; i++)
		    fprintf(testp1,"\nx%d",i);

	    fprintf(testp1,"\nGeneral");
	    for(i=nc+nb; i< n; i++)
		    fprintf(testp1,"\nx%d",i);

	
	    fprintf(testp1,"\nEND");


		fclose(testp1);
    }
}



/***********************************************************************/





/**************************************************************/
/*This block is used for M&R RNG*/
float rand_(int seed)
{long zi, lowprd, hi31;

  zi      = zrng[seed];
  lowprd  = (zi & 65535)*MULT1;
  hi31    = (zi >> 16)*MULT1 + (lowprd >> 16);
  zi      = ((lowprd & 65535) - MODLUS) + ((hi31 & 32767) << 16) + (hi31 >> 15);
  if (zi<0) zi += MODLUS;
  lowprd  = (zi & 65535)*MULT2;
  hi31    = (zi >> 16)*MULT2 + (lowprd >> 16);
  zi      = ((lowprd & 65535)-MODLUS)+((hi31 & 32767) << 16) + (hi31 >> 15);
  if (zi<0) zi += MODLUS;
  zrng[seed] = zi;
  return ((zi >> 7 | 1) + 1)/ 16777216.0;
  }

/*******************************************************************/

