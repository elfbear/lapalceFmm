/*
  test_adap_fmm.c: demonstrates the use of parallel adaptive fmm-laplace and 
  fmm-yukawa packages. 
  Copyright (c) 2012 Bo Zhang, Jingfang Huang, Nikos P. Pitsianis, Xiaobai Sun

  This program is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public Licenses as published by 
  the Free Software Foundation, either version 3 of the Licenses, or any 
  later version. 

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see http://www.gnu.org/licenses/.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <math.h>
#include "adap_fmm.h"
#include "test_adap_fmm.h"

int main(int argc, char **argv)
{
  int nparts, s, accuracy, distribution;

  double beta, *ploc, *pcharge, *pot, *field, *dxx, *dyy, *dzz, *dxy, *dxz, *dyz, elapsed;

  struct timeval tic, toc;

  test_parser(argc, argv, &beta, &nparts, &s, &accuracy, &distribution);

  test_init(distribution, nparts, accuracy, s, beta, &ploc, &pcharge, &pot, &field, 
            &dxx, &dyy, &dzz, &dxy, &dxz, &dyz);

  printf("\n\tPROGRESS\n");
  printf("======================================================\n");
  gettimeofday(&tic, 0);

  adap_fmm_init(accuracy, nparts);

  gettimeofday(&toc, 0);
  elapsed = (double) ((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
  printf("t(FMM_INIT)          | %20.4e\n", elapsed);

  gettimeofday(&tic, 0);

  adap_fmm_graph(nparts, s, beta, ploc, pcharge);
  
  gettimeofday(&toc, 0);
  elapsed = (double) ((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
  printf("t(FMM_GRAPH)         | %20.4e\n", elapsed);

  gettimeofday(&tic,0);

  adap_fmm_compute();
  
  gettimeofday(&toc, 0);
  elapsed = (double) ((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
  printf("t(FMM_COMPUTE)       | %20.4e\n", elapsed);

  gettimeofday(&tic, 0);

  adap_fmm_post(pot, field, dxx, dyy, dzz, dxy, dxz, dyz);
  
  gettimeofday(&toc, 0);
  elapsed = (double) ((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
  printf("t(FMM_POST)          | %20.4e\n", elapsed);

  gettimeofday(&tic, 0);

  adap_fmm_clean ();

  gettimeofday(&toc, 0);
  elapsed = (double) ((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
  printf("t(FMM_CLEAN)         | %20.4e\n", elapsed);
  printf("======================================================\n\n");
 
  test_verify(beta, nparts, ploc, pcharge, pot, field, accuracy, dxx, dyy, dzz, dxy, dxz, dyz);

  test_clean(ploc, pcharge, pot, field, dxx, dyy, dzz, dxy, dxz, dyz);

  return 0;
}

void test_parser(int argc, char ** argv, double *beta, int *nparts, int *s, 
		 int *accuracy, int *distribution)
{
  // Setup default demo parameters: The kernel function is expressed as 
  // k(r) = exp(-beta*r)/r. When beta is zero, the demo is for Laplace kernel.
  // When beta is a positive number, the demo is for Yukawa kernel. 

  // Three types of distributions are provided for the test: (1) particles uniformly
  // distributed inside a box, referred to as CUBE; (2) particles uniformly distributed 
  // over a spherical surface, referred to as SPHERE; and (3) the restriction of type
  // (2) in the first octant, referred to as OCTANT. 

  *beta = 0; 
  *nparts = 1000;
  *s = 80;
  *accuracy = 3;
  *distribution = 1;

  // Parse command line flags when necessary. 
  int i = 0;
  if ( argc > 1 ) {
    for ( i = 1; i < argc; i++ ) {
      if ( argv[i][0] == '-' ) {
	switch ( argv[i][1] ) {
	case 'b':
	  *beta = atof(argv[++i]);
	  break;
	case 'n':
	  *nparts = atoi(argv[++i]);
	  break;
	case 's':
	  *s = atoi(argv[++i]);
	  break;
	case 'a':
	  *accuracy = atoi(argv[++i]);
	  break;
	case 'd':
	  *distribution = atoi(argv[++i]);
	  break;
	default:
	  break;
	}
      }
    }
  }
  return;
}

void test_init(const int distribution, const int nparts, const int accuracy, 
	       const int s, const double beta, double **ploc, double **pcharge, 
	       double **pot, double **field, double **dxx, double **dyy, double **dzz,
	       double **dxy, double **dxz, double **dyz)
{
  // The function completes two tasks: (1) Allocate memory to hold particle locations
  // and the charges carried by them, and to hold the computed potential and field 
  // result; (2) Generate the input data as specified by the distribution variable. 

  char* datatype[] = {"", "  CUBE", "SPHERE", "OCTANT"};

  // Allocate memory to hold particle information and output results. 
  (*ploc) = (double *)calloc(nparts*3, sizeof(double));
  (*pcharge) = (double *)calloc(nparts, sizeof(double));
  (*pot) = (double *)calloc(nparts, sizeof(double));
  (*field) = (double *)calloc(nparts*3, sizeof(double));
  (*dxx) = (double *)calloc(nparts, sizeof(double));
  (*dyy) = (double*)calloc(nparts, sizeof(double));
  (*dzz) = (double *)calloc(nparts, sizeof(double));
  (*dxy) = (double *)calloc(nparts, sizeof(double));
  (*dxz) = (double *)calloc(nparts, sizeof(double));
  (*dyz) = (double *)calloc(nparts, sizeof(double));  
  

  int allocFailure = (*ploc==0) || (*pcharge==0) || (*pot==0) || (*field==0) ||
  (*dxx == 0) || (*dyy == 0) ||(*dzz == 0) ||(*dxy == 0) || (*dxz == 0) || (*dyz == 0);
  if ( allocFailure ) {
    printf("Error in %s, line %d: unable to allocate memory\n",  __FILE__, __LINE__);
    exit(-1);
  }

  // Generate input data set for the demo
  test_data(nparts, distribution, *ploc, *pcharge);

  // Print summary of the demo
  if ( beta > 0 ) {
    printf("\n\tDEMO FMM-YUKAWA RUN\n"
	   "\n\tSETUP\n\n"
	   "======================================================\n"
	   "# OF PARTICLE        | %20d\n"
	   "DISTRIBUTION         | \t\t     %s\n"
	   "S                    | %20d\n"
	   "ACCURACY             | %20d\n"
	   "======================================================\n", 
	   nparts, datatype[distribution], s, accuracy);
  } else if  ( beta == 0 ) {
    printf("\n\tDEMO FMM-LAPLACE RUN\n"
	   "\n\tSETUP\n\n"
	   "======================================================\n"
	   "# OF PARTICLE        | %20d\n"
	   "DISTRIBUTION         | \t\t     %s\n"
	   "S                    | %20d\n"
	   "ACCURACY             | %20d\n"
	   "======================================================\n", 
	   nparts, datatype[distribution], s, accuracy);
  }
  return;
}


void test_verify(const double beta, const int nparts, const double *ploc, 
		 const double *pcharge, const double *pot, const double *field, 
		 const int accuracy,const double *dxx, const double *dyy, 
		 const double *dzz, const double *dxy, const double *dxz, 
		 const double *dyz)
{
  // The verify function computes the potential and field result using pairwise
  // interaction at the first "ncheck" generated particle locations, and compare
  // the accurate results with the ones returned from FMM. 

  const int ncheck = nparts;// (NDIRECT < nparts? NDIRECT: nparts);
  double *dpot = (double *)calloc(ncheck, sizeof(double));
  double *dfield = (double *)calloc(ncheck*3, sizeof(double));
  double *ddxx = (double *)calloc(ncheck, sizeof(double));
  double *ddyy = (double *)calloc(ncheck, sizeof(double));
  double *ddzz = (double *)calloc(ncheck, sizeof(double));
  double *ddxy = (double *)calloc(ncheck, sizeof(double));
  double *ddxz = (double *)calloc(ncheck, sizeof(double));
  double *ddyz = (double *)calloc(ncheck, sizeof(double));
  
  int allocFailure = (dpot==0) || (dfield==0) || (ddxx==0) || (ddyy==0)|| (ddzz==0) ||
                     (ddxy==0) || (ddxz==0) || (ddyz==0);
  if ( allocFailure ) {
    printf("Error in %s, line %d: unable to allocate memory\n", __FILE__, __LINE__);
    exit(-1);
  }

  // Compute accurate results using pairwise interaction. 
  int i;

  if ( beta == 0 ) {
    cilk_for ( i = 0; i < ncheck; i++ ) 
      lapdirect(nparts, ploc, pcharge, i, &dpot[i], &dfield[3*i], 
	        &ddxx[i], &ddyy[i], &ddzz[i], &ddxy[i], &ddxz[i], &ddyz[i]);
  } else {
    cilk_for ( i = 0; i < ncheck; i++)
      yukdirect(beta, nparts, ploc, pcharge, i, &dpot[i], &dfield[3*i],
	   &ddxx[i], &ddyy[i], &ddzz[i], &ddxy[i], &ddxz[i], &ddyz[i]);
  }

  // Compute the numerical error, including L2 and L-infty norms.
  double salg = 0; 
  double salg2 = 0; 
  double salg3 = 0;
  double salg4 = 0;
  double salg5 =0;
  double salg6 = 0;
  double salg7 = 0;
  double salg8 = 0;
  
  double stot = 0; 
  double stot2 = 0; 
  double stot3 = 0;
  double stot4=0;
  double stot5=0;
  double stot6 =0;
  double stot7 = 0;
  double stot8 =0;
  double errmax = 0; 
  double salg21=0, salg22=0, salg23=0;
  double stot21 =0, stot22=0, stot23=0;
  for ( i = 0; i < ncheck; i++ ) { 
  
    salg += pow(dpot[i] - pot[i], 2);
    stot += pow(pot[i], 2); 
	
    int j = 3*i;      
    //salg2 += pow(dfield[j] - field[j], 2) + pow(dfield[j+1] - field[j+1], 2) +
    //         pow(dfield[j+2] - field[j+2], 2);      
    //stot2 += pow(dfield[j], 2 ) + pow(dfield[j+1], 2) + pow (dfield[j+2], 2);

    salg21 += pow(dfield[j] - field[j], 2);
    salg22 += pow(dfield[j+1] - field[j+1], 2);
    salg23 += pow(dfield[j+2] - field[j+2], 2);
    
    stot21 += pow(dfield[j], 2 );
    stot22 += pow(dfield[j+1], 2 );
    stot23 += pow(dfield[j+2], 2 );

    salg3 += pow(ddxx[i]-dxx[i], 2) ; // dxx
    salg4 += pow(ddyy[i]-dyy[i],2);   // dyy
    salg5 += pow(ddzz[i]-dzz[i],2);   // dzz
    salg6 += pow(ddxy[i]-dxy[i], 2); // dxy 
    salg7 += pow(ddxz[i]-dxz[i],2);
    salg8 += pow(ddyz[i]-dyz[i],2);
	
    stot3 += pow(ddxx[i], 2);
    stot4 += pow(ddyy[i], 2);
    stot5 += pow(ddzz[i], 2);
    stot6 += pow(ddxy[i], 2);// 
    stot7 += pow(ddxz[i], 2);
    stot8 += pow(ddyz[i], 2);
	
    
    if ( fabs(dpot[i] - pot[i]) >= errmax )
      errmax = fabs(dpot[i] - pot[i]);      
  }

  printf("\n\tERROR ANALYSIS\n\n"
	 "======================================================\n"
	 "ALGORITHM L2 ERROR   | %20.4e\n"
	 "ALGORITHM MAX ERROR  | %20.4e\n"
	 "L2 ERROR OF FIELD    | %20.4e\n  %20.2e\n %20.2e\n"
	 "L2 ERROR OF 2nd Order Derivative | %20.4e\n"
	 "norm(dxx-ddxx) |%20.4e\n"
	 "norm(ddxx)| %20.4e\n"
	 "======================================================\n",	
	 sqrt(salg/stot), errmax, sqrt(salg21/stot21),sqrt(salg22/stot22), sqrt(salg23/stot23), sqrt(salg3/stot3),
	 sqrt(salg3), sqrt(stot3) );
  printf("norm(dyy-ddyy) = %20.4e\n norm(dzz-ddzz) =%20.4e\n norm(dxy-ddxy) =%20.4e\n norm(dxz-ddxz) = %20.4e\n norm(dyz-ddyz) = %20.4e\n",
           sqrt(salg4), sqrt(salg5), sqrt(salg6),sqrt(salg7),sqrt(salg8) );
  printf(" norm(ddyy) = %20.4e\n  norm(ddzz) = %20.4e\n  norm(ddxy) = %20.4e\n ",
           sqrt(stot4), sqrt(stot5), sqrt(stot6) );
  

  int IfPass = (sqrt(salg/stot) <= pow(10.0, -accuracy))&& 
    (sqrt(salg2/stot2) <= pow(10.0, -accuracy+1)) &&
	(sqrt(salg3/stot3) <= pow(10.0, -accuracy+2)) &&(sqrt(salg4/stot4) <= pow(10.0, -accuracy+2))
	&&(sqrt(salg5/stot5) <= pow(10.0, -accuracy+2))&&(sqrt(salg6/stot6) <= pow(10.0, -accuracy+2));

  printf("norm(dxx-ddxx)/norm(ddxx) =  %20.4e\n", sqrt(salg3/stot3));
  printf("norm(dyy-ddyy)/norm(dyy)  =  %20.4e\n", sqrt(salg4/stot4));
  printf("norm(dzz-ddzz)/norm(dzz)  =  %20.4e\n", sqrt(salg5/stot5));
  printf("norm(dxy-ddxy)/norm(dxy)  =  %20.4e\n", sqrt(salg6/stot6));
  printf("norm(dxz-ddxz)/norm(dxz)  =  %20.4e\n", sqrt(salg7/stot7));
  printf("norm(dyz-ddyz)/norm(dyz) =  %20.4e\n", sqrt(salg8/stot8));

  free(dpot);
  free(dfield);
  free(ddxx);
  free(ddyy);
  free(ddzz);
  free(ddxy);
  free(ddxz);
  free(ddyz);
  return;
}

void test_clean(double *ploc, double *pcharge, double *pot, double *field,
		  double *dxx, double *dyy, double *dzz, double *dxy, double *dxz, double *dyz )
{
  free(ploc);
  free(pcharge);
  free(pot);
  free(field);
  free(dxx);
  free(dyy);
  free(dzz);
  free(dxy);
  free(dxz);
  free(dyz);
}

void test_data(int nparts, int distribution, double *ploc, double *pcharge)
{
  if ( distribution == 1 ) {
    // generate data randomly distributed in a unit cube
    int i; 
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
      pcharge[i] = 1.0;//(double) rand()/RAND_MAX - 0.5;
      ploc[j] = 0.1*((double) rand()/RAND_MAX - 0.5);
      ploc[j+1] = (double) rand()/RAND_MAX - 0.5;
      ploc[j+2] = (double) rand()/RAND_MAX - 0.5;
    }
  } else if ( distribution == 2 ) {
    // generate data randomly distributed over a spherical surface
    int i;
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
      pcharge[i] = (double) rand()/RAND_MAX - 0.5;
      double theta = ((double) rand()/RAND_MAX - 0.5)*M_PI;
      double phi = ((double) rand()/RAND_MAX - 0.5)*M_PI;
      ploc[j] = sin(theta)*sin(phi);
      ploc[j+1] = sin(theta)*cos(phi);
      ploc[j+2] = cos(theta);
    }
  } else if ( distribution == 3 ) {
    // generate data randomly distributed over a spherical surface that is in the first octant
    int i;
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
      pcharge[i] = (double) rand()/RAND_MAX - 0.5;
      double theta = ((double) rand()/RAND_MAX)*M_PI_2;
      double phi = ((double) rand()/RAND_MAX)*M_PI_2;
      ploc[j] = sin(theta)*sin(phi);
      ploc[j+1] = sin(theta)*cos(phi);
      ploc[j+2] = cos(theta);
    }
  }

  return;
}

void lapdirect(const int nparts, const double *ploc, const double *pcharge, 
	       const int i, double *pot, double *field, double *dxx, double *dyy, 
	       double *dzz, double *dxy, double *dxz, double *dyz)
{
  int j, j3, i3;
  double rx, ry, rz, rr, rdis, rmul;
  *pot = 0;
  field[0] = 0;
  field[1] = 0;
  field[2] = 0;
  i3 = i*3;
  
  *dxx = 0; *dyy = 0; *dzz = 0;
  *dxy = 0; *dxz = 0; *dyz = 0;

  for ( j = 0; j < i; j++ ) {
    j3 = j*3;
    rx = ploc[i3]   - ploc[j3];
    ry = ploc[i3+1] - ploc[j3+1];
    rz = ploc[i3+2] - ploc[j3+2];
    rr = rx*rx+ry*ry+rz*rz;
    rdis = sqrt(rr);

    *pot += pcharge[j]/rdis;
    rmul = pcharge[j]/(rdis*rr);
    field[0] += rmul*rx;
    field[1] += rmul*ry;
    field[2] += rmul*rz;
	
    *dxx += rmul * ( 2*rx*rx - ry*ry - rz*rz )/rr;
    *dyy += rmul * ( 2*ry*ry - rx*rx - rz*rz )/rr;
    *dzz += rmul * ( 2*rz*rz - rx*rx - ry*ry )/rr;
    *dxy += rmul * ( 3*rx*ry )/rr;
    *dxz += rmul * ( 3*rx*rz )/rr;
    *dyz += rmul * ( 3*ry*rz )/rr; 
  }

  for ( j = i+1; j < nparts; j++ ) {
    j3 = j*3;
    rx = ploc[i3]   - ploc[j3];
    ry = ploc[i3+1] - ploc[j3+1];
    rz = ploc[i3+2] - ploc[j3+2];
    rr = rx*rx+ry*ry+rz*rz;
    rdis = sqrt(rr);

    *pot += pcharge[j]/rdis;
    rmul = pcharge[j]/(rdis*rr);
    field[0] += rmul*rx;
    field[1] += rmul*ry;
    field[2] += rmul*rz;
  
    *dxx += rmul * ( 2*rx*rx - ry*ry - rz*rz )/rr;
    *dyy += rmul * ( 2*ry*ry - rx*rx - rz*rz )/rr;
    *dzz += rmul * ( 2*rz*rz - rx*rx - ry*ry )/rr;
    *dxy += rmul * ( 3*rx*ry )/rr;
    *dxz += rmul * ( 3*rx*rz )/rr;
    *dyz += rmul * ( 3*ry*rz )/rr;
  }
}

void yukdirect(const double beta, const int nparts, const double *ploc, 
	       const double *pcharge, const int i, double *pot, double *field,
		   double *dxx, double *dyy, double *dzz, double *dxy, double *dxz, double *dyz)
{
  int j, j3, i3;
  double rx, ry, rz, rr, rdis, expr, term1, term2;
  
  *pot = 0;
  field[0] = 0;
  field[1] = 0;
  field[2] = 0;
  
  *dxx = 0; *dyy = 0; *dzz = 0;
  *dxy = 0; *dxz = 0; *dyz = 0;
  
  i3 = i*3;
  for ( j = 0; j < i; j++ ) {
    j3 = 3*j;
    rx = ploc[i3] - ploc[j3];
    ry = ploc[i3+1] - ploc[j3+1];
    rz = ploc[i3+2] - ploc[j3+2];
    rr = rx*rx + ry*ry + rz*rz;
    rdis = sqrt(rr)*beta;
    expr = exp(-rdis);
    term1 = pcharge[j]/rdis*expr*M_PI_2;
    *pot += term1;
    
    term2 = -term1*(1+rdis)/rr;
    field[0] += rx*term2;
    field[1] += ry*term2;
    field[2] += rz*term2;
  }
  
  for ( j = i+1; j < nparts; j++ ) {
    j3 = 3*j;
    rx = ploc[i3] - ploc[j3];
    ry = ploc[i3+1] - ploc[j3+1];
    rz = ploc[i3+2] - ploc[j3+2];
    rr = rx*rx + ry*ry + rz*rz;
    rdis = sqrt(rr)*beta;
    expr = exp(-rdis);
    term1 = pcharge[j]/rdis*expr*M_PI_2;
    *pot += term1;
    
    term2 = -term1*(1+rdis)/rr;
    field[0] += rx*term2;
    field[1] += ry*term2;
    field[2] += rz*term2;
  }
}

