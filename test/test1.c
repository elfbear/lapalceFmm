#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <math.h>
#include "adap_fmm.h"
#include "RPY.h"

int main(int argc, char **argv)
{
  int nparts, s, accuracy, distribution; 
  double *ploc, *pcharge, *pot, *field, *charge, *RPY;
  double *dxx, *dyy, *dzz, *dxy, *dyz, *dxz;
  double ha, pi, c0, c1, cmu, c2;
  double beta, elapsed;
  struct timeval tic, toc;
  
  charge = (double *)calloc(nparts, sizeof(double));
  RPY    = (double *)calloc(3*nparts, sizeof(double));
  
  // ha - raids of the particle
  ha = 0.1;
  pi = 4*atan(1.0); 
  
  cmu = 1.0; 
  c0 = cmu/(6.0*pi*ha);
  c1 = cmu/(8.0*pi);
  c2 = cmu*ha*ha/(12.0*pi);

  test_parser(argc, argv, &beta, &nparts, &s, &accuracy, &distribution);
  test_init(distribution, nparts, accuracy, s, beta, &ploc, &pcharge, &pot, 
   &field, &dxx, &dyy, &dzz, &dxy, &dxz, &dyz);
  
  printf("\n\tPROGRESS\n");
  printf("======================================================\n");
  gettimeofday(&tic, 0);
//*************************************************
// step 1: First three FMM calling
//*************************************************
  // FMM Initialization
  
  adap_fmm_init(accuracy, nparts);
  
  gettimeofday(&toc, 0);
  elapsed = (double) ((toc.tv_usec - tic.tv_usec)/1.0e6 + toc.tv_sec - tic.tv_sec);
  printf("t(FMM_INIT)          | %20.4e\n", elapsed);
 
  int i, m;
  for(i = 0; i < 3; i++){ 
    for(m = 0; m < nparts; m++){
       charge[m] = pcharge[3*m + i];  //charge[i*nparts + m]; //July 9th
	}       
	printf("i = %d", i);
    	// FMM GRAPH
	adap_fmm_graph(nparts, s, beta, ploc, charge);
        printf("graph finished"); 
    	// FMMS compute 
	adap_fmm_compute();
  
    	// FMM POST
	adap_fmm_post(pot, field, dxx, dyy, dzz, dxy, dxz, dyz);
	
	// FMM CLEAN
	adap_fmm_clean();
	
	// Update for RPY
	int i1 = (i+1) % 3; 
	int i2 = (i+2) % 3;
	
	for(m = 0; m < nparts; m++)
	{   	  
		RPY[3*m+i] += c1* pot[m];            // pot[i * nparts + m] +=  c1* ppot[m];
		   
	    RPY[3*m + i]  += (-c1) * ploc[3*m + i] * field[3*m + i];    //pot[i *nparts + m] += (-c1) * ans;
	    RPY[3*m + i1] += (-c1) * ploc[3*m + i] * field[3*m + i1];   //pot[i1*nparts + m] += (-c1) * ans1; 
	    RPY[3*m + i2] += (-c1) * ploc[3*m + i] * field[3*m + i2];    //pot[i2*nparts + m] += (-c1) * ans2; 
	
		if(i==0){
	
			RPY[3*m ]    -= c2 * dxx[m];
			RPY[3*m + 1] -= c2 * dxy[m];
			RPY[3*m + 2] -= c2 * dxz[m];
		
		}else if(i==1){
	
			RPY[3*m ]    -= c2 * dxy[m];
			RPY[3*m + 1] -= c2 * dyy[m];
			RPY[3*m + 2] -= c2 * dyz[m];
	
		}else{
	
			RPY[3*m ]    -= c2 * dxz[m];
			RPY[3*m + 1] -= c2 * dyz[m];
			RPY[3*m + 2] -= c2 * dzz[m];
		}
	}
 }
 //*************************************************
// step 2: The 4th FMM calling
//**************************************************

    for(m = 0; m < nparts; m++) {
	
	    charge[m] = 0;
		pot[m] = 0;
		field[3*m] = 0;
		field[3*m+1] = 0;
		field[3*m+2] = 0;
		dxx[m] = 0;
		dyy[m] = 0;
		dzz[m] = 0.0;
		dxy[m] = 0.0;
		dxz[m] = 0.0;
		dyz[m] = 0.0;
		
	    for(i = 0; i < 3; i++){
		
           charge[m]  += ploc[3*m + i] *pcharge[3*m + i];  //charge[i*nparts + m];
	    }   
	   
	   charge[m] = c1 * charge[m];
	   
	}
    // FMM GRAPH
	adap_fmm_graph(nparts, s, beta, ploc, pcharge);
 
    // FMMS compute 
	adap_fmm_compute();
  
    // FMM POST
	adap_fmm_post(pot, field, dxx, dyy, dzz, dxy, dxz, dyz);
	
	// FMM CLEAN
	adap_fmm_clean(); 
    
    // Update RPY	
	for(m = 0; m < nparts; m++){
    	for(i=0; i<3; i++)
	    {   	  
		   pot[3*m+i] += field[m+i];  		
	    }
    }
	test_clean(ploc, pcharge, charge, pot, RPY,field, dxx, dyy, dzz, dxy, dxz, dyz);

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
  *s = 80;    // s - max number of particles per partition box
  *accuracy = 3;  // accuracy - number of digits of accuracy (3 or 6 supported) 
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
  (*ploc)    = (double *)calloc(3*nparts, sizeof(double));
  (*pcharge)  = (double *)calloc(3*nparts, sizeof(double));
  (*pot)     = (double *)calloc(nparts, sizeof(double));
  
  (*field)   = (double *)calloc(3*nparts, sizeof(double));
  (*dxx) = (double *)calloc(nparts, sizeof(double));
  (*dyy) = (double *)calloc(nparts, sizeof(double));
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

void test_data(int nparts, int distribution, double *ploc, double *pcharge)
{
  if ( distribution == 1 ) {
    // generate data randomly distributed in a unit cube
    int i;
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
      ploc[j]   = ((double) rand()/RAND_MAX - 0.5);
      ploc[j+1] = ((double) rand()/RAND_MAX - 0.5);
      ploc[j+2] = ((double) rand()/RAND_MAX - 0.5);
	  
	  // July 9th 
	  pcharge[j]   = (double) rand()/RAND_MAX;
	  pcharge[j+1] = (double) rand()/RAND_MAX;
	  pcharge[j+2] = (double) rand()/RAND_MAX;
    }
	
	/*for (i=0; i<3; i++){
		for(int j=0; j<nparts; j++)
			charge[i*nparts+j] = (double) rand()/RAND_MAX;
	}
	*/ // July 9th 
	
  } else if ( distribution == 2 ) {
    // generate data randomly distributed over a spherical surface
    int i;
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
      double theta = ((double) rand()/RAND_MAX - 0.5)*M_PI;
      double phi   = ((double) rand()/RAND_MAX - 0.5)*M_PI;
      
	  ploc[j]   = sin(theta)*sin(phi);
      ploc[j+1] = sin(theta)*cos(phi);
      ploc[j+2] = cos(theta);
	  // July 9th 
	  pcharge[j]   = (double) rand()/RAND_MAX;
	  pcharge[j+1] = (double) rand()/RAND_MAX;
	  pcharge[j+2] = (double) rand()/RAND_MAX;
    }
	
	/*for (i=0; i<3; i++){
		for(int j=0; j<nparts; j++)
			charge[i*nparts+j] = (double) rand()/RAND_MAX;
	} */ // July 9th 
  } else if ( distribution == 3 ) {
    // generate data randomly distributed over a spherical surface that is in the first octant
    int i;
    for ( i = 0; i < nparts; i++ ) {
      int j = 3*i;
	  double theta = ((double) rand()/RAND_MAX)*M_PI_2;
      double phi = ((double) rand()/RAND_MAX)*M_PI_2;
      ploc[j]      = sin(theta)*sin(phi);
      ploc[j+1]    = sin(theta)*cos(phi);
      ploc[j+2]    = cos(theta);
	  
	  // July 9th 
	  pcharge[j]   = (double) rand()/RAND_MAX;
	  pcharge[j+1] = (double) rand()/RAND_MAX;
	  pcharge[j+2] = (double) rand()/RAND_MAX;
	  
    }
	
	/*for (i=0; i<3; i++){
		for(int j=0; j<nparts; j++)
			charge[i*nparts+j] = (double) rand()/RAND_MAX;
	} 	*/ // July 9th 
  }
  
  return;
}

void test_clean(double *ploc, double *pcharge, double *charge, double *pot, double *RPY, double *field,
		  double *dxx, double *dyy, double *dzz, double *dxy, double *dxz, double *dyz )
{
  free(ploc);
  free(pcharge);
  free(charge);
  free(pot);
  free(RPY);
  free(field);
  free(dxx);
  free(dyy);
  free(dzz);
  free(dxy);
  free(dxz);
  free(dyz);
}



