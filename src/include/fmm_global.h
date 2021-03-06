#include <sys/time.h>
#include <stdio.h>
#include "fmm_ds.h"

int PTERMS, NLAMBS, PGSZ, NLEV, NBOXES, *NUMFOUR, *NUMPHYS, 
  *PERM, *CONTENT, *LEAFBOX, NPARTS, PTERMS2, NEXPMAX;

double *FMMLOC, *FMMCHARGE, *FMMPOT, *FMMFIELD, *FMMPOTN, *FMMFIELDN, 
  *FMMDXX,  *FMMDYY,  *FMMDZZ,  *FMMDXY,  *FMMDXZ,  *FMMDYZ, 
  *FMMDXXN, *FMMDYYN, *FMMDZZN, *FMMDXYN, *FMMDXZN, *FMMDYZN, 
  SIZE,  *WHTS, *RLAMS, *RDPLUS, *RDMINUS, *RDSQ3, *RDMSQ3,
  *RLSC, *CARRAY, *ZS;

dcomplex *MPOLE, *LOCAL, *LEXPU, *LEXPD, *LEXPN, *LEXPS, *LEXPE, 
  *LEXPW, *XS, *YS, *FEXPE, *FEXPO, *FEXPBACK;

fmmbox *BOXES;

fmmlist *LIST;



