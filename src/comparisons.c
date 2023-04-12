#include <stdio.h>
#include "Rinternals.h"

/* Chris Wallace <chris.wallace@cimr.cam.ac.uk> */
/* GPL */

/* SEXP countmatches(SEXP Rx, SEXP Ry) { */

/*   int nprotect=0; */
/*   SEXP Rdim, Rcounts, Rxx, Ryy; */

/*   PROTECT(Rdim = getAttrib(Rx, R_DimSymbol)); */
/*   nprotect++; */
/*   int nx = INTEGER(Rdim)[0]; */
/*   int mx = INTEGER(Rdim)[1]; */

/*   PROTECT(Rdim = getAttrib(Ry, R_DimSymbol)); */
/*   nprotect++; */
/*   int ny = INTEGER(Rdim)[0]; */
/*   int my = INTEGER(Rdim)[1]; */

/*   if(mx != my) */
/*     error("x and y need equal number of columns"); */

/*   // pointers to x, y */
/*   PROTECT(Rxx = coerceVector(Rx, INTEGER)); */
/*   PROTECT(Ryy = coerceVector(Ry, INTEGER)); */
/*   nprotect+=2; */
/*   int *xx = CHAR(Rxx); */
/*   int *yy = CHAR(Ryy); */

/*   // return matrix */
/*   PROTECT(Rcounts = allocMatrix(INTSXP, nx, ny)); */
/*   nprotect++; */
/*   int *counts = INTEGER(Rcounts); */
  
/*   int i=0, j=0, ii=0, jj=0, k=0; */

/*   int ij=0; */
/*   for(i=0; i<nx; i++) { // index rows of x */
/*     for(j=0; j<ny; j++) { // index rows of y, ij indexes counts */
/*       for(k=0, ii=i, jj=j; k<mx; k++, ii+=nx, jj+=ny) { // index elements of each row in x, y */
/* 	/\* int xx=(int) x[ii]; *\/ */
/* 	/\* int yy=(int) y[jj]; *\/ */
/* 	if( xx[ii]!=0 && yy[jj]!=0 && xx[ii]==yy[jj] ) { */
/* 	  //	  Rprintf("%d %d %d %d : [%d, %d] : %d == %d\n",i,j,k,ij,ii,jj,xx,yy); */
/* 	  counts[ij]++; */
/* 	} */
/*       } */
/*       // Rprintf("%d %d\n",ij,counts[ij]); */
/*       ij++; */
/*     } */
/*   } */
	
/*   UNPROTECT(nprotect); */
/*   return(Rcounts); */

/* } */

SEXP countdiffs(SEXP Rx, SEXP Ry, SEXP maxDiff, SEXP Rtype, SEXP Rquick, SEXP pBar) {

  int nprotect=0;
  SEXP Rdimx, Rdimy, Rcounts;

  SEXP utilsPackage, percentComplete;
  PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
  PROTECT(percentComplete = allocVector(INTSXP, 1));
  nprotect+=2;
  int *rPercentComplete = INTEGER(percentComplete);
 
  PROTECT(Rdimx = getAttrib(Rx, R_DimSymbol));
  nprotect++;
  int nx = INTEGER(Rdimx)[0];
  int mx = INTEGER(Rdimx)[1];

  PROTECT(Rdimy = getAttrib(Ry, R_DimSymbol));
  nprotect++;
  int ny = INTEGER(Rdimy)[0];
  int my = INTEGER(Rdimy)[1];

  if(mx != my)
    error("x and y need equal number of columns");

  // maximum number of mismatches allowed
  int maxdiff = INTEGER(maxDiff)[0];
   
  // type of things to count 
  int type = INTEGER(Rtype)[0];

  // be quick by assuming <=1 match per sample
  int quick = INTEGER(Rquick)[0];

  // pointers to x, y
  unsigned char *x = RAW(Rx);
  unsigned char *y = RAW(Ry);

  // return indices of samples with < tol mismatches, number of mismatches
  int MM = nx;
  if(MM > ny)
    MM = ny;
  if(quick == 0)
    MM=MM*4; // worst case: at maximum, each sample in x or y may have four matches in y or x
  int xindex[MM], yindex[MM], mismatch[MM], total[MM];
  
  int i=0, j=0, ii=0, jj=0, k=0;

  // record matches
  int xflag[nx], yflag[ny];
  for(i=0; i<nx; i++)
    xflag[i]=0;
  for(i=0; i<ny; i++)
    yflag[i]=0;

  int ij=0;
  for(i=0; i<nx; i++) { // index rows of x
    if(xflag[i]==1) // already matched
      continue;
    for(j=0; j<ny; j++) { // index rows of y, ij indexes counts
      if(yflag[j]==1)
	continue;
      int nonzero = 0, different=0;
      
      for(k=0, ii=i, jj=j; k<mx; k++, ii+=nx, jj+=ny) { // index elements of each row in x, y
	int xx=(int) x[ii];
	int yy=(int) y[jj];
	if( xx!=0 && yy!=0) {
	  nonzero++;
	  if((type==0 && xx!=yy) || (type==1 && xx==2 && yy!=2) || (type==1 && xx!=2 && yy==2)) {
	    different++;
	    if(different == maxdiff)
	      break;
	  }
	}
      }

      if(different == maxdiff) // different samples
	continue;

      // low mismatch - store
      /* fprintf(stderr, "i:%i  j:%i, ij:%i\n", i, j, ij); */
      xindex[ij] = i+1; // switch to 1-based
      yindex[ij] = j+1; // switch to 1-based
      mismatch[ij] = different;
      total[ij] = nonzero;
      if(quick==1) {
	xflag[i] = 1;
	yflag[j] = 1;
      }
      ij++;
    }
   *rPercentComplete = i; //this value increments
   eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
  }

  // trim Rcount
  PROTECT(Rcounts = allocMatrix(INTSXP, ij, 4));
  nprotect++;
  int *pRcounts = INTEGER(Rcounts);
  for(i=0; i<ij; i++) {
    pRcounts[i] = xindex[i];
    pRcounts[i+ij] = yindex[i];
    pRcounts[i+2*ij] = mismatch[i];
    pRcounts[i+3*ij] = total[i];
  }
	
  UNPROTECT(nprotect);
  return(Rcounts);

}


/* SEXP countmatches(SEXP x, SEXP y, SEXP maxDiff) { */

/*   int nprotect=0; */

/*   PROTECT(SEXP Rdim = getAttrib(x, R_DimSymbol)); */
/*   nprotect++; */
/*   int nx = INTEGER(Rdim)[0]; */
/*   int mx = INTEGER(Rdim)[1]; */

/*   PROTECT(SEXP Rdim = getAttrib(y, R_DimSymbol)); */
/*   nprotect++; */
/*   int ny = INTEGER(Rdim)[0]; */
/*   int my = INTEGER(Rdim)[1]; */

/*   if(mx != my) */
/*     error("x and y need equal number of columns"); */

/*   int maxdiff = INTEGER(maxDiff)[0]; */

/*   // return matrix */
/*   PROTECT(Rcounts = allocMatrix(INTSXP, nx * ny)); */
/*   nprotect++; */
/*   counts = INTEGER(Rcounts); */
  
/*   int i=0, j=0, ii=0, jj=0, k=0; */

/*   int ij=0; */
/*   for(i=0; i<nx; i++) { // index rows of x */
/*     for(j=0; j<ny; j++) { // index rows of y, ij indexes counts */
/*       for(k=0, ii=i, jj=j; k<mx; k++, ii+=nx, jj+=ny) { // index elements of each row in x, y */
/* 	int xx=(int) x[ii]; */
/* 	int yy=(int) y[jj]; */
/* 	if( xx!=0 && yy!=0 && xx!=yy ) { */
/* 	  //	  Rprintf("%d %d %d %d : [%d, %d] : %d == %d\n",i,j,k,ij,ii,jj,xx,yy); */
/* 	  counts[ij]++; */
/* 	  if(counts[ij]==maxdiff) */
/* 	    break; */
/* 	} */
/*       } */
/*       // Rprintf("%d %d\n",ij,counts[ij]); */
/*       ij++; */
/*     } */
/*   } */
	
/*   UNPROTECT(nprotect); */
/*   return(Rcounts); */
    
/* } */
