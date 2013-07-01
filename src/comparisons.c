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

SEXP countdiffs(SEXP Rx, SEXP Ry, SEXP maxDiff) {

  int nprotect=0;
  SEXP Rdim, Rcounts;

  PROTECT(Rdim = getAttrib(Rx, R_DimSymbol));
  nprotect++;
  int nx = INTEGER(Rdim)[0];
  int mx = INTEGER(Rdim)[1];

  PROTECT(Rdim = getAttrib(Ry, R_DimSymbol));
  nprotect++;
  int ny = INTEGER(Rdim)[0];
  int my = INTEGER(Rdim)[1];

  if(mx != my)
    error("x and y need equal number of columns");

  // maximum number of mismatches allowed
  int maxdiff = INTEGER(maxDiff)[0];
   
   // pointers to x, y
  unsigned char *x = RAW(Rx);
  unsigned char *y = RAW(Ry);

  // return indices of samples with < tol mismatches, number of mismatches
  int M = nx;
  if(M < ny)
    M = ny;
  int counts[M*8];
  
  int i=0, j=0, ii=0, jj=0, k=0;

  int ij=0;
  for(i=0; i<nx; i++) { // index rows of x
    for(j=0; j<ny; j++) { // index rows of y, ij indexes counts
      int nonzero = 0, mismatch=0;
      
      for(k=0, ii=i, jj=j; k<mx; k++, ii+=nx, jj+=ny) { // index elements of each row in x, y
	int xx=(int) x[ii];
	int yy=(int) y[jj];
	if( xx!=0 && yy!=0) {
	  nonzero++;
	  if(xx!=yy ) {
	    mismatch++;
	    if(mismatch == maxdiff)
	      break;
	  }
	}
      }

      if(mismatch == maxdiff) // different samples
	continue;

      // low mismatch - store
      counts[ij] = i;
      counts[ij + M] = j;
      counts[ij + 2*M] = mismatch;
      counts[ij + 3*M] = nonzero;
      ij++;
    }
  }

  // trim Rcount
  PROTECT(Rcounts = allocMatrix(INTSXP, ij, 4));
  nprotect++;
  int *pRcounts = INTEGER(Rcounts);
  for(i=0; i<ij; i++)
    for(j=0; j<4; j++)
      pRcounts[i + j*ij] = counts[i + j*M];
	
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
