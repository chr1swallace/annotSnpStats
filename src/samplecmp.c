#include <stdio.h>
#include <R.h>

/* Chris Wallace <chris.wallace@cimr.cam.ac.uk> */
/* GPL */

void samplecmp(char *x, char *y, int *counts, int *Nx, int *Mx, int *Ny, int *My) {
  int nx = *Nx;
  int mx = *Mx;
  int ny = *Ny;
  int my = *My;
  int i=0, j=0, ii=0, jj=0, k=0;

  if(mx != my)
    error("x and y need equal number of columns");

  int ij=0;
  for(i=0; i<nx; i++) { // index rows of x
    for(j=0; j<ny; j++) { // index rows of y, ij indexes counts
      for(k=0, ii=i, jj=j; k<mx; k++, ii+=nx, jj+=ny) { // index elements of each row in x, y
	int xx=(int) x[ii];
	int yy=(int) y[jj];
	if( xx!=0 && yy!=0 && xx==yy ) {
	  //	  Rprintf("%d %d %d %d : [%d, %d] : %d == %d\n",i,j,k,ij,ii,jj,xx,yy);
	  counts[ij]++;
	}
      }
      // Rprintf("%d %d\n",ij,counts[ij]);
      ij++;
    }
  }
	
  return;

}

void samplediff(char *x, char *y, int *counts, int *maxDiff, int *Nx, int *Mx, int *Ny, int *My) {
  int maxdiff = *maxDiff;
  int nx = *Nx;
  int mx = *Mx;
  int ny = *Ny;
  int my = *My;  
  int i=0, j=0, ii=0, jj=0, k=0;

  if(mx != my)
    error("x and y need equal number of columns");

  int ij=0;
  for(i=0; i<nx; i++) { // index rows of x
    for(j=0; j<ny; j++) { // index rows of y, ij indexes counts
      for(k=0, ii=i, jj=j; k<mx; k++, ii+=nx, jj+=ny) { // index elements of each row in x, y
	int xx=(int) x[ii];
	int yy=(int) y[jj];
	if( xx!=0 && yy!=0 && xx!=yy ) {
	  //	  Rprintf("%d %d %d %d : [%d, %d] : %d == %d\n",i,j,k,ij,ii,jj,xx,yy);
	  counts[ij]++;
	  if(counts[ij] > mx)
	    break;
	}
      }
      // Rprintf("%d %d\n",ij,counts[ij]);
      ij++;
    }
  }
	
  return;

}
