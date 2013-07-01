#include <stdio.h>
#include <R.h>

/* Chris Wallace <chris.wallace@cimr.cam.ac.uk> */
/* GPL */

void samplecmp(char *x, char *y, int *nhet, int *Nx, int *Mx, int *Ny, int *My) {
  int nx = *Nx;
  int mx = *Mx;
  int ny = *Ny;
  int my = *My;
  int i=0, j=0, ii=0, jj=0, k=0;

  if(mx != my)
    error("x and y need equal number of columns");

  nhet[1]=4;

  int ij=0;
  for(i=0; i<nx; i++) { // index rows of x
    for(j=0; j<ny; j++, ij++) { // index rows of y, ij indexes nhet
      for(k=0, ii=i, jj=j; k<mx; k++, ii+=nx, jj+=ny) { // index elements of each row in x, y
   	if( (int) x[ii]==2 && (int) y[jj]==2 ) {
	  nhet[ij]++;
	}
      }
    }
  }
	
  return;

}
