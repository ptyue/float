/* write out data that can be read by the visualisation package */

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

void wreslt_(int *it, double *t, int *nbrigid, double *xpos, double *xcv,
             double *ycv, double *strsx, double *strsy,
             int *nnode, int *nvert, int *nelem, int *inod,
             double *x, double *y,
             int *nbd, int *ibdnod, int *nic, int *ic,
             int *nbound, int *nside, int *bdseg, double *pwall,
             int *ncpvelo, int *ncpelas,
             double *u, double *p, double *strm, double *vort,
             double *selas, char *outfile)

{
  register int   i, j, k, n;
  int    nnvar;
  float *rwork;
  FILE  *fileout;

  /* allocate work space */
  nnvar = 2*(*nvert);
  if ( nnvar < 4*(*nbd) ) nnvar = 4*(*nbd);
  if ( nnvar < 3*(*nbrigid) ) nnvar = 3*(*nbrigid);
  rwork = (float *)calloc(nnvar,sizeof(float));
  if (rwork == NULL)
  { printf("from WRESLT: rwork == NULL"); exit(0);
  }
/*
   printf("%6d%6d%6d%6d%6d\n",*it,*nvert,*nnode,*nbd,*nelem);
   printf("%6d%6d%6d\n",*ncpelas,*nbrigid);
*/

  /* open output file */
/*  printf("%s\n", outfile); exit(0); */
  fileout = fopen(outfile, "w");

  /* write data to file */
  fwrite(it, sizeof(int),  1, fileout);

  rwork[0] = (float)*t;
  fwrite(&rwork[0],sizeof(float), 1, fileout);

  fwrite(nvert, sizeof(int),  1, fileout);
  fwrite(nelem, sizeof(int),  1, fileout);

  for (i = 0; i < *nvert; i++) 
  { j = i*2;
    rwork[j]   =  (float)x[i];
    rwork[j+1] =  (float)y[i];
  }
    fwrite(&rwork[0],sizeof(float),2*(*nvert), fileout);

  for (i = 0; i < *nelem; i++) 
  { j = i*6;
    fwrite(&inod[j],sizeof(int), 3, fileout);
  }

  nnvar = 5 + *ncpelas;
  fwrite(&nnvar, sizeof(int),  1,     fileout);

  for (i = 0; i < *nvert; i++)
  { j = i*2;
    rwork[i]        =  (float)u[j];
    rwork[*nvert+i] =  (float)u[j+1];
  }
  fwrite(&rwork[0],     sizeof(float),*nvert, fileout);
  fwrite(&rwork[*nvert],sizeof(float),*nvert, fileout);

  for (i = 0; i < *nvert; i++) rwork[i] =  (float)p[i];
  fwrite(&rwork[0],sizeof(float),*nvert, fileout);

  for (i = 0; i < *nvert; i++) rwork[i] =  (float)strm[i];
  fwrite(&rwork[0],  sizeof(float),*nvert, fileout);

  for (i = 0; i < *nvert; i++) rwork[i] =  (float)vort[i];
  fwrite(&rwork[0],  sizeof(float),*nvert, fileout);

  if ( *ncpelas>0 )
    for (j=0; j<*ncpelas; j++)
    { for (i = 0, k=j; i < *nvert; i++, k+=*ncpelas)
        rwork[i] =  (float)selas[k];
      fwrite(&rwork[0],sizeof(float),*nvert, fileout);
    }

  fwrite(nnode, sizeof(int),  1,     fileout);
  fwrite(nbd,   sizeof(int),  1,     fileout);
  fwrite(nic,sizeof(int),  1,     fileout);
  fwrite(&ibdnod[0], sizeof(int), *nbd,      fileout);
  fwrite(&ic[0], sizeof(int), *nic+1, fileout);
  fwrite(nbrigid, sizeof(int), 1, fileout);

  for (i = 0 ; i < 3*(*nbrigid); i++) rwork[i] = (float)xpos[i];
  fwrite(&rwork[0],sizeof(float),3*(*nbrigid), fileout);

  fwrite(nbound, sizeof(int), 1,  fileout);
  fwrite(&nside[0],sizeof(int), *nbound, fileout);

  fwrite(bdseg, sizeof(int), 1,  fileout);
  for (i = 0 ; i < 5*(*bdseg); i++) rwork[i] = (float)pwall[i];
  fwrite(&rwork[0],sizeof(float),5*(*bdseg), fileout);

  fwrite(nbd, sizeof(int), 1,  fileout);
  for (n=0 ; n < *nbd; n++)
    { j=n*4; 
      rwork[j]   = (float)xcv[n];
      rwork[j+1] = (float)ycv[n];
      rwork[j+2] = (float)strsx[n];
      rwork[j+3] = (float)strsy[n];
    }
  fwrite(&rwork[0],sizeof(float), 4*(*nbd), fileout);

  fclose(fileout);

  /* free work space */
  free(rwork);

}
