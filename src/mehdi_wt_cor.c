/* mehdi_wt_cor.c */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/* Simple EXPressions from R
  x is a vector, y is a matrix, w is a matrix
  all have the same number of rows, 
  y and w have the same number of columns,
  the output correlation vector has the 
  variables from y as columns */

SEXP mehdi_wt_cor(SEXP x, SEXP y, SEXP w){
  /* i rows, j columns, k index of respondent
  npair = number of complete observations
  ncx = number of columns of x
  ncy = number of columns of y */

  int i, j, k, n, ncx, ncy, df;
  double sumx, sumy, sumxy, meanx, meany, sdx, sdy, sumx2, sumy2, cor, stat, npair;
  SEXP mat;
//   SEXP stat;
//   SEXP df;

  /*  get dimensions of data */
  n = INTEGER(getAttrib(x, R_DimSymbol))[0];
  ncx = INTEGER(getAttrib(x, R_DimSymbol))[1] + 2;
  ncy = INTEGER(getAttrib(y, R_DimSymbol))[1];

  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(y = coerceVector(y, REALSXP));
  PROTECT(w = coerceVector(w, REALSXP));
  double *xx = REAL(x), *yy = REAL(y), *ww = REAL(w);

  /* allocate space for output matrix */
  PROTECT(mat = allocMatrix(REALSXP, ncx, ncy));
//   PROTECT(stat = allocMatrix(REALSXP, ncx, ncy));
//   PROTECT(df = allocMatrix(INTSXP, ncx, ncy));

  /* for each i and j calculate each correlation and save in matrix mat */
//   for (i = 0; i < ncx; i++){
  for (i = 0; i < 1; i++){
    for (j = 0; j < ncy; j++){
      npair = sumx = sumy = sumxy = sumx2 = sumy2 = df = 0;
      
      /* for loop to calculate sums, the for loop is needed as we have a check for pairwise missing: ISNAN  */
      for (k = 0; k < n; k++){
        if (!ISNAN(xx[k+i*n]) && !ISNAN(yy[k+j*n]) && !ISNAN(ww[k+j*n])){
          npair += ww[k+j*n];
          sumx += ww[k+j*n]*xx[k+i*n];
          sumy += ww[k+j*n]*yy[k+j*n];
          sumxy += ww[k+j*n]*xx[k+i*n]*yy[k+j*n];
        }
      }

      meanx = sumx / npair;
      meany = sumy / npair;

      for (k = 0; k < n; k++){
        if (!ISNAN(xx[k+i*n]) && !ISNAN(yy[k+j*n]) && !ISNAN(ww[k+j*n])){
          sumx2 += ww[k+j*n]*xx[k+i*n]*xx[k+i*n] - 2*ww[k+j*n]*xx[k+i*n]*meanx + ww[k+j*n]*meanx*meanx;
          sumy2 += ww[k+j*n]*yy[k+j*n]*yy[k+j*n] - 2*ww[k+j*n]*yy[k+j*n]*meany + ww[k+j*n]*meany*meany;
	  df += 1;
        }
      }

      sdx = sqrt(sumx2/(npair-1));
      sdy = sqrt(sumy2/(npair-1));

      /*  final correlation */
      cor = (sumxy - npair*meanx*meany)/((npair-1)*sdx*sdy);
      stat = sqrt(df - 2) * cor / sqrt(1 - cor * cor);

      REAL(mat)[i+j*ncx] = cor;
      REAL(mat)[i+1+j*ncx] = stat;
      if((df - 2) > 0) {
	REAL(mat)[i+2+j*ncx] = df - 2;
      } else { REAL(mat)[i+2+j*ncx] = stat; }
    }
  }

  UNPROTECT(4);

  return mat;
}

// stat = sqrt(df) * r / sqrt(1 - r^2)
// df = n - 2L
// p = 2 * min(pt(STATISTIC, df), pt(STATISTIC, df, lower.tail=FALSE)))

