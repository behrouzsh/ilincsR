
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <R.h> // to show errors in R


double calcMean (double *x, int n);
double calcStdev (double *x, double mu, int n);
double calcCov(double *x, double *y, int n, double xmu, double ymu);        

void rCorrelationWrapper2 ( double *X, int *dim, double *mu, double *sd, int *RowRange, int *ColRange, double *corr) {

    int i, j, n = dim[0], p = dim[1];
    int RowStart = RowRange[0], RowEnd = RowRange[1], ColStart = ColRange[0], ColEnd = ColRange[1];
    double xyCov;

    Rprintf("\n p: %d, %d <= row < %d, %d <= col < %d", p, RowStart, RowEnd, ColStart, ColEnd);

    if(RowStart==ColStart && RowEnd==ColEnd){
        for(i=RowStart; i<RowEnd; i++){
            for(j=i; j<ColEnd; j++){
                Rprintf("\n i: %d, j: %d, p: %d", i, j, p);
                xyCov = calcCov(X + i*n, X + j*n, n, mu[i], mu[j]);
                *(corr + j*p + i) = xyCov/(sd[i]*sd[j]);
            }
        }
    } else {
        for(i=RowStart; i<RowEnd; i++){
            for (j=ColStart; j<ColEnd; j++){
                xyCov = calcCov(X + i*n, X + j*n, n, mu[i], mu[j]);
                *(corr + j*p + i) = xyCov/(sd[i]*sd[j]);
            }
        }
    }
}


// function to calculate mean 

double calcMean (double *x, int n){
    double s = 0;
    int i;
    for(i=0; i<n; i++){     
        s = s + *(x+i);
    }
    return(s/n);
}

// function to calculate standard devation

double calcStdev (double *x, double mu, int n){
    double t, sd = 0;
    int i;

    for (i=0; i<n; i++){
        t = *(x + i) - mu;
        sd = sd + t*t;
    }    
    return(sqrt(sd/(n-1)));
}


// function to calculate covariance

double calcCov(double *x, double *y, int n, double xmu, double ymu){
    double s = 0;
    int i;

    for(i=0; i<n; i++){
        s = s + (*(x+i)-xmu)*(*(y+i)-ymu);
    }
    return(s/(n-1));
}

//############################################################

  #include <stdio.h>
  #include <math.h>

  int main() {
        int x[100], y[100], xy[100], xsquare[100], ysquare[100];
        int i, n, xsum, ysum, xysum, xsqr_sum, ysqr_sum;
        float coeff, num, deno;

        xsum = ysum = xysum = xsqr_sum = ysqr_sum = 0;

        /* get the number of entries from the user */
        printf("Enter the value for n:");
        scanf("%d", &n);

        /* get the values for x and y  from the user */
        printf("Enter the value for x and y:\n");
        for (i = 0; i < n; i++) {
                printf("x[%d] and y[%d]: ", i, i);
                scanf("%d%d", &x[i], &y[i]);
        }

        /* find the needed data to manipulate correlation coeff */
        for (i = 0; i < n; i++) {
                xy[i] = x[i] * y[i];
                xsquare[i] = x[i] * x[i];
                ysquare[i] = y[i] * y[i];
                xsum = xsum + x[i];
                ysum = ysum + y[i];
                xysum = xysum + xy[i];
                xsqr_sum = xsqr_sum + xsquare[i];
                ysqr_sum = ysqr_sum + ysquare[i];
        }

        num = 1.0 * ((n * xysum) - (xsum * ysum));
        deno = 1.0 * ((n * xsqr_sum - xsum * xsum)* (n * ysqr_sum - ysum * ysum));

        /* calculate correlation coefficient */
        coeff = num / sqrt(deno);

        /* print the result */
        printf("Correlation Coefficient : %.4f\n", coeff);
        return 0;
  }

  
  
//############################################################
// Program to find correlation coefficient in c++
#include<bits/stdc++.h>
 
using namespace std;
 
// function that returns correlation coefficient.
float correlationCoefficient(int X[], int Y[], int n)
{
 
    int sum_X = 0, sum_Y = 0, sum_XY = 0;
    int squareSum_X = 0, squareSum_Y = 0;
 
    for (int i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X = sum_X + X[i];
 
        // sum of elements of array Y.
        sum_Y = sum_Y + Y[i];
 
        // sum of X[i] * Y[i].
        sum_XY = sum_XY + X[i] * Y[i];
 
        // sum of square of array elements.
        squareSum_X = squareSum_X + X[i] * X[i];
        squareSum_Y = squareSum_Y + Y[i] * Y[i];
    }
 
    // use formula for calculating correlation coefficient.
    float corr = (float)(n * sum_XY - sum_X * sum_Y) 
                  / sqrt((n * squareSum_X - sum_X * sum_X) 
                      * (n * squareSum_Y - sum_Y * sum_Y));
 
    return corr;
}
 
// Driver function
int main()
{
 
    int X[] = {15, 18, 21, 24, 27};
    int Y[] = {25, 25, 27, 31, 32};
 
    //Find the size of array.
    int n = sizeof(X)/sizeof(X[0]);
 
    //Function call to correlationCoefficient.
    cout<<correlationCoefficient(X, Y, n);
 
    return 0;
}
