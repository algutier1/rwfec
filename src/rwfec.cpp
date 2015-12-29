#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/


// [[Rcpp::export]]
List rcpp_hello() {
  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}

NumericVector correlateCpp(NumericVector a, NumericVector b) {
  int na = a.size(), nb = b.size();
  int nab = na + nb -1;
  NumericVector xab(nab);
  for (int i =0; i < na; i++)
    for (int j = 0; j < nb; j++)
      xab[i+j] += a[i] * b[j];
  return xab;
}

// Convolution Example
//          inf
//y(n) =    sum[x(k)h(n-k)]
//       k=-inf
// [[Rcpp::export]]
NumericVector rwconv(NumericVector h, NumericVector x) {
  int nh = h.size(), nx = x.size();
  int ny = nh + nx - 1;
  NumericVector y(ny);
  for (int n = 0; n < ny; n++) {
    for (int k = 0; k <= n; k++)  { 
   /*  printf("n = %d, k = %d, n-k = %d, x[%d] = %f, h[%d]=%f, y[%d] = %f \n",n, k, n-k, k, x[k], n-k, h[n-k], n, y[n]);*/
      y[n] += x[k]*h[n-k];
    }
  if (y[n] < 1e-10 ) y[n] = 0;
  }
  return y;
}


/* sinc function http://www.inside-r.org/node/175318
package phonTools */
/* error shown below, no matching function for call to sin,
however, this seems to compile just fine. I believe because
it uses Rcpp sugar, not completely sure, which RStudio does
not seem to know about*/
// [[Rcpp::export]]
NumericVector sinc(NumericVector x) {
  int nx = x.size();
  NumericVector y(nx);
  for (int n = 0; n < nx; n++)
    if (x[n]==0) y[n] = 1;
    else y[n] = sin(x[n])/x[n];
    return y;
}
