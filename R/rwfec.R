# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#' Hello, World!
#'
#' Prints 'Hello, world!'
#' @examples
#' hello()
#' @export
hello <- function() {
  print("Hello, world!")
}


#' Convolutional Coder
#' This convcoder documentation
#' 
#' @param x bits received vector of bits (0's and 1's)
#' @param G matrix of n code words n rows, K length code words
#' @return ldsfjsldaf
#' @examples
#' # Example coder encodes the sequence 1101. The figure illustrates 
#' # a rate 1/2 (k=1, n=2), constraint length K = 3, convolutional
#' # code. Convolutional code is described by the codes g1 = (1,1,0) 
#' # and g2 = (1,1,1). Coded output is 1,1,0,1,1,1,1,0,1,0,0,1
#' #
#' #                        (1,1,1) n = 6 5 4 3 2 1
#' #                  -------> + -->    1 1 0 0 0 1
#' #                 |         ^
#' #                 |     ____|
#' # n = 4 3 2 1     |    |               n =   12 11 10 9 8 7 6 5 4 3 2 1
#' #   = 1 0 1 1  -> i -> x -> x2         y  =  1  1  0  1 0 0 1 0 1 0 1 1
#' #                 |         |
#' #                 |         v         n = 6 5 4 3 2 1 
#' #                  -------> + -------->   1 0 0 1 1 1
#' #                        (1,0,1)
#' #     g1 = 1, 1, 0
#' #     g2 = 1, 0, 1           
#' # Example 
#'   g1=c(1,1,1)
#'   g2=c(1,0,1)
#'   x=c(1,1,0,1)  # t= 0 1 2 3 ...
#'   G=rbind(g1,g2)
#'   y=convcoder(y,G)
#' @export
 convcoder <- function(x,G) {
   y=0
   if (is.matrix(G)) {
     Gdim=dim(G)
     yy=matrix(rep(0,length(x)+Gdim[2]-1),Gdim[1],length(x)+Gdim[2]-1)
     for (i in 1:Gdim[1]) {
       yi= rwconv(x,G[i,]) %% 2
       yy[i,] = yi
     }
   } else {
     cat(sprintf("G must be a matrix \n"))
   }
   y = as.vector(yy)
   N = length(y) - Gdim[1]*(Gdim[2]-1)
       # the last  Gdim[1]*(Gdim[2]-1) entries are equivalent to "zero flushing" the coder so remove
   if ( length(y) > N )  y = y[1:N]
   return(y)
 }

