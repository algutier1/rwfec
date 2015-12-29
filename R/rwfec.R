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
#' 
#' Example coder to aid the discussion.
#' 
#'                        (1,1,0) t = 6 5 4 3 2 1
#'                  -------> + -->    0 1 1 1 0 1
#'                 |         ^
#'                 |     ____|
#' t = 4 3 2 1     |    |               t =   12 11 10 9 8 7 6 5 4 3 2 1
#'     1 0 1 1  -> i -> x -> x2         y  =   1  0  0 1 0 1 1 1 1 0 1 1
#'          |         |
#'          |         v         t = 6 5 4 3 2 1 
#'           -------> + -------->   1 0 0 1 1 1
#'                 (1,0,1)
#'                
#'     G1 = 1, 1, 0
#'     G2 = 1, 0, 1           
#'      
#'    The figure illustrates a rate 1/2 (k=1, n=2), constraint length K = 3, convolutional
#'    code. 
#'   
#' @param x bits received vector of bits (0's and 1's)
#' @param G matrix of n code words n rows, K length code words
#' @return Returns a complex vector of QPSK symbols. If Ns > 1
#' then the returned signal is shaped with pulse shape, p.
#' @examples
#'   g1=c(1,1,0)
#'   g2=c(1,0,1)
#'  # t= 0 1 2 3 ...
#'   x=c(1,1,0,1)
#'   G=rbind(g1,g2)
#'   y=rwfec(x,G)
#' @export
 convcoder <- function(x,G) {
   y=0
   if (is.matrix(G)) {
     cat(sprintf("hello \n"))
     Gdim=dim(G)
     yy=matrix(rep(0,length(x)+Gdim[2]-1),Gdim[1],length(x)+Gdim[2]-1)
     print(G)
     print(yy)
     for (i in 1:Gdim[1]) {
       print(i)
       print(G[i,])
       yyy= rwconv(x,G[i,]) %% 2
       yy[i,] = yyy
       print(yy[i,])
     }
   } else {
     cat(sprintf("G must be a matrix \n"))
   }
   y = as.vector(yy)
   return(y)
 }

