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

// Convolution
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

// Viterbi decoder
//          inf
//y(n) =    sum[x(k)h(n-k)]
//       k=-inf
// [[Rcpp::export]]
NumericVector rwviterbi(NumericVector x, NumericMatrix G, NumericVector tb,  NumericVector v) {
 
  int  Nx = x.size();   /* input vector size */
  int  Ngr = G.nrow();  /* Generator matrix rows (outputs per input bit) */
  int  Ngc = G.ncol();  /* Generator matrix cols (K-1) */
  int  Nv = v.size();   /* Verbose flag */
  int  Ntailbyte = tb.size(); /* tail byte flag */
  int K=Ngc; /* constraint length */
  int Ntb=4*(K-1); /* N tail byte length */
  int Nd = 5; /* trellis decoding depth */
  int Nstates = (int)pow(2,K-1);
  NumericVector trout((int)(Nx/Ngr)); /*  trellis output */
  int oidx=0, Nout=Nx/Ngr;  /* output idx and output cardinality */
   
   
   /* outline of rwviterbi
    * 
    *   Check inputs
    *     iterate over trellis 
    *        buttrefly branch computations (iterate over states)
    *           state1,state2 ---> nextstate1
    *           state1,state2 ---> nextstate2
    *        add tailbyting bits
    *        traceback
    *      final traceback
    *      print verbose information
    */

   /* basic check, input parameters are of 
       correct cardinality and proceed */
   if (Ngc > 1 && Ngr > 1 && Ngr <= 16 && Nx > 1 && Nv == 1 && Ntailbyte==1) {
       /* ok, proceed */
    

    int Ny= Nx+Ntb*tb(0)*Ngr;
    int Nouttmp = Ny/Ngr;
    
    NumericVector y(Ny);  /* copy input x to y and operate on y */
    NumericVector trouttmp(Nouttmp); /* trellis out tmp */
    
    
    /* copy input x to y and add extra space for tail byteing */
    if (v(0) != 0)  printf("Nout = %d, Nouttmp = %d\n", Nout, Nouttmp);
    for (int i = 0; i < Nx ; i++) y(i)=x(i);
    if (tb(0)!=1 ) { if (v(0) != 0 ) printf("tb flag (tail byte) = %d\n", (int)tb(0));
      } else { 
      if (v(0)!= 0) printf("tail byte = 1\n");
      for (int i = Nx; i < Ny; i+=Ngr)
          for (int k=0; k< Ngr; k++) {
              y(i+k)=0;
         }
    }



     /*trellis variables */
     int state1, state2, nextstate1, nextstate2; /* butterfly variables */
     int trellispm   [Nd][Nstates] ;  /* trellis path metric by [depth] and [state] */
     int trellisin   [Nd][Nstates] ;  /* trellis in to each state by [depth] and [state] */  
     int trellisprestate [Nd][Nstates]; /* trellis prev state by [depth] and [state] */
     int trlidx = 0, trlidxp1 = 1; /* trellis indexes */
     int ix = 0;
     int bmtmp1 =0, bmtmp2=0, bm=0;
     
     /*trace back parameters: trbidx, trbstate, trbprevstate */
     int trbidx, trbstate=0, trbprevstate=0;
     
   
   /**********************/
   /*build state map     */
   /*  - next states       */
   /*  - and output tables */
   /**********************/
  
   unsigned int b,a;
   a=0x1000;
   b=0x0000;

    int state_out [Nstates][2][Ngr];
    int state_next [Nstates][2];
    for (int i=0; i < Nstates; i++) {
      for (int j=0; j < 2; j++) {   /* input 0 j = 0, input 1 j= 1 */
           state_next[i][j] = (2*i + j) % Nstates;
        /*   printf(" state = %d, input = %d, state_next = %d, state_prev[%d][%d] = %d\n", i, j,state_next[i][j],  (2*i+j) % Nstates,j, i ); */
          for (int k =0 ; k < Ngr ; k++) {    /* code G[k] = 1 ... Ngr */
            state_out[i][j][k] = j*G(k,0);
          /*   printf("%d * %f %d\n",j,G(k,0), state_out[i][j][k] );  */
            for (int l = 1 ; l < K ; l++)  {  /* code word sum( statebit[l]*G[k][l]) */
              state_out[i][j][k] += (((int)(i/pow(2,l-1))) %2)* G(k,l);
          /*     printf("%d * %f %d \n",(((int)(i/pow(2,l-1))) %2),G(k,l),  state_out[i][j][k]); */
            } 
            state_out[i][j][k] = state_out[i][j][k] % 2;
          }
           
       /*   printf("state = %d, input = %d, state_out_1 = %d, state_out_2 = %d, state_next = %d\n",
                 i,j,state_out[i][j][0], state_out[i][j][1],state_next[i][j] ); */
      
      }
    } /* done building state map */
   
      
    /*initialize trellis parameters
     *   path metrics
     *   previous state
     *   trellis input
     */
  
    for (int i=0; i < Nd; i++) {
      for (int j=0; j < Nstates; j++) {
        trellispm[i][j]=0;
        trellisprestate[i][j]=0;
        trellisin[i][j]=0;
      }
    }
    
    /* initialize output to 0 
    * trout and trouttmp
    */
    for (int i = 0; i < Nout ; i ++) {
      trout[i]=0;
    }
      /* Nouttmp will be longer than Nout if doing tail byting */
    for (int i = 0; i < Nouttmp ; i++) {
      trouttmp[i]=0;
    }
    
    

    /**** Trellis/CirculaR Buffer ****/
    /*
       At each step of the trellis, i, 
          include prevstate, nextate  "butteryfly" computation
          trellis indexes ...
             trlidx   ...  trellis index
             trlidxp1  ... trellis index plus 1
             ix += Ngr   ... increment by Ngr (output bits per input bit) 
         see pic below
         prevstate and nextstate metrics including 
           pathmetric, 
           prevstate  to nextstate corresponding to min path metrick 
           input associated with the prev state to nextstate min path metric
     
        once reach end of input stream, (if tb==1) add tail byte inputs and iterate 
          through trellis until done with tail bytes inputs
    
    example                                 |<-- traceback when i >= Nd-1
    Nd = 5 (trellis depth)                  v
    i =          0       1      2     3     4    5    6   ...
       i = 0  prev     next      
              state    state       each step through trellis looks at prev --> next state transition
          trlidx=0  tridxp1=1
     
                 i = 1 prev    next
                       state   state
                   trlidx=1  tridxp1=2
     
                        |<--first input stored here as transition from initial state   
                        v
  trlidx/p1  =  0/1    1/2    2/3    3/4    4/5    5/0
               <--- trellis circular buffer----------->      
                s1  x  s1  x  s1  x  s1  x  s1  x   s1  x
                s2  x  s2  x  s2  x  s2  x  s2  x   s2  x
                s3  x  s3  x  s3  x  s3  x  s3  x   s3  x
                s4  x  s4  x  s4  x  s4  x  s4  x   s4  x
    */

      /*******************************************************************/
      /* step through trellis ... iterate over inputs and manage trellis */ 
      /*                          butterfly computations                 */
      /*   one iteration for each of Nx/Ngr code bit inputs              */
      /*******************************************************************/
      
      for (int i = 0 ; i < Nouttmp; i++ ){
          
              /**** print if verbose (v0) == 1 is on  ****/
        if (v(0)!=0) {
        printf("\n");
          printf("   trlidx = %d                  trlidxp1 = %d\n", trlidx, trlidxp1);
        }
        
              /**** print if verbose (v0) == 1 is on  ****/
              /* print received sequence  */
         if(v(0)!=0) {
          for (int k =0; k < Ngr; k++) {
            printf("x[%d] = %d, ",ix+k,(int)y[ix+k]);
          }
          printf("\n");
         }
          
          /* butteryfly computations as interate through states */
          
          /*******************************************/
          /* branch metrics...butteryfly comutations */
          /* iterate over inputs, manage trellis     */ 
          /*   iterate from state1 to Nstates/2      */
          /*     input = 1 or 0                      */
          /*     state1 = j                          */
          /*     state2 = Nstates/2 + j              */
          /*     nextstate = 2j                      */
          /*     nextstate = 2j + 1                  */
          /*                                         */  
          /*          0/xx                           */
          /* state1 ------> nextstate1               */
          /*         \  / 0/xx                       */
          /*          /\                             */
          /*        /   \ 1/xx                       */
          /* state2 ----- > nextstate2               */
          /*           1/xx                          */
          /* note:                                   */
          /*   nextstates have 2 inputs              */
          /*   compute bm for each                   */
          /*   select the path with least pm         */
          /*                                         */
          /*******************************************/
                   /**** print if verbose (v0) == 1 is on  ****/
          if (v(0)!=0) {
             printf("state/in/out   -->  input  pmetric, bmetric, prestate, nextstate\n");
          }
          
              /* iterate over states ... butteryfly comps */
          for (int j = 0; j < Nstates/2 ; j++) {  
            
            state1 = j;
            state2 = Nstates/2+j;
            nextstate1 = 2*j;
            nextstate2 = 2*j+1;
            
            /* print if verbose on */
            /* state1 -> nextate1 */
            if (v(0)!=0) {
              printf(" %d/0/",state1);
              for (int k = 0 ; k < Ngr; k++) {
                printf("%d",state_out[state1][0][k]);
              }
            }
            /* print if verbose */
            /* state2 -> nextstate1 */
            if (v(0)!=0) {
              printf(" %d/0/",state2);
              for (int k = 0 ; k < Ngr; k++) {
                printf("%d",state_out[state2][0][k]);
              }
            }
           
           
            /**** nextstate 1 ****/
            /* branch metrics, note input 0 */
            bmtmp1=0;
            bmtmp2=0;
            for (int k=0; k < Ngr; k++) {
              if (y[ix+k] != state_out[state1][0][k]) bmtmp1+=1;
              if (y[ix+k] != state_out[state2][0][k]) bmtmp2+=1;
            } 
            
            if ( bmtmp1 < bmtmp2 ) {
              trellisprestate[trlidxp1][nextstate1]=state1;
              trellispm[trlidxp1][nextstate1]=trellispm[trlidx][state1]+bmtmp1;
              trellisin[trlidxp1][nextstate1]=0;
            } else {
              trellisprestate[trlidxp1][nextstate1]=state2;
              trellispm[trlidxp1][nextstate1]=trellispm[trlidx][state2]+bmtmp2;
              trellisin[trlidxp1][nextstate1]=0;
            }
             
            bm=fmin(bmtmp1,bmtmp2);
            
            /* print if verbose on */
            if (v(0)!=0) {
               printf("     --> %d      %d     %d        %d          %d\n",trellisin[trlidxp1][nextstate1],trellispm[trlidxp1][nextstate1], bm, trellisprestate[trlidxp1][nextstate1], nextstate1); 
            }
            
            /* print if verbose on */
            /*branch metrics state1 -> nextstate2 */
            if (v(0)!=0) {
              printf(" %d/1/",state1);
              for (int k = 0 ; k < Ngr; k++) {
                printf("%d",state_out[state1][1][k]);
              }
            }
            
            /* print if verbose on */
            /*branch metrics state2 -> nextstate2 */
            if (v(0)!=0) {
               printf(" %d/1/",state2);
               for (int k = 0 ; k < Ngr; k++) {
                 printf("%d",state_out[state2][1][k]);
              }
            }
            
            
            /**** nextstate 2 *****/
            /* branch metrics, note input 1 */
            bmtmp1=0;
            bmtmp2=0;
            for (int k=0; k < Ngr; k++) {
              if (y[ix+k] != state_out[state1][1][k]) bmtmp1+=1;
              if (y[ix+k] != state_out[state2][1][k]) bmtmp2+=1;
            }   
            
            if ( bmtmp1 < bmtmp2 ) {
              trellisprestate[trlidxp1][nextstate2]=state1;
              trellispm[trlidxp1][nextstate2]=trellispm[trlidx][state1]+bmtmp1;
              trellisin[trlidxp1][nextstate2]=1;
            } else {
              trellisprestate[trlidxp1][nextstate2]=state2;
              trellispm[trlidxp1][nextstate2]=trellispm[trlidx][state2]+bmtmp2;
              trellisin[trlidxp1][nextstate2]=1;
            }
            
            bm = fmin(bmtmp1,bmtmp2);  
            if (v(0)!=0 ) {
               printf("     --> %d      %d     %d        %d          %d\n",trellisin[trlidxp1][nextstate2],trellispm[trlidxp1][nextstate2],bm,trellisprestate[trlidxp1][nextstate2], nextstate2); 
            }
          } /*end iterate over states ... butteryfly computations */

          
          
          
          
          /**** trace back through trellis ... get trellis output ****/
          /* i.e., estimate intput bit 
           trace back if Nd is reached */

          
          if (i >= Nd - 1 ) { /*start traceback at i = Nd-1, update trbidx (pointer) and trbstate */
            for (int itrb = trlidxp1;  itrb > trlidxp1-Nd ; itrb-- ) { /* trace through trellis, at end of each iteration trbstate, trbidx are correct */
              trbidx=itrb;
              if ( trbidx < 0 ) trbidx=Nd-(-itrb); /* trbidx will become negative as wrap around trellis */
              trbstate=trbprevstate;  /*shift trbstate to trbprevstate */     
              if (itrb == trlidxp1 ) { /* first time find state w/ min path metric */
                for (int itrb2 = 0; itrb2 < Nstates; itrb2++) { 
                  if (trellispm[trbidx][itrb2] < trellispm[trbidx][trbstate])  { 
                     /* printf("  itrb2=%d, trellispm[trbidx][itrb2] = %d, trellispm[trbidx][trbstate] = %d\n",itrb2, trellispm[trbidx][itrb2],trellispm[trbidx][trbstate]); */
                     trbstate = itrb2;
                     /* printf("  new trbstate = %d\n", trbstate); */
                  }
                }
              }
              trbprevstate=trellisprestate[trbidx][trbstate]; /* get the prevstate corresponding to trbstate */
              /* printf("  trbidx = %d, trbstate = %d, trbprevstate = %d\n",trbidx,trbstate,trbprevstate); */
            }  /* trbstate and trbidx (at Nd) are accurate */
                                            /* trellis out is the input at trellis trbidx and trbstate 
                                                first output will correspond to trellis index 1
                                                recall first time through trellis first input is stored
                                                at trellis index of 1 and trellis index of 0 is the prevstate (initial state)*/
            trouttmp[oidx]=trellisin[trbidx][trbstate]; 
            if (v(0)!=0) {
               printf(" i = %d, trlidxp1 = %d, itrb = %d, trbidx = %d, oidx = %d, trbstate = %d, trouttmp[%d] = %d\n",i, trlidxp1, trbidx, trbidx, oidx, trbstate, oidx, (int)trouttmp[oidx]); 
            }
            oidx++;
          }   /* end traceback */
          
         ix+= Ngr;  /* increment trellis input index by Ngr */
         
         
         /*** add tail byte bits ***/
         /* add tail byte y values
          *   note, if tb != 1 then ix is never greather than Nx, no tail byting
          *   if tb = = 1 then once input reaches ix = Nx start add tail byting inputs
          *   effectively add zero's to the input stream to clear 
          *   the coder state. Add zero inputs from the minimum path metric state.
          *   
          *   
          */                                   
         if (ix >= Nx && ix <= Ny-Ngr && tb(0) == 1)    {
            if (v(0)!=0) printf("tail byte bits\n");
            /* go back to ix-Ngr   .... find state with lowest path metric
             y(ix,ix+1... ix+Ngr)=  input=0:state/output
             print ix ... state .. y */
            int trbstate=0;
            
            for (int j = 0; j < Nstates; j ++) { 
              if (trellispm[trlidxp1][j] < trellispm[trlidxp1][trbstate])  { 
                /* printf("  j = %d, trellispm[trlidxp1][j] = %d, trellispm[trlidxp1][trbstate] = %d\n",j, trellispm[trlidxo1][j],trellispm[trlidxp1][trbstate]); */
                trbstate = j;
                if (v(0)!=0) printf("  new trbstate = %d\n", trbstate);
              }
            }
           if (v(0)!=0) printf("  ix = %d, min pm state = %d",ix, trbstate);
            for (int k = 0; k < Ngr ; k++) {
              y(ix+k) = state_out[trbstate][0][k];
              if (v(0)!=0) printf("  y[%d] = %d",ix+k, (int)y(ix+k));
            }
            if(v(0)!=0) printf("\n");
         }                               
                                            
         trlidx = (trlidx+1) % Nd;
         trlidxp1 = (trlidxp1+1) % Nd;
         
         
         
      } /* end iterate over trellis */


      /***** final trace back ***/
      /* final traceback ... get output bits less than Nd deep */  
         /*start traceback at i = Nd-1, update trbidx (pointer) and trbstate */   
      trbidx=trlidxp1;
      trbprevstate=0;
      trbstate=0;
      for  (int i = Nouttmp-1; i >= oidx; i-- ) { 
           if ( trbidx < 0 ) trbidx=Nd-1; /* trbidx will become negative as wrap around trellis */
           trbstate=trbprevstate;  /*shift trbstate to trbprevstate */  
           if (v(0)!=0) { printf(" i = %d, trbidx =%d\n",i,trbidx); }
           if (i == Nouttmp-1 ) { /* first time find state w/ min path metric */
             for (int itrb2 = 0; itrb2 < Nstates; itrb2++) { 
               if (trellispm[trbidx][itrb2] < trellispm[trbidx][trbstate])  { 
                 /* printf("  itrb2=%d, trellispm[trbidx][itrb2] = %d, trellispm[trbidx][trbstate] = %d\n",itrb2, trellispm[trbidx][itrb2],trellispm[trbidx][trbstate]); */
                 trbstate = itrb2;
                 if (v(0)!=0) { printf("  min metric trbstate = %d\n", trbstate); }
               }
             }
           }
          trbprevstate=trellisprestate[trbidx][trbstate];   /* get the prevstate corresponding to trbstate */
          /* printf("  trbidx = %d, trbstate = %d, trbprevstate = %d\n",trbidx,trbstate,trbprevstate); */

        /* trellis out is the input at trellis trbidx and trbstate 
           first output will correspond to trellis index 1
           recall first time through trellis first input is stored
           at trellis index of 1 and trellis index of 0 is the prevstate (initial state)*/
          if (v(0)!=0) {
             printf(" i = %d, trbidx = %d, trbstate = %d, trouttmp = %d\n", i, trbidx, trbstate, (int)trouttmp[i]);
            }
          trouttmp[i]=trellisin[trbidx][trbstate];
        if (v(0)!=0) { printf(" i = %d, trlidxp1 = %d, itrb = %d, trbidx = %d, oidx = %d, trbstate = %d, trout[%d] = %d\n",i, trlidxp1, trbidx, trbidx, oidx, trbstate, oidx, (int)trout[oidx]); }
         trbidx--;
        }   /* end final traceback */
      
      /***** copy trouttmp to trout ****/
      for (int i =0 ; i < Nout ; i++) trout(i)=trouttmp(i);
      
    /*********************************************/
    /* Verbose (v !=0), print path metrics       */    
    /*  -  state map and other parameters        */
    /*********************************************/


    if (v(0) != 0 ) { 
      
      printf("Nstates = %d \n",Nstates);
      
      printf("Nd = %d\n",Nd);
     
      printf("Constraint length K = %d\n",K);
      
      printf("Code generator matrix G, rows = %d, cols =%d\n", Ngr, Ngc);
      
      if (tb(0)==1)  printf("x +tail byte zeros = ");
      else printf("x = ");
      for (int i = 0; i < Ny; i++) {
        printf("%d",(int)y[i]);
      }
      printf("\n");
      
      printf("Nx = %d, Nout = %d\n",Nx,Nout);
      
    
      for (int i = 0; i < Ngr ; i++ ) {
        for (int j = 0; j < Ngc ; j++) {
          printf(" %d", (int)G(i,j));
        }
        printf("\n");
      }
      
      
      /*print b MSB(left) to LSB(right) */
      printf(" b = ");
      for (int i = 15; i >= 0; i-- ) {
        printf("%d", (b >> i) & 1);
      }
      printf("\n");
      
      
      /*************/
      /* state map */
      /*************/
      printf("input state nextstate output\n");
      printf("----- ----- ---------  ------\n");
      for (int i = 0; i < Nstates ; i ++ ) {
        for (int j = 0; j < 2; j++) {
          printf("  %d     %d       %d         ",j, i, state_next[i][j]); 
          for (int k = 0; k< Ngr; k++ ) {
            printf("%d ",state_out[i][j][k]);
          }
          printf("\n");
        }
      }    


    }  /**** done with verbose, print path metrics ****/ 
      

    
  } else {
     printf("Use Error\n");
  } /* end check for input consistency */

  return trout;
} /****** end rwviterbi *****/

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
