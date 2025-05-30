/* ========================================================================== */
/*                                                                            */
/*   Functional Local Linear Regression                                       */
/*   Functional Kernel Smoothing                                              */
/*   (c) 2019 Stanislav Nagy                                                  */
/*   29/07/2019, nagy@karlin.mff.cuni.cz                                      */
/*                                                                            */
/*   Description                                                              */
/*   Fast implementation of the Nadaraya - Watson kernel smoother, and the    */
/*   local linear smoother procedure for functional predictor and scalar      */
/*   response. Both the version with a chosen bandwidth, and the version with */
/*   leave-one-out cross-validation for bandwdith selection are provided.     */
/*   Based on BLAS and LAPACK routines.                                       */
/*                                                                            */
/*   Ferraty, F. and Nagy, S. (2019):                                         */ 
/*      Scalar-on-function local linear regression and beyond                 */
/*                                                                            */
/* ========================================================================== */

#include <R.h>
#include <R_ext/Lapack.h>

// Solves A x = b with the solution x clobbering b (LAPACK)
// first call DPOTRF and then DPOTRS

void matsolve(double *a, double *b, int *nrowb, int *ncolb, int *info)
{
    F77_CALL(dpotrf)("L", nrowb, a, nrowb, info);
    if (*info != 0)
        warning("Cholesky decomposition failed");
    F77_CALL(dpotrs)("L", nrowb, ncolb, a, nrowb, b, nrowb, info);
    if (*info != 0)
        warning("solution failed");
}

double kern(double t, int k) {
   // positive kernel functions
   // k is the switch between kernels
   double result;
   double pi = atan(1)*4;

   result = 0;
   if(k==1){ // uniform kernel
   if ((t < 0)||(t > 1))
      result = 0;
   else
      result = 1;
      }
   if(k==2){ // triangular kernel
   if ((t < 0)||(t > 1))
      result = 0;
   else
      result = (1-t)*2;
      } 
   if(k==3){ // Epanechnikov kernel
   if ((t < 0)||(t > 1))
      result = 0;
   else
      result = (1-pow(t,2))*3/2;
      }
   if(k==4){ // biweight kernel
   if ((t < 0)||(t > 1))
      result = 0;
   else
      result = pow(1-pow(t,2),2)*15/8;
      }
   if(k==5){ // triweight kernel
   if ((t < 0)||(t > 1))
      result = 0;
   else
      result = pow(1-pow(t,2),3)*35/16;
      }
   if(k==6){ // Gaussian kernel
   if (t < 0)
      result = 0;
   else
      result = sqrt(2/pi)*exp(-pow(t,2)/2);
      }
   return result; 
}

void ksm(double *D, double *Y, double *h, int *n, int *m, int *kernI, double *res)
{
  int i, j;
  double den;
  
  for(j=0;j<*m;j++){
    res[j] = 0;
    den = 0;
    if(h[j]==-1){ // bandwidth -1 is the code for h=infinity, i.e. just the mean of Y
      for(i=0;i<*n;i++){
       res[j] = res[j] + Y[i];
       den = den + 1;
      }} else {
      for(i=0;i<*n;i++){
       res[j] = res[j] + Y[i]*kern(D[j + i*(*m)]/h[j], kernI[0]);
       den = den + kern(D[j + i*(*m)]/h[j], kernI[0]);
      }
    }
    res[j] = res[j]/den;
  }
}

void ksm_cv(double *D, double *Y, int *n, double *H, int *nH, double *CV, int *nCV, int *kernI)
{
  int i,j,k,hi;
  double Di[*n-1], Yi[*n-1];
  double h, Yt; //, CVmin; //hmin;
  int n1 = *n - 1, one = 1;

  // for each bandwidth h leave-one-out cross-validation
  // CVmin = pow(10,100);
  // hmin = 0;
  for(hi=0;hi<*nH;hi++){
    CV[hi] = 0;
    for(i=0;i<*nCV;i++){
      k = 0;
      for(j=0;j<*n;j++){
        if(i!=j){
          Di[k] = D[j * (*n) + i];
          Yi[k] = Y[j];
          k++;
          }
      }
      h = H[hi * (*n) + i];
      ksm(Di, Yi, &h, &n1, &one, kernI, &Yt);
      CV[hi] += pow(Yt-Y[i],2);
    }
    CV[hi] = CV[hi]/(*nCV);
    //if(CV[hi]<CVmin){
    //  CVmin = CV[hi];
      // hmin = h;
    // }
  }
  // cross-validated bandwidth estimator 
  // ksm(Dnew, Y, &hmin, n, &one, kernI, res);
}

void ksm_der(double *C, double *D, double *Dnew, double *Y, double *h1, double *h2, int *n, int *nTest, int *J, int *kernI, double *res)
// kernel-based derivative, Hall et al. (2009)
// C[n][J]; D[n][n]; Dnew[nTest][n]; Y[n]; h1[nTest]; h2[nTest]
{
  int i,j,i1,i2;
  double term, num, den;
  //double Q[*n][*n][*J];
  double ksi[*n][*n][*J];
  double Yii[*n][*n];
  
  for(i1=0;i1<*n;i1++){
    for(i2=i1;i2<*n;i2++){
      Yii[i1][i2] = Y[i1] - Y[i2];
      Yii[i2][i1] = -Yii[i1][i2];
      for(j=0;j<*J;j++){
        // Q[i1][i2][j] = kern((1 - pow(C[i1 + j*(*n)]-C[i2 + j*(*n)],2)/pow(D[i1 + i2*(*n)],2))/h2[0], kernI[0]);
        //Q[i2][i1][j] = Q[i1][i2][j];
        ksi[i1][i2][j] = C[i1 + j*(*n)]-C[i2 + j*(*n)];
        ksi[i2][i1][j] = -ksi[i1][i2][j]; 
       }
     }
   } 
   for(i=0;i<*nTest;i++){
    if((h1[i]>0) & (h2[i]>0)){
      for(j=0;j<*J;j++){
        num = 0.0;
        den = 0.0;
        for(i1=0;i1<*n;i1++){ 
          for(i2=0;i2<*n;i2++){
            if(ksi[i1][i2][j]>0){
              term = kern(Dnew[i+i1*(*nTest)]/h1[i],kernI[0])*kern(Dnew[i+i2*(*nTest)]/h1[i],kernI[0])*kern((1 - pow(ksi[i1][i2][j]/D[i1 + i2*(*n)],2))/h2[i], kernI[0]);
              num = num + Yii[i1][i2]*term;
              den = den + ksi[i1][i2][j]*term;
            }
          }  
        } 
      res[i + j*(*nTest)] = num/den;
      }
    } else { // if the bandwidth is set to Infty, derivatives are all zero (constant regression)
      for(j=0;j<*J;j++) res[i + j*(*nTest)] = 0.0;
    }                        
   }
}


#include <R_ext/BLAS.h>

void matvecmult(double *a, double *b, int *nrow, int *ncol, double *result)
{
  double zero = 0.0;
  double one = 1.0;
  int ione = 1;
  F77_CALL(dgemv)("n", nrow, ncol, &one, a, nrow, b, &ione, &zero, result,
           &ione);
}

void matmatmult(double *a, double *b, int *nrowa, int *ncola, int *ncolb,
                double *c)
{
  double one = 1.0;
  double zero = 0.0;
  F77_CALL(dgemm)("n", "n", nrowa, ncolb, ncola, &one, a, nrowa, b, ncola,
           &zero, c, nrowa);
}

void llsm_single(double *C, double *Cnew, double *D, double *Y, double *h, int *n, int *m, int *J, int *kernI, double *res)
{
  int i, j, k, ione = 1, info = 0;
  double C2[*J * *n], M1[*J * *J], C3[*n * *J];
  double C1[(*n) * (*J)], restemp[*J];
    
  for(k=0;k<*m;k++){ // for each function
    // create matrix with rows coefficients corresponding to X - x
    for(i=0;i<*n;i++){
      C1[0 + i] = 1; // C1[i,1] = 1
      for(j=1;j<*J;j++){
        C1[j * (*n) + i] = C[j * (*n) + i] - Cnew[k + j*(*m)]; // C1[i,j] = C[i,j] - Cnew[k,j], i>1
        }}
    for(j=0;j<*J;j++) restemp[j] = 0.0;
    for(i=0;i<*n;i++)
    for(j=0;j<*J;j++){
      // C2 is C1 transposed, and cut at jj
      // C3 is C1, cut at jj
    if(h[k]==-1){
      C2[i * (*J) + j] = C1[j * (*n) + i] * 1;
      } else C2[i * (*J) + j] = C1[j * (*n) + i] * kern(D[k+i*(*m)]/h[k],kernI[0]);
        C3[j * (*n) + i] = C1[j * (*n) + i];
      }
    // matrix multiplication
    matmatmult(C2, C3, J, n, J, M1);
    matvecmult(C2, Y, J, n, restemp);
    matsolve(M1, restemp, J, &ione, &info);
    if(info!=0) restemp[0] = R_PosInf;
    for(j=0;j<*J;j++) res[k+j*(*m)] = restemp[j];
    }
}

void llsm(double *C, double *Cnew, double *D, double *Y, double *h, int *n, int *J, int *kernI, double *res)
{
  int i, j, ione = 1, info = 0;
  int ii, jj;
  double C1[(*n) * (*J)];

  // create matrix with rows coefficients corresponding to X - x
  for(i=0;i<*n;i++){
    C1[0 + i] = 1;
    for(j=1;j<*J;j++){
      C1[j * (*n) + i] = C[j * (*n) + i] - Cnew[j];
      }}
  // check for possible problems with matrix K (if close to non-regular matsolve does not work)
  //for(i=0;i<*n;i++)
  //  if(kern(D[i]/h[0])<pow(10,-75)) info = 1;
  // create matrix C %*% diag(kern(D/h))
  for(jj=0;jj<*J;jj++){ // for all sizes of matrices from 1 to J
    int jjb = jj+1;
    double C2[jjb * *n], M1[jjb * jjb], C3[*n * jjb], tres[jjb];
     //if(info==0){
      for(i=0;i<*n;i++)
        for(j=0;j<jjb;j++){
          // C2 is C1 transposed, and cut at jj
          // C3 is C1, cut at jj
          if(h[0]==-1){ // bandwidth -1 is the code for h=infinity, i.e. the usual linear regression
            C2[i * (jjb) + j] = C1[j * (*n) + i] * 1;
            } else C2[i * (jjb) + j] = C1[j * (*n) + i] * kern(D[i]/h[0],kernI[0]);
          C3[j * (*n) + i] = C1[j * (*n) + i];
        }
      // matrix multiplication
      matmatmult(C2, C3, &jjb, n, &jjb, M1);
      matvecmult(C2, Y, &jjb, n, tres);
      matsolve(M1, tres, &jjb, &ione, &info);
      if(info==0){
        for(ii=0;ii<*J;ii++) if(ii<jjb) res[ii*(*J) + jj] = tres[ii]; else res[ii*(*J)+jj] = R_NaReal;
      } else res[0*(*J) + jj] = R_PosInf;
    //} else res[0*(*J) + jj] = R_PosInf; // Inf if we have an ill-conditioned matrix
  }
}

void llsm_cv(double *C, double *D, double *Y, int *n, int *J, double *H, int *Hind, int *nH, double *CV, double *CVB, double *CVD, 
     int *nCV, int *kernI, double *res)
     // H is a matrix n*nH, each function has its own set of bandwidths
     //, double *Yt, double *Ci, double *Yi, double *Di, double *Cinew)
             //double *Cnew, double *Dnew
{
  int i,j,k,l,hi;
  double Di[*n-1], Ci[(*n-1) * (*J)],  Yi[*n-1], Cinew[*J];
  double h, h2, CVmin[*J], Yt[*J * (*J)]; //hmin[*J] 
  int n1 = *n - 1;
  // for cross-validation with derivatives
  int k2,k3,kdif,ndif;
  double Ydif, CVDmin[*J], CVt;
  
  // for each bandwidth h leave-one-out cross-validation
  for(j=0;j<*J;j++){
    CVmin[j] = R_PosInf;
    CVDmin[j] = R_PosInf;
    // hmin[j] = 0; 
  }
  for(hi=0;hi<*nH;hi++){
    for(j=0;j<*J;j++){
      CV[j*(*nH)+hi] = 0;
      CVB[j*(*nH)+hi] = 0;
      CVD[j*(*nH)+hi] = 0;
    }
    for(i=0;i<*nCV;i++){ 
      // cross-validation for first nCV functions only in a leave-one-out sense
      for(j=0;j<*J;j++){
        Cinew[j] = C[j * (*n) + i];
      }
      k = 0;
      for(j=0;j<*n;j++){
        if(j!=i){
          Di[k] = D[j * (*n) + i];
          Yi[k] = Y[j];
          for(l=0;l<*J;l++){
            Ci[l * (*n-1) + k] = C[l * (*n) + j];
          }
          k++;
        }
      }
      h = H[hi * (*n) + i];
      if(h==-1) h2 = R_PosInf; else h2 = h; 
      llsm(Ci,Cinew,Di,Yi,&h,&n1,J,kernI,Yt);  // Yt is a matrix J times J
      // first, we cross-validate only for the intercept term (the others represent the derivatives)
      for(k=0;k<*J;k++){
        CV[k*(*nH) + hi] += pow(Yt[0+k]-Y[i],2);  // MSE
        CVB[k*(*nH) + hi] += Yt[0+k]-Y[i];        // bias
        // cross-validation for derivative now
        if(Hind[hi * (*n) + i] == 1){ // if cross-validation for derivative also for this h should be performed
          CVt = 0.0;
          ndif = 0;
          for(k2=0;k2<(*n-1);k2++){
              for(k3=k2+1;k3<(*n);k3++){  // for all X_{k2} and X_{k3} with k2<k3
                if((k2!=i)&&(k3!=i)&&(D[k2*(*n)+k3]>0)&&(D[i*(*n)+k2]<h2)&&(D[i*(*n)+k3]<h2)){  // the last two conditions replaced by weights in CVt below
                // if ||X_{k2}-X_{k3}||>0 and ||X_{k2}-X_i||<h and ||X_{k3}-X_i||<h
                  Ydif = 0.0;
                  if(k>0) { for(kdif=1;kdif<=k;kdif++){ // not the intercept term, only derivatives (kdif>=1)
                    Ydif += Yt[kdif*(*J)+k]*(C[kdif*(*n)+k2] - C[kdif*(*n)+k3]);
                    }  // estimate the derivative of m at X_i in direction (X_{k2}-X_{k3})/||X_{k2}-X_{k3}||
                    } else {Ydif = 0.0;} // in the intercept only model (usual linear model), estimated derivative is zero
                  CVt += pow((Y[k2]-Y[k3])/(D[k2*(*n)+k3]) - Ydif/(D[k2*(*n)+k3]),2);
                  ndif++;
                  }
              }
          }
        } else {
        CVt = 0.0;
        ndif = 0;
        }
        CVD[k*(*nH)+hi] += CVt/ndif; // gives Inf if ndif = 0 (no couple of points in the neighbourhood)
      }
    }
    for(j=0;j<*J;j++){
      CV[j*(*nH) + hi] = CV[j*(*nH) + hi]/(*nCV);
      CVB[j*(*nH) + hi] = CV[j*(*nH) + hi] - pow(CVB[j*(*nH) + hi]/(*nCV),2); // variance (MSE - bias^2)
      CVD[j*(*nH) + hi] = CVD[j*(*nH) + hi]/(*nCV);
      if(CV[j*(*nH)+hi]<R_PosInf) {if(CV[j*(*nH) + hi]<CVmin[j]){
        CVmin[j] = CV[j*(*nH) + hi];
        }}
      if(CVD[j*(*nH)+hi]<R_PosInf) {if(CVD[j*(*nH) + hi]<CVDmin[j]){
        CVDmin[j] = CVD[j*(*nH) + hi];
        // res[0]=1.1;
        }}
      /*if(CV[j*(*nH) + hi]<CVmin[j]){
        CVmin[j] = CV[j*(*nH) + hi];
        }
      if(CVD[j*(*nH) + hi]<CVDmin[j]){
        CVDmin[j] = CVD[j*(*nH) + hi];
        }*/
    }
  }
}

///
/// Leave procedures
///

void llsm_leave(double *C, double *Cnew, int *Cint, double *D, double *Y, double *h, int *n, int *J, int *kernI, double *res)
{
  int i, j, ione = 1, info = 0;
  int ii, jj;
  double C1[(*n) * (*J)];
  double cnst;

  // create matrix with rows coefficients corresponding to X - x
  for(i=0;i<*n;i++){
    C1[0 + i] = 1;
    for(j=1;j<*J;j++){
      C1[j * (*n) + i] = C[j * (*n) + i] - Cnew[j];
      }}
  // check for possible problems with matrix K (if close to non-regular matsolve does not work)
  //for(i=0;i<*n;i++)
  //  if(kern(D[i]/h[0])<pow(10,-75)) info = 1;
  // create matrix C %*% diag(kern(D/h))
  for(jj=0;jj<*J;jj++){ // for all sizes of matrices from 1 to J
    int jjb = jj+1;
    double C2[jjb * *n], M1[jjb * jjb], C3[*n * jjb], tres[jjb];
     //if(info==0){
      for(i=0;i<*n;i++){       
        if(i==(*Cint)) cnst = 0.0; else cnst = 1.0;  // if i-th function is the Cint-th function this contributes 0
        for(j=0;j<jjb;j++){
          // C2 is C1 transposed, and cut at jj
          // C3 is C1, cut at jj
          if(h[0]==-1){ // bandwidth -1 is the code for h=infinity, i.e. the usual linear regression
            C2[i * (jjb) + j] = C1[j * (*n) + i] * 1;
            } else C2[i * (jjb) + j] = C1[j * (*n) + i] * kern(D[i]/h[0],kernI[0]) * cnst;
          C3[j * (*n) + i] = C1[j * (*n) + i];
        }
        }
      // matrix multiplication
      matmatmult(C2, C3, &jjb, n, &jjb, M1);
      matvecmult(C2, Y, &jjb, n, tres);
      matsolve(M1, tres, &jjb, &ione, &info);
      if(info==0){
        for(ii=0;ii<*J;ii++) if(ii<jjb) res[ii*(*J) + jj] = tres[ii]; else res[ii*(*J)+jj] = R_NaReal;
      } else res[0*(*J) + jj] = R_PosInf;
    //} else res[0*(*J) + jj] = R_PosInf; // Inf if we have an ill-conditioned matrix
  }
}

void llsm_single_leave(double *C, double *Cnew, double *D, double *Y, double *h, int *n, int *m, int *J, int *kernI, double *res)
// just as llsm_single, but with left out j-th observation from Cnew when evaluating prediction for the j-th observation from C
// to be used in the leave-one-out cross-validation when Cnew is the same as C
{
  int i, j, k, ione = 1, info = 0;
  double C2[*J * *n], M1[*J * *J], C3[*n * *J];
  double C1[(*n) * (*J)], restemp[*J];
  double cnst;
    
  for(k=0;k<*m;k++){ // for each function
    // create matrix with rows coefficients corresponding to X - x
    for(i=0;i<*n;i++){
      C1[0 + i] = 1;
      for(j=1;j<*J;j++){
        C1[j * (*n) + i] = C[j * (*n) + i] - Cnew[k + j*(*m)];
        }}
    for(j=0;j<*J;j++) restemp[j] = 0.0;
    for(i=0;i<*n;i++){
      if(i==k) cnst = 0.0; else cnst = 1.0;  // if i-th function is the same as k-th function this contributes 0
    for(j=0;j<*J;j++){
      // C2 is C1 transposed, and cut at jj
      // C3 is C1, cut at jj
    if(h[k]==-1){
      C2[i * (*J) + j] = C1[j * (*n) + i] * 1;
      } else C2[i * (*J) + j] = C1[j * (*n) + i] * kern(D[k+i*(*m)]/h[k],kernI[0]) * cnst;
        C3[j * (*n) + i] = C1[j * (*n) + i];
      }
      }
    // matrix multiplication
    matmatmult(C2, C3, J, n, J, M1);
    matvecmult(C2, Y, J, n, restemp);
    matsolve(M1, restemp, J, &ione, &info);
    if(info!=0) restemp[0] = R_PosInf;
    for(j=0;j<*J;j++) res[k+j*(*m)] = restemp[j];
    }
}

void llsm_cv_leave(double *C, double *D, double *Y, int *n, int *J, double *H, int *Hind, int *nH, double *CV, double *CVB, double *CVD, 
     int *nCV, int *kernI, double *res)
     // H is a matrix n*nH, each function has its own set of bandwidths
     //, double *Yt, double *Ci, double *Yi, double *Di, double *Cinew)
             //double *Cnew, double *Dnew
{
  int i,j,k,l,hi;
  double Di[*n-1], Ci[(*n-1) * (*J)],  Yi[*n-1], Cinew[*J];
  double h, h2, CVmin[*J], Yt[*J * (*J)]; //hmin[*J] 
  int n1 = *n - 1;
  // for cross-validation with derivatives
  int k2,k3,kdif,ndif;
  double Ydif, CVDmin[*J], CVt;
  
  // for each bandwidth h leave-one-out cross-validation
  for(j=0;j<*J;j++){
    CVmin[j] = R_PosInf;
    CVDmin[j] = R_PosInf;
    // hmin[j] = 0; 
  }
  for(hi=0;hi<*nH;hi++){
    for(j=0;j<*J;j++){
      CV[j*(*nH)+hi] = 0;
      CVB[j*(*nH)+hi] = 0;
      CVD[j*(*nH)+hi] = 0;
    }
    for(i=0;i<*nCV;i++){ 
      // cross-validation for first nCV functions only in a leave-one-out sense
      for(j=0;j<*J;j++){
        Cinew[j] = C[j * (*n) + i];
      }
      k = 0;
      for(j=0;j<*n;j++){
        if(j!=i){
          Di[k] = D[j * (*n) + i];
          Yi[k] = Y[j];
          for(l=0;l<*J;l++){
            Ci[l * (*n-1) + k] = C[l * (*n) + j];
          }
          k++;
        }
      }
      h = H[hi * (*n) + i];
      if(h==-1) h2 = R_PosInf; else h2 = h; 
      llsm_leave(Ci,Cinew,&i,Di,Yi,&h,&n1,J,kernI,Yt);  // Yt is a matrix J times J
      // first, we cross-validate only for the intercept term (the others represent the derivatives)
      for(k=0;k<*J;k++){
        CV[k*(*nH) + hi] += pow(Yt[0+k]-Y[i],2);  // MSE
        CVB[k*(*nH) + hi] += Yt[0+k]-Y[i];        // bias
        // cross-validation for derivative now
        if(Hind[hi * (*n) + i] == 1){ // if cross-validation for derivative also for this h should be performed
          CVt = 0.0;
          ndif = 0;
          for(k2=0;k2<(*n-1);k2++){
              for(k3=k2+1;k3<(*n);k3++){  // for all X_{k2} and X_{k3} with k2<k3
                if((k2!=i)&&(k3!=i)&&(D[k2*(*n)+k3]>0)&&(D[i*(*n)+k2]<h2)&&(D[i*(*n)+k3]<h2)){  // the last two conditions replaced by weights in CVt below
                // if ||X_{k2}-X_{k3}||>0 and ||X_{k2}-X_i||<h and ||X_{k3}-X_i||<h
                  Ydif = 0.0;
                  if(k>0) { for(kdif=1;kdif<=k;kdif++){ // not the intercept term, only derivatives (kdif>=1)
                    Ydif += Yt[kdif*(*J)+k]*(C[kdif*(*n)+k2] - C[kdif*(*n)+k3]);
                    }  // estimate the derivative of m at X_i in direction (X_{k2}-X_{k3})/||X_{k2}-X_{k3}||
                    } else {Ydif = 0.0;} // in the intercept only model (usual linear model), estimated derivative is zero
                  CVt += pow((Y[k2]-Y[k3])/(D[k2*(*n)+k3]) - Ydif/(D[k2*(*n)+k3]),2);
                  ndif++;
                  }
              }
          }
        } else {
        CVt = 0.0;
        ndif = 0;
        }
        CVD[k*(*nH)+hi] += CVt/ndif; // gives Inf if ndif = 0 (no couple of points in the neighbourhood)
      }
    }
    for(j=0;j<*J;j++){
      CV[j*(*nH) + hi] = CV[j*(*nH) + hi]/(*nCV);
      CVB[j*(*nH) + hi] = CV[j*(*nH) + hi] - pow(CVB[j*(*nH) + hi]/(*nCV),2); // variance (MSE - bias^2)
      CVD[j*(*nH) + hi] = CVD[j*(*nH) + hi]/(*nCV);
      if(CV[j*(*nH)+hi]<R_PosInf) {if(CV[j*(*nH) + hi]<CVmin[j]){
        CVmin[j] = CV[j*(*nH) + hi];
        }}
      if(CVD[j*(*nH)+hi]<R_PosInf) {if(CVD[j*(*nH) + hi]<CVDmin[j]){
        CVDmin[j] = CVD[j*(*nH) + hi];
        // res[0]=1.1;
        }}
      /*if(CV[j*(*nH) + hi]<CVmin[j]){
        CVmin[j] = CV[j*(*nH) + hi];
        }
      if(CVD[j*(*nH) + hi]<CVDmin[j]){
        CVDmin[j] = CVD[j*(*nH) + hi];
        }*/
    }
  }
}


/// single version of llsm_cv (only for one J)

void llsm_cv_single(double *C, double *D, double *Y, int *n, int *J, double *H, int *nH, double *CV, double *CVB, 
     int *nCV, int *kernI, double *res)//, //)
     // H is a matrix n*nH, each function has its own set of bandwidths
     // double *Yt, double *Ci, double *Yi, double *Di, double *Cinew)
             //double *Cnew, double *Dnew
{
  int i,j,k,l,hi;
  double Di[*n-1], Ci[(*n-1) * (*J)],  Yi[*n-1], Cinew[*J];
  double h, Yt[*J];
  int n1 = *n - 1, one = 1;
  
  // for each bandwidth h leave-one-out cross-validation
  for(hi=0;hi<*nH;hi++){
    for(j=0;j<*J;j++){
      CV[j*(*nH)+hi] = 0;
      CVB[j*(*nH)+hi] = 0;
    }
    for(i=0;i<*nCV;i++){ 
      // cross-validation for first nCV functions only in a leave-one-out sense
      for(j=0;j<*J;j++){
        Cinew[j] = C[j * (*n) + i];
      }
      k = 0;
      for(j=0;j<*n;j++){
        if(j!=i){
          Di[k] = D[j * (*n) + i];
          Yi[k] = Y[j];
          for(l=0;l<*J;l++){
            Ci[l * (*n-1) + k] = C[l * (*n) + j];  // Ci[k,l] = C[j,l]
          } 
          k++;
        }
      }  
      h = H[hi * (*n) + i];
      llsm_single(Ci,Cinew,Di,Yi,&h,&n1,&one,J,kernI,Yt);  // only Yt[0] is used, but Yt has J elements
      // cross-validate
      k=*J-1;
      CV[k*(*nH) + hi] += pow(Yt[0]-Y[i],2);  // MSE
      CVB[k*(*nH) + hi] += Yt[0]-Y[i];        // bias
    }  
    j=*J-1;
    CV[j*(*nH) + hi] = CV[j*(*nH) + hi]/(*nCV);
    CVB[j*(*nH) + hi] = CV[j*(*nH) + hi] - pow(CVB[j*(*nH) + hi]/(*nCV),2); // variance (MSE - bias^2)
  } 
}

void lqsm_single(double *C, double *Cnew, double *D, double *Y, double *h, int *n, int *m, int *kernI, double *res)
// local quadratic regression, for derivatives in the Mueller-Yao procedure based on kNN estimation
{
  int i, j, k, ithree = 3, ione = 1, info = 0;
  double C2[3 * *n], M1[3*3], C3[*n * 3];
  double C1[(*n) * 3], restemp[3];
    
  for(k=0;k<*m;k++){ // for each function
    // create matrix with rows coefficients corresponding to X - x
    for(i=0;i<*n;i++){
      C1[0 + i] = 1; // C1[i,1] = 1
      C1[1 * (*n) + i] = C[1 * (*n) + i] - Cnew[k + 1*(*m)]; // C1[i,2] = C[i,2] - Cnew[k,2]
      C1[2 * (*n) + i] = pow(C1[1 * (*n) + i],2); // C1[i,3] = C1[i,2]^2
      }
    for(j=0;j<3;j++) restemp[j] = 0.0;
    for(i=0;i<*n;i++)
    for(j=0;j<3;j++){
      // C2 is C1 transposed, and cut at jj
      // C3 is C1, cut at jj
    if(h[k]==-1){
      C2[i * 3 + j] = C1[j * (*n) + i] * 1;
      } else C2[i * 3 + j] = C1[j * (*n) + i] * kern(D[k+i*(*m)]/h[k],kernI[0]);
        C3[j * (*n) + i] = C1[j * (*n) + i];
      }
    // matrix multiplication
    matmatmult(C2, C3, &ithree, n, &ithree, M1);
    matvecmult(C2, Y, &ithree, n, restemp);
    matsolve(M1, restemp, &ithree, &ione, &info);
    if(info!=0) restemp[0] = R_PosInf;
    for(j=0;j<3;j++) res[k+j*(*m)] = restemp[j];
    }
}

void lqsm_cv_single(double *C, double *D, double *Y, int *n, int *J, double *H, int *nH, double *CV, double *CVB, 
     int *nCV, int *kernI, double *res)//, //)
     // H is a matrix n*nH, each function has its own set of bandwidths
     // double *Yt, double *Ci, double *Yi, double *Di, double *Cinew)
             //double *Cnew, double *Dnew
{
  int i,j,k,l,hi;
  double Di[*n-1], Ci[(*n-1) * (*J)],  Yi[*n-1], Cinew[*J];
  double h, Yt[*J];
  int n1 = *n - 1, one = 1;
  
  // for each bandwidth h leave-one-out cross-validation
  for(hi=0;hi<*nH;hi++){
    for(j=0;j<*J;j++){
      CV[j*(*nH)+hi] = 0;
      CVB[j*(*nH)+hi] = 0;
    }
    for(i=0;i<*nCV;i++){ 
      // cross-validation for first nCV functions only in a leave-one-out sense
      for(j=0;j<*J;j++){
        Cinew[j] = C[j * (*n) + i];
      }
      k = 0;
      for(j=0;j<*n;j++){
        if(j!=i){
          Di[k] = D[j * (*n) + i];
          Yi[k] = Y[j];
          for(l=0;l<*J;l++){
            Ci[l * (*n-1) + k] = C[l * (*n) + j];  // Ci[k,l] = C[j,l]
          } 
          k++;
        }
      }  
      h = H[hi * (*n) + i];
      lqsm_single(Ci,Cinew,Di,Yi,&h,&n1,&one,kernI,Yt);  // only Yt[0] is used, but Yt has J elements
      // cross-validate
      k=*J-1;
      CV[k*(*nH) + hi] += pow(Yt[0]-Y[i],2);  // MSE
      CVB[k*(*nH) + hi] += Yt[0]-Y[i];        // bias
    }  
    j=*J-1;
    CV[j*(*nH) + hi] = CV[j*(*nH) + hi]/(*nCV);
    CVB[j*(*nH) + hi] = CV[j*(*nH) + hi] - pow(CVB[j*(*nH) + hi]/(*nCV),2); // variance (MSE - bias^2)
  } 
}