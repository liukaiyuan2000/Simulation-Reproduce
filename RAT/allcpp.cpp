#include <RcppArmadillo.h>
//using namespace arma;
//using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

arma::vec vecToRanks(arma::vec x) {
  int n=x.size();
  arma::uvec xSortedInds =  sort_index(x);
  arma::vec ranks = arma::linspace(0, n-1, n);
  arma::vec repeats = arma::linspace(0, n-1, n);
  int i,k = 1;
  double last = x[xSortedInds[0]];
  for (i = 0; i < n; i++) {
    double current = x[xSortedInds[i]];
    if (current != last) {
      k++;
      last = current;
    }
    repeats[i] = k;
  }
  int m=xSortedInds.size();
  for (int i = 0; i < m; i++) {
    ranks[xSortedInds[i]] = repeats[i];
  }
  return ranks;
}

// [[Rcpp::export]]
arma::mat outerminus(arma::vec X) {
  int n=X.size();
  arma::mat A=repmat(X,1,n); 
  arma::mat C=repmat(X.t(),n,1);  
  return A-C;
}


//[[Rcpp::export]]
arma::mat bubbleSortX(arma::vec x,arma::vec y) {
  int n=x.size();
  arma::mat  A = arma::ones<arma::mat>(n,2);
  int i,j;
  double tmpx,tmpy;
  for(i = n - 1; i > 0; i--) {
    for(j = 0; j < i; j++) {
      if(x(j+1)<x(j)){
        tmpx = x(j+1);
        x(j+1) = x(j);
        x(j) = tmpx;
        tmpy = y(j+1);
        y(j+1) = y(j);
        y(j) = tmpy;
      }
    }
    A(i,0)=x(i);
    A(i,1)=y(i);
    
  }
  
  A(0,0)=x(0);
  A(0,1)=y(0);
  return(A);
} 

//[[Rcpp::export]]
double mtau(arma::vec x,arma::vec y) {
  int n=x.size();
  arma::mat  A = bubbleSortX(x,y);
  double c=0;
  for(int i=0;i<n;i++){
    double Ri=0;
    for(int j=0;j<i;j++){
      if(A(j,1)<A(i,1))
        Ri=Ri+1;
    }
    c+=Ri;
  }
  return( 4*c/n/(n-1)-1.0);
} 


// [[Rcpp::export]]
arma::vec taucpp(arma::mat X){
  arma::vec TN = arma::zeros<arma::vec>(2);
  int p=X.n_cols, n=X.n_rows;
  arma::mat Sigma2=arma::zeros(p,p);
  for(int i=0;i<(p-1);i++){
    for(int j=(i+1);j<p;j++){
      Sigma2(i,j)= pow(mtau(X.col(i),X.col(j)),2);
    }
  }
  double Tn=0;
  Tn=sum(sum( Sigma2 ));
  double mu=p*(p-1)*(2*n+5)/9.0/n/(n-1);
  double sigmanp=4.0*p*(p-1)*(n-2)*(100*pow(n,3)+492*pow(n,2)+731*n+279)/2025.0/(pow(n,3)*pow(n-1,3));
  TN(0)=(Tn-mu)/sqrt(sigmanp);
  sigmanp=2.0*(2*n+5)/9.0/n/(n-1);
  TN(1)=max(max(Sigma2))/(sigmanp);
  return(TN);
}



// [[Rcpp::export]]
arma::vec allmain(arma::mat X){
  arma::vec TN = arma::zeros<arma::vec>(6);
  int n=X.n_rows;
  int p=X.n_cols;
  arma::mat Sigma=cor(X);
  Sigma.diag().zeros();

  //compute the statistic of rho
  arma::mat Sigma1=arma::zeros(p,p);
  arma::mat Sigma2=arma::zeros(p,p);
  double mu=1.0*p*(p-1)/(n-1)/2.0;
  arma::vec index1,index2;
  for(int i=0;i<(p-1);i++){
    for(int j=(i+1);j<p;j++){
      index1=vecToRanks(X.col(i));
      index2=vecToRanks(X.col(j));
      Sigma2(i,j)=1.0-6.0*sum(pow(index1-index2,2.0))/(pow(n,3.0)-n);
      Sigma1(i,j)=1.0-3.0*sum(abs(index1-index2))/(pow(n,2.0)-1.0); 
    }
  }
  
  
  
  //compute the statistic of sum Mao's rho(2017 CSTM)
  mu=1.0*p*(p-1)/(n-1)/2.0;
  double estd1=p*(p-1)*(25.0*pow(n,3)-57.0*pow(n,2)-40.0*n+108);
  double estd2=25.0*pow(n-1,3)*n*(n+1);
  double estd=sqrt(estd1/estd2);
  TN(0)=(sum(sum(pow(Sigma2,2.0)))-mu)/estd;
  
  //compute the statistic based on sum Mao's tau
  arma::vec Tn=taucpp(X);
  TN(1)=Tn(0);
  
  //compute the statistic based on sum Spearman's footrule  
  mu=1.0*p*(p-1)*(2*pow(n,2)+7)/(n+1)/pow(n-1,2)/10.0;
  estd1=1.0*p*(p-1)*(28*pow(n,5)-14*pow(n,4)+172*pow(n,3)+142*pow(n,2)+1159*n+1726);
  estd2=1.0*175*pow(n+1,3)*pow(n-1,4);
  estd=sqrt(estd1/estd2);
  TN(2)=(sum(sum(pow(Sigma1,2.0)))-mu)/estd;
  
  //compute the statistic of max rho
  estd=1.0/(n-1);
  TN(3)=max(max(pow(Sigma2,2.0)))/estd;
  
  //compute the statistic based on max tau  
  TN(4)=Tn(1);
  
  //compute the statistic based on max footrule  
  estd1=2.0*pow(n,2)+7.0;
  estd2=5.0*(n+1)*pow(n-1,2);
  double estd0=(estd1/estd2);
  TN(5)=max(max(pow(Sigma1,2.0)))/estd0;
  return(TN);
}



