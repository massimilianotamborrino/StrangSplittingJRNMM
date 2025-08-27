#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//------------------------------------------------------------------------------
// Author: Massimiliano Tamborrino
// Date: 2024-12-09
//
// Description: Strang splitting method for path simulation of the stochastic
//              multi-population JRNMM, proposed in the paper:
//
//             Network inference via approximate Bayesian computation.
//             Illustration on a stochastic multi-population neural mass model, by S. Ditlevsen, M. Tamborrino, I. Tubikanec, Ann. Appl. Stat. 2025
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Matrix-Vector Multiplication
//
// Input: matrix mat, vector vec
// Output: product of mat and vec
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector mv_mult_JRNMM_(NumericMatrix mat, NumericVector vec)
{
  NumericVector ret(mat.nrow());
  double temp=0;

  for(int i = 0; i < mat.nrow(); i++)
  {
    for(int j = 0; j < vec.size(); j++)
    {
      temp = temp + mat(i,j) * vec[j];
    }
    ret[i]=temp;
    temp=0;
  }
  return ret;
};

//------------------------------------------------------------------------------
//Sigmoid function of the JRNMM
//
// Input:
// x             value at which the sigmoid function should be evaluated
// vmax, v0, r   model parameters of the JRNMM
//
// Output:       value of sigmoid function evaluated in x
//------------------------------------------------------------------------------
// [[Rcpp::export]]
double sigmoid_JRNMM_Cpp_(double x, double vmax, double v0, double r)
{
  double ret=vmax/(1+exp(r*(v0-x)));
  return ret;
};

//------------------------------------------------------------------------------
// Diagonal entering into the block matrix of the exponential matrix of the linear SDE in the splitting scheme
// Exponential matrix of the underlying linear SDE in the JRNMM
//
// Input:
// dGamma vector of diagonal entries of the Gamma matrix
// h step size used for path simulation
//
// Output: See above

// [[Rcpp::export]]
NumericMatrix exponential_matrix_JRNMM_(NumericVector dGamma, double h){
  NumericVector dGh(dGamma.size()),expdG(dGamma.size());
  dGh=  dGamma*h;
  expdG=exp(-dGh);
  NumericMatrix M(4,dGamma.size());
  M(0,_)= expdG*(1+dGh);;
  M(1,_)= expdG*h;
  M(2,_)= -expdG*pow(dGh,2)/h;
  M(3,_)= expdG*(1-dGh);
  return M;
}

//------------------------------------------------------------------------------
// Covariance matrix of the linear SDE in the splitting scheme
// Input:
// dGamma vector of diagonal entries of the Gamma matrix
// h step size used for path simulation
// dSigma vector of diagonal entries of the Sigma matrix
// Output: See above
// [[Rcpp::export]]
NumericMatrix covariance_matrix_JRNMM_(NumericVector dGamma,NumericVector dSigma, double h){
  int n=dGamma.size();
  NumericVector dSigma2(n),dGh(n),exp2dGh(n),factor(n);
  dSigma2=pow(dSigma,2);
  dGh= dGamma*h;
  exp2dGh=exp(-2*dGh);factor=1-exp2dGh*pow(dGh,2);
  NumericMatrix COV(3,n);
  COV(0,_)= 0.25*pow(dGamma,-3)*dSigma2*(factor-exp2dGh*pow(1+dGh,2));
  COV(1,_)= 0.5*dSigma2*exp2dGh*pow(h,2);
  COV(2,_)= 0.25*1/dGamma*dSigma2*(factor-exp2dGh*pow(1-dGh,2));
  NumericMatrix COV2(2*n,2*n);
  for(int i=0;i<n;++i){
    COV2(i,i)= COV(0,i);
    COV2(i+n,i+n)= COV(2,i);
    COV2(i,i+n)=COV(1,i);
    COV2(i+n,i)=COV(1,i);}
  return(COV2);
}


//------------------------------------------------------------------------------
//linear SDE: exact solution of 1st subequation of splitting procedure for the JRNMM
//
// Input:
// vec       current value X^[1](ti) of the solution of the 1st subequation
// dm        exponential matrix appearing in the exact solution of the 1st subequation
// xi        normally distributed random vector xi appearing in the exact solution of the 1st subequation
//
// Output:   next value X^[1](ti+1) of the solution of the 1st subequation
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector linear_JRNMM_Cpp_(NumericVector vec, NumericMatrix dm, NumericVector xi)
{
  NumericVector ret=mv_mult_JRNMM_(dm,vec)+xi;
  return ret;
};

//------------------------------------------------------------------------------
//nonlinear ODE: exact solution of 2nd subequation of splitting procedure for the JRNMM
//
// Input:
// N        number of populations in the JRNMM
// vec      current value X^[2](ti) of the solution of the 2nd subequation
// h        step size for path simulation
// Theta    matrix containing continuous model parameters of the JRNMM
// Rho      matrix containing the coupling direction parameters of the JRNMM
// K        matrix containing the coupling strength parameters of the JRNMM
//
// Output:  next value X^[2](ti+1) of the solution of the 2nd subequation
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector nonlinear_JRNMM_Cpp_(int N, NumericVector vec, double h, NumericMatrix Theta, NumericMatrix Rho, NumericMatrix K)
{
  NumericVector help(vec.size());
  double sum=0;
  int l=(6*N)/2;

  for(int k = 0; k < N; k++)
  {
    help(3*k+l)=Theta(k,0)*Theta(k,2)*sigmoid_JRNMM_Cpp_(vec(3*k+1)-vec(3*k+2),Theta(k,8),Theta(k,6),Theta(k,7));
    help(3*k+(l+2))=Theta(k,1)*Theta(k,3)*0.25*Theta(k,4)*sigmoid_JRNMM_Cpp_(0.25*Theta(k,4)*vec(3*k),Theta(k,8),Theta(k,6),Theta(k,7));

    sum=0;
    for(int j = 0; j < N; j++)
    {
      if(j!=k){
        sum=sum+Rho(j,k)*K(j,k)*vec(3*j);
      }
    }
    help(3*k+(l+1))=Theta(k,0)*Theta(k,2)*(Theta(k,5)+0.8*Theta(k,4)*sigmoid_JRNMM_Cpp_(Theta(k,4)*vec(3*k),Theta(k,8),Theta(k,6),Theta(k,7))+sum);

  }

  NumericVector ret=vec+h*help;
  return ret;
};

//------------------------------------------------------------------------------
//Strang splitting method for the JRNMM
//
// Input:
// N        number of populations in the JRNMM
// grid     time grid for path simulation
// h        step size for path simulation
// startv   starting value X0 for path simulation
// dm       exponential matrix appearing in the exact solution of the 1st subequation
// meanVec  zero vector used to generate the random vector xi appearing in the exact solution of the 1st subequation
// covMat   covariance matrix (after Cholesky decomposition) appearing in the exact solution of the 1st subequation
// Theta    matrix containing continuous model parameters of the JRNMM
// Rho      matrix containing the coupling direction parameters of the JRNMM
// K        matrix containing the coupling strength parameters of the JRNMM
//
// Output:  path of the stochastic N-population JRNMM
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix JRNMM_Splitting_Cpp_(int N_i, NumericVector grid_i, double h_i, NumericVector startv, NumericMatrix dm_i, NumericVector meanVec_i, NumericMatrix covMat_i, NumericMatrix Theta_i, NumericMatrix Rho_i, NumericMatrix K_i)
{
  int N=N_i;
  double h=h_i;
  NumericVector start=startv;
  NumericVector grid=grid_i;
  int iter=grid.size();
  int dim=6*N;

  NumericMatrix dm=dm_i;
  NumericVector meanVec=meanVec_i;
  NumericMatrix covMat=covMat_i;
  NumericMatrix Theta=Theta_i;
  NumericMatrix Rho=Rho_i;
  NumericMatrix K=K_i;

  NumericMatrix sol(dim,iter);
  sol(_, 0)=start;
  NumericVector newv=start;

  Function f_rmvn("rmvn");
  NumericMatrix xi=f_rmvn(iter,meanVec,covMat,Named("isChol", true));
  NumericVector xi_rel;

  for(int i=1;i<iter;i++)
  {
    xi_rel=xi(i,_);
    newv=nonlinear_JRNMM_Cpp_(N,newv,h/2,Theta,Rho,K);
    newv=linear_JRNMM_Cpp_(newv,dm,xi_rel);
    newv=nonlinear_JRNMM_Cpp_(N,newv,h/2,Theta,Rho,K);
    sol(_,i)=newv;
  }

  NumericMatrix ret=sol;

  return sol;
};


// [[Rcpp::export]]
NumericVector meanHO_(int N, double h,NumericMatrix expM,NumericVector x){
  NumericVector vec(6*N);
  for(int i=0; i<3*N;++i){
    vec[i]=expM(0,i)*x[i]+expM(1,i)*x[3*N+i];
    vec[3*N+i]=expM(2,i)*x[i]+expM(3,i)*x[3*N+i];
  }
  return vec;
}

//------------------------------------------------------------------------------
// Faster code to simulate all components of a stochastic N-population JRNMM taking
// advantage of block-matrix decomposition

// Input:
// N        number of populations in the JRNMM
// grid     time grid for path simulation
// h        step size for path simulation
// startv   starting value X0 for path simulation
// dGamma   vector of diagonal entries of the Gamma matrix
// dSigma   vector of diagonal entries of the Sigma matrix
// Theta    matrix containing continuous model parameters of the JRNMM
// Rho      matrix containing the coupling direction parameters of the JRNMM
// K        matrix containing the coupling strength parameters of the JRNMM
//
// Output:  path of the stochastic N-population JRNMM
// [[Rcpp::export]]
NumericMatrix fastJRNMM_Splitting_Cpp_(int N, NumericVector grid, double h, NumericVector start, NumericVector dGamma,NumericVector dSigma, NumericMatrix Theta, NumericMatrix Rho, NumericMatrix K)
{
  int iter=grid.size();
  int dim=6*N;
  NumericMatrix sol(dim,iter);
  sol(_, 0)=start;
  NumericVector newv=start;
  NumericVector meanVec(dim);

  NumericMatrix expM(4,3*N);
  expM=exponential_matrix_JRNMM_(dGamma,h);

  NumericMatrix Cov(dim,dim);
  Cov=covariance_matrix_JRNMM_(dGamma,dSigma, h);

  Function f_rmvn("rmvn");
  NumericMatrix xi=f_rmvn(iter,meanVec,Cov);

  for(int i=1;i<iter;i++)
  {
    newv=nonlinear_JRNMM_Cpp_(N,newv,h/2,Theta,Rho,K);
    newv=meanHO_(N,h,expM,newv)+xi(i,_);
    newv=nonlinear_JRNMM_Cpp_(N,newv,h/2,Theta,Rho,K);
    sol(_,i)=newv;
  }
  return sol;
};

//------------------------------------------------------------------------------
// Function returning the observed components from a system of N populations of JRNMMs.
// Input:
// N    number of populations in the JRNMM
// sol  paths of the 6N coordinates of the N populations JRNMMs
// Output:
// N observed paths from the system

// [[Rcpp::export]]
NumericMatrix observedJRNMM_(int N, NumericMatrix sol)
{
  int iter=sol.cols();
  NumericMatrix Y(N,iter);
  int index_l, index_r;
  for(int j=0;j<N;j++){
    index_l=3*j+1;
    index_r=3*j+2;
    Y(j,_)=sol(index_l,_)-sol(index_r,_);
  }
  return(Y);}

//------------------------------------------------------------------------------
// Function to compute the K matrix given N, L, c
// Input:
// N number of populations of neurons
// L strength parameter
// c [0,1] parameter
// Output: Kmatrix given L and c

// [[Rcpp::export]]
NumericMatrix KmatrixgivenLc_(int N, double L, double c)
{
  double K_value;
  NumericMatrix K(N,N);
  for(int j=0;j<N-1;++j){
    K(j,j)=R_PosInf;
    for(int k=j+1;k<N;++k){
      K_value=pow(c,abs(k-j)-1)*L;
      K(j,k)=K_value;
      K(k,j)=K_value;
    }
    K(N-1,N-1)=R_PosInf;
  }
  return K;
};
