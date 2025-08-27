#' @useDynLib StrangSplittingJRNMM
#' @importFrom Rcpp sourceCpp evalCpp
NULL

#'@rdname JRNMM_Splitting
#'@title JRNMM_Splitting
#'@description Strang splitting method for path simulation of the
#' stochastic multi-population JRNMM
#'@param N number of populations of neurons
#'@param grid time points at which to simulate the process
#'@param h step size used for path simulation
#'@param startv starting point x0
#'@param dm exponential matrix coming from the linear SDE in the splitting scheme
#'@param meanVec 6N dimensional vector of 0s
#'@param covMat covariance matrix coming from the linear SDE in the splitting scheme
#'@param Theta  vector of continuous parameters
#'@param Rho matrix of {0,1} discrete parameters
#'@param K Matrix of strength parameters
#'@return path of the stochastic N-population JRNMM
#'@export
JRNMM_Splitting <- function(N, grid, h, startv, dm, meanVec, covMat, Theta, Rho, K){
  return(JRNMM_Splitting_Cpp_(N, grid, h, startv, dm, meanVec, covMat, Theta, Rho, K))
}

#'@rdname fast_JRNMM_Splitting
#'@title fast_JRNMM_Splitting
#'@description (faster version of the) Strang splitting method for path simulation of the
#' stochastic multi-population JRNMM
#'@param N number of populations of neurons
#'@param grid time points at which to simulate the process
#'@param h step size used for path simulation
#'@param startv starting point x0
#'@param dGamma vector of diagonal entries of the Gamma matrix
#'@param dSigma vector of diagonal entries of the Sigma matrix
#'@param Theta  vector of continuous parameters
#'@param Rho matrix of {0,1} discrete parameters
#'@param K Matrix of strength parameters
#'@return path of the stochastic N-population JRNMM
#'@export
fast_JRNMM_Splitting <- function(N, grid, h, startv, dGamma,dSigma, Theta, Rho, K){
  return(fastJRNMM_Splitting_Cpp_(N, grid, h, startv, dGamma,dSigma, Theta, Rho, K))
}

#'@rdname observedJRNMM
#'@title observedJRNMM
#'@description Strang splitting method for path simulation of the observed components of the
#' stochastic multi-population JRNMM.
#'@param N number of populations of neurons
#'@param grid time points at which to simulate the process
#'@param h step size used for path simulation
#'@param startv starting point x0
#'@param dGamma vector of diagonal entries of the Gamma matrix
#'@param dSigma vector of diagonal entries of the Sigma matrix
#'@param Theta  vector of continuous parameters
#'@param Rho matrix of {0,1} discrete parameters
#'@param K Matrix of strength parameters
#'@return observed components of the stochastic N-population JRNMM
#'@export
observedJRNMM <- function(N, grid, h, startv, dGamma,dSigma, Theta, Rho, K){
  Y<-fast_JRNMM_Splitting(N, grid, h, startv, dGamma,dSigma, Theta, Rho, K)
  return(observedJRNMM_(N, Y))
}

#'@rdname exponential_matrix_JRNMM
#'@title exponential_matrix_JRNMM
#'@description Calculation of the diagonal entering into the block matrix of
# the exponential matrix of the linear SDE in the splitting scheme
#'@param dGamma vector of diagonal entries of the Gamma matrix
#'@param h step size used for path simulation
#'@return See description
#'@export
exponential_matrix_JRNMM <- function(dGamma, h){
  return(exponential_matrix_JRNMM_(dGamma, h))
}

#'@rdname covariance_matrix_JRNMM
#'@title covariance_matrix_JRNMM
#'@description Covariance matrix of the linear SDE in the splitting scheme
#'@param dGamma vector of diagonal entries of the Gamma matrix
#'@param h step size used for path simulation
#'@param dSigma vector of diagonal entries of the Sigma matrix
#'@return See description
#'@export
covariance_matrix_JRNMM <- function(dGamma,dSigma, h){
  return(covariance_matrix_JRNMM_(dGamma,dSigma, h))
}


#'@rdname mv_mult_JRNMM
#'@title mv_mult_JRNMM
#'@description Matrix-Vector Multiplication
#'@param mat matrix
#'@param vet vector
#'@return See description
#'@export
mv_mult_JRNMM <- function(mat,vet){
  return(mv_mult_JRNMM_(mat,vet))
}

#'@rdname sigmoid_JRNMM
#'@title sigmoid_JRNMM
#'@description Sigmoid function of the JRNMM
#'@param x value at which the sigmoid function should be evaluated
#'@param vmax value of the JRNMM
#'@param v0  value of the JRNMM
#'@param r value of the JRNMM
#'@return See description
#'@export
sigmoid_JRNMM <- function(x,vmax,v0,r){
  return(sigmoid_JRNMM_Cpp_(x,vmax,v0,r))
}


#'@rdname linear_JRNMM
#'@title linear_JRNMM
#'@description linear SDE: exact solution of 1st subequation of splitting procedure for the JRNMM
#'@param vec vector coming from the non-linear ODE from the Strang splitting scheme
#'@param dm exponential matrix coming from the linear SDE in the splitting scheme
#'@param xi gaussian increments
#'@return See description
#'@export
linear_JRNMM <- function(vec,dm,xi){
  return(linear_JRNMM_Cpp_(vec,dm,xi))
}


#'@rdname nonlinear_JRNMM
#'@title nonlinear_JRNMM
#'@description nonlinear ODE: exact solution of 2nd subequation of splitting procedure for the JRNMM
#'@param N number of populations of neurons
#'@param vec starting point before applying the solution to the non-linear ODE from the Strang splitting scheme
#'@param h step size used for path simulation
#'@param Theta  vector of continuous parameters
#'@param Rho matrix of {0,1} discrete parameters
#'@param K Matrix of strength parameters
#'@return See description
#'@export
nonlinear_JRNMM <- function(N,vec,h,Theta,Rho,K){
  return(nonlinear_JRNMM_Cpp_(N,vec,h,Theta,Rho,K))
}

#'@rdname KmatrixgivenLc
#'@title KmatrixgivenLc
#'@description Code to calculate the K matrix given L and c
#'@param N number of populations of neurons
#'@param L strength parameter
#'@param c [0,1] parameter
#'@return Kmatrix given L and c
#'@export
KmatrixgivenLc <- function(N,L,c){
  return(KmatrixgivenLc_(N,L,c))
}
