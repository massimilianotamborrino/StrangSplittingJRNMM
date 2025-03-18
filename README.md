# StrangSplittingJRNMM

R package for the Strang numerical splitting scheme of N populations of the Jansen-and-Rit neural mass models (JRNMM) proposed in Algorithm 2 in 
[1] S. Ditlevsen, M. Tamborrino, I. Tubikanec. Network Inference via Approximate Bayesian Computation. Illustration on a stochastic multi-population Neural Mass Model. Preprint at ArXiv: 2306.15787v2 https://arxiv.org/abs/2306.15787 

The R-package is written and maintained by Massimiliano Tamborrino (firstname dot secondname at warwick.ac.uk).

# What can you find in the package
In this package, we provide the codes for the simulation of trajectories from N populations of the JRNMM. Each population is a 6-dimensional SDE, so the process \(X=(Q,P)^T\) is 6N dimensional, with 3N-dimensional components \(Q=(Q^1,\ldots, Q^N)^T=(X_1^1,X_2^1,X_3^1,\ldots, X_1^N,X_2^N,X_3^N)^T)\) and \(P=(P^1,\ldots, P^N)^T=(X_4^1,X_5^1,X_6^1,\ldots, X_4^N,X_5^N,X_6^N)^T)\). 

The main routine is "fast_JRNMM_Splitting.R", which simulates trajectories of \(X(t), t\in[0,T]\) using the block structure of the model, yielding faster simulations compared to the alternative "JRNMM_Splitting" routine. The "fast_JRNMM_Splitting.R" routines takes as user input the number of populations "N", the time-discretisation grid "grid" (obtained by discretising \([0,T]\) by a time-step \(h\)), the time step "h", the starting vector "startv" representing the initial position of the process at time 0, the vector of diagonal entries of the matrix Gamma, "dGamma", the vector of diagonal entries of the matrix Sigma, "dSigma", the vector of continuous parameter "Theta", the matrix of {0,1} connection parameters "Rho" and the matrix of strength parameters (without diagonal) "K". We refer to Appendix A of [1] for more details on the Gamma and Sigma matrices. 

The "JRNMM_Splitting" routine has the same inputs as "fast_JRNMM_Splitting.R" with the only difference of passing the entire diagonal matrices Gamma and Sigma instead of only the diagonal entries. This routine is not taking advantage of the fact that these matrices are diagonal, leading to higher runtimes. This routine is left here for illustration purposes, but the use of "fast_JRNMM_Splitting.R" is recommended.

# How to install the package
* Tools/Install packages/ select the source folder
*To update The simplest way is to do it via devtools, using devtools::install_github("massimilianotamborrino/StrangSplittingJRNMM")

# Output of "fast_JRNMM_Splitting.R" and "JRNMM_Splitting.R"
Both routines return a 6NxM matrix, where the number of rows, 6N, corresponds to the 6N components of the N JRNMM populations, with M being the number of discrete-time points where the trajectories are evaluated, e.g., using equidistant points in \([0,T]\), we will have \(t_i=ih, i=0,\ldots, M, M=T/h+1\). Remember that the first 3N components (rows) corresponds to \(Q\) and the second 3N to \(P\). For inference purposes, the available observations are \(Y(t)=(Y^1(t),\ldots, Y^N(t))^T=(X_2^1(t)-X_3^1(t),\ldots,X_2^N(t)-X_3^N(t)\) which can be immediately determined from the first 3N components of the output \(X(t)\).
