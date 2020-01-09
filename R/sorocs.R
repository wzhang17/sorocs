#' A Bayesian nonparametric Dirichlet process mixtures to Eestimate Receiver Operating Characteristic (ROC) Surface model
#' @import MASS MCMCpack mvtnorm stats
#' @description A Bayesian nonparametric Dirichlet process mixtures to estimate the receiver operating characteristic (ROC) 
#' surfaces and the associated volume under the surface (VUS), a summary measure similar to the area under the 
#' curve measure for ROC curves. To model distributions flexibly, including their skewness and multi-modality
#' characteristics a Bayesian nonparametric Dirichlet process mixtures was used. Between-setting correlations is handled
#' by dependent Dirichlet process mixtures that have bivariate distributions with nonzero
#' correlations as their bases. To  accommodate ordering constraints, the stochastic ordering in the
#' context of mixture distributions was adopted.
#'
#' @param Yvariable1 Dependent variable at setting 1
#' @param Yvariable2 Dependent variable at setting 2
#' @param Xvariable1 independent variable at setting 1
#' @param Xvariable2 independent variable at setting 2
#' @param nsim Number of simulations
#' @param nburn Burn in number
#' @param gridY a regular sequence spanning the range of Y variable
#' @param gam0 Initial value for the test score distributions (e.g., a priori information between different disease populations for a single test or between multiple correlated tests)
#' @param gam1 Initial value for the test score distributions
#' @param lamb0 Initial value forthe test score distributions
#' @param lamb1 Initial value for the test score distributions
#' @param alpha1 fixed values of the precision parameters of the Dirichlet process
#' @param alpha2 fixed values of the precision parameters of the Dirichlet process
#' @param alpha3 fixed values of the precision parameters of the Dirichlet process
#' @param lambda1 fixed values of the precision parameters of the Dirichlet process
#' @param lambda2 fixed values of the precision parameters of the Dirichlet process
#' @param H trucation level number for Dirichlet process prior trucation approximation
#' @param L trucation level number for Dirichlet process prior trucation approximation 
#' @param mu1 fixed values of the bivariate normal parameters of the Dirichlet process
#' @param mu2 fixed values of the bivariate normal parameters of the Dirichlet process
#' @param mu3 fixed values of the bivariate normal parameters of the Dirichlet process
#' @param m1 fixed values of the bivariate normal parameters of the Dirichlet process
#' @param m2 fixed values of the bivariate normal parameters of the Dirichlet process
#' @param m3 fixed values of the bivariate normal parameters of the Dirichlet process
#' @param A1 Initial values of the bivariate normal parameters of the Dirichlet process
#' @param A2 Initial values of the bivariate normal parameters of the Dirichlet process
#' @param A3 Initial values of the bivariate normal parameters of the Dirichlet process
#' @param Sig1 Initial values of the inverse Wishart distribution parameters of the Dirichlet process
#' @param Sig2 Initial values of the inverse Wishart distribution parameters of the Dirichlet process
#' @param Sig3 Initial values of the inverse Wishart distribution parameters of the Dirichlet process

#' @param nu Initial values of the inverse Wishart distribution parameters of the Dirichlet process
#' @param C0 Initial values of the inverse Wishart distribution parameters of the Dirichlet process
#' @param a1 Initial shape values of the inverse-gamma base distributions for the Dirichlet process
#' @param a2 Initial shape values of the inverse-gamma base distributions for the Dirichlet process
#' @param b1 Initial scale values of the inverse-gamma base distributions for the Dirichlet process
#' @param b2 Initial scale values of the inverse-gamma base distributions for the Dirichlet process
#' @keywords Dirichlet process mixtures
#' @export
#' @return A list of posterior estimates
#' @examples
#' library(MASS)
#' library(MCMCpack)
#' library(mvtnorm)
#' data(asrm)
#' try1 <- sorocs:::sorocs(Yvariable1 =asrm$logREscoremean2, Yvariable2=asrm$logREscoremean1, 
#' gridY=seq(0,5,by=0.5), Xvariable1=asrm$TN12/asrm$JN12 , Xvariable2 =asrm$TNN12/asrm$JNN12)

sorocs <-function(## iterations
  nsim=4,
  nburn=2,
  Yvariable1 ,   # Setting 1
  Yvariable2  ,    # Setting 2
  gridY = seq(0,5,by=0.05),
  ## make X1 and X2
  Xvariable1,
  Xvariable2,
  gam0=-4.6 ,
  gam1=9.2 ,
  lamb0=-4.6 ,
  lamb1=9.2,
  ## initial values
  H=30,       
  L=30 ,        
  ## fixed values
  alpha1=1 ,
  alpha2=1 ,
  alpha3=1 ,
  
  lambda1=1 ,
  lambda2=1 ,
  
  mu1=matrix(c(0.5,0.5),2,1),
  mu2=matrix(c(1,1),2,1) ,
  mu3=matrix(c(3,3),2,1) ,
  
  m1=c(0,0) ,
  m2=c(0,0) ,
  m3=c(0,0) ,
  
  A1=10*diag(2) ,
  A2=10*diag(2) ,
  A3=10*diag(2) ,
  
  Sig1=matrix(c(1,0.5,0.5,1),2,2) ,
  Sig2=matrix(c(1,0.5,0.5,1),2,2) ,
  Sig3=matrix(c(1,0.5,0.5,1),2,2) ,
  
  nu=6 ,
  C0=10*diag(2) ,
  
  a1=2 ,
  a2=2 ,
  
  b1=0.1 ,
  b2=0.1 
  
){
  force(gridY)

  ngrid=length(gridY)
  iteration=1:((nsim-nburn)/10)
  N=length(Yvariable1)
#  gam0=-4.6
#  gam1=9.2

#  lamb0=-4.6
#  lamb1=9.2

  ## make initial diseased population based on IEs' diagnosis
  pb0=1/(1+exp(gam0+gam1*Xvariable1))
  pb1=1/(1+exp(lamb0+lamb1*Xvariable2))*exp(gam0+gam1*Xvariable1)/(1+exp(gam0+gam1*Xvariable1))
  pb2=exp(lamb0+lamb1*Xvariable2)/(1+exp(lamb0+lamb1*Xvariable2))*exp(gam0+gam1*Xvariable1)/(1+exp(gam0+gam1*Xvariable1))

  pb1[is.na(pb1)==T]=(1-pb0[is.na(pb1)==T])/2
  pb2[is.na(pb2)==T]=(1-pb0[is.na(pb2)==T])/2

  di=rep(0,N)
  for (i in 1:N){
    di[i]=sample(c(0,1,2),1,replace=TRUE,prob=c(pb0[i],pb1[i],pb2[i]))
  }


  ## initial values
#  H=30
#  L=30

  p_h1=rep(1,H)
  p_h2=rep(1,H)
  p_h3=rep(1,H)

  p_h1=p_h1/sum(p_h1)
  p_h2=p_h2/sum(p_h2)
  p_h3=p_h3/sum(p_h3)

  q_l1=rep(1,L)
  q_l2=rep(1,L)

  q_l1=q_l1/sum(q_l1)
  q_l2=q_l2/sum(q_l2)

  V_h1=rep(0.5,H)
  V_h2=rep(0.5,H)
  V_h3=rep(0.5,H)

  W_l1=rep(0.5,L)
  W_l2=rep(0.5,L)

  theta_i11=rep(1,N)
  theta_i12=rep(1,N)
  theta_i13=rep(1,N)
  theta_i21=rep(1,N)
  theta_i22=rep(1,N)
  theta_i23=rep(1,N)

  sigma_i1=rep(0.1,2*N)
  sigma_i2=rep(0.1,N)

  theta_star11=rep(1,H)
  theta_star12=rep(1,H)
  theta_star13=rep(1,H)
  theta_star21=rep(1,H)
  theta_star22=rep(1,H)
  theta_star23=rep(1,H)

  theta_star1=rbind(theta_star11,theta_star21)
  theta_star2=rbind(theta_star12,theta_star22)
  theta_star3=rbind(theta_star13,theta_star23)

  sigma_star1=rep(0.1,L)
  sigma_star2=rep(0.1,L)


  ## fixed values
#  alpha1=1
#  alpha2=1
#  alpha3=1

#  lambda1=1
#  lambda2=1

#  mu1=matrix(c(0.5,0.5),2,1)
#  mu2=matrix(c(1,1),2,1)
#  mu3=matrix(c(3,3),2,1)

#  m1=c(0,0)
#  m2=c(0,0)
#  m3=c(0,0)

#  A1=10*diag(2)
#  A2=10*diag(2)
#  A3=10*diag(2)

#  Sig1=matrix(c(1,0.5,0.5,1),2,2)
#  Sig2=matrix(c(1,0.5,0.5,1),2,2)
#  Sig3=matrix(c(1,0.5,0.5,1),2,2)

#  nu=6
#  C0=10*diag(2)

#  a1=2
#  a2=2

#  b1=0.1
#  b2=0.1


  ## define output files
  Y11_pred_pdf=matrix(1,nsim,ngrid)
  Y11_pred_cdf=matrix(1,nsim,ngrid)
  Y21_pred_pdf=matrix(1,nsim,ngrid)
  Y21_pred_cdf=matrix(1,nsim,ngrid)

  Y12_pred_pdf=matrix(1,nsim,ngrid)
  Y12_pred_cdf=matrix(1,nsim,ngrid)
  Y22_pred_pdf=matrix(1,nsim,ngrid)
  Y22_pred_cdf=matrix(1,nsim,ngrid)

  Y13_pred_pdf=matrix(1,nsim,ngrid)
  Y13_pred_cdf=matrix(1,nsim,ngrid)
  Y23_pred_pdf=matrix(1,nsim,ngrid)
  Y23_pred_cdf=matrix(1,nsim,ngrid)

  Y11_pdf=array(1,dim=c(H,L,ngrid))
  Y11_cdf=array(1,dim=c(H,L,ngrid))
  Y21_pdf=array(1,dim=c(H,L,L,ngrid))
  Y21_cdf=array(1,dim=c(H,L,L,ngrid))
  Y12_pdf=array(1,dim=c(H,H,L,ngrid))
  Y12_cdf=array(1,dim=c(H,H,L,ngrid))
  Y22_pdf=array(1,dim=c(H,H,L,L,ngrid))
  Y22_cdf=array(1,dim=c(H,H,L,L,ngrid))
  Y13_pdf=array(1,dim=c(H,H,H,L,ngrid))
  Y13_cdf=array(1,dim=c(H,H,H,L,ngrid))
  Y23_pdf=array(1,dim=c(H,H,L,L,ngrid))
  Y23_cdf=array(1,dim=c(H,H,L,L,ngrid))
  Y23_pdf2=matrix(1,H,ngrid)
  Y23_cdf2=matrix(1,H,ngrid)

  VUS1_temp=rep(1,nsim)
  VUS2_temp=rep(1,nsim)

  VUS1_out=rep(1,(nsim-nburn))
  VUS2_out=rep(1,(nsim-nburn))

  ## extra preallocation
  iset=c(1:N)

  ss1=matrix(c(1,0.5,0.5,1),2,2)
  ss2=matrix(c(1,0.5,0.5,1),2,2)
  ss3=matrix(c(1,0.5,0.5,1),2,2)
  ss4=0.1
  ss5=0.1

  p1=rep(1,H)
  p2=rep(1,H)
  p3=rep(1,H)
  p4=rep(1,L)
  p5=rep(1,L)

  Sig11_T=rep(1,(nsim-nburn)/10)
  Sig12_T=rep(1,(nsim-nburn)/10)
  Sig13_T=rep(1,(nsim-nburn)/10)
  Sig21_T=rep(1,(nsim-nburn)/10)
  Sig22_T=rep(1,(nsim-nburn)/10)
  Sig23_T=rep(1,(nsim-nburn)/10)
  Sig31_T=rep(1,(nsim-nburn)/10)
  Sig32_T=rep(1,(nsim-nburn)/10)
  Sig33_T=rep(1,(nsim-nburn)/10)

  theta_star11_T=rep(1,(nsim-nburn)/10)
  theta_star12_T=rep(1,(nsim-nburn)/10)
  theta_star13_T=rep(1,(nsim-nburn)/10)
  theta_star21_T=rep(1,(nsim-nburn)/10)
  theta_star22_T=rep(1,(nsim-nburn)/10)
  theta_star23_T=rep(1,(nsim-nburn)/10)

  sigma_star1_T=rep(1,(nsim-nburn)/10)
  sigma_star2_T=rep(1,(nsim-nburn)/10)

  theta_star11_out=rbind(theta_star11,matrix(1,nsim,H))
  theta_star12_out=rbind(theta_star12,matrix(1,nsim,H))
  theta_star13_out=rbind(theta_star13,matrix(1,nsim,H))
  theta_star21_out=rbind(theta_star21,matrix(1,nsim,H))
  theta_star22_out=rbind(theta_star22,matrix(1,nsim,H))
  theta_star23_out=rbind(theta_star23,matrix(1,nsim,H))

  sigma_star1_out=rbind(sigma_star1,matrix(1,nsim,L))
  sigma_star2_out=rbind(sigma_star2,matrix(1,nsim,L))

  V_h1_T=rep(1,(nsim-nburn)/10)
  V_h2_T=rep(1,(nsim-nburn)/10)
  V_h3_T=rep(1,(nsim-nburn)/10)

  di_out=rbind(di,matrix(1,nsim,N))

  corr1_out=rep(0.5,nsim)
  corr2_out=rep(0.5,nsim)
  corr3_out=rep(0.5,nsim)

  ## make Metropolis-Hastings function for theta_star11 and theta_star21
  MH_th1=function(theta1,theta2,h){
    ist1=iset[S_i1==h & di==0]
    ist2=iset[S_i1==h & di==1]
    ist3=iset[S_i1==h & di==2]

    part0=dmvnorm(c(theta1,theta2),mu1,Sig1)

    part1=prod(dnorm(Yvariable1[ist1],theta1,sqrt(sigma_i1[ist1])))
    part2=prod(dnorm(Yvariable2[ist1],theta2,sqrt(pmax(sigma_i1[(N+ist1)],sigma_i2[ist1]))))

    part3=prod(dnorm(Yvariable1[ist2],pmax(theta1,theta_i12[ist2]),sqrt(sigma_i1[ist2])))
    part4=prod(dnorm(Yvariable2[ist2],pmax(theta2,theta_i22[ist2]),sqrt(pmax(sigma_i1[(N+ist2)],sigma_i2[ist2]))))

    part5=prod(dnorm(Yvariable1[ist3],pmax(theta1,theta_i12[ist3],theta_i13[ist3]),sqrt(sigma_i1[ist3])))
    part6=prod(dnorm(Yvariable2[ist3],pmax(theta2,theta_i22[ist3],theta_i23[ist3]),sqrt(pmax(sigma_i1[(N+ist3)],sigma_i2[ist3]))))

    part1[is.na(part1)==TRUE]=1
    part2[is.na(part2)==TRUE]=1
    part3[is.na(part3)==TRUE]=1
    part4[is.na(part4)==TRUE]=1
    part5[is.na(part5)==TRUE]=1
    part6[is.na(part6)==TRUE]=1

    fn=part0*part1*part2*part3*part4*part5*part6

    return(fn)
  }

  ## make Metropolis-Hastings function for theta_star12 and theta_star22
  MH_th2=function(theta1,theta2,h){
    ist2=iset[S_i2==h & di==1]
    ist3=iset[S_i2==h & di==2]

    part0=dmvnorm(c(theta1,theta2),mu2,Sig2)

    part1=prod(dnorm(Yvariable1[ist2],pmax(theta_i11[ist2],theta1),sqrt(sigma_i1[ist2])))
    part2=prod(dnorm(Yvariable2[ist2],pmax(theta_i21[ist2],theta2),sqrt(pmax(sigma_i1[(N+ist2)],sigma_i2[ist2]))))

    part3=prod(dnorm(Yvariable1[ist3],pmax(theta_i11[ist3],theta1,theta_i13[ist3]),sqrt(sigma_i1[ist3])))
    part4=prod(dnorm(Yvariable2[ist3],pmax(theta_i21[ist3],theta2,theta_i23[ist3]),sqrt(pmax(sigma_i1[(N+ist3)],sigma_i2[ist3]))))

    part1[is.na(part1)==TRUE]=1
    part2[is.na(part2)==TRUE]=1
    part3[is.na(part3)==TRUE]=1
    part4[is.na(part4)==TRUE]=1

    fn=part0*part1*part2*part3*part4

    return(fn)
  }

  ## make Metropolis-Hastings function for theta_star13 and theta_star23
  MH_th3=function(theta1,theta2,h){
    ist3=iset[S_i3==h & di==2]

    part0=dmvnorm(c(theta1,theta2),mu3,Sig3)

    part1=prod(dnorm(Yvariable1[ist3],pmax(theta_i11[ist3],theta_i12[ist3],theta1),sqrt(sigma_i1[ist3])))
    part2=prod(dnorm(Yvariable2[ist3],pmax(theta_i21[ist3],theta_i22[ist3],theta2),sqrt(pmax(sigma_i1[(N+ist3)],sigma_i2[ist3]))))

    part1[is.na(part1)==TRUE]=1
    part2[is.na(part2)==TRUE]=1

    fn=part0*part1*part2

    return(fn)
  }

  ## make Metropolis-Hastings function for sigma_star1
  MH_sig1=function(sigma,l){
    ist1=iset[T_i1[1:N]==l & di==0]
    ist2=iset[T_i1[(N+1):(2*N)]==l & di==0]
    ist3=iset[T_i1[1:N]==l & di==1]
    ist4=iset[T_i1[(N+1):(2*N)]==l & di==1]
    ist5=iset[T_i1[1:N]==l & di==2]
    ist6=iset[T_i1[(N+1):(2*N)]==l & di==2]

    part0=dinvgamma(sigma,shape=a1,scale=b1)

    part1=prod(dnorm(Yvariable1[ist1],theta_i11[ist1],sqrt(sigma)))
    part2=prod(dnorm(Yvariable2[ist2],theta_i21[ist2],sqrt(pmax(sigma,sigma_i2[ist2]))))

    part3=prod(dnorm(Yvariable1[ist3],pmax(theta_i11[ist3],theta_i12[ist3]),sqrt(sigma)))
    part4=prod(dnorm(Yvariable2[ist4],pmax(theta_i21[ist4],theta_i22[ist4]),sqrt(pmax(sigma,sigma_i2[ist4]))))

    part5=prod(dnorm(Yvariable1[ist5],pmax(theta_i11[ist5],theta_i12[ist5],theta_i13[ist5]),sqrt(sigma)))
    part6=prod(dnorm(Yvariable2[ist6],pmax(theta_i21[ist6],theta_i22[ist6],theta_i23[ist6]),sqrt(pmax(sigma,sigma_i2[ist6]))))

    part1[is.na(part1)==TRUE]=1
    part2[is.na(part2)==TRUE]=1
    part3[is.na(part3)==TRUE]=1
    part4[is.na(part4)==TRUE]=1
    part5[is.na(part5)==TRUE]=1
    part6[is.na(part6)==TRUE]=1

    fn=part0*part1*part2*part3*part4*part5*part6

    return(fn)
  }

  ## make Metropolis-Hastings function for sigma_star2
  MH_sig2=function(sigma,l){
    ist1=iset[T_i2==l & di==0]
    ist2=iset[T_i2==l & di==1]
    ist3=iset[T_i2==l & di==2]

    part0=dinvgamma(sigma,shape=a2,scale=b2)

    part1=prod(dnorm(Yvariable2[ist1],theta_i21[ist1],sqrt(pmax(sigma_i1[(N+ist1)],sigma))))
    part2=prod(dnorm(Yvariable2[ist2],pmax(theta_i21[ist2],theta_i22[ist2]),sqrt(pmax(sigma_i1[(N+ist2)],sigma))))
    part3=prod(dnorm(Yvariable2[ist3],pmax(theta_i21[ist3],theta_i22[ist3],theta_i23[ist3]),sqrt(pmax(sigma_i1[(N+ist3)],sigma))))

    part1[is.na(part1)==TRUE]=1
    part2[is.na(part2)==TRUE]=1
    part3[is.na(part3)==TRUE]=1

    fn=part0*part1*part2*part3

    return(fn)
  }

  ## For easier calculation for predictive densities

  dnormfun <<- function(muvalue, sigmavalue, aNumber = gridY[i] ){
    force(aNumber)
    ansvalue=dnorm(aNumber,muvalue,sqrt(sigmavalue))
    return(ansvalue)
  }


  pnormfun <<-  function(muvalue, sigmavalue, aNumber = gridY[i] ){
    force(aNumber)
    returnans=pnorm(aNumber,muvalue,sqrt(sigmavalue))
    return(returnans)
  }



  #####################
  ###### MCMC
  #####################
  ptm=proc.time()
  for (gt in 1:nsim){

    ## update indicator S_i1 of theta_i11 and theta_i21
    S_i1=rep(0,N)
    for (i in iset[di==0]){                                        # for disease=1
      pr1=dnorm(Yvariable1[i],theta_star11,sqrt(sigma_i1[i]))
      pr2=dnorm(Yvariable2[i],theta_star21,sqrt(pmax(sigma_i1[N+i],sigma_i2[i])))
      pr1[is.na(pr1)==TRUE]=1
      pr2[is.na(pr2)==TRUE]=1

      prob1=p_h1*pr1*pr2
      prob1=prob1/sum(prob1)
      S_i1[i]=sample(1:H,1,replace=TRUE,prob=prob1)
    }

    for (i in iset[di==1]){                                        # for disease=2
      pr1=dnorm(Yvariable1[i],pmax(theta_star11,theta_i12[i]),sqrt(sigma_i1[i]))
      pr2=dnorm(Yvariable2[i],pmax(theta_star21,theta_i22[i]),sqrt(pmax(sigma_i1[N+i],sigma_i2[i])))
      pr1[is.na(pr1)==TRUE]=1
      pr2[is.na(pr2)==TRUE]=1

      prob1=p_h1*pr1*pr2
      prob1=prob1/sum(prob1)
      S_i1[i]=sample(1:H,1,replace=TRUE,prob=prob1)
    }

    for (i in iset[di==2]){                                        # for disease=3
      pr1=dnorm(Yvariable1[i],pmax(theta_star11,theta_i12[i],theta_i13[i]),sqrt(sigma_i1[i]))
      pr2=dnorm(Yvariable2[i],pmax(theta_star21,theta_i22[i],theta_i23[i]),sqrt(pmax(sigma_i1[N+i],sigma_i2[i])))
      pr1[is.na(pr1)==TRUE]=1
      pr2[is.na(pr2)==TRUE]=1

      prob1=p_h1*pr1*pr2
      prob1=prob1/sum(prob1)
      S_i1[i]=sample(1:H,1,replace=TRUE,prob=prob1)
    }

    ## update indicator S_i2 of theta_i12 and theta_i22
    S_i2=rep(0,N)
    for (i in iset[di==1]){                                        # for disease=2
      pr1=dnorm(Yvariable1[i],pmax(theta_i11[i],theta_star12),sqrt(sigma_i1[i]))
      pr2=dnorm(Yvariable2[i],pmax(theta_i21[i],theta_star22),sqrt(pmax(sigma_i1[N+i],sigma_i2[i])))
      pr1[is.na(pr1)==TRUE]=1
      pr2[is.na(pr2)==TRUE]=1

      prob2=p_h2*pr1*pr2
      prob2=prob2/sum(prob2)
      S_i2[i]=sample(1:H,1,replace=TRUE,prob=prob2)
    }

    for (i in iset[di==2]){                                        # for disease=3
      pr1=dnorm(Yvariable1[i],pmax(theta_i11[i],theta_star12,theta_i13[i]),sqrt(sigma_i1[i]))
      pr2=dnorm(Yvariable2[i],pmax(theta_i21[i],theta_star22,theta_i23[i]),sqrt(pmax(sigma_i1[N+i],sigma_i2[i])))
      pr1[is.na(pr1)==TRUE]=1
      pr2[is.na(pr2)==TRUE]=1

      prob2=p_h2*pr1*pr2
      prob2=prob2/sum(prob2)
      S_i2[i]=sample(1:H,1,replace=TRUE,prob=prob2)
    }

    ## update indicator S_i3 of theta_i13 and theta_i23
    S_i3=rep(0,N)
    for (i in iset[di==2]){                                        # for disease=3
      pr1=dnorm(Yvariable1[i],pmax(theta_i11[i],theta_i12[i],theta_star13),sqrt(sigma_i1[i]))
      pr2=dnorm(Yvariable2[i],pmax(theta_i21[i],theta_i22[i],theta_star23),sqrt(pmax(sigma_i1[N+i],sigma_i2[i])))
      pr1[is.na(pr1)==TRUE]=1
      pr2[is.na(pr2)==TRUE]=1

      prob3=p_h3*pr1*pr2
      prob3=prob3/sum(prob3)
      S_i3[i]=sample(1:H,1,replace=TRUE,prob=prob3)
    }

    ## update indicator T_i1 of sigma_i1
    T_i1=rep(0,2*N)
    for (i in iset[di==0 & is.na(Yvariable1)==FALSE]){                                        # for disease=1 at setting 1
      prob4=q_l1*dnorm(Yvariable1[i],theta_i11[i],sqrt(sigma_star1))
      prob4=prob4/sum(prob4)
      T_i1[i]=sample(1:L,1,replace=TRUE,prob=prob4)
    }
    for (i in iset[di==0 & is.na(Yvariable2)==FALSE]){                                        # for disease=1 at setting 2
      prob5=q_l1*dnorm(Yvariable2[i],theta_i21[i],sqrt(pmax(sigma_star1,sigma_i2[i])))
      prob5=prob5/sum(prob5)
      T_i1[(N+i)]=sample(1:L,1,replace=TRUE,prob=prob5)
    }

    for (i in iset[di==1 & is.na(Yvariable1)==FALSE]){                                        # for disease=2 at setting 1
      prob6=q_l1*dnorm(Yvariable1[i],pmax(theta_i11[i],theta_i12[i]),sqrt(sigma_star1))
      prob6=prob6/sum(prob6)
      T_i1[i]=sample(1:L,1,replace=TRUE,prob=prob6)
    }
    for (i in iset[di==1 & is.na(Yvariable2)==FALSE]){                                        # for disease=2 at setting 2
      prob7=q_l1*dnorm(Yvariable2[i],pmax(theta_i21[i],theta_i22[i]),sqrt(pmax(sigma_star1,sigma_i2[i])))
      prob7=prob7/sum(prob7)
      T_i1[(N+i)]=sample(1:L,1,replace=TRUE,prob=prob7)
    }

    for (i in iset[di==2 & is.na(Yvariable1)==FALSE]){                                        # for disease=3 at setting 1
      prob8=q_l1*dnorm(Yvariable1[i],pmax(theta_i11[i],theta_i12[i],theta_i13[i]),sqrt(sigma_star1))
      prob8=prob8/sum(prob8)
      T_i1[i]=sample(1:L,1,replace=TRUE,prob=prob8)
    }
    for (i in iset[di==2 & is.na(Yvariable2)==FALSE]){                                        # for disease=3 at setting 2
      prob9=q_l1*dnorm(Yvariable2[i],pmax(theta_i21[i],theta_i22[i],theta_i23[i]),sqrt(pmax(sigma_star1,sigma_i2[i])))
      prob9=prob9/sum(prob9)
      T_i1[(N+i)]=sample(1:L,1,replace=TRUE,prob=prob9)
    }

    ## update indicator T_i2 of sigma_i2
    T_i2=rep(0,N)
    for (i in iset[di==0 & is.na(Yvariable2)==FALSE]){                                        # for disease=1 at setting 2
      prob10=q_l2*dnorm(Yvariable2[i],theta_i21[i],sqrt(pmax(sigma_i1[(N+i)],sigma_star2)))
      prob10=prob10/sum(prob10)
      T_i2[i]=sample(1:L,1,replace=TRUE,prob=prob10)
    }

    for (i in iset[di==1 & is.na(Yvariable2)==FALSE]){                                        # for disease=2 at setting 2
      prob11=q_l2*dnorm(Yvariable2[i],pmax(theta_i21[i],theta_i22[i]),sqrt(pmax(sigma_i1[(N+i)],sigma_star2)))
      prob11=prob11/sum(prob11)
      T_i2[i]=sample(1:L,1,replace=TRUE,prob=prob11)
    }

    for (i in iset[di==2 & is.na(Yvariable2)==FALSE]){                                        # for disease=3 at setting 2
      prob12=q_l2*dnorm(Yvariable2[i],pmax(theta_i21[i],theta_i22[i],theta_i23[i]),sqrt(pmax(sigma_i1[(N+i)],sigma_star2)))
      prob12=prob12/sum(prob12)
      T_i2[i]=sample(1:L,1,replace=TRUE,prob=prob12)
    }



    ## update theta_star11 and theta_star21
    for (h in 1:H){
      if (gt!=1 & gt%%100==1){
        p1[h]=sum(diff(theta_star11_out[(gt-100):(gt-1),h])!=0)/100
        if (p1[h]<0.34){ss1=ss1/sqrt(2)}
        else if (p1[h]>0.54){ss1=ss1*sqrt(2)}
        else {ss1=ss1}
      }
      else {ss1=ss1}

      uu1=runif(1,0,1)
      th1=rmvnorm(1,c(theta_star11[h],theta_star21[h]),ss1)
      RR1=MH_th1(th1[1],th1[2],h)/MH_th1(theta_star11[h],theta_star21[h],h)
      if (uu1<=RR1&is.na(RR1)=="FALSE"){theta_star11[h]=th1[1]}&{theta_star21[h]=th1[2]}
      else {theta_star11[h]=theta_star11[h]}&{theta_star21[h]=theta_star21[h]}
    }


    ## update theta_star12 and theta_star22
    for (h in 1:H){
      if (gt!=1 & gt%%100==1){
        p2[h]=sum(diff(theta_star12_out[(gt-100):(gt-1),h])!=0)/100
        if (p2[h]<0.34){ss2=ss2/sqrt(2)}
        else if (p2[h]>0.54){ss2=ss2*sqrt(2)}
        else {ss2=ss2}
      }
      else {ss2=ss2}

      uu2=runif(1,0,1)
      th2=rmvnorm(1,c(theta_star12[h],theta_star22[h]),ss2)
      RR2=MH_th2(th2[1],th2[2],h)/MH_th2(theta_star12[h],theta_star22[h],h)
      if (uu2<=RR2&is.na(RR2)=="FALSE"){theta_star12[h]=th2[1]}&{theta_star22[h]=th2[2]}
      else {theta_star12[h]=theta_star12[h]}&{theta_star22[h]=theta_star22[h]}
    }


    ## update theta_star13 and theta_star23
    for (h in 1:H){
      if (gt!=1 & gt%%100==1){
        p3[h]=sum(diff(theta_star13_out[(gt-100):(gt-1),h])!=0)/100
        if (p3[h]<0.34){ss3=ss3/sqrt(2)}
        else if (p3[h]>0.54){ss3=ss3*sqrt(2)}
        else {ss3=ss3}
      }
      else {ss3=ss3}

      uu3=runif(1,0,1)
      th3=rmvnorm(1,c(theta_star13[h],theta_star23[h]),ss3)
      RR3=MH_th3(th3[1],th3[2],h)/MH_th3(theta_star13[h],theta_star23[h],h)
      if (uu3<=RR3&is.na(RR3)=="FALSE"){theta_star13[h]=th3[1]}&{theta_star23[h]=th3[2]}
      else {theta_star13[h]=theta_star13[h]}&{theta_star23[h]=theta_star23[h]}
    }

    theta_star1=rbind(theta_star11,theta_star21)
    theta_star2=rbind(theta_star12,theta_star22)
    theta_star3=rbind(theta_star13,theta_star23)

    ## update sigma_star1
    for (l in 1:L){
      if (gt!=1 & gt%%100==1){
        p4[l]=sum(diff(sigma_star1_out[(gt-100):(gt-1),l])!=0)/100
        if (p4[l]<0.34){ss4=ss4/sqrt(2)}
        else if (p4[l]>0.54){ss4=ss4*sqrt(2)}
        else {ss4=ss4}
      }
      else {ss4=ss4}

      uu4=runif(1,0,1)
      sig1=rnorm(1,sigma_star1[l],ss4)
      if (sig1>0){
        RR4=MH_sig1(sig1,l)/MH_sig1(sigma_star1[l],l)
        if (uu4<=RR4&is.na(RR4)=="FALSE"){sigma_star1[l]=sig1}
        else {sigma_star1[l]=sigma_star1[l]}
      }
      else {sigma_star1[l]=sigma_star1[l]}
    }


    ## update sigma_star2
    for (l in 1:L){
      if (gt!=1 & gt%%100==1){
        p5[l]=sum(diff(sigma_star2_out[(gt-100):(gt-1),l])!=0)/100
        if (p5[l]<0.34){ss5=ss5/sqrt(2)}
        else if (p5[l]>0.54){ss5=ss5*sqrt(2)}
        else {ss5=ss5}
      }
      else {ss5=ss5}

      uu5=runif(1,0,1)
      sig2=rnorm(1,sigma_star2[l],ss5)
      if (sig2>0){
        RR5=MH_sig2(sig2,l)/MH_sig2(sigma_star2[l],l)
        if (uu5<=RR5&is.na(RR5)=="FALSE"){sigma_star2[l]=sig2}
        else {sigma_star2[l]=sigma_star2[l]}
      }
      else {sigma_star2[l]=sigma_star2[l]}
    }



    ## update theta_i and sigma_i
    for (i in 1:N){
      theta_i11[i]=theta_star11[S_i1[i]]
      theta_i21[i]=theta_star21[S_i1[i]]

      if (S_i2[i]==0){theta_i12[i]=theta_i12[i]}&{theta_i22[i]=theta_i22[i]}
      else {theta_i12[i]=theta_star12[S_i2[i]]}&{theta_i22[i]=theta_star22[S_i2[i]]}

      if (S_i3[i]==0){theta_i13[i]=theta_i13[i]}&{theta_i23[i]=theta_i23[i]}
      else {theta_i13[i]=theta_star13[S_i3[i]]}&{theta_i23[i]=theta_star23[S_i3[i]]}

      if (T_i2[i]==0){sigma_i2[i]=sigma_i2[i]} else {sigma_i2[i]=sigma_star2[T_i2[i]]}
    }

    for (i in 1:(2*N)){
      if (T_i1[i]==0){sigma_i1[i]=sigma_i1[i]} else {sigma_i1[i]=sigma_star1[T_i1[i]]}
    }


    ## update mu
    mu1=t(rmvnorm(1,solve(solve(A1)+H*solve(Sig1))%*%(solve(A1)%*%m1+solve(Sig1)%*%rowSums(theta_star1)),solve(solve(A1)+H*solve(Sig1))))
    mu2=t(rmvnorm(1,solve(solve(A2)+H*solve(Sig2))%*%(solve(A2)%*%m2+solve(Sig2)%*%rowSums(theta_star2)),solve(solve(A2)+H*solve(Sig2))))
    mu3=t(rmvnorm(1,solve(solve(A3)+H*solve(Sig3))%*%(solve(A3)%*%m3+solve(Sig3)%*%rowSums(theta_star3)),solve(solve(A3)+H*solve(Sig3))))


    # update Sig
    Sig1=riwish(nu+H,C0+(theta_star1-mu1%*%rep(1,H))%*%t((theta_star1-mu1%*%rep(1,H))))
    Sig2=riwish(nu+H,C0+(theta_star2-mu2%*%rep(1,H))%*%t((theta_star2-mu2%*%rep(1,H))))
    Sig3=riwish(nu+H,C0+(theta_star3-mu3%*%rep(1,H))%*%t((theta_star3-mu3%*%rep(1,H))))



    ## update V_h1, V_h2, V_h3
    for (h in 1:H){
      n_st1=sum(S_i1==h)
      n_st2=sum(S_i1>h)
      V_h1[h]=rbeta(1,1+n_st1,alpha1+n_st2)

      n_st3=sum(S_i2==h)
      n_st4=sum(S_i2>h)
      V_h2[h]=rbeta(1,1+n_st3,alpha2+n_st4)

      n_st5=sum(S_i3==h)
      n_st6=sum(S_i3>h)
      V_h3[h]=rbeta(1,1+n_st5,alpha3+n_st6)
    }


    ## update W_l1, W_l2
    for (l in 1:L){
      n_st7=sum(T_i1==l)
      n_st8=sum(T_i1>l)
      W_l1[l]=rbeta(1,1+n_st7,lambda1+n_st8)

      n_st9=sum(T_i2==l)
      n_st10=sum(T_i2>l)
      W_l2[l]=rbeta(1,1+n_st9,lambda2+n_st10)
    }


    ## update p_h1, p_h2, p_h3
    cp1=1
    cp2=1
    cp3=1
    for (h in 1:H){
      p_h1[h]=V_h1[h]*cp1
      cp1=cp1*(1-V_h1[h])

      p_h2[h]=V_h2[h]*cp2
      cp2=cp2*(1-V_h2[h])

      p_h3[h]=V_h3[h]*cp3
      cp3=cp3*(1-V_h3[h])
    }
    p_h1=p_h1/sum(p_h1)
    p_h2=p_h2/sum(p_h2)
    p_h3=p_h3/sum(p_h3)


    ## update q_l1, q_l2
    cp4=1
    cp5=1
    for (l in 1:L){
      q_l1[l]=W_l1[l]*cp4
      cp4=cp4*(1-W_l1[l])

      q_l2[l]=W_l2[l]*cp5
      cp5=cp5*(1-W_l2[l])
    }
    q_l1=q_l1/sum(q_l1)
    q_l2=q_l2/sum(q_l2)


    ## update di
    dis01=dnorm(Yvariable1,theta_i11,sqrt(sigma_i1[1:N]))
    dis02=dnorm(Yvariable2,theta_i21,sqrt(pmax(sigma_i1[(N+1):(2*N)],sigma_i2)))

    dis11=dnorm(Yvariable1,pmax(theta_i11,theta_i12),sqrt(sigma_i1[1:N]))
    dis12=dnorm(Yvariable2,pmax(theta_i21,theta_i22),sqrt(pmax(sigma_i1[(N+1):(2*N)],sigma_i2)))

    dis21=dnorm(Yvariable1,pmax(theta_i11,theta_i12,theta_i13),sqrt(sigma_i1[1:N]))
    dis22=dnorm(Yvariable2,pmax(theta_i21,theta_i22,theta_i23),sqrt(pmax(sigma_i1[(N+1):(2*N)],sigma_i2)))

    dis01[is.na(dis01)==TRUE]=1
    dis02[is.na(dis02)==TRUE]=1
    dis11[is.na(dis11)==TRUE]=1
    dis12[is.na(dis12)==TRUE]=1
    dis21[is.na(dis21)==TRUE]=1
    dis22[is.na(dis22)==TRUE]=1

    prd0=dis01*dis02*(1/(1+exp(gam0+gam1*Xvariable1)))
    prd1=dis11*dis12*(1/(1+exp(lamb0+lamb1*Xvariable2))*exp(gam0+gam1*Xvariable1)/(1+exp(gam0+gam1*Xvariable1)))
    prd2=dis21*dis22*(exp(lamb0+lamb1*Xvariable2)/(1+exp(lamb0+lamb1*Xvariable2))*exp(gam0+gam1*Xvariable1)/(1+exp(gam0+gam1*Xvariable1)))

    prd1[is.na(prd1)==TRUE]=0
    prd2[is.na(prd2)==TRUE]=0

    pd0=prd0/(prd0+prd1+prd2)
    pd1=prd1/(prd0+prd1+prd2)
    pd2=prd2/(prd0+prd1+prd2)

    for (i in 1:N){
      di[i]=sample(c(0,1,2),1,replace=TRUE,prob=c(pd0[i],pd1[i],pd2[i]))
    }


    ## predictive pdf and cdf over Y11
    #rm(dnormfun)
    for (i in 1:ngrid){
      Y11_pdf[,,i]=kronecker(p_h1,q_l1)*kronecker(theta_star11,sigma_star1,FUN="dnormfun")
      Y11_cdf[,,i]=kronecker(p_h1,q_l1)*kronecker(theta_star11,sigma_star1,FUN="pnormfun")
    }
    Y11_pred_pdf[gt,]=colSums(Y11_pdf,dims=2)
    Y11_pred_cdf[gt,]=colSums(Y11_cdf,dims=2)


    ## predictive pdf and cdf over Y21
    for (i in 1:ngrid){
      Y21_pdf[,,,i]=kronecker(p_h1,kronecker(q_l1,q_l2))*kronecker(theta_star21,kronecker(sigma_star1,sigma_star2,FUN="pmax"),FUN="dnormfun")
      Y21_cdf[,,,i]=kronecker(p_h1,kronecker(q_l1,q_l2))*kronecker(theta_star21,kronecker(sigma_star1,sigma_star2,FUN="pmax"),FUN="pnormfun")
    }
    Y21_pred_pdf[gt,]=colSums(Y21_pdf,dims=3)
    Y21_pred_cdf[gt,]=colSums(Y21_cdf,dims=3)


    ## predictive pdf and cdf over Y12
    for (i in 1:ngrid){
      Y12_pdf[,,,i]=kronecker(p_h1,kronecker(p_h2,q_l1))*kronecker(kronecker(theta_star11,theta_star12,FUN="pmax"),sigma_star1,FUN="dnormfun")
      Y12_cdf[,,,i]=kronecker(p_h1,kronecker(p_h2,q_l1))*kronecker(kronecker(theta_star11,theta_star12,FUN="pmax"),sigma_star1,FUN="pnormfun")
    }
    Y12_pred_pdf[gt,]=colSums(Y12_pdf,dims=3)
    Y12_pred_cdf[gt,]=colSums(Y12_cdf,dims=3)


    ## predictive pdf and cdf over Y22
    for (i in 1:ngrid){
      Y22_pdf[,,,,i]=kronecker(p_h1,kronecker(p_h2,kronecker(q_l1,q_l2)))*kronecker(kronecker(theta_star21,theta_star22,FUN="pmax"),kronecker(sigma_star1,sigma_star2,FUN="pmax"),FUN="dnormfun")
      Y22_cdf[,,,,i]=kronecker(p_h1,kronecker(p_h2,kronecker(q_l1,q_l2)))*kronecker(kronecker(theta_star21,theta_star22,FUN="pmax"),kronecker(sigma_star1,sigma_star2,FUN="pmax"),FUN="pnormfun")
    }
    Y22_pred_pdf[gt,]=colSums(Y22_pdf,dims=4)
    Y22_pred_cdf[gt,]=colSums(Y22_cdf,dims=4)


    ## predictive pdf and cdf over Y13
    for (i in 1:ngrid){
      Y13_pdf[,,,,i]=kronecker(p_h1,kronecker(p_h2,kronecker(p_h3,q_l1)))*kronecker(kronecker(theta_star11,kronecker(theta_star12,theta_star13,FUN="pmax"),FUN="pmax"),sigma_star1,FUN="dnormfun")
      Y13_cdf[,,,,i]=kronecker(p_h1,kronecker(p_h2,kronecker(p_h3,q_l1)))*kronecker(kronecker(theta_star11,kronecker(theta_star12,theta_star13,FUN="pmax"),FUN="pmax"),sigma_star1,FUN="pnormfun")
    }
    Y13_pred_pdf[gt,]=colSums(Y13_pdf,dims=4)
    Y13_pred_cdf[gt,]=colSums(Y13_cdf,dims=4)


    ## predictive pdf and cdf over Y23
    for (h1 in 1:H){
      for (i in 1:ngrid){
        Y23_pdf[,,,,i]=kronecker(p_h1[h1],kronecker(p_h2,kronecker(p_h3,kronecker(q_l1,q_l2))))*kronecker(kronecker(theta_star21[h1],kronecker(theta_star22,theta_star23,FUN="pmax"),FUN="pmax"),kronecker(sigma_star1,sigma_star2,FUN="pmax"),FUN="dnormfun")
        Y23_cdf[,,,,i]=kronecker(p_h1[h1],kronecker(p_h2,kronecker(p_h3,kronecker(q_l1,q_l2))))*kronecker(kronecker(theta_star21[h1],kronecker(theta_star22,theta_star23,FUN="pmax"),FUN="pmax"),kronecker(sigma_star1,sigma_star2,FUN="pmax"),FUN="pnormfun")
      }
      Y23_pdf2[h1,]=colSums(Y23_pdf,dims=4)
      Y23_cdf2[h1,]=colSums(Y23_cdf,dims=4)
    }
    Y23_pred_pdf[gt,]=colSums(Y23_pdf2)
    Y23_pred_cdf[gt,]=colSums(Y23_cdf2)



    ## calculate VUS
    PrY11=sample(gridY,200,replace=TRUE,prob=Y11_pred_pdf[gt,])   #randomly selected person for d=1 at setting 1
    PrY12=sample(gridY,200,replace=TRUE,prob=Y12_pred_pdf[gt,])   #randomly selected person for d=2 at setting 1
    PrY13=sample(gridY,200,replace=TRUE,prob=Y13_pred_pdf[gt,])   #randomly selected person for d=3 at setting 1

    PrY21=sample(gridY,200,replace=TRUE,prob=Y21_pred_pdf[gt,])   #randomly selected person for d=1 at setting 2
    PrY22=sample(gridY,200,replace=TRUE,prob=Y22_pred_pdf[gt,])   #randomly selected person for d=2 at setting 2
    PrY23=sample(gridY,200,replace=TRUE,prob=Y23_pred_pdf[gt,])   #randomly selected person for d=3 at setting 2

    VUS1_temp[gt]=sum(PrY11<PrY12 & PrY12<PrY13)/200
    VUS2_temp[gt]=sum(PrY21<PrY22 & PrY22<PrY23)/200

    ## make output files
    theta_star11_out[gt+1,]=theta_star11
    theta_star12_out[gt+1,]=theta_star12
    theta_star13_out[gt+1,]=theta_star13
    theta_star21_out[gt+1,]=theta_star21
    theta_star22_out[gt+1,]=theta_star22
    theta_star23_out[gt+1,]=theta_star23

    sigma_star1_out[gt+1,]=sigma_star1
    sigma_star2_out[gt+1,]=sigma_star2

    di_out[gt+1,]=di

    corr1_out[gt]=Sig1[1,2]/(sqrt(Sig1[1,1])*sqrt(Sig1[2,2]))
    corr2_out[gt]=Sig2[1,2]/(sqrt(Sig2[1,1])*sqrt(Sig2[2,2]))
    corr3_out[gt]=Sig3[1,2]/(sqrt(Sig3[1,1])*sqrt(Sig3[2,2]))

    # save data at every 10th iteration
    if (gt>nburn & gt%%10==0){
      theta_star11_T[(gt-nburn)/10]=theta_star11[0.3*H]
      theta_star12_T[(gt-nburn)/10]=theta_star12[0.5*H]
      theta_star13_T[(gt-nburn)/10]=theta_star13[0.7*H]
      theta_star21_T[(gt-nburn)/10]=theta_star21[0.8*H]
      theta_star22_T[(gt-nburn)/10]=theta_star22[0.5*H]
      theta_star23_T[(gt-nburn)/10]=theta_star23[0.2*H]

      sigma_star1_T[(gt-nburn)/10]=sigma_star1[0.3*L]
      sigma_star2_T[(gt-nburn)/10]=sigma_star2[0.7*L]

      V_h1_T[(gt-nburn)/10]=V_h1[0.3*H]
      V_h2_T[(gt-nburn)/10]=V_h2[0.8*H]
      V_h3_T[(gt-nburn)/10]=V_h3[0.5*H]

      Sig11_T[(gt-nburn)/10]=Sig1[1,1]
      Sig12_T[(gt-nburn)/10]=Sig1[1,2]
      Sig13_T[(gt-nburn)/10]=Sig1[2,2]

      Sig21_T[(gt-nburn)/10]=Sig2[1,1]
      Sig22_T[(gt-nburn)/10]=Sig2[1,2]
      Sig23_T[(gt-nburn)/10]=Sig2[2,2]

      Sig31_T[(gt-nburn)/10]=Sig3[1,1]
      Sig32_T[(gt-nburn)/10]=Sig3[1,2]
      Sig33_T[(gt-nburn)/10]=Sig3[2,2]
    }

    if (gt%%1000==0){
      print(gt)
    }
  }

  ## make posterior means and 95% credible intervals
  Y11_pred_pdfm=colMeans(Y21_pred_pdf[(nburn+1):nsim,])
  Y11_pred_pdfci=apply(Y21_pred_pdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))
  Y11_pred_cdfm=colMeans(Y21_pred_cdf[(nburn+1):nsim,])
  Y11_pred_cdfci=apply(Y21_pred_cdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))

  Y12_pred_pdfm=colMeans(Y22_pred_pdf[(nburn+1):nsim,])
  Y12_pred_pdfci=apply(Y22_pred_pdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))
  Y12_pred_cdfm=colMeans(Y22_pred_cdf[(nburn+1):nsim,])
  Y12_pred_cdfci=apply(Y22_pred_cdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))

  Y13_pred_pdfm=colMeans(Y23_pred_pdf[(nburn+1):nsim,])
  Y13_pred_pdfci=apply(Y23_pred_pdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))
  Y13_pred_cdfm=colMeans(Y23_pred_cdf[(nburn+1):nsim,])
  Y13_pred_cdfci=apply(Y23_pred_cdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))

  Y21_pred_pdfm=colMeans(Y11_pred_pdf[(nburn+1):nsim,])
  Y21_pred_pdfci=apply(Y11_pred_pdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))
  Y21_pred_cdfm=colMeans(Y11_pred_cdf[(nburn+1):nsim,])
  Y21_pred_cdfci=apply(Y11_pred_cdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))

  Y22_pred_pdfm=colMeans(Y12_pred_pdf[(nburn+1):nsim,])
  Y22_pred_pdfci=apply(Y12_pred_pdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))
  Y22_pred_cdfm=colMeans(Y12_pred_cdf[(nburn+1):nsim,])
  Y22_pred_cdfci=apply(Y12_pred_cdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))

  Y23_pred_pdfm=colMeans(Y13_pred_pdf[(nburn+1):nsim,])
  Y23_pred_pdfci=apply(Y13_pred_pdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))
  Y23_pred_cdfm=colMeans(Y13_pred_cdf[(nburn+1):nsim,])
  Y23_pred_cdfci=apply(Y13_pred_cdf[(nburn+1):nsim,],2,quantile,c(0.025,0.975))

  di_temp=matrix(0,3,N)
  di_post=rep(0,N)
  for (i in 1:N){
    di_temp[1,i]=sum(di_out[(nburn+1):nsim,i]==0)
    di_temp[2,i]=sum(di_out[(nburn+1):nsim,i]==1)
    di_temp[3,i]=sum(di_out[(nburn+1):nsim,i]==2)
  }

  for (i in 1:N){
    if (di_temp[1,i]>(di_temp[2,i]+di_temp[3,i])){di_post[i]=0}
    else if (di_temp[1,i]<(di_temp[2,i]+di_temp[3,i])){
      if (di_temp[2,i]>di_temp[3,i]){di_post[i]=1}
      else if (di_temp[2,i]<di_temp[3,i]){di_post[i]=2}
    }
  }

  corr1_post=mean(corr1_out[(nburn+1):nsim])
  corr2_post=mean(corr2_out[(nburn+1):nsim])
  corr3_post=mean(corr3_out[(nburn+1):nsim])

  corr1_ci=quantile(corr1_out[(nburn+1):nsim],probs=c(0.025,0.975))
  corr2_ci=quantile(corr2_out[(nburn+1):nsim],probs=c(0.025,0.975))
  corr3_ci=quantile(corr3_out[(nburn+1):nsim],probs=c(0.025,0.975))

  ## calculate VUS
  VUS1m=mean(VUS2_temp[(nburn+1):nsim])
  VUS1ci=quantile(VUS2_temp[(nburn+1):nsim],probs=c(0.025,0.975))
  VUS2m=mean(VUS1_temp[(nburn+1):nsim])
  VUS2ci=quantile(VUS1_temp[(nburn+1):nsim],probs=c(0.025,0.975))

  VUS1_out=VUS2_temp[(nburn+1):nsim]
  VUS2_out=VUS1_temp[(nburn+1):nsim]

  proc.time()-ptm
  ## return posterior means and 95% credible intervals
  ROC_return <- list( Y11_pred_pdfm, Y11_pred_pdfci, Y11_pred_cdfm, Y11_pred_cdfci,
                      Y12_pred_pdfm, Y12_pred_pdfci, Y12_pred_cdfm, Y12_pred_cdfci,
                      Y13_pred_pdfm, Y13_pred_pdfci, Y13_pred_cdfm, Y13_pred_cdfci,
                      Y21_pred_pdfm, Y21_pred_pdfci, Y21_pred_cdfm, Y21_pred_cdfci,
                      Y22_pred_pdfm, Y22_pred_pdfci, Y22_pred_cdfm, Y22_pred_cdfci,
                      Y23_pred_pdfm, Y23_pred_pdfci, Y23_pred_cdfm, Y23_pred_cdfci,
                      di_post, corr1_post, corr2_post, corr3_post, corr1_ci, corr2_ci,
                      corr3_ci,
                      ## calculate VUS
                      VUS1m, VUS1ci, VUS2m, VUS2ci, VUS1_out, VUS2_out
  )
}






