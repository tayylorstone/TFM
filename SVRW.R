# SVRW.r
# adapted from https://people.smp.uq.edu.au/DirkKroese/statbook/
SVRW<-function(ystar,h,omega2h,Vh){
  T = length(h);
  ## parameters for the Gaussian mixture
  pi = c(0.0073, .10556, .00002, .04395, .34001, .24566, .2575);
  mui = c(-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518,-1.08819) - 1.2704; 
  sig2i = c(5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261);
  sigi = sqrt(sig2i);
  ## sample S from a 7-point distrete distribution
  temprand = matrix(runif(T));
  q = matrix(rep(pi,each=T),T) * dnorm(matrix(rep(ystar,each=7),byrow=TRUE,,7),matrix(rep(h,each=7),byrow=TRUE,,7)
                                       +matrix(rep(mui,each=T),T),matrix(rep(sigi,each=T),T));
  q = q/matrix(rep(rowSums(q),each=7),byrow=TRUE,,7);
  S = matrix(7 - rowSums(matrix(rep(temprand,each=7),byrow=TRUE,,7)<t(apply(q,1,cumsum)))+1);
  ## sample h
  ####<Making a sparse matrix>####
  Sp=diag(T);                    #
  d=matrix(rep(0,T),1);          #
  sparse=rbind(d,Sp);            #
  sparse<-sparse[-(T+1),]        #
  ################################
  H = diag(T) - sparse;
  invOmegah = diag(T)*c(1/Vh, 1/omega2h*rep(1,T-1));
  d = matrix(mui[S]); invSigystar = diag(T)*c(1/sig2i[S]);
  Kh = t(H) %*% invOmegah %*% H + invSigystar;
  Ch = t(chol(Kh));
  hhat = solve(Kh,(invSigystar %*% (ystar-d)));
  h = hhat + solve(t(Ch),matrix(rnorm(T)));
  result = list(h,S);
  return (result)
}