# adapted from https://people.smp.uq.edu.au/DirkKroese/statbook/
gamrand<-function(alpha,lambda){
  if (alpha>1)
  {
    d=alpha-1/3;
    c=1/sqrt(9*d); 
    flag=1;
    while (flag){
      Z=rnorm(1);
      if (Z>(-1/c)){
        V=(1+c*Z)^3;
        U=runif(1);
        flag=log(U)>(0.5*Z^2+d-d*V+d*log(V));  
      }
    }
    x=d*V/lambda;
  }
  else{
    x=gamrand(alpha+1,lambda);
    x=x*runif(1)^(1/alpha);
  }
};