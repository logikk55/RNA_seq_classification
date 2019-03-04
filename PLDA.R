## Poisson Linear Discriminant Anallysis classifier
# This implementation is largely based on the PoiClaClu package
# Own implementation for the purpose of understanding it better

PLDA <- function(x,x_test, y, beta = 1, prior=NULL, method = c('mle','deseq','quantile')){
  # x : n x p matrix, where n indicates observations and p indicates features
  # x_test : m x p matrix for testing, where m indicates observations and p indicates features
  # y : class labels of each observation in x
  # beta : prior for Gamma(beta, beta) distribution
  # prior : prior probabilities for classes y. Vector of lenght unique(y)
  # method: method to compute size factor
  
  Xi_dot = rowSums(x)
  Xdot_j = colSums(x)
  Xdotdot = sum(x)
  
  Nhat = outer(Xi_dot, Xdot_j, '*')/Xdotdot
  
  classes = sort(unique(y))
  
  dhat <- matrix(0,nrow=length(classes),ncol=length(Xdot_j))
  for(i in 1:length(classes)){
    numerator = colSums(x[y==classes[i],])+beta
    denominator = colSums(Nhat[y==classes[i],])+beta
    division = numerator/denominator
    colnames(division) = NULL
    dhat[i,] = division
  }
 
  shat = compute_shat(x, x_test, method)  
  ghat = colSums(x)
  
  p = matrix(0,ncol = length(classes),nrow = nrow(x_test))
  for (i in 1:length(classes)){
    p[,i] = rowSums(x_test*1/log(dhat[i,]))-rowSums(outer(shat, ghat, '*')*1/dhat[i,])+log(prior[i])
    #p[,i] = exp(log_p) 
    #p[,i] = rowSums(scale(x_test, center=FALSE, scale=(1/log(dhat[i,])))) - rowSums(scale(outer(shat, ghat, '*'), center = FALSE, scale = (1/dhat[i,]))) + log(prior[i])
  }
  yhat = classes[apply(p,1,which.max)]
  return(yhat) #maybe good to return something else too
}


compute_shat <- function(x_train,x_test,method = c('mle','deseq','quatile')){
  # NOTE! For now only mle possible! Implementation for deseq and quantile is in progress
  # Compute size factor based on defined method
  # 'mle'
  # 'deseq'
  # 'quantile'
 if (method == 'mle')
   shat <- rowSums(x_test)/sum(x_train)
 return(shat)
}
