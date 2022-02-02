f_df_lambda <- function(mm, df){
  d2 = svd(mm)$d^2
  ff <- function(la)  (df-sum(sapply(d2,function(x) x/(x+la ))))^2
  guess= (mean(sqrt(d2))^2)/(df/length(d2)) - mean(sqrt(d2))^2
  opt = optim(par = guess, ff, method = "Brent", lower = 0, upper = 1e8)
  
  structure(opt$par, d = c(dfh = sum(sapply(d2,function(x) x/(x+opt$par))), df = df))
}

# get deviance of ridge glmnet given the df
f_df_glmnet <- function(X, y, df, family = "binomial"){
  L = f_df_lambda(X, df)
  gn <- glmnet(x = X, y = y, family = family, intercept = F,
               alpha = 0, lambda = L/nrow(X), standardize = F)
  #browser()
  structure(gn$nulldev*(1-gn$dev.ratio), null_dev = gn$nulldev, L = c(L/nrow(X)), d = attr(L,"d"))
}

# p value of regularized model using LRT 
# could be problematic be careful for some cases ! 
p_df_glmnet <- function(X, y, df, family = "binomial"){
  re = f_df_glmnet(X,y,df, family)
  d_dev = c(attr(re,"null_dev")-re)
  structure( -expm1( pchisq( d_dev,df = df, log.p = T) ), dev = re)
}