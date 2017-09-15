qhat <- function( obs, dist, isBin=FALSE ){
  if (length(obs)>1) qhat.vector(obs,dist,isBin)
  else qhat.vector(obs,dist,isBin)
}

qhat.scalar <- function( obs, dist, isBin ){
  if (isBin) return( sum(obs==dist)/length(dist) )
  else return( sum(obs>dist)/length(dist) )
}

qhat.vector <- function( obs, dist, isBin ){
  if (length(obs)!=ncol(dist)) return("error: length(obs) must equal ncol(dist)")
  else return( sapply(1:length(obs),
                      function(x){ qhat.scalar(obs[x],dist[,x],isBin) }) )
}

ksVec <- function( paramMat1, paramMat2, type="statistic" ){
  sapply( 1:ncol(paramMat1), function(i){
    as.numeric(ks.test(paramMat1[,i],paramMat2[,i])[type]) })
}

predscore <- function( obsdata,
                       model,
                       preddata = obsdata,
                       params = NA,
                       partitions = nrow(obsdata),
                       n.predict = nrow(obsdata),
                       type="stan_lm"){
  startTime <- Sys.time()
  
  ## Set up 
  yname <- deparse(attr(model$terms,"variables")[[2]])
  n <- nrow(obsdata)
  isBin <- FALSE
  if (partitions==nrow(obsdata) & !is.na(params[1])){
    partitions <- 10
    print("Defaulting to 10 partitions...") }
  partSize <- floor(n/partitions)
  
  ## Check for improper inputs
  if (partitions>n){ 
    stop("Must have fewer than n partitions") }
  if ((!is.na(params[1])) & (partSize < 25)){ 
    stop("When comparing parameters, partitions should contain at
         least 25 observations to allow for estimation.") }
  if ((!is.na(params[1])) & sum(!(params %in% names(model$stanfit)))>0){
    stop("Items in params must be names of model parameters") }
  
  ## Prep data partitions
  perm <- sample.int(n)
  idxsets <- lapply(1:partitions,function(k){
                (1 + (k-1)*partSize):(k*partSize) })
  if (partitions*partSize != n){
    idxsets[[partitions]] <- c(idxsets[[partitions]],
                             perm[(k*partSize+1):n]) }

  ## Cross-validation
  if (identical(preddata,obsdata)){
    results <- data.frame( set=rep(1:partitions,each=partSize),
                           q=rep(NA,partitions*partSize) )
    if (!is.na(params[1])){
      results <- as.data.frame(matrix(NA,nrow=partitions,
                      ncol=length(params),dimnames=list(NULL,params)))
      resultsPval <- results
    }
    for (k in 1:partitions){
      testset_k <- obsdata[perm[idxsets[[k]]],]
      trainset_k <- obsdata[perm[(1:n)[
                    !(1:n %in% idxsets[[k]])]],]
      fit_k <- stan_lm( formula(model),
                 data=trainset_k,
                 prior=NULL, chains=1,
                 iter=2*max(n.predict,1000),
                 show_messages=FALSE )
      if (is.na(params[1])){
        pred_k <- posterior_predict(fit_k,
                 newdata=testset_k,
                 draws=n.predict)
        results$q[idxsets[[k]]] <- qhat( testset_k[,yname], pred_k, isBin )
      } else {
        n.compare <- min( n.predict,
                          nrow(model$stanfit),
                          nrow(fit_k$stanfit) )
        pred_k <- as.matrix(fit_k$stanfit)[sample.int(
                    nrow(fit_k$stanfit),n.compare),params]
        pred_orig <- as.matrix(model$stanfit)[sample.int(
                    nrow(model$stanfit),n.compare),params]
        results[k,] <- ksVec(pred_k,pred_orig)
        resultsPval[k,] <- ksVec(pred_k,pred_orig,"p.value")
      }
  }} else {
  ## Validation
    if (is.na(params[1])){
      pred <- posterior_predict( model, newdata=preddata,
                                draws=n.predict )
      results <- data.frame( q=qhat(preddata[,yname], pred, isBin) )
    } else {
      fit_pred <- stan_lm( formula(model),
                       data=preddata,
                       prior=NULL, chains=1,
                       iter=2*max(n.predict,1000),
                       show_messages=FALSE )
      n.compare <- min( n.predict,
                        nrow(model$stanfit),
                        nrow(fit_pred$stanfit) )
      pred_orig <- as.matrix(model$stanfit)[sample.int(
                        nrow(model$stanfit),n.compare),params]
      pred <- as.matrix(fit_pred$stanfit)[sample.int(
                        nrow(fit_pred$stanfit),n.compare),params]
      results <- as.matrix( cbind( ksVec(pred_orig,pred),
                                   ksVec(pred_orig,pred,"p.value") ),
                            dimnames=list( c("ks.stat","p.value"),
                                           params ) )
  }}
  
  print(paste("Prediction scoring complete.",
              round(as.numeric(Sys.time()-startTime)/60,2),
              "minutes elapsed."))
  if (identical(preddata,obsdata) & !is.na(params[1])){ 
    return(list(results=results,pvalues=resultsPval)) 
  } else { return(results) }
}


plot.predscore <- function( q1, q2=NULL,
                            main1="cross-validation",
                            main2="validation"){
  plotcol <- 2; if (length(q2)>0){ plotcol <- 2 }
  par(mfrow=c(2,plotcol))
  
  ## Raw Quantiles (Uniform)
  hist(q1,breaks=25,main="main1")
  if (length(q2)>0){ hist(q2,breaks=25,main="main2") }
  
  plot(ecdf(q1),main="Empirical CDFs",xlab="q",ylab="Fn(q)",
       do.points=FALSE,verticals=TRUE)
  if (length(q2)>0){ plot(ecdf(q2),do.points=FALSE,verticals=TRUE,
        col="blue",add=TRUE) }
  plot(ecdf(runif(500)),do.points=FALSE,verticals=TRUE,
        col="red",add=TRUE)
  if (length(q2)>0){ legend("topleft",c(main1,main2,"uniform"),
      col=c("black","blue","red"), lty=1)
  } else { legend("topleft",c(main1,"uniform"),
      col=c("black","red"), lty=1) }
  
  ## Transformed Quantiles (Normal)
  x1 <- qnorm(q1); x1[x1==-Inf] <- NA; x1[x1==Inf] <- NA
  if (length(q2)>0){
    x2 <- qnorm(q2); x2[x2==-Inf] <- NA; x2[x2==Inf] <- NA }
  
  qqnorm(x1, main=main1); qqline(x1,col="red")
  legend("topleft","Normal pdf",lty=1,col="red",cex=.8)
  if (length(q2)>0){ qqnorm(x2, main=main2); qqline(x2,col="red") }
  
  plot(ecdf(x1),main="Empirical CDFs",xlab="q",ylab="Fn(q)",
       do.points=FALSE,verticals=TRUE)
  if (length(q2)>0){ plot(ecdf(x2),col="blue",do.points=FALSE,
                          verticals=TRUE,add=TRUE) }
  plot(ecdf(rnorm(1000)),col="red",do.points=FALSE,
       verticals=TRUE,add=TRUE)
  if (length(q2)>0){ legend("topleft",c(main1,main2,"normal"),
         col=c("black","blue","red"), lty=1) 
  } else { legend("topleft",c(main1,"normal"),
         col=c("black","red")) }
}
