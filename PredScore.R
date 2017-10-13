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

fitModel <- function(model,dataset,n=NA){
  if ("stanreg" %in% class(model)){
    stan_lm( formula(model),
             data=dataset,
             prior=NULL, chains=1,
             iter=2*max(n,1000),
             show_messages=FALSE )
  } else if (class(model)=="lm"){
    lm( formula(model),
        data=dataset )
  }
}

getPred <- function(model,dataset,n=NA){
  if ("stanreg" %in% class(model)){
    posterior_predict(model,
                      newdata=dataset,
                      draws=n)
  } else if (class(model)=="lm"){
    predict(model,dataset)
  }
}

predscore <- function( obsdata,
                       model,
                       preddata = obsdata,
                       params = NA,
                       partitions = nrow(obsdata),
                       n.predict = nrow(obsdata) ){
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
    if (!is.na(params[1])){ ## if comparing model parameters
      results <- as.data.frame(matrix(NA,nrow=partitions,
                      ncol=length(params),dimnames=list(NULL,params)))
      resultsPval <- results
    }
    for (k in 1:partitions){
      testset_k <- obsdata[perm[idxsets[[k]]],]
      trainset_k <- obsdata[perm[(1:n)[
                    !(1:n %in% idxsets[[k]])]],]
      fit_k <- fitModel( model, trainset_k, n.predict )
      if (is.na(params[1])){  ## comparing at outcome-level
        pred_k <- getPred(fit_k,testset_k,n.predict) 
        results$q[idxsets[[k]]] <- qhat( testset_k[,yname], pred_k, isBin )
      } else { ## comparing coefficients
        if ("stanreg" %in% class(model)){
          n.compare <- min( n.predict,
                          nrow(model$stanfit),
                          nrow(fit_k$stanfit) )
          pred_k <- as.matrix(fit_k$stanfit)[sample.int(
                    nrow(fit_k$stanfit),n.compare),params]
          pred_orig <- as.matrix(model$stanfit)[sample.int(
                    nrow(model$stanfit),n.compare),params]
        } else {
          pred_k <- coef(fit_k)[params]
          pred_orig <- coef(model)[params]
        }
        results[k,] <- ksVec(pred_k,pred_orig)
        resultsPval[k,] <- ksVec(pred_k,pred_orig,"p.value")
      }
  }} else {
  ## Validation
    if (is.na(params[1])){  ## comparing outcomes
      pred <- getPred(model,preddata,n.predict)
      results <- data.frame( q=qhat(preddata[,yname], pred, isBin) )
    } else { ## comparing parameters
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
      results <- as.data.frame(matrix( cbind( ksVec(pred_orig,pred),
                                   ksVec(pred_orig,pred,"p.value") ),
                            nrow=2, ncol=2,
                            dimnames=list( params,
                                           c("ks.stat","p.value") )))
    }}
  
  print(paste("Prediction scoring complete.",
              round(as.numeric(Sys.time()-startTime)/60,2),
              "minutes elapsed."))
  if (identical(preddata,obsdata) & !is.na(params[1])){ 
    return(list(results=results,pvalues=resultsPval)) 
  } else { return(results) }
}



predscore.nonBayes <- function( obsdata,
                       model,
                       preddata = obsdata,
                       params = NA,
                       partitions = nrow(obsdata),
                       n.predict = 25 ){
  startTime <- Sys.time()
  
  ## Set up 
  yname <- deparse(attr(model$terms,"variables")[[2]])
  n <- nrow(obsdata)
  isBin <- FALSE
  partSize <- floor(n/partitions)
  pred_kMat <- NULL
  
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
    for (k in 1:partitions){
      testset_k <- obsdata[perm[idxsets[[k]]],]
      train_obs_k <- perm[(1:n)[!(1:n %in% idxsets[[k]])]]
      for (r in 1:n.predict){
        trainset_kr <- obsdata[sample( train_obs_k,
                            length(train_obs_k), replace=TRUE ),]
        fit_kr <- fitModel( model, trainset_kr )
        pred_kr <- getPred( fit_kr, testset_k )
        pred_kMat <- rbind( pred_kMat, pred_kr )
      }
      results$q[idxsets[[k]]] <- qhat(testset_k[,yname],pred_kMat,isBin)
    }

  ## Validation
  } else {  
    pred_vMat <- getPred(model,preddata)
    for (r in 1:(n.predict-1)){
      trainset_r <- obsdata[sample( 1:n, n, replace=TRUE ),]
      fit_r <- fitModel(model,trainset_r)
      pred_r <- getPred(fit_r,preddata)
      pred_vMat <- rbind( pred_vMat, pred_r )
    }
    results <- data.frame( #id=1:length(preddata),
                           q=qhat(preddata[,yname],pred_vMat,isBin) )
  }
  
  print(paste("Prediction scoring complete.",
              round(as.numeric(Sys.time()-startTime)/60,2),
              "minutes elapsed."))
  return(results)
}





plot.predscore <- function( q1, q2=NULL,
                            main1="cross-validation",
                            main2="validation"){
  plotcol <- 1; if (length(q2)>0){ plotcol <- 2 }
  par(mfrow=c(2,plotcol))
  
  ## Raw Quantiles (Uniform)
  hist(q1,breaks=25,main=main1)
  if (length(q2)>0){ hist(q2,breaks=25,main=main2) }
  
  plot(ecdf(q1),main="Empirical CDFs",xlab="q",ylab="Fn(q)",
       do.points=FALSE,verticals=TRUE)
  if (length(q2)>0){ plot(ecdf(q2),do.points=FALSE,verticals=TRUE,
        col="blue",add=TRUE) }
  plot(ecdf(runif(500)),do.points=FALSE,verticals=TRUE,
        col="red",add=TRUE)
  if (length(q2)>0){ legend("topleft",c(main1,main2,"uniform"),
      col=c("black","blue","red"), lty=1, cex=.5)
  } else { legend("topleft",c(main1,"uniform"),
      col=c("black","red"), lty=1, cex=.5) }
  
  ## Transformed Quantiles (Normal)
  x1 <- qnorm(q1); x1[x1==-Inf] <- NA; x1[x1==Inf] <- NA
  if (length(q2)>0){
    x2 <- qnorm(q2); x2[x2==-Inf] <- NA; x2[x2==Inf] <- NA }
  
  #qqnorm(x1, main=main1); qqline(x1,col="red")
  #legend("topleft","Normal pdf",lty=1,col="red",cex=.8)
  #if (length(q2)>0){ qqnorm(x2, main=main2); qqline(x2,col="red") }
  
  plot(ecdf(x1),main="Empirical CDFs",xlab="q",ylab="Fn(q)",
       do.points=FALSE,verticals=TRUE)
  if (length(q2)>0){ plot(ecdf(x2),col="blue",do.points=FALSE,
                          verticals=TRUE,add=TRUE) }
  plot(ecdf(rnorm(1000)),col="red",do.points=FALSE,
       verticals=TRUE,add=TRUE)
  if (length(q2)>0){ legend("topleft",c(main1,main2,"normal"),
         col=c("black","blue","red"), lty=1, cex=.5) 
  } else { legend("topleft",c(main1,"normal"),
         col=c("black","red"), cex=.5) }
}
