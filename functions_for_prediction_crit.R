################################################################################
########## FUNCTION FOR GENERATING INITIAL DESIGN ##############################
################################################################################
dinicial <- function(comb)
{
  ind <- matrix(sample(seq(1:nrow(comb)),N,replace=TRUE),N,1)
  ind <- as.matrix(ind[order(ind)])
  di <- comb[ind,]
  d <- prod(round(eigen(t(di[,-1])%*%di[,-1], symmetric=TRUE, 
                        only.values=TRUE)$values, 6))
  list(di=di,d=d)
}
################################################################################
############ CONTROLS FOR NON-SINGULAR INITIAL DESIGN ##########################
################################################################################
SampleD <- function(comb)
{
  d <- 0
  while(d<10^(-6))
  {
    dini <- dinicial(comb)
    d <- dini$d
  }
  X <- as.matrix(dini$di)
  list(X=X,d=d)
}
################################################################################
########## FUNCTION FOR CODING THE LEVELS IN CUBIC REGION ######################
################################################################################
codifica <- function(x)
{
  xc <- (x-(max(x)+min(x))/2)/((max(x)-min(x))/2)
}
################################################################################
########### FUNCTION FOR CODING THE LEVELS IN SPHERICAL REGION ################
################################################################################
Sphcand <- function(t)
{ 
  candc <- t[apply(t^2,1,sum)==K,]
  for(k in (K-1):1)
  {
    candc <- rbind(candc,sqrt(K/k)*t[apply(t^2,1,sum)==k,])
  }
  candc <- rbind(candc,t[apply(t^2,1,sum)==0,])
  
  list(candc=candc)
}
################################################################################
################ FUNCTION TO COUNT TREATMENT NUMBER #############################
################################################################################
TreatLabels <- function(tcand)
{
  n <- nrow(tcand)
  Treat <- matrix(c(1,rep(0,n-1)),n,1)
  Label <- 1
  for (i in 2:n)
  {
    Label <- Label+1
    Treat[i] <- Label
    for (j in 1:(i-1))
    {
      if(min(as.numeric(tcand[i,]==tcand[j,]))==1) Treat[i] <- Treat[j]
    }
    if(Treat[i] < Label) Label <- Label-1
  }
  Treat <- as.matrix(cbind(Treat,tcand))
  list(Treat=Treat)
}

################################################################################
############ FUNCTION FOR EXPANDING THE X MATRIX ###############################
################################################################################
matrixX <- function(x,Terms)
{
  Nr <- nrow(x)
  X <- matrix(c(rep(1,Nr)),Nr,1)
  contador <- 0
  
  for(i in 1:K)
  {
    contador <- contador+1
    if(Terms[i]>0){X <- cbind(X,x[,i])}
  }     
  for(i in 1:K)
  {
    contador <- contador+1
    if(Terms[contador]>0){X <- cbind(X,(x[,i])^2)}
  }
  if(K>1)
  {
    for(i in 1:(K-1))
    {
      for(j in (i+1):K)
      {
        contador <- contador+1
        if(Terms[contador]>0){X <- cbind(X,x[,i]*x[,j])}
      }
    }
  }
  X <- matrix(X,nrow=Nr)
  list(X=X)
}

################################################################################
########### FUNCTION FOR MOMENTS MATRIX ########################################
################################################################################

mat.I <- function(Npar,K,Cubic){
  Num <-Npar+1
  S<-matrix(c(rep(0,Num)),Num,Num)
  S[1,1]<-1 #intercept
  
  if(Cubic=='Y'){
    sig1  <- 1/3
    sig2  <- 1/3
    sigd3 <- 1/5
    sig3  <- 1/9
    sig4  <- 1/9
  } 
  else{
    radius.sph<-sqrt(K)
    sig1  <- (radius.sph^2)/K
    sig2  <- (radius.sph^2)/K
    sigd3 <- (3*(radius.sph^4))/(K*(K+2))
    sig3  <- (radius.sph^4)/(K*(K+2))
    sig4  <- sig3
  }
  
  #1st row e 1st column we have sig1 and so on
  for (j in (K+2):(2*K+1))
    S[1,j]<- sig1
  for (i in (K+2):(2*K+1))
    S[i,1]<- sig1
  
  for(i in 2:(K+2))  
    S[i,i]<-sig2
  
  for (i in (K+1):(2*K)+1)
    for (j in (K+2) : (2*K+1) )
      ifelse(i==j, S[i,j]<-sigd3, S[i,j]<-sig3)
  
  ### interactions
  for (i in (2*K+2):Num)
    S[i,i]<-sig4
  
  list(S=S)
}

################################################################################
######### FUNCTION FOR CALCULATING THE CRITERIA: D, A or L, I and I_D ##########
################################################################################
criterion <- function(X1,M,W,moment,Cubic,vol,Kap)
{
  det <- prod(round(eigen(M, symmetric=TRUE, only.values=TRUE)$values, 6))
  critD <- (det/N)^(1/Npar)
  if(critD>0)
  {
    critL  <- 1
    critI  <- 1
    critI_D <- 1
    if(max(Kap[3],Kap[4],Kap[6],Kap[7],Kap[8],Kap[9])>0)
    { 
      Minv <- solve(M)
      critL <- 1/c(W%*%matrix(diag(Minv)[-1],nc=1))
      critI <- vol/sum(diag(moment%*%Minv))
      critI_D <- vol/sum(diag((Minv[-1,-1])%*%moment[-1,-1]))
    }
  }else{
        critD   <- 0
        critL  <- 0
        critI  <- 0
        critI_D <- 0
       }
   list(critD=critD,critL=critL,critI=critI,critI_D=critI_D)
}
################################################################################
########## FUNCTION TO EXCHANGE ROWS - Algoritmo point exchange ################
################################################################################
swap <- function(cand,D,crita,w,moment,Npar,Cubic,vol,Kap)
{

  busca <- 0
  n <- nrow(cand)
  for(i in 1:n)
  {
    if(D[1,1]!=cand[i,1])
    {  
      Xc <- D      
      Xc[1,] <- cand[i,]
      M <- t(Xc[,-1])%*%Xc[,-1]
      criteC  <- criterion(Xc,M,W,moment,Cubic,vol,Kap)
      criteD  <- criteC$critD
      criteL  <- criteC$critL
      criteI <- criteC$critI
      criteI_D  <- criteC$critI_D
      if(criteD>0)
      {
       if(max(Kap[2],Kap[3],Kap[7],Kap[9])>0)
        {
          critDP <- 0
          critLP <- 0
          critIP <- 0
          critI_DP <-0
          df2 <- N-nlevels(as.factor(Xc[,1]))
          critDF <- (N-df2)
          if(df2>0)
          {
            critDP  <- criteD/qf(prob1,Npar,df2)
            critLP  <- criteL/qf(prob2,1,df2)
            critIP  <- criteI/qf(prob3,1,df2)
            critI_DP <- criteI_D/qf(prob3,1,df2)
          }
        }else{
          critDP  <- 1
          critLP  <- 1
          critDF   <- 1
          critIP  <- 1
          critI_DP <- 1
          critDF   <- 1 
        }
        critc <- (criteD^Kap[1])*(critDP^Kap[2])*(criteL^Kap[4])*(critLP^Kap[3])*(critDF^Kap[5])*(criteI^Kap[6])*(critIP^Kap[7])*(criteI_D^Kap[8])*(critI_DP^Kap[9])
      }else {critc <- 0}
      
      if(critc>crita)
      {
        D <- Xc
        crita <- critc
        busca <- 1
      } 
    }
  }
  for(l in 2:N) 
  {
    if(D[l,1]!=D[(l-1),1])
    {  
      Xc <- D 
      for(i in 1:n)
      {
        if(D[l,1]!=cand[i,1])
        {
          Xc[l,] <- cand[i,]   
          M <- t(Xc[,-1])%*%Xc[,-1]
          criteC  <- criterion(Xc,M,W,moment,Cubic,vol,Kap)
          criteD  <- criteC$critD
          criteL  <- criteC$critL
          criteI <- criteC$critI
          criteI_DP  <- criteC$critI_D
          if(criteD>0)
          {
            if(max(Kap[2],Kap[3],Kap[7],Kap[9])>0)
            {
              critDP <- 0
              critLP <- 0
              critIP <- 0
              critI_DP <-0
              df2 <- N-nlevels(as.factor(Xc[,1]))
              critDF <- (N-df2)
              if(df2>0)
              {
                critDP  <- criteD/qf(prob1,Npar,df2)
                critLP  <- criteL/qf(prob2,1,df2)
                critIP  <- criteI/qf(prob3,1,df2)
                critI_DP <- criteI_D/qf(prob3,1,df2)
              }
            }else{
              critDP  <- 1
              critLP  <- 1
              critDF   <- 1
              critIP  <- 1
              critI_DP <- 1
              critDF   <- 1 
            }
            critc <- (criteD^Kap[1])*(critDP^Kap[2])*(criteL^Kap[4])*(critLP^Kap[3])*(critDF^Kap[5])*(criteI^Kap[6])*(critIP^Kap[7])*(criteI_D^Kap[8])*(critI_DP^Kap[9])
          }else {critc <- 0}
          
          if(critc>crita)
          {
            D <- Xc
            crita <- critc
            busca <- 1
          } 
        }
      }
    } 
  }
  list(D=D,crita=crita,busca=busca)
}

##################################################################################
######### FUNCTION TO DRIVE THE SEARCH ###########################################
##################################################################################
SearchTreat <- function()
{
  inicio <- Sys.time()
  Kap <- Kappa/sum(Kappa)
#  Kap0 <- Kap[1];Kap1 <- Kap[2];Kap2 <- Kap[3];Kap3 <- Kap[4];
#  Kap4 <- Kap[5];Kap5 <- Kap[6];Kap6 <- Kap[7];Kap7 <- Kap[8];Kap8 <- Kap[9];
  cand <-as.matrix(expand.grid(Levels))
  candc <- apply(cand,2,codifica)
  if(Cubic=='N')  candc <- as.matrix(Sphcand(candc)$candc)
  candc <- TreatLabels(candc)$Treat
  cand <- cbind(candc[,1],matrixX(candc[,-1],Terms)$X)
  moment <- as.matrix(mat.I(Npar,K,Cubic)$S)
  vol <- ifelse(Cubic=="Y",(2^K),((pi^(K/2))/gamma(K/2+1))*sqrt(K)^K)
  critall <- matrix(0,nr=Ntries)
  for(m in 1:Ntries) 
  {
    X <- SampleD(cand)$X
    M <- t(X[,-1])%*%X[,-1]
    criteC  <- criterion(X,M,W,moment,Cubic,vol,Kap)
    criteD  <- criteC$critD
    criteL  <- criteC$critL
    criteI  <- criteC$critI
    criteI_D <- criteC$critI_D
    critDF  <- 0
    if(criteD>0)
    {
      if(sum(Kap[2],Kap[3],Kap[5],Kap[7],Kap[9])>0)
      {
        critPED <- 0
        critPEL <- 0
        critPEI <- 0
        critPEI_D <- 0 
        critDF  <- nlevels(as.factor(X[,1]))
        df2 <- N-critDF
        if(df2>0)
        {
          critDP <- criteD/qf(prob1,Npar,df2)
          critLP <- criteL/qf(prob2,1,df2)
          critIP <- criteI/qf(prob3,1,df2)
          critI_DP <- criteI_D/qf(prob3,1,df2)
        }   
      }else{
            critDP  <- 1
            critLP  <- 1
            critIP  <- 1
            critI_DP <- 1
            critDF  <- 1
        }
        crite <- (criteD^Kap[1])*(critDP^Kap[2])*(criteL^Kap[4])*(critLP^Kap[3])*(critDF^Kap[5])*(criteI^Kap[6])*(critIP^Kap[7])*(criteI_D^Kap[8])*(critI_DP^Kap[9])   
    }else{crite <- 0}          
    busca <- 1
    while(busca==1)
    {
      Xs <- swap(cand,X,crite,W,moment,Npar,Cubic,vol,Kap)
      X <- Xs$D
      crite <- Xs$crita
      busca <- Xs$busca
    }
    critall[m] <- crite
    if(m==1)
    {
      critopt <- crite
      Xopt <- X
    }
    else
    {
      if(crite > critopt)
      {
        Xopt <- X
        critopt <- crite
      }
    }
    
#    list(critopt=critopt,Xopt.d=Xopt)
  }
  
  tfim <- Sys.time()
  time <- tfim-inicio
  
  ############################################################################################
  ##########################  Optimal Design #################################################
  ############################################################################################
  
  critDF <- nlevels(as.factor(Xopt[,1]))
  dfbest <- N-nlevels(as.factor(Xopt[,1]))
  LoF <- critDF-(Npar+1)
  M.opt <- t(Xopt[,-1])%*%Xopt[,-1]
  Kapall <- rep(1,length(Kap))
  copt <- criterion(Xopt,M.opt,W,moment,Cubic,vol,Kapall)
  D   <- copt$critD
  DP <- ifelse(dfbest>0,D/qf(prob1,Npar,dfbest),0)
  L   <- copt$critL
  LP <- ifelse(dfbest>0,L/qf(prob2,1,dfbest),0)
  Minv <- solve(M.opt)
  V   <- diag(Minv)
  I   <- copt$critI
  IP  <- ifelse(dfbest>0,I/qf(prob3,1,dfbest),0)
  I_D  <- copt$critI_D
  I_DP <- ifelse(dfbest>0,I_D/qf(prob3,1,dfbest),0)
  # CHECKING
  critFINAL <-(D^Kap[1])*(DP^Kap[2])*(LP^Kap[3])*(L^Kap[4])*(critDF^Kap[5])*(I^Kap[6])*(IP^Kap[7])*(I_D^Kap[8])*(I_DP^Kap[9])
  criteria.opt <-cbind(compound.crit=critopt,critFINAL=critFINAL,D=D,DP=DP,A=L,
                        AP=LP,I=I,IP=IP,I_D=I_D,I_DP=I_DP,N=N,df_PE=dfbest,df_LoF=LoF)
  
  XBest <- Xopt[,3:(K+2)]
  # PRINTING THE RESULTS
  print(time)
  print(list("Weights actually used for A type criterion"=W))
  print(list("Priorities (kappa) in the Compound Criteria"=Kap))
  print("Criterion values for each Try")
  stem(critall)
  print(list("Criterion values, all parts"=criteria.opt))
  print(list("Variances of each term for the Optimum Design"=V))
  print("Optimal Design")
  print(XBest)
  list(time=time,Kappa=Kap,Weights=W,XBest=XBest,Xopt=Xopt,critall=critall,criteria.opt=criteria.opt,V=V)
  }



Efficiency <- function(d)
{
  X1 <- TreatLabels(d)$Treat
  critDF <- nlevels(as.factor(X1[,1]))
  dfPE <- N-critDF
  X <- matrixX(X1[,-1],Terms)$X
  X1 <- cbind(X1[,1],X)
  M <- t(X)%*%X
  moment <- as.matrix(mat.I(Npar,K,Cubic)$S)
  vol <- ifelse(Cubic=="Y",(2^K),((pi^(K/2))/gamma(K/2+1))*sqrt(K)^K)
  Allcrit <- criterion(X1,M,W,moment,Cubic,vol)
  critD <- Allcrit$crit
  critA <- Allcrit$critL
  critI <- Allcrit$critI
  critId <- Allcrit$critId
  critDP <- ifelse(dfPE,critD/qf(prob1,Npar,dfPE),0)
  critAP <- ifelse(dfPE,critA/qf(prob2,1,dfPE),0)
  critIP <- ifelse(dfPE,critI/qf(prob3,1,dfPE),0)
  critIDP <- ifelse(dfPE,critId/qf(prob3,1,dfPE),0)
  critBias <- Allcrit$critBias
  crit <- c(critD=critD,critDP=critDP,critA=critA,critAP=critAP,critI=critI,critIP=critIP,critID=critId,critIDP=critIDP,critBias=critBias,dfPE=dfPE)
  return(crit)
}

