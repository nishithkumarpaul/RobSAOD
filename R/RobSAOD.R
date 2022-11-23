RobSAOD <-
function(x,ngroup=2,noGrpMem=c(ncol(x)/2,ncol(x)/2)){
if (ngroup!=length(noGrpMem))
stop("Please Enter the Individual Column Number of Each Group and Check the Number of group")
x<-as.matrix(x)
#####################################################
###     Gradual obs                               ###
#####################################################
grwdat<-NULL
nCol<-0
for (i in 1:ngroup){
ncol1<-nCol+1
nCol<-nCol+noGrpMem[i]
grwdat[[i]]<-x[,ncol1:nCol]
}

#####################################################
###     Weight calculation                        ###
#####################################################
nl<-length(noGrpMem)
nCol<-0

  med<-NULL
  mad<-NULL
  dW<-NULL
for(i in 1:nl){
  ncol1<-nCol+1
  nCol<-nCol+noGrpMem[i]

  med[[i]]<-abs((x[,ncol1:nCol]-apply(x[,ncol1:nCol],1,median)))
  mad[[i]]<-(apply(x[,ncol1:nCol],1,mad))

  dW[[i]]<-(1.96)^2/((med[[i]])/mad[[i]])^2
}
dataWeight<-NULL
 for(i in 1:nl){
   dataWeight<-cbind(dataWeight,dW[[i]])
 }
nx<-length(dataWeight[,1])
cx<-length(dataWeight[1,])

weightDat<-matrix(rnorm(nx*cx), nrow=nx,ncol=cx)

   for( f in 1:nx)
     {
      for( g in 1:cx)
        {
           weightDat[f,g]=min(1,dataWeight[f,g])
        }
      }

#################################################
### put weight 1 if its weight is greater .2 ####
#################################################

weightDat1<-matrix(rnorm(nx*cx), nrow=nx,ncol=cx)

   for( p in 1:nx)
     {
      for( q in 1:cx)
        {
         if(weightDat[p,q]>=.223)
            weightDat1[p,q]<-1
         else
           weightDat1[p,q]<-weightDat[p,q]
         }
      }

w<-weightDat1
#########################################################
###  Weighted Mean calculation each group and polled ####
#########################################################
nCol<-0
  xwm<-NULL
  xwms<-NULL
for(i in 1:nl){
  ncol1<-nCol+1
  nCol<-nCol+noGrpMem[i]


  for(j in 1:(nrow(x[,ncol1:nCol]))){
  xwms[j]<-wmean(x[j,ncol1:nCol],w[j,ncol1:nCol])
  }
  xwm[[i]]<-xwms
}
xw<-NULL
for(j in 1:(nrow(x))){
xw[j]<-wmean(x[j,],w[j,])
}
#xw<-unlist(xw1)
#x11<-unlist(x1b)
#x12<-unlist(x2b)
#x13<-unlist(x3b)
#xwm<-unlist(xw)
##########################################
###  rwi calculation                  ####
##########################################
snk<-sum(noGrpMem)
pnk<-prod(noGrpMem)
xnk<-0
xnk1<-NULL
for(i in 1:nl){
xnk1<-xnk+noGrpMem[i]*(xwm[[i]]-xw)^2
xnk<-xnk1
}

rwi<-sqrt((snk/pnk)*xnk)

#############################################
#### Calculation of  weighted variance   ####
#############################################

nCol<-0
  xwv<-NULL
  xwvs<-NULL
for(i in 1:nl){
  ncol1<-nCol+1
  nCol<-nCol+noGrpMem[i]


  for(j in 1:(nrow(x[,ncol1:nCol]))){
  xwvs[j]<-wVarSs(x[j,ncol1:nCol],w[j,ncol1:nCol],xwm[[i]][j])
  }
 xwv[[i]]<-xwvs
}
proCoeff<-sum(1/noGrpMem)/sum(noGrpMem-1)

Sxwv<-0
for(i in 1:nl){
Sxwv<-Sxwv+xwv[[i]]
}
#############################################
#### Calculation swi                     ####
#############################################
swi<-sqrt(proCoeff*Sxwv)

sNot<-quantile(swi, probs = seq(0, 1, 0.05))

########################################
###   Find coeficient of variation    ##
########################################

coVarw<-NULL
for(kt in 1:length(sNot))
{
dIw<-rwi/(swi+sNot[kt])
coVarw[kt]<-(sd(dIw)/abs(mean(dIw)))*100
}
dExw<-rwi/(swi+sNot[which(coVarw==min(coVarw))])
if (nl<3){
pv1<-dt(dExw, df=(sum(noGrpMem)-nl))
pv<-pv1
  } else pv<-dt(dExw, df=sum(noGrpMem)-nl)/(2*(nl))
return(pv)
}
