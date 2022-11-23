wVarSs <-
function(x,w=c(rep(1,length(x))),center=0){
vwtx<-0
swt<-sum(w)
for(i in 1:length(x)){
vwtx<-vwtx+(x[i]*w[i]-center)^2
}
wvn<-vwtx
return(wvn)
}
