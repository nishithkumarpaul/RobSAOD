wmean <-
function(x,w=c(rep(1,length(x)))){
swtx<-0
for(i in 1:length(x)){
swtx<-swtx+x[i]*w[i]
}
swt<-sum(w)
wmn<-swtx/swt
return(wmn)
}
