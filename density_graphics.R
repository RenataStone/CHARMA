
# density function
dchen=function(y, lambda, mu)
{
   log(0.5)/(1-exp(mu^lambda))*lambda*y^(lambda-1)*exp(log(0.5)/(1-exp(mu^lambda))*(1-exp(y^lambda))+y^lambda)
}
## ==============================================================

mar_b<-2.5
mar_e<-2.5
mar_c<-0.5
mar_d<-0.5
dist_text<-1.3
dist_tick<-0.5


x = seq(0,20, by = 0.001)
y = dchen(x, 0.8, 1)

name_lambda<-paste("lfixochen",".pdf",sep="")
pdf(file = nome2,width = 3.5, height = 3.8,family = "Times") 
par(mar=c(2.8, 2.7, 1.2, 1)) 
par(mgp=c(1.7, 0.45, 0))
plot(x, y, type="l", col=16, lty=1,lwd=1, xlab="y",
     ylab="f(y)", ylim = c(0,0.5))
curve(dchen(x, 0.8, 3),col=2,lty=2,lwd=1,add=T)
curve(dchen(x,0.8, 6),col=4,lty=4,lwd=1,add=T)
curve(dchen(x, 0.8, 10),col=1,lty=3,lwd=1,add=T)
legend(x=5.5, y=0.47, 
       c(expression(plain(lambda)==0.8~"," ~ plain(mu)==1),
         expression(plain(lambda)==0.8~"," ~ plain(mu)==3), 
         expression(plain(lambda)==0.8~"," ~ plain(mu)==6),
         expression(plain(lambda)==0.8~","~ plain(mu)==10)),
       lty=c(1,2,4,3),lwd=c(1,1,1,1),bty="n",col=c(16,2,4,1))
box();axis(1) 
dev.off() 


###

x = seq(0,20, by = 0.001)
y = dchen(x, 0.5, 5)

name_mu<-paste("mufixochen",".pdf",sep="")
pdf(file = nome2,width = 3.5, height = 3.8,family = "Times") 
par(mar=c(2.8, 2.7, 1.2, 1)) 
par(mgp=c(1.7, 0.45, 0))
plot(x, y, type="l", col=16, lty=1,lwd=1, xlab="y",
     ylab="f(y)", ylim = c(0,0.5))
curve(dchen(x, 0.7, 5),col=2,lty=2,lwd=1,add=T)
curve(dchen(x,0.9, 5),col=4,lty=4,lwd=1,add=T)
curve(dchen(x, 1.1, 5),col=1,lty=3,lwd=1,add=T)
legend(x=6.5, y=0.47, 
       c(expression(plain(lambda)==0.5~"," ~ plain(mu)==5),
         expression(plain(lambda)==0.7~"," ~ plain(mu)==5), 
         expression(plain(lambda)==0.9~"," ~ plain(mu)==5),
         expression(plain(lambda)==1.1~","~ plain(mu)==5)),
       lty=c(1,2,4,3),lwd=c(1,1,1,1),bty="n",col=c(16,2,4,1))
box();axis(1) 
dev.off() 

