load('mcmc_samples.RData')

posterior.full = samples.full



pdf('params.pdf',width=6.5,height=6.5)

par(mar=c(5,5,1,1),
    oma=c(0,0,1.5,1.5))

layout(matrix(1:8,4,2))

param.names = c(
  'wAlb SF shape',
  'wAlb DF shape',
  'WT SF shape',
  'WT DF shape',
  'wAlb SF rate',
  'wAlb DF rate',
  'WT SF rate',
  'WT DF rate')

for(ii in 1:8){
  d.armstrong = stats::density(posterior.armstrong[,ii])
  d.full = stats::density(posterior.full[,ii])
  plot(-100,-100,
       xlim = range(c(d.armstrong$x,d.full$x)),
       ylim = range(c(d.armstrong$y,d.full$y)),
       xlab = param.names[ii],
       ylab = 'Density', las = 1)
  polygon(
    c(d.armstrong$x,d.armstrong$x[1]),
    c(d.armstrong$y,d.armstrong$y[1]),
    col=rgb(0,0,0,0.4))
  if(ii%%2){
    polygon(
      c(d.full$x,d.full$x[1]),
      c(d.full$y,d.full$y[1]),
      col=rgb(0,0,1,0.4))
    legend('topright',legend=c('Prior','Posterior'),
           fill=c(rgb(0,0,0,0.4),rgb(0,0,1,0.4)),bty='n')
  } else {
    polygon(
      c(d.full$x,d.full$x[1]),
      c(d.full$y,d.full$y[1]),
      col=rgb(1,0,0,0.4))
    legend('topright',legend=c('Prior','Posterior'),
           fill=c(rgb(0,0,0,0.4),rgb(1,0,0,0.4)),bty='n')
  }
}

dev.off()



pdf('model_fit.pdf',width=6.5,height=6.5)

par(mar=c(5,5,1,1),
    oma=c(0,0,1.5,1.5))

layout(matrix(1:4,2,2))

t = seq(min(d$Day.post.infection)-1,max(d$Day.post.infection)+1,by=0.1)

plot(t,pgamma(t,params.map[1],params.map[4+1]),type='l',
     xlim=c(4.9,10.1),ylim=c(0,1),xaxs='i',yaxs='i',
     xlab='Days post infection',ylab='% disseminated',las=1)
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(1,4+1))
quantile50 = numeric()
for(ii in 1:nrow(post)){
  lines(t,pgamma(t,post[ii,1],post[ii,2]),col=rgb(0,0,1,0.15))
  quantile50[ii] = uniroot(
    function(x){pgamma(x,post[ii,1],post[ii,2])-0.5},interval=c(0,100))$root
}
quantile50 = quantile(quantile50,c(0.025,0.5,0.975))
# 2.5%      50%    97.5% 
# 7.717196 8.378418 9.010905 
arrows(max(quantile50),0.5,4.9,0.5,col='gray',length=0.075)
arrows(quantile50,rep(0.5,3),quantile50,rep(0,3),col='gray',length=0.075)
lines(t,pgamma(t,params.map[1],params.map[4+1]),col='cyan',lwd=3)
points(
  d$Day.post.infection[par.ind==1],
  d$Positive[par.ind==1]/d$Total[par.ind==1],
  pch=19)
segments(
  x0 = d$Day.post.infection[par.ind==1],
  y0 = qbeta(0.025,1+d$Positive[par.ind==1],1+d$Negative[par.ind==1]),
  x1 = d$Day.post.infection[par.ind==1],
  y1 = qbeta(0.975,1+d$Positive[par.ind==1],1+d$Negative[par.ind==1]))
mtext('SF',3,line=1)

plot(t,pgamma(t,params.map[3],params.map[4+3]),type='l',
     xlim=c(4.9,10.1),ylim=c(0,1),xaxs='i',yaxs='i',
     xlab='Days post infection',ylab='% disseminated',las=1)
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(3,4+3))
for(ii in 1:nrow(post)){
  lines(t,pgamma(t,post[ii,1],post[ii,2]),col=rgb(0,0,1,0.15))
}
lines(t,pgamma(t,params.map[3],params.map[4+3]),col='cyan',lwd=3)
points(
  d$Day.post.infection[par.ind==3],
  d$Positive[par.ind==3]/d$Total[par.ind==3],
  pch=19)
segments(
  x0 = d$Day.post.infection[par.ind==3],
  y0 = qbeta(0.025,1+d$Positive[par.ind==3],1+d$Negative[par.ind==3]),
  x1 = d$Day.post.infection[par.ind==3],
  y1 = qbeta(0.975,1+d$Positive[par.ind==3],1+d$Negative[par.ind==3]))

plot(t,pgamma(t,params.map[2],params.map[4+2]),type='l',
     xlim=c(4.9,10.1),ylim=c(0,1),xaxs='i',yaxs='i',
     xlab='Days post infection',ylab='% disseminated',las=1)
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(2,4+2))
quantile50 = numeric()
for(ii in 1:nrow(post)){
  lines(t,pgamma(t,post[ii,1],post[ii,2]),col=rgb(1,0,0,0.15))
  quantile50[ii] = uniroot(
    function(x){pgamma(x,post[ii,1],post[ii,2])-0.5},interval=c(0,100))$root
}
quantile50 = quantile(quantile50,c(0.025,0.5,0.975))
# 2.5%      50%    97.5% 
# 6.029424 6.858753 7.620342 
arrows(max(quantile50),0.5,4.9,0.5,col='gray',length=0.075)
arrows(quantile50,rep(0.5,3),quantile50,rep(0,3),col='gray',length=0.075)
lines(t,pgamma(t,params.map[2],params.map[4+2]),col='yellow',lwd=3)
points(
  d$Day.post.infection[par.ind==2],
  d$Positive[par.ind==2]/d$Total[par.ind==2],
  pch=19)
segments(
  x0 = d$Day.post.infection[par.ind==2],
  y0 = qbeta(0.025,1+d$Positive[par.ind==2],1+d$Negative[par.ind==2]),
  x1 = d$Day.post.infection[par.ind==2],
  y1 = qbeta(0.975,1+d$Positive[par.ind==2],1+d$Negative[par.ind==2]))
mtext('DF',3,line=1)
mtext('wAlb',4,line=1)

plot(t,pgamma(t,params.map[4],params.map[4+4]),type='l',
     xlim=c(4.9,10.1),ylim=c(0,1),xaxs='i',yaxs='i',
     xlab='Days post infection',ylab='% disseminated',las=1)
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(4,4+4))
for(ii in 1:nrow(post)){
  lines(t,pgamma(t,post[ii,1],post[ii,2]),col=rgb(1,0,0,0.15))
}
lines(t,pgamma(t,params.map[4],params.map[4+4]),col='yellow',lwd=3)
points(
  d$Day.post.infection[par.ind==4],
  d$Positive[par.ind==4]/d$Total[par.ind==4],
  pch=19)
segments(
  x0 = d$Day.post.infection[par.ind==4],
  y0 = qbeta(0.025,1+d$Positive[par.ind==4],1+d$Negative[par.ind==4]),
  x1 = d$Day.post.infection[par.ind==4],
  y1 = qbeta(0.975,1+d$Positive[par.ind==4],1+d$Negative[par.ind==4]))
mtext('WT',4,line=1)

dev.off()



pdf('eip_predictions.pdf',width=6.5,height=6.5)

par(mar=c(5,5,1,1),
    oma=c(0,0,1.5,1.5))

layout(matrix(1:4,2,2))

t = seq(0,21,by=0.1)

plot(t,dgamma(t,shape=params.map[1],rate=params.map[4+1]),type='l',
     ylim=c(0,0.2),
     xlab='EIP',ylab='Probability density',las=1)
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(1,4+1))
for(ii in 1:nrow(post)){
  lines(t,dgamma(t,shape=post[ii,1],rate=post[ii,2]),col=rgb(0,0,1,0.15))
}
lines(t,dgamma(t,shape=params.map[1],rate=params.map[4+1]),col='cyan',lwd=3)
mtext('SF',3,line=1)

plot(t,dgamma(t,shape=params.map[3],rate=params.map[4+3]),type='l',
     ylim=c(0,0.4),
     xlab='EIP',ylab='Probability density',las=1)
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(3,4+3))
for(ii in 1:nrow(post)){
  lines(t,dgamma(t,shape=post[ii,1],rate=post[ii,2]),col=rgb(0,0,1,0.15))
}
lines(t,dgamma(t,shape=params.map[3],rate=params.map[4+3]),col='cyan',lwd=3)

plot(t,dgamma(t,shape=params.map[2],rate=params.map[4+2]),type='l',
     ylim=c(0,0.2),
     xlab='EIP',ylab='Probability density',las=1)
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(2,4+2))
for(ii in 1:nrow(post)){
  lines(t,dgamma(t,shape=post[ii,1],rate=post[ii,2]),col=rgb(1,0,0,0.15))
}
lines(t,dgamma(t,shape=params.map[2],rate=params.map[4+2]),col='yellow',lwd=3)
mtext('DF',3,line=1)
mtext('wAlb',4,line=1)

plot(t,dgamma(t,shape=params.map[4],rate=params.map[4+4]),type='l',
     ylim=c(0,0.4),
     xlab='EIP',ylab='Probability density',las=1)
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(4,4+4))
for(ii in 1:nrow(post)){
  lines(t,dgamma(t,shape=post[ii,1],rate=post[ii,2]),col=rgb(1,0,0,0.15))
}
lines(t,dgamma(t,params.map[4],params.map[4+4]),col='yellow',lwd=3)
mtext('WT',4,line=1)

dev.off()




pdf('pr_surv_eip.pdf',width=6.5,height=3.25)

par(mar=c(3.5,2,2,1),
    oma=c(0,2.5,0,0),
    xpd=TRUE)

layout(matrix(1:3,1,3))

lifespan = 4

post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(1,4+1))
pr.surv.eip_1_4 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(2,4+2))
pr.surv.eip_2_4 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(3,4+3))
pr.surv.eip_3_4 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(4,4+4))
pr.surv.eip_4_4 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})

lifespan = 7

post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(1,4+1))
pr.surv.eip_1_7 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(2,4+2))
pr.surv.eip_2_7 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(3,4+3))
pr.surv.eip_3_7 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(4,4+4))
pr.surv.eip_4_7 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})

lifespan = 10

post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(1,4+1))
pr.surv.eip_1_10 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(2,4+2))
pr.surv.eip_2_10 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(3,4+3))
pr.surv.eip_3_10 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})
post = getSample(out.full, start = 5e4, end = NULL, thin = 5e2, whichParameters = c(4,4+4))
pr.surv.eip_4_10 = sapply(
  1:nrow(post),
  function(ii){
    integrate(f=function(n){exp(-n/lifespan)*dgamma(n,shape=post[ii,1],rate=post[ii,2])},lower=0,upper=28)$value})

boxplot(
  pr.surv.eip_1_4,pr.surv.eip_2_4,pr.surv.eip_3_4,pr.surv.eip_4_4,
  col=c('blue','red','blue','red'),
  names=c('SF','DF','SF','DF'),notch=T,
  ylim=c(0,1),yaxs='i',las=1,border=c('blue','red','blue','red'),
  ylab='')
mtext('wAlb',side=1,at=1.5,line=2.25,cex=0.7)
mtext('WT',side=1,at=3.5,line=2.25,cex=0.7)
legend('topleft',col=rgb(c(0,1),c(0,0),c(1,0),c(1,1)),legend=c('SF','DF'),pch=15,bty='n')
mtext('Average lifespan = 4 d',line=0.5)
mtext('Probability of surviving the EIP',2,line=3,cex=1)

boxplot(
  pr.surv.eip_1_7,pr.surv.eip_2_7,pr.surv.eip_3_7,pr.surv.eip_4_7,
  col=c('blue','red','blue','red'),
  names=c('SF','DF','SF','DF'),notch=T,
  ylim=c(0,1),yaxs='i',las=1,border=c('blue','red','blue','red'),
  ylab='')
mtext('wAlb',side=1,at=1.5,line=2.25,cex=0.7)
mtext('WT',side=1,at=3.5,line=2.25,cex=0.7)
legend('topleft',col=rgb(c(0,1),c(0,0),c(1,0),c(1,1)),legend=c('SF','DF'),pch=15,bty='n')
mtext('Average lifespan = 7 d',line=0.5)

boxplot(
  pr.surv.eip_1_10,pr.surv.eip_2_10,pr.surv.eip_3_10,pr.surv.eip_4_10,
  col=c('blue','red','blue','red'),
  names=c('SF','DF','SF','DF'),notch=T,
  ylim=c(0,1),yaxs='i',las=1,border=c('blue','red','blue','red'),
  ylab='')
mtext('wAlb',side=1,at=1.5,line=2.25,cex=0.7)
mtext('WT',side=1,at=3.5,line=2.25,cex=0.7)
legend('topleft',col=rgb(c(0,1),c(0,0),c(1,0),c(1,1)),legend=c('SF','DF'),pch=15,bty='n')
mtext('Average lifespan = 10 d',line=0.5)

dev.off()



pdf('relative_R0_wAlb.pdf',width=3.25,height=3.25)

par(mar=c(4.5,4.5,0.5,0.5),
    oma=c(0,0,0,0),
    xpd=FALSE)

boxplot(
  pr.surv.eip_2_4/pr.surv.eip_1_4,
  pr.surv.eip_2_7/pr.surv.eip_1_7,
  pr.surv.eip_2_10/pr.surv.eip_1_10,
  col='red',
  names=c(4,7,10),notch=T,
  las=1,ylim=c(0,3.3),yaxs='i',
  xlab='Average lifespan',
  ylab=expression('R'['0,DF']*' / R'['0,SF']*' for wAlb'))
abline(h=1,lty=2)

dev.off()

quantile(pr.surv.eip_2_4/pr.surv.eip_1_4,c(0.025,0.5,0.975))
# 2.5%      50%    97.5% 
# 1.339913 1.630907 2.100147 
quantile(pr.surv.eip_2_7/pr.surv.eip_1_7,c(0.025,0.5,0.975))
# 2.5%      50%    97.5% 
# 1.131652 1.258486 1.464587  
quantile(pr.surv.eip_2_10/pr.surv.eip_1_10,c(0.025,0.5,0.975))
# 2.5%      50%    97.5% 
# 1.066245 1.155353 1.284370  



pdf('relative_R0_WT.pdf',width=3.25,height=3.25)

par(mar=c(4.5,4.5,0.5,0.5),
    oma=c(0,0,0,0),
    xpd=FALSE)

boxplot(
  pr.surv.eip_4_4/pr.surv.eip_3_4,
  pr.surv.eip_4_7/pr.surv.eip_3_7,
  pr.surv.eip_4_10/pr.surv.eip_3_10,
  col='red',
  names=c(4,7,10),notch=T,
  las=1,ylim=c(0,3.3),yaxs='i',
  xlab='Average lifespan',
  ylab=expression('R'['0,DF']*' / R'['0,SF']*' for WT'))
abline(h=1,lty=2)

dev.off()



pdf('relative_R0_ratio.pdf',width=3.25,height=3.25)

par(mar=c(4.5,4.5,0.5,0.5),
    oma=c(0,0,0,0),
    xpd=FALSE)

boxplot(
  (pr.surv.eip_2_4/pr.surv.eip_1_4)/(pr.surv.eip_4_4/pr.surv.eip_3_4),
  (pr.surv.eip_2_7/pr.surv.eip_1_7)/(pr.surv.eip_4_7/pr.surv.eip_3_7),
  (pr.surv.eip_2_10/pr.surv.eip_1_10)/(pr.surv.eip_4_10/pr.surv.eip_3_10),
  col='red',
  names=c(4,7,10),notch=T,
  las=1,ylim=c(0,2),yaxs='i',
  xlab='Average lifespan',
  ylab=expression('Ratio of R'['0,DF']*' / R'['0,SF']*' for wAlb:WT'))
abline(h=1,lty=2)

dev.off()



pdf('pr_surv_eip_oddsratios.pdf',width=3.25,height=3.25)

par(mar=c(4.5,4.5,0.5,0.5),
    oma=c(0,0,0,0),
    xpd=FALSE)

boxplot(
  pr.surv.eip_1_4/(1-pr.surv.eip_1_4)*(1-pr.surv.eip_3_4)/pr.surv.eip_3_4,
  pr.surv.eip_1_7/(1-pr.surv.eip_1_7)*(1-pr.surv.eip_3_7)/pr.surv.eip_3_7,
  pr.surv.eip_1_10/(1-pr.surv.eip_1_10)*(1-pr.surv.eip_3_10)/pr.surv.eip_3_10,
  pr.surv.eip_2_4/(1-pr.surv.eip_2_4)*(1-pr.surv.eip_4_4)/pr.surv.eip_4_4,
  pr.surv.eip_2_7/(1-pr.surv.eip_2_7)*(1-pr.surv.eip_4_7)/pr.surv.eip_4_7,
  pr.surv.eip_2_10/(1-pr.surv.eip_2_10)*(1-pr.surv.eip_4_10)/pr.surv.eip_4_10,
  col=c(rep('blue',3),rep('red',3)),
  names=rep(c(4,7,10),2),notch=T,
  las=1,
  xlab='Average lifespan',
  ylab='OR of surviving EIP (wAlb:WT)')
abline(h=1,lty=2)
legend('topright',legend=c('SF','DF'),
       fill=c('blue','red'),bty='n')

dev.off()



pdf('relative_R0_wAlb_WT.pdf',width=3.25,height=3.25)

par(mar=c(4.5,4.5,0.5,0.5),
    oma=c(0,0,0,0),
    xpd=FALSE)

boxplot(
  pr.surv.eip_1_4/pr.surv.eip_3_4,
  pr.surv.eip_1_7/pr.surv.eip_3_7,
  pr.surv.eip_1_10/pr.surv.eip_3_10,
  pr.surv.eip_2_4/pr.surv.eip_4_4,
  pr.surv.eip_2_7/pr.surv.eip_4_7,
  pr.surv.eip_2_10/pr.surv.eip_4_10,
  col=c(rep('blue',3),rep('red',3)),
  names=rep(c(4,7,10),2),notch=T,
  las=1,ylim=c(0,1),yaxs='i',
  xlab='Average lifespan',
  ylab=expression('R'['0,wAlb']*' / R'['0,WT']))
legend('topright',legend=c('SF','DF'),
       fill=c('blue','red'),bty='n')

dev.off()

