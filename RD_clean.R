library(bmsR)
library(qpcR)

aces = read.csv("~/Dropbox (University of Oregon)/ACEs/ACE_trial_data.csv")

ace_sub = read.csv("~/Dropbox (University of Oregon)/ACEs/ACE_dataset.csv")
S2 = unique(ace_sub$Participant.ID)
aces = aces[aces$Participant.ID%in%ace_sub$Participant.ID,]
S = unique(aces$Participant.ID)
aces$N=240
aces$LL=aces$leave_TRUE
aces$s2_u=0
aces$X1=aces$X2=aces$T1=aces$T2=0
aces$Apple_return=as.numeric(aces$Apple_return)
aces$Lag1=as.numeric(aces$Lag1)
aces$Lag2=as.numeric(aces$Lag2)

safelog = function(x){
  x[x==0]=1e-100
  y=log(x)
  return(y)
}
mvt = function(pars,dat){
  alpha = pars[1]; beta = pars[2]
  dat$Apple_return[is.na(dat$Apple_return)]=1e-6
  fg=bg=c();fg[1]=bg[1]=dat$Apple_return[1]
  cost = dat$choice_time
  for(k in 2:nrow(dat)){
    bg[k] = ((1-alpha)^cost[k])*(dat$Apple_return[k]/cost[k])+(1-(1-alpha)^cost[k])*bg[k-1]
    fg[k] = mean(c(dat$Lag1[k],dat$Lag2[k]),na.rm=T)/
      mean(c(dat$choice_time[k],dat$choice_time[k-1]),na.rm=T)
  }
  P =pnorm(beta*(bg-fg))
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
R1 = function(pars,dat){
  beta = pars[1]; alpha = pars[2]; lapse = 0
  r = cbind(dat$X2,dat$X1)
  r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/abs(r)*beta
  k = s2_e/dat$s2_u
  V = r/(1+(k*t))
  P = (1-lapse)*pnorm(alpha*(V[,1]-V[,2]))+lapse*0.5
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}
r1=MVT=list()
ntrials = c()
for(s in 1:length(S)){
  ix = aces[aces$Participant.ID==S[s],]
  ix = ix[!is.na(ix$Response),]
  ix$N = nrow(ix)
  ntrials[s] = nrow(ix)
  ix$T1=3;ix$T2=ifelse(ix$Environment==1,6,12)
  for(n in 1:nrow(ix)){
    ix$X1[n] = mean(c(ix$Lag1[n],ix$Lag2[n]),na.rm=T)
    ix$X2[n] = 10
    if(n==1) ix$X1[n] = 10
    ix$s2_u[n] = max(1e-5,var(c(ix$X1[1:n],ix$X2[1:n]),na.rm=T))
  }
  r1[[s]] = nlminb(start = rtnorm(2,0,1,0,1), objective = R1, lower = c(1e-6,1e-6), upper = c(50,100), dat = ix)
  MVT[[s]] = nlminb(start = rtnorm(2,0,1,0,1), objective = mvt, lower = rep(1e-6,2), upper = c(1,50), dat = ix)
  print(s)
}

bic = matrix(0,length(r1),2)
pars=matrix(0,length(r1),2)
for(i in 1:length(r1)){
  bic[i,1] = 2*r1[[i]]$objective+log(ntrials[i])*2
  bic[i,2] = 2*MVT[[i]]$objective+log(ntrials[i])*2
  
  pars[i,1] = r1[[i]]$par[1]
  pars[i,2] = r1[[i]]$par[2]
}

m = akaike.weights(bic)$weights
m = matrix(m,length(r1),2,byrow=F)
m[m==0]=1e-250
pxp = VB_bms(log(m))
pxp
colMeans(bic)


df = data.frame(beta = pars[,1],#beta
                alpha = pars[,2],#alpha
                subs = S[1:nrow(m)])
df = df[df$subs%in%S2,]
df = df[match(S2,df$subs),]
df$richscore = ace_sub$tot_score_rich
df$poorscore = ace_sub$tot_score_poor
df$acescore = ace_sub$ACE_Score
df$acegroup = as.factor(ace_sub$ACE_Group)
df$risk = ace_sub$Risk_all
a1 = aov(richscore~acescore*beta+acescore*alpha,df)
summary(a1)
cor.test(df$richscore,df$beta);cor.test(df$richscore,df$alpha)

a2 = aov(poorscore~acescore*beta+acescore*alpha,df)
summary(a2)
cor.test(df$poorscore,df$beta);cor.test(df$poorscore,df$alpha)

cor.test(df$acescore,df$beta)
cor.test(df$acescore,df$alpha)

a1 = aov(poorscore~acegroup*beta+acegroup*alpha,df)
summary(a1)


################################################################# PLOTS#################################################################
################################################################# PLOTS#################################################################
################################################################# PLOTS#################################################################
################################################################# PLOTS#################################################################
################################################################# PLOTS#################################################################
################################################################# PLOTS#################################################################
################################################################# PLOTS#################################################################
################################################################# PLOTS#################################################################
################################################################# PLOTS#################################################################
bic_df = data.frame(BIC=c(bic[,1],bic[,2]),
                    sub = rep(1:145,2),
                    model = rep(c("RI","MVT"),each=nrow(bic)))
ggplot(bic_df,aes(x=model,y=BIC,fill=model))+
  geom_violin(width=0.5)+
  geom_boxplot(color="black",width=0.15,alpha=0.2,position=position_dodge(width=0.5))+
  geom_point(position=position_jitterdodge(.1,jitter.height = .1),color="black",shape=21,size=2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c("#4AB4E0","#E18B2B"),labels=c("exemplar","prototype"))+
  ylab("BIC")+
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=25),
        axis.title=element_text(size=30,face="bold"))+
  theme(legend.title=element_text(size=30),
        legend.text = element_text(size=30))+
  theme(legend.position = "none")

pxp_df = data.frame(PXP=pxp$pxp,
                    model = c("RI","MVT"))
ggplot(pxp_df,aes(x=model,y=PXP,fill=model))+
  geom_bar(stat="identity")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c("#4AB4E0","#E18B2B"),labels=c("exemplar","prototype"))+
  ylab("PXP")+
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=25),
        axis.title=element_text(size=30,face="bold"))+
  theme(legend.title=element_text(size=30),
        legend.text = element_text(size=30))+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,1.1),expand= c(0,0))

pars_df = data.frame(beta = c(pars[,1],pars[,1]),
                     alpha = c(pars[,2],pars[,2]),
                     scores = c(df$richscore,df$poorscore),
                     environment = rep(c("rich","poor"),each=nrow(df)),
                     sub = rep(1:145,2))

ggplot()+
  geom_smooth(method="lm",aes(x=df$beta,y=df$richscore),col="#17940C",size=2)+
  geom_smooth(method="lm",aes(x=df$beta,y=df$poorscore),col="#E01292",size=2)+
  geom_point(aes(x=df$beta,y=df$richscore),shape=21,fill="#17940C",col="black",size=3)+
  geom_point(aes(x=df$beta,y=df$poorscore),shape=21,fill="#E01292",col="black",size=3)+
  ylab("Apples Accumulated")+xlab("Beta")+
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=25),
        axis.title=element_text(size=30,face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot(df)+
  geom_smooth(method="lm",aes(x=df$alpha,y=df$richscore),col="#17940C",size=2)+
  geom_smooth(method="lm",aes(x=df$alpha,y=df$poorscore),col="#E01292",size=2)+
  geom_point(aes(x=df$alpha,y=df$richscore),shape=21,fill="#17940C",col="black",size=3)+
  geom_point(aes(x=df$alpha,y=df$poorscore),shape=21,fill="#E01292",col="black",size=3)+
  ylab("Apples Accumulated")+xlab("Alpha")+
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=25),
        axis.title=element_text(size=30,face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggplot(df)+
  geom_smooth(method="lm",aes(x=df$beta,y=df$acescore),col="#E09577",size=3)+
  geom_smooth(method="lm",aes(x=df$alpha,y=df$acescore),col="#4B9488",size=3)+
  ylab("ACEs")+xlab("parameter estimate")+
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=25),
        axis.title=element_text(size=30,face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_jitter(aes(x=df$beta,y=df$acescore),fill="#E09577",size=3,shape=21,col="black")+
  geom_jitter(aes(x=df$alpha,y=df$acescore),fill="#4B9488",size=3,shape=21,col="black")

dfpars = data.frame(par.est = c(df$beta,df$alpha),
                    sub = rep(1:145,2),
                    ACEs = rep(df$acegroup,2),
                    par = rep(c("beta","alpha"),each=nrow(df)))
ggplot(dfpars,aes(x=ACEs,y=par.est,fill=ACEs))+
  geom_bar(stat="summary",col="black")+
  geom_errorbar(stat="summary",width=0.2)+
  scale_fill_manual(values=c("white","gray"))+
  facet_wrap(~par)+
  ylab("Parameter Estimate")+xlab("ACEs Group (0=Low#; 1=High#)")+
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=25),
        axis.title=element_text(size=30,face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_continuous(limits = c(0,50),expand= c(0,0))+
  theme(strip.text.x = element_text(
      size = 30,face = "bold.italic"
    ),strip.text.y = element_text(
      size = 30,face = "bold.italic"
    )
  )+theme(legend.position = "none")
