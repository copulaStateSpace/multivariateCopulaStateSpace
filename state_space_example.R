
library(rstan)
library(VineCopula)

#function to estimate the mode
estimate_mode=function(x) {
  tryCatch( { d=density(x)
  return(d$x[which.max(d$y)])
  }, error = function(e) {print(e)
    return(mean(x))})
}

#simulate from the factor latent Dvine model
state_factor_cop_sim = function(T,fam_obs,tau_obs,par2_obs,fam_lat,tau_lat,par2_lat){
  d=length(tau_obs)
  U=matrix(nrow=T, ncol=d)
  if(length(fam_obs)==1)
  {
    fam_obs=rep(fam_obs,d)
  }
  par_obs=BiCopTau2Par(family=fam_obs,tau=tau_obs)
  par_lat=BiCopTau2Par(family=fam_lat,tau=tau_lat)
  v=numeric(T)
  v[1]=runif(1)
  for(t in 2:T)
  {
    v[t]=BiCopCondSim(1,cond.val=v[t-1],cond.var=2,family=fam_lat,par=par_lat,par2=par2_lat)
  }
  for(j in 1:d)
  {
    U[,j]=BiCopCondSim(T,cond.val=v,cond.var=2,family=fam_obs[j],par=par_obs[j],par2=par2_obs[j])
  }
  return(list(U=U,T=T,v.true=v))
}

#build stan model
stan_m = stan_model(file='state_space.stan')


sample_faccopss=function(U, iter=200, burnin=100, fam_lat=100, fam_obs=100, chains=1, cores=1, para="C", lags=1, get_U_sim=TRUE, adapt_delta=0.8)
{
  T=dim(U)[1]
  d=dim(U)[2]
  n_lags=length(lags)
  Z_student=qt(U, df=4) #STAN has no t quantile function, to allow for a student t copula with 4 df we transform the udata
  lengths=numeric(d) #stores how many non NA values each dimension has
  indices1=matrix(0,nrow=T, ncol=d) #stores the non NA indices for each dimension
  for(j in 1:d)
  {
    no_NA=which(!is.na(U[,j]))
    lengths[j]=length(no_NA)
    indices1[1:lengths[j], j]=no_NA
  }
  indices2=indices1+1 # non NA indices shifted by one, needed inside the STAN program
  Z_stan=Z_student
  Z_stan[is.na(Z_stan)]=1 #STAN doesnt accept NA values, just set to any values, they are not used inside the STAN program
  data_stan=list(T=T, d=d, lengths=lengths, U=Z_stan, fam_obs=fam_obs,fam_lat=fam_lat, indices1=indices1, indices2=indices2)
  init_list=list(list(tau_obs1=0.6, tau_obs=as.array(rep(0.1,d-1)), tau_lat=0.7, v_cont=rep(0.5,T+1)))
  
  fit=sampling(stan_m,iter=iter,warmup=burnin,chains=chains,cores=cores,data=data_stan, init=init_list, control=list(adapt_delta=adapt_delta))

  mcmc_iters=extract(fit, permute=FALSE)
  tau_array=mcmc_iters[,,1:d]
  tau_lat_array=mcmc_iters[,,d+1]
  v_array=pt(mcmc_iters[,,(d+3):(d+3+T-1)], df=4)
  len=iter-burnin
  U_sim=array(dim=c(T,d,len))
  f_samples=matrix(nrow=len, ncol=d+1)
  m1=matrix(nrow=len, ncol=4)

  #sample linking copula families
  for(j in 1:d)
  {
    no_NA=which(!is.na(U[,j]))
    for(i in 1:len)
    {
      if(fam_obs==100)
      {
        if(tau_array[i,j]>0)
        {
          m1[i,]=cbind(sum(log(BiCopPDF(v_array[i,no_NA],U[no_NA,j], par=(BiCopTau2Par(family=1, tau=tau_array[i,j])), family=1))),
                       sum(log(BiCopPDF(v_array[i,no_NA],U[no_NA,j], par=(BiCopTau2Par(family=2, tau=tau_array[i,j])), family=2, par2=4))), 
                       sum(log(BiCopPDF(v_array[i,no_NA],U[no_NA,j], par=(BiCopTau2Par(family=3, tau=min(tau_array[i,j], 0.95))), family=3))),
                       sum(log(BiCopPDF(v_array[i,no_NA],U[no_NA,j], par=(BiCopTau2Par(family=4, tau=min(tau_array[i,j], 0.95))), family=4))))
        }else{
          m1[i,]=cbind(sum(log(BiCopPDF(v_array[i,no_NA],U[no_NA,j], par=(BiCopTau2Par(family=1, tau=tau_array[i,j])), family=1))),
                       sum(log(BiCopPDF(v_array[i,no_NA],U[no_NA,j], par=(BiCopTau2Par(family=2, tau=tau_array[i,j])), family=2, par2=4))),
                       sum(log(BiCopPDF(v_array[i,no_NA],U[no_NA,j], par=(BiCopTau2Par(family=23, tau=max(tau_array[i,j], -0.95))), family=23))),
                       sum(log(BiCopPDF(v_array[i,no_NA],U[no_NA,j], par=(BiCopTau2Par(family=24, tau=max(tau_array[i,j], -0.95))), family=24))))
        }
        f_samples[i,j]=sample(x=1:4, size=1, prob=exp(m1[i,]-max(m1[i,])))
      }else{
        f_samples[i,j]=fam_obs
      }
    }
  }
  
  #simulate from the predictive distribution
  fam_vec=numeric(len)
  if(get_U_sim)
  {
    for(j in 1:d)
    {
      for(t in 1:T)
      {
        fam_vec=f_samples[,j]
        ind=which((tau_array[,j]<0) & (f_samples[,j] %in% c(3,4)))
        fam_vec[ind]=fam_vec[ind]+20
        U_sim[t,j,]=BiCopCondSim(N=len, cond.val= v_array[,t], cond.var=2, family=fam_vec, par= BiCopTau2Par(tau=pmax(pmin(tau_array[,j],0.98),-0.98), family=fam_vec), par2=4)
      }
    }
  }
  
  #sample latent copula family
  m2=matrix(nrow=len, ncol=4)
  for(i in 1:len)
  {
    if(fam_lat==100)
    {
      if(tau_lat_array[i]>=0)
      {
        m2[i,]=cbind(sum(log(BiCopPDF(v_array[i,2:T],v_array[i,1:(T-1)], par=(BiCopTau2Par(family=1, tau=tau_lat_array[i])), family=1))),
                     sum(log(BiCopPDF(v_array[i,2:T],v_array[i,1:(T-1)], par=(BiCopTau2Par(family=2, tau=tau_lat_array[i])), family=2, par2=4))), 
                     sum(log(BiCopPDF(v_array[i,2:T],v_array[i,1:(T-1)], par=(BiCopTau2Par(family=3, tau=min(tau_lat_array[i], 0.95))), family=3))),
                     sum(log(BiCopPDF(v_array[i,2:T],v_array[i,1:(T-1)], par=(BiCopTau2Par(family=4, tau=min(tau_lat_array[i], 0.95))), family=4))))
        
      }else{
        m2[i,]=cbind(sum(log(BiCopPDF(v_array[i,2:T],v_array[i,1:(T-1)], par=(BiCopTau2Par(family=1, tau=tau_lat_array[i])), family=1))), 
                     sum(log(BiCopPDF(v_array[i,2:T],v_array[i,1:(T-1)], par=(BiCopTau2Par(family=2, tau=tau_lat_array[i])), family=2, par2=4))),
                     sum(log(BiCopPDF(v_array[i,2:T],v_array[i,1:(T-1)], par=(BiCopTau2Par(family=23, tau=max(tau_lat_array[i], -0.95))), family=23))),
                     sum(log(BiCopPDF(v_array[i,2:T],v_array[i,1:(T-1)], par=(BiCopTau2Par(family=24, tau=max(tau_lat_array[i], -0.95))), family=24))))
        
      }
      f_samples[i,d+1]=sample(x=1:4, size=1, prob=exp(m2[i,]-max(m2[i,])))
    }else{
      f_samples[i,d+1]=fam_lat
    }
  }
  v_mode=apply(v_array, 2, estimate_mode)
  v_q5=apply(v_array, 2, quantile, prob=0.05)
  v_q95=apply(v_array, 2, quantile, prob=0.95)
  tau_mode=c(apply(tau_array, 2, estimate_mode), estimate_mode(tau_lat_array))
  tau_q5=c(apply(tau_array, 2, quantile, prob=0.05), quantile(tau_lat_array, prob=0.05))
  tau_q95=c(apply(tau_array, 2, quantile, prob=0.95), quantile(tau_lat_array, prob=0.95))
  
  #select families as marginal posterior mode estimates
  f_sel=numeric(d+1)
  for(j in 1:(d+1))
  {
    f_sel[j]=as.numeric(names(sort(table(f_samples[,j]),decreasing=TRUE)[1]))
  }
  
  return(list(fit=fit,v_mode=v_mode, v_q5=v_q5, v_q95=v_q95, tau_mode=tau_mode,tau_q5=tau_q5, tau_q95=tau_q95,f_samples=f_samples,f_sel=f_sel, U_sim=U_sim))
}


d=8

set.seed(123)
sim=state_factor_cop_sim(T=500,fam_obs=c(1,1,2,2,3,3,4,4),tau_obs=rep(0.6,d),par2_obs=rep(4,d),fam_lat=1,tau_lat=0.7,par2_lat=4)
fit=sample_faccopss(sim$U, iter=500, burnin=250, fam_obs=100,fam_lat=100, para="NC")

fit$f_sel
fit$tau_mode

traceplot(fit$fit)

#plot estimated states on the Z-scale
plot(qnorm(fit$v_mode), type="l")
lines(qnorm(sim$v.true), col="red", lty=2)
lines(qnorm(fit$v_q5), col="grey")
lines(qnorm(fit$v_q95), col="grey")
