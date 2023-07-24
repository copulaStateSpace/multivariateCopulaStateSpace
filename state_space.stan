functions
{

  real bicop_vec_wt_lpdf(vector z1, vector z2, real tau, int family,int T)
  {
    real ll=0;
    real t1;
    real t2;
    real theta;
    real u_rot;
    real v_rot;
    real ll1=0;
    real ll2=0;
    real ll3=0;
    real ll13=0;
    real ll4=0;
    real normalizer;
    real u[T];
    real v[T];
    for(t in 1:T)
    {
      u[t]=student_t_cdf(z1[t],4, 0, 1);
      v[t]=student_t_cdf(z2[t],4,0,1);
    }
    if(family == 1 || family==100)
    {
      theta = sin(tau*pi()/2);
      for(t in 1:T)
      {
        t1 = inv_Phi(u[t]); 
        t2 = inv_Phi(v[t]);
        ll1 += log(1.0/sqrt(1.0-pow(theta,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*theta*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(theta,2.0)))));
      }
      ll=ll1;
    }
    
    
    if(family == 2 || family==100)
    {
      theta = sin(tau*pi()/2);
      for(t in 1:T)
      {
        ll2 += -log(2*pi())-1.0/2.0*log(1-pow(theta,2)) - student_t_lpdf(z1[t]|4,0,1) - student_t_lpdf(z2[t]|4,0,1) - (4.0+2.0)/2.0 * log(1+( pow(z1[t],2) + pow(z2[t],2) - 2.0 * theta *z1[t]*z2[t])/(4.0*(1-pow(theta,2))));
      }
      ll=ll2;
    }
    
    
    
    if(family == 3 || family==100)
    {
      if(tau>=0)
      {
        theta = 2 * tau/(1 - tau);
        for(t in 1:T)
        {
          ll3+=log1p(theta)-(1.0+theta)*log(u[t]*v[t])-(2.0+1.0/(theta))*log(pow(u[t],-theta)+pow(v[t],-theta)-1.0);
        }
      }
      if(tau<0)
      {
        theta = -(2 * tau)/(1 + tau);
        for(t in 1:T)
        {
          u_rot=1-u[t];
          ll3+=log1p(theta)-(1.0+theta)*log(u_rot*v[t])-(2.0+1.0/(theta))*log(pow(u_rot,-theta)+pow(v[t],-theta)-1.0);
        }
      }
      ll=ll3;
    }
    
    
    if(family == 13|| family==100)
    {
      if(tau>=0)
      {
        theta = 2 * tau/(1 - tau);
        for(t in 1:T)
        {
          u_rot = 1 - u[t];
          v_rot = 1 - v[t];
          ll13+=log1p(theta)-(1.0+theta)*log(u_rot*v_rot)-(2.0+1.0/(theta))*log(pow(u_rot,-theta)+pow(v_rot,-theta)-1.0);
        }
      }
      if(tau<0)
      {
        theta = -(2 * tau)/(1 + tau);
        for(t in 1:T)
        {
          u_rot=u[t];
          v_rot = 1 - v[t];
          ll13+=log1p(theta)-(1.0+theta)*log(u_rot*v_rot)-(2.0+1.0/(theta))*log(pow(u_rot,-theta)+pow(v_rot,-theta)-1.0);
        }
      }
      ll=ll13;
    }
    
    
    
    
    if(family == 4|| family==100)
    {
      if(tau>=0)
      {
        theta=1/(1-tau);
        for(t in 1:T)
        {
          t1 = pow(-log(u[t]),theta)+pow(-log(v[t]),theta);
          ll4+= -pow(t1,1.0/(theta))+(2.0/(theta)-2.0)*log(t1)+(theta-1.0)*log(log(u[t])*log(v[t]))-log(u[t]*v[t])+log1p((theta-1.0)*pow(t1,-1.0/(theta)));
        }
      }
      if(tau<0)
      {
        theta=1/(1+tau);
        for(t in 1:T)
        {
          u_rot=1-u[t];
          t1 = pow(-log(u_rot),theta)+pow(-log(v[t]),theta);
          ll4+= -pow(t1,1.0/(theta))+(2.0/(theta)-2.0)*log(t1)+(theta-1.0)*log(log(u_rot)*log(v[t]))-log(u_rot*v[t])+log1p((theta-1.0)*pow(t1,-1.0/(theta)));
        }
      }
      ll=ll4;
    }
    

    if(family==100)
    {
      normalizer=ll4;
      ll=log((exp(ll1-normalizer)+exp(ll2-normalizer)+exp(ll3-normalizer)+exp(ll4-normalizer))/4.0)+normalizer;
      
      return ll;
    }

    return ll;
  }
 }


data
{
  int T;
  int d;
  int lengths[d];
  matrix[T,d] U;
  int fam_obs;
  int fam_lat;
  int indices1[T,d];
  int indices2[T,d];
}

parameters
{
  real<lower=0.3,upper=0.95> tau_obs1;
  vector<lower=-0.95,upper=0.95>[d-1] tau_obs;
  real<lower=-0.95,upper=0.95> tau_lat;
  vector[T+1] v_cont;
}

model
{
  tau_obs1 ~ beta(10,1.5);
  v_cont[1]~student_t(4,0,1);
  target += bicop_vec_wt_lpdf(v_cont[2:(T+1)]|v_cont[1:T],tau_lat,fam_lat,T)+student_t_lpdf(v_cont[2:(T+1)]|4,0,1);
  U[indices1[1:lengths[1],1],1] ~bicop_vec_wt(v_cont[indices2[1:lengths[1],1]],tau_obs1,fam_obs,lengths[1]);
  for(j in 2:d)
  { 
    U[indices1[1:lengths[j],j],j] ~bicop_vec_wt(v_cont[indices2[1:lengths[j],j]],tau_obs[j-1],fam_obs, lengths[j]);
  }
} 
