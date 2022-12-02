pvh_ftn = function(beta, phi, v, dist)
{
  source("hlik_vi_ftn.R")
  
  if(dist == "GEV")
  {
    hlik = 0
    
    theta = c(beta, phi)
    
    sig = phi[1]
    zeta = phi[2]
    alpha = phi[3]
    
    for(i in 1:q)
    {
      hlik = hlik + hlik_vi_ftn(v[i], i, theta, dist)
    }
    
    Hvv = rep(0, q)
    
    for(i in 1:q)
    {
      for(j in 1:n_vec[i])
      {
        if(n_vec[i] > 1)
        {
          muij = X[[i]][j, ]%*%beta + v[i]
          
          Mij = (1 + zeta*(y_list[[i]][j] - muij)/sig)^(-1/zeta)
          
          Hvv[i] = Hvv[i] - delta_list[[i]][j]*(1 + zeta)*( sig^(-2) )*( Mij^(1 + zeta) )*( zeta*Mij^(zeta-1) - Mij^zeta ) - (1 - delta_list[[i]][j])*( sig^(-2) )*( Mij^(1 + zeta) )*( (1+zeta)*(Mij^(zeta))/(exp(Mij) - 1) - ( Mij^(1 + zeta) )*exp(Mij)/( (exp(Mij) - 1)^2 ) )
        }else{
          muij = X[[i]]%*%beta + v[i]
          
          Mij = (1 + zeta*(y_list[[i]][j] - muij)/sig)^(-1/zeta)
          
          Hvv[i] = Hvv[i] - delta_list[[i]][j]*(1 + zeta)*( sig^(-2) )*( Mij^(1 + zeta) )*( zeta*Mij^(zeta-1) - Mij^zeta ) - (1 - delta_list[[i]][j])*( sig^(-2) )*( Mij^(1 + zeta) )*( (1+zeta)*(Mij^(zeta))/(exp(Mij) - 1) - ( Mij^(1 + zeta) )*exp(Mij)/( (exp(Mij) - 1)^2 ) )
        }
      }
    }
    
    Hvv = Hvv + 1/alpha
    
    pvh = hlik - 0.5*sum(log(Hvv))
  }else if(dist == "N")
  {
    hlik = 0
    
    theta = c(beta, phi)
    lambda = phi[1]
    alpha = phi[2]
    
    for(i in 1:q)
    {
      hlik = hlik + hlik_vi_ftn(v[i], i, theta, dist)
    }
    
    Hvv = rep(0, q)
    
    for(i in 1:q)
    {
      for(j in 1:n_vec[i])
      {
        if(n_vec[i] > 1)
        {
          muij = as.numeric( X[[i]][j, ]%*%beta + v[i] )
          
          zij = as.numeric( (y_list[[i]][j] - muij)/sqrt(lambda) )
          
          Nij = dnorm(zij)
          
          Pij = pnorm(zij)
          
          Hvv[i] = Hvv[i] + delta_list[[i]][j]/sqrt(lambda) - ( ( (1-delta_list[[i]][j])/lambda )/(1 - Pij)^2 )*( zij*Nij*(1-Pij) - (Nij^2) )
          
        }else{
          muij = as.numeric( X[[i]]%*%beta + v[i] )
          
          zij = as.numeric( (y_list[[i]] - muij)/sqrt(lambda) )
          
          Nij = dnorm(zij)
          
          Pij = pnorm(zij)
          
          Hvv[i] = Hvv[i] + delta_list[[i]]/sqrt(lambda) - ( ( (1-delta_list[[i]])/lambda )/(1 - Pij)^2 )*( zij*Nij*(1-Pij) - (Nij^2) )
        }
      }
    }
    
    Hvv = Hvv + 1/alpha
    
    pvh = hlik - 0.5*sum(log(Hvv))
  }
  

  return(-pvh)
}