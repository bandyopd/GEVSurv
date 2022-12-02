hlik_vi_ftn = function(v_i, i, theta, dist)
{
  if(dist == "GEV")
  {
    p = length(theta) - 3
    
    beta = theta[1:p]
    
    sig = theta[p+1]
    
    zeta = theta[p+2]
    
    alpha = theta[p+3]
    
    hlik_i = 0
    
    for(j in 1:n_vec[i])
    {
      if(n_vec[i] > 1)
      {
        muij = X[[i]][j, ]%*%beta + v_i
        Mij = (1 + zeta*(y_list[[i]][j] - muij)/sig)^(-1/zeta)
        hlik_i = hlik_i + delta_list[[i]][j]*( -log(sig) + (1 + zeta)*log(Mij) - Mij) + (1 - delta_list[[i]][j])*log(1 - exp(-Mij))
      }else{
        muij = X[[i]]%*%beta + v_i
        Mij = (1 + zeta*(y_list[[i]][j] - muij)/sig)^(-1/zeta)
        hlik_i = hlik_i + delta_list[[i]][j]*( -log(sig) + (1 + zeta)*log(Mij) - Mij) + (1 - delta_list[[i]][j])*log(1 - exp(-Mij))
      }
    }
    
    hlik_i = hlik_i - 0.5*log(2*pi*alpha) - 0.5*(v_i^(2))/alpha
  }else if(dist == "N")
  {
    p = length(theta) - 2
    
    beta = theta[1:p]
    
    lambda = theta[p+1]
    
    alpha = theta[p+2]
    
    hlik_i = 0
    
    for(j in 1:n_vec[i])
    {
      if(n_vec[i] > 1)
      {
        muij = as.numeric( X[[i]][j, ]%*%beta + v_i)
        zij = as.numeric( (y_list[[i]][j] - muij)/sqrt(lambda) )
        Nij = dnorm(zij)
        Pij = pnorm(zij)
        hlik_i = hlik_i + delta_list[[i]][j]*(log(Nij) - log(sqrt(lambda))) + (1 - delta_list[[i]][j])*log(1 - Pij)
      }else{
        muij = as.numeric( X[[i]]%*%beta + v_i)
        zij = as.numeric( (y_list[[i]] - muij)/sqrt(lambda) )
        Nij = dnorm(zij)
        Pij = pnorm(zij)
        hlik_i = hlik_i + delta_list[[i]]*(log(Nij) - log(sqrt(lambda))) + (1 - delta_list[[i]])*log(1 - Pij)
      }
    }
    
    hlik_i = hlik_i - 0.5*log(2*pi*alpha) - 0.5*(v_i^(2))/alpha
  }
  
  return(hlik_i)
  
}
