beta_se_ftn = function(theta, v, dist)
{
  if(dist == "GEV")
  {
    beta = theta[1:p]
    
    sig = theta[p+1]
    
    zeta = theta[p+2]
    
    alpha = theta[p+3]
    
    Hbetabeta = matrix(0, nrow = p, ncol = p)
    
    Hbetav = matrix(0, nrow = p, ncol = q)
    
    Hvv = rep(0, q)
    
    for(i in 1:q)
    {
      if(n_vec[i] > 1)
      {
        for(j in 1:n_vec[i])
        {
          muij = X[[i]][j, ]%*%beta + v[i]
          
          Mij = (1 + zeta*(y_list[[i]][j] - muij)/sig)^(-1/zeta)
          
          Hbetabeta = Hbetabeta - as.numeric( delta_list[[i]][j]*(1 + zeta)*( sig^(-2) )*( zeta*Mij^(2*zeta) - Mij^(1 + 2*zeta) ) )*X[[i]][j, ]%*%t(X[[i]][j, ]) - as.numeric( (1 - delta_list[[i]][j])*( sig^(-2) )*( (1+zeta)*(Mij^(1 + 2*zeta))/(exp(Mij) - 1) - ( Mij^(2 + 2*zeta) )*exp(Mij)/( (exp(Mij) - 1)^2 ) ) )*X[[i]][j, ]%*%t(X[[i]][j, ])
          
          Hbetav[, i] = Hbetav[, i] - as.numeric( delta_list[[i]][j]*(1 + zeta)*( sig^(-2) )*( zeta*Mij^(2*zeta) - Mij^(1 + 2*zeta) ) )*X[[i]][j, ]  - as.numeric( (1 - delta_list[[i]][j])*( sig^(-2) )*( (1+zeta)*(Mij^(1 + 2*zeta))/(exp(Mij) - 1) - ( Mij^(2 + 2*zeta) )*exp(Mij)/( (exp(Mij) - 1)^2 ) ) )*X[[i]][j, ]
          
          Hvv[i] = Hvv[i] - delta_list[[i]][j]*(1 + zeta)*( sig^(-2) )*( Mij^(1 + zeta) )*( zeta*Mij^(zeta-1) - Mij^zeta ) - (1 - delta_list[[i]][j])*( sig^(-2) )*( Mij^(1 + zeta) )*( (1+zeta)*(Mij^(zeta))/(exp(Mij) - 1) - ( Mij^(1 + zeta) )*exp(Mij)/( (exp(Mij) - 1)^2 ) )
        }
      }else{
        muij = X[[i]]%*%beta + v[i]
        
        Mij = (1 + zeta*(y_list[[i]] - muij)/sig)^(-1/zeta)
        
        Hbetabeta = - as.numeric( delta_list[[i]]*(1 + zeta)*( sig^(-2) )*( zeta*Mij^(2*zeta) - Mij^(1 + 2*zeta) ) )*X[[i]]%*%t(X[[i]]) - as.numeric( (1 - delta_list[[i]])*( sig^(-2) )*( (1+zeta)*(Mij^(1 + 2*zeta))/(exp(Mij) - 1) - ( Mij^(2 + 2*zeta) )*exp(Mij)/( (exp(Mij) - 1)^2 ) ) )*X[[i]]%*%t(X[[i]])
        
        Hbetav[, i] = - as.numeric( delta_list[[i]]*(1 + zeta)*( sig^(-2) )*( zeta*Mij^(2*zeta) - Mij^(1 + 2*zeta) ) )*X[[i]]  - as.numeric( (1 - delta_list[[i]])*( sig^(-2) )*( (1+zeta)*(Mij^(1 + 2*zeta))/(exp(Mij) - 1) - ( Mij^(2 + 2*zeta) )*exp(Mij)/( (exp(Mij) - 1)^2 ) ) )*X[[i]]
        
        Hvv[i] = - delta_list[[i]]*(1 + zeta)*( sig^(-2) )*( Mij^(1 + zeta) )*( zeta*Mij^(zeta-1) - Mij^zeta ) - (1 - delta_list[[i]])*( sig^(-2) )*( Mij^(1 + zeta) )*( (1+zeta)*(Mij^(zeta))/(exp(Mij) - 1) - ( Mij^(1 + zeta) )*exp(Mij)/( (exp(Mij) - 1)^2 ) )
      }
    }
    
    Hvv = Hvv + 1/alpha
    
    beta_se = sqrt(diag( solve( Hbetabeta - Hbetav%*%diag(1/Hvv)%*%t(Hbetav) ) ) )
  }else if(dist == "N")
  {
    beta = theta[1:p]
    
    lambda = theta[p+1]
    alpha = theta[p+2]
    
    Hbetabeta = matrix(0, nrow = p, ncol = p)
    
    Hbetav = matrix(0, nrow = p, ncol = q)
    
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
          
          Hbetabeta = Hbetabeta + as.numeric(delta_list[[i]][j]/sqrt(lambda))*X[[i]][j, ]%*%t(X[[i]][j, ]) - ( ( ( 1-delta_list[[i]][j])/lambda )/(1 - Pij)^2 )*( zij*Nij*(1-Pij) - (Nij^2) )*X[[i]][j, ]%*%t(X[[i]][j, ])
          
          Hbetav[, i] = Hbetav[, i] + delta_list[[i]][j]*X[[i]][j, ]/sqrt(lambda) - ( ( ( 1-delta_list[[i]][j])/lambda )/(1 - Pij)^2 )*( zij*Nij*(1-Pij) - (Nij^2) )*X[[i]][j, ]
          
          Hvv[i] = Hvv[i] + delta_list[[i]][j]/sqrt(lambda) - ( ( (1-delta_list[[i]][j])/lambda )/(1 - Pij)^2 )*( zij*Nij*(1-Pij) - (Nij^2) )
          
        }else{
          muij = as.numeric( X[[i]]%*%beta + v[i] )
          
          zij = as.numeric( (y_list[[i]] - muij)/sqrt(lambda) )
          
          Nij = dnorm(zij)
          
          Pij = pnorm(zij)
          
          Hbetabeta = Hbetabeta + delta_list[[i]]*X[[i]]%*%t(X[[i]])/sqrt(lambda) - ( ( ( 1-delta_list[[i]])/lambda )/(1 - Pij)^2 )*( zij*Nij*(1-Pij) - (Nij^2) )*X[[i]]%*%t(X[[i]])
          
          Hbetav[, i] = Hbetav[, i] + delta_list[[i]]*X[[i]]/sqrt(lambda) - ( ( ( 1-delta_list[[i]])/lambda )/(1 - Pij)^2 )*( zij*Nij*(1-Pij) - (Nij^2) )*X[[i]]
          
          Hvv[i] = Hvv[i] + delta_list[[i]]/sqrt(lambda) - ( ( (1-delta_list[[i]])/lambda )/(1 - Pij)^2 )*( zij*Nij*(1-Pij) - (Nij^2) )
        }
      }
    }
    
    Hvv = Hvv + 1/alpha
    
    beta_se = sqrt(diag( solve( Hbetabeta - Hbetav%*%diag(1/Hvv)%*%t(Hbetav) ) ) )
  }
  
  return(beta_se)
}

