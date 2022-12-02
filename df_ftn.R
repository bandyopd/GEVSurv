df_ftn = function(theta, v, dist)
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
    
    Hvv2 = Hvv + 1/alpha
    
    H11 = cbind(Hbetabeta, Hbetav)
    H12 = cbind(t(Hbetav), diag(Hvv))
    
    H122 = cbind(t(Hbetav), diag(Hvv2))
    
    H1 = rbind(H11, H12)  #H_ast in paper
    H2 = rbind(H11, H122) #H in paper
    
    df = sum(diag(solve(H2, H1)))
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
    
    Hvv2 = Hvv + 1/alpha
    
    H11 = cbind(Hbetabeta, Hbetav)
    H12 = cbind(t(Hbetav), diag(Hvv))
    
    H122 = cbind(t(Hbetav), diag(Hvv2))
    
    H1 = rbind(H11, H12)  #H_ast in paper
    H2 = rbind(H11, H122) #H in paper
    
    # HH = solve(H2)
    # 
    # HHH = HH%*%H1
    df = sum(diag(solve(H2, H1)))
  }
  
  return(df)
}

