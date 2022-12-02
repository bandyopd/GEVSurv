#################################################################################
## R code to run the h-likelihood estimation procedure for fitting the GEV model
## Written by Jeongseop Han; Edited by Dipankar Bandyopadhyay on 12/01/2022
#################################################################################


rm(list = ls())

setwd("C:\\Users\\dbandyop\\OneDrive - VCUHealth\\Work\\IL-DO-HA\\Code\\GEV_code_final\\GEV_code_final\\")  # put your working directory here

# Source various .R files 

#optimize v
source("hlik_vi_ftn.R")

#optimize beta
source("pvh_ftn.R")

#optimize phi
source("ppsih_ftn.R")

#se for beta
source("beta_se_ftn.R")


####################
#### TEETH data ####
####################


teeth = read.csv("TeethGEV.csv")

X_mat = model.matrix(~ Mobility + BOP + Plaque + PDmean + CALmean + Crown + Filled + Decayed + Dfsites + Age + Gender + Diabetes + Tobacco + Molar, data = teeth)

gp = factor(teeth$id)

n_vec = table(gp)
N=nrow(X_mat)
p=ncol(X_mat)
q=length(n_vec)
y = log(teeth$time)
delta = teeth$delta

tol = 1e-3

X = list()

X[[1]] = X_mat[1:n_vec[1], ]

for(i in 2:q)
{
  if(n_vec[i] > 1)
  {
    X[[i]] = X_mat[ (sum(n_vec[1:(i-1)]) + 1):( sum(n_vec[1:i]) ),  ]
  }else{
    X[[i]] = X_mat[(sum(n_vec[1:(i-1)]) + 1), ]
  }
}

y_list = list()

y_list[[1]] = y[1:n_vec[1]]

for(i in 2:q)
{
  y_list[[i]] = y[ (sum(n_vec[1:(i-1)]) + 1):( sum(n_vec[1:i]) ) ]
}

delta_list = list()

delta_list[[1]] = delta[1:n_vec[1]]

for(i in 2:q)
{
  delta_list[[i]] = delta[ (sum(n_vec[1:(i-1)]) + 1):( sum(n_vec[1:i]) ) ]
}

beta_temp = c(2, -0.618, -0.004, 0.001, -0.192, -0.277, 0.342, 0.796, -0.591, 0.262, -0.008, -0.069, -0.280, -1.148, -0.262)

sig_temp = 1.35

zeta_temp = -0.5

alpha_temp = 2

phi_temp = c(sig_temp, zeta_temp, alpha_temp)

theta_temp = c(beta_temp, phi_temp)

kk = 1

while(1)
{
  v_new = rep(0, q)
  
  for(i in 1:q)
  {
    res_v = optimize(f = hlik_vi_ftn, interval = c(-10, 10), maximum = TRUE, i=i, theta = theta_temp, dist = "GEV")
    
    v_new[i] = res_v$maximum
  }
  
  res_beta = optim(par = beta_temp, fn = pvh_ftn, phi = phi_temp, v = v_new, dist = "GEV")
  
  beta_new = res_beta$par
  
  psi_new = c(beta_new, v_new)
  
  res_phi = optim(par = phi_temp, fn = ppsih_ftn, psi = psi_new, dist = "GEV")
  
  phi_new = res_phi$par
  
  theta_new = c(beta_new, phi_new)
  
  err = max(abs(theta_new[-1] - theta_temp[-1]))
  
  if(round(err, 3) == 0)
  {
    theta_hat = theta_new
    v_hat = v_new
    break
  }else if(kk > 100)
  {
    theta_hat = theta_new
    v_hat = v_new
    break
  }else{
    theta_temp = theta_new
    print(round(theta_temp, 3))
    kk = kk + 1
  }
}

beta.hat = theta_hat[1:p]
sig.hat = theta_hat[p+1]
zeta.hat = theta_hat[p+2]
alpha.hat = theta_hat[p+3]

se.beta.hat = beta_se_ftn(theta_hat, v_hat, dist = "GEV")

beta.conf.low = beta.hat - 1.96*se.beta.hat
beta.conf.upp = beta.hat - 1.96*se.beta.hat

2*res_beta$value + 2*length(theta_hat)  #mAIC

source("ell1_i_ftn.R")

ell1 = 0

for(i in 1:q)
{
  ell1 = ell1 + ell1_i_ftn(v_hat[i], i, theta_hat, dist = "GEV")
}

source("df_ftn.R")

-2*ell1 + 2*df_ftn(theta_hat, v_hat, dist = "GEV")  #cAIC
