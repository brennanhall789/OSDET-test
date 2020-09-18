library(foreach)
library(doParallel)
library(beepr)
library(plotrix)

##### Simulate #####
func_file = 'osdet.R'
source(func_file)

### get multiple replicants of M monte carlo simulations of events { S_n > n*beta }
#### for increasing number of N


#### Gaussian ##### 

## simulate for varying values of n.
set.seed(232)
M = 1e6
N = seq(from=1, to = 40, by = 1)
cores=12
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  system.time(osdet <- foreach(k = N) %dopar% { 
  	replicate(M, gauss_osdet(n = k, b = 1, l = 50, mu = 0, sigma2 = 1))
  })
  beep(5)
  stopCluster(cl)


###############


avg_osdet <- array(as.numeric(unlist(osdet)), dim=c(M, length(N)))
  avg_osdet <- apply(avg_osdet, 2, sum)/M
  logAvg <- log(avg_osdet)
m2_osdet <- array(as.double(unlist(osdet)), dim=c(M, length(N)))
  m2_osdet <- apply(m2_osdet, 2, function(x) sum(x^2)/M)
  var_osdet <- m2_osdet - avg_osdet^2
  
logCI = data.frame(low = logAvg - 1.645*sqrt((var_osdet)/N), 
                   upr = logAvg + 1.645*sqrt((var_osdet)/N))

re_osdet = var_osdet/(avg_osdet^2)

plot(N,avg_osdet, type='l',main=expression(paste("Estimate of P(",S[n],"> n",beta,")")), 
     ylab = expression(paste("E[",hat(Q)[K],"]")), xlab="n")
#plotCI(avg_osdet, y=NULL, ui = CI$upr, li = CI$low, pch=20, slty=1,scol="red",add=T)

plot(N,logAvg, type='l',main=expression(paste("OSDET Log estimate of P(",S[n],"> n",beta,")")), 
     ylab = expression(paste("log",hat(alpha)[n])), xlab="n")
# plotCI(logAvg, y=NULL, ui = logCI$upr, li = logCI$low, pch=20, slty=1,scol="red",add=T)
# legend("topright", c("95% Confidence Intervals"), col=2, lty=1)


plot(N, re_osdet, type = 'l', main="Squared CV for OSDET", 
     ylab = expression(paste(cv,"(",Q[K],")")), xlab="n")

#######
## variance ratio: p(1-p)/sample variance
# consider for n=5, P(S_n/n > 1)
# Sn/n ~ N(0,1)
p=array(length(N))
for(i in 1:length(N)){
  p[i] = 1-pnorm(1, mean=0, sd=1)
}
var_ratio=p*(1-p)/var_osdet

df = data.frame(Estimate = c(1e-2, 1e-3, 1e-4, 1e-5), 
                OSDET_Var_Ratio = var_ratio, OET_Var_Ratio=)

##########
## simulate for varying values of MC runs, M.####
set.seed(232)

M = seq(from=100, to = 3000, by = 100)
N = 10

cl <- makeCluster(12)
registerDoParallel(cl)

output_gauss2 <- foreach(m = M) %dopar% {
  replicate(m, gauss_osdet(n = N, b = 1, l = 50, mu = 0, sigma2 = 1))
}
beep(2)
stopCluster(cl)

# estimates
mean_out_gauss2 <- lapply(output_gauss2, mean)



####### GAMMA #######

## simulate for varying values of n.
set.seed(232)
M1 = 100
N1 = seq(from=100, to = 1000, by = 100)

cl <- makeCluster(10)
registerDoParallel(cl)

output_gam1 <- foreach(k = 1:length(N1)) %dopar% {
  replicate(M1, gam_osdet(n = N1[k], b = 1, lambda = 10, scale = 0.5, shape = 0.5))
}
beep(2)
stopCluster(cl)

# estimates
mean_out_gam1 <- lapply(output_gam1, mean)

#########
## simulate for varying values of MC runs, M.
set.seed(232)
M2 = seq(from=100, to = 5000, by = 100)
N2 = 10

cl <- makeCluster(10)
registerDoParallel(cl)

output_gam2 <- foreach(m = M2) %dopar% {
  replicate(m, gam_osdet(n = N2, b = 2, lambda = 20, scale = 0.5, shape = 0.5))
}
beep(2)
stopCluster(cl)

# estimates
mean_out_gam2 <- lapply(output_gam2, mean)

