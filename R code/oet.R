###############################
oet <- function(n, b, s = 0, mu, sigma2){

  theta = (b-mu)/sigma2
  mu1 = mu + theta*sigma2
  
  s = sum(rnorm(n))
  
  L = exp(theta*s - n*gauss_phi(theta, mu = mu1, sigma2 = sigma2))
  sbar = s*L
  
Ind <- 0+(sbar > n*b)
Y = L*Ind
return(Y)
}

#############################

set.seed(232)
Mo = 1e6
No = seq(from=1, to = 40, by = 1)
cores=12

  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  out <- foreach(k = No) %dopar% { 
    replicate(Mo, oet(n = k, b = 1, mu = 0, sigma2 = 1))
  }
  beep(5)
  stopCluster(cl)
  
###########################

avg_oet <- array(as.numeric(unlist(out)), dim=c(Mo, length(No)))
  avg_oet <- apply(avg_oet, 2, sum)/Mo
  logAvg_oet <- log(avg_oet)
m2_oet <- array(as.double(unlist(out)), dim=c(Mo, length(No)))
m2_oet <- apply(m2_oet, 2, function(x) sum(x^2))/Mo
  var_oet <- m2_oet - avg_oet^2

re_oet = var_oet/(avg_oet^2)

logCI_oet = data.frame(low = logAvg_oet - 1.645*sqrt(var_oet/No), 
                   upr = logAvg_oet + 1.645*sqrt(var_oet/No))


plot(No, logAvg_oet, type='l',
     main=expression(paste("OET Log estimate of P(",S[n],"> n",beta,")")), 
     ylab = expression(paste("log",hat(alpha)[n])), xlab="n")
# plotCI(logAvg_oet, y=NULL, ui = logCI_oet$upr, 
#        li = logCI_oet$low, pch=20, slty=1,scol="red",add=T)
# legend("topright", c("95% Confidence Intervals"), col=2, lty=1)

plot(No, re_oet, type = 'l', main="Squared CV for Standard OET", 
     ylab = expression(paste(cv,"(",Q[K],")")), xlab="n")

pe=array(length(No))
for(i in 1:length(No)){
  pe[i] = 1-pnorm(1, mean=0, sd=1)
}
var_ratioe=pe*(1-pe)/var_oet
