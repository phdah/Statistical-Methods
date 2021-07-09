######### Home assignment 1 ######### 
# Part 2
# Philip Dahlqvist-Sjöberg
# 29-05-2020
# NOTE: Run all these libraries first.
library(dplyr) # Main data handler
library(ICAOD) # Locally function
library(matrixcalc) # svd.inverse
library(nloptr)

# Function to write matrices to Latex code
write_matex <- function(x) {
  begin <- "\\begin{bmatrix}"
  end <- "\\end{bmatrix}"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  writeLines(c(begin, X, end))
}

###### 1 ###### 
# a)
# Standardized information matrix function
sINF_matrix <- function(design){
  x <- design[1,]
  w <- design[2,]
  M <- matrix(0,2,2)
  for(i in 1:length(x)){
  z1 <- c(1,x[i]^2)
  z2 <- c(x[i]^2, x[i]^4)
  M <- M + (w[i] * matrix(c(z1, z2), 2, byrow = T))
  }
  return(M)
}

# Design
design <- matrix(c(0,3,6,1/3,1/3,1/3), ncol = 3, byrow = T)

# Run stand. information matrix
M <- sINF_matrix(design)
write_matex(M)
# b)
# Region 
grid <- seq(0,6, by = 0.1)

# Standardized predictor variance of the design
sVar_matrix <- function(x, M){
  f_x <- matrix(c(1, x^2), 2)
  S <- t(f_x) %*% solve(M) %*% f_x
  return(S)
}

plot(grid, sapply(grid, sVar_matrix, M), type = "l", xlab = "x", 
     ylab = expression(paste('d(x,', xi,beta, ')')))

# c)
# The figure in b), shows that the deisgn is not optimal. For a optimal
# design, all three x-values would be maximum in the plot.

locally(formula =~ b0 + b1*x^2, predvars = "x", parvars = c("b0","b1"),
        lx = 0 ,ux = 6, inipars = c(1, 0),iter = 100, k = 2) # points >= parameters
# Two points at c(0,6) are optimal with equal weight. 

# Plot same as b), but with optimal design.
design_star <- matrix(c(0,6, 0.5, 0.5), 2, byrow = T)
plot(grid, sapply(grid, sVar_matrix, M = sINF_matrix(design_star)), type = "l", xlab ="x", 
     ylab = expression(paste('d(x,', xi,beta, ')')))

# d)
RE_D <- function(design, design_star, p){
  (det(solve(sINF_matrix(design_star))) / det(solve(sINF_matrix(design))))^(1/p)
}

# How many %subject non-D-optimal design need to be same as efficient as D optimal design
RE_D_value <- RE_D(design, design_star, p = 2)
round(((1/RE_D_value)-1)*100,1)

# 17.7% more subjects needed to get same efficiency as optimal design

# e)
#criterion function for D-optimality
crit = function(w, x){
  Mi = matrix(0,2,2)
  for(i in 1:length(x)){
    z1 <- c(1,x[i]^2)
    z2 <- c(x[i]^2, x[i]^4)
    z <- matrix(c(z1, z2), 2, byrow = T)
    Mi = Mi+w[i]*z
  }
  sum(1/(eigen(Mi)$values))
}


# We use this function to search optimum value of w while x values are fixed,
# which minimize the criterion function.

critw = function(w, x){
  g = w^2/sum(w^2)  
  crit(g, x) 
}


# We use this function to search optimum values of x and w simulatneously.
# which minimize the criterion function.

critwx = function(wx){
  r = length(wx)/2
  c = length(wx)/r
  mat = matrix(wx,r,c,byrow =FALSE)
  x = mat[,1]
  w = mat[,2]
  w = w^2/sum(w^2)
  crit(w,x) 
} 

# Simplify design by removing values with low weight. and average the value wich
# has less difference(0.001) and sum-up their weight.

Consolidate = function(x,w,ex,ew){
  de = data.frame(cbind(x,w))
  de = de[which(de$w>ew),]
  x = de[,1]
  
  k = dim(de)[1]
  g = rep(0,k)
  for(i in 1:k){
    
    if(g[i]==0){
      g[i] = max(g)+1
      for(j in (i+1):k) {
        if(i<k)
          if(abs(x[i]-x[j])<=ex)
          {g[j]=g[i]} 
      }}}
  de$g = g
  f1 = aggregate(x~g,de,mean)[,2]
  f2 = aggregate(w~g,de,sum)[,2]
  newx = data.frame(x=f1,w=f2) 
  newx
}

# Using initial design now we minimize varience for both x and w 

final = function(xop, wop, iter, accuracy){
  i=1;improvement=1
  while((i<iter) &  (improvement>accuracy)){
    i=i+1  
    last = cbind(xop,wop)
    n = length(xop)
    
    wxinit = c(xop,sqrt(wop))
    
    opdes = neldermead(x0=wxinit,fn=critwx, lower=c(rep(xl,n),rep(0,n)),upper= c(rep(xu,n),rep(1,n)))$par
    
    r = length(opdes )/2
    c = length(opdes )/r
    mat = matrix(opdes,r,c,byrow =FALSE)
    x = mat[,1]
    w = mat[,2]
    w = w^2/sum(w^2)  
    
    #again simplify the sedign
    sdes = Consolidate(x, w, 0.01, 0.001)
    xop = sdes$x
    wop = sdes$w
    
    if(length(xop)!=nrow(last)){
      improvement=999
    }  
    else
      improvement=max(abs(last-cbind(xop,wop)))
  }
  rr = data.frame(x=xop,w=wop,var=crit(wop,xop),improvement=improvement)
  rr[order(rr$x),]
}

# Finally utilization of all above functions

# Design space for x 

xl=0
xu=6


# Starting values of design point with weight

x=runif(41,xl,xu);
winit = rep(1,length(x)) / length(x);
winit = sqrt(winit);

# Fixed value of x and find optimal value for w by maximizing function critw

wop = neldermead(x0=winit,fn=critw,x=x)$par

# Covert wop between 0 and 1 such that its sum is 1
wop = (wop^2)/sum(wop^2);
initialdes = data.frame(x,w=wop,var=crit(wop,x))

# Simplify the design

sdes = Consolidate(x,wop,0,0.0001)
xop = sdes$x
wop = sdes$w
siminitialdes = data.frame(x=xop,w=wop,var=crit(wop,xop))

# Final

res_final <- final(xop,wop,iter=100,accuracy=0.000000001)
res_final

# A-optimal design
Aopt <-function(x, w, b0, b1, fimfunc){
  sum(diag(svd.inverse(fimfunc(x = x, w = w, b0 = b0, b1 = b1))))
}

res1<-locally(formula = ~ b0 + b1*x^2, predvars = "x", 
              parvars = c("b0", "b1"),
              lx = 0, ux = 6, crtfunc = Aopt, inipars =c(1,2),k = 2,iter = 100)
res1
# Two points, c(~0,6), in both both methods. 

###### 2 ######
# a)
data <- data.frame(matrix(c(1.6907 , 59,  6,
                  1.7242,  60 , 13,
                  1.7552 , 62 , 18,
                  1.7842,  56 , 28,
                  1.8113,  63 , 52,
                  1.8369,  59 , 53,
                  1.8610,  62 , 61,
                  1.8839,  60 , 60), ncol = 3, byrow = T))
names(data) <- c("Dose", "Beetles", "Killed")

# Extend the observed values to each have their own row in dataframe
data_2a = matrix(c(rep(1.6907,59), rep(1.7242, 60), rep(1.7552, 62), rep(1.7842, 56), rep(1.8113, 63),
                   rep(1.8369, 59), rep(1.8610, 62), rep(1.8839, 60), rep(0, 53), rep(1, 6), rep(0, 47), 
                   rep(1, 13), rep(0, 44), rep(1, 18), rep(0, 28), rep(1, 28), rep(0, 11), rep(1, 52),
                   rep(0, 6), rep(1, 53), rep(0, 1), rep(1, 61), rep(1, 60)), byrow = F, ncol = 2)

data_2a = data.frame(data_2a)

colnames(data_2a) = c("Dose", "Killed")

# Fit logistic regression
fit_logist <- glm(Killed~Dose + I(Dose^2), data = data_2a, family = binomial(link = "logit"))

# Save parameters
b0_log <- fit_logist$coefficients[[1]]
b1_log <- fit_logist$coefficients[[2]]
b2_log <- fit_logist$coefficients[[3]]

# b)
#criterion function for D-optimality, logistic model
crit = function(w, x, b0, b1, b2){
  p1 = 1+exp(-(b0+b1*x+b2*x^2))
  e = 1/p1
  v = e*(1-e)
  Mi = matrix(0,3,3)
  for(i in 1:length(x)){
    z = c(1,x[i], x[i]^2)
    z1 = z%*%t(z)
    Mi = Mi+w[i]*v[i]*z1
  }
  -log(det(Mi))
}


# We use this function to search optimum value of w while x values are fixed,
# which minimize the criterion function.

critw = function(w, x, b0, b1, b2){
  g = w^2/sum(w^2)  
  crit(g, x, b0, b1, b2) 
}


# We use this function to search optimum values of x and w simulatneously.
# which minimize the criterion function.

critwx = function(wx, b0, b1, b2){
  r = length(wx)/2
  c = length(wx)/r
  mat = matrix(wx,r,c,byrow =FALSE)
  x = mat[,1]
  w = mat[,2]
  w = w^2/sum(w^2)
  crit(w,x,b0,b1,b2) 
} 

# Simplify design by removing values with low weight. and average the value wich
# has less difference(0.001) and sum-up their weight.

Consolidate = function(x,w,ex,ew){
  de = data.frame(cbind(x,w))
  de = de[which(de$w>ew),]
  x = de[,1]
  
  k = dim(de)[1]
  g = rep(0,k)
  for(i in 1:k){
    
    if(g[i]==0){
      g[i] = max(g)+1
      for(j in (i+1):k) {
        if(i<k)
          if(abs(x[i]-x[j])<=ex)
          {g[j]=g[i]} 
      }}}
  de$g = g
  f1 = aggregate(x~g,de,mean)[,2]
  f2 = aggregate(w~g,de,sum)[,2]
  newx = data.frame(x=f1,w=f2) 
  newx
}

# Using initial design now we minimize varience for both x and w 

final = function(xop, wop, b0, b1, b2, iter, accuracy){
  i=1;improvement=1
  while((i<iter) &  (improvement>accuracy)){
    i=i+1  
    last = cbind(xop,wop)
    n = length(xop)
    
    wxinit = c(xop,sqrt(wop))
    
    opdes = neldermead(x0=wxinit,fn=critwx, b0 = b0, b1 = b1, b2 = b2,lower=c(rep(xl,n),rep(0,n)),upper= c(rep(xu,n),rep(1,n)))$par
    
    r = length(opdes )/2
    c = length(opdes )/r
    mat = matrix(opdes,r,c,byrow =FALSE)
    x = mat[,1]
    w = mat[,2]
    w = w^2/sum(w^2)  
    
    #again simplify the sedign
    sdes = Consolidate(x, w, 0.01, 0.001)
    xop = sdes$x
    wop = sdes$w
    
    if(length(xop)!=nrow(last)){
      improvement=999
    }  
    else
      improvement=max(abs(last-cbind(xop,wop)))
  }
  rr = data.frame(x=xop,w=wop,var=crit(wop,xop,b0,b1,b2),improvement=improvement)
  rr[order(rr$x),]
}

# Finally utilization of all above functions

# Design space for x 

xl=1.5
xu=2

b0 = b0_log; b1 = b1_log; b2 = b2_log

# Starting values of design point with weight

x=runif(41,xl,xu);
winit = rep(1,length(x)) / length(x);
winit = sqrt(winit);

# Fixed value of x and find optimal value for w by maximizing function critw

wop = neldermead(x0=winit,fn=critw,x=x,b0 = b0,b1=b1, b2 = b2)$par

# Covert wop between 0 and 1 such that its sum is 1
wop = (wop^2)/sum(wop^2);
initialdes = data.frame(x,w=wop,var=crit(wop,x,b0,b1, b2))

# Simplify the design

sdes = Consolidate(x,wop,0,0.0001)
xop = sdes$x
wop = sdes$w
siminitialdes = data.frame(x=xop,w=wop,var=crit(wop,xop,b0, b1, b2))

# Final

res_final <- final(xop,wop,b0_log,b1_log, b2_log,iter=100,accuracy=0.000000001)
res_final

# D-optimal design locally
res1 <- locally(formula = ~ exp(b0+b1*x+b2*x^2)/(1+exp(b0+b1*x+b2*x^2)), predvars = "x", 
              parvars = c("b0", "b1","b2"), 
              lx =1.5, ux =2, inipars =c(b0_log,b1_log,b2_log), k=4, iter=100)
res1
# Same result for both methods

# c)
# RE_d test for different parameter input
# Standardized information matrix function
crit_RE = function(w, x, b0, b1, b2){
  p1 = 1+exp(-(b0+b1*x+b2*x^2))
  e = 1/p1
  v = e*(1-e)
  Mi = matrix(0,3,3)
  for(i in 1:length(x)){
    z = c(1,x[i], x[i]^2)
    z1 = z%*%t(z)
    Mi = Mi+w[i]*v[i]*z1
  }
  Mi
}

# Optimal design from b)
w <- c(0.3114205, 0.1885823, 0.1885764, 0.3114207)
x <- c(1.517310, 1.590779, 1.737716, 1.811185)

# Changes of the parameter values
e1 <- -1
e2 <- -1
e3 <- -1

# Calculate RE_d
RE_D_value <- (det(solve(crit_RE(w,x,b0,b1,b2)))/det(solve(crit_RE(w,x,b0+e1,b1+e2,b2+e3))))^(1/3)
RE_D_value
N_obs <- ((1/RE_D_value)-1)*100
round(c(RE_D_value, N_obs),4)
# Determinant of cov-matrix is a scalar that tells us how much variance the design has. 
# For values >1, the glm parameters give more variance than altered parameters. 

# d)
# Bayesian- Sequentiall design, see pdf

# e)
# E-optimal criterion
crit = function(w, x, b0, b1, b2){
  p1 = 1+exp(-(b0+b1*x+b2*x^2))
  e = 1/p1
  v = e*(1-e)
  Mi = matrix(0,3,3)
  for(i in 1:length(x)){
    z = c(1,x[i], x[i]^2)
    z1 = z%*%t(z)
    Mi = Mi+w[i]*v[i]*z1
  }
  u <- eigen(Mi)$values
  -min(u)
}

# Design space for x 

xl=1.5
xu=2

b0 = b0_log; b1 = b1_log; b2 = b2_log

# Starting values of design point with weight

x=runif(41,xl,xu);
winit = rep(1,length(x)) / length(x);
winit = sqrt(winit);

# Fixed value of x and find optimal value for w by maximizing function critw

wop = neldermead(x0=winit,fn=critw,x=x,b0 = b0,b1=b1, b2 = b2)$par

# Covert wop between 0 and 1 such that its sum is 1
wop = (wop^2)/sum(wop^2);
initialdes = data.frame(x,w=wop,var=crit(wop,x,b0,b1, b2))

# Simplify the design

sdes = Consolidate(x,wop,0,0.0001)
xop = sdes$x
wop = sdes$w
siminitialdes = data.frame(x=xop,w=wop,var=crit(wop,xop,b0, b1, b2))

res_final <- final(xop,wop,b0,b1, b2,iter=100,accuracy=0.000000001)
res_final 

# E-optimal design, locally
Eopt <-function(x,w,b0,b1,b2, fimfunc){
  e1 <- eigen(svd.inverse(fimfunc(x = x, w = w, b0=b0_log,b1=b1_log,b2=b2_log)))$values
  max(e1)
}

res3<-locally(formula = ~ exp(b0+b1*x+b2*x^2)/(1+exp(b0+b1*x+b2*x^2)), predvars = "x",
              parvars = c("b0", "b1","b2"), k = 3,
              lx = 1.5, ux = 2, inipars = c(b0_log, b1_log, b2_log), iter = 100)
res3

###### 3 ######
# a)
# D-optimal design E_max model
res1 <- locally(formula = ~E_0 + (E_M*x/(ED_50 + x)), predvars = "x", 
              parvars = c("E_0", "E_M","ED_50"),
              lx = 0, ux = 100, inipars = c(22, 11.2, 70), k = 3, iter = 100)
res1


# b)
# 1
# Standardized information matrix up to (Mi), also see pdf
#criterion function for D-optimality 3PL (parameter c=t due to package overright)
crit<-function(w,x,a,b, t){
  p1<-t+(1-t)/(1+exp(-a*(x-b)))
  v<-p1*(1-p1)
  Mi<-matrix(0,3,3)
  for(i in 1:length(x)){
    za <- -((-x[i]+b)/(t*exp(-a*(x[i]-b))+1))
    zb <- -((a)/(t*exp(-a*(x[i]-b))+1))
    zt <- 1/(exp(a*(x[i]-b))+t)+1/(1-t)
    z <- c(za, zb, zt)
    z1<-z%*%t(z)
    Mi<-Mi+w[i]*v[i]*z1
  }
  
# for other optimality criterion function change here 
-log(det(Mi))
}

# 2

# We use this function to search optimum value of w while x values are fixed,
# which minimize the criterion function.

critw = function(w, x, a, b, t){
  g = w^2/sum(w^2)  
  crit(g, x, a, b, t) 
}


# We use this function to search optimum values of x and w simulatneously.
# which minimize the criterion function.

critwx = function(wx, a, b, t){
  r = length(wx)/2
  c1 = length(wx)/r
  mat = matrix(wx,r,c1,byrow =FALSE)
  x = mat[,1]
  w = mat[,2]
  w = w^2/sum(w^2)
  crit(w,x, a, b, t) 
} 

# Simplify design by removing values with low weight. and average the value wich
# has less difference(0.001) and sum-up their weight.

Consolidate = function(x,w,ex,ew){
  de = data.frame(cbind(x,w))
  de = de[which(de$w>ew),]
  x = de[,1]
  
  k = dim(de)[1]
  g = rep(0,k)
  for(i in 1:k){
    
    if(g[i]==0){
      g[i] = max(g)+1
      for(j in (i+1):k) {
        if(i<k)
          if(abs(x[i]-x[j])<=ex)
          {g[j]=g[i]} 
      }}}
  de$g = g
  f1 = aggregate(x~g,de,mean)[,2]
  f2 = aggregate(w~g,de,sum)[,2]
  newx = data.frame(x=f1,w=f2) 
  newx
}

# Using initial design now we minimize varience for both x and w 

final = function(xop, wop, a, b, t, iter, accuracy){
  i=1;improvement=1
  while((i<iter) &  (improvement>accuracy)){
    i=i+1  
    last = cbind(xop,wop)
    n = length(xop)
    
    wxinit = c(xop,sqrt(wop))
    
    opdes = neldermead(x0=wxinit,fn=critwx,  a=a, b=b, t=t,lower=c(rep(xl,n),rep(0,n)),upper= c(rep(xu,n),rep(1,n)))$par
    
    r = length(opdes )/2
    c1 = length(opdes )/r
    mat = matrix(opdes,r,c1,byrow =FALSE)
    x = mat[,1]
    w = mat[,2]
    w = w^2/sum(w^2)  
    
    #again simplify the sedign
    sdes = Consolidate(x, w, 0.01, 0.001)
    xop = sdes$x
    wop = sdes$w
    
    if(length(xop)!=nrow(last)){
      improvement=999
    }  
    else
      improvement=max(abs(last-cbind(xop,wop)))
  }
  rr = data.frame(x=xop,w=wop,var=crit(wop,xop, a, b, t),improvement=improvement)
  rr[order(rr$x),]
}

# Finally utilization of all above functions

# Design space for x 

xl=-3
xu=3

a=0.5; b=1; t=0.05

# Starting values of design point with weight

x=runif(41,xl,xu);
winit = rep(1,length(x)) / length(x);
winit = sqrt(winit);

# Fixed value of x and find optimal value for w by maximizing function critw

wop = neldermead(x0=winit,fn=critw, x = x, a = a, b = b, t = t)$par

# Covert wop between 0 and 1 such that its sum is 1
wop = (wop^2)/sum(wop^2);
initialdes = data.frame(x,w=wop,var=crit(wop,x, a, b, c))

# Simplify the design

sdes = Consolidate(x,wop,0,0.0001)
xop = sdes$x
wop = sdes$w
siminitialdes = data.frame(x=xop,w=wop,var=crit(wop,xop, a, b, c))

# Final

res_final <- final(xop,wop, a, b, t,iter=100,accuracy=0.000000001)
res_final
write_matex(round(res_final[,1:2],4))

# D-optimal design, locally 3PL
res2 <- locally(formula = ~ c + (1-c)/(1+exp(-a*(x-b))), predvars = "x", 
                parvars = c("a", "b","c"), k = 3,
                lx = -3, ux = 3, inipars = c(0.5, 1, 0.05), iter = 100)
res2

############