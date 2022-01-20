### 13.2.2 ###
# install.packages("Rglpk")
library(Rglpk)
LP_solver = function(c, cstr=list(), trace=FALSE){
  Aeq = Reduce(rbind, cstr[names(cstr) %in% "Aeq"])
  aeq = Reduce(c, cstr[names(cstr) %in% "aeq"])
  A = Reduce(rbind, cstr[names(cstr) %in% "A"])
  a = Reduce(c, cstr[names(cstr) %in% "a"])

  sol = Rglpk_solve_LP(
    obj=c, 
    mat=rbind(Aeq, A), 
    dir=c(rep("==", nrow(Aeq)), rep(">=", nrow(A))), 
    rhs=c(aeq, a), 
    verbose=trace
  )

  status = sol$status
  if(status){
    solution = rep(NA, length(c))
  } else {
    solution = sol$solution
  }
  list(solution=solution, status=status)
}

### 13.2.3 ###
QP_solver = function(c, Q, cstr=list(), trace=FALSE){
  Aeq = Reduce(rbind, cstr[names(cstr) %in% "Aeq"])
  aeq = Reduce(c, cstr[names(cstr) %in% "aeq"])
  A = Reduce(rbind, cstr[names(cstr) %in% "A"])
  a = Reduce(c, cstr[names(cstr) %in% "a"])

  sol = try(solve.QP(
    Dmat=Q, 
    dvec=-2*c, 
    Amat=t(rbind(Aeq, A)), 
    bvec=c(aeq, a), 
    meq=nrow(Aeq)
  ), silent=TRUE)
  if(trace){
    cat(sol)
  }
  if(inherits(sol, "try-error")){
    list(solution=rep(NA, length(c)), status=1)
  } else {
    list(solution=sol$solution, status=0)
  }
}

### 13.4 ###
x = apply(EuStockMarkets, 2, function(x) x[-1] / x[-length(x)] - 1)
pftPerf = function(x, w, W0=1000){
  W0 * cumprod(c(1, 1 + x %*% w))
}
nc = ncol(x)
w = rep(1/nc, nc)
plot(pftPerf(x, w), type="l", main="Port")
grid()

### 13.5.1 ###
targetReturn = function(x, target){
  list(Aeq=rbind(colMeans(x)), aeq=target)
}

fullInvest = function(x, target){
  list(Aeq=matrix(1, nrow=1, ncol=ncol(x)), aeq=1)
}

longOnly = function(x){
  list(A=diag(1, ncol(x)), a=rep(0, ncol(x)))
}

# groupBudget = ...

### 13.5.2 ###
MV_QP = function(
  x, target, 
  Sigma=cov(x), ..., 
  cstr=c(fullInvest(x), targetReturn(x, target), longOnly(x), ...), 
  trace=FALSE
  ) {
  
  # quadratic coefficients
  size = ncol(x)
  c = rep(0, size)
  Q = Sigma

  # optimization
  sol = QP_solver(c, Q, cstr, trace)

  # extract weight
  weights = sol$solution
  names(weights) = colnames(x)
  weights
}

w = MV_QP(x, mean(x))
sum(w)
barplot(w, ylim = c(0, 1), las = 2, main = "Test MV Implementation")

### 13.5.3 ###
# install.packages("robustbase")
library(robustbase)
w1 = MV_QP(x, mean(x), Sigma = cov(x))
w2 = MV_QP(x, mean(x), Sigma = covMcd(x)$cov)
par(mfrow=c(1, 2))
barplot(w1, ylim = c(0, 1), las = 2, main = "MV with classical cov")
barplot(w2, ylim = c(0, 1), las = 2, main = "MV with robust cov")

### 13.5.4 ###
w = MV_QP(x, cstr=c(fullInvest(x), longOnly(x)))
barplot(w, ylim = c(0, 1), las = 2, main = "Minimum Variance Portfolio")

### 13.5.5 ###
CVaR_LP = function(
  x, target, alpha=0.95,...,
  cstr=c(fullInvest(x), targetReturn(x, target), longOnly(x)), 
  trace=FALSE
){
  # number of scenarios
  J = nrow(x)

  # number of assets
  size = ncol(x)

  # objective coefficients
  c_weights = rep(0, size)
  c_VaR = 1
  c_Scenarios = rep(1 / ((1 - alpha) * J), J)
  c = c(c_weights, c_VaR, c_Scenarios)

  # extract values from constraint to extend them with CVaR constraints
  Aeq = Reduce(rbind, cstr[names(cstr) %in% "Aeq"])
  aeq = Reduce(c, cstr[names(cstr) %in% "aeq"])
  A = Reduce(rbind, cstr[names(cstr) %in% "A"])
  a = Reduce(c, cstr[names(cstr) %in% "a"])

  # build first two blocks of the constraint matrix
  M1 = cbind(Aeq, simple_triplet_zero_matrix(nrow(Aeq), J+1))
  M2 = cbind(A, simple_triplet_zero_matrix(nrow(A), J+1))

  # identity matrix and vector of zeros
  I = simple_triplet_diag_matrix(1, J)

  # block CVaR constraint (y x +alpha+z_j>=0)
  M3 = cbind(x, rep(1, J), I)

  # block CVaR constraint (z_j >= 0)
  M4 = cbind(simple_triplet_zero_matrix(J, size+1), I)

  # vector of zeros used for the rhs of M3 and M4
  zeros = rep(0, J)

  # combine constraints
  cstr = list(
    Aeq = M1, 
    aeq = aeq, 
    A = rbind(M2, M3, M4), 
    a = c(a, zeros, zeros)
  )

  # optimization
  sol = LP_solver(c, cstr, trace=trace)

  # extract weights
  weights = sol$solution[1:size]
  names(weights) = colnames(x)

  # extract VaR and CVaR
  VaR = sol$solution[size+1]
  CVaR = c(c%*%sol$solution)

  attr(weights, "risk") = c(VaR=VaR, CVaR=CVaR)
  
  weights
}

CVaR_LP(x, mean(x), 0.5)

### 13.6 ###
mu = apply(x, 2, mean)
reward = seq(from=min(mu), to=max(mu), length.out=300)
sigma = apply(x, 2, sd)

Sigma = cov(x)
riskCov = sapply(reward, function(targetReturn){
  w = MV_QP(x, targetReturn, Sigma)
  sd(c(x%*%w))
})
xlim = range(c(sigma, riskCov), na.rm=TRUE)
ylim = range(mu)
plot(
  riskCov, reward, type="l", xlim=xlim, ylim=ylim, 
  xlab="Risk", ylab="Reyard", main="Efficient Frontier"
)
points(sigma, mu, col="steelblue", pch=19, cex=0.8)
text(sigma, mu, labels=colnames(x), pos=2, cex=0.8)
