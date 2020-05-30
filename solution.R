#rough
sr = 1000
for(i in seq(1,144,1)){
  print(sr)
  sr = exp(log(sr) + dt2*(r-q - 0.5*(vol[37,389]^2)) + vol[37,389]*rnorm(1,mean=0,sd=dt2))
}

print(max(sr-1100,0)*exp(-r*3))


x = c()
y = c()
z = c()

for(i in seq(1,181,1)){
  for(j in seq(1,1000,1)){
    x = c(x,i)
    y = c(y,j)
    z = c(z,localvol[i,j])
  }
}

plot(x,y,z)






vol = read.csv("qc_vols2.csv")

#Part 1
s0 = 1000
r = 0.05
q = 0.02
mu = r-q
k = seq(50,2750,(2700/999))
t = c(1/180,seq(1,180,1)/12)
call = function(s0, r, q, k, t, sigma) {
  f = s0*exp((r-q)*t)
  d1 = (log(f/k) + 0.5*t*sigma^2) / (sigma*sqrt(t))
  d2 = d1 - sigma*sqrt(t)
  c = exp(-r*t)*(f*pnorm(d1) - k*pnorm(d2))
  c
}

#print(call(1000,r,q,1100,3,vol[37,389]))
#print(call(1000,r,q,1100,3,localv(3,1100)))

callop = matrix(,nrow=181,ncol=1000)
for(i in seq(1,181,1)){
  for(j in seq(1,1000,1)){
    callop[i,j] = call(s0, r, q, k[j], t[i], vol[i,j])
  }
}
result1 = sum(callop)/(s0*1000*181)

#Part2
dk = 2700/999
dt = 1/12
dsigk = matrix(,nrow=181,ncol=1000)
for(i in seq(1,181,1)){
  for(j in seq(1,1000,1)){
    if(j==1){
      dsigk[i,1] = (vol[i,2] - vol[i,1])/dk
    } else if(j==1000){
      dsigk[i,1000] = (vol[i,1000] - vol[i,999])/dk
    } else {
      dsigk[i,j] = (vol[i,j+1] - vol[i,j-1])/(2*dk)
    }
  }
}

dsigt = matrix(,nrow=181,ncol=1000)
for(i in seq(1,181,1)){
  for(j in seq(1,1000,1)){
    if(i==1){
      dsigt[1,j] = (vol[2,j] - vol[1,j])/dt
    } else if(i==181){
      dsigt[181,j] = (vol[181,j] - vol[180,j])/dt
    } else {
      dsigt[i,j] = (vol[i+1,j] - vol[i-1,j])/(2*dt)
    }
  }
}

d2sig = matrix(,nrow=181,ncol=1000)
for(i in seq(1,181,1)){
  for(j in seq(1,1000,1)){
    if(j==1){
      d2sig[i,1] = 0
    } else if(j==1000){
      d2sig[i,1000] = 0
    } else {
      d2sig[i,j] = (vol[i,j+1] + vol[i,j-1] - 2*vol[i,j])/(dk*dk)
    }
  }
}

lv = function(s0,k,sigma,t,d1sigmak, d1sigmat, d2sigma,r,q) {
  f = s0*exp((r-q)*t)
  x = (log(f/k) + (sigma^2)*t/2)/(sigma*sqrt(t))
  c = sqrt((sigma^2 + 2*sigma*t*(d1sigmat + mu*k*d1sigmak))/(1+2*k*x*sqrt(t)*d1sigmak + k*k*t*sigma*d2sigma + k*k*x*(x-sigma*sqrt(t))*t*(d1sigmak)^2))
  c
}

localvol = matrix(,nrow=181,ncol=1000)
for(i in seq(1,181,1)){
  for(j in seq(1,1000,1)){
    localvol[i,j] = lv(s0,k[j],vol[i,j],t[i],dsigk[i,j],dsigt[i,j],d2sig[i,j],r,q)
  }
}

for(i in seq(1,181,1)){
  for(j in seq(1,1000,1)){
    if(is.nan(localvol[i,j])==TRUE){
      localvol[i,j] = (localvol[i-1,j] + localvol[i-2,j] + localvol[i-3,j] + localvol[i-4,j])/4 
    }
  }
}

sim = 10000
dt1 = 1/52

localv = function(t1,s1) {
  vol1 = 0
  
  tindex = 0
  if(t1<=t[1]){
    tindex = 1
  } else if(t1>=t[181]){
    tindex = 181
  } else {
    tindex = floor((t1 - t[2])/dt) + 2
  }
  
  if((s1-50)%%dk == 0){
    vol1 = localvol[tindex, ((s1-50)/dk) + 1]
  } else {
    if(s1<=k[1]){
      vol1 = localvol[tindex, 1]
    } else if(s1>=k[1000]) {
      vol1 = localvol[tindex, 1000]
    } else {
      k1 = ceiling((s1-50)/dk)
      l1 = localvol[tindex,k1]
      l2 = localvol[tindex,k1+1]
      vol1 = l1 + (l2-l1)*(s1-k[k1])/(dk) 
    }
  }
  vol1
}

payoff = matrix(,nrow=181,ncol=1000)
for(i in seq(1,181,1)){
  for(j in seq(1,1000,1)){
    priceS = 0
    for(iter in seq(1,sim,1)) {
      n = ceiling(t[i]/dt1)
      s = s0
      for(iter1 in seq(1,n,1)) {
        s = exp(log(s) + dt1*(mu - 0.5*localv(iter1*dt1,s)^2) + localv(iter1*dt1,s)*rnorm(1,mean=0,sd=dt1))
      }
      priceS = priceS + max(s-k[j],0)
    }
    priceS = priceS/sim
    payoff[i,j] = priceS*exp(-r*t[i])
  }
}

result2 = sum(payoff)/(s0*1000*181)


#rough
vol4 = c()
k4 = c()
priceS = 0
sim=10000
for(iter in seq(1,sim,1)) {
  n = ceiling(t[61]/dt1)
  s = s0
  for(iter1 in seq(1,n,1)) {
    if(iter==1){
      vol4 = c(vol4, localv(iter1*dt1,s))
      k4 = c(k4,s)
    }
    s = exp(log(s) + dt1*(mu - 0.5*localv(iter1*dt1,s)^2) + sqrt(dt1)*localv(iter1*dt1,s)*rnorm(1,mean=0,sd=1))
  }
  priceS = priceS + max(s-k[350],0)
}
priceS = priceS/sim
print(priceS*exp(-r*t[61]))
plot(k4,vol4)
v4a = rep(vol[61,350],length(k4))
lines(k4,v4a)


#Part3
kb = 1.1*s0
ll = 0.4*s0
ul = 1.5*s0
nsim = 20000
dt2 = 1/48

tplot = seq(1/48,3,1/48)
lvolplot = c()
volplot = rep(vol[37,389],144)
upb = 0
downb = 0

barrieropt = matrix(,nrow=nsim, ncol=144)
for(i in seq(1,nsim,1)){
  s=s0
  for(j in seq(1,144,1)){
    reqvol = localv(j*dt2,s)
    if(i==1){
      lvolplot = c(lvolplot,reqvol)
    }
    s = exp(log(s) + dt2*(mu - 0.5*reqvol^2) + sqrt(dt2)*reqvol*rnorm(1,mean=0,sd=1))
    barrieropt[i,j] = s
    if(j%%4==0){
      if(s>=ul || s<=ll){
        if(s>=ul){
          upb = upb + 1
        }
        if(s<=ll){
          downb = downb + 1
        }
        barrieropt[i,144] = 0
        break
      }
    }
  }
}

priceB = 0
for(i in seq(1,20000,1)){
  priceB = priceB + max(barrieropt[i,144] - kb,0)
}

price = exp(-r*3)*priceB/nsim

plot(tplot,lvolplot,type="l",col="red")
lines(tplot,volplot,col="green")

#Part4
kindex = 350
tindex = 61
lvol = vol[tindex,kindex]
hvol = vol[tindex,kindex+1]
sigvol = lvol + (s0 - k[kindex])*(hvol - lvol)/(dk)

ML = sigvol*sqrt(t[tindex])*(-5)
MU = sigvol*sqrt(t[tindex])*(5)

Smax = exp(log(s0) + MU)
Smin = exp(log(s0) + ML)

m=400
n=500

tfd = seq(0,t[tindex],t[tindex]/n)

xfd = c()
for(i in seq(0,m+1,1)){
  xfd = c(xfd, (ML + i*(MU-ML)/(m+1)))
}

v = matrix(,nrow=n+1, ncol=m+2)
l = matrix(,nrow=n+1, ncol=m+2)
c = matrix(,nrow=n+1, ncol=m+2)
u = matrix(,nrow=n+1, ncol=m+2)
L = matrix(,nrow=n+1, ncol=m+2)
C = matrix(,nrow=n+1, ncol=m+2)
U = matrix(,nrow=n+1, ncol=m+2)
theta = 0.5

dx = (MU - ML)/(m+1)
#dx=1
dt3 = t[tindex]/n

for(i in seq(1,n+1,1)){
  for(j in seq(1,m+2,1)){
    a = (mu - 0.5*localv(tfd[i],s0*exp(xfd[j]))^2)
    #a = mu
    c1 = 0.5*localv(tfd[i],s0*exp(xfd[j]))^2
    f = -r
    l[i,j] = -theta*dt3*(-0.5*(a/dx) + c1/(dx*dx))
    c[i,j] = 1-theta*dt3*(-2*c1/(dx*dx) + f)
    u[i,j] = -theta*dt3*(0.5*(a/dx) + c1/(dx*dx))
    L[i,j] = (1-theta)*dt3*(-0.5*(a/dx) + c1/(dx*dx))
    C[i,j] = 1 + (1-theta)*dt3*(-2*c1/(dx*dx) + f)
    U[i,j] = (1-theta)*dt3*(0.5*(a/dx) + c1/(dx*dx))
  }
}

for(j in seq(1,m+2,1)){
  v[n+1,j] = max(s0*exp(xfd[j]) - k[kindex], 0)
}

for(i in seq(1,n+1,1)){
  v[i,1] = 0
  v[i,m+2] = s0*exp(-q*(t[tindex]-tfd[i])) - k[kindex]*exp(-r*(t[tindex]-tfd[i]))
}

for(i in seq(n,1,-1)){
  A = matrix(0, nrow=m,ncol=m)
  b = c()
  for(j in seq(2,m+1,1)){
    if(j==2){
      A[j-1,j-1] = c[i,j]
      A[j-1,j] = u[i,j]
      b = c(b, L[i+1,j]*v[i+1,j-1] + C[i+1,j]*v[i+1,j] + U[i+1,j]*v[i+1,j+1] + (1-theta)*v[i+1,1]*L[i+1,2] + theta*v[i,1]*l[i+1,2])
    } else if(j==m+1){
      A[j-1,j-1] = c[i,j]
      A[j-1,j-2] = l[i,j]
      b = c(b, L[i+1,j]*v[i+1,j-1] + C[i+1,j]*v[i+1,j] + U[i+1,j]*v[i+1,j+1] + (1-theta)*v[i+1,m+2]*U[i+1,m+1] + theta*v[i,m+2]*u[i+1,m+1] - u[i,j]*v[i,j+1])
    } else {
      A[j-1,j-1] = c[i,j]
      A[j-1,j-2] = l[i,j]
      A[j-1,j] = u[i,j]
      b = c(b, L[i+1,j]*v[i+1,j-1] + C[i+1,j]*v[i+1,j] + U[i+1,j]*v[i+1,j+1])
    }
  }
  soln = solve(A,b)
  for(j in seq(2,m+1,1)){
    v[i,j] = soln[j-1] 
  }
}


#rough
kindex = 350
tindex = 61
lvol = vol[tindex,kindex]
hvol = vol[tindex,kindex+1]
sigvol = lvol + (s0 - k[kindex])*(hvol - lvol)/(dk)

ML = sigvol*sqrt(t[tindex])*(-5)
MU = sigvol*sqrt(t[tindex])*(5)

Smax = exp(log(s0) + MU)
Smin = exp(log(s0) + ML)

m=400
n=500

tfd = seq(0,t[tindex],t[tindex]/n)

xfd = c()
for(i in seq(0,m+1,1)){
  xfd = c(xfd, (Smin + i*(Smax-Smin)/(m+1)))
}


Ac = matrix(,nrow=n+1,ncol=m+2)
Bc = matrix(,nrow=n+1,ncol=m+2)
Cc = matrix(,nrow=n+1,ncol=m+2)
v = matrix(,nrow=n+1,ncol=m+2)

for(i in seq(1,n+1,1)){
  for(j in seq(1,m+2,1)){
    sigmac = localv(tfd[i],xfd[j])
    Ac[i,j] = 0.5*mu*j*dt3 - 0.5*sigmac*sigmac*j*j*dt3
    Bc[i,j] = 1 + sigmac*sigmac*j*j*dt3 + r*dt3
    Cc[i,j] = -0.5*mu*j*dt3 - 0.5*sigmac*sigmac**j*j*dt3
  }
}


for(j in seq(1,m+2,1)){
  v[n+1,j] = max(xfd[j] - k[kindex], 0)
}

for(i in seq(1,n+1,1)){
  v[i,1] = 0
  v[i,m+2] = xfd[m+2]*exp(-q*(t[tindex]-tfd[i])) - k[kindex]*exp(-r*(t[tindex]-tfd[i]))
}

for(i in seq(n,1,-1)){
  A = matrix(0, nrow=m,ncol=m)
  b = c()
  for(j in seq(2,m+1,1)){
    if(j==2){
      A[j-1,j-1] = Bc[i,j]
      A[j-1,j] = Cc[i,j]
      b = c(b,v[i+1,j])
    } else if(j==m+1){
      A[j-1,j-1] = Bc[i,j]
      A[j-1,j-2] = Ac[i,j]
      b = c(b, v[i+1,j]-Cc[i,j]*v[i,j+1])
    } else {
      A[j-1,j-1] = Bc[i,j]
      A[j-1,j-2] = Ac[i,j]
      A[j-1,j] = Cc[i,j]
      b = c(b, v[i+1,j])
    }
  }
  soln = solve(A,b)
  for(j in seq(2,m+1,1)){
    v[i,j] = soln[j-1] 
  }
}
