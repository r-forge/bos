# EM algorithm for a mixture of d-variate BOS model (with cond. indep.).
clustMultiBOS <- function (x,m,k,ntrials=1) {
# input:
#   x   ordinal data set [n,d]
#   m   nb of modalities [d]
#   k   number of components
#   ntrials nb of random trials for EM OPTIONAL
# output:
#   par_ml parameters (structure)
#           $pmix   proportions [k]
#           $mu  true modality {1,...,m} [k,d]
#           $pknow   P(sait) [0,1] [k,d]
#   ml      ml value
#   bic     bic value
#   tik   conditional probabilities [n,k]

n = nrow(x)
d = ncol(x)
ml = -Inf
par=list()

# mu to explore at each imbricated EM (intra-component EM)
tabmu = matrix(0,d,max(m))
for (id in 1:d) {tabmu[id,1:m[id]] = 1:m[id]}

for (itrial in 1:ntrials){
  mlold = -Inf
  iter = 1
  # mu, pknow and pmix at random
  par$mu = matrix(0,k,d)
  par$pknow = matrix(runif(k*d),k,d)
  for (ik in 1:k){
    for (id in 1:d){
       par$mu[ik,id] = floor(runif(1,1,m[id])+1)    
  }}
  #library(gtools)
  par$pmix = rdirichlet(1,rep(1,k)) 

  nostop = 1

  while (nostop){
    # E step
    # ------
    # compute all fikj
    cat('*')
    fikj = array(0,dim=c(n,k,d))
    for (id in 1:d){
      for (ik in 1:k){
        pallx = matrix(0,m[id],1)
        for (i in 1:m[id]){  
          fikj[x[,id]==i,ik,id] = pej(i,m[id],m[id],par$mu[ik,id],par$pknow[ik,id])
    }}}
  
    # compute all fik
    fik = matrix(par$pmix,n,k,byrow=TRUE) * apply(fikj,1:2,prod)
    # compute all fi
    fi = matrix(rowSums(fik),n,1)
    # compute all tik
    tik = condprobnolog(fik,fi)
    # compute ml
    mlnew = sum(log(fi))
  
    # M step
    # ------
    # each component is estimated independently
    for (ik in 1:k){
      tmp=dordiem(x,m,tabmu,par$pknow[ik,],tik[,ik])
      par$mu[ik,]=tmp$mu_ml
      par$pknow[ik,]=tmp$p_ml
    }
    nk = matrix(colSums(tik),k,1)
    par$pmix = t(nk / sum(nk))

    # stop?
    if (abs(mlnew-mlold) < 1e-3){
      nostop = 0
      if (mlnew > ml) {
        ml = mlnew
        par_ml = par
    }}

    mlold = mlnew
    iter = iter+1

  } # nostop
        
} # ntrials

bic = ml - 0.5*(k-1+k*d)*log(n)

res=list(par_ml=par_ml,ml=ml,bic=bic,tik=tik)

return(res)
}

###########################################

condprobnolog <- function(f,sumf){
# Conditional probabilities of belonging to clusters (version without log).
# input:
#   f       components density [n,k]
#	  sumf 	  mixture density [n]
# output:
# 	t	conditional probabilities [n,k]
  n = nrow(f)
  k = ncol(f)  
  t = f / (rep(sumf,k,bycol=TRUE)) # t [n,k]
  return(t)
}

# ====== Estimation of the BOS model via EM algoritm =====
ordiem <- function (x,m,tabmu0,tabp0,w=NULL) {
# input:
#   x   ordinal data set [n,1]
#   m   nb of modalities
#   tabmu0  all starting mu0 to try
#   tabp0 all starting p0 to try for each starting mu0
#   w   weights [n,1] OPTIONAL
# output:
#   mu_ml  ml true modality {1,...,m}
#   p_ml   ml P(sait) [0,1]
#   ml     ml value
  
  n = length(x)
# all weights t to 1 of missing 
  if (is.null(w)) {w = rep(1,n)}
  p_ml = NaN
  mu_ml = NaN
  
  ntot = sum(w) # sample size
  ml = -Inf
  
  # try all mu values
  for (mu in tabmu0){
    mlold = -Inf
    for (p in tabp0){
      nostop = 1
      while (nostop){
      # E step
  
  # first: compute px for each modality
  pallx = rep(0,m)
   
  for (i in 1:m) {pallx[i] = pej(i,m,m,mu,p)}  
  
  # second: affect each modality value to the corresponding units
  px = rep(0,n)
  for (i in 1:m)  {px[x==i] = max(1e-300,pallx[i])}
  #print(px)
  mlnew = sum(w*log(px)) # log-likelihood
  #print(log(max(1e-300,px)))
  
  #print(w*log(max(1e-300,px)))
  
  # first: compute pxz1 for each modality
  pallxz1 = matrix(0,m,m-1)
  for (i in 1:m){
    for (j in 1:(m-1)){
      z1tozmm1 = matrix(0,1,m-1)
      z1tozmm1[j] = 1
      pallxz1[i,j] = max(1e-300,pej(i,m,m,mu,p,z1tozmm1)) / max(1e-300,pallx[i])
    }
  }

  # second: affect each modality value to the corresponding units           
  pxz1 = matrix(0,n,m-1)

  for (i in 1:m){
  whereisi = (x==i)
  pxz1[whereisi,] = matrix(1,sum(whereisi),1)%*%pallxz1[i,]              
  }
   
  pxz1 = (as.matrix(w)%*%matrix(1,1,m-1)) * pxz1;
  
  # M step
  # ------
  
  p = sum(pxz1) / (ntot*(m-1))
  
  
# print("-----")  
# cat('ancien maxL :',mlold,' nouveau maxL :',mlnew,'\n')
# cat('mu :',mu,'p :',p,'\n')

if (abs(mlnew-mlold)/ntot < 1e-3) {
  nostop = 0
  if (mlnew > ml)  {
    ml = mlnew
    p_ml = p
    mu_ml = mu
  }
  } 
  mlold = mlnew
  
}#while
  
}#p
}#mu
  return(list(mu_ml=mu_ml,p_ml=p_ml,ml=ml))
}


############################################################################# 

# ====== EM algorithm for a d-variate BOS model. =====
dordiem <- function (x,m,tabmu0,p0,w=NULL) {
# input:
#   x   ordinal data set [n,d]
#   m   nb of modalities [d]
#   tabmu0  all starting mu0 to try cell(d,1)
#   p0  starting p0 to try for each starting mu0 [d]
#   w   weights [n] OPTIONAL
# output:
#   mu_ml  ml true modality {1,...,m} [d]
#   p_ml   ml P(sait) [0,1] [d]
#   ml     ml value

n = nrow(x)
d = ncol(x)
  
#  all weights w to 1 of missing 
if (is.null(w)) w = rep(1,n)

ntot = sum(w)
ml = -Inf

# each dimension is estimated independently
mu_ml = rep(0,d)
p_ml = rep(0,d)
ml = 0
for (id in 1:d){
  res=  ordiem(x[,id],m[id],tabmu0[id,],p0[id],w)
  mu_ml[id]=res$mu_ml
  p_ml[id]=res$p_ml
  mlid = res$ml
  ml = ml + mlid
}
return(list(mu_ml=mu_ml,p_ml=p_ml,ml=ml))
}


############################################################################# 
pej <- function(ej,j,m,mu,p,z1tozjm1=NULL){
# Proba ej (unconditional) with a hierarchical random dichotomic model.
# input:
#   ej  interval ej+1 [binf,bsup]
#   j   dichotomical step
#   m   nb of modalities
#   mu  true modality {1,...,m}
#   p   P(sait) [0,1]
#   z1tozjm1    values {0,1} of z_1 to z_{j-1} OPTIONAL
#               1: z is known and is 1
#               2: z is unknown and is 0
# output:
#   pej

# j == 1
if (j == 1){ 
proba = 1
return(proba)
}

# if ej is not an interval then build [ej,ej]
if (length(ej) == 1) {ej = c(ej,ej)}

# create zjm1toz1
if (is.null(z1tozjm1)){
z1tozjm1 = matrix(0,1,j-1) # 0 at place j <=> unknown value for z_j
}
z1tozjm2 = z1tozjm1[1:(length(z1tozjm1)-1)];
zjm1 = z1tozjm1[length(z1tozjm1)];


# version nouvelle plus rapide (exple de gain : 17* pour m=5, 374* pour m=6)
# ----------------------------
# idee de la mehode : on coupe les branches de proba nulle, cad ou ejm1 pas
# inclus dans ej

# j > 1
if (zjm1) { # zjm1 is known   ####### ATTENTION : pas sur que ca marche sous R
proba = 0
tabint = allej(j-1,m)
if (!is.matrix(tabint)) {tabint = t(as.matrix(allej(j-1,m)))}
nbtabint = length(tabint[,1])
for (i in 1:nbtabint){
ejm1 = tabint[i,] # ej minus 1
if ((ej[1] >= ejm1[1]) && (ej[2] <= ejm1[2])){ # pour accelerer, il faut verifier ejm1 est inclus dans ej !
  proba = proba + pejp1zj1_ej(ej,ejm1,mu,p) * pej(ejm1,j-1,m,mu,p,z1tozjm2);
}
}
}
else # zjm1 is unknown    
{proba = 0
tabint = allej(j-1,m)
if (!is.matrix(tabint)) {tabint = t(as.matrix(allej(j-1,m)))}
nbtabint = length(tabint[,1])
for (i in 1:nbtabint){
ejm1 = tabint[i,]# ej minus 1
if ((ej[1] >= ejm1[1]) && (ej[2] <= ejm1[2])){ # pour accelerer, il faut verifier ejm1 est inclus dans ej !
   proba = proba + pejp1_ej(ej,ejm1,mu,p) * pej(ejm1,j-1,m,mu,p,z1tozjm2)
}
}
} # zjm1
return(proba)
}

############################################################################# 
allej <- function(j,m){
# tabint = allej(j,m)
# All possible intervals ej with a hierarchical random dichotomic model.
# input:
#  j   dichotomical step
#  m   nb of modalities
# output:
#  tabint  table of intervals [binf,bsup]

# cas j = 1
if (j == 1) {
tabint = c(1,m)
return(tabint)
}

if (j >m) {
  tabint = c()
  return(tabint)
}

# cas j > 1
tabint = c()
for (sizeej in 1:(m-j+1)){
for (binf in 1:(m-sizeej+1)){
  bsup = binf + sizeej-1
  tabint = rbind(tabint,c(binf,bsup))
}
}
return(tabint)
}

############################################################################# 
pejp1zj1_ej <-function(ejp1,ej,mu,p){
# Proba (ej+1,zj=1) cond. to ej with a hierarchical random dichotomic model.
# input:
#   ejp1 interval ej+1 [binf,bsup]
#   ej  interval ej [binf,bsup]
#   mu  true modality {1,...,m}
#   p   P(sait) [0,1]
# output:
#   pejp1_yjej

proba = 0
for (yj in ej[1]:ej[2]){
proba = proba + pejp1zj1_yjej(ejp1,yj,ej,mu,p)
}
proba = proba /  (ej[2]-ej[1]+1)
return(proba)
}

############################################################################# 
pejp1zj1_yjej<-function(ejp1,yj,ej,mu,p){
# Proba (ej+1,zj=1) cond. to ej and yj with a hierarchical random dichotomic model.
# input:
#   ejp1 interval ej+1 [binf,bsup]
#   yj  integer in ej
#   ej  interval ej [binf,bsup]
#   mu  true modality {1,...,m}
#   p   P(sait) [0,1]
# output:
#   pejp1_yjej

ejminus = c(ej[1],(yj-1))
ejequal = c(yj,yj)
ejplus = c((yj+1),ej[2])

# pejp1_yjejzj1
if (ejminus[1] > ejminus[2]){
dmuejminus = Inf
}
else{
dmuejminus = min(abs(mu-ejminus))
}

if (ejplus[1] > ejplus[2]){
dmuejplus = Inf
}
else
{
  dmuejplus = min(abs(mu-ejplus))
}

dmuejequal = min(abs(mu-ejequal))

dmuejp1 = min(abs(mu-ejp1))


if (dmuejp1 == min(c(dmuejminus,dmuejequal,dmuejplus)) && (all(ejp1==ejminus) |  all(ejp1==ejequal) |  all(ejp1==ejplus)))
{pejp1_yjejzj1 = 1}
else
{pejp1_yjejzj1 = 0}

# pejp1zj1_yjej
proba = p*pejp1_yjejzj1
return(proba)
}

############################################################################# 

pejp1_ej<-function(ejp1,ej,mu,p){

# Proba ej+1 cond. to ej with a hierarchical random dichotomic model.
# input:
#   ejp1 interval ej+1 [binf,bsup]
#   ej  interval ej [binf,bsup]
#   mu  true modality {1,...,m}
#   p   P(sait) [0,1]
# output:
#   pejp1_yjej

proba = 0;

# attention : on suppose toujours ejp1 inclus dans ej (pas verifie ici)

if (ejp1[2]== ejp1[1]){ # |ejp1|=1
  if ((ejp1[2] < ej[2]) && (ejp1[1] > ej[1])){ # ejp1 ne touche aucun bord
    allyj = ejp1[1]
  }
  else{
    if (ejp1[2] < ej[2]){ # ejp1 ne touche pas le bord droit
      allyj = c(ejp1[1],ejp1[1]+1)}
    else{  # ejp1 ne touche pas le bord gauche
      allyj = c(ejp1[1]-1,ejp1[1])
    }
  }
}
else{ # |ejp1|>1
  if (ejp1[2] < ej[2]){ # ejp1 ne touche pas le bord droit
    allyj = ejp1[2]+1}
  else { # ejp1 ne touche pas le bord gauche
    allyj = ejp1[1]-1
  }
}

for (yj in allyj){proba = proba + pejp1_yjej(ejp1,yj,ej,mu,p) * pyj_ej(yj,ej)}

return(proba)
}
  
############################################################################# 
pejp1_yjej<-function(ejp1,yj,ej,mu,p){
 
  # Proba ej+1 cond. to ej and yj with a hierarchical random dichotomic model.
  # input:
  #   ejp1 interval ej+1 [binf,bsup]
  #   yj  integer in ej
  #   ej  interval ej [binf,bsup]
  #   mu  true modality {1,...,m}
  #   p   P(sait) [0,1]
  # output:
  #   pejp1_yjej
  
  ejminus = c(ej[1],(yj-1))
  ejequal = c(yj,yj)
  ejplus = c((yj+1),ej[2])
  
  # pejp1_yjejzj0
  if (all(ejp1==ejminus) |  all(ejp1==ejequal) |  all(ejp1==ejplus)){
  pejp1_yjejzj0 = (ejp1[2]-ejp1[1]+1) / (ej[2]-ej[1]+1)
  }
  else{
    pejp1_yjejzj0 = 0
  }
  
  # pejp1_yjejzj1
  if (ejminus[1] > ejminus[2]){
  dmuejminus = Inf
  }
  else
  {
    dmuejminus = min(abs(mu-ejminus))
  }
  
  if (ejplus[1] > ejplus[2]){
    dmuejplus = Inf
  }
  else{
    dmuejplus = min(abs(mu-ejplus))
  }
  
  dmuejequal = min(abs(mu-ejequal))
  
  dmuejp1 = min(abs(mu-ejp1))
  
  if (dmuejp1 == min(c(dmuejminus,dmuejequal,dmuejplus)) && (all(ejp1==ejminus) |  all(ejp1==ejequal) |  all(ejp1==ejplus))){
  pejp1_yjejzj1 = 1
  }
  else{
    pejp1_yjejzj1 = 0
  }
  
  # pejp1_yjej
  proba = p*pejp1_yjejzj1 + (1-p)*pejp1_yjejzj0;
  
  return(proba)
}

############################################################################# 
pyj_ej<-function(yj,ej){
  
# proba = pyj_ej(yj,ej)
# Proba yj cond. to ej with a hierarchical random dichotomic model.
# input:
#   yj  integer
#   ej  interval [binf,bsup]
# output:
#   pyj_ej
  
  if ((ej[1] <= yj) && (yj <= ej[2])){
  proba = 1  / (ej[2]-ej[1]+1);
  }
  else{
    proba = 0
  }
 return(proba)
}
  
  
