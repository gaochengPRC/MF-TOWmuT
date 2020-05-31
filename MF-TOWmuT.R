## R Program for TOWmuT/MF-TOWmuT 
## Created by Cheng Gao
## Date: 05/31/2020
## Email: chenggao@mtu.edu
## Department of Mathematical Sciences
## Michigan Technological University
## Copyright reserved. 
## Claim: This program can only be used for academic research. 
## Note: Further details are available only upon request, please feel free to contact me.
## Ref: Cheng Gao et al, MF-TOWmuT: Testing an Optimally Weighted Combination of Common and Rare Variants with Multiple Traits Using
## Family Data (submitted to Genetic Epidemiology)

#----------------------------------------------------------------------------------
#----------------------- return test statistic ------------------------------------
## num.row: number of individuals 
## oc: 0, without covariates; otherwise, with covariates
## sh: if sh == i, means the i-th repetition
## geno.dat/pheno.dat: centralized genotype/phenotype data matrix
## z1/z2: two covariates, you can easily adapt the code to any number of covariates
#----------------------------------------------------------------------------------

get.stat <- function(oc,sh,geno.dat,pheno.dat,z1,z2,num.row) {
  
  if(oc == 0){
    geno.resid <- geno.dat
    pheno.resid <- pheno.dat
  } else{
    
    geno.resid <- resid(lm(geno.dat ~ cbind(z1, z2)))
    pheno.resid <- resid(lm(pheno.dat ~ cbind(z1, z2)))
  }
   
  if (sh == 0) {  # sh = 0 means original data
    
    pheno.resid <-  pheno.resid
    
  } else {
    
    shuf <- sample(1:num.row, num.row, replace = FALSE)
    
    pheno.resid <-  pheno.resid[shuf, ]
  }
  
  geno.resid.ind <- geno.resid
  pheno.resid.ind <- pheno.resid
  
  ## MF-TOWMuT
  geno.resid <- sig.star.inv.sqrt %*% geno.resid # difference
  A.inv <- ginv(t(geno.resid) %*% geno.resid / num.row )  
  B <- geno.resid %*% A.inv %*% t(geno.resid) # not change in permutation 
  pheno.resid.shuf <- sig.star.inv.sqrt %*% pheno.resid # difference
  
  Y.var.inv <- ginv(t(pheno.resid.shuf) %*% pheno.resid.shuf) # note: when phenotypes are shuffled, Y.var.inv will change !!
  DD <- Y.var.inv %*% t(pheno.resid.shuf) %*% B %*% pheno.resid.shuf
  test.stat <- eigen(DD)$values[1]
  
  ## TOWMuT
  A.inv.ind <-  ginv(t(geno.resid.ind) %*% geno.resid.ind / num.row  )
  B.ind <- geno.resid.ind %*% A.inv.ind %*% t(geno.resid.ind) # not change in permutation 
  pheno.resid.ind.shuf <- pheno.resid.ind
  Y.var.inv.ind <- ginv(t(pheno.resid.ind.shuf) %*% pheno.resid.ind.shuf) # note: when phenotypes are shuffled, Y.var.inv will change !!
  DD.ind <- Y.var.inv.ind %*% t(pheno.resid.ind.shuf) %*% B.ind %*% pheno.resid.ind.shuf
  test.stat.ind <- eigen(DD.ind)$values[1]
  
  return(list(test.stat.ind = test.stat.ind, test.stat = test.stat))
  
}

#-----------------------------------------------------------------------
#---------------- return p-value ---------------------------------------
## sig.star.inv.sqrt: square root of inverse of kinship matrix
## num.permu: the number of permutations
#----------------------------------------------------------------------- 

get.p <- function(oc,geno.dat,pheno.dat,z1,z2,sig.star.inv.sqrt) {
  
  # obtain the test statistic
  tv <- lapply(0:num.permu, FUN = function(i) {matrix(unlist(get.stat(oc,i,geno.dat,pheno.dat,z1,z2,num.row)),ncol=2)})
  tv <- do.call('rbind',tv)
  tv.ind <- tv[,1]
  p.value.ind <- mean(tv.ind[-1] > tv.ind[1]) # calculate p-value for TOWMuT
  
  tv.fam <- tv[,2]
  p.value.fam <- mean(tv.fam[-1] > tv.fam[1]) # calculate p-value for MF-TOWMuT
  
  return(list(p.value.ind = p.value.ind, p.value.fam = p.value.fam))
}









