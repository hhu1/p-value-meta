brown = function(m) {
  # implementation of method 4 (Brown method) in: Accuracy Evaluation of the Unified P-Value from Combining Correlated P-Values
  # the input is a matrix, whose columns correspond to interested phenotypes, and rows corresponds to genes. The values in the matrix are the p-value of association between the gene and each phenotype, based on any rare-variant association test. 
  
  m[is.na(m)] = 1
  m[m==0] = 1e-20
  
  num_genes = nrow(m)
  num_pts = ncol(m)
  
  cor_m = cor(log(m))
  cov_m = 3.263 * cor_m + 0.710 * (cor_m^2) + 0.027* (cor_m^3)
  
  sum_term = 0 # this is the common sum term in c and f calculation
  for (i in 1:(num_pts-1)) {
    for (j in (i+1):num_pts) {
      sum_term = sum_term + cov_m[i,j]
    } 
  }
    
  c = (2* num_pts + sum_term) / (2* num_pts)
  f = 4 * num_pts^2 / (2*num_pts + sum_term)
  
  ps = numeric()
  for (i in 1:num_genes) {
    tau = sum(-2* log(m[i,]))
    ps[i] = pchisq(tau/c, f, lower.tail=F)
  }
  return(ps) # returning a vector of p-values, one for each gene
}

rvat_chi2 = function(m) {
  # the input is a matrix, whose columns correspond to interested phenotypes, and rows corresponds to genes. The values in the matrix are the p-value of association between the gene and each phenotype, based on any rare-variant association test. 
  m[is.na(m)] = 1
  m[m==0] = 1e-16
  
  #signs_m =  ((m>0)-0.5)*2
  #abs_m = abs(m)
  #matrix = qnorm(abs_m/2, lower.tail=F) * signs_m
  matrix= qnorm(m, lower.tail=F)
  
  num_genes = nrow(m)
  num_pts = ncol(m)
  V_inv = solve(cor(matrix))
  ps = numeric()
  
  for (i in 1:nrow(matrix)) {
    t = matrix[i,]
    test_stat = (t(t) %*% V_inv %*% t)[1,1]
    p_value = pchisq(test_stat, num_pts, lower.tail =F )
    ps[i] = p_value
  }
  return(ps)
  # return a vector of p-values corresponding to each gene in the input matrix (rows of the matrix). The p-value tests for the association between each gene and the multi-traits.
}

# example: calculate the meta-analysis p-values for 300 genes across 8 traits. 
m =matrix(runif(300*8), nrow=300)
p1 = rvat_chi2(m)
p2 = brown(m)
