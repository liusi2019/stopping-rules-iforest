## file for general case with parameter setting for target one
## hard coded parameters:
## alpha, 1 - the confidence level, this part is actually hard coded in the code
## %% fix the k value to be 5 here







## data generation done
## add three modifications -- done with all three


#####################################################
##########################preceding part
##################################################
### set parameters
evaluation <- function(a){ 
  library(MASS)
  rho = 0.05
  alpha = 0.05
  #proportion = 0.01
  
  ### set the structure of data 
  
  vsize <- c(100, 1000, 10000)
  vpro <- c(5, 10, 20, 30)
  #start from 10 instead because the results for 5 is already out there
  #vpro <- c(10, 20, 30)
  vdist <- c(3, 5, 8, 10)
  vd <- c(10, 20, 40, 100)
  veta <- c(0.10, 0.20)
  vtype <- c(0, 1)
  
  #---------------------------------------------------------------------------
  for(ind_round in 1:14){
    if(ind_round%in%(1:3)){
      n = vsize[ind_round]
      j = 1
      dist = vdist[2]
      d = vd[1]
      eta = veta[2]
      a_type = vtype[1]
    }else if(ind_round%in%(4:6)){
      n = vsize[1]
      j = round(ind_round - 2)
      dist = vdist[2]
      d = vd[1]
      eta = veta[2]
      a_type = vtype[1]
    }else if(ind_round%in%(7:9)){
      n = vsize[1]
      j = 1
      if(ind_round == 7){
        dist = vdist[ind_round - 6]
      }else{
        dist = vdist[ind_round - 5]
      }
      d = vd[1]
      eta = veta[2]
      a_type = vtype[1]
    }else if(ind_round%in%(10:12)){
      n = vsize[1]
      j = 1
      dist = vdist[2]
      d = vd[ind_round - 8]
      eta = veta[2]
      a_type = vtype[1]
    }else if(ind_round == 13){
      n = vsize[1]
      j = 1
      dist = vdist[2]
      d = vd[1]
      eta = veta[ind_round - 12]
      a_type = vtype[1]
    }else if(ind_round == 14){
      n = vsize[1]
      j = 1
      dist = vdist[2]
      d = vd[1]
      eta = veta[2]
      a_type = vtype[ind_round - 12]
    }


  #---------------------------------------------------------------------------
  
  

  a_dim = round(d * eta)
  
  proportion = vpro[j]/100
  #k = round(n*proportion)
  #using fixed k = 5
  k = 5
  
  rseed <- as.numeric(unlist(read.csv("general_random_seed.csv", header = FALSE)))
  set.seed(rseed[a])
  
  dat <- matrix(ncol = d, nrow = 2*n)
  vmat = matrix(0, ncol = d, nrow = d)
  diag(vmat) = 1
  if (a_type == 0){
  for(i in (1:round(n*proportion))){
    center = rep(0, d)
    center[sample(d, a_dim,replace = F)] = dist
    dat[i,] = mvrnorm(1, center, vmat)
  } 
  }else{
    a_index = sample(d, a_dim, replace = F)
    for(i in (1:round(n*proportion))){
      center = rep(0, d)
      center[a_index] = dist
      dat[i,] = mvrnorm(1, center, vmat)
    } 
  }
  
  nnrow <- round(n * (2 - proportion))
  rvec <- rnorm(d * nnrow, 0, 1)
  
  dat[(n*proportion+1):(2*n),] <- base::matrix(rvec, nrow = round(n*(2-proportion)), ncol = d)
  data2 <- dat[(1:n),]
  data1 <- dat[((n+1):(2*n)),]
  

  nametag <- paste("pro",vpro[j],"_n",n,"_d", d, "_dist", dist, "_eta", eta, "_type", a_type, sep = "")
  #write.csv(data1, file = paste("data1_",vpro[j],"_",n,"_", a, ".csv", sep = ""), row.names = FALSE)
  write.csv(data2, file = paste("synthetic/permutation/", nametag, "/data2_pro",vpro[j],"_n",n,"_d", d, "_dist", dist, "_eta", eta, "_type", a_type, "_", a, ".csv", sep = ""), row.names = FALSE)
  ### get c(n)----this is the c(n) used in the isolation forest paper for approximating the average path length
  c_n <- function(n){
    H = log(n - 1) + 0.5772156649
    return(2*H - 2*(n-1)/n)
  }
  ### function for anomaly score----given average depth h, calculate the corresponding anormaly score
  ans <- function(h, c_n){
    return(2^{-(h/c_n)})
  }
  ### grow isolation trees
  ### here average depths from 10000 trees are used as the "ground truth". 
  ### may not be able to do this if we have a larger data set, because iforest dies a lot when trying to calculate 10000 trees on 10000 points.
  
  #vseed <- rep(0, 50)
  #for (i in 1:50){
  #  vseed[i] = sample(0:32767, 1)    
  #}
  
  time_list = rep("b", 50)
  seed_list = rep(0, 50)

  op <- options(digits.secs = 30)
  total_seed <- as.numeric(unlist(read.csv("random_seed_nature.csv", header = FALSE)))
  seed_to_use <- total_seed[((a-1)* 50 +1) : (a * 50) ]
  curr_seed = seed_to_use[1]
  seed_list[1] = curr_seed
  time_list[1] = Sys.time()
  system(paste('./iforest','-i', paste('synthetic/permutation/', nametag, '/data2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = "") ,'-o',paste('synthetic/permutation/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""),'-t 200 -s 0 -p -x', paste('synthetic/permutation/', nametag, '/data2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""),'-z', curr_seed, '-g'), wait = TRUE)
  depth2 <- read.csv(paste('synthetic/permutation/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), header = FALSE)
  depth2 <- depth2[, 2:201]
  system(paste('rm','-f',paste('synthetic/permutation/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = "")))
  
  for (i in 2:50){
    curr_seed = seed_to_use[i]
    seed_list[i] = curr_seed    
    time_list[i] = Sys.time()
    system(paste('./iforest','-i', paste('synthetic/permutation/', nametag, '/data2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = "") ,'-o',paste('synthetic/permutation/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""),'-t 200 -s 0 -p -x', paste('synthetic/permutation/', nametag, '/data2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), '-z', curr_seed, '-g'), wait = TRUE)
    depth2_mid <- read.csv(paste('synthetic/permutation/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), header = FALSE)
    depth2_mid <- depth2_mid[, 2:201]
    depth2 = cbind(depth2, depth2_mid)
    system(paste('rm','-f',paste('synthetic/permutation/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = "")))
  }
  
  #system(paste('rm','-f',file = paste("data1_",vpro[j],"_",n,"_", a, ".csv", sep = ""), row.names = FALSE))
  #system(paste('rm','-f',file = paste("m_truth_toy/data2_",vpro[j],"_",n,"_", a, ".csv", sep = ""), row.names = FALSE))
  
  ### order the rows(each row stands for depths for all points from one tree) according to the anomaly score, in the descreasing order
  
  avgdepth2 <- apply(depth2, 1, mean)
  depth2$score <- ans(avgdepth2, c_n(n))
  
  #write.csv(depth2, file = paste("m_truth_toy_00000/big_depth_",vpro[j],"_",n,"_", a, ".csv", sep = ""), row.names = FALSE)
  depth_order = depth2[order(-depth2$score),]
  
  ################################################################
  ##### Here comes formulas and stopping rules
  ################################################################
  ### construct a data frame to store the result
  #simdat<-data.frame(array(dim=c(1,17)))
  #names(simdat)<-c('targetsize','nt_tCI', 'pre_tCI','npa_tCI','nps_tCI', 'time_tCI', 
  #                 'nt_Ho','pre_Ho','npa_Ho','nps_Ho','time_Ho', 'nt_EBB', 'pre_EBB','npa_EBB', 'nps_EBB','time_EBB','truth')
  
  simdat<-data.frame(array(dim=c(1,18)))
  names(simdat)<-c('targetsize','nt_tCI', 'pre_tCI','npa_tCI','nps_tCI', 'time_tCI', 
                   'nt_Ho','pre_Ho','npa_Ho','nps_Ho','time_Ho', 'nt_EBB', 'pre_EBB','npa_EBB', 'nps_EBB','time_EBB','truth', 'AUC')
  
  library(pROC)
  
  label = c(rep(1, round(n* proportion)), rep(0, round(n*(1 - proportion))))
  ROC1 = roc(label, ans(avgdepth2, c_n(n)), direction = "<")
  AUC1 = auc(ROC1)
  simdat$AUC = AUC1
  
  ### set the relaxd set that we target on
  names <- row.names(depth_order[which(depth_order$score >= (1-rho)*depth_order$score[k]),]) 
  ### record the size of the relaxed target set
  simdat$targetsize <- length(names)
  ###########################################################################
  ######next we have code for the 4 other confidence interval methods
  ######one trick here is to make use of the fact that in R, if you take a subset of 
  ######rows from a dataframe, the rows in the new dataframe will keep the original rownames.
  ########################################################
  #####here comes normal approach
  ##### ininitialize S_left
  ## insert the timing part
  ptm <- proc.time()
  S_left = seq(1,n, by = 1)
  t = 0
  ##### do the loop from g = 10 to g = 100
  ##### which corresponds to the number of trees from 10^2 = 100 to 100^2 = 10000
  for(g in 1:100){
    h = 100 * g
    ##### select the depths from the first h trees
    depth_h = depth2[,1:(h)]
    ##### from the result above, select the rows corresponding to S_left
    depth_h = depth_h[S_left,]
    ##### calculate the average depths
    depth_h$avg = apply(depth_h[,1:h], 1, mean)
    ##### initialize upper bounds, lower bounds and half widths to 0
    depth_h$upper = rep(0, nrow(depth_h))
    depth_h$lower = rep(0, nrow(depth_h))
    depth_h$halfwidth = rep(0, nrow(depth_h))
    ##### calculate half width for each point to construct CI
    ##### need to adjust for t--the number of CIs we have already calculated
    for(i in 1:nrow(depth_h)){ 
      ## here we use "finite" equal Bonferroni corrections for arounds from 10 to 100
      #depth_h$halfwidth[i] = qt(1-0.025/(nrow(depth_h)* 91), h-1)*sd(depth_h[i,1:h])/sqrt(h)
      depth_h$halfwidth[i] = qt(1-0.025/(nrow(depth_h)* 100), h-1)*sd(depth_h[i,1:h])/sqrt(h)
    }
    ##### construct CI by calculating upper bound and lower bound
    ##### for calculating lower bound, if average - halfwidth <0, just take 0
    depth_h$upper = depth_h$avg + depth_h$halfwidth
    depth_h$lower = pmax(depth_h$avg - depth_h$halfwidth, rep(0, nrow(depth_h)))
    ##### the following steps are for getting D_lower and D_upper
    ##### the following two lines will order the row names of the points by the upper bound or lower bound
    order_upper = order(depth_h$upper)
    order_lower = order(depth_h$lower)
    ##### get D_lower
    ##### --> this part(two lines) is changed, the orignial version is, select the points in S_left with upper bound greater than the value of the (k-1)th smallest lower bound
    ##### note that, we change the previous version to the new version
    ##### -->d_low_k_1 = depth_h$lower[order_lower[k-1]]
    ##### -->D_lower = depth_h$upper> d_low_k_1
    ##### now we select the points in S_left with upper bound greater than or equal to the value of the kth smallest lower bound
    d_low_k_1 = depth_h$lower[order_lower[k]]
    D_lower = depth_h$upper>= d_low_k_1
    ##### get D_upper
    ##### --> this part (two lines) is changed, the original version is, select the points in S_left with lower bound smaller than the value of the (k+1)th smallest upper bound
    ##### --> also is the (size(S_left) - k)th largest upper bound
    ##### -->d_up_k_plus_1 = depth_h$upper[order_upper[k+1]]
    ##### -->D_upper = depth_h$lower < d_up_k_plus_1
    ##### now we select the points in S_left with lower bound smaller than or equal to the value of the (k)th smallest upper bound
    d_up_k_1 = depth_h$upper[order_upper[k]]
    D_upper = depth_h$lower <= d_up_k_1
    ##### get C, which is the intersection of D_lower and D_upper
    C = D_lower&D_upper
    ##### calculate the anomaly score of the smallest lower bound among points in C
    s_upper_d_k = ans(min(depth_h$lower[C]),c_n(n))
    ##### calculate the anomaly score of upper bound for all points in S_left
    depth_h$ans = ans(depth_h$upper,c_n(n))
    ##### choose the points to output in A_hat
    A_hat<-row.names(depth_h[depth_h$ans>= (1-rho)*s_upper_d_k,])
    ##### get the kth smallest upper bound
    d_upper_k = depth_h$upper[order_upper[k]]
    ##### get what shoud be left for the new S_left after this round
    ##### get rid of all the points who have lower bound greater than the kth smallest upper bound
    S_left =row.names(depth_h[depth_h$lower<= d_upper_k,])
    ## what to output if we stop because the size of A_hat is no less than k
    if(length(A_hat)>=k){
      ## the number of trees
      simdat$nt_tCI = h
      ## the size of A_hat we output
      simdat$npa_tCI = length(A_hat)
      ## this nps is the size of S_left, if we stop because the size of S_left shrinks to no greater than k
      ## and 0 if we stop because the size of A_hat reaches no less than k
      simdat$nps_tCI = 0
      ## the precision of our output
      simdat$pre_tCI= as.numeric(all(A_hat %in% names))
      break
    }
    ## what to output if we stop because the size of S_left shrinks to no greater than k
    if(length(S_left)<=k){
      simdat$nt_tCI = h
      simdat$npa_tCI = 0
      simdat$nps_tCI = length(S_left)
      simdat$pre_tCI= as.numeric(all(S_left %in% names)) 
      break
    }
    ### this is an additional stopping rule to gaurantee the algorithm will stop and give output
    ### if we achieve neither length(A_hat)>= k nor length(S_left)<= k before we have g^2 = 100^2 = 10000 trees
    if(g == 100){
      simdat$nt_tCI = 10000
      ### I use 10000 for npa and nps in this case, just to indicate why we stop
      simdat$npa_tCI = 10000
      simdat$nps_tCI = 10000
      depth_add = depth2[,1:(h)]
      depth_add$avg = apply(depth_add[,1:h], 1, mean)
      depth_add_order = depth_add[order(depth_add$avg),]
      ### in this case we output the first k points according to average depths from 10000 trees
      #add_names = row.names(depth_add_order[1:k,])
      add_names = row.names(depth_add_order[depth_add_order$avg<= depth_add_order$avg[k],])
      simdat$pre_tCI = as.numeric(all(add_names %in% names)) 
    }
  }
  simdat$time_tCI = (proc.time() - ptm)[['elapsed']]
  
  ########################################################
  #####here comes Hoeffding bound approach
  ##### first get 100 trees to calculate the range for each point
  ptm = proc.time()
  depth_100 = depth2[,1:100]
  maxdepth = apply(depth_100, 1, max)
  S_left = seq(1,n, by = 1)
  t = 0
  for(g in 1:100){
    h = 100 * g
    depth_h = depth2[,1:(h)]
    ### put the maximum depths from above into the data frame
    depth_h$maxdepth = maxdepth
    depth_h = depth_h[S_left,]
    depth_h$avg = apply(depth_h[,1:h], 1, mean)
    depth_h$upper = rep(0, nrow(depth_h))
    depth_h$lower = rep(0, nrow(depth_h))
    depth_h$halfwidth = rep(0, nrow(depth_h))
    ### need to adjust for t here
    for(i in 1:nrow(depth_h)){ 
      #depth_h$halfwidth[i] = depth_h$maxdepth[i] * sqrt(log(2*nrow(depth_h)* 91/0.05)/(2*h))
      depth_h$halfwidth[i] = depth_h$maxdepth[i] * sqrt(log(2*nrow(depth_h)* 100/0.05)/(2*h))
    }
    depth_h$upper = depth_h$avg + depth_h$halfwidth
    depth_h$lower = pmax(depth_h$avg - depth_h$halfwidth, rep(0, nrow(depth_h)))
    ## get D_lower and D_upper
    order_upper = order(depth_h$upper)
    order_lower = order(depth_h$lower)
    ##### get D_lower
    ##### --> this part(two lines) is changed, the orignial version is, select the points in S_left with upper bound greater than the value of the (k-1)th smallest lower bound
    ##### note that, we change the previous version to the new version
    ##### -->d_low_k_1 = depth_h$lower[order_lower[k-1]]
    ##### -->D_lower = depth_h$upper> d_low_k_1
    ##### now we select the points in S_left with upper bound greater than or equal to the value of the kth smallest lower bound
    d_low_k_1 = depth_h$lower[order_lower[k]]
    D_lower = depth_h$upper>= d_low_k_1
    ##### get D_upper
    ##### --> this part (two lines) is changed, the original version is, select the points in S_left with lower bound smaller than the value of the (k+1)th smallest upper bound
    ##### --> also is the (size(S_left) - k)th largest upper bound
    ##### -->d_up_k_plus_1 = depth_h$upper[order_upper[k+1]]
    ##### -->D_upper = depth_h$lower < d_up_k_plus_1
    ##### now we select the points in S_left with lower bound smaller than or equal to the value of the (k)th smallest upper bound
    d_up_k_1 = depth_h$upper[order_upper[k]]
    D_upper = depth_h$lower <= d_up_k_1
    #######get C
    C = D_lower&D_upper
    ##get s_upper_d_k
    s_upper_d_k = ans(min(depth_h$lower[C]),c_n(n))
    ##get the set to output, A_hat
    depth_h$ans = ans(depth_h$upper,c_n(n))
    A_hat<-row.names(depth_h[depth_h$ans>= (1-rho)*s_upper_d_k,])
    ##get d_upper_k
    d_upper_k = depth_h$upper[order_upper[k]]
    ##get what shoud be left after one round, the new S_left
    S_left =row.names(depth_h[depth_h$lower<= d_upper_k,])
    ##output
    if(length(A_hat)>=k){
      simdat$nt_Ho = h
      simdat$npa_Ho = length(A_hat)
      simdat$nps_Ho = 0
      simdat$pre_Ho= as.numeric(all(A_hat %in% names))
      break
    }
    if(length(S_left)<=k){
      simdat$nt_Ho = h
      simdat$npa_Ho = 0
      simdat$nps_Ho = length(S_left)
      simdat$pre_Ho= as.numeric(all(S_left %in% names)) 
      break
    }
    ### tadditional stopping rule to gaurantee the algorithm will stop and give output
    if(g == 100){
      simdat$nt_Ho = 10000
      simdat$npa_Ho = 10000
      simdat$nps_Ho = 10000
      depth_add = depth2[,1:(h)]
      depth_add$avg = apply(depth_add[,1:h], 1, mean)
      depth_add_order = depth_add[order(depth_add$avg),]
      #add_names = row.names(depth_add_order[1:k,])
      add_names = row.names(depth_add_order[depth_add_order$avg<= depth_add_order$avg[k],])
      simdat$pre_Ho = as.numeric(all(add_names %in% names)) 
    }
  }
  simdat$time_Ho = (proc.time() - ptm)[['elapsed']]
  
  ########################################################
  #####here comes Bernstein bound approach
  ptm = proc.time()
  depth_100 = depth2[,1:100]
  maxdepth = apply(depth_100, 1, max)
  S_left = seq(1,n, by = 1)
  t = 0
  for(g in 1:100){
    h = 100 * g
    depth_h = depth2[,1:(h)]
    depth_h$maxdepth = maxdepth
    depth_h = depth_h[S_left,]
    depth_h$avg = apply(depth_h[,1:h], 1, mean)
    depth_h$upper = rep(0, nrow(depth_h))
    depth_h$lower = rep(0, nrow(depth_h))
    depth_h$halfwidth = rep(0, nrow(depth_h))
    ## need to adjust for t here
    for(i in 1:nrow(depth_h)){ 
      #depth_h$halfwidth[i] = depth_h$maxdepth[i]*(sqrt(2*(var(as.numeric(depth_h[i,1:h]))/(depth_h$maxdepth[i]^2))*log(4*nrow(depth_h)* 91/0.05)/h) + 7*log(4*nrow(depth_h)* 91/0.05)/(3*(h-1)))
      depth_h$halfwidth[i] = depth_h$maxdepth[i]*(sqrt(2*((stats::var(as.numeric(depth_h[i,1:h])))/(depth_h$maxdepth[i]^2))*log(4*nrow(depth_h)* 100/0.05)/h) + 7*log(4*nrow(depth_h)* 100/0.05)/(3*(h-1)))
    }
    depth_h$upper = depth_h$avg + depth_h$halfwidth
    depth_h$lower = pmax(depth_h$avg - depth_h$halfwidth, rep(0, nrow(depth_h)))
    ## get D_lower and D_upper
    order_upper = order(depth_h$upper)
    order_lower = order(depth_h$lower)
    ##### get D_lower
    ##### --> this part(two lines) is changed, the orignial version is, select the points in S_left with upper bound greater than the value of the (k-1)th smallest lower bound
    ##### note that, we change the previous version to the new version
    ##### -->d_low_k_1 = depth_h$lower[order_lower[k-1]]
    ##### -->D_lower = depth_h$upper> d_low_k_1
    ##### now we select the points in S_left with upper bound greater than or equal to the value of the kth smallest lower bound
    d_low_k_1 = depth_h$lower[order_lower[k]]
    D_lower = depth_h$upper>= d_low_k_1
    ##### get D_upper
    ##### --> this part (two lines) is changed, the original version is, select the points in S_left with lower bound smaller than the value of the (k+1)th smallest upper bound
    ##### --> also is the (size(S_left) - k)th largest upper bound
    ##### -->d_up_k_plus_1 = depth_h$upper[order_upper[k+1]]
    ##### -->D_upper = depth_h$lower < d_up_k_plus_1
    ##### now we select the points in S_left with lower bound smaller than or equal to the value of the (k)th smallest upper bound
    d_up_k_1 = depth_h$upper[order_upper[k]]
    D_upper = depth_h$lower <= d_up_k_1
    #######get C
    C = D_lower&D_upper
    ##get s_upper_d_k
    s_upper_d_k = ans(min(depth_h$lower[C]),c_n(n))
    ##get the set to output, A_hat
    depth_h$ans = ans(depth_h$upper,c_n(n))
    A_hat<-row.names(depth_h[depth_h$ans>= (1-rho)*s_upper_d_k,])
    ##get d_upper_k
    d_upper_k = depth_h$upper[order_upper[k]]
    ##get what shoud be left after one round, the new S_left
    S_left =row.names(depth_h[depth_h$lower<= d_upper_k,])
    ##output
    if(length(A_hat)>=k){
      simdat$nt_EBB = h
      simdat$npa_EBB = length(A_hat)
      simdat$nps_EBB = 0
      simdat$pre_EBB= as.numeric(all(A_hat %in% names))
      break
    }
    if(length(S_left)<=k){
      simdat$nt_EBB = h
      simdat$npa_EBB = 0
      simdat$nps_EBB = length(S_left)
      simdat$pre_EBB= as.numeric(all(S_left %in% names)) 
      break
    }
    ### additional stopping rule to gaurantee the algorithm will stop and give output
    if(g == 100){
      simdat$nt_EBB = 10000
      simdat$npa_EBB = 10000
      simdat$nps_EBB = 10000
      depth_add = depth2[,1:(h)]
      depth_add$avg = apply(depth_add[,1:h], 1, mean)
      depth_add_order = depth_add[order(depth_add$avg),]
      #add_names = row.names(depth_add_order[1:k,])
      add_names = row.names(depth_add_order[depth_add_order$avg<= depth_add_order$avg[k],])
      simdat$pre_EBB = as.numeric(all(add_names %in% names)) 
    }
  }
  simdat$time_EBB = (proc.time() - ptm)[['elapsed']]  
  
  check_list = rep(0, 10000)
  for (t in 2:10000){
    depth_med = depth2[,1:t]
    depth_med$avg = apply(depth_med, 1, mean)
    kth_value = sort(depth_med$avg, partial = k)[k]
    names_med = row.names(depth_med[which(depth_med$avg <= kth_value),])
    check_list[t] = as.numeric(all(names_med %in% names))
  }
  
  simdat$truth = tail(which(check_list ==0),1) + 1
  write.csv(as.matrix(t(simdat)), file = paste('synthetic/permutation/', nametag, '/simdat_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), row.names = FALSE)
  write.csv(time_list, file = paste('synthetic/permutation/', nametag, '/timelist_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), row.names = FALSE)
  write.csv(seed_list, file = paste('synthetic/permutation/', nametag, '/seedlist_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), row.names = FALSE)
}
}

args = commandArgs(trailingOnly = TRUE)
a = as.numeric(args[1])
evaluation(a)
  
  
  
  

