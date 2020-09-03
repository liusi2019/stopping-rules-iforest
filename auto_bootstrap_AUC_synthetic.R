

## auto bootstrap based method for synthetic data sets


## things to change
# (1) parameter setting and data generating -- done
# (2) growing trees -- done
# (3) add the ground truth part -- done
# (4) modify the output -- done

# concerns: 
# algorithm using 2 * threshold
# ground truth using 2 * k?


## file for general case with parameter setting for target two

## (1) this version uses 1/100 Bonferroni correction
## (2) modified the output when the algorithm doesn't stop until 10000, 
##     from using "add_names = row.names(depth_add_order[1:k,])"
##     to         "add_names = row.names(depth_add_order[depth_add_order$avg<= depth_add_order$avg[k],])"



## add the 10 000 condition


## automate the bootstrap method on benchmark data sets

## things to change
## (1) all input and output routes
## (2) the data frame holding outputs
## (3) Connect data resource and algorithm 



### set parameters
evaluation <- function(a){ 
  
  B = 10000
  delta = 0.05
  # modify the following line to reduce computation ----------------------> done
  #sub_delta = delta/200
  sub_delta = delta/100
  
  library(MASS)
  library(pROC)

  ## -------------------------------------------------------------------------------------
  vsize <- c(100, 1000, 10000)
  vpro <- c(5, 10, 20, 30)
  #start from 10 instead because the results for 5 is already out there
  #vpro <- c(10, 20, 30)
  vdist <- c(3, 5, 8, 10)
  vd <- c(10, 20, 40, 100)
  veta <- c(0.10, 0.20)
  vtype <- c(0, 1)
  
  ####### ---------------------------------------------------------------------------
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
    
    proportion = vpro[j]/100
    #fix k to be 5
    #k = round(n*proportion)
    k = 5
    
    a_dim = round(d * eta)
  
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
  
  
    nametag <- paste("bt_", "pro", vpro[j],"_n",n,"_d", d, "_dist", dist, "_eta", eta, "_type", a_type, sep = "")
    write.csv(data2, file = paste("synthetic/bootstrap/", nametag, "/data2_pro",vpro[j],"_n",n,"_d", d, "_dist", dist, "_eta", eta, "_type", a_type, "_", a, ".csv", sep = ""), row.names = FALSE)
  
    ### get c(n)----this is the c(n) used in the isolation forest paper for approximating the average path length
    c_n <- function(n){
      H = log(n - 1) + 0.5772156649
      return(2*H - 2*(n-1)/n)
    }
    ### function for anomaly score----given average depth h, calculate the corresponding anormaly score
    ans <- function(h, c_n){
      return(2^{-(h/c_n)})
    }
  
    get_A_hat_b <- function(trees, n, k){
      boots_trees <- trees[, sample(1:ncol(trees), replace = TRUE)]
      avg_d <- apply(boots_trees, 1, mean)
      score_d <- ans(avg_d, c_n(n))
      A_hat_b <- which(score_d >= score_d[(order(-score_d))[k]])
      return(A_hat_b)
    }


    simdat<-data.frame(array(dim=c(1,10)))
    names(simdat)<-c('targetsize','nt_bt', 'stop_bt', 'size_bt', 'threshold', 'success_bt', 'recall_bt', 'time_bt', 'truth', "AUC")
    simdat$stop_bt = 0
  
    time_list = rep("b", 50)
    seed_list = rep(0, 50)
  
    op <- options(digits.secs = 30)
    total_seed <- as.numeric(unlist(read.csv("random_seed_nature.csv", header = FALSE)))
    seed_to_use <- total_seed[((a-1)* 50 +1) : (a * 50) ]
    curr_seed = seed_to_use[1]
    seed_list[1] = curr_seed
    time_list[1] = Sys.time()
    system(paste('./iforest','-i', paste('synthetic/bootstrap/', nametag, '/data2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = "") ,'-o',paste('synthetic/bootstrap/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""),'-t 200 -s 0 -p -x', paste('synthetic/bootstrap/', nametag, '/data2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""),'-z', curr_seed, '-g'), wait = TRUE)
    depth <- read.csv(paste('synthetic/bootstrap/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), header = FALSE)
    depth <- depth[, 2:201]
    system(paste('rm','-f',paste('synthetic/bootstrap/',  nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = "")))
  
    for (i in 2:50){
      curr_seed = seed_to_use[i]
      seed_list[i] = curr_seed    
      time_list[i] = Sys.time()
      system(paste('./iforest','-i', paste('synthetic/bootstrap/', nametag, '/data2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = "") ,'-o',paste('synthetic/bootstrap/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""),'-t 200 -s 0 -p -x', paste('synthetic/bootstrap/', nametag, '/data2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), '-z', curr_seed, '-g'), wait = TRUE)
      depth_mid <- read.csv(paste('synthetic/bootstrap/', nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), header = FALSE)
      depth_mid <- depth_mid[, 2:201]
      depth = cbind(depth, depth_mid)
      system(paste('rm','-f',paste('synthetic/bootstrap/',  nametag, '/depth2_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = "")))
    }
  
    avg_d <- apply(depth, 1, mean)
    score_d <- ans(avg_d, c_n(n))
    target_names = which(score_d >= score_d[(order(-score_d))[k]])
    simdat$targetsize = length(target_names)
  
    label = c(rep(1, round(n* proportion)), rep(0, round(n*(1 - proportion))))  
    ROC1 = roc(label, ans(avg_d, c_n(n)), direction = "<")
    AUC1 = auc(ROC1)
    simdat$AUC = AUC1
    
    
    ptm = proc.time()
  # need to change it back to 1:100   ---------------------------------------------------------> done
    # for (w in 1:200){
    for (w in 1:100){
         # trees <- depth[, 1: (w * 50)]
         trees <- depth[, 1: (w * 100)]
         boots_target = lapply(1:B, function(x) get_A_hat_b(trees, n, k))
         boots_length = unlist(lapply(boots_target, function(x) length(x)))
         threshold <- mean(boots_length)
         is_it_in_union <- rep(1, n) ## point
         freq_in_replicates <- rep(0, n) ## point
         A_union_ind <- rep(1, B) ## this one is to keep track of which simulated target sets are we still considering now
         for (h in 1:n){
         freq_in_replicates[h] <- mean(unlist(lapply(boots_target, function(x) as.numeric(h %in% x))))
         }
         is_it_in_union[which(freq_in_replicates == 0)] = 0
         for (l in 1:n){
           all_freq =  freq_in_replicates * is_it_in_union
           ind_set_to_delete <- which(all_freq == min(all_freq[all_freq >0])) # indexes of all the points with smallest positive frequency
           ind_to_delete <- ind_set_to_delete[sample(length(ind_set_to_delete), 1)] # randomly pick one of those points to work on
           ## we need to use boots_target in the following line to keep track of the index of the set in the bootstrap resamples
           pool_set_to_delete <- which(unlist(lapply(boots_target, function(x) as.numeric(ind_to_delete %in% x))) > 0) # which sets does this point belong to
           pool_set_to_delete <- intersect(pool_set_to_delete, which(A_union_ind != 0))
           set_to_delete <- pool_set_to_delete[sample(length(pool_set_to_delete), 1)] # randomly choose one of the sets to delete
           A_union_ind[set_to_delete] = 0
           new_boots_target = boots_target[A_union_ind != 0]
           new_freq_in_replicates = rep(1, n)
           for (h in 1:n){
             new_freq_in_replicates[h] <- mean(unlist(lapply(new_boots_target, function(x) as.numeric(h %in% x))))
           }
           is_it_in_union[which(new_freq_in_replicates == 0)] = 0
           if (sum(A_union_ind) == ceiling((1 - sub_delta) * B)){ # should use B instead of explicit 10000
            break 
           }
         }
         A_hat = unique(unlist(new_boots_target))
         if (length(A_hat) <= 2 * threshold){
           simdat$stop_bt = 1
           break
         }
         freq_in_replicates = new_freq_in_replicates 
    }
    simdat$time_bt = (proc.time() - ptm)[['elapsed']]  

    if (simdat$stop_bt == 0){
      A_hat = which(score_d >= score_d[(order(-score_d))[2*k]])
    }
    simdat$success_bt = mean(target_names %in% A_hat) == 1
    simdat$recall_bt = mean(target_names %in% A_hat)
    simdat$size_bt = length(A_hat)
    simdat$threshold = threshold
    #simdat$nt_bt = 50 * w
    simdat$nt_bt = 100 * w
    
    check_list = rep(0, 10000)
    ## need to change it back to 10000 ----------------------------------->
    for (t in 2:10000){
      depth_med = depth[,1:t]
      depth_med$avg = apply(depth_med, 1, mean)
      two_kth_value = sort(depth_med$avg, partial = (2*k))[2*k]
      names_med = row.names(depth_med[which(depth_med$avg <= two_kth_value),])
      check_list[t] = as.numeric(all(target_names %in% names_med))
    }
    simdat$truth = tail(which(check_list ==0),1) + 1
    
    write.csv(as.matrix(t(simdat)), file = paste('synthetic/bootstrap/', nametag, '/simdat_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), row.names = FALSE)
    write.csv(time_list, file = paste('synthetic/bootstrap/', nametag, '/timelist_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), row.names = FALSE)
    write.csv(seed_list, file = paste('synthetic/bootstrap/', nametag, '/seedlist_pro',vpro[j],'_n',n,'_d', d, '_dist', dist, '_eta', eta, '_type', a_type, '_', a, '.csv', sep = ""), row.names = FALSE)
  }
}

args = commandArgs(trailingOnly = TRUE)
a = as.numeric(args[1])
evaluation(a)