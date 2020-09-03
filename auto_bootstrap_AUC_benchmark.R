
## this version correct the update frequency issue in boostrap 

## landsat
## LR
## OCR
## pageb
## shuttle
## covertype

## construct folders




## (1) this version uses 1/100 Bonferroni correction
## (2) modified the output when the algorithm doesn't stop until 10000, 
##     from using "add_names = row.names(depth_add_order[1:k,])"
##     to         "add_names = row.names(depth_add_order[depth_add_order$avg<= depth_add_order$avg[k],])"

## automate the premuation method on benchmark data sets


## things to change
## (1) all input and output routes
## (2) the data frame holding outputs
## (3) Connect data resource and algorithm 


## -----------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------

### set parameters
evaluation <- function(a){ 
  library(pROC)
  B = 10000
  delta = 0.05
  # modify the following line to reduce computation ----------------------> done
  #sub_delta = delta/200
  sub_delta = delta/100

  n = 1000
  # hard code k to 5 instead
  #k = round(n*proportion)
  vpro = c(0.005)
  j = 1
  k = 5
  
  
  vector_name = c("landsat", "LR", "OCR", "pageb", "shuttle", "covertype")

  
  for (v in 1:length(vector_name)){
 #   proportion = vpro[j]/100
    data_name = vector_name[v]
    rseed <- as.numeric(unlist(read.csv("general_random_seed.csv", header = FALSE)))
    set.seed(rseed[a])
    
    ##------------------------------------------------------
    if (data_name == "landsat"){
      data_trn=read.csv('sat.trn',sep="", header = FALSE)
      data_tst=read.csv('sat.tst',sep="", header = FALSE)
      names(data_tst)=names(data_trn)
      data=rbind(data_trn,data_tst)
      n_col=ncol(data)
      o_col=n_col
      known_classes=c(1,7)
      big_nominal = subset(data,data[[o_col]] %in% known_classes)
      nominal = big_nominal[sample(nrow(big_nominal), round(n * (1 - vpro[j]))), ]
      big_anomaly = subset(data, !data[[o_col]] %in% known_classes)
      anomaly = big_anomaly[sample(nrow(big_anomaly), round(n * vpro[j])), ]
      data2 = rbind(anomaly, nominal)
      data2 = data2[-o_col]
    }else if(data_name == "LR"){
      data=read.csv('letter-recognition.data', header = FALSE)
      n_col=ncol(data)
      o_col=1
      data[[1]]=as.numeric(data[[1]])
      known_classes=c(1,3,5,7,9,11,13,15,17,19,21,23,26,25)
      big_nominal = subset(data,data[[o_col]] %in% known_classes)
      nominal = big_nominal[sample(nrow(big_nominal), round(n * (1 - vpro[j]))), ]
      big_anomaly = subset(data, !data[[o_col]] %in% known_classes)
      anomaly = big_anomaly[sample(nrow(big_anomaly), round(n * vpro[j])), ]
      data2 = rbind(anomaly, nominal)
      data2 = data2[-o_col]
    }else if(data_name == "OCR"){
      data_trn=read.csv('optdigits.tra', header = FALSE)
      data_tes=read.csv('optdigits.tes', header = FALSE)
      names(data_tes)=names(data_trn)
      data=rbind(data_trn,data_tes)
      n_col=ncol(data)
      o_col=n_col
      known_classes=c(1,3,4,5,7)
      big_nominal = subset(data,data[[o_col]] %in% known_classes)
      nominal = big_nominal[sample(nrow(big_nominal), round(n * (1 - vpro[j]))), ]
      big_anomaly = subset(data, !data[[o_col]] %in% known_classes)
      anomaly = big_anomaly[sample(nrow(big_anomaly), round(n * vpro[j])), ]
      data2 = rbind(anomaly, nominal)
      data2 = data2[-o_col]
    }else if(data_name == "pageb"){
      data=read.csv("page-blocks.data",sep="", header = FALSE)
      n_col=ncol(data)
      o_col=n_col
      known_classes=c(1,5)
      big_nominal = subset(data,data[[o_col]] %in% known_classes)
      nominal = big_nominal[sample(nrow(big_nominal), round(n * (1 - vpro[j]))), ]
      big_anomaly = subset(data, !data[[o_col]] %in% known_classes)
      anomaly = big_anomaly[sample(nrow(big_anomaly), round(n * vpro[j])), ]
      data2 = rbind(anomaly, nominal)
      data2 = data2[-o_col]
    }else if(data_name == "shuttle"){
      data_trn=read.csv('shuttle.trn',sep="", header = FALSE)
      data_tst=read.csv('shuttle.tst',sep="", header = FALSE)
      names(data_tst)=names(data_trn)
      data=rbind(data_trn,data_tst)
      o_col=ncol(data)#output class column
      ##known_class=c(1)
      known_classes=c(1,4)
      big_nominal = subset(data,data[[o_col]] %in% known_classes)
      nominal = big_nominal[sample(nrow(big_nominal), round(n * (1 - vpro[j]))), ]
      big_anomaly = subset(data, !data[[o_col]] %in% known_classes)
      anomaly = big_anomaly[sample(nrow(big_anomaly), round(n * vpro[j])), ]
      data2 = rbind(anomaly, nominal)
      data2 = data2[-o_col]
    }else if(data_name == "covertype"){
      ## only keep the numerical features
      data=read.csv("covtype.data",sep=',', header= FALSE)
      o_col=ncol(data)#output class column
      known_classes=c(1,2, 3, 7)
      big_nominal = subset(data,data[[o_col]] %in% known_classes)
      nominal = big_nominal[sample(nrow(big_nominal), round(n * (1 - vpro[j]))), ]
      big_anomaly = subset(data, !data[[o_col]] %in% known_classes)
      anomaly = big_anomaly[sample(nrow(big_anomaly), round(n * vpro[j])), ]
      data2 = rbind(anomaly, nominal)
      data2 = data2[, 1:10]
    }
  

    nametag <- paste("auto_bootstrap_AUC_correct_", data_name, "_pro",vpro[j],"_n",n, sep = "")
    write.csv(data2, file = paste("benchmark/", nametag, "/data2_pro",vpro[j],"_n",n, "_", a, ".csv", sep = ""), row.names = FALSE)
   
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
    system(paste('./iforest','-i', paste("benchmark/", nametag, "/data2_pro",vpro[j],"_n",n, "_", a, ".csv", sep = "") ,'-o',paste('benchmark/', nametag, '/depth2_pro',vpro[j],'_n',n, '_', a, '.csv', sep = ""),'-t 200 -s 0 -p -x', paste('benchmark/', nametag, '/data2_pro',vpro[j],'_n',n, '_', a, '.csv', sep = ""),'-z', curr_seed, '-g'), wait = TRUE)
    depth <- read.csv(paste('benchmark/', nametag, '/depth2_pro',vpro[j],'_n',n,'_', a, '.csv', sep = ""), header = FALSE)
    depth <- depth[, 2:201]
    system(paste('rm','-f',paste('benchmark/', nametag, '/depth2_pro',vpro[j],'_n',n, '_', a, '.csv', sep = "")))
    
      
    for (i in 2:50){
      curr_seed = seed_to_use[i]
      seed_list[i] = curr_seed    
      time_list[i] = Sys.time()
      system(paste('./iforest','-i', paste('benchmark/', nametag, '/data2_pro',vpro[j],'_n',n,'_', a, '.csv', sep = "") ,'-o',paste('benchmark/', nametag, '/depth2_pro',vpro[j],'_n',n,'_', a, '.csv', sep = ""),'-t 200 -s 0 -p -x', paste('benchmark/', nametag, '/data2_pro',vpro[j],'_n',n, '_', a, '.csv', sep = ""), '-z', curr_seed, '-g'), wait = TRUE)
      depth_mid <- read.csv(paste('benchmark/', nametag, '/depth2_pro',vpro[j],'_n',n, '_', a, '.csv', sep = ""), header = FALSE)
      depth_mid <- depth_mid[, 2:201]
      depth = cbind(depth, depth_mid)
      system(paste('rm','-f',paste('benchmark/', nametag, '/depth2_pro',vpro[j],'_n',n, '_', a, '.csv', sep = "")))
    }
    
    avg_d <- apply(depth, 1, mean)
    score_d <- ans(avg_d, c_n(n))
    target_names = which(score_d >= score_d[(order(-score_d))[k]])
    simdat$targetsize = length(target_names)
    
    label = c(rep(1, round(n*vpro[j])), rep(0, round(n * (1 - vpro[j]))))
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
    
    simdat$success_bt = mean(target_names %in% A_hat) == 1
    simdat$recall_bt = mean(target_names %in% A_hat)
    simdat$size_bt = length(A_hat)
    simdat$threshold = threshold
    #simdat$nt_bt = 50 * w
    simdat$nt_bt = 100 * w
    
    check_list = rep(0, 10000)

    for (t in 2:10000){
      depth_med = depth[,1:t]
      depth_med$avg = apply(depth_med, 1, mean)
      two_kth_value = sort(depth_med$avg, partial = (2*k))[2*k]
      names_med = row.names(depth_med[which(depth_med$avg <= two_kth_value),])
      check_list[t] = as.numeric(all(target_names %in% names_med))
    }

    simdat$truth = tail(which(check_list ==0),1) + 1
    write.csv(as.matrix(t(simdat)), file = paste('benchmark/', nametag, '/simdat_pro',vpro[j],'_n',n, '_', a, '.csv', sep = ""), row.names = FALSE)
    write.csv(time_list, file = paste('benchmark/', nametag, '/timelist_pro',vpro[j],'_n',n, '_', a, '.csv', sep = ""), row.names = FALSE)
    write.csv(seed_list, file = paste('benchmark/', nametag, '/seedlist_pro',vpro[j],'_n',n, '_', a, '.csv', sep = ""), row.names = FALSE)    
    
  }
}

args = commandArgs(trailingOnly = TRUE)
a = as.numeric(args[1])
evaluation(a)








