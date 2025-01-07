library(gtools)
library(mvtnorm)
library(ggplot2)
library(patchwork)

#########################################

# For displaying results

heat_tile <- function(edge_probs, title = 'Edge inclusion probabilities'){
  
  # rownames(edge_probs) = colnames(edge_probs) = names(leukemia)
  
  edge_probs_melt <- reshape2::melt(edge_probs)
  
  heatmap <- ggplot(data=edge_probs_melt, aes(factor(Var2), factor(Var1), fill= value)) +
    geom_tile(color = "white",
              lwd = 1,
              linetype = 1) +
    geom_text(aes(label = round(value,2)), color = "white", size = 4.5)+
    scale_fill_viridis_c(option='G', direction = -1) +
    scale_x_discrete(name = '', expand = c(0, 0))+
    scale_y_discrete(name = '', expand = c(0, 0))+
    guides(fill = guide_colourbar(barwidth = 1,
                                  barheight = 10,
                                  title = "probs"))+
    ggtitle(title)+
    theme(plot.title = element_text(size = 17, face = "bold",hjust = 0.5),
          axis.title.y=element_text(angle=0, vjust=0.6),
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 12))
  
  heatmap
  
}

symmetric_matrix <- function(matrix) {
  n <- nrow(matrix)
  sym_matrix <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (matrix[i, j] == 1 || matrix[j, i] == 1) {
        sym_matrix[i, j] <- 1
        sym_matrix[j, i] <- 1
      }
    }
  }
  
  return(sym_matrix)
}

#########################################

# From BCDAG package

rDAG = function(q, w){
  DAG = matrix(0, q, q); colnames(DAG) = rownames(DAG) = 1:q
  DAG[lower.tri(DAG)] = stats::rbinom(n = q*(q-1)/2, size = 1, prob = w)
  return(DAG)
}

rDAGWishart <- function(n, DAG, a, U) {
  
  q <- ncol(DAG)
  ajs <- sapply(1:q, function(j) a+sum(DAG[,j]==1)-q+1)
  # this is the a star defined from the theory which 
  # assures compatibility in terms of Markov properties 
  
  L.array <- array(0, dim = c(q,q,n))
  D.array <- array(0, dim = c(q,q,n))
  
  for (i in 1:n) {
    params <- lapply(1:q, function(j) rnodeDAGWishart(j, DAG, ajs[j], U))
    
    sigmas <- sapply(1:q, function(x) params[[x]]$sigmaj)
    L <- lapply(1:q, function(x) params[[x]]$Lj)
    # I simply take out sigmas and L from params which already contains them 
    
    D.array[,,i] <- diag(sigmas) # I fill one matrix
    for (j in 1:q) {
      whc <- which(DAG[,j] == 1) # I identify the parents 
      L.array[whc,j,i] <- as.numeric(L[[j]])
      
      # The Cholesky factorization (L) represents the lower triangular matrix 
      # used to decompose the covariance matrix. For each parent node identified 
      # (whc), the Cholesky factorization values from the L list (which contains 
      # Cholesky factorizations for each node in the DAG) are assigned to the 
      # corresponding positions in the L.array for the current sample i. This 
      # ensures that the conditional distribution of each node is constructed 
      # accurately based on its parents' contributions.
    }
    diag(L.array[,,i]) <- 1
  }
  
  if (n == 1) {
    D.array <- D.array[,,1]
    L.array <- L.array[,,1]
  }
  
  return(list(D = D.array, L = L.array))
}

rnodeDAGWishart <- function(node, DAG, aj, U) {
  
  # the function aims to sample parameters related to the covariance matrix 
  # associated with the conditional distribution of the specified node given 
  # its parents in the DAG. 
  
  q <- ncol(data)
  n <- nrow(data)
  
  j <- node
  pa <- pa(j, DAG)
  
  out <- list(sigmaj = 0, Lj = 0)
  
  if (length(pa) == 0) {
    U_jj <- U[j,j]
    out$sigmaj <- stats::rgamma(1, shape = aj/2, rate = U_jj/2)^-1
  } else {
    U_paj.j <- U[pa,j]
    invU_papa <- chol2inv(chol(U[pa,pa]))
    U_jj <- U[j,j] - t(U_paj.j)%*%invU_papa%*%U_paj.j
    
    out$sigmaj <- stats::rgamma(1, shape = aj/2, rate = U_jj/2)^-1
    out$Lj <- mvtnorm::rmvnorm(1, -invU_papa%*%U_paj.j, out$sigmaj*invU_papa)
  }
  
  return(out)
}

pa <- function(node, DAG) {
  pa <- which(DAG[,node] != 0)
  return(pa)
}

operation <- function(op, A, nodes) {
  x <- nodes[1]
  y <- nodes[2]
  
  if(op == 1) {
    A[x,y] = 1
    return(A)
  }
  
  if(op == 2) {
    A[x,y] = 0
    return(A)
  }
  
  if(op == 3) {
    A[x,y] = 0
    A[y,x] = 1
    return(A)
  }
}

DW_nodelml <- function(node, DAG, tXX, n, a, U) {
  # slide 19
  j <- node
  pa <- pa(j, DAG)
  # for that column being different from zero, which(DAG[,node] != 0)
  q <- ncol(tXX) # number  of parameters
  
  a.star <- (a+length(pa)-q+1)
  # pag. 6 it guarantees compatibility among prior distributions 
  # for Markov equivalent DAGs
  
  
  Upost <- U + tXX
  # updating
  
  if (length(pa) == 0) {
    # isolated node:
    U_jj <- U[j,j]
    Upost_jj <- Upost[j,j]
    
    prior.normcost <- -lgamma(a.star/2) + a.star/2*log(U_jj/2)
    # l stands for logarithm 
    post.normcost <- -lgamma(a.star/2 + n/2) + (a.star/2 + n/2)*log(Upost_jj/2)
    # in logarithmic scale
    # p.8 I have all the formula 
    
  } else {
    U_paj.j <- U[pa,j]
    U_jj <- U[j,j] - t(U_paj.j)%*%chol2inv(chol(U[pa,pa]))%*%U_paj.j
    # chol2inv -> Inverse from Choleski 
    Upost_paj.j <- Upost[pa,j]
    Upost_jj <- Upost[j,j] - t(Upost_paj.j)%*%chol2inv(chol(Upost[pa,pa]))%*%Upost_paj.j
    
    prior.normcost <- -lgamma(a.star/2) + a.star/2*log(U_jj/2) + 0.5*log(det(as.matrix(U[pa,pa])))
    post.normcost <- -lgamma(a.star/2 + n/2) + (a.star/2 + n/2)*log(Upost_jj/2) + 0.5*log(det(as.matrix(Upost[pa,pa])))
    
  }
  
  nodelml <- -n/2*log(2*pi) + prior.normcost - post.normcost
  
  return(nodelml)
}

#########################################

# BDMCMC implementation

softmax <- function(x, scale_factor = 1) {
  x_scaled <- x / scale_factor
  max_x_scaled <- max(x_scaled)
  exp_x_scaled <- exp(x_scaled - max_x_scaled)
  return(exp_x_scaled / sum(exp_x_scaled))
}

time_dag_update <- function(A, tXX, n, a, U, w, return_DAG = F){
  
  q <- ncol(A)
  A_na <- A
  diag(A_na) <- NA
  
  id_set = c()
  dd_set = c()
  
  ## set of nodes for id
  set_id = which(A_na == 0, TRUE)
  # If I want to add I need nothing at that location
  if(length(set_id) != 0){
    id_set = cbind(1, set_id)
  }
  len_id <- nrow(id_set)
  
  ## set of nodes for dd
  set_dd = which(A_na == 1, TRUE)
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }
  
  O <- rbind(id_set, dd_set)
  nrows <- nrow(O)
  # need to consider somehow if a valid DAG or not!!
  
  # start_time <- Sys.time()
  
  log.prior.ratio_ins <-log(w / (1 - w))
  log.prior.ratio_del <- log((1 - w) / w)
  # to be given as input so that I don't need to recompute each time 
  
  probs <- sapply(1:nrows, function(i) {
    
    A_new <- operation(O[i, 1], A, O[i, 2:3])
    
    if (gRbase::is.DAG(A_new)) {
      log_prior_ratio <- ifelse(i <= len_id, log.prior.ratio_ins, log.prior.ratio_del)
      
      log_new_ml <- DW_nodelml(O[i, 3], A_new, tXX, n, a, U)
      log_old_ml <- DW_nodelml(O[i, 3], A, tXX, n, a, U)
      # cat(log_new_ml, end = ' ')
      # cat(log_old_ml, end = ' ')
      # cat(log_new_ml - log_old_ml + log_prior_ratio, end = ' ')
      return(log_new_ml - log_old_ml + log_prior_ratio)
    } else {
      return(NA)
    }
  })
  valid_indices <- which(!is.na(probs))
  probs_valid <- probs[valid_indices]
  # plot(probs_valid/sum(probs_valid))
  # valid_probs <- exp(probs_valid/sum(probs_valid))
  # print(valid_probs)
  # plot(exp(valid_probs))
  range_adjusted_log_weights <- range(probs_valid - max(probs_valid))
  scale_factor <- abs(min(range_adjusted_log_weights)) / 10
  valid_probs <- softmax(probs_valid, scale_factor)
  # plot(valid_probs)
  # print(valid_probs)
  valid_probs[is.infinite(valid_probs)] <- 1e16
  # par(mfrow=c(1,3))
  # plot(probs_valid)
  # plot(valid_probs)
  # plot(probs_valid/sum(probs_valid))
  index <- sample(valid_indices, 1, prob = valid_probs)
  # print(paste(length(valid_probs),which.max(valid_probs),O[index,1]))
  
  A_next <- operation(O[index, 1], A, O[index, 2:3])
  # DAG <- A_next # !!!!!!!!!!!
  # par(mfrow=c(1,2))
  # plot(probs[valid_indices])
  # plot(exp(probs[valid_indices]), log = 'y')
  # print(exp(probs[valid_indices]))
  
  # print(paste('sum is',sum(probs[valid_indices])))
  # boxplot(valid_probs)
  exp_probs_valid <- exp(probs_valid)
  exp_probs_valid[is.infinite(exp_probs_valid)] <- 1e16
  exp_probs_valid[exp_probs_valid == 0] <- 1e-16
  exp_probs_valid[exp_probs_valid < 1e-7] <- 1e-7
  wait_time <- 1/(sum(exp_probs_valid))
  if (wait_time > 400){
    print(paste('-- HERE IS BIGGER THAN 400',wait_time,sum(probs_valid),sum(exp_probs_valid)))
    print(wait_time)
    print(exp_probs_valid, end= ' ')
    print(probs_valid, end= ' ')
  }
  if (wait_time > 1e-07 & wait_time < 1e-06){
    print(paste('-- HERE IS lower THAN 1',wait_time,sum(probs_valid),sum(exp_probs_valid)))
    # print(paste(wait_time),sum(probs_valid))
    print(exp_probs_valid, end= ' ')
    print(probs_valid, end= ' ')
    # print(sum(probs_valid))
  }
  # print(wait_time)
  # print(paste('wait_time', 1/sum(exp(probs[valid_indices]))))
  if (return_DAG == T){
    return(list(A_next = A_next, wait_time = wait_time))
  } else {
    return(wait_time)
  }
}

local_moves <- function(DAG, tXX, n, a, U, w){
  
  k = 0
  waiting_time=c()
  dags <- list()
  
  for (i in 1:200){
    dags[[i]] <- DAG
    # if (i %in% 1:15){
    #   print(heat_tile(true, title='True DAG')+heat_tile(A))
    #   Sys.sleep(1)
    # }
    
    next_move <- time_dag_update(A = DAG, tXX, n, a, U, w, return_DAG = T)
    A_next <- next_move$A_next 
    wait_time <- next_move$wait_time
    DAG <- A_next
    # heat_tile(true)+heat_tile(symmetric_matrix(true))+heat_tile(DAG)
    
    # print(paste(exp(wait_time),' - ', sum(!(true == DAG))))
    # print(exp(wait_time))
    
    # before <- DAG
    
    waiting_time[i] <- wait_time 
    if (i >= 10) {
      if (i %% 10 == 0){
        if (waiting_time[i] <= waiting_time[i - 9]){
          break
        }
      }
    }
    # if (i == 1){
    #   times[[j]][1]=wait_time
    # } else {
    #   times[[j]] <- c(times[[j]],wait_time)
    # }
    if (i > 4){
      if (waiting_time[i-2]==wait_time) {
        k = k + 1
        # print(paste('max_time changed in',max_time))
      } else  {
        k = 0
      }
    }
    if (k > 5){
      print('k>5')
      break
    }
    
  }
  
  # final_times <- tail(waiting_time,length(waiting_time)-4)
  # final_dag <- dags[[length(waiting_time)-length(final_times)+which.max(tail(waiting_time,length(final_times)))]]
  # final_dags <- list()
  # plot(waiting_time)
  # for (i in 1:(length(dags)-4)){
  #   final_dags[[i]] <- dags[[i+4]]
  # }
  final_dag <- dags[[which.max(waiting_time)]]
  return(list(DAG = final_dag, length = length(waiting_time), max_time = max(waiting_time),
              waiting_times = waiting_time, dags = dags))
}

global_move <- function(A, tXX, n, a, U, w, q, n_op) {
  
  q <- ncol(A)
  A_na <- A
  diag(A_na) <- NA
  
  id_set = c()
  dd_set = c()
  
  ## set of nodes for id
  set_id = which(A_na == 0, TRUE)
  # If I want to add I need nothing at that location
  if(length(set_id) != 0){
    id_set = cbind(1, set_id)
  }
  len_id <- nrow(id_set)
  
  ## set of nodes for dd
  set_dd = which(A_na == 1, TRUE)
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }
  
  O <- rbind(id_set, dd_set)
  
  comb_indices <- combinations(nrow(O), n_op)
  
  combs <- matrix(NA, nrow(comb_indices), n_op*3)
  
  for (i in 1:nrow(comb_indices)){
    for (j in 1:n_op){
      combs[i,j] <- O[comb_indices[i,j],1]
    }
    k <-  n_op + 1
    for (j in rep(1:n_op)){ 
      combs[i,k] <- O[comb_indices[i,j],2]
      k = k+1
      combs[i,k] <- O[comb_indices[i,j],3]
      k = k+1
    }
  }
  
  log.prior.ratio_ins <-log(w / (1 - w))
  log.prior.ratio_del <- log((1 - w) / w)
  
  log_old_mls <- numeric(q)
  for (i in 1:q){
    log_old_mls[i] <- DW_nodelml(i, A, tXX, n, a, U)
  }
  
  # filtered_indices <- !apply(combs[, 1:n_op], 1, function(x) sum(x == 2) >= 2)
  # n_samples <- min(10000, nrow(filtered_combs))
  # sampled_rows <- sample(nrow(filtered_combs), n_samples, replace = TRUE)
  # rows_with_2 <- apply(combs[, 1:n_op], 1, function(x) sum(x == 2) >= 2)
  # indices <- c(sampled_rows, rows_with_2)
  # indices <- sort(indices)
  
  filtered_combs <- combs[!apply(combs[, 1:3], 1, function(x) sum(x == 2) >= 2), ]
  n_samples <- min(10000, nrow(filtered_combs))
  sampled_rows <- filtered_combs[sample(nrow(filtered_combs), n_samples, replace = TRUE), ]
  rows_with_2 <- combs[apply(combs[, 1:3], 1, function(x) sum(x == 2) >= 2), ]
  combined <- rbind(sampled_rows, rows_with_2)
  n_rows <- nrow(combined)
  
  # indices <- sample(1:nrow(combs),n, replace = F)
  # indices <- sort(indices)
  
  sampled_dags <- list()
  probs <- numeric(n_rows)
  for (i in 1:n_rows){
    # probs <- sapply(1:n, function(i) {
    
    A_new <- A
    k <- n_op +1
    for (j in 1:(n_op)){
      A_new <- operation(combined[i,j], A_new, combined[i, k:(k+1)])
      if (k < (n_op)*3){
        k <- k + 2
      } else {
        break
      }
    }
    sampled_dags[[i]] <- A_new
    
    if (gRbase::is.DAG(A_new)) {
      
      log_prior_ratios <- sapply(1:n_op, function(j) ifelse(combined[i,j]==1, log.prior.ratio_ins, log.prior.ratio_del))
      log_prior_ratio <- sum(log_prior_ratios)
      
      log_new_mls <- sapply(seq(n_op+2,(n_op+1)+(n_op)*2,2), function(k) {
        DW_nodelml(combined[i, k], A_new, tXX, n, a, U)})
      log_new_ml <- sum(log_new_mls)
      
      log_old_ml <- sapply(seq(n_op+2,(n_op+1)+(n_op)*2,2), function(j) {
        log_old_mls[combined[i,j]]})
      log_old_ml_tot <- sum(log_old_ml)
      
      probs[i] <- (log_new_ml - log_old_ml_tot + log_prior_ratio)
    } else {
      probs[i] <- NA
    }
  }
  
  valid_indices <- which(!is.na(probs))
  probs_valid <- probs[valid_indices]
  range_adjusted_log_weights <- range(probs_valid - max(probs_valid))
  scale_factor <- abs(min(range_adjusted_log_weights)) / 10
  valid_probs <- softmax(probs_valid, scale_factor)
  # plot(valid_probs)
  # index <- sample(valid_indices, 1, prob = valid_probs)
  # print(paste(length(valid_probs),which.max(valid_probs),O[index,1]))
  d <- 10
  sorted_valid_indices <- order(valid_probs, decreasing = TRUE)[1:d]
  sorted_original_indices <- valid_indices[sorted_valid_indices]
  DAG_for_local <- lapply(sorted_original_indices, function(x) sampled_dags[[x]])
  
  times_after_global <- sapply(1:d, function(i) {
    time_dag_update(DAG_for_local[[i]], tXX, n, a, U, w, return_DAG = F)})
  max_time_global <- max(times_after_global)
  which_max_time_global <- which.max(times_after_global)
  max_DAG_global <- DAG_for_local[[which_max_time_global]]
  
  max_times <- numeric(d)
  length_times <- numeric(d)
  final_dags <- list()
  dags_sequence <- list()
  times_sequence <- list()
  for (i in 1:d){
    #' need to get the DAG to input the local function 
    #' compute local update, get the time 
    #' decide on the maximum local 
    #' return that update to the outer two functions 
    #' see how it is going, from there understand the situation 
    #' implement parallelization and see if I have acceptable results 
    #' early stopping very important, find a criterion, start with 1e-5 e.g.
    
    dag <- DAG_for_local[[i]]
    # this i is like valid_indices, thus the subset
    after_local <- local_moves(DAG = dag, tXX, n, a, U, w)
    final_dags[[i]] <- after_local$DAG
    max_times[i] <- after_local$max_time
    length_times[i] <- after_local$length
    dags_sequence[[i]] <- after_local$dags
    times_sequence[[i]] <- after_local$waiting_times
    
  }
  max_time_local <- max(max_times)
  print(paste('MAX TIME LOCAL', max_time_local, 'MAX TIME GLOBAL', max_time_global))
  selected_move <- which.max(max_times)
  selected_DAG <- final_dags[[selected_move]]
  selected_length <- length_times[selected_move]
  selected_dags <- dags_sequence[[selected_move]]
  selected_waiting_times <- times_sequence[[selected_move]]
  
  if (max_time_local > max_time_global){
    print('-Done local!!')
    print(selected_waiting_times, end = ' ')
    # print('INSIDE GLOBAL CHECK COMPARISON')
    # for (i in 1:length(selected_waiting_times)){
    #   up <- time_dag_update(A = selected_dags[[i]], tXX, n, a, U, w, return_DAG = F)
    #   cat(up, end = ' ')
    #   # print(sum(selected_dags[[i]]==1), end = ' ')
    # }
    return(list(DAG = selected_DAG, length = selected_length, 
                dags_sequence = selected_dags, times_sequence = selected_waiting_times))
  } else {
    print(paste('Inside when local is not done',length(max_time_global)))
    print(max_time_global)
    return(list(DAG = max_DAG_global, length = 1, 
                times_sequence = max_time_global, dags_sequence = max_DAG_global))
  }
  
}

local_and_global_steps <- function(A, tXX, n, a, U, w, number_of_moves, previous_wait_time){
  
  #' local moves 
  #' sample how many operations for the global
  #' if likelihood and prior is better then I go with other locals 
  #' better means above a certain threshold (to be assessed)
  #' if not better I try to change the n. of operations, if nothing works I do stop 
  #' If in general I've done more then a certain number of steps I stop 
  
  options <- c(2:3)
  n_op_vec <- c()
  # print(heat_tile(true, title='True DAG')+heat_tile(A))
  # Sys.sleep(1)
  n_op <- 2
  print(paste('2- n_op sampled is',n_op))
  n_op_vec <- c(n_op_vec, n_op)
  # cat(paste('Vector is', paste(n_op_vec, collapse = ',')), "\n")
  after_global <- global_move(A = A, tXX, n, a, U, w, q, n_op)
  # print('INSIDE STEP CHECK COMPARISON')
  # print(after_global$times_sequence, end = '')
  # for (i in 1:length(after_global$times_sequence)){
  #   up <- time_dag_update(A = after_global$dags_sequence[[i]], tXX, n, a, U, w, return_DAG = F)
  #   cat(up, end = '  ')
  #   # print(sum(after_global$dags_sequence[[i]]==1), end = ' ')
  # }
  number_of_moves <- number_of_moves +  after_global$length
  # print(DAG_after_global)
  wait_time_global <- max(after_global$times_sequence)
  print(paste('IN STEP WAIT TIME GLOBAL', wait_time_global))
  print(paste('after first global n.wrong',sum(!(true == after_global$DAG)),'wait time is', wait_time_global))
  
  if (wait_time_global > previous_wait_time){
    print(paste('comparing current with previous',wait_time_global,previous_wait_time))
    print(heat_tile(true, title='True DAG')+heat_tile(symmetric_matrix(true), title = 'True with symmetric edges')+heat_tile(after_global$DAG, title = 'selected DAG'))
    Sys.sleep(1)
    return(list(better = TRUE, number_of_moves = number_of_moves, DAG = after_global$DAG, 
                time = wait_time_global, dags_sequence = after_global$dags_sequence, 
                times_sequence = after_global$times_sequence))
  } else {
    while (TRUE){
      if (length(n_op_vec) == 2){
        print('Not found better DAG configuration with any global move')
        return(list(better = FALSE, time = NULL, DAG = A, number_of_moves = 0))
      }
      print('SAMPLING OTHER OP')
      # remaining_options <- setdiff(options, n_op_vec)
      # if (length(remaining_options) == 1) {
      #   n_op <- remaining_options
      # } else {
      #   n_op <- sample(remaining_options, 1)
      # }
      n_op = 3
      print(paste('n_op sampled is',n_op))
      n_op_vec <- c(n_op_vec, n_op)
      # cat(paste('Vector is', paste(n_op_vec, collapse = ',')), "\n")
      after_global <- global_move(A = A, tXX, n, a, U, w, q, n_op)
      # print(after_local$DAG == DAG_after_global)
      wait_time_global <- max(after_global$times_sequence)
      print(paste('after ANOTHER global n. wrong',sum(!(true == after_global$DAG)),'wait time is', wait_time_global))
      # print(sum(!(true == DAG_after_global)))
      if (wait_time_global > previous_wait_time){
        print(paste('comparing current with previous',wait_time_global,previous_wait_time))
        print(heat_tile(true, title='True DAG')+heat_tile(symmetric_matrix(true), title = 'True with symmetric edges')+heat_tile(after_global$DAG, title = 'selected DAG'))
        Sys.sleep(1)
        return(list(better = TRUE, number_of_moves = number_of_moves, DAG = after_global$DAG, 
                    time = wait_time_global, dags_sequence = after_global$dags_sequence, 
                    times_sequence = after_global$times_sequence))
      }
      else{
        print(paste('OFF: comparing current with previous',wait_time_global,previous_wait_time))
      }
    }
  }
  
}

birth_death <- function (tXX, n, a, U, w, max_moves = 1000){
  verbose = T
  
  Graphs <- array(NA, dim=c(q,q, max_moves))
  L <- array(NA, dim=c(q,q, max_moves))
  D <- array(NA, dim=c(q,q, max_moves))
  waiting_times <- numeric(max_moves)
  
  initial <- rDAG(q,w)
  after_local <- local_moves(initial, tXX, n, a, U, w)
  print(paste('1 - Did first local'))
  k = 1
  for (i in 1:after_local$length){
    Graphs[,,k] <- after_local$dags[[k]]
    waiting_times[k] <- after_local$waiting_times[k]
    print(waiting_times[k])
    # mmg <- Graphs[,,k]
    # up <- time_dag_update(A = mmg, tXX, n, a, U, w, return_DAG = T)
    # print(paste('COMPARING SAVED OBJECTS:',waiting_times[k], up$wait_time))
    k = k + 1
  }
  first_wait_time <- after_local$max_time
  # here need to fill with locals
  step <- list(number_of_moves = 0, DAG = after_local$DAG, time = first_wait_time)
  # here need to fill with global/locals
  print(paste('initial random dag', sum(!(initial == true)), 'wait time', first_wait_time))
  
  while (TRUE){
    
    # print(step$DAG)
    if (!is.null(step$time)) {
      previous_time <- step$time
      print(paste('PREVIOUS TIME IN STEP', previous_time))
    }
    step <- local_and_global_steps(A = step$DAG, tXX, n, a, U, w, 
                                   number_of_moves = step$number_of_moves, previous_wait_time = previous_time)
    
    if (step$better == FALSE){
      after_local <- local_moves(step$DAG, tXX, n, a, U, w)
      wait_time <- after_local$max_time
      print(paste('LAST LOCAL TIME', wait_time))
      if (wait_time > previous_time){
        print('LAST LOCAL ACCEPTED')
        for (i in 1:after_local$length){
          Graphs[,,k] <- after_local$dags[[i]]
          waiting_times[k] <- after_local$waiting_times[i]
          print(waiting_times[k])
          k = k + 1
          step$time <- after_local$max_time
          step$DAG <- after_local$DAG
          step$number_of_moves <- step$number_of_moves + after_local$length 
        } 
      }else {
        break
      }
    } else if (step$number_of_moves >= max_moves){
      print('exceeded number of  moves')
      break
    } else {
      print( '-- Global move accepted')
      print(paste(length(step$times_sequence), length(step$dags_sequence)))
      for (i in 1:length(step$times_sequence)){
        if (is.list(step$dags_sequence)) {
          Graphs[,,k] <- step$dags_sequence[[i]]
        } else {
          Graphs[,,k] <- step$dags_sequence
        }
        if (is.vector(step$times_sequence)) {
          waiting_times[k] <- step$times_sequence[i]
        } else {
          waiting_times[k] <- step$times_sequence
        }
        # print(waiting_times[k])
        # up <- time_dag_update(A = Graphs[,,k], tXX, n, a, U, w, return_DAG = T)
        # print(paste('COMPARING SAVED OBJECTS:',waiting_times[k], up$wait_time))
        k = k + 1
      }
      # here need to fill with global/locals
    }
  }
  # time <- time_dag_update(step$DAG, tXX, n, a, U, w, return_DAG = F)
  # return(list(DAG = step$DAG, time = time))
  if (verbose == TRUE) {
    cat("\nSampling parameters...")
    pb <- utils::txtProgressBar(min = 2, max = (k-1), style = 3)
  }
  L <- array(NA, dim = c(q,q,(k-1)))
  D <- array(NA, dim = c(q,q,(k-1)))
  for (i in 1:(k-1)) {
    postparams <- rDAGWishart(1, Graphs[,,i], a+n, U+tXX)
    L[,,i] <- postparams$L
    D[,,i] <- postparams$D
    if (verbose == TRUE) {
      utils::setTxtProgressBar(pb, i)
      close(pb)
    }
  }
  Graphs <- Graphs[,,1:(k-1)]
  waiting_times <- waiting_times[1:(k-1)]
  
  return(list(Graphs = Graphs, waiting_times = waiting_times, 
              L = L, D = D))
  
}

######################################

data_creation <- function(n = 1000, w = 0.2, random = F, boot = F){
  
  if (random == T) { 
    DAG <- rDAG(q,w)
    true = DAG
  } else {
    DAG <- matrix(0,8,8)
    nodes <- c(2,1,3,1,7,2,8,2,7,3,5,4,6,4)
    nod_mat <- matrix(nodes, length(nodes)/2,2, byrow = T)
    for (i in 1:nrow(nod_mat)){
      DAG[nod_mat[i,1],nod_mat[i,2]] = 1
    }
    true = DAG
  }
  q = nrow(DAG)
  w = w
  outDL = rDAGWishart(n = 1, DAG = DAG, a = q, U = diag(1, q))
  L = outDL$L; D = outDL$D
  Sigma = solve(t(L))%*%D%*%solve(L)
  X = mvtnorm::rmvnorm(n = n, sigma = Sigma)
  n <- nrow(X)
  if (boot == TRUE) {
    if (length(X)<800){
      remaining_samples <- 500-nrow(X)
      indices <- sample(1:nrow(X), remaining_samples, replace = TRUE)
      samples <- X[indices,]
      X <- rbind(X, samples)
    }
  }
  data <- X
  X <- scale(data, scale = FALSE)
  tXX <- crossprod(X)
  waiting_time <- c()
  
  assign("DAG", DAG, envir = .GlobalEnv)
  assign("q", q, envir = .GlobalEnv)
  assign("n", n, envir = .GlobalEnv)
  assign("w", w, envir = .GlobalEnv)
  assign("outDL", outDL, envir = .GlobalEnv)
  assign("L", L, envir = .GlobalEnv)
  assign("D", D, envir = .GlobalEnv)
  assign("Sigma", Sigma, envir = .GlobalEnv)
  assign("X", X, envir = .GlobalEnv)
  assign("data", data, envir = .GlobalEnv)
  assign("tXX", tXX, envir = .GlobalEnv)
  assign("waiting_time", waiting_time, envir = .GlobalEnv)
  assign("true", true, envir = .GlobalEnv)
}

diff_time <- function(start_time = start_time ){
  end_time <- Sys.time()
  
  # Calculate the elapsed time
  elapsed_time <- end_time - start_time
  
  return(elapsed_time)
}

######################################

q <- 8 ; w <- 0.2
data_creation(n = 1000, w = 0.2,  random = F, boot = F)
result = birth_death(tXX, n, a = q, U = diag(1,q)/n, w, max_moves = 400)
final_dags <- result$Graphs
wait <- result$waiting_times
mmg_index <- which.max(wait)
print(heat_tile(true, title='True DAG')+heat_tile(symmetric_matrix(true), title = 'True with symmetric edges')+heat_tile(final_dags[,,mmg_index], title = 'selected DAG'))
mmg <- result$Graphs[,,mmg_index]
causalDisco::shd(true, mmg)

##########################################

## PARALLELIZATION 

data_creation(n = 1000, w = 0.2,  random = F, boot = F)
heat_tile(true)

nn <- 20
final_dags <- list()
final_times <- list()
for (i in 1:nn){
  result = birth_death(tXX, n, a = q, U = diag(1,q)/n, w, max_moves = 400)
  # final_dags[[i]] <- result$DAG
  print(result$waiting_times)
  final_times[[i]] <- result$waiting_times
  print(i)
}

selected_index <- which.max(final_times)
print(heat_tile(true, title='True DAG')+heat_tile(symmetric_matrix(true), title = 'True with symmetric edges')+heat_tile(final_dags[[selected_index]], title = 'selected DAG'))

####################################

## SIMULATIONS 

dir1 <- "C:/D-drive-15734/Sara/Università/Magistrale/Secondo anno/Tesi/computational part/BDMCMC/BDMCMC runs/simulations/n = 50"
dir2 <- "C:/D-drive-15734/Sara/Università/Magistrale/Secondo anno/Tesi/computational part/BDMCMC/BDMCMC runs/simulations/n = 100"
dir3 <- "C:/D-drive-15734/Sara/Università/Magistrale/Secondo anno/Tesi/computational part/BDMCMC/BDMCMC runs/simulations/n = 200"
dir4 <- "C:/D-drive-15734/Sara/Università/Magistrale/Secondo anno/Tesi/computational part/BDMCMC/BDMCMC runs/simulations/n = 500"
dir5 <- "C:/D-drive-15734/Sara/Università/Magistrale/Secondo anno/Tesi/computational part/BDMCMC/BDMCMC runs/simulations/n = 1000"
dir5 <- "C:/D-drive-15734/Sara/Università/Magistrale/Secondo anno/Tesi/computational part/BDMCMC/BDMCMC runs/simulations/trials"
dirs <- c(dir1, dir2, dir3, dir4, dir5)

setwd("C:/D-drive-15734/Sara/Università/Magistrale/Secondo anno/Tesi/computational part/Adaptive Metropolis Hastings/Adaptive runs/simulations")
data_50 = readRDS('data_n_50.RDS')
data_100 = readRDS('data_n_100.RDS')
data_200 = readRDS('data_n_200.RDS')
data_500 = readRDS('data_n_500.RDS')
outDL = readRDS('original_parameters.RDS')
L = outDL$L; D = outDL$D
Sigma = solve(t(L))%*%D%*%solve(L)
set.seed(123)
X = mvtnorm::rmvnorm(n = 1000, sigma = Sigma)
x1 <- scale(X)
data <- list(data_50,data_100,data_200,data_500,x1)

numerosity <- c(50,100,200,500,1000)

for (j in 1:length(numerosity)){ # numerosity 
  
  X <- data[[j]]
  tXX <- crossprod(X)
  n <- numerosity[j]
  setwd(dirs[j])
  timing <- numeric(50)
  
  for (iter in 1:50){
    start_time <- Sys.time()
    
    result = birth_death(tXX, n, a = q, U = diag(1,q)/n, w, max_moves = 400)
    timing_sim <- diff_time(start_time)
    print(paste(timing_sim, iter, numerosity[j]))
    timing[iter] <- timing_sim
    
    saveRDS(result, file=paste0('out_sim_n_',n,'_iter_',iter,'.RDS'))
    rm(result)
    gc()
    
  }
  saveRDS(timing, 'timing.RDS')
  
}




######################
# what I have 

plot(outbd$waiting_times)
wait <- outbd$waiting_times
length(wait)
mm <- which.max(tail(wait,10))
third_dim <- dim(outbd$Graphs)[3]
mmg_index <- third_dim - 10 + mm
mmg <- outbd$Graphs[,,mmg_index]
dag <- readRDS("C:/D-drive-15734/Sara/Università/Magistrale/Secondo anno/Tesi/computational part/Adaptive Metropolis Hastings/Adaptive runs/simulations/original_dag.RDS")

ssh <- numeric(50)
for (i in 1:50){
  outbd <- readRDS(paste0('out_sim_n_100_iter_',i,'.RDS'))
  wait <- outbd$waiting_times
  mm <- which.max(wait)
  print(wait[mm])
  mmg <- outbd$Graphs[,,mm]
  ssh[i] <- causalDisco::shd(dag, mmg)
  # print(heat_tile(true, title='True DAG')+heat_tile(mmg, title = 'selected DAG'))
  # Sys.sleep(3)
  
}

boxplot(ssh)

x <- local_moves(DAG, tXX, n, a, U, w)


for (i in 29:34){
  mmg <- outbd$Graphs[,,i]
  print(heat_tile(true, title='True DAG')+heat_tile(symmetric_matrix(true), title = 'True with symmetric edges')+heat_tile(mmg, title = 'selected DAG'))
  Sys.sleep(2)
  
}
