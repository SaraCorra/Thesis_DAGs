
set_creation <- function(DAG){
  
  A <- DAG
  q <- ncol(A)
  A_na <- A
  diag(A_na) <- NA
  
  id_set = c()
  dd_set = c()
  rd_set = c()
  
  ## set of nodes for id
  set_id = which(A_na == 0, TRUE)
  # If I want to add I need nothing at that location
  if(length(set_id) != 0){
    id_set = cbind(1, set_id)
  }
  
  ## set of nodes for dd
  set_dd = which(A_na == 1, TRUE)
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }
  
  ## set of nodes for rd
  set_rd = which(A_na == 1, TRUE)
  if(length(set_rd != 0)){
    rd_set = cbind(3, set_rd)
  }
  
  O = rbind(id_set, dd_set, rd_set)
  
  return(O)
}


learn_DAG <- function(S, burn,
                      data, a, U, w,
                      fast = FALSE, save.memory = FALSE, collapse = FALSE,
                      verbose = TRUE, adaptive_mcmc=FALSE, adapt_factor = 0.1) {
  
  input <- as.list(environment())
  
  ## Input check
  
  data_check <- sum(is.na(data)) == 0
  S.burn_check <- is.numeric(c(S,burn)) & length(S) == 1 & length(burn) == 1
  S.burn_check <- if (S.burn_check) {
    S.burn_check & (S %% 1 == 0) & (burn %% 1 == 0) # verify if is.integer() can be used
  } else {
    S.burn_check
  }
  a_check <- is.numeric(a) & (length(a) == 1) & (a > ncol(data) - 1)
  w_check <- is.numeric(w) & (length(w) == 1) & (w <= 1) & (w >= 0)
  U_check <- is.numeric(U) & (dim(U)[1] == dim(U)[2]) & (prod(eigen(U)$values) > 0) & isSymmetric(U)
  U.data_check <- dim(U)[1] == ncol(data)
  
  if (data_check == FALSE) {
    stop("Data must not contain NAs")
  }
  if (S.burn_check == FALSE) {
    stop("S and burn must be integer numbers")
  }
  if (a_check == FALSE) {
    stop("a must be at least equal to the number of variables")
  }
  if (w_check == FALSE) {
    stop("w must be a number between 0 and 1")
  }
  if (U_check == FALSE) {
    stop("U must be a squared symmetric positive definite matrix")
  }
  if (U.data_check == FALSE) {
    stop("U must be a squared spd matrix with dimensions equal to the number of variables")
  }
  
  n.iter <- input$burn + input$S
  X <- scale(data, scale = FALSE)
  tXX <- crossprod(X)
  
  n <- dim(data)[1]
  q <- dim(data)[2]
  
  ## Initialize arrays or vectors depending on save.memory
  if (save.memory == TRUE) {
    # if I save memory I just need to initialize a vector
    Graphs <- vector("double", n.iter) # zero times n.iter
    L <- vector("double", n.iter)
    D <- vector("double", n.iter)
  } else {
    Graphs <- array(0, dim = c(q,q,n.iter))
    L <- array(0, dim = c(q,q,n.iter))
    D <- array(0, dim = c(q,q,n.iter))
    # these are the (q,q,s)
  }
  
  currentDAG <- matrix(0, ncol = q, nrow = q)
  current_node <- 1
  current_table <- set_creation(currentDAG)
  current_table <- cbind(current_table, rep(1/nrow(current_table), nrow(current_table)))
  adapt_probs <- matrix(1/q^2, q, q)
  alphavec <- numeric(n.iter)
  
  ## iterations
  
  if (save.memory == FALSE) {
    type = "collapsed"
    if (verbose == TRUE) {
      cat("Sampling DAGs...")
      pb <- utils::txtProgressBar(min = 2, max = n.iter, style = 3)
    }
    
    for (i in 1:n.iter) {
      
      if (input$adaptive_mcmc & i!=1){
        table <- current_table
        adapt_probs <- prop$adapt_probs
      } else if (input$adaptive_mcmc & i==1){
        table <- current_table
      }else {
        table <- NULL
        adapt_probs <- NULL
      }
      
      prop <- propose_DAG(
        currentDAG,
        input = input,
        iter = i,
        proposal_prob = table, 
        adapt_probs = adapt_probs
      )
      
      if (input$adaptive_mcmc){
        proposed_table <- prop$proposed_proposalprob
      } else {
        proposed_table <- NULL
      }
      ratio_eval <- acceptreject_DAG(tXX, n,currentDAG, prop$proposedDAG,
                                     prop$op.node, prop$op.type, a, U, w,
                                     prop$current.opcard,
                                     prop$proposed.opcard,
                                     current_nodeindex = current_node,
                                     proposed_nodeindex = prop$proposed_nodeindex, # always exists
                                     current_proptable = current_table, 
                                     proposed_proptable = proposed_table,
                                     adaptive_mcmc = input$adaptive_mcmc)
      current_node <- ratio_eval$current_nodeindex
      if (input$adaptive_mcmc){
        current_table <- ratio_eval$current_proptable
      } else {
        current_table <- NULL
      }
      
      if (ratio_eval$is.accepted == TRUE) {
        currentDAG <- prop$proposedDAG
      }
      
      alphavec[i] <- ratio_eval$alpha
      Graphs[,,i] <- currentDAG
      if (verbose == TRUE) {
        utils::setTxtProgressBar(pb, i)
        close(pb)
      }
    }
    if (collapse == FALSE) {
      type = "complete"
      if (verbose == TRUE) {
        cat("\nSampling parameters...")
        pb <- utils::txtProgressBar(min = 2, max = n.iter, style = 3)
      }
      for (i in 1:n.iter) {
        postparams <- rDAGWishart(1, Graphs[,,i], a+n, U+tXX)
        L[,,i] <- postparams$L
        D[,,i] <- postparams$D
        if (verbose == TRUE) {
          utils::setTxtProgressBar(pb, i)
          close(pb)
        }
      }
    }
    Graphs <- Graphs[,,(burn+1):n.iter]
    L <- L[,,(burn+1):n.iter]
    D <- D[,,(burn+1):n.iter]
    # subsets the Graphs array to retain only the sampled DAGs after the 
    # burn-in period.
    
  } else { # save memory is true 
    
    type = "compressed and collapsed"
    if (verbose == TRUE) {
      cat("Sampling DAGs...")
      pb <- utils::txtProgressBar(min = 2, max = n.iter, style = 3)
    }
    for (i in 1:n.iter) {
      
      # start_time <- Sys.time()
      
      if (input$adaptive_mcmc & i!=1){
        table <- current_table
        adapt_probs <- prop$adapt_probs
      } else if (input$adaptive_mcmc & i==1){
        table <- current_table
      }else {
        table <- NULL
        adapt_probs <- NULL
      }
      
      prop <- propose_DAG(
        currentDAG,
        input = input,
        iter = i,
        proposal_prob = table, 
        adapt_probs = adapt_probs
      )
      
      # print(paste('end proposal', diff_time(start_time)))
      
      # start_time <- Sys.time()
      
      if (input$adaptive_mcmc){
        proposed_table <- prop$proposed_proposalprob
      } else {
        proposed_table <- NULL
      }
      
      ratio_eval <- acceptreject_DAG(tXX, n,currentDAG, prop$proposedDAG,
                                     prop$op.node, prop$op.type, a, U, w,
                                     prop$current.opcard,
                                     prop$proposed.opcard,
                                     current_nodeindex = current_node,
                                     proposed_nodeindex = prop$proposed_nodeindex, # always exists
                                     current_proptable = current_table, 
                                     proposed_proptable = proposed_table,
                                     adaptive_mcmc = input$adaptive_mcmc)
      current_node <- ratio_eval$current_nodeindex
      if (input$adaptive_mcmc){
        current_table <- ratio_eval$current_proptable
      } else {
        current_table <- NULL
      }
      
      # print(paste('end ratio', diff_time(start_time)))
      
      if (ratio_eval$is.accepted == TRUE) {
        currentDAG <- prop$proposedDAG
      }
      
      alphavec[i] <- ratio_eval$alpha
      Graphs[i] <- bd_encode(currentDAG)
      # here I also encode!!!  
      if (verbose == TRUE) {
        utils::setTxtProgressBar(pb, i)
        close(pb)
      }
    }
    if (collapse == FALSE) {
      type = "compressed"
      if (verbose == TRUE) {
        cat("\nSampling parameters...")
        pb <- utils::txtProgressBar(min = 2, max = n.iter, style = 3)
      }
      for (i in 1:n.iter) {
        postparams <- rDAGWishart(1, bd_decode(Graphs[i]), a+n, U+tXX)
        L[i] <- bd_encode(postparams$L)
        D[i] <- bd_encode(postparams$D)
        if (verbose == TRUE) {
          utils::setTxtProgressBar(pb, i)
          close(pb)
        }
      }
    }
    Graphs <- utils::tail(Graphs, S)
    # here I get rid of the burnin 
    L <- utils::tail(L, S)
    D <- utils::tail(D, S)
    alphavec <- utils::tail(alphavec, S)
  }
  
  if (collapse == FALSE) {
    out <- new_bcdag(list(Graphs = Graphs, L = L, D = D), input = input, type = type)
  } else {
    out <- new_bcdag(list(Graphs = Graphs), input = input, type = type)
  }
  return(list(out, alpha = alphavec))
}


#################################################

propose_DAG <- function(DAG, input, iter, proposal_prob, adapt_probs) {
  
  A <- DAG
  q <- ncol(A)
  
  if (!input$adaptive_mcmc){ proposal_prob <- set_creation(A)}
  
  # Sample one random operator and verify it produces a DAG
  
  if (input$fast == FALSE) {
    # I take all the valid operations and randomly select one 
    proposed.opcardvec <- vector(length = nrow(proposal_prob))
    # I will have the proposed DAGs
    for (i in 1:nrow(proposal_prob)) {
      proposed.opcardvec[i] <- gRbase::is.DAG(operation(proposal_prob[i,1], DAG, proposal_prob[i,2:3]))
    }
    # fill in with the TRUE or FALSE
    proposed.opcard <- sum(proposed.opcardvec)
    # how many true proposed DAGs I have
    if (input$adaptive_mcmc) { 
      i <- sample(which(proposed.opcardvec), 1, 
                  prob = proposal_prob[which(proposed.opcardvec),4])
      # random sample one of the true DAGs
    } else { # adaptive_mcmc=FALSE
      i <- sample(which(proposed.opcardvec), 1)
    }
    A_next <- operation(proposal_prob[i,1], A, proposal_prob[i,2:3])
    # do the actual operation on the selected one 
    current.opcard <- get_opcard(A_next)
    # how many proper operations I can perform starting from the new DAG
  } else { # fast=TRUE
    repeat {  # this is a sort of while 
      if (input$adaptive_mcmc) { 
        i <- sample(nrow(proposal_prob), 1, prob = proposal_prob[,4])
        # random sample one from the total then check
      } else { # adaptive_mcmc == F
        i <- sample(nrow(proposal_prob), 1)
      }
      # randomly sample one operation in proportion to the cardinality 
      A_next <- operation(proposal_prob[i,1], A, proposal_prob[i,2:3])
      verify <- gRbase::is.DAG(A_next)
      
      if (verify == TRUE) {
        break
      }
    }
    proposed.opcard <- nrow(proposal_prob) # I don't understand why not getting 
    # the ones defining an actual DAG
    current.opcard <- nrow(proposal_prob)
    # shouldn't it relate to the new selected DAG somehow?
    
  } # this close the else statement
  
  op.type <- proposal_prob[i,1]
  if (op.type == 3) {
    op.node <- proposal_prob[i,-1]
  } else {
    op.node <- proposal_prob[i,3]
    # Therefore, only the child node (the node to which the edge points) is 
    # relevant for determining op.node in the context of an edge reversal 
    # operation. The parent node remains the same in this operation.
  }
  
  if (input$adaptive_mcmc){
    
    proposed_nodeindex <- i
    
    proposal_prob <- set_creation(A_next)
    adapt_probs <- adapt_probs*((iter-1)/(iter)) + A_next/(iter)
    post_prob <- sapply(1:nrow(proposal_prob), function(x) adapt_probs[proposal_prob[x,2],proposal_prob[x,3]])
    # post_prob <- adapt_probs[proposal_prob[, 2], proposal_prob[, 3]]
    proposal_probabilities <- ifelse(proposal_prob[,1] == 1, (post_prob+input$adapt_factor)/(1-post_prob+input$adapt_factor),
                                     ifelse(proposal_prob[,1] == 2, (1-post_prob+input$adapt_factor)/(post_prob+input$adapt_factor),
                                            ifelse(proposal_prob[,1] == 3, 1, NA)))
    proposal_prob <- cbind(proposal_prob, proposal_probabilities)
    
  }
  
  if (input$adaptive_mcmc){
    
    return(list(proposedDAG = A_next, op.type = op.type, op.node = op.node,
                current.opcard = current.opcard,
                proposed.opcard = proposed.opcard, 
                proposed_nodeindex = proposed_nodeindex,
                proposed_proposalprob = proposal_prob,
                adapt_probs = adapt_probs))
  } else {
    
    return(list(proposedDAG = A_next, op.type = op.type, op.node = op.node,
                current.opcard = current.opcard,
                proposed.opcard = proposed.opcard))
    
  }
  
}

#################################################

acceptreject_DAG <- function(tXX, n, currentDAG, proposedDAG, node, op.type,
                             a, U, w, current.opcard, proposed.opcard,
                             current_nodeindex, proposed_nodeindex, 
                             current_proptable, proposed_proptable, 
                             adaptive_mcmc) {
  
  logprior.ratios <- c(log(w/(1-w)), log((1-w)/w), log(1))
  logprior.ratio <- logprior.ratios[op.type]
  # select the ratio depending on the operation type (delete, remove, invert)
  
  if (adaptive_mcmc){
    current_prop <- current_proptable[current_nodeindex,4]
    proposed_prop <- proposed_proptable[proposed_nodeindex,4]
    logproposal.ratio <- log(current_prop) - log(proposed_prop)
  } else {
    logproposal.ratio <- log(current.opcard) - log(proposed.opcard)
    # this is the correction for symmetry in uniform case 
  }
  
  if (op.type != 3) {
    current_lml <- DW_nodelml(node, currentDAG, tXX, n, a, U)
    proposed_lml <- DW_nodelml(node, proposedDAG, tXX, n, a, U)
  } else {
    current_lml <- DW_nodelml(node[1], currentDAG, tXX, n, a, U) +
      DW_nodelml(node[2], currentDAG, tXX, n, a, U)
    proposed_lml <- DW_nodelml(node[1], proposedDAG, tXX, n, a, U) +
      DW_nodelml(node[2], proposedDAG, tXX, n, a, U)
  }
  
  
  acp.ratio <- min(0, proposed_lml - current_lml + logprior.ratio +
                     logproposal.ratio)
  
  is.accepted <- log(stats::runif(1)) < acp.ratio
  if (is.accepted){
    current_nodeindex = proposed_nodeindex
    current_proptable <- proposed_proptable
  } 
  # this is the problem 
  # need to recreate the function with the burnin 
  
  return(list(is.accepted = is.accepted, current_proptable = current_proptable,
              current_nodeindex = current_nodeindex, alpha = acp.ratio))
}

#######################################

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

#########################

bd_encode <- function(matrix, separator = ";") {
  paste(matrix, collapse = separator)
}

#########################

bd_decode <- function(string, separator = ";") {
  vec4mat <- as.numeric(strsplit(string, separator)[[1]])
  q <- length(vec4mat)
  matrix(vec4mat, ncol = sqrt(q))
}

##########################

new_bcdag <- function(x = list(), input = list(), type = "complete") {
  stopifnot(is.list(x))
  stopifnot(is.list(input))
  type <- match.arg(type, c("complete", "compressed", "collapsed", "compressed and collapsed"))
  
  structure(x,
            class = "bcdag",
            type = type,
            input = input)
}

#######################################

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

##################################

get_opcard <- function(DAG) {
  A <- DAG
  q <- ncol(A)
  A_na <- A
  diag(A_na) <- NA
  
  # Define the set of possible operations!
  
  id_set = c()
  dd_set = c()
  rd_set = c()
  
  ## set of nodes for id
  set_id = which(A_na == 0, TRUE)
  if(length(set_id) != 0){
    id_set = cbind(1, set_id)
  }
  
  ## set of nodes for dd
  set_dd = which(A_na == 1, TRUE)
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }
  
  ## set of nodes for rd
  set_rd = which(A_na == 1, TRUE)
  if(length(set_rd != 0)){
    rd_set = cbind(3, set_rd)
  }
  
  O = rbind(id_set, dd_set, rd_set)
  op.cardvec <- vector(length = nrow(O))
  for (i in 1:nrow(O)) {
    op.cardvec[i] <- gRbase::is.DAG(operation(O[i,1], DAG, O[i,2:3]))
  }
  op.card <- sum(op.cardvec)
  return(op.card)
}

################################

DW_nodelml <- function(node, DAG, tXX, n, a, U) {
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

################################

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

############################

pa <- function(node, DAG) {
  pa <- which(DAG[,node] != 0)
  return(pa)
}

###################################

  