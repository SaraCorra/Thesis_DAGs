# Import all functions from the Scripts
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Adaptive_metropolis.R")
source("BDMCMC.R")

data_creation(n = 1000, w = 0.2,  random = F, boot = F)
# BDMCMC simulation:
result = birth_death(tXX, n, a = q, U = diag(1,q)/n, w, max_moves = 400)
final_dags <- result$Graphs
wait <- result$waiting_times
mmg_index <- which.max(wait)
print(heat_tile(true, title='True DAG')+heat_tile(symmetric_matrix(true), title = 'True with symmetric edges')+heat_tile(final_dags[,,mmg_index], title = 'selected DAG'))

# Adaptive simulation:

result = learn_DAG(S = 60000, burn = 5000, a = q, U = diag(1,q)/n, 
                   data = X, w = 0.2, fast = TRUE, 
                   save.memory = T, collapse = T, adaptive_mcmc = T,
                   adapt_factor = 0.5)


