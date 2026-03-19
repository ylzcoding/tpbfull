library(rstan)

rstan_options(auto_write = TRUE)

cat("Compiling Stan model on this machine architecture...\n")

model_hc <- stan_model(file = "tpb_with_cauchy.stan")
model_gamma <- stan_model(file = "tpb_with_gamma.stan")

saveRDS(model_hc, file = "tpb_with_cauchy.rds")
saveRDS(model_gamma, file = "tpb_with_gamma.rds")

cat("Success! Model saved \n")