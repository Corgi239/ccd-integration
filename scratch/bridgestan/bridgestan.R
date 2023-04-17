# https://github.com/roualdes/bridgestan/tree/main/R
# https://roualdes.github.io/bridgestan/languages/r.html
devtools::install_github("https://github.com/roualdes/bridgestan", subdir="R")

library(bridgestan)

# Clone bridgestan repo with
# `git clone --recurse-submodules https://github.com/roualdes/bridgestan.git` 

# From the bridgestan directory, run 
# `make ~/github/ccd-integration/scratch/eight_schools_noncentered_model.so`
# to generate the .so file
model <- StanModel$new(
    "scratch/eight_schools_noncentered_model.so", 
    "scratch/eight_schools.json", 1234
)

model$log_density_gradient(rep(0., 10))
model$log_density(rep(0., 10))
