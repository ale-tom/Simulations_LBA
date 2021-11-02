pacman::p_load('tidyverse','magrittr','rtdists')



set.seed(123)

## LBA: recovers parameters
rt1 <- rLBA(500, A=0.5, b=c(1.2,1), t0 = c(0.5,0.2), mean_v=c(2.4, 1.6), sd_v=c(1,1.2))


head(rt1)

prop.table(table(rt1$response))

lattice::densityplot(~rt, rt1, group = response, auto.key=TRUE, plot.points=FALSE, weights = rep(1/nrow(rt1), nrow(rt1)), ylab = "Density")

objective_fun <- function(par, rt, response, distribution = "norm") {
  # simple parameters
  spar <- par[!grepl("[12]$", names(par))]  
  
  # distribution parameters:
  dist_par_names <- unique(sub("[12]$", "", grep("[12]$" ,names(par), value = TRUE)))
  dist_par <- vector("list", length = length(dist_par_names))
  names(dist_par) <- dist_par_names
  for (i in dist_par_names) dist_par[[i]] <- as.list(unname(par[grep(i, names(par))]))
  dist_par$sd_v <- c(1, dist_par$sd_v) 
  
  # get summed log-likelihood:
  d <- do.call(dLBA, args = c(rt=list(rt), response=list(response), spar, dist_par, 
                              distribution=distribution, silent=TRUE))
  if (any(d == 0)) return(1e6)
  else return(-sum(log(d)))
}

objective_fun(c(A=0.5, b=1, t0=0.5, mean_v1=2.4, mean_v2=1.6, sd_v1=1.2), 
              rt=rt1$rt, response=rt1$response)

init_par <- c(runif(3, 0, 0.5), runif(3, 0.5, 2))
names(init_par) <- c("A", "b", "t0", "mean_v1", "mean_v2", "sd_v2")
nlminb(objective_fun, start = init_par, rt=rt1$rt, response=rt1$response, lower = 0)
