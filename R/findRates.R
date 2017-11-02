#' Use diversitree to find ML speciation and extinction rates
#'
#' Makes an attempt to automatically calculate lambda and mu using diversitree given
#' starting values for the optimization. 
#'
#' @param tree A phylogeny in ape format.
#' @param prop.complete The proportion of all species in the tree that are currently
#' contained in it. For instance, if a tree contains 9 species, but there should be
#' 10 species in the tree, set sampling to 0.9.
#' @param ini.lambda Initial speciation value for the optimization. Defaults to 1.
#' @param ini.mu Initial extinction value for the optimization. Defaults to 0.1.
#'
#' @details Add details.
#' 
#' @return A named numeric vector with the estimated speciation and extinction rates.
#'
#' @export
#'
#' @references ETM unpublished
#'
#' @examples
#' #load a molecular tree up
#' data(bird.families)
#'
#' tree <- multi2di(bird.families)
#'
#' findRates(tree, 0.95)

findRates <- function(tree, prop.complete, ini.lambda=1, ini.mu=0.1)
{
	iniLik <- make.bd(tree, sampling.f=prop.complete)
	results <- find.mle(iniLik, method="subplex", c(ini.lambda,ini.mu))
	results$par
}