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
#' @param rate.estimate Whether to use 'laser', 'ape', or 'diversitree' to
#' calculate diversification rates. The latter is the only of those that can account for
#' missing taxa, but in my experience returns systematically biased values.
#'
#' @details Add details.
#' 
#' @return A named numeric vector with the estimated speciation and extinction rates.
#'
#' @export
#'
#' @importFrom diversitree make.bd find.mle
#' @importFrom ape is.binary.tree multi2di birthdeath
#' @importFrom laser bd
#' @importFrom phytools bd
#'
#' @references ETM unpublished
#'
#' @examples
#' #load a molecular tree up
#' data(bird.families)
#'
#' tree <- multi2di(bird.families)
#'
#' findRates(tree=tree, prop.complete=0.95, rate.estimate="diversitree")

findRates <- function(tree, prop.complete, ini.lambda=1, ini.mu=0.1, rate.estimate)
{
	if(rate.estimate=="laser")
	{
		#pull the branching times, pass to laser
		temp <- laser::bd(TreeSim::getx(tree))

		#make some quick calculations of lambda and mu
		tempLambda <- temp$r - (temp$r * temp$a)
		tempMu <- temp$r * temp$a

		finalResults <- c(tempLambda, tempMu)
		names(finalResults) <- c("lambda","mu")
	}

	else if(rate.estimate=="ape")
	{
		finalResults <- phytools::bd(ape::birthdeath(tree))
		names(finalResults) <- c("lambda","mu")
	}

	else if(rate.estimate=="diversitree")
	{
		#if the tree is not binary then this function will fail. if that's the case,
		#convert polytomies to bifurcating branches of length 0. note that this does
		#not mean the actual tree ends up with polytomies
		if(!ape::is.binary.tree(tree))
		{
			tree <- ape::multi2di(tree)
		}
		iniLik <- diversitree::make.bd(tree, sampling.f=prop.complete)
		results <- diversitree::find.mle(iniLik, method="subplex", c(ini.lambda,ini.mu))
		finalResults <- results$par
	}

	else
	{
		stop("rate.estimate must be set to one of 'laser', 'ape', or 'diversitree'")
	}
	finalResults
}
