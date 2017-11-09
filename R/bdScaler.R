#' Generate missing speciation events
#'
#' Place speciation event in a speciation/extinction rate-informed position
#'
#' @param tree A phylogeny in ape format.
#' @param lambda The speciation rate.
#' @param mu The exctinction rate.
#' @param min.age Minimum age of the missing speciation event to be inferred. Note that
#' this is the age counting back from the tips, such that 0 is the present. Thus, min.age
#' must always be less than max.age.
#' @param max.age Maximum age of the missing speciation event to be inferred.
#'
#' @details This function assumes that a single species is missing, then takes the
#' provided lambda and mu values, and the range of time over which the speciation event
#' occurred, and generates the missing time. 
#' 
#' @return A missing speciation time.
#'
#' @export
#'
#' @importFrom TreeSim getx corsim
#'
#' @references ETM unpublished
#'
#' @examples
#' #load a molecular tree up
#' data(bird.families)
#'
#' bdScaler(bird.families, lambda=1, mu=0, min.age=0, max.age=10)

bdScaler <- function(tree, lambda, mu, min.age, max.age)
{
	saveRDS(tree, "sliced.RDS")

	#generate the vector of branching times in the input tree with the getx function
	originalTimes <- TreeSim::getx(tree)

	#simulate the missing branching time. having big problems with this sometimes
	#being deeper than it should be, so wrap this in a while statement. 
	acceptable <- 0
	while(acceptable == 0)
	{
		temp <- TreeSim::corsim(x=originalTimes, lambda, mu, missing=1,
			tyoung=min.age, told=max.age)

		#figure out which of the new branching times was not one of the original times
		result <- temp[!(temp %in% originalTimes)]

		if(result <= max(originalTimes) & result >= 0)
		{
			acceptable <- 1
		}
	}

	#not sure why this has a name, but weird formatting, so remove to make it nicer
	names(result) <- NULL
	result
}
