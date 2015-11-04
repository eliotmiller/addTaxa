#' Generate a speciation event
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
#' @references ETM unpublished
#'
#' @examples
#' #load a molecular tree up
#' data(bird.families)
#'
#' bdScaler(bird.families, lambda=1, mu=0, min.age=0, max.age=10)

bdScaler <- function(tree, lambda, mu, min.age, max.age)
{
	#generate the vector of branching times in the input tree with the getx function
	times <- getx(tree)
	
	#use corsim to simulate a missing speciation event in the appropriate age range
	timesPlus <- corsim(x=times, lambda=lambda, mu=mu, missing=1, told=max.age,
		tyoung=min.age)

	#figure out which of the new branching times is different than the input, use that
	output <- setdiff(timesPlus, times)

	return(output)
}
