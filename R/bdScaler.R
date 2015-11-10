#' Generate missing speciation events
#'
#' Place speciation event in a speciation/extinction rate-informed position
#'
#' @param tree A phylogeny in ape format.
#' @param lambda The speciation rate.
#' @param mu The exctinction rate.
#' @param missing Number of missing speciation events to be inferred.
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

bdScaler <- function(tree, lambda, mu, missing, min.age, max.age)
{
	#generate the vector of branching times in the input tree with the getx function
	originalTimes <- getx(tree)
	
	#create a vector of length zero that are the possible speciation times
	possTimes <- c()

	#create a simple placeholder here to tally the number of while loops. if it goes
	#through 10 loops and still has not found a solution, have it use a uniform
	#distribution to generate the value
	placeH <- 1
	
	while(length(possTimes) == 0)
	{
		#use corsim to simulate the missing speciation events in the appropriate age range
		allTimes <- corsim(x=originalTimes, lambda=lambda, mu=mu, missing=missing,
			told=0, tyoung=0)

		#figure out which of the new branching times is different than the input
		newTimes <- setdiff(allTimes, originalTimes)

		#subset these new times to those that are between min and max age. if the length
		#of possTimes is not zero it should break out of the while loop
		possTimes <- newTimes[newTimes > min.age & newTimes < max.age]
		
		if(placeH > 10)
		{
			possTimes <- runif(n=1, min=min.age, max=max.age)
		}
		
		#add one to placeH
		placeH <- placeH + 1
	}
	
	#sample a single new branching time from that. it turns out after much time spent
	#troubleshooting that if you sample a vector of length 1 weird things happen
	if(length(possTimes) == 1)
	{
		output <- possTimes
	}
	else
	{
		output <- sample(possTimes, size=1)
	}

	return(output)
}
