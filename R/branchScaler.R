#' Scale a vector of branch lengths
#'
#' Simple utility function to scale tree-wide distribution of branch lengths
#'
#' @param input.vector A vector of numbers such as the length of all the branches in
#' a phylogeny.
#' @param max.age The maximum value you want to occur in the vector.
#'
#' @details This function simply divides the input vector, e.g. a distribution of all the
#' edge lengths in a tree, by the maximum length in that vector, then multiplies the
#' resulting values by the "max.age" argument. 
#' 
#' @return A scaled vector.
#'
#' @export
#'
#' @references Mast et al. 2015. Paraphyly changes understanding of timing and tempo of 
#' diversification in subtribe Hakeinae (Proteaceae), a giant Australian plant radiation.
#' American Journal of Botany.
#'
#' @examples
#' #load a molecular tree up
#' data(bird.families)
#'
#' branchScaler(bird.families$edge.length, 2)

branchScaler <- function(input.vector, max.age)
{
	scaled <- input.vector/max(input.vector)

	output.vector <- scaled * max.age

	return(output.vector)
}
