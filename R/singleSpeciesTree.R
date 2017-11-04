#' Create a single-species tree
#'
#' Generate a single-species tree with appropriate stem branch length
#'
#' @param tree An ape-style phylogenetic tree.
#' @param tip A character vector corresponding to one of the species in the phylogeny.
#' 
#' @details Ape cannot readily deal with single-species trees. This functions will extract
#' a single-species tree with appropriate branch length to the most recent common
#' ancestor. The tree cannot be plotted, but it can be used for analyses.
#'
#' @return A single-species phylo object.
#'
#' @export
#'
#' @importFrom ape read.tree
#'
#' @references Mast et al. 2015. Paraphyly changes understanding of timing and tempo of 
#' diversification in subtribe Hakeinae (Proteaceae), a giant Australian plant radiation.
#' American Journal of Botany.
#'
#' @examples
#' #load a molecular tree up
#' data(bird.families)
#'
#' #extract a single-species tree
#' tree <- singleSpeciesTree(bird.families, "Rheidae")

singleSpeciesTree <- function(tree, tip)
{
	#determine the edge length to that tip. set up a dummy frame for this
	temp <- data.frame(id = 1:length(tree$tip.label), species=tree$tip.label)
	tipID <- temp$id[temp$species == tip]
	#set up another dummy frame to determine what the right edge length is leading to tip
	temp2 <- data.frame(tree$edge)
	row.names(temp2) <- 1:dim(temp2)[1]
	#determine the right edge length
	theEdge <- as.numeric(row.names(temp2[temp2[,2]==tipID,]))
	#get the actual	length of the edge
	edgeLength <- tree$edge.length[theEdge]
	#create a silly Newick file in words
	inWords <- paste("(", tip, ":", edgeLength, ")", ";", sep="")
	singleSpTree <- read.tree(text = inWords)
	return(singleSpTree)
}
