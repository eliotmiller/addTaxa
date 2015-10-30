#' Identify the node that subtends each named clade
#'
#' Given a phylogeny and a data frame of species and named clades, identify the nodes that
#' subtend those named clades. 
#'
#' @param tree An ape-style phylogenetic tree
#' @param clade.membership An optional data frame with first column = "species", second
#' column = "clade". These are named, monophyletic clades to which the species in the
#' input phylogeny can belong. Not every species in the input phylogeny needs to be in the
#' clade membership data frame, but missing species added to species not in this frame
#' will not be included in the output clade membership data frame. They will of course
#' still be included in the output phylogenies.
#' 
#' @details Given a phylogeny and a data frame of species and named clades, identify the 
#' nodes that subtend those named clades. If a clade is composed of a single species it
#' will return the tip number of the clade
#'
#' @return A list of named elements. Each element is a single number (a node). Names
#' correspond to the name of the clade. 
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
#' #create a data frame of all taxa from the phylogeny, and make up clade memberships
#' #for each. note that the names on this data frame differ from "groupings" in other
#' #functions
#' dummy.frame <- data.frame(species=bird.families$tip.label, 
#' clade=c(rep("nonPasserine", 95), rep("suboscine", 9), rep("basalOscine", 13), 
#' rep("oscine", 20)))
#'
#' #use the function
#' getMRCAs(bird.families, dummy.frame)

getMRCAs <- function(tree, clade.membership)
{
	clades <- unique(clade.membership$clade)
	
	results <- list()
	
	for(i in 1:length(clades))
	{
		#important that this is a character vector
		cladeSpp <- as.character(clade.membership$species[clade.membership$clade == clades[[i]]])
		#the function getMRCA does not work for single species clades, so need to just get
		#node of this species
		if(length(cladeSpp) == 1)
		{
			#identify the node number of the single species
			results[[i]] <- which(tree$tip.label==cladeSpp)
		}
		else
		{
			#find the MRCA
			results[[i]] <- getMRCA(tree, cladeSpp)
		}
	}
	names(results) <- clades
	results
}
