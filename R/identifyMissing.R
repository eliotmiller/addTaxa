#' Identify taxa that are not included in phylogeny
#'
#' Given a data frame of species and taxonomic assignments, and an accepted phylogeny with
#' some of those species in it, identify which species are not included in the phylogeny.
#'
#' @param tree An ape-style phylogenetic tree
#' @param groupings A data frame with two columns, "species" and "group". Missing species,
#' to be added, are taken as those that do not match a value in the tip labels of tree.
#' @param print.to.screen Logical. Will print a list of the missing taxa to screen.
#' 
#' @details Utility function to identify taxa in the groupings data frame that are not
#' included in the input phylogeny.
#'
#' @return A data frame of missing taxa with two columns, species and group. 
#'
#' @export
#'
#' @references Eliot Miller unpublished
#'
#' @examples
#' #load a molecular tree up
#' data(bird.families)
#'
#' #create a data frame of all taxa from the phylogeny, and make up species groups
#' #for each.
#' dummy.frame <- data.frame(species=bird.families$tip.label, 
#' group=c(rep("nonPasserine", 95), rep("suboscine", 9), rep("basalOscine", 13), 
#' rep("oscine", 20)))
#' 
#' #now make up a few passerine taxa that are missing and add these into the dummy frame
#' toAdd <- data.frame(species=c("Icteria", "Yuhina", "Pityriasis", "Macgregoria"), 
#' group=c(rep("oscine", 2), rep("basalOscine", 2)))
#' groupsDF <- rbind(dummy.frame, toAdd)
#'
#' #use the function to see what was missing from the original tree as compared with the 
#' #new data frame
#' identifyMissing(bird.families, groupsDF, print.to.screen=TRUE)

identifyMissing <- function(tree, groupings, print.to.screen)
{
	#set the species and group to character vectors
	groupings$species <- as.character(groupings$species)
	groupings$group <- as.character(groupings$group)

	#subset the groupings data frame to those spp we only have taxonomic information for
	results <- groupings[!(groupings$species %in% tree$tip.label), ]

	if(print.to.screen==TRUE)
	{
		print("You are missing phylogenetic information for these species:", quote=FALSE)
		print(results$species, quote=FALSE)
	}
	
	results
}
