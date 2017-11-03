#' Identify taxa that are not included in phylogeny
#'
#' Given a data frame of species and taxonomic assignments, and an accepted phylogeny with
#' some of those species in it, identify which species are not included in the phylogeny.
#'
#' @param tree An ape-style phylogenetic tree
#' @param groupings A data frame with two columns, "species" and "group". Missing species,
#' to be added, are taken as those that do not match a value in the tip labels of tree.
#' 
#' @details Utility function to identify taxa in the groupings data frame that are not
#' included in the input phylogeny.
#'
#' @return A data frame of missing taxa with two columns, species and group. 
#'
#' @export
#'
#' @references Mast et al. 2015. Paraphyly changes understanding of timing and tempo of 
#' diversification in subtribe Hakeinae (Proteaceae), a giant Australian plant radiation.
#' American Journal of Botany.
#'
#' @examples
#' data(chelonia)
#' tree <- chelonia$phy
#'
#' #some species in this tree are identified to subspecies. drop those
#' temp <- lapply(strsplit(tree$tip.label, "_"), length)
#' names(temp) <- tree$tip.label
#' temp <- temp[temp==2]
#' tree <- drop.tip(tree, setdiff(tree$tip.label, names(temp)))
#'
#' #create an example groupings data frame.
#' groupsDF <- data.frame(species=tree$tip.label)
#' groupsDF$group <- unlist(lapply(strsplit(tree$tip.label, "_"), "[", 1))
#'
#' #use the function to drop 100 species (that is over half the species)
#' example <- tipDropper(tree, groupsDF, 100)
#'
#' #identify missing species
#' identifyMissing(example, groupsDF)

identifyMissing <- function(tree, groupings)
{
	#set the species and group to character vectors
	groupings$species <- as.character(groupings$species)
	groupings$group <- as.character(groupings$group)

	#subset the groupings data frame to those spp we only have taxonomic information for
	results <- groupings[!(groupings$species %in% tree$tip.label), ]
	
	results
}
