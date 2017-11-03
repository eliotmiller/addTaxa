#' Randomly drop tips from a phylogeny
#'
#' Given a phylogeny and a data frame of species and the taxonomic group to which they
#' belong, will drop requested number of species, with some caveats.
#'
#' @param tree An ape-style phylogenetic tree
#' @param groupings A data frame with two columns, "species" and "group". Missing species,
#' to be added, are taken as those that do not match a value in the tip labels of tree.
#' @param no.to.drop The number of tips to drop from the tree.
#' 
#' @details This function drops tips from a phylogeny but ensures that no genera are
#' entirely removed from the tree. No check is made to ensure that the number of tips
#' to drop can be can be accomplished without dropping below the threshold of leaving
#' one species per genus.
#'
#' @return A pruned phylogeny.
#'
#' @export
#'
#' @examples
#' #load an example tree up
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
#' #use the function to drop 20 species
#' example <- tipDropper(tree, groupsDF, 20)

tipDropper <- function(tree, groupings, no.to.drop)
{
	for(i in 1:no.to.drop)
	{
		#this figures out how many species are in each genus in the tree
		genusVector <- unlist(lapply(strsplit(tree$tip.label, "_"), "[", 1))

		#count them up
		sppCount <- plyr::count(genusVector)

		#drop any genera that only have one species in the tree in them
		genera <- sppCount$x[sppCount$freq > 1]

		#sample a genus
		genus <- sample(genera, 1)

		#sample a species in that genus
		spp <- groupings$species[groupings$group == genus]
		sp <- sample(spp, 1)

		#drop the tip from the tree
		tree <- drop.tip(tree, as.character(sp))

		#drop the tip from groupings to make sure you don't hit it again
		groupings <- groupings[groupings$species != sp,]
	}

	tree
}
