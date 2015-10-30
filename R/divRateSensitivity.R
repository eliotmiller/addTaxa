#' Calculate sensitivity of absolute diversification rates to missing taxa
#'
#' Given a data frame of species and taxonomic assignments, an accepted phylogeny with
#' some of those species in it, and a number of named, monophyletic clades, will add
#' missing species next to a taxonomic relatives and quantify sensitivity of
#' clade-specific absolute diversification rates.
#'
#' @param tree An ape-style phylogenetic tree
#' @param groupings A data frame with two columns, "species" and "group". Missing species,
#' to be added, are taken as those that do not match a value in the tip labels of tree.
#' @param from.node Whether species should be added in a "polytomy" with, "crown" (more 
#' recently diverged), "stem" (previously diverged), or "randomly" crown-wards or
#' stem-wards from the tip they are bound to.
#' @param clade.membership An optional data frame with first column = "species", second
#' column = "clade". These are named, monophyletic clades to which the species in the
#' input phylogeny can belong. Not every species in the input phylogeny needs to be in the
#' clade membership data frame, but missing species added to species not in this frame
#' will not be included in the output clade membership data frame. They will of course
#' still be included in the output phylogenies.
#' @param crown.can.move Logical. If TRUE, and if missing taxa are to be added stemwards,
#' this will allow the age of the crown group to potentially shift back in time. If FALSE,
#' and if missing taxa were to be added stemwards, this will prevent taxa from being added
#' below the crown group. It will force these to be bound crownwards if the "basal"
#' taxon in the clade is selected as the species to be bound to. Note that the argument
#' name is somewhat misleading in that the crown ages can still shift when missing taxa
#' are added. Specifically, if missing taxa are bound to single-species clades, then the
#' crown age will "shift" forward in time (the single-taxon clade did not actually have a 
#' crown age).
#' @param calc.from.crown Whether the crown or stem method of bd.ms from geiger will be
#' used.
#' @param epsilon The extinction rate parameter for bd.ms
#' @param no.trees The number of complete trees across which the user would like to test
#' the diversification rate sensitivity. 
#' 
#' @details Adds missing taxa according to all rules and details outlined in
#' randomlyAddTaxa(). During this process, it tabulates to which named clades missing
#' species are added, and calculates the diversification rate of each clade after each
#' complete tree is created. In this way, the sensitivity of these absolute
#' diversification rates to missing taxa can be assessed.
#'
#' @return A dataframe where each row corresponds to a complete tree, and each column
#' corresponds to the calculated absolute diversification rate of a named clade.
#'
#' @export
#'
#' @references Mast et al. 2015. Paraphyly changes understanding of timing and tempo of 
#' diversification in subtribe Hakeinae (Proteaceae), a giant Australian plant radiation.
#' American Journal of Botany.
#'
#' @examples
#' #load a molecular tree up. resolve polytomies
#' data(bird.families)
#' bird.families <- multi2di(bird.families)
#'
#' #create a data frame of all taxa from the phylogeny, and make up species groups
#' #for each.
#' dummyGroups <- data.frame(species=bird.families$tip.label, 
#' group=c(rep("nonPasserine", 95), rep("suboscine", 9), rep("basalOscine", 13), 
#' rep("oscine", 20)))
#' 
#' #now make up a few passerine taxa that are missing and add these into the dummy frame
#' toAdd <- data.frame(species=c("Icteria", "Yuhina", "Pityriasis", "Macgregoria"), 
#' group=c(rep("oscine", 2), rep("basalOscine", 2)))
#' groupsDF <- rbind(dummyGroups, toAdd)
#'
#' #these groups were actually monophyletic. but make a slightly more detailed clade
#' #membership frame to see how one would be used. note that it doesn't include the
#' #missing taxa, and taxonomy does not exactly follow modern understanding
#' cladesDF <- data.frame(species=bird.families$tip.label, 
#' clade=c(rep("nonPasserine", 95), rep("suboscine", 9), rep("basalOscineOther", 6),
#' rep("basalOscineCore", 7), rep("oscineBase", 14), rep("oscineDerived", 6)))
#'
#' #examples of changing the from.node argument. you can plot or write these trees out to
#' #better see what the differences are
#' test <- divRateSensitivity(tree=bird.families, groupings=groupsDF, from.node="crown", 
#' 	clade.membership=cladesDF, crown.can.move=TRUE, calc.from.crown=TRUE, epsilon=0,
#'	no.trees=10)

divRateSensitivity <- function(tree, groupings, from.node, clade.membership,
	crown.can.move, calc.from.crown, epsilon, no.trees)
{
	#generate your complete trees and updated tables of to which clade each sp belongs
	treesNtables <- randomlyAddTaxa(tree=tree, groupings=groupings, 
		from.node=from.node, clade.membership=clade.membership, crown.can.move, 
		no.trees=no.trees, print.to.screen=FALSE)
	#lapply your getMRCAs function over these complete trees and tables. you will get a
	#list of lists of most recent common ancestors for each named clade
	mrcas <- lapply(1:length(treesNtables[[1]]), function(x) 
		getMRCAs(tree=treesNtables[[1]][[x]], clade.membership=treesNtables[[2]][[x]]))
	#now use these actual observed taxa numbers per clade to calculate
	#diversification rates per each named clade. make a list of lists of pruned trees.
	#if calculating from the stem is chosen, then even single species clades can get
	#a rate, though bd.ms will throw a warning about how it used a crown method (fairly
	#sure that is not how best to interpret what it did, so just suppressWarnings())
	
	prunedTrees <- list()
	
	results <- list()
	
	if(calc.from.crown==FALSE)
	{
		for(i in 1:no.trees)
		{
			prunedTrees[[i]] <- extractClade(treesNtables[[1]][[i]], mrcas[[i]], 
				root.edge=1)
		}
		for(i in 1:no.trees)
		{
			#need to correctly add a suppress warnings here
			results[[i]] <- lapply(prunedTrees[[i]], bd.ms, epsilon=epsilon, 
				crown=calc.from.crown)
		}
	}
	
	else if(calc.from.crown==TRUE)
	{
		for(i in 1:no.trees)
		{
			prunedTrees[[i]] <- extractClade(treesNtables[[1]][[i]], mrcas[[i]], 
				root.edge=0)
		}
		for(i in 1:no.trees)
		{
			results[[i]] <- lapply(prunedTrees[[i]], bd.ms, epsilon=epsilon, 
				crown=calc.from.crown)
			#note that there is no such thing as a crown calculation for a clade with only
			#a single species (bd.ms assumed these single species clades with stems should
			#be calculated "stem"-wise). Set these values to NA. First identify the clades
			#that have only 1 species in them. Here you subset the names of the pruned
			#trees to those trees that return true to the length of tip label being 
			#less than or equal to 1 
			problemClades <- names(prunedTrees[[i]])[lapply(prunedTrees[[i]], function(x)
				length(x$tip.label)) <= 1]
			results[[i]][problemClades] <- NA
		}
	}

	results <- lapply(results, as.data.frame)
	results <- Reduce(rbind, results)
	results
}
