#' Identify species to which missing taxa cannot be bound stemwards
#'
#' Given a tree and an optional data frame of clade membership, identifies species to
#' which missing taxa cannot be bound stemwards. 
#'
#' @param tree An ape-style phylogenetic tree
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
#' @param print.to.screen Logical. Will print a list of the missing taxa to screen.
#' 
#' @details It is difficult to add a missing species below the root. This would "happen" 
#' (it would actually throw an error) if the lineage that is sister to all others was a 
#' single species. Similarly, if named clades contain only a single species, then adding
#' something stemwards would render the clade polyphyletic. Finally, if named clades have
#' a species that is sister to the rest of the species in the clade, adding taxa stemwards
#' from this species would change the crown age of the clade, which may be undesirable.
#' This function identifies such species. If no named clades are provided, only potential 
#' root issues are assessed.
#'
#' @return A vector of taxa to which missing taxa cannot be added stemwards. Returns NA
# if no such taxa exist in the phylogeny.
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
#' #define a named vector of node numbers, including a root
#' temp <- c()
#' temp[1] <- length(bird.families$tip.label) + 1
#' names(temp) <- "root"
#'
#' #create a data frame of all taxa from the phylogeny, and make up clade memberships
#' #for each. note that the names on this data frame differ from "groupings" in other
#' #functions
#' dummy.frame <- data.frame(species=bird.families$tip.label, 
#' clade=c(rep("nonPasserine", 95), rep("suboscine", 9), rep("basalOscine", 13), 
#' rep("oscine", 20)))
#'
#' #note effect of crown.can.move argument. because Acanthisittidae is sister
#' #to rest of suboscines, missing taxa cannot be bound stemwards from this tip or it will
#' #change the age of the suboscine clade. 
#' dangerList(bird.families, dummy.frame, crown.can.move=TRUE, print.to.screen=FALSE)
#' dangerList(bird.families, dummy.frame, crown.can.move=FALSE, print.to.screen=FALSE)

dangerList <- function(tree, clade.membership, crown.can.move, print.to.screen)
{
	#define the nodes we need to check out. set up empty vector and ID root node. 
	possible <- c()
	possible[1] <- length(tree$tip.label) + 1
	names(possible) <- "root"
	#then if provided identify the MRCA of the named clades
	if(!missing(clade.membership))
	{
		others <- getMRCAs(tree, clade.membership)
		#getMRCAs returns a list, so unlist and append to root node
		possible <- append(possible, unlist(others))
	}

	#use the isSingleton function to check whether these nodes are trouble spots
	singletonResults <- isSingleton(tree, possible, crown.can.move)
	
	#subset the possible nodes to those that lead to singletons
	problems <- possible[singletonResults]
	
	#if any of the problems are single taxon clades, get both the node numbers of the
	#descendants of the problem nodes, and the node number of the problem tips
	if(any(problems <= length(tree$tip.label)))
	{
		singles <- problems[problems <= length(tree$tip.label)]
		descendantNodes <- c(tree$edge[tree$edge[,1] %in% problems,][,2], singles)
	}
	
	else
	{
		descendantNodes <- tree$edge[tree$edge[,1] %in% problems,][,2]
	}
	
	#subset those to just tips (still node numbers)
	descendantNodes <- descendantNodes[descendantNodes <= length(tree$tip.label)]
	
	#subset those to unique nodes (duplication can happen with tips and roots being same
	descendantNodes <- unique(descendantNodes)

	#determine who the actual taxa are, and print to screen if called for
	if(print.to.screen==TRUE & length(descendantNodes) >=1)
	{
		#figure out who those species are (assuming there are any)
		results <- tree$tip.label[descendantNodes]
		print("Will never add missing species to the stem of this/these taxon/taxa:", 
			quote=FALSE)
		print(results, quote=FALSE)
	}
	
	else if(print.to.screen==FALSE & length(descendantNodes) >=1)
	{
		#figure out who those species are (assuming there are any)
		results <- tree$tip.label[descendantNodes]
	}

	#assume that descendant nodes will evaluate to length < 1 if there are not any	
	else
	{
		results <- "NA"
	}
	
	results
}
