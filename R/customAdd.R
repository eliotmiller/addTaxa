#' Bind species to hypothesized relatives in a phylogeny
#'
#' Given a data frame of species and taxonomic assignments, and an accepted phylogeny with
#' some of those species in it, will add the missing species in next to a taxonomic 
#' relative. Optionally, can keep track of which named clades species are bound to.
#'
#' @param tree An ape-style phylogenetic tree.
#' @param addition.statements See details.
#' @param no.trees The number of desired final trees with all missing species from 
#' groupings added.
#' 
#' @details All the gory stuff.
#'
#' @return A list with two elements: (1) multiPhylo object with number of trees as 
#' determined by no.trees, and, if clade.membership was provided, (2) a list of data 
#' frames summarizing the addition results.
#'
#' @export
#'
#' @importFrom stats runif
#' @importFrom phytools bind.tip
#' @importFrom ape is.binary.tree is.monophyletic getMRCA branching.times
#' @importFrom paleotree timeSliceTree
#' @importFrom R.utils withTimeout
#'
#' @references Mast et al. 2015. Paraphyly changes understanding of timing and tempo of 
#' diversification in subtribe Hakeinae (Proteaceae), a giant Australian plant radiation.
#' American Journal of Botany.

customAdd <- function(tree, addition.statements, no.trees)
{
	#set up a blank list and set aside the orig tree and DF to reload below
	trees <- list()
	addition.results <- list()
	origTree <- tree

	#use identify missing to figure out which taxa are not in the input tree
	toAdd <- addition.statements[!(addition.statements$species %in% tree$tip.label), ]

	#begin outer loop where you aggregate complete trees
	for(i in 1:no.trees)
	{
	  #shuffle the toAdd table so that each new tree more efficiently explores
	  #possible tree space. i feel like this needs to be nested within loop i in the 
	  #addTaxa function. should check that sometime.
	  toAdd <- toAdd[order(sample(toAdd$species)),]
	  
	  #set up a results file to save into
	  results <- data.frame(missing.sp=toAdd$species, bind.node=0,
	                        bumped.back=0, respected.exclusion=0, sister.taxa="",
	                        warnings="")
	  
	  #begin inner loop where you build up individual trees
		for(j in 1:dim(toAdd)[1])
		{
		  #run the additionPrep function
		  print(toAdd$species[j])
		  prepResults <- additionPrep(tree=tree, addition.statements=addition.statements,
		                              missing.sp=toAdd$species[j])
		  
			#we are doing midpoint branch lengths for all addition. first 
		  #identify the node that subtends the selected node to bind to
			parent <- tree$edge[,1][tree$edge[,2]==prepResults$bind.node]

			#find the distance between these two nodes. first
			#set up a temporary matrix and give it row names. this allows you to pull out
			#the index of the edge in question, and subset the edge.lengths based on that
			#index, to get needed branch lengths later
			tempMatrix <- tree$edge
			rownames(tempMatrix) <- 1:dim(tempMatrix)[1]

			#use the matrix to get the index needed below
			nodeIndex <- rownames(tempMatrix)[tempMatrix[,1]==parent 
				& tempMatrix[,2]==prepResults$bind.node]
			nodeIndex <- as.numeric(nodeIndex)

			#define the distance to bind as half distance to the parent node
			bindDist <- tree$edge.length[nodeIndex]/2

			#bind the new tip in
			tree <- phytools::bind.tip(tree=tree, tip.label=toAdd[j,"species"],
			                           where=prepResults$bind.node, position=bindDist)
			
      #log the info in the results
			results$bind.node[j] <- prepResults$bind.node
			results$bumped.back[j] <- prepResults$bumped.back
			results$respected.exclusion[j] <- prepResults$respected.exclusion
			results$sister.taxa[j] <- prepResults$respected.exclusion
			results$warnings[j] <- prepResults$warnings
		}

		#save the complete tree as an element in the list object
		trees[[i]] <- tree
		
		#save the results out too
		addition.results[[i]] <- results

		#set tree back to original
		tree <- origTree
	}
	
	#set the class of trees to multiPhylo, bind with clade.tables and return
	class(trees) <- "multiPhylo"
	finalResults <- list(trees=trees, addition.results=addition.results)
	finalResults
}
