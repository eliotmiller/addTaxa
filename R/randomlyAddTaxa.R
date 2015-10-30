#' Bind species to hypothesized relatives in a phylogeny
#'
#' Given a data frame of species and taxonomic assignments, and an accepted phylogeny with
#' some of those species in it, will add the missing species in next to a taxonomic 
#' relative. Optionally, can keep track of which named clades species are bound to.
#'
#' @param tree An ape-style phylogenetic tree
#' @param groupings A data frame with two columns, "species" and "group". Missing species,
#' to be added, are taken as those that do not match a value in the tip labels of tree.
#' @param from.node Whether species should be added in a "polytomy" with, "crown" (more 
#' recently diverged), "stem" (previously diverged), or "randomly" crown-wards or
#' stem-wards from the tip they are bound to.
#' @param branch.position Once the algorithm has selected to which branch a species will be
#' bound, this determines where on the branch that happens. Currently there are three
#' options. The first, "midpoint", simple breaks the branch in half and binds the species
#' there. The second, "uniform", will sample from a uniform distribution with a minimum of
#' zero and a maximum of the full branch length. The third, "bd", is currently being
#' developed and is not fully operational. For now, all it does is take the distribution
#' of branches from the full tree and scale it to the branch in question, then sample from
#' that. Thus, if there are lots of long branches on average in the phylogeny, the species
#' being bound will tend to be bound closer to the subtending node. Since this is just
#' based on the tree-wide distribution of branch lengths as opposed to the lengths in that
#' "portion" of the tree, this is not currently recommended. Note that options 2 and 3
#' (particularly "bd") mean that species can conceivably be added at the maximum distance
#' possible from the tip, i.e. at the node subtending them. This would create branch
#' lengths of zero. To avoid this, the argument "optional.offset" can be set to some small
#' value, e.g. 1e-9.
#' @param optional.offset Value to avoid polytomies (see above).
#' @param no.trees The number of desired final trees with all missing species from 
#' groupings added.
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
#' @param file.name Optional name for output RData file. If not provided then the object
#' will be stored in memory instead of saved to workspace. Either way, it is calculated in
#' memory. The function could probably be sped up by having this save directly to file,
#' appending each tree instead of first saving them all to memory. The saved object can be
#' loaded with load(), and will be imported with the name "taxaAddResults".
#' 
#' @details Given a data frame of two columns, "species" and "group", will take a species
#' that is absent from the phylogeny and bind it to a randomly selected taxonomic 
#' relative. Note that the current implementation of the function iteratively builds up
#' the tree, meaning that a species being bound into the tree can be bound to a
#' species that was in the input phylogeny, or to a species that was added during a
#' previous iteration of the function. Previously, added species could only be bound to
#' species that were in the input phylogeny. This new feature slows the function down,
#' but it results in more realistic-looking, less ladderized trees. Additionally, previous
#' iterations of this function also worked from top to bottom of the list of missing taxa.
#' This potentially results in a biased exploration of possible tree space. The function
#' now samples missing taxa and updates the missing list accordingly to allow a less
#' biased exploration of possible complete trees. Four distinct methods
#' of adding new taxa are possible. With "polytomy", the new species is simply assigned to
#' a tip. Both species end up with branch lengths to their most recent common ancestor of
#' 0. With "crown", the new species is bound in at half the distance between the species
#' it is being bound to and that species' original parent node. With "stem", the new
#' species is bound to the parent node of the species it was selected to be bound to. The
#' new species is bound in at half the distance between the parent node and grandparent
#' node. With "randomly", each new species is randomly bound either crown-wards or
#' stem-wards, following the descriptions above. Note that even if "stem" or "randomly" 
#' are chosen, if a species is to be bound to the sister species to all others (the most
#' "basal" species in the phylogeny), it will automatically be bound stem-wards. With all
#' options, if the input tree is ultrametric, the output trees should remain so. 
#' Currently, no effort is made to ensure that the taxonomic groups of the missing species
#' are actually to be found in the species in the input tree.
#'
#' @return A list with two elements: (1) multiPhylo object with number of trees as 
#' determined by no.trees, and, if clade.membership was provided, (2) a list of data 
#' frames summarizing which named clade each species in the complete phylogeny belongs to
#' (if any).
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
#' crown <- randomlyAddTaxa(tree=bird.families, groupings=groupsDF, from.node="crown", 
#' 	no.trees=10, clade.membership=cladesDF, crown.can.move=TRUE, print.to.screen=FALSE)
#' stem <- randomlyAddTaxa(tree=bird.families, groupings=groupsDF, from.node="stem", 
#' 	no.trees=10, clade.membership=cladesDF, crown.can.move=TRUE, print.to.screen=FALSE)
#' polytomy <- randomlyAddTaxa(tree=bird.families, groupings=groupsDF,
#'	from.node="polytomy", no.trees=10, clade.membership=cladesDF, crown.can.move=TRUE, 
#'	print.to.screen=FALSE)

randomlyAddTaxa <- function(tree, groupings, from.node="randomly",
	branch.position="midpoint", optional.offset=0, no.trees, clade.membership,
	crown.can.move=TRUE, print.to.screen=FALSE, file.name)
{
	#add a line to throw an error if from.node is not properly specified
	if(from.node != "randomly" & from.node != "crown" & from.node != "stem" 
		& from.node != "polytomy")
	{
		stop("from.node must be set to one of 'randomly', 'crown', 'stem', or 'polytomy'")
	}

	#not sure what would happen if provided with a polytomous tree, so throw an error here
	if(is.binary.tree(tree) == FALSE)
	{
		stop("Input tree is not binary. Cannot currently account for this.")
	}

	#use identify missing to figure out which taxa are not in the input tree
	possGroupings <- identifyMissing(tree, groupings, print.to.screen)

	#add a check to ensure that all missing species are included in a taxonomic group that
	#is included in the input tree
	inTree <- groupings[!(groupings$species %in% possGroupings$species),]
	inTree <- unique(inTree$group)

	if(length(setdiff(unique(possGroupings$group), inTree)) > 0)
	{
		stop("All missing species must be part of a taxonomic group in the input tree")
	}

	#use dangerList() to figure out who cannot be bound to stemwards
	problemTaxa <- dangerList(tree, clade.membership, crown.can.move, print.to.screen)

	#set to character clade memberships, if provided
	if(!missing(clade.membership))
	{
		clade.membership$species <- as.character(clade.membership$species)
		clade.membership$clade <- as.character(clade.membership$clade)
	}

	#set up an empty list to save the complete trees into
	random.trees <- list()
	
	#set up an empty list to save the optionally complete clade.membership data frames
	clade.memberships <- list()
	
	for(i in 1:no.trees)
	{
		#shuffle the possGroupings table so that each new tree more efficiently explores
		#possible tree space
		possGroupings <- possGroupings[order(sample(possGroupings$species)),]

		#define a new tree. this is so that it starts with the input tree not the final
		#tree at the end of each for loop
		new.tree <- tree
		
		#define a temp clade membership table so it starts with the input at end of each.
		#also set up an empty table to keep track of where species are bound in according
		#to the optional clade.membership input
		if(!missing(clade.membership))
		{
			clade.membership.temp <- clade.membership
			tabulations <- matrix(nrow=dim(possGroupings)[1], ncol=2)
			tabulations <- as.data.frame(tabulations)
			names(tabulations) <- c("species", "clade")

		}
		else
		{
			clade.memberships[[i]] <- "No clade membership input provided"
		}

		for(j in 1:dim(possGroupings)[1])
		{
			#if from.node was set to "randomly", define whether the current species to be
			#added will be added above or below the node
			if(from.node=="randomly")
			{
				chooseFrom <- c("crown","stem")
				from.nodeUsed <- sample(chooseFrom,1)
			}
			
			else
			{
				from.nodeUsed <- from.node
			}

			#subset the groupings to those we have phylogenetic information for.
			#note that because new.tree gets updated after each iteration,
			#this gets updated after each species is added in, so an unsequenced
			#species can get bound to what was also an unsequenced sp on a previous loop
			grouped <- groupings[groupings$species %in% new.tree$tip.label, ]

			#subset the species that are grouped to those whose group matches
			#the species you are adding in. 
			bindingToList <- grouped$species[grouped$group == possGroupings$group[j]]

			#check whether there is more than one species in the list (the species group).
			#if so, set big enough to 1. later down the line, this means we will be ok
			#with adding the missing species stemwards. if there is a single species only
			#set big enough to 0. in this case, adding the missing species stemward would
			#render the species group polyphyletic. while it is entirely possible that a
			#species group could still be polyphyletic and have multiple taxa in input
			#phylogeny, our adding of missing taxa will not be responsible for making the
			#species group polyphyletic (in other words we do not preclude making the 
			#species group "more" polyphyletic)
			
			if(length(bindingToList) > 1)
			{
				bigEnough <- 1
			}
			else
			{
				bigEnough <- 0
			}

			#randomly choose one of the species that is currently in the phylogeny and is
			#in the same taxonomic group as the species being added.
			#identify which "node" that randomly chosen species is
			bindingToSpecies <- sample(bindingToList, 1)
	
			#add a line to automatically switch from.node back to "crown" if the species
			#being bound to is on the danger list
			if(bindingToSpecies %in% problemTaxa 
				& (from.node=="stem" | from.node=="randomly"))
			{
				from.nodeUsed <- "crown"
			}

			#this identifies the node of the bindingToSpecies
			bindingTo <- which(new.tree$tip.label==bindingToSpecies)
			
			#add an internal if statement here that will keep track of, if provided, the
			#clade membership of where missing species are added
			if(!missing(clade.membership))
			{
				#identify the clade, if there is one, of the species that will be bound to
				boundClade <- clade.membership.temp$clade[clade.membership.temp$species==
					bindingToSpecies]
				
				#if bindingToSpecies is not in clade membership temp then boundClade
				#will have a length of 0
				if(length(boundClade)==1)
				{
					#set the first column of the given row of the tabulations table to the 
					#species being bound
					tabulations[j,1] <- possGroupings$species[j]
					
					#set the second column to the clade identity of the species to which
					#it is being bound. 
					tabulations[j,2] <- boundClade
					
					#if you are tabulating things this way, then you need to update the
					#clade membership input so that we can identify clades at the end of
					#the process
					clade.membership.temp <- rbind(clade.membership.temp, tabulations[j,])
				}
			}

			#identify the node that subtends the species you selected
			parent <- new.tree$edge[,1][new.tree$edge[,2]==bindingTo]

			#identify the node that subtends the parent node
			grandparent <- new.tree$edge[,1][new.tree$edge[,2]==parent]

			#set up a temporary matrix and give it row names. this allows you to pull out
			#the index of the edge in question, and subset the edge.lengths based on that
			#index, to get needed branch lengths later
			tempMatrix <- new.tree$edge
			rownames(tempMatrix) <- 1:dim(tempMatrix)[1]

			#use the matrix to get the index needed below
			parentIndex <- rownames(tempMatrix)[tempMatrix[,1]==parent 
				& tempMatrix[,2]==bindingTo]
			grandparentIndex <- rownames(tempMatrix)[tempMatrix[,1]==grandparent 
				& tempMatrix[,2]==parent]
			parentIndex <- as.numeric(parentIndex)
			grandparentIndex <- as.numeric(grandparentIndex)
			
			#find the edge length between the tip species you randomly selected and
			#its parent node, and the edge length between the parent and grandparent node
			parentDistance <- new.tree$edge.length[parentIndex]
			grandparentDistance <- new.tree$edge.length[grandparentIndex]

			#use bind.tip to add in the tip. if only 1 taxon is present in the phylogeny
			#for the species group in question, the new species has to be bound directly
			#to this to maintain monophyly			
			if(bigEnough == 0)
			{
				if(branch.position=="midpoint")
				{
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j], 
						where=bindingTo, position=parentDistance/2)
				}
				else if(branch.position=="bd")
				{
					#scale the distribution of original branching times to the age from
					#the parent node to the present. to avoid polytomies, it is useful to
					#add a very small constant here. sample from the new distribution
					newPositions <- branchScaler(input.vector=tree$edge.length,
						max.age=parentDistance)
					newPosition <- sample(newPositions, 1)
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j], 
						where=bindingTo, position=parentDistance-newPosition
						 + optional.offset)
				}
				else if(branch.position=="uniform")
				{
					#scale the distribution of original branching times to the age from
					#the parent node to the present. to avoid polytomies, it is useful to
					#add a very small constant here. sample from the new distribution
					newPositions <- branchScaler(input.vector=tree$edge.length,
						max.age=parentDistance)
					newPosition <- sample(newPositions, 1)
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j], 
						where=bindingTo,
						position=parentDistance-runif(n=1, min=0, max=parentDistance)
							+ optional.offset)
				}
				else
				{
					stop("Ensure arguments 'branch.position' and 'from.node' are properly specified")
				}				
			}

			#bind new species directly to randomly selected species if "crown" is selected
			#give the new species a branch length of half the original distance of the
			#selected species and its parent node
			else if(bigEnough == 1 & from.nodeUsed == "crown")
			{
				if(branch.position=="midpoint")
				{
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j], 
						where=bindingTo, position=parentDistance/2)
				}
				else if(branch.position=="bd")
				{
					#scale the distribution of original branching times to the age from
					#the parent node to the present, then sample from it
					newPositions <- branchScaler(input.vector=tree$edge.length,
						max.age=parentDistance)
					newPosition <- sample(newPositions, 1)
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j], 
						where=bindingTo, position=parentDistance-newPosition
							+ optional.offset)
				}
				else if(branch.position=="uniform")
				{
					#scale the distribution of original branching times to the age from
					#the parent node to the present. to avoid polytomies, it is useful to
					#add a very small constant here. sample from the new distribution
					newPositions <- branchScaler(input.vector=tree$edge.length,
						max.age=parentDistance)
					newPosition <- sample(newPositions, 1)
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j], 
						where=bindingTo, 
						position=parentDistance-runif(n=1, min=0, max=parentDistance)
							+ optional.offset)
				}
				else
				{
					stop("Ensure arguments 'branch.position' and 'from.node' are properly specified")
				}				
			}
			
			#if more than one species is present in group, and if "polytomy" is selected
			#bind the new species to the parent node and do not assign it any position			
			else if(bigEnough == 1 & from.nodeUsed == "polytomy")
			{
				new.tree <- bind.tip(tree=new.tree, tip.label=possGroupings$species[j],
					where=bindingTo)
			}

			#bind new species to parent node if "stem" is selected. give it an additional
			#branch length of half the distance to the grand parent node (phytools will
			#go ahead and make it ultrametric) if "midpoint" is selected, else let "bd"
			#find it
			else if(bigEnough == 1 & from.nodeUsed == "stem")
			{
				if(branch.position=="midpoint")
				{
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j],
						where=parent, position=grandparentDistance/2)
					workAround <- write.tree(new.tree)
					new.tree <- read.newick(text=workAround)
				}
				else if(branch.position=="bd")
				{
					#scale the distribution of original branching times to the age between
					#the grandparent and parent nodes, then sample from it
					newPositions <- branchScaler(input.vector=tree$edge.length,
						max.age=grandparentDistance)
					newPosition <- sample(newPositions, 1)
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j], 
						where=parent, position=grandparentDistance-newPosition
							+ optional.offset)
					workAround <- write.tree(new.tree)
					new.tree <- read.newick(text=workAround)
				}
				else if(branch.position=="uniform")
				{
					#scale the distribution of original branching times to the age from
					#the parent node to the present. to avoid polytomies, it is useful to
					#add a very small constant here. sample from the new distribution
					newPositions <- branchScaler(input.vector=tree$edge.length,
						max.age=parentDistance)
					newPosition <- sample(newPositions, 1)
					new.tree <- bind.tip(tree=new.tree,
						tip.label=possGroupings$species[j], 
						where=bindingTo,
						position=parentDistance-runif(n=1, min=0, max=parentDistance)
							+ optional.offset)
				}
				else
				{
					stop("Ensure arguments 'branch.position' and 'from.node' are properly specified")
				}				
			}
		}		
		#add the last version of new tree in as a new element in the list of random, 
		#final trees
		random.trees[[i]] <- new.tree
		
		if(!missing(clade.membership))
		{
			clade.memberships[[i]] <- clade.membership.temp
		}
	}

	class(random.trees) <- "multiPhylo"
	taxaAddResults <- list("trees"=random.trees, "clade.memberships"=clade.memberships)

	if(missing(file.name))
	{
		return(taxaAddResults)
	}
	else
	{
		save(taxaAddResults, file=file.name)
		print("Results saved to working directory")
	}
}
