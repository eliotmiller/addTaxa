#' Bind species to hypothesized relatives in a phylogeny
#'
#' Given a data frame of species and taxonomic assignments, and an accepted phylogeny with
#' some of those species in it, will add the missing species in next to a taxonomic 
#' relative. Optionally, can keep track of which named clades species are bound to.
#'
#' @param tree An ape-style phylogenetic tree.
#' @param groupings A data frame with two columns, "species" and "group". Missing species,
#' to be added, are taken as those that do not match a value in the tip labels of tree.
#' @param branch.position
#' Once the algorithm has selected to which node a species will
#' be bound, this determines where on the branch that happens. Note that this refers to
#' the position below (stemwards) from the node the tip will be bound.Currently there are
#' four options. The first, "polytomy", creates a polytomy at the node to which the new
#' tip is added (need to check but it might actually make a dichotomous branch of length 0).
#' The second, "midpoint", simply splits the difference between the node and its parent
#' and binds the new tip in there.
#' The third, "uniform", will sample from a uniform distribution with a minimum of
#' zero and a maximum of the full distance between the node and its parent. Note that
#' this means that species can conceivably be added at the maximum or minimum
#' distance possible from the tip, which would create a polytomy (or branch of length 0?).
#' I have removed checks to account and deal with this--if you need them back let me know. 
#' The fourth option is
#' "bd." This uses the corsim function from the TreeSim package to simulate the missing
#' speciation events according to speciation (lambda) and extinction (mu) values calculated
#' internally by addTaxa using diversitree.
#' @param ini.lambda Initial speciation value for the "bd" optimization, if that option
#' of branch position is chosen. Defaults to 1.
#' @param ini.mu Initial extinction value for the "bd" optimization, if that option
#' of branch position is chosen. Defaults to 0.1.
#' @param no.trees The number of desired final trees with all missing species from 
#' groupings added.
#' @param clade.membership An optional data frame with first column = "species", second
#' column = "clade". These are named, monophyletic clades to which the species in the
#' input phylogeny can belong. Not every species in the input phylogeny needs to be in the
#' clade membership data frame, but missing species added to species not in this frame
#' will not be included in the output clade membership data frame. Such species will
#' still be included in the output phylogenies. Note that these clades need to be
#' mutually exclusive. For example, clade 2 cannot be contained within clade 1.
#' @param crown.can.move Logical. If TRUE, and if missing taxa are to be added stemwards,
#' this will allow the age of the crown group to potentially shift back in time. If FALSE,
#' and if missing taxa were to be added stemwards, this will prevent taxa from being added
#' below the crown group. It will force these to be bound crownwards if the crown node
#' in the clade is selected to be bound to. Note that the argument
#' name is somewhat misleading in that the crown ages can still shift when missing taxa
#' are added. Specifically, if missing taxa are bound to what were initially
#' single-species clades in the input tree, then the
#' crown age will "shift" forward in time (the single-taxon clade did not actually have a 
#' crown age). If crown.can.move is set to FALSE, then after one taxon is added to such a
#' single-species clade, the crown age of that clade then becomes fixed and will not move.
#' 
#' @details Given a data frame of two columns, which *must* be named "species" and
#' "group", will take a species
#' that is absent from the phylogeny and bind it to a randomly selected taxonomic 
#' relative. The algorithm works as follows. First, species are identified that are in
#' the groupings data frame but are not present in the tree. The order of these missing
#' species is then randomized. One of these missing species (A) is selected, and a
#' species (B) from the tree that is in that species' group is identified. If B is the
#' only species in the tree in that group, A is bound to B at a distance below the tip
#' determined by the branch.position argument. If the group of A+B contains additional
#' species in the tree, the function then checks whether those species are monophyletic.
#' If so, the function identifies all possible valid positions within the group to which
#' A could be added. The root of the phylogeny and, if crown.can.move is set to FALSE,
#' the crown node of the group, are excluded from consideration. If the species group is
#' not monophyletic, the function bumps one node down (stemwards) in the tree towards the
#' root and checks whether the species that descend from it all belong to the same group
#' as B. The function continues this process until it finds a node that leads to species
#' in multiple species groups, or hits the root. Then, all possible positions upstream
#' (crownwards) from the deepest node encountered are tabulated, one is randomly selected,
#' and A is bound accordingly. This process is repeated iteratively until the tree contains
#' all species in the groupings data frame. There are a four options for how far below A
#' will be added to whichever node is ultimately selected for binding. These are:
#' polytomy, midpoint, uniform, and birth-death model, described in the branch.position
#' argument above. Additionally, the function can take a clade membership data frame,
#' which must contain the columns "species" and "clade". When missing species are added
#' into these clades, the data frame is updated accordingly, which facilitates
#' calculations of the sensitivity of diversification rate later.
#'
#' @return A list with two elements: (1) multiPhylo object with number of trees as 
#' determined by no.trees, and, if clade.membership was provided, (2) a list of data 
#' frames summarizing which named clade each species in the complete phylogeny belongs to
#' (if any).
#'
#' @export
#'
#' @importFrom stats runif
#' @importFrom phytools bind.tip
#' @importFrom ape is.binary.tree is.monophyletic getMRCA branching.times
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
#' #use the function to drop 100 species (there were 194 in tree)
#' example <- tipDropper(tree, groupsDF, 100)
#'
#' #add those missing species back in
#' newTrees <- addTaxa(tree=example, groupings=groupsDF, branch.position="bd",
#'   no.trees=1)

addTaxa <- function(tree, groupings, branch.position="midpoint", 
	ini.lambda=1, ini.mu=0, no.trees, clade.membership, crown.can.move=TRUE)
{
	#set up a blank list and set aside the orig tree and DF to reload below
	trees <- list()
	clade.tables <- list()
	origTree <- tree
	if(!missing(clade.membership))
	{
		origMembership <- clade.membership
	}

	#not sure what would happen if provided with a polytomous tree, so throw an 
	#error here
	if(is.binary.tree(tree) == FALSE)
	{
		stop("Input tree is not binary. Cannot currently account for this.")
	}

	#force groupings DF to character if it isn't. do same for clade.membership
	#if it is provided
	groupings[,1] <- as.character(groupings[,1])
	groupings[,2] <- as.character(groupings[,2])

	if(!missing(clade.membership))
	{
		clade.membership[,1] <- as.character(clade.membership[,1])
		clade.membership[,2] <- as.character(clade.membership[,2])

		#insert a check here to ensure all the named clades are actually monophyletic
		#split the clade.membership table on clade. if any are not monophyletic, stop
		splitUp <- split(clade.membership, clade.membership$clade)
		mono <- unlist(lapply(splitUp, function(x) is.monophyletic(tree, x$species)))
		if(!all(mono))
		{
			stop("All named clades in clade.membership table are not monophyletic")
		}

		#insert another check to make sure the clades are not nested
		if(length(clade.membership$species) != length(unique(clade.membership$species)))
		{
			stop("The clades are not mutually exclusive")
		}

		#insert another check to make sure all spp in clade.membership are in the
		#tree. you get really confusing error messages if they are not
		if(length(setdiff(clade.membership$species, tree$tip.label)) > 0)
		{
			stop("clade.membership contains species that are not in tree")
		}
	}

	#use identify missing to figure out which taxa are not in the input tree
	toAdd <- identifyMissing(tree, groupings)

	#shuffle the toAdd table so that each new tree more efficiently explores
	#possible tree space. 
	toAdd <- toAdd[order(sample(toAdd$species)),]

	#add a check to ensure that all missing species are included in a taxonomic
	#group that
	#is included in the input tree
	inTree <- groupings[!(groupings$species %in% toAdd$species),]
	inTree <- unique(inTree$group)

	if(length(setdiff(unique(toAdd$group), inTree)) > 0)
	{
		stop("All missing species must be part of a taxonomic group in the input tree")
	}

	#use the findRates function to estimate speciation and extinction rates for the
	#original tree. do not update these rates as the tree is built up (you used to
	#do this).
	if(branch.position=="bd")
	{
		rates <- findRates(origTree,
			prop.complete=length(origTree$tip.label)/dim(groupings)[1],
			ini.lambda=ini.lambda, ini.mu=ini.mu)
	}

	#begin outer loop where you aggregate complete trees
	for(i in 1:no.trees)
	{
		#begin inner loop where you build up individual trees
		for(j in 1:dim(toAdd)[1])
		{
			#if the clade.membership table is provided, and if crown.can.move is FALSE,
			#find the MRCA of all the named clades. continually update this so that as
			#species are added to named clades, these nodes change also. add the root
			#node to it, to simplify some things below
			if(!missing(clade.membership) & crown.can.move==FALSE)
			{
				problemNodes <- unlist(getMRCAs(tree, clade.membership))
				rootNode <- length(tree$tip.label) + 1
				names(rootNode) <- "root"
				problemNodes <- c(problemNodes, rootNode)
			}

			#if those conditions are not true, just set problemNodes to the root
			else
			{
				problemNodes <- length(tree$tip.label) + 1
				names(problemNodes) <- "root"
			}

			#find all tips in the species group that tip i belongs to
			relatedSpp <- groupings$species[groupings$group == toAdd$group[j]]

			#this included species that potentially are not yet in the tree. drop any
			#species that are not yet in the tree
			relatedSpp <- relatedSpp[relatedSpp %in% tree$tip.label]

			#if relatedSpp is of length 1, just bind directly to it and bounce down to
			#bottom to determine binding distance. a few important things to notice
			#here. first, getMRCAs will have returned the MRCA of any single-taxon
			#clades as the tip node of the taxon. this would mean that that node gets
			#excluded as a problem node and cannot be bound to. however, we do not
			#assess problem nodes in this step. the steps below where we do should not
			#be relevant because they shouldn't be single taxon clades. if on one
			#iteration we bind a species to a single species taxon, and if 
			#crown.can.move is set to FALSE, the new MRCA of the clade is going to
			#become the crown and THEN it can no longer move.
			if(length(relatedSpp) == 1)
			{
				#identify the node that represents that species
				bindTo <- which(tree$tip.label==relatedSpp)
			}

			#if there are > 1 relatedSpp, figure out whether it is a monophyletic group
			else if(is.monophyletic(tree, relatedSpp))
			{
				#if true, find the MRCA of the clade of related spp
				crownNode <- getMRCA(tree, relatedSpp)

				#now find all nodes that descend from the crownNode. previously, you were
				#only returning tips if branch.position was set to polytomy, but actually
				#as long as the species' group is monophyletic it's ok to add new species
				#into that group in polytomies with internal nodes. the tips=FALSE arg
				#returns tips and internal nodes; tips=TRUE just returns tips
				dNodes <- geiger:::.get.descendants.of.node(node=crownNode,
					phy=tree, tips=FALSE)
				
				#confirm that crownNode is not a part of the set of problem nodes
				#defined above. if it is, only consider nodes that descend from it.
				if(length(problemNodes[problemNodes == crownNode])==1)
				{
					allNodes <- dNodes
				}
				
				#otherwise add the crown node to consideration as a place to bind the tip
				else
				{
					#importantly, note that you do not have to figure out who the
					#parent node of crownNode is, because bind.tip's position arg
					#will, if given a positive value, bind stemward from crownNode
					allNodes <- c(crownNode, dNodes)
				}

				#REALLY REALLY IMPORTANT. sample produces really odd results if sampling
				#from a numeric vector of length 1. need to account for that with this
				#overcomplicated if else statement. then randomly sample allNodes for 1
				if(length(allNodes)==1)
				{
					bindTo <- allNodes
				}
				else
				{
					bindTo <- sample(allNodes, 1)
				}
			}

			#if relatedSpp are not monophyletic, go into a more costly procedure to
			#figure out where you can bind the tip in.
			else
			{
				#first randomly sample a species from relatedSpp
				bindingToSp <- sample(relatedSpp, 1)

				#identify the node that represents that species
				bindingTo <- which(tree$tip.label==bindingToSp)

				#identify the node that subtends the node in question. at this point
				#you shouldn't have to worry about going "below" the root, since there
				#will always be one node stemwards from a tip. after that we check to
				#make sure we don't go past the root
				parent <- tree$edge[,1][tree$edge[,2]==bindingTo]

				#go into a while loop where keep bumping down the tree until you hit
				#a node that doesn't lead to a set of monophyletic ancestors
				#according to the species groups.
				keepGoing <- TRUE

				while(keepGoing)
				{
					#check whether the parent node is part of problem nodes. if it is,
					#set keepGoing to FALSE so it will just bind to the bindingTo sp
					if(length(problemNodes[problemNodes == bindingTo])==1)
					{
						keepGoing <- FALSE
					}

					#figure out which tip nodes descend from this parent node
					toExamineNodes <- geiger:::.get.descendants.of.node(node=parent,
						phy=tree, tips=TRUE)

					#figure out which species these nodes represent
					toExamine <- tree$tip.label[toExamineNodes]

					#figure out whether all these species belong to the same group.
					#this will pull a vector of the group to which each species belongs
					groupToExamine <- groupings$group[groupings$species %in% toExamine]

					#if there is one unique element in this vector, then the node is
					#monophyletic. you do not need to check if it is the root, did so above
					if(length(unique(groupToExamine)) == 1)
					{
						#identify the set of possible nodes that descend from parent
						dNodes <- geiger:::.get.descendants.of.node(node=parent,
							phy=tree, tips=FALSE)
						allNodes <- c(parent, dNodes)

						#set bindingTo to be parent, then set parent to be the node one
						#below what it was
						bindingTo <- parent
						parent <- tree$edge[,1][tree$edge[,2]==parent]
					}

					#if the vector is not monophyletic, no need to keepGoing
					else
					{
						keepGoing <- FALSE
					}
				}

				#hopefully the previous while loop bumped the bindingTo node up to whatever
				#the most stemwards monophyletic node was encompassing the randomly sampled
				#species. if bindingTo is a species, then there are no descendants, and need
				#to bind directly to that tip. otherwise, identify the nodes that descend
				#from bindingTo and randomly sample one to actually bind to.		
				if(bindingTo <= length(tree$tip.label))
				{
					allNodes <- bindingTo
				}

				else
				{
					#identify descendants
					dNodes <- geiger:::.get.descendants.of.node(node=bindingTo,
						phy=tree, tips=FALSE)

					#at this point there should be no way that bindingTo is
					#the root, so add it to all descendant nodes for sampling from below
					allNodes <- c(bindingTo, dNodes)
				}		

				#REALLY REALLY IMPORTANT. sample produces really odd results if sampling
				#from a numeric vector of length 1. need to account for that with this
				#overcomplicated if else statement. then randomly sample allNodes for 1
				if(length(allNodes)==1)
				{
					bindTo <- allNodes
				}
				else
				{
					bindTo <- sample(allNodes, 1)
				}
			}

			#calculate position based on the chosen method
			if(branch.position=="polytomy")
			{
				bindDist <- 0
			}

			else if(branch.position=="midpoint")
			{
				#identify the node that subtends the selected node to bind to
				parent <- tree$edge[,1][tree$edge[,2]==bindTo]

				#find the distance between these two nodes. first
				#set up a temporary matrix and give it row names. this allows you to pull out
				#the index of the edge in question, and subset the edge.lengths based on that
				#index, to get needed branch lengths later
				tempMatrix <- tree$edge
				rownames(tempMatrix) <- 1:dim(tempMatrix)[1]

				#use the matrix to get the index needed below
				nodeIndex <- rownames(tempMatrix)[tempMatrix[,1]==parent 
					& tempMatrix[,2]==bindTo]
				nodeIndex <- as.numeric(nodeIndex)

				#define the distance to bind as half distance to the parent node
				bindDist <- tree$edge.length[nodeIndex]/2
			}

			else if(branch.position=="uniform")
			{
				#identify the node that subtends the selected node to bind to
				parent <- tree$edge[,1][tree$edge[,2]==bindTo]

				#find the distance between these two nodes. first
				#set up a temporary matrix and give it row names. this allows you to pull out
				#the index of the edge in question, and subset the edge.lengths based on that
				#index, to get needed branch lengths later
				tempMatrix <- tree$edge
				rownames(tempMatrix) <- 1:dim(tempMatrix)[1]

				#use the matrix to get the index needed below
				nodeIndex <- rownames(tempMatrix)[tempMatrix[,1]==parent 
					& tempMatrix[,2]==bindTo]
				nodeIndex <- as.numeric(nodeIndex)

				#define bindDist as a uniform value between 0 and max. could create polytomies
				#with code like this, in theory. used to add an offset here to account for that,
				#but seems like over thinking it. add later if necessary
				bindDist <- stats::runif(n=1, min=0, max=tree$edge.length[nodeIndex])
			}

			else if(branch.position=="bd")
			{
				#identify the node that subtends the selected node to bind to
				parent <- tree$edge[,1][tree$edge[,2]==bindTo]

				#find the ages of these nodes. if bindTo is a tip, then its age is 0. create
				#vector of branching times because you need it either way for age of parent
				ages <- branching.times(tree)

				if(bindTo <= length(tree$tip.label))
				{
					bindToAge <- 0
				}

				else
				{
					bindToAge <- ages[names(ages)==bindTo]
				}

				parentAge <- ages[names(ages)==parent]

				#calculate the age of the missing speciation event. this can very occasionally
				#throw errors for reasons that I don't understand. wrap it up in a tryCatch
				#and use the midpoint method if it throws an error. i think actually it might
				#return as 'numeric(0)', so try adding that logical below.				
				missingAge <- try(bdScaler(tree=tree,
					lambda=rates["lambda"], mu=rates["mu"],
					min.age=bindToAge, max.age=parentAge), silent=TRUE)

				if(class(missingAge)=="try-error" | length(missingAge) == 0)
				{
					#identify the node that subtends the selected node to bind to
					parent <- tree$edge[,1][tree$edge[,2]==bindTo]

					#find the distance between these two nodes. first
					#set up a temporary matrix and give it row names. this allows you to pull out
					#the index of the edge in question, and subset the edge.lengths based on that
					#index, to get needed branch lengths later
					tempMatrix <- tree$edge
					rownames(tempMatrix) <- 1:dim(tempMatrix)[1]

					#use the matrix to get the index needed below
					nodeIndex <- rownames(tempMatrix)[tempMatrix[,1]==parent 
						& tempMatrix[,2]==bindTo]
					nodeIndex <- as.numeric(nodeIndex)

					#define the distance to bind as half distance to the parent node
					bindDist <- tree$edge.length[nodeIndex]/2

					warning("used 'midpoint' branch.position argument for addition of a taxon")
				}

				#if missingAge is not a try-error, it should be class numeric. find
				#the distance below bindTo to bind tip in.
				else
				{
					bindDist <- missingAge-bindToAge

					#add a check here that if bindDist is below the parent node
					#it sets the age to the parent node. you would think this shouldn't be
					#possible, but TreeSim seems willing to return some weird things, so also
					#include another check that if bindDist is a negative number, it goes to 0.
					if(bindDist > parentAge-bindToAge)
					{
						bindDist <- parentAge-bindToAge
					}

					if(bindDist < 0)
					{
						bindDist <- 0
					}
				}
			}

			else
			{
				stop("branch.position must be set to one of 'polytomy', 'midpoint', 'uniform', or 'bd'")
			}

			#if clade membership table was provided, figure out which clade the bindTo
			#node belongs to, if any, and update the table accordingly
			if(!missing(clade.membership))
			{
				#this seems like it could be costly, but I can't think of a better
				#way to do it. for every
				#named clade, find the MRCA, then find all nodes that descend from
				#each of those. then figure out if bindTo is represented in any of
				#those sets. find the MRCAs of all named clades
				mrcas <- unlist(getMRCAs(tree, clade.membership))

				#identify all descendants. notice that you have getMRCAs
				#programmed so that it returns the tip of single taxon clades.
				#however, there is no such thing as the descendants of a tip.
				#deal with that below. one good thing is that if the clades
				#are monophyletic and non-overlapping, which you checked for at start
				#there is no way two clades can share a descendant.
				allDesc <- lapply(mrcas, function(x)
					geiger:::.get.descendants.of.node(node=x, phy=tree, tips=FALSE))

				#if any of these named clades return null for descendants, it
				#should mean that they are a single-taxon clade. replace those
				#descendants with the single tips
				allDesc[unlist(lapply(allDesc, is.null))] <- mrcas[unlist(lapply(allDesc, is.null))]

				#subset those descendants to any (should only be one) that match
				anyMatches <- lapply(allDesc, function(x) x[x==bindTo])

				#subset the names of the clades to the hopefully single node that 
				#represents the crown node of a clade and is considered an ancestor
				#of bindTo
				clade <- names(unlist(anyMatches))

				#if the length of clade is 0, then the node we are binding the tip to
				#was not in a named clade. in that case, do not bind it into the DF,
				#as do not need to keep track of it
				if(length(clade) == 1)
				{
					toBind <- data.frame(species=toAdd[j,1], clade=clade)
					clade.membership <- rbind(clade.membership, toBind)
				}
			}

			#finally, bind the new tip in
			tree <- phytools::bind.tip(tree=tree, tip.label=toAdd[j,1], where=bindTo,
				position=bindDist)
		}

		#if branch.position was bd, and the tree has polytomies, throw a warning
		if(branch.position=="bd" & !ape::is.binary.tree(tree))
		{
			warning("used multi2di to calculate diversification rates; polytomies present in tree")
		}

		#save the complete tree as an element in the list object. do same for
		#clade membership if it exists
		trees[[i]] <- tree
	
		if(!missing(clade.membership))
		{
			clade.tables[[i]] <- clade.membership
		}
		
		else
		{
			clade.tables[[i]] <- "clade.membership table not provided"
		}
		
		#set tree back to original tree and original clade.membership
		tree <- origTree

		if(!missing(clade.membership))
		{
			clade.membership <- origMembership
		}
	}
	
	#set the class of trees to multiPhylo, bind with clade.tables and return
	class(trees) <- "multiPhylo"
	results <- list(trees=trees, clade.tables=clade.tables)
	results
}
