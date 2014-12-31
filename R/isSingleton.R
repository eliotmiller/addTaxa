#' Determine whether a node leads to a tip
#'
#' Takes a vector of node numbers and determines whether it leads to a tip or to an
#' internal node (plus more).
#'
#' @param tree Phylogeny in ape-format
#' @param nodes Vector of node numbers, with an optional node named root.
#' @param crown.can.move True or false. Related to whether you want to allow missing taxa
#' to be added stemwards from the input crown data of the clade in question. See details.
#'
#' @details This function determines what descends from a node. In the context of the
#' package, these will generally be nodes that subtend named clades. It takes a vector of
#' node numbers, with an optional node named "root". If no such root is provided,
#' then it will not evaluate whether there would be problems with binding taxa to the
#' species that is sister to all others (if one exists). If a tip is provided it 
#' returns TRUE. If an internal node that leads to other internal nodes is provided, it
#' returns FALSE. This is because any such node will always have taxa bound crownwards
#' from it; missing taxa are bound to species, not nodes, so adding stemwards to a species
#' descending from this node will mean the species is still crownwards from the node in
#' question. Finally, if an internal node is provided that leads to at least one tip,
#' the function returns TRUE if crown.can.move=FALSE, else returns FALSE. 
#'
#' @return A vector of TRUE/FALSE corresponding to the rules detailed above.
#'
#' @export
#'
#' @references Eliot Miller unpublished
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
#' #use the function getMRCAs() to determine the nodes subtending these named clades
#' nodes <- getMRCAs(bird.families, dummy.frame)
#' #unlist the results and append to the root node defined above
#' nodes <- append(temp, nodes)
#'
#' #note effect of crown.can.move argument. because Acanthisittidae is sister
#' #to rest of suboscines, missing taxa cannot be bound stemwards from this tip or it will
#' #change the age of the suboscine clade. 
#' isSingleton(bird.families, nodes, crown.can.move=TRUE)
#' isSingleton(bird.families, nodes, crown.can.move=FALSE)

isSingleton <- function(tree, nodes, crown.can.move)
{
	nodes <- unlist(nodes)

	#define a vector that will be set to TRUE if a given node leads to problem tips that
	#should not be bound stemwards (i.e. you do not ever want taxa bound below this node)
	singletons <- c()
	
	for(i in 1:length(nodes))
	{
		#define the nodes that descend from the node in question
		nodeDescendants <- tree$edge[tree$edge[,1] == nodes[i],][,2]
		
		#this will evaluate to TRUE if a node that is just a tip is provided
		if(nodes[i] <= length(tree$tip.label))
		{
			singletons[i] <- TRUE
		}

		#this will evaluate to TRUE if the root leads directly to at least one tip, i.e.
		#there is at least 1 nodeDescendants with node less than or equal to the node 
		#numbers of the tips
		
		else if(names(nodes[i]) == "root" 
			& sum(nodeDescendants <= length(tree$tip.label)) >= 1)
		{
			singletons[i] <- TRUE
		}
		
		#if one of nodeDescendants is a tip AND if crown.can.move=FALSE,
		#this will evaluate to TRUE, meaning taxa will not be bound stemwards from the
		#node subtending the clade in question
		else if(sum(nodeDescendants <= length(tree$tip.label)) >= 1 
			& crown.can.move==FALSE)
		{
			singletons[i] <- TRUE
		}
		
		#if none of the above conditions are met, return FALSE
		else
		{
			singletons[i] <- FALSE
		}
	}
	singletons
}
