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
#' speciation events according to speciation (lambda) and extinction (mu) values passed
#' to addTaxa. 
#' @param ini.lambda Initial speciation value for the "bd" optimization, if that option
#' of branch position is chosen. Defaults to 1.
#' @param ini.mu Initial extinction value for the "bd" optimization, if that option
#' of branch position is chosen. Defaults to 0.1.
#' @param bd.type Whether to use a 'local' or a 'global' birth-death estimation. The former
#' is slower but tends to perform better, particularly if there are shifts in diversification
#' rates in the tree. The latter runs faster.
#' @param local.cutoff If using the local form of the bd method, this argument specifies
#' at what point the local neighborhood is deemed sufficiently large to estimate a local
#' diversification rate. For example, if this argument is set to 10, a local estimate is
#' not derived until the species being added is being bound into a monophyletic clade of
#' at least 10 tips belonging to that species' group. The method seems to perform better
#' when this cutoff is set low, but smaller numbers increase run time. The cutoff must be
#' 3 or higher, as the birth-death estimation will always fail for two taxa. It can also
#' fail or timeout (see below) for larger cutoffs--when it does, it uses the midpoint
#' branch.position method and moves onto the next taxon.
#' @param timeout The amount of time to sample in search of a new branch position that
#' accords with the local birth death estimation and that falls within the age of the
#' branch to which the new species is being bound. Larger values may increase accuracy,
#' but definitely increase run time. Default is 2 seconds. 
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
#' @param calc.from.crown Whether the crown or stem method of bd.ms from geiger will be
#' used.
#' @param epsilon The extinction rate parameter for bd.ms
#' 
#' @details Adds missing taxa according to all rules and details outlined in
#' addTaxa(). During this process, it tabulates to which named clades missing
#' species are added, and calculates the diversification rate of each clade after each
#' complete tree is created. In this way, the sensitivity of these absolute
#' diversification rates to missing taxa can be assessed.
#'
#' @return A dataframe where each row corresponds to a complete tree, and each column
#' corresponds to the calculated absolute diversification rate of a named clade.
#'
#' @export
#'
#' @importFrom geiger bd.ms
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
#' #drop 20 species from the tree
#' example <- tipDropper(tree, groupsDF, 20)
#'
#' #create a dummy data frame to pay attention to diversification rates
#' #of two arbitrary clades
#' graptemysEtc <- data.frame(
#'   species=c(example$tip.label[grep("Graptemys", example$tip.label)],
#'   example$tip.label[grep("Trachemys", example$tip.label)],
#'   example$tip.label[grep("Malaclemys", example$tip.label)]),
#'   clade="graptemys")
#' gopherus <- data.frame(
#'   species=example$tip.label[grep("Gopherus", example$tip.label)],
#'   clade="gopherus")
#' cladesDF <- rbind(graptemysEtc, gopherus)
#'
#' #use the divRates function to examine the sensitivity of diversification rates
#' #to the addTaxa function
#' sensitivity <- divRates(tree=example, groupings=groupsDF, branch.position="bd",
#'   bd.type="global", no.trees=3, clade.membership=cladesDF, crown.can.move=TRUE,
#'   calc.from.crown=TRUE, epsilon=0.1)

divRates <- function(tree, groupings, branch.position, ini.lambda=1,
	ini.mu=0.1, bd.type, local.cutoff, timeout, no.trees, clade.membership,
	crown.can.move, calc.from.crown, epsilon)
{
	#generate your complete trees and updated tables of to which clade each sp belongs
	treesNtables <- addTaxa(tree=tree, groupings=groupings, 
		branch.position=branch.position, ini.lambda=ini.lambda, ini.mu=ini.mu,
		bd.type=bd.type, local.cutoff=local.cutoff, timeout=timeout,
		no.trees=no.trees, clade.membership=clade.membership, crown.can.move)
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
