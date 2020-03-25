# addTaxa
## An R package for adding missing taxa to phylogenies

The original implementations of these functions were detailed in [our American Journal of Botany paper](http://www.amjbot.org/content/102/10/1634) on the timing and tempo of the Australian Hakeinae. The version of addTaxa released in conjunction with that paper is available [here](https://github.com/eliotmiller/addTaxa/releases/tag/v0.1). The current version of the package contains an updated algorithm for adding missing taxa, and it contains new options for scaling branch lengths when taxa are added in. The new version runs more smoothly than the original version, but it has not been extensively user-tested. Create a GitHub issue and let me know if you find anything weird.

#### Why should I use addTaxa?
addTaxa allows users to add missing taxa, e.g., species, to phylogenies based on user-provided taxonomic information. It can also track additions to focal clades within the original phylogeny, permitting one to estimate the sensitivity of absolute diversification rate calculations for these clades as a function of phylogenetic uncertainty.

#### How does the new algorithm work? 
Given a data frame of two columns, which *must* be named "species" and "group", the function will take a species that is absent from the phylogeny and bind it to a randomly selected taxonomic relative. First, species are identified that are in the groupings data frame but are not present in the tree. One of these missing species (A) is selected, and a species (B) from the tree that is in that species' group is identified. If B is the only species in the tree in that group, A is bound to B at a distance below the tip determined by the branch.position argument. If the group of A+B contains additional species which are in the tree, the function then checks whether those species are monophyletic. If so, the function identifies all possible valid positions within the group to which A could be added. The root of the phylogeny and, if crown.can.move is set to FALSE, the crown node of the group, are excluded from consideration. If the species group is not monophyletic, the function bumps one node down (stemwards) in the tree towards the root and checks whether the species that descend from it all belong to the same group as B. The function continues this process until it finds a node that leads to species in multiple species groups, or hits the root. Then, all possible positions upstream (crownwards) from the deepest (oldest) acceptable are tabulated, one is randomly selected, and A is bound accordingly. This process is repeated iteratively until the tree contains all species in the groupings data frame. There are four options for how far below A will be added to whichever node is ultimately selected for binding. These are: polytomy, midpoint, uniform, and a birth-death model. "polytomy" creates a polytomy at the node to which the new tip is added (need to check but it might actually make a dichotomous branch of length 0). "midpoint", simply splits the difference between the node and its parent and binds the new tip in there. "uniform", will sample from a uniform distribution with a minimum of zero and a maximum of the full distance between the node and its parent. "bd" uses the corsim function from the [TreeSim package](https://cran.r-project.org/web/packages/TreeSim/index.html) to simulate the missing speciation events according to speciation (lambda) and extinction (mu) values calculated internally by addTaxa using [diversitree](https://cran.r-project.org/web/packages/diversitree/index.html). Additionally, the function can take a clade membership data frame, which must contain the columns "species" and "clade". When missing species are added into these clades, the data frame is updated accordingly, which facilitates calculations of the sensitivity of diversification rate later.

#### How do I install addTaxa?
addTaxa is currently available only from GitHub.

```r
library(devtools)
install_github("eliotmiller/addTaxa")
library(addTaxa)
```

#### How do I add missing species to a phylogeny?

```r
#load the turtle phylogeny from the geiger package
data(chelonia)
tree <- chelonia$phy

#some species in this tree are identified to subspecies. drop those
temp <- lapply(strsplit(tree$tip.label, "_"), length)
names(temp) <- tree$tip.label
temp <- temp[temp==2]
tree <- drop.tip(tree, setdiff(tree$tip.label, names(temp)))

#create an example groupings data frame. note that the function requires the first
#column to be named "species" and the second column to be named "group"
groupsDF <- data.frame(species=tree$tip.label)
groupsDF$group <- unlist(lapply(strsplit(tree$tip.label, "_"), "[", 1))

#the tipDropper function will randomly drop some species from the tree, ensuring
#that no genera are lost entirely. use the function to drop 100 of the 194 species
#in the phylogeny
example <- tipDropper(tree, groupsDF, 100)

#add those missing species back in, using diversitree to estimate speciation and
#extinction rate, then using those estimates to generate believable branch lengths
#with TreeSim.
newTrees <- addTaxa(tree=example, groupings=groupsDF, branch.position="bd",
  bd.type="global", no.trees=1)
```

#### How do I keep track of diversification rate sensitivity?
```r
#redefine example, but only drop 20 species from the tree this time
example <- tipDropper(tree, groupsDF, 20)

#create a dummy data frame to pay attention to diversification rates
#of two arbitrary (BUT MONOPHYLETIC) clades.
graptemysEtc <- data.frame(
  species=c(example$tip.label[grep("Graptemys", example$tip.label)],
  example$tip.label[grep("Trachemys", example$tip.label)],
  example$tip.label[grep("Malaclemys", example$tip.label)]),
  clade="graptemys")
gopherus <- data.frame(
  species=example$tip.label[grep("Gopherus", example$tip.label)],
  clade="gopherus")
cladesDF <- rbind(graptemysEtc, gopherus)

#use the divRates function to examine the sensitivity of diversification rates
#to the addTaxa function. allow the crown node of the named clades to change
#age as taxa are added, use the crown-form of the Magallon and Sanderson absolute
#diversification rate calculation, and set the extinction rate for that calculation
#to 0.1
sensitivity <- divRates(tree=example, groupings=groupsDF, branch.position="bd",
  bd.type="global", no.trees=10, clade.membership=cladesDF, crown.can.move=TRUE,
  calc.from.crown=TRUE, epsilon=0.1)

#depending on which species were removed and where those species were added back
#in, it's likely there was no variation in diversification rate in this example.
 ```

