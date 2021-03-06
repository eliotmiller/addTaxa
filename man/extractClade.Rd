% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractClade.R
\name{extractClade}
\alias{extractClade}
\title{Extract multiple subclades from a phylogeny based on node numbers}
\usage{
extractClade(tree, nodes, root.edge = 0)
}
\arguments{
\item{tree}{An ape-style phylogenetic tree.}

\item{nodes}{A named list or vector of node numbers that subtend the clade in question
to extract. Can extract and create single species trees with appropriate branch length.
The input vector or list must have names!}

\item{root.edge}{Default is 0, meaning that no stem is left below the node from which
the subclade is extracted. If root.edge is 1, then a root edge (stem) is also included
with the extracted clade.}
}
\value{
A named list of trees.
}
\description{
Given a tree and a vector or list of nodes, this function extracts subclades and
returns a list of smaller trees.
}
\details{
Given a named list or vector of node numbers that subtend one or more clades,
extracts those clades and returns a list of trees. The nodes provided must have names,
though the tree itself does not need to have named nodes. This is because the output
trees are provided with these names, which are used in further downstream analyses. 
If the node provided relates to a tip, then the tree that is returned contains the 
stem length from the ancestor of that species to the tip in question.
}
\examples{
#load a molecular tree up
data(bird.families)

#define a named vector of node numbers, including a root
temp <- c()
temp[1] <- length(bird.families$tip.label) + 1
names(temp) <- "root"

#create a data frame of all taxa from the phylogeny, and make up clade memberships
#for each. note that the names on this data frame differ from "groupings" in other
#functions
dummy.frame <- data.frame(species=bird.families$tip.label, 
clade=c(rep("nonPasserine", 95), rep("suboscine", 9), rep("basalOscine", 13), 
rep("oscine", 20)))

#use the function getMRCAs() to determine the nodes subtending these named clades
nodes <- getMRCAs(bird.families, dummy.frame)
#unlist the results and append to the root node defined above
nodes <- append(temp, nodes)

#use the function. note the effect of including a root edge or not. also note that 
#because non-passerines are not monophyletic, the "subtree" corresponding to that
#"clade" is the whole tree
temp1 <- extractClade(bird.families, nodes)
temp2 <- extractClade(bird.families, nodes, root.edge=1)
plot(temp1$oscine, root.edge=TRUE)
plot(temp2$oscine, root.edge=TRUE)
}
\references{
Mast et al. 2015. Paraphyly changes understanding of timing and tempo of 
diversification in subtribe Hakeinae (Proteaceae), a giant Australian plant radiation.
American Journal of Botany.
}
