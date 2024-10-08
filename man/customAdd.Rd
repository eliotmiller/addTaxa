% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/customAdd.R
\name{customAdd}
\alias{customAdd}
\title{Bind species to hypothesized relatives in a phylogeny}
\usage{
customAdd(tree, addition.statements, no.trees)
}
\arguments{
\item{tree}{An ape-style phylogenetic tree.}

\item{addition.statements}{See details.}

\item{no.trees}{The number of desired final trees with all missing species from 
groupings added.}
}
\value{
A list with two elements: (1) multiPhylo object with number of trees as 
determined by no.trees, and, if clade.membership was provided, (2) a list of data 
frames summarizing the addition results.
}
\description{
Given a data frame of species and taxonomic assignments, and an accepted phylogeny with
some of those species in it, will add the missing species in next to a taxonomic 
relative. Optionally, can keep track of which named clades species are bound to.
}
\details{
All the gory stuff.
}
\references{
Mast et al. 2015. Paraphyly changes understanding of timing and tempo of 
diversification in subtribe Hakeinae (Proteaceae), a giant Australian plant radiation.
American Journal of Botany.
}
