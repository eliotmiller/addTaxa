#' Bind species to hypothesized relatives in a phylogeny
#'
#' Given a data frame of species and taxonomic assignments, and a phylogeny with
#' some of those species in it, will determine the allowed positions for a
#' particular missing taxon.
#'
#' @param tree An ape-style phylogenetic tree.
#' @param addition.statements A data frame with some specific formatting.
#' @param missing.sp The name of the missing species.
#' 
#' @details The add.to column should list a unique genus, a unique family, or a species
#' group, where each entry in the group is a valid taxon separated by a comma. The 
#' add.to column will very infrequently have a semi-colon in it. If it does,
#' what it means is that we should attempt to sub out each taxon for a clade, but that
#' if that clade doesn't exist yet, we should remove it for now. Unless one of the taxa
#' separated by a semi-colon is a species group, in which case we should add to what exists
#' of that group, if anything. This could happen if, for example,
#' two species are missing from a genus, and we know that the species are sister, 
#' but we don't know where within the genus they go. Importantly, the narrower taxonomic
#' scope taxon should always be listed second here, e.g. if you want the final arrangement
#' for the three species in GenusA to be ((sp1,sp2),sp3), and only sp3 is in there at the
#' moment, then both addition statements for sp1 and sp2 should be written as
#' GenusA; GenusA_sp1, GenusA_sp2. The algorithm will first check whether either species
#' in that group is present. If not, it will bind the taxon in to something within GenusA,
#' which will be sp3. If so, it will bind the taxon in to whichever of sp1 or sp2 is present.
#'
#' @return A vector of valid nodes to which to bind missing.sp. Something about warnings
#' and errors too.
#'
#' @export
#'
#' @importFrom phytools bind.tip
#' @importFrom ape is.binary.tree is.monophyletic getMRCA branching.times
#'
#' @references Mast et al. 2015. Paraphyly changes understanding of timing and tempo of 
#' diversification in subtribe Hakeinae (Proteaceae), a giant Australian plant radiation.
#' American Journal of Botany.

additionPrep <- function(tree, addition.statements, missing.sp)
{
  #commenting out for now. not positive this is needed
  #throw an error if the tree is not fully dichotomous
  #if(is.binary(tree) == FALSE)
  #{
  #  stop("Input tree is not binary. Cannot currently account for this.")
  #}
  
  #the function definitely fails with a weird error if the tree is not
  #ultrametric. check and stop if so
  if(is.ultrametric(tree) == FALSE)
  {
    stop("Input tree is not ultrametric.")
  }
  
  #confirm that all the relevant columns are in the addition.statements
  allCols <- c("species","family","genus","add.to","do.not.break")
  
  if(all(allCols %in% names(addition.statements)) != TRUE)
  {
    stop("You are either missing columns or one is misspelled")
  }
  
  #identify the root node here, for later use
  rootNode <- length(tree$tip.label) + 1
  
  #set up some variables to use for messaging reports.
  #this will change to 1 if you were unable to use the second element in add.to
  bumpedBack <- -1
  
  #this will change to 1 if you respected some do.not.break statements
  respected <- 0
  
  #define blank character vector for the species you will add missing.sp as sister to
  sisterTaxa <- ""
  
  #define a blank character vector for warning messages
  warnMessages <- ""
  
  #pull the add.to column for that species and parse
  toParse <- addition.statements$add.to[addition.statements$species==missing.sp]
  
  #split it on semi-colons
  splitUp <- as.list(strsplit(toParse, "; ")[[1]])
  
  # convert all taxon statements to species groups.
  for(i in 1:length(splitUp))
  {
    # see what matches to genus or to family
    theFamily <- unique(addition.statements$family[addition.statements$family==splitUp[[i]]])
    theGenus <- unique(addition.statements$genus[addition.statements$genus==splitUp[[i]]])
    
    # if splitUp[[i]] is a family, sub in those species to the statement
    if(length(theFamily)==1 & length(theGenus)==0)
    {
      tempTaxon <- theFamily
      allFamilySpp <- addition.statements$species[addition.statements$family==tempTaxon]
      splitUp[[i]] <- allFamilySpp
    }
    
    # if it's a genus, handle accordingly
    else if(length(theFamily)==0 & length(theGenus)==1)
    {
      tempTaxon <- theGenus
      allGenusSpp <- addition.statements$species[addition.statements$genus==tempTaxon]
      splitUp[[i]] <- allGenusSpp
    }
    
    # haven't encountered this but just in case
    else if(length(theFamily)==1 & length(theGenus)==1)
    {
      #throw an error if both of these are 1
      stop("Genus and family names cannot overlap")
    }
    
    # otherwise assume it already is a species group and handle accordingly
    else
    {
      parsedSpp <- unlist(strsplit(splitUp[[i]], ", "))
      
      # add a check to ensure everything is correctly spelled if species are being added directly
      checks <- addition.statements$species[addition.statements$species %in% parsedSpp]
      
      # if these are not of equal length, there's a misspelling or something
      if(length(checks) != length(parsedSpp))
      {
        print(missing.sp)
        stop("One of your species is misspelled")
      }
      
      # set it into the right place
      splitUp[[i]] <- parsedSpp
    }
  }
  
  # now figure out which statement we'll use. start at end of statements, which
  # we assume is most specific/precise, and work from there. reverse these and
  # start at the end
  targets <- rev(splitUp)
  for(i in 1:length(targets))
  {
    # figure out if any of the species is already in the tree
    spp <- tree$tip.label[tree$tip.label %in% targets[[i]]]
    
    # log one for each iteration you bump up. you start at -1, so just start now
    bumpedBack <- bumpedBack + 1
    
    # now break out of the loop if spp contains anything
    if(length(spp) > 0)
    {
      break
    }
  }
  
  # if spp doesn't exist yet (it would have been created if any of splitUp
  # contained a taxon in the tree), throw an error
  if(!exists("spp"))
  {
    stop("None of those species are in the tree yet")
  }
  
  #store a warning if spp are not monophyletic
  if(!is.monophyletic(tree, spp))
  {
    warnMessages <- "The identified addition taxon is not monophyletic"
  }

  # now process your addition statement. if there is only a single species here,
  # you need to get MRCA in a different way
  if(length(spp)==1)
  {
    allNodes <- which(tree$tip.label==spp)
  }
  
  #grab the MRCA otherwise
  else
  {
    crownNode <- getMRCA(tree, spp)
    
    #confirm the crown node isn't the root. that would be bad
    if(crownNode == rootNode)
    {
      stop("Your crown node is your root node")
    }
    
    #otherwise you are going to consider the MRCA plus all descendant nodes as
    #potential binding sites. we are ignoring issues of existing non-monophyly here (we
    #threw a warning above but that was all)
    else
    {
      #tips=FALSE also returns internal nodes
      dNodes <- geiger:::.get.descendants.of.node(node=crownNode,
                                                  phy=tree, tips=FALSE)
      allNodes <- c(crownNode, dNodes)
    }
  }
  
  #now dig into your do.not.break statements. i am not going to try and do anything
  #like, e.g., respect whether one taxon doesn't exist yet. simply split on semi-colons,
  #then loop through every clade and remove nodes from consideration in allNodes
  toParse <- addition.statements$do.not.break[addition.statements$species==missing.sp]
  
  #split it on semi-colons
  exclusions <- as.list(strsplit(toParse, "; ")[[1]])
  
  # if there are no exclusions, skip the next steps
  if(length(exclusions) == 0)
  {
    TRUE
  }
  
  # otherwise process the exclusions
  else
  {
    # convert all taxon statements to species groups.
    for(i in 1:length(exclusions))
    {
      # see what matches to genus or to family
      theFamily <- unique(addition.statements$family[addition.statements$family==exclusions[[i]]])
      theGenus <- unique(addition.statements$genus[addition.statements$genus==exclusions[[i]]])
      
      # if exclusions[[i]] is a family, sub in those species to the statement
      if(length(theFamily)==1 & length(theGenus)==0)
      {
        tempTaxon <- theFamily
        allFamilySpp <- addition.statements$species[addition.statements$family==tempTaxon]
        exclusions[[i]] <- allFamilySpp
      }
      
      # if it's a genus, handle accordingly
      else if(length(theFamily)==0 & length(theGenus)==1)
      {
        tempTaxon <- theGenus
        allGenusSpp <- addition.statements$species[addition.statements$genus==tempTaxon]
        exclusions[[i]] <- allGenusSpp
      }
      
      # haven't encountered this but just in case
      else if(length(theFamily)==1 & length(theGenus)==1)
      {
        #throw an error if both of these are 1
        stop("Genus and family names cannot overlap")
      }
      
      # otherwise assume it already is a species group and handle accordingly
      else
      {
        parsedSpp <- unlist(strsplit(exclusions[[i]], ", "))
  
        # add a check to ensure everything is correctly spelled if species are being added directly
        checks <- addition.statements$species[addition.statements$species %in% parsedSpp]
        
        # if these are not of equal length, there's a misspelling or something
        if(length(checks) != length(parsedSpp))
        {
          print(missing.sp)
          stop("One of your species is misspelled")
        }
        
        #set into right place
        exclusions[[i]] <- parsedSpp
      }
    }

    #set the exclusion indicator
    respected <- 1
    
    for(i in 1:length(exclusions))
    {
      # subset exclusions to those species that are actually in the tree
      spp <- tree$tip.label[tree$tip.label %in% exclusions[[i]]]
      
      #it is possible that the species you are trying not to break up haven't been added to the
      #tree yet. if length of spp == 0, bump to next exclusion
      if(length(spp)==0)
      {
        next()
      }
      
      #if there is only a single species here, you need to get MRCA in a different way than below
      else if(length(spp)==1)
      {
        excludeNodes <- which(tree$tip.label==spp)
      }
      
      #grab the MRCA otherwise
      else
      {
        crownNode <- getMRCA(tree, spp)
        
        #confirm the crown node isn't the root. that would be bad
        if(crownNode == rootNode)
        {
          stop("Your crown node is your root node. You're in trouble.")
        }
        
        #otherwise you are going to consider the MRCA plus all descendant nodes as
        #potential binding sites. we are ignoring issues of non-monophyly here (we
        #threw a warning above but that was all)
        else
        {
          #tips=FALSE also returns internal nodes
          dNodes <- geiger:::.get.descendants.of.node(node=crownNode,
                                                      phy=tree, tips=FALSE)
          
          #I don't believe we want to bind in the MRCA here, but need to check
          excludeNodes <- dNodes
        }
      }
      
      #i am only partially sure the logic here is correct, but implementing this for now.
      #the reason: if for example, you want to add
      #a species as sister to a clade, but at the moment that clade is only represented by
      #a single species, then excludeNodes will equal allNodes, and there will be no valid
      #placement for missing.sp. In my head this isn't necessary, since the species you do
      #want added to that clade should get added correctly later, but part of this depends
      #on the add statements I guess. For example, if the add statement for Genus 1 sp A was
      #add to Genus 1 but do not break monophyly of Genus 1, and only Genus 1 sp B was
      #present, this would allow sp A to be added sister to sp B. Then if the add statement
      #for Genus 1 sp C was add to Genus 1, then it could be added anywhere wrt to sp B and A.
      #the solution would be to make the add statement for sp C be spC+spA.
      if(length(allNodes)==1 & length(excludeNodes)==1)
      {
        if(allNodes==excludeNodes)
        {
          #for lack of a better solution, just set excludeNode to rootNode since we know
          #we never want to add there
          excludeNodes <- rootNode
        }
      }
      
      #redefine allNodes to exclude those in excludeNodes
      allNodes <- allNodes[!(allNodes %in% excludeNodes)]
    }
  }
  
  #confirm that allNodes is at least one long. 
  if(length(allNodes)==1)
  {
    bindTo <- allNodes
  }
  else if(length(allNodes) > 1)
  {
    bindTo <- sample(allNodes, 1)
  }
  
  else
  {
    stop("No valid positions exist for binding this taxon")
  }

  #figure out if bindTo is an internal node or a tip. set sisterTaxa accordingly
  if(bindTo >= rootNode)
  {
    sisterTaxa <- extract.clade(tree, bindTo)$tip.label
  }
  
  else
  {
    sisterTaxa <- tree$tip.label[bindTo]
  }
  
  #create a results object and return
  results <- list(bind.node=bindTo,
                  bumped.back=bumpedBack,
                  respected.exclusion=respected,
                  sister.taxa=paste(sisterTaxa, collapse=", "),
                  warnings=warnMessages)
  results
}
