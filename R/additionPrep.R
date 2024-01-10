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
  #throw an error if the tree is fully dichotomous
  if(is.binary(tree) == FALSE)
  {
    stop("Input tree is not binary. Cannot currently account for this.")
  }
  
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
  bumpedBack <- 0
  
  #this will change to 1 if you respected some do.not.break statements
  respected <- 0
  
  #define blank character vector for the species you will add missing.sp as sister to
  sisterTaxa <- ""
  
  #define a blank character vector for warning messages
  warnMessages <- ""
  
  #pull the add.to column for that species and parse
  toParse <- addition.statements$add.to[addition.statements$species==missing.sp]
  
  #split it on semi-colons
  splitUp <- unlist(strsplit(toParse, "; "))
  
  #if the length is 2, check whether the second taxon/species group exists, and
  #focus on that if so. lengths of > 2 will be ignored for now. if the groups don't
  #exist, or length equals 1, focus on the first element
  if(length(splitUp)==2)
  {
    theFamily <- unique(addition.statements$family[addition.statements$family==splitUp[2]])
    theGenus <- unique(addition.statements$genus[addition.statements$genus==splitUp[2]])
    
    #if one of these exists, set tempTaxon to be that thing
    if(length(theFamily)==1 & length(theGenus)==0)
    {
      tempTaxon <- theFamily
      
      #check whether that family is in the tree yet
      allFamilySpp <- addition.statements$species[addition.statements$family==tempTaxon]
      
      #if at least one species from the family is in the tree, set taxon to be that family
      if(sum(tree$tip.label %in% allFamilySpp) > 0)
      {
        taxon <- tempTaxon
      }
      
      #otherwise set taxon to be the first element, and note you did it
      else
      {
        taxon <- splitUp[1]
        bumpedBack <- 1
      }
    }
    
    else if(length(theFamily)==0 & length(theGenus)==1)
    {
      tempTaxon <- theGenus
      
      #check whether that genus is in the tree yet
      allGenusSpp <- addition.statements$species[addition.statements$genus==tempTaxon]
      
      #if at least one species from the genus is in the tree, set taxon to be that genus
      if(sum(tree$tip.label %in% allGenusSpp) > 0)
      {
        taxon <- tempTaxon
      }
      
      #otherwise set taxon to be the first element, and note you did it
      else
      {
        taxon <- splitUp[1]
        bumpedBack <- 1
      }
    }

    else if(length(theFamily)==1 & length(theGenus)==1)
    {
      #throw an error if both of these are 1
      stop("Genus and family names cannot overlap")
    }
        
    #otherwise assume it's a species group and try there
    else
    {
      parsedSpp <- unlist(strsplit(splitUp[2], ", "))
      
      #if at least one species from the group is in the tree, set taxon to be that group
      if(sum(tree$tip.label %in% parsedSpp) > 0)
      {
        taxon <- parsedSpp
      }
      
      #otherwise set taxon to be the first element, and note you did it
      else
      {
        taxon <- splitUp[1]
        bumpedBack <- 1
      }
    }
  }
  
  #if length != 2, just set taxon to be the first element
  else
  {
    taxon <- splitUp[1]
  }

  #we need to re-figure out if the second element is a genus, family, or species group,
  #or if there was never a second element, figure it out for the first time
  theFamily <- unique(addition.statements$family[addition.statements$family==taxon])
  theGenus <- unique(addition.statements$genus[addition.statements$genus==taxon])
  
  #if there is only one of these things find all the species that belong to that
  #taxon and check whether they are monophyletic
  if(length(theFamily)==1 & length(theGenus)==0)
  {
    #this will be all the species in the whole family
    allSpp <- addition.statements$species[addition.statements$family==theFamily]
    
    #subset to those that are actually in the family
    spp <- allSpp[allSpp %in% tree$tip.label]
    
    #store a warning if spp are not monophyletic
    if(!is.monophyletic(tree, spp))
    {
      warnMessages <- "The identified addition taxon is not monophyletic"
    }
  }
  
  #do the same for a genus here
  else if(length(theFamily)==0 & length(theGenus)==1)
  {
    #this will be all the species in the whole genus
    allSpp <- addition.statements$species[addition.statements$genus==theGenus]
    
    #subset to those that are actually in the genus
    spp <- allSpp[allSpp %in% tree$tip.label]
    
    #store a warning if spp are not monophyletic
    if(!is.monophyletic(tree, spp))
    {
      warnMessages <- "The identified addition taxon is not monophyletic"
    }
  }
  
  #if both of these are length = 1, something is wrong
  else if(length(theFamily)==1 & length(theGenus)==1)
  {
    stop("Genus and family names cannot overlap")
  }
  
  #otherwise assume this is a species group. confirm all species are valid first
  else
  {
    #remember that you will have passed it a single character vector split up by commas
    #with spaces. parse this into a longer vector
    parsedSpp <- unlist(strsplit(taxon, ", "))
    
    allSpp <- addition.statements$species[addition.statements$species %in% parsedSpp]
    
    #if these are not of equal length, there's a misspelling or something
    if(length(allSpp) != length(parsedSpp))
    {
      print(missing.sp)
      stop("One of your species is misspelled")
    }
    
    #subset to species in the tree
    spp <- allSpp[allSpp %in% tree$tip.label]
  }
  
  #if the length of spp == 0, then there are no species in the tree from that group.
  #and now we are out of luck
  if(length(spp) == 0)
  {
    stop("None of those species are in the tree yet")
  }
  
  #otherwise identify the valid node(s) you can bind missing.sp to
  else
  {
    #if there is only a single species here, you need to get MRCA in a different way
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
      #potential binding sites. we are ignoring issues of non-monophyly here (we
      #threw a warning above but that was all)
      else
      {
        #tips=FALSE also returns internal nodes
        dNodes <- geiger:::.get.descendants.of.node(node=crownNode,
                                                    phy=tree, tips=FALSE)
        allNodes <- c(crownNode, dNodes)
      }
    }
  }
  
  #now dig into your do.not.break statements. i am not going to try and do anything
  #like, e.g., respect whether one taxon doesn't exist yet. simply split on semi-colons,
  #then loop through every clade and remove nodes from consideration in allNodes
  toParse <- addition.statements$do.not.break[addition.statements$species==missing.sp]
  
  #split it on semi-colons
  exclusions <- unlist(strsplit(toParse, "; "))
  
  #if there are no exclusions, skip this step
  if(length(exclusions) == 0)
  {
    TRUE
  }
  
  #otherwise loop over them
  else
  {
    #set the exclusion indicator
    respected <- 1
    
    for(i in 1:length(exclusions))
    {
      #figure out if that element is a genus, family, or species group
      theFamily <- unique(addition.statements$family[addition.statements$family==exclusions[i]])
      theGenus <- unique(addition.statements$genus[addition.statements$genus==exclusions[i]])
      
      #if there is only one of these things find all the species that belong to that
      #taxon and check whether they are monophyletic
      if(length(theFamily)==1 & length(theGenus)==0)
      {
        #this will be all the species in the whole family
        allSpp <- addition.statements$species[addition.statements$family==theFamily]
        
        #subset to those that are actually in the family
        spp <- allSpp[allSpp %in% tree$tip.label]
        
        #append a warning if spp are not monophyletic
        if(!is.monophyletic(tree, spp))
        {
          tempWarning <- "The exclusion taxon is not monophyletic"
          warnMessages <- paste(warnMessages, tempWarning, sep=". ")
        }
      }
      
      #do the same for a genus here
      else if(length(theFamily)==0 & length(theGenus)==1)
      {
        #this will be all the species in the whole genus
        allSpp <- addition.statements$species[addition.statements$genus==theGenus]
        
        #subset to those that are actually in the tree
        spp <- allSpp[allSpp %in% tree$tip.label]
        
        #append a warning if spp are not monophyletic
        if(!is.monophyletic(tree, spp))
        {
          tempWarning <- "The exclusion taxon is not monophyletic"
          warnMessages <- paste(warnMessages, tempWarning, sep=". ")
        }
      }
      
      #if both of these are length = 1, something is wrong
      else if(length(theFamily)==1 & length(theGenus)==1)
      {
        stop("Genus and family names cannot overlap")
      }
      
      #otherwise assume this is a species group. confirm all species are valid first
      else
      {
        #remember that you will have passed it a single character vector split up by commas
        #with spaces. parse this into a longer vector
        parsedSpp <- unlist(strsplit(exclusions[i], ", "))
        
        allSpp <- addition.statements$species[addition.statements$species %in% parsedSpp]
        
        #if these are not of equal length, there's a misspelling or something
        if(length(allSpp) != length(parsedSpp))
        {
          print(missing.sp)
          stop("One of your species is misspelled")
        }
        
        #subset to species in the tree
        spp <- allSpp[allSpp %in% tree$tip.label]
      }
      
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
  
  #confirm that allNodes is at least one long. the way you sample needs to change
  #here too. if it's not at least one long, there are no valid nodes and you're hosed
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
