#Phyloscanning for Paleo-eFLOWER (Hervé Sauquet, Apr 2017, last updated Apr 2020) modified to apply to any nexus matrix and tree. 

#Clear buffer of all objects
	rm(list=ls())

#Load packages
	library(ape) #framework for handling phylogenetic trees
	library(phangorn) #used for MP ancestral state reconstruction (ancestral.pars function)
	library(geiger) #required for function tips
	library(plyr) #required for calculating step frequencies using the count function
	library(Claddis)
	
plotDoylegram = function(treename, matrixname, fossiltaxa,){	
#Load custom functions
 sapply(list.files(pattern="[.]R$", path="./Functions/", full.names=TRUE), source);
#Input tree
  tree <- read.nexus(treename)
	tree$edge.length <- rep(1,length(tree$edge[,1]))
#Identify fossil taxa
#	fossiltaxa <- 'Mesodescolea'
	nrfossils <- length(fossiltaxa)

#Read character matrix
	charmatrix <- read_nexus_matrix(matrixname, equalize_weights = FALSE)
	pmatrix <- charmatrix$matrix_1$matrix
	pmatrix[is.na(pmatrix)] <- '?' 
	
#Prepare matrix in the strange format required by phangorn
	#Build a custom so-called 'contrast' matrix to define polymorphic states
		#! this may see, very complicated, but there is no other way: phangorn does not seem to have a native system to read polymorphic data for the custom "USER" type (e.g., morphological), unlike corHMM
		#! the code below assumed that polymorphic data have been provided as requested by corHMM (e.g., 0/1 = "0&1")
		#! the code below will only work for characters with symbols in "0123456789" (more work needed for allowing non-numeric symbols)
		#! the alternative is to transform all polymorphic data into missing data, but this is not the exact same inference; using contrasts is much better and replicates exactly what Mesquite does
		nrstates <- 	max(charmatrix$matrix_1$maximum_values) ++ 1
		contrast <- diag(nrstates)
		levels <- levels(factor(pmatrix))
		symbols <- levels(factor(as.numeric(levels)))
		dimnames(contrast) <- list(symbols, symbols)
		for (i in 1:length(levels)) {
			if (is.na(as.numeric(levels)[i])==TRUE) { #NA values introduced by this conversion correspond to polymorphic and missing data codes
				if (levels[i]=="?") {
					newcontrast <- rep(1,nrstates)
					contrast <- rbind(contrast, newcontrast)
				} else {
					atomicsymbols <- unlist(strsplit(levels[i], "&"))
					newcontrast <- rep(0,nrstates)
					contrast <- rbind(contrast, newcontrast)
					contrast[nrow(contrast), atomicsymbols] <- rep(1, length(atomicsymbols))
				}
				rownames(contrast)[nrow(contrast)]=levels[i]
			}
		}
	pdata <- phyDat(pmatrix, type="USER", contrast=contrast)

#Calculate tree length with Fitch (unordered) parsimony using ancestral.pars (phangorn package)
	tree$edge.length[tree$edge.length <= 0] <- 1e-05 #quick fix to the issue of negative branch lengths in some MCC trees (allows the code to run properly, but produces suboptimal graphical results)
	nrsteps <- parsimony(tree, pdata, method="fitch")

#Prepare tree stuff for fossil analyses
	tips <- tree$tip.label #Prepare a vector with the complete list of tips in "tree" ("tips")
	nrtips <- length(tree$tip.label)
	nredges <- nrtips + tree$Nnode - 1 #total number of branches (edges) in tree: number of terminal branches (tips) + number of internal branches (= number of internal nodes - 1); this would normally be equal to 2*nrtips-2 unless there are polytomies

#Loop through fossil taxa
	# fossiltaxon <- fossiltaxa[sample(1:nrfossils,1)]
	for (f in 1: nrfossils)		{
			fossiltaxon <- fossiltaxa[f]
			fossilchars <- pmatrix[fossiltaxon,]
			nrmissing <- length(fossilchars[fossilchars=="?"])
			nrchars <- ncol(pmatrix) - 1 - nrmissing
			singlefossiltree <- read.tree(text=paste("(", fossiltaxon, ":300);")) #note arbitrary branch length
			#Loop through all edges of the tree and get a parsimony score for each of them
				edgescores <- vector(mode="integer", length=0)
				for (i in 1:nredges)
					{
				  print(i)
						edgeendnode <- tree$edge[[i,2]] #edge and node numbers have nothing in common, unfortunately, hence the conversion (method borrowed from asrAndmore.R)
						treewithfossil <- bind.tree(tree, singlefossiltree, where=edgeendnode, position=tree$edge.length[i]/2)
						edgescores[i] <- parsimony(treewithfossil, pdata, method="fitch")
					}
				scorespread <- max(edgescores) - min(edgescores) + 1
				colorgradient <- c("black", "green", colorRampPalette(c("yellow", "grey"))(scorespread-2))
				fossilTree(shape="long", tree, fossiltaxon, nrchars, edgescores, colorgradient,  subfolder="", graphpar=c(80,8,20,0.05))
				fossilTree(shape="fan", tree, fossiltaxon, nrchars, edgescores, colorgradient,  subfolder="", graphpar=c(80,8,20,0.05))
		}
}