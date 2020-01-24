#Phyloscanning for Paleo-eFLOWER (Hervé Sauquet, Apr 2017)

#Clear buffer of all objects
	rm(list=ls())

#Load packages
	library(ape) #framework for handling phylogenetic trees
	library(corHMM) #used for ML ancestral state reconstruction (rayDISC function)
	library(phangorn) #used for MP ancestral state reconstruction (ancestral.pars function)
	library(geiger) #required for function tips
	library(plotrix) #required for plotting tables
	# library(xlsx) #required for writing summary results in multi-sheet Excel workbooks
	library(numbers) #required for calculating Bell numbers
	library(boa) #required for calculating 95% HPDs
	library(plyr) #required for calculating state frequencies using the count function

#Load custom functions
	sapply(list.files(pattern="[.]R$", path="./Functions/", full.names=TRUE), source);

#Read data
	dsversion <- "V7"
	series <- "A"
	suffix <- ""
	treefilename <- paste("eFLOWER_",series,"_MCC_updated.phy",sep="")
	usrdata <- readData(
		datafolder="./Data/", 
		treefile=treefilename, 
		maptreefile="", 
		datafile="Paleo-eFLOWER_2020-01-22.csv", 
		charnamefile="Paleo-eFLOWER_2020-01-22_characters.csv", 
		ingroupdef=c("Amborella_trichopoda","Spathiostemon_javensis"), 
		cladedefsfile="CladeDefinitions_eFLOWER.csv",
		trimtotree=FALSE
	)

#Identify fossil taxa
	fossiltaxa_index <- grep("X_",usrdata$charmatrix[,1])
	fossiltaxa <- as.character(usrdata$charmatrix[fossiltaxa_index,1])
	nrfossils <- length(fossiltaxa)
	#Report on number of fossil taxa
	charmatrix <- usrdata$charmatrix
#	charmatrix <- usrdata$charmatrix[-fossiltaxa,]

#Prepare matrix in the strange format required by phangorn
	pmatrix <- as.matrix(charmatrix[,-1])
	rownames(pmatrix) <- charmatrix[,1]
	#Build a custom so-called 'contrast' matrix to define polymorphic states
		#! this may see, very complicated, but there is no other way: phangorn does not seem to have a native system to read polymorphic data for the custom "USER" type (e.g., morphological), unlike corHMM
		#! the code below assumed that polymorphic data have been provided as requested by corHMM (e.g., 0/1 = "0&1")
		#! the code below will only work for characters with symbols in "0123456789" (more work needed for allowing non-numeric symbols)
		#! the alternative is to transform all polymorphic data into missing data, but this is not the exact same inference; using contrasts is much better and replicates exactly what Mesquite does
		nrstates <- 5
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
	tree <- usrdata$treeUP
	nrsteps <- parsimony(tree, pdata, method="fitch")
	#Report on default tree length

#Prepare tree stuff for fossil analyses
	tips <- tree$tip.label #Prepare a vector with the complete list of tips in "tree" ("tips")
	nrtips <- length(tree$tip.label)
	nredges <- nrtips + tree$Nnode - 1 #total number of branches (edges) in tree: number of terminal branches (tips) + number of internal branches (= number of internal nodes - 1); this would normally be equal to 2*nrtips-2 unless there are polytomies

#Loop through fossil taxa
	# fossiltaxon <- fossiltaxa[sample(1:nrfossils,1)]
	for (f in 1: nrfossils)
		{
			fossiltaxon <- fossiltaxa[f]
			fossilchars <- charmatrix[f,]
			nrmissing <- length(fossilchars[fossilchars=="?"])
			nrchars <- ncol(charmatrix) - 1 - nrmissing
			singlefossiltree <- read.tree(text=paste("(", fossiltaxon, ":300);")) #note arbitrary branch length
			#Loop through all edges of the tree and get a parsimony score for each of them
				edgescores <- vector(mode="integer", length=0)
				# edgeco <- rep("grey",nredges)
				# plot(tree, edge.color=edgeco, show.tip.label=FALSE)
				for (i in 1:nredges)
					{
						edgeendnode <- tree$edge[[i,2]] #edge and node numbers have nothing in common, unfortunately, hence the conversion (method borrowed from asrAndmore.R)
						treewithfossil <- bind.tree(tree, singlefossiltree, where=edgeendnode, position=tree$edge.length[i]/2)
						edgescores[i] <- parsimony(treewithfossil, pdata, method="fitch")
						# edgeco[i] <- "blue"
					}
				scorespread <- max(edgescores) - min(edgescores) + 1
				# colorgradient <- c(topo.colors(3), rep("grey", scorespread-3))
				# colorgradient <- c("red", "yellow", colorRampPalette(c("blue", "cyan"))(scorespread-2))
				# colorgradient <- c("red", "gold1", colorRampPalette(c("blue", "cyan"))(scorespread-2))
				# colorgradient <- c("red", "blue", colorRampPalette(c("yellow", "grey"))(scorespread-2))
				# colorgradient <- c("red", "green", colorRampPalette(c("yellow", "grey"))(scorespread-2))
				# colorgradient <- c("black", "blue", colorRampPalette(c("yellow", "grey"))(scorespread-2))
				colorgradient <- c("black", "green", colorRampPalette(c("yellow", "grey"))(scorespread-2))
				# colorgradient <- c("black", "green", colorRampPalette(c("yellow", "white"))(scorespread-2))
				# colorgradient <- rainbow(scorespread)
				fossilTree(shape="long", tree, fossiltaxon, nrchars, edgescores, colorgradient, usrdata$mapclades, usrdata$mapnodesUP, subfolder="", graphpar=c(80,8,20,0.05))
				fossilTree(shape="fan", tree, fossiltaxon, nrchars, edgescores, colorgradient, usrdata$mapclades, usrdata$mapnodesUP, subfolder="", graphpar=c(80,8,20,0.05))
				# pie(rep(1, scorespread), col=colorgradient, labels=min(edgescores):max(edgescores))
				# edgeco <- colorgradient[edgescores-min(edgescores)+1]
				# plot(tree, edge.color=edgeco, show.tip.label=FALSE)
				# hist(edgescores, labels=TRUE, col=colorgradient)
		}
