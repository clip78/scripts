# Load the required libraries
library(SSOAP)		#To use WikiPathways web service
library(XML)
library(GEOquery)	#To download GEO experiments
library(GSEABase)	#Gene Set data structures
library(PGSEA)		#Parametric Gene Set Enrichment
library(limma)		#Used for mining the GSEA results

# Create a SOAPServer instance for the web service
srv = SOAPServer("http://www.wikipathways.org/wpi/webservice/webservice.php");

# Get a list of all pathways from WikiPathways
reply = .SOAP(
	srv, "listPathways", 
	action=I("listPathways"), handlers=NULL
)
# Parse the xml document
doc = xmlParse(reply$content, asText=TRUE)

# Extract the pathway info
pathwayNodes = xmlElementsByTagName(xmlRoot(doc), "pathways", TRUE)
pathways = lapply(pathwayNodes, function(n) {
	children = xmlChildren(n, addNames= TRUE)
	if(xmlValue(children$species) == "Homo sapiens") {
		p = list()
		p[["id"]] = xmlValue(children$id)
		p[["name"]] = xmlValue(children$name)
		p[["species"]] = xmlValue(children$species)
		p[["url"]] = xmlValue(children$url)
		return(p)
	} else {
		return() # Skip non-human pathways
	}
})

# Remove NULL entries (non-human pathways)
pathways = pathways[!sapply(pathways, is.null)]

# A function that downloads a GeneSet for a pathway
createGS = function(p) {
	print(p[["id"]])
	# Download the gene list (translated to Entrez ids)
	reply = .SOAP(srv, "getXrefList", pwId = p[["id"]], code="L", 
						action=I("getXrefList"), handlers=NULL)
	doc = xmlParse(reply$content, asText=TRUE)
	# Find the xref nodes with an xpath query
	resultNodes = xmlElementsByTagName(xmlRoot(doc), "xrefs", TRUE)
	# Extract the gene ids
	geneIds = sapply(resultNodes, xmlValue)
	if(length(geneIds) > 0) { # Skip empty lists
		# Create a GeneSet object
		geneSet = GeneSet(geneIds, geneIdType=EntrezIdentifier(), 
			setName=paste(p[["id"]], " (", p[["name"]], ")", sep="")
		)
		return(geneSet)
	}
}
geneSets = lapply(pathways, createGS) #Apply the createGS function on all pathways
geneSets = geneSets[!sapply(geneSets, is.null)] #Remove empty sets
geneSetCollection = GeneSetCollection(geneSets)

# Save the downloaded gene sets, so we can use them for later calculations
save(geneSetCollection, file="/home/thomas/code/webservice-examples/Hs_genesets.Rd")

# Get an expression set from GEO
gse = getGEO("GSE7023", GSEMatrix = TRUE)
eset = gse[[1]]

# Reformat the phenotype labels
subtype = gsub("subtype: ", "", phenoData(eset)$characteristics_ch1)
subtype = gsub("\\.", "_", subtype)

# Run the PGSEA calculation (reference sample is 'NO')
pg = PGSEA(eset, geneSetCollection, ref = which(subtype == "NO"), p.value=0.005)

# Plot the result matrix
smcPlot(
	pg, factor(subtype), scale=c(-15,15), show.grid = TRUE, margins = c(1, 1, 4, 10), 		col = .rwb, r.cex = 0.35
)

# Find the pathways most different between P1 and P2B
pgNF <- PGSEA(eset, geneSetCollection, ref = which(subtype == "NO"), p.value = NA)
design <- model.matrix(~-1 + factor(subtype))
colnames(design) <- names(table(subtype))
fit <- lmFit(pgNF, design)
contrast.matrix <- makeContrasts(P2B - P1, levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

# Plot the top 25
smcPlot(pg[as.numeric(rownames(topTable(fit, n = 30, resort.by = "logFC"))),], factor(subtype, levels = c("P1", "P2B")), col = .rwb, scale = c(-15, 15), margins = c(1, 1, 4, 10), show.grid = TRUE, r.cex = 0.75)
