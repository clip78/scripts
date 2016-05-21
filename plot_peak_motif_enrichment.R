library(Biostrings)
library(MotifDb)
library(R.utils)
library(gplots)
library(RColorBrewer)

#config
verbose=F

####################
# generate_enrichment_map
# - coverage - dataframe representing output from a bedtools coverage (coverageBed) in histogram mode on a set of peaks
# - window - size of the window of enrichment
generate_enrichment_map<-function(coverage,window=2001) {
  names=as.character(unique(coverage$V8))
  enrichment_map=matrix(nrow=length(names),ncol=window,data=0)
  
  for(i in 1:(length(names)-1)) {
    subset=coverage[coverage$V8==names[i],]
    scores=as.numeric(as.character(subset$V4))
    print(length(scores))
    if(length(scores)==window)
      enrichment_map[i,1:length(scores)]=scores
    else {
      print(names[i])
    }
  }
  enrichment_map
}

####################
# plot enrichment from a pre-generated coverage file


macs_peaks_enrichment<-read.table('SEG1default_summits_slop_coverage.bed')
enrichment_map=generate_enrichment_map(macs_peaks_enrichment)

macs_peaks<-read.table('SEG1default_summits.bed')
ordered_dna_map<-enrichment_map[order(macs_peaks$V5),]

pal<-colorRampPalette(rev(brewer.pal('YlOrRd',n=9)))
heatmap.2(log2(ordered_dna_map+1),Rowv=NA,Colv=NA,trace="none",density="none",col=pal,
          labRow=NA,labCol=NA,ylab="Peaks",xlab="Genome position (2000bp window around summit)",
          main="Log read count in SEG1 peaks (MACS1.4)")

