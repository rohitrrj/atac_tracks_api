#* Return the Track
#* @param distflank Zoom out number
#* @param Genes Gene name / locus
#* @param Track_list Named list of URLs for BigWig files
#* @param Track_cols Sequence of colors for the tracks
#* @post /atacTracks
#* @get /atacTracks
#* @serializer rds
atacTracks = function(req){
  source("getGeneIDs.R")
  
  # post body
  body <- jsonlite::fromJSON(req$postBody)

  distflank<-body$distflank
  Genes<-body$Genes
  Track_list<-body$Track_list
  Track_cols<-body$Track_cols
  
  distflank<-as.numeric(distflank)
  if(grepl(":", Genes)){
    message("This is region")
    Genes_toplot_gr<-parse2GRanges(Genes)
    Genes_toplot_gr<-Genes_toplot_gr + distflank
    ids <- getGeneIDsFromTxDb_RRJ(Genes_toplot_gr, TxDb.Hsapiens.UCSC.hg19.knownGene)
    symbols <- mget(ids, org.Hs.egSYMBOL)
    if(length(symbols) > 0){
      genes <- geneTrack(ids, TxDb.Hsapiens.UCSC.hg19.knownGene, 
                         symbols, asList=FALSE)
      auto <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC",group="auto")
      genes@dat<-genes@dat[genes@dat@seqnames %in% auto,]
      Genes_toplot<-genes
      setTrackStyleParam(Genes_toplot, "ylabpos", 'upstream')
      setTrackStyleParam(Genes_toplot, "color", 'black')
      # setTrackYaxisParam(Genes_toplot, "gp", list(col = "black", lty = "solid", lwd = 3, fontsize = 12))
      setTrackYaxisParam(Genes_toplot, "gp", list(col = "black"))
      eval(parse(text=(paste("obj<-list(","Locus","=Genes_toplot)",sep=""))))
    }
    
  } else {
    message("This is gene")
    Genes_toplot <- Genes
    message(Genes_toplot)
    entrezIDforGenes_toplot <- get(Genes_toplot, org.Hs.egSYMBOL2EG)
    Genes_toplot_gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene, single.strand.genes.only=FALSE)[entrezIDforGenes_toplot]
    Genes_toplot_gr <- keepStandardChromosomes(Genes_toplot_gr,pruning.mode="tidy")
    seqlevels(Genes_toplot_gr) = as.character(unique(seqnames(Genes_toplot_gr)))
    Genes_toplot_gr<-unlist(Genes_toplot_gr) + distflank
    ids <- getGeneIDsFromTxDb_RRJ(Genes_toplot_gr, TxDb.Hsapiens.UCSC.hg19.knownGene)
    symbols <- mget(ids, org.Hs.egSYMBOL)
    genes <- geneTrack(ids, TxDb.Hsapiens.UCSC.hg19.knownGene, 
                       symbols, asList=FALSE)
    auto <- extractSeqlevelsByGroup(species="Homo sapiens", style="UCSC",group="auto")
    genes@dat<-genes@dat[genes@dat@seqnames %in% auto,]
    Genes_toplot<-genes
    setTrackStyleParam(Genes_toplot, "ylabpos", 'upstream')
    setTrackStyleParam(Genes_toplot, "color", 'black')
    # setTrackYaxisParam(Genes_toplot, "gp", list(col = "black", lty = "solid", lwd = 3, fontsize = 12))
    setTrackYaxisParam(Genes_toplot, "gp", list(col = "black"))
    eval(parse(text=(paste("obj<-list(","Locus","=Genes_toplot)",sep=""))))
  }
  
  seqlevelsStyle(Genes_toplot_gr) <- "UCSC"
  
  AllSamples<-Track_list
  j=1
  for(i in AllSamples){
    message(parse(text=(paste0(names(AllSamples)[j]," <- importScore(file = \"",i,"\", format=\"BigWig\",ranges = Genes_toplot_gr)"))))
    eval(parse(text=(paste0(names(AllSamples)[j]," <- importScore(file = \"",i,"\", format=\"BigWig\",ranges = Genes_toplot_gr)"))))
    message(parse(text=(paste0("setTrackStyleParam(",names(AllSamples)[j],", \"color\", c(\"",Track_cols[j],"\",\"",Track_cols[j],"\"))"))))
    eval(parse(text=(paste0("setTrackStyleParam(",names(AllSamples)[j],", \"color\", c(\"",Track_cols[j],"\",\"",Track_cols[j],"\"))"))))
    j = j+1
  }
  
  AllSamplesobj_score<-NULL
  for(i in 1:length(AllSamples)){
    if(i == length(AllSamples)){
      AllSamplesobj_score<-paste0(AllSamplesobj_score,parse(text=(names(AllSamples)[i])),"$dat$score")
    }
    else{
      AllSamplesobj_score<-paste0(AllSamplesobj_score,parse(text=(names(AllSamples)[i])),"$dat$score",",")
    }
  }
  
  y_max<-0
  eval(parse(text=(paste0("y_max<-ceiling(max(c(",AllSamplesobj_score,")))"))))
  
  for(i in names(AllSamples)){
    eval(parse(text=(paste0("setTrackStyleParam(",i,", \"ylim\", c(0,y_max))"))))
  }
  
  AllSamplesobj<-NULL
  for(i in 1:length(AllSamples)){
    if(i == length(AllSamples)){
      AllSamplesobj<-paste0(AllSamplesobj,parse(text=(names(AllSamples)[i])))
    }
    else{
      AllSamplesobj<-paste0(AllSamplesobj,parse(text=(names(AllSamples)[i])),",")
    }
  }
  
  if(length(symbols) > 0){
    t <- try(Refseq_Genes <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                               org.Hs.eg.db,
                                               gr=Genes_toplot_gr))
    if ("try-error" %in% class(t)){
      eval(parse(text=(paste0("optSty <- optimizeStyle(trackList(obj[1],",AllSamplesobj,"), theme=NULL)"))))
      # optSty <- optimizeStyle(trackList(HC_1, HC_2, HC_3, HC_4, CAD_1, CAD_2, CAD_3, CAD_4), theme=NULL)
      trackList <- optSty$tracks
      viewerStyle <- optSty$style
      # vp <- viewTracks(trackList, gr=Genes_toplot_gr, viewerStyle=viewerStyle)
    }else {
      Refseq_Genes_names<-c()
      for (i in 1:length(Refseq_Genes)){
        setTrackStyleParam(Refseq_Genes[[i]], "ylabpos", "upstream")
        setTrackStyleParam(Refseq_Genes[[i]], "ylabgp", list(cex=.6))
        setTrackStyleParam(Refseq_Genes[[i]], "color", 'black')
        setTrackYaxisParam(Refseq_Genes[[i]], "gp", list(col = "black", lty = "solid", lwd = 3, fontsize = 32))
        Refseq_Genes_names<-c(Refseq_Genes_names,paste0(Refseq_Genes[[i]]$dat$symbol,"::",Refseq_Genes[[i]]$dat$transcript)[1])
      }
      names(Refseq_Genes)<-Refseq_Genes_names
      
      # eval(parse(text=(paste("obj<-list(",Genes,"=Genes_toplot,Transcripts=Refseq_Genes)",sep=""))))
    }
    eval(parse(text=(paste("obj<-list(","Transcripts=Refseq_Genes)",sep=""))))
    
    # eval(parse(text=(paste("obj<-list(",Genes,"=Genes_toplot)",sep=""))))
    
    eval(parse(text=(paste0("optSty <- optimizeStyle(trackList(obj[1],",AllSamplesobj,"), theme=NULL)"))))
    trackList <- optSty$tracks
    viewerStyle <- optSty$style
    # vp <- viewTracks(trackList, gr=Genes_toplot_gr, viewerStyle=viewerStyle)
  } else {
    # eval(parse(text=(paste("obj<-list(","Transcripts=Refseq_Genes)",sep=""))))
    eval(parse(text=(paste0("optSty <- optimizeStyle(trackList(",AllSamplesobj,"), theme=NULL)"))))
    trackList <- optSty$tracks
    viewerStyle <- optSty$style
    setTrackViewerStyleParam(viewerStyle, "xaxis", TRUE)
    setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .01, .01))
    # vp <- viewTracks(trackList, gr=Genes_toplot_gr, viewerStyle=viewerStyle)
  }
  track_args<-list(
    tracks = trackList,
    geneRegion = Genes_toplot_gr,
    view = viewerStyle
  )
}

print('app.R running')