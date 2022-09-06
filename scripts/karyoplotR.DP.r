#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=3)
library(karyoploteR)

#args = c("hifi_HG00514/RFMix1KG5ContinentPops/", "tmp.pdf")

if (length(args) == 2) {
  list_inputs <- list.files(args[1], full.names = T, pattern="msp.tsv")
} else {
  list_inputs <- list.files(args[1], full.names = T, pattern=paste(args[3], "msp.tsv", sep="."))
}

df_allchr <- NULL
for (i in seq(1,length(list_inputs))){
  header <- readLines(list_inputs[i],n=1)
  l_header <- strsplit(header,"\t| ")[[1]]
  mapPop2ID <- as.data.frame(do.call(rbind, sapply(l_header[ 3: length(l_header)], function(x) strsplit(x,"="))), row.names = F)

#  mapping_pop_color = data.frame(popid=c("AFR","EUR","EAS","AMR","SAS"), col=c('#ff7f00','#377eb8','#e41a1c','#e41a1c','#377eb8'))
  mapping_pop_color = data.frame(popid=c("AFR","EUR","EAS","AMR","SAS","NDL","DNS"), col=c('#ff7f00','#377eb8','#4daf4a','#e41a1c','#984ea3','black','#ffff33'))
  mapping_pop_color$pop <- mapPop2ID[ match(mapping_pop_color$popid, mapPop2ID$V1), 2]
  mapping_pop_color$col <- as.character(mapping_pop_color$col)
  
  
  dat <- read.delim(list_inputs[i], skip=1)
  dat$spos <- dat$spos+1
  dat$len <- dat$epos - dat$spos
  dat$density <- dat$n.snps / dat$len

  dat$hap0_anc <- mapping_pop_color[ match(dat[,7], mapping_pop_color$pop), "popid"]
  dat$hap1_anc <- mapping_pop_color[ match(dat[,8], mapping_pop_color$pop), "popid"]
  dat$hap0_col <- mapping_pop_color[ match(dat[,7], mapping_pop_color$pop), "col"]
  dat$hap1_col <- mapping_pop_color[ match(dat[,8], mapping_pop_color$pop), "col"]

  df = data.frame(chr=dat$X.chm, start=dat$spos, end=dat$epos, hap0_col=dat$hap0_col, hap1_col=dat$hap1_col, 
                  hap0_anc=dat$hap0_anc, hap1_anc=dat$hap1_anc, length=dat$len, n.snps=dat$n.snps, density=dat$density)
  df$hap0_col <- as.character(df$hap0_col)
  df$hap1_col <- as.character(df$hap1_col)

  df_allchr <- rbind(df_allchr, df)
    
}

df_allchr <- df_allchr[ df_allchr$density > 0.0001,]

############ David P's ideogram plotting function ###############
#' Function to plot assignments of haplotypes to 1000G population
#' 
#'
#' @param data.tab Table of assignments of H1 assembly to the reference haplotypes.
#' @param genome A \pkg{biovizBase} reference genome id. default: 'hg38'
#' @param chromosomes User defined set of chromosomes to plot.
#' @param title A \code{character} to use as a title of the plot.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
#'
plotChromPaints <- function(data.tab, genome = 'hg38', title=NULL) {
  ## Required libraries
  library(ggnewscale)
  library(ggbio)
  library(biovizBase)
  library(ggplot2)
  ## Load the data
#  data.df <- utils::read.table(data.tab, header = TRUE, stringsAsFactors = FALSE, comment.char = "&")
  data.df <- data.tab
  colnames(data.df) <- gsub(colnames(data.df), pattern = 'chr|chromosome', replacement = 'seqnames')
  
  ## Get coressponding ideogram from the database
  suppressMessages( hg38IdeogramCyto <- biovizBase::getIdeogram(genome, cytobands = TRUE) )
  
  ## Remove chromosomes that are not defined in chromosomes
  tmp <- levels(df_allchr$chr)
  chromosomes <- tmp[order(as.numeric(gsub("[^0-9]+","",tmp)))]

  if (!is.null(chromosomes)) {
    hg38IdeogramCyto <- keepSeqlevels(hg38IdeogramCyto, value = chromosomes, pruning.mode = 'coarse')
    seqlevels(hg38IdeogramCyto) <- chromosomes
    
    data.df <-  data.df[data.df$seqnames %in% chromosomes,]
    data.df$seqnames <- factor(data.df$seqnames, levels=chromosomes)
  }  
  
  ## Set the ggplot theme for final chromosome ideogram
  theme_vertical <- theme(legend.position ="top",
                          axis.line = element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.text.y=element_text(size=12),
                          strip.background = element_blank(),
                          strip.text.y = element_text(angle = 180),
                          strip.text.x = element_text(size = 12, margin = margin(0,2,0,2, "cm")),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank(),
                          panel.spacing.x=unit(0, "lines"),
                          legend.text=element_text(size=10))
  
  
  ## Convert ideogram bands stored in GRanges object into the data.frame
  ideo.df <- as.data.frame(hg38IdeogramCyto)
  
  ## Get chromosome breaks and labels
  max.len <- signif(max(ideo.df $end), digits = 2)
  breaks <- seq(from = 0, to = max.len, length.out = 6)
  labels <- breaks / 1000000
  labels <- paste0(labels, 'Mb')
  chr.num <- length(unique(ideo.df $seqnames))
  
  ## Plot the ideogram using ggplot2
  ideo <- ggplot() + 
    geom_rect(data=ideo.df, aes(ymin=start, ymax=end, xmin=0.25, xmax=0.75, fill=gieStain), color='black', show.legend=FALSE) + 
    scale_fill_giemsa() +
    facet_grid(. ~ seqnames, switch = 'x') + 
    scale_y_continuous(breaks = breaks, labels = labels) +
    theme_vertical +
    xlab("") +
    ylab("")
  
  ## Prepare data for plotting
  plt.df <- reshape2::melt(data.df, measure.vars = c('hap0_col', 'hap1_col'))
  levels(plt.df$variable) <- c('hap0', 'hap1')
  plt.df$xmin <- 0
  plt.df$xmax <- 0
  plt.df$xmin[plt.df$variable == 'hap0'] <- 0.9
  plt.df$xmax[plt.df$variable == 'hap0'] <- 1.75
  plt.df$xmin[plt.df$variable == 'hap1'] <- -0.75
  plt.df$xmax[plt.df$variable == 'hap1'] <- 0.1
  
  cols <- levels(factor(plt.df$value))
  popID <- NULL
  for (ii in seq(1,length(cols))){
    popID <- c(popID, as.character(plt.df[ plt.df$value == cols[ii], "hap0_anc"][1]))
  }
  ## Convert user supplied GRanges object into the data.frame
  plt <- ideo + 
    new_scale("fill") + 
    geom_rect(data=plt.df, aes(ymin=start, ymax=end, xmin=xmin, xmax=xmax, fill=plt.df$value)) + 
      scale_fill_manual(values=cols, labels=popID, name="") + xlim(-2,2)
  
  ## Decrese spacing between facets  
  plt <- plt + theme(panel.spacing.x=unit(0.5, "lines") , panel.spacing.y=unit(1,"lines"))
  
  ## Add title
  if (!is.null(title)) {
    plt <- plt + ggtitle(title)
  }
  return(plt)
}

################### End of David P's plotting function ########################


pdf(args[2], width=12, height=9)
plotChromPaints(data.tab = df_allchr, genome = 'hg38', title = "")
#plotChromPaints(data.tab = df_allchr, genome = 'hg38', chromosomes = paste0('chr', c(1:22)), title = "")
dev.off()
