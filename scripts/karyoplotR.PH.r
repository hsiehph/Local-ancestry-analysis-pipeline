#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=3)
library(karyoploteR)

#args = c("hifi_HG00514/RFMix1KG5ContinentPops/", "tmp")

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

#  mapping_pop_color = data.frame(popid=c("AFR","EUR","EAS","AMR","SAS"), col=c('#ffff33','#377eb8','#e41a1c','#e41a1c','#377eb8'))
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

cols <- levels(factor(df_allchr$hap0_col))
popID <- NULL
for (ii in seq(1,length(cols))){
  popID <- c(popID, as.character(df_allchr[ df_allchr$hap0_col == cols[ii], "hap0_anc"][1]))
}


tmp <- levels(df_allchr$chr)
chromosomes <- tmp[order(as.numeric(gsub("[^0-9]+","",tmp)))]

num_chrs <- length(chromosomes)
plot.params <- getDefaultPlotParams(plot.type=1)
plot.params$ideogramheight <- 30
  
if (num_chrs < 12){
  pdf(paste(args[2], ".PH.pdf",sep=""), height=9, width=12)
  kp <- plotKaryotype(chromosomes=chromosomes, plot.params = plot.params)
  kpPlotRegions(kp, data=df_allchr[ df_allchr$density>0.0001,], col=df_allchr[ df_allchr$density>0.0001,"hap0_col"], border=NULL, r0=0, r1=0.3)
  kpAddLabels(kp, labels="Hap0", r0=0, r1=0.3)
  kpPlotRegions(kp, data=df_allchr[ df_allchr$density>0.0001,], col=df_allchr[ df_allchr$density>0.0001,"hap1_col"], border=NULL, r0=0.4, r1=0.7)
  kpAddLabels(kp, labels="Hap1", r0=0.4, r1=0.7)
  legend("topright", legend=popID, fill=cols, bty="n", horiz=TRUE)
  dev.off()
} else{
  part1END_chrs <- round(num_chrs/2)
  
  pdf(paste(args[2],".part1.PH.pdf",sep=""), height=9, width=12)
  kp <- plotKaryotype(chromosomes=chromosomes[1:part1END_chrs], plot.params = plot.params)
  part1_df_allchr <- df_allchr[ df_allchr$chr %in% chromosomes[1:part1END_chrs] & df_allchr$density>0.0001, ]
  kpPlotRegions(kp, data=part1_df_allchr, col=part1_df_allchr$hap0_col, border=NULL, r0=0, r1=0.3)
  kpAddLabels(kp, labels="Hap0", r0=0, r1=0.3)
  kpPlotRegions(kp, data=part1_df_allchr, col=part1_df_allchr$hap1_col, border=NULL, r0=0.4, r1=0.7)
  kpAddLabels(kp, labels="Hap1", r0=0.4, r1=0.7)
  legend("bottomright", legend=popID, fill=cols, bty="n", horiz=TRUE)
  dev.off()
  
  pdf(paste(args[2],".part2.PH.pdf",sep=""), height=9, width=12)
  kp <- plotKaryotype(chromosomes=chromosomes[(part1END_chrs+1):num_chrs], plot.params = plot.params)
  part2_df_allchr <- df_allchr[ df_allchr$chr %in% chromosomes[ (part1END_chrs+1):num_chrs] & df_allchr$density>0.0001, ]
  kpPlotRegions(kp, data=part2_df_allchr, col=part2_df_allchr$hap0_col, border=NULL, r0=0, r1=0.3)
  kpAddLabels(kp, labels="Hap0", r0=0, r1=0.3)
  kpPlotRegions(kp, data=part2_df_allchr, col=part2_df_allchr$hap1_col, border=NULL, r0=0.4, r1=0.7)
  kpAddLabels(kp, labels="Hap1", r0=0.4, r1=0.7)
  legend("bottomright", legend=popID, fill=cols, bty="n", horiz=TRUE)
  dev.off()
} 

