options(max.print = 1000)
options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressMessages(library(Rsamtools))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
# suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(data.table))

# bam file
# f <- "~/DropboxPartners/ShareefMethylationAnalysis/k562/test.bam"
# f <- "190610-190719.ociaml100-B.human.sort.bam"
f <- commandArgs(trailingOnly=TRUE)[1]
bamFile <- BamFile(f)

# stats accros all chromosomes
r <- list(mapped=0, illegal=0, mapq=0, pos=0, isize=0, duplicate=0)
r.pos <- r.cov.plus <- r.cov.minus <- list()

# loop over each chromosome
chr <- scanBamHeader(bamFile)$targets
chr <- chr[names(chr)!='chrL']

for(i in names(chr)) {
  # i <- "chr1"
  message("- processing ", i)
  
  # find MspI cut sites
  #chrL corresponds to lambda DNA
  message("  finding MspI cut sites: ", appendLF = FALSE)
  if(i == "chrL") {
    m <- matchPattern("CCGG", import("NC_001416.1.fasta")[[1]])
  } else {
    m <- matchPattern("CCGG", Hsapiens[[i]])
  }
  m <- start(m) + 2  # center of CCGG
  message(length(m), " sites found")
  
  # collect names of reads that are kept
  r.cov.plus[[i]] <- r.cov.minus[[i]] <- matrix(0, length(m), 10)  # 10 columns for downsampling analysis
  r.pos[[i]] <- data.frame(chr=i, pos=m)
  
  
  # pass 1, collect valid read ids on plus strand
  message("  scanning plus strand: ", appendLF = FALSE)
  bam.param <- ScanBamParam(which = GRanges(i, IRanges(1, chr[i])),
                            what = c("qname", "flag", "mapq", "rname", "pos", "isize", "cigar", "qual"),
                            flag = scanBamFlag(isUnmappedQuery = FALSE, isPaired = TRUE, isFirstMateRead = TRUE, isMinusStrand = FALSE))
  open(bamFile)
  
  # load data
  d <- as.data.table(scanBam(bamFile, param = bam.param)[[1]])
  if(nrow(d)==0) break
  
  # calculate average base quality and sort (needed later)
  d$qual <- sapply(as(d$qual, "IntegerList"), sum)
  d <- d[order(-qual, qname)]

  # number of reads
  r$mapped <- r$mapped + nrow(d)
  
  # illegal mapping
  d <- d[!bamFlagTest(flag, "isNotPassingQualityControls")]
  r$illegal <- r$illegal + nrow(d)
  
  # mapping quality
  d <- d[mapq > 0]
  r$mapq <- r$mapq + nrow(d)
  
  # mapping position, R1 needs to be at MspI cut site, correcting for soft clip
  x <- ifelse(substr(d$cigar, 2, 2) == "S", as.numeric(substr(d$cigar, 1, 1)), 0)
  d$pos <- d$pos - x
  d$isize <- d$isize + x
  d <- d[is.element(pos, m - 1)]
  r$pos <- r$pos + nrow(d)
  
  # insert size
  d <- d[isize >= 20 & isize <= 600]
  r$isize <- r$isize + nrow(d)
  
  # pcr duplicate / downsampling analysis. use hexamer sequence as UMI. take read with highest average base quality (reads are still sorted).
  # because reads are read in chunks of 1M, very few reads between chunks are not filtered even though they should be.
  set.seed(123)
  d$ran <- ceiling(runif(nrow(d), 0, 10))  # random number of downsampling
  d$umi <- paste0(d$pos, ":", d$isize, ":", substr(d$qname, nchar(d$qname)-5, 1000))
  for(j in seq(10)) {
    dd <- d[ran<=j]
    dd.split <- split(dd$qname, dd$umi)
    dd.split[lengths(dd.split)>1] <- lapply(dd.split[lengths(dd.split)>1], head, 1)  # only take read with highest quality if more than one
    dd <- dd[is.element(qname, unlist(dd.split))]
    r.cov.plus[[i]][, j] <- as.numeric(table(factor(match(dd$pos+1, m), seq(length(m)))))
  }
  r$duplicate <- r$duplicate + nrow(dd)
  
  # record the reads to keep
  r.qname.plus <- dd$qname

  message(length(r.qname.plus), " valid read pairs found")
  close(bamFile)
  
  rm(list=c("d", "dd"))
  invisible(gc())
  
  
  # pass 2, collect valid read ids from minus strand
  message("  scanning minus strand: ", appendLF = FALSE)
  bam.param <- ScanBamParam(which = GRanges(i, IRanges(1, chr[i])),
                            what = c("qname", "flag", "mapq", "rname", "pos", "isize", "cigar", "qual"),
                            flag = scanBamFlag(isUnmappedQuery = FALSE, isPaired = TRUE, isFirstMateRead = TRUE, isMinusStrand = TRUE))
  open(bamFile)
  
  # load data
  d <- as.data.table(scanBam(bamFile, param = bam.param)[[1]])
  if(nrow(d)==0) break
  
  # calculate average base quality and sort
  d$qual <- sapply(as(d$qual, "IntegerList"), sum)
  d <- d[order(-qual, qname)]
  
  # number of reads
  r$mapped <- r$mapped + nrow(d)
  
  # illegal mapping
  d <- d[!bamFlagTest(flag, "isNotPassingQualityControls")]
  r$illegal <- r$illegal + nrow(d)
  
  # mapping quality
  d <- d[mapq > 0]
  r$mapq <- r$mapq + nrow(d)
  
  # mapping position, R1 needs to be at MspI cut site, correcting for soft clip
  x <- ifelse(substr(d$cigar, nchar(d$cigar), nchar(d$cigar)) == "S" & is.element(substr(d$cigar, nchar(d$cigar)-2, nchar(d$cigar)-2), LETTERS),
              as.numeric(substr(d$cigar, nchar(d$cigar)-1, nchar(d$cigar)-1)), 0)
  d$pos <- d$pos + cigarWidthAlongReferenceSpace(d$cigar) + x
  d$isize <- d$isize - x
  d <- d[is.element(pos, m + 1)]
  r$pos <- r$pos + nrow(d)
  
  # insert size
  d <- d[isize <= -20 & isize >= -600]
  r$isize <- r$isize + nrow(d)
  
  # pcr duplicate / downsampling analysis. use hexamer sequence as UMI. take read with highest average base quality (reads are still sorted).
  # because reads are read in chunks of 1M, very few reads between chunks are not filtered even though they should be.
  set.seed(123)
  d$ran <- ceiling(runif(nrow(d), 0, 10))  # random number of downsampling
  d$umi <- paste0(d$pos, ":", -d$isize, ":", substr(d$qname, nchar(d$qname)-5, 1000))
  for(j in seq(10)) {
    dd <- d[ran<=j]
    dd.split <- split(dd$qname, dd$umi)
    dd.split[lengths(dd.split)>1] <- lapply(dd.split[lengths(dd.split)>1], head, 1)  # only take read with highest quality if more than one
    dd <- dd[is.element(qname, unlist(dd.split))]
    r.cov.minus[[i]][, j] <- as.numeric(table(factor(match(dd$pos-1, m), seq(length(m)))))
  }
  r$duplicate <- r$duplicate + nrow(dd)
  
  # record the reads to keep
  r.qname.minus <- dd$qname
  
  message(length(r.qname.minus), " valid read pairs found")
  close(bamFile)
  
  rm(list=c("d", "dd"))
  invisible(gc())
  
  
  # pass 3, filter bam file for valid read ids
  message("  filtering bam file: ", appendLF = FALSE)
  r.qname <- unique(sort(c(r.qname.plus, r.qname.minus)))
  
  filterBam(file = bamFile,
            destination = paste0(f, ".", i, ".filter.bam"),
            param = ScanBamParam(which = GRanges(i, IRanges(1, chr[i])),
                                 what = "qname"),
            filter = FilterRules(list(f = function(x) is.element(x$qname, r.qname))),
            indexDestination = FALSE)
  
  message(length(r.qname), " read pairs written")
}

# write final stats
message("- merging bam files: ", appendLF = FALSE)
write.table(t(unlist(r)), row.names=tail(strsplit(f, "/")[[1]], 1), file=paste0(f, ".filter.stats"), sep="\t", quote=FALSE)

write.table(data.frame(do.call(rbind, r.pos), strand="+", do.call(rbind, r.cov.plus)), file=paste0(f, ".filter.cov_plus"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(do.call(rbind, r.pos), strand="-", do.call(rbind, r.cov.minus)), file=paste0(f, ".filter.cov_minus"), sep="\t", quote=FALSE, row.names=FALSE)

# merge bam
mergeBam(files = paste0(f, ".", names(chr), ".filter.bam"), destination = paste0(f, ".filter.bam"), indexDestination = TRUE, overwrite = TRUE)
file.remove(paste0(f, ".", names(chr), ".filter.bam"))
message("done")



