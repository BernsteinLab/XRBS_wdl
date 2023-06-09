options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressMessages(library(ShortRead))


## expected barcodes
a <- commandArgs(trailingOnly=TRUE)[2]
# a <- "2"

a.split <- lapply(strsplit(strsplit(a, ",")[[1]], "-"), as.numeric)
a.split[lengths(a.split)==2] <- lapply(a.split[lengths(a.split)==2], function(x) seq(x[1], x[2]))
# a.split

a.unlist <- unique(sort(LETTERS[unlist(a.split)]))
# a.unlist


## known barcode sequences
# load
bc0 <- read.table("/seq/epiprod02/sgaldon/Alba/barcodes24.txt", header = TRUE)
bc0$letter <- LETTERS[seq(nrow(bc0))]
# bc0

# only expected barcodes
bc0 <- bc0[is.element(bc0$letter, a.unlist), ]
# bc0

# generate all barcode sequences with 1mm (2mm if only looking for single barcode)
per <- function(x) unlist(lapply(lapply(seq(nchar(x)), function(i) paste0(substr(x, 1, i-1), c("A", "C", "G", "T"), substr(x, i+1, 100))), setdiff, x))
bc0.per <- lapply(bc0$seq, per)
if(length(bc0.per)==1) bc0.per <- list(setdiff(unique(sort(unlist(lapply(bc0.per[[1]], per)))), bc0$seq))  # mismatch of mismatch = 2mm
# table(adist(bc1$seq, bc0$seq)[, 1])
bc1 <- bc0[rep(seq(length(bc0.per)), lengths(bc0.per)), ]
bc1$seq <- unlist(bc0.per)
# bc1

# remove sequences that occur multiple times
bc01 <- rbind(bc0, bc1)
bc01 <- bc01[!is.element(bc01$seq, names(which(table(bc01$seq) > 1))), ]
bc01 <- bc01[order(bc01$seq), ]
# bc01

# create vector
bc <- bc01$letter
names(bc) <- bc01$seq
bc <- bc[is.element(bc, a.unlist)]
# bc


## read and write fastq files at the same time
f1 <- commandArgs(trailingOnly=TRUE)[1]
# f1 <- "biotin-1M_R1.fastq.gz"
f2 <- sub("_R1", "_R2", f1)
f <- strsplit(f1, "_R1")[[1]][1]

# matrix to collect stats per barcode
n <- matrix(0, length(a.unlist)+1, 2, dimnames = list(paste0(f, "-", c(a.unlist, "NA")), c("all", "ygg")))

# loop
strm1 <- FastqStreamer(f1, n=1E6)  # 1M reads by default
strm2 <- FastqStreamer(f2, n=1E6)
repeat {
  fq1 <- yield(strm1)
  if(length(fq1) == 0) break
  fq2 <- yield(strm2)
  
  # match to barcodes
  fq2.eight <- as.vector(subseq(sread(fq2), 1, 8))
  fq2.bc <- bc[match(fq2.eight, names(bc))]
  fq2.ygg <- is.element(as.vector(subseq(sread(fq2), 9, 11)), c("CGG", "TGG"))  # sequence following barcode (informative CpG in MspI cut site)
  
  # count barcode occurence
  n[, 1] <- n[, 1] + table(factor(fq2.bc, a.unlist), useNA="a")
  n[, 2] <- n[, 2] + table(factor(fq2.bc[fq2.ygg], a.unlist), useNA="a")
  
  # shorten R2 and add sample barcode to id
  fq1@id <- fq2@id <- BStringSet(paste0(sub(" .:N:0:", ":", as.vector(fq2@id)), ":", fq2.eight, "+", as.vector(subseq(sread(fq1), 1, 6))))
  fq1 <- narrow(fq1, 7, width(fq1))  # remove hexamer sequence
  fq2 <- narrow(fq2, 9, width(fq2))  # remove barcode sequence

  # write uncompressed fastq
  for(i in a.unlist) {
    writeFastq(fq1[which(fq2.bc == i)], file = paste0(f, "-", i, "_R2.fastq"), mode = "a", compress = FALSE)  # append, switch read 1 and read 2
    writeFastq(fq2[which(fq2.bc == i)], file = paste0(f, "-", i, "_R1.fastq"), mode = "a", compress = FALSE)
  }
  writeFastq(fq1[which(is.na(fq2.bc))], file = paste0(f, "-NA_R2.fastq"), mode = "a", compress = FALSE)
  writeFastq(fq2[which(is.na(fq2.bc))], file = paste0(f, "-NA_R1.fastq"), mode = "a", compress = FALSE)
  invisible(gc())
}
close(strm1)
close(strm2)

# write summary
write.table(n, file = paste0(f, ".bcsplit.txt"), sep="\t", quote=FALSE)


