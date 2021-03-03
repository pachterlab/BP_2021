# Get custom intron exon index with the 1kb at 3' end
library(GenomicRanges)
library(GenomicFeatures)
library(biomartr)
library(plyranges)
library(Biostrings)
library(BUSpaRse)

gtf <- getGTF(db = "ensembl", "Gallus gallus")
gn <- getGenome(db = "ensembl", "Gallus gallus")
gtf <- plyranges::read_gff(gtf)
gn <- readDNAStringSet(gn)
names(gn) <- str_remove(names(gn), " dna:.*")

# Strategy: Create a pseudo-exon by the 3' end of each transcript, and get union of this flanking region and the exons.

tx <- gtf[gtf$type == "transcript"]
fls <- flank(tx, 1000, start = FALSE)
fls <- split(fls, fls$transcript_id)

tr2g <- tr2g_ensembl("Gallus gallus", other_attrs = "description", use_transcript_version = FALSE, use_gene_version = FALSE)
tr2g$gene_name[tr2g$gene_name == ""] <- tr2g$gene[tr2g$gene_name == ""]
write_tsv(tr2g, "./output/BronnerLab_Sha/tr2g_gallus.tsv")

exons <- gtf[gtf$type == "exon"]
exons <- split(exons, exons$transcript_id)
exons_fl <- union(exons, fls)
exons_fl <- revElements(exons_fl, any(strand(exons_fl) == "-"))
all.equal(lengths(exons), lengths(exons_fl))

exons_fl <- unlist(exons_fl)
exons_fl$type <- "exon"
exons_fl$transcript_id <- names(exons_fl)
exons_fl$gene_id <- tr2g$gene[match(names(exons_fl), tr2g$transcript)]
names(exons_fl) <- NULL

# Get rid of negative ranges
start(exons_fl)[start(exons_fl) < 1] <- 1
end(exons_fl)[end(exons_fl) < 1] <- 1

# Get rid of ranges running outside the chromosome
ws <- setNames(width(gn), names(gn))
max_range <- ws[as.vector(seqnames(exons_fl))]
start(exons_fl)[start(exons_fl) > max_range] <- max_range[start(exons_fl) > max_range]
end(exons_fl)[end(exons_fl) > max_range] <- max_range[end(exons_fl) > max_range]

chroms_use <- seqlevels(exons_fl)[1:35]
exons_fl <- keepSeqlevels(exons_fl, chroms_use, pruning.mode = "coarse")
gn <- gn[chroms_use]
isCircular(gn) <- isCircular(exons_fl) <- setNames(chroms_use == "MT", chroms_use)

# I'll use gene IDs here, to avoid duplicates
get_velocity_files(exons_fl, L = 30, Genome = gn, out_path = "./output/BronnerLab_Sha",
                   transcript_version = NULL, gene_version = NULL)

# Write transcriptome fasta file for non-RNA-velocity application
exons_tx <- split(exons_fl, exons_fl$transcript_id)
tx <- extractTranscriptSeqs(gn, exons_tx)
writeXStringSet(tx, "./output/BronnerLab_Sha/gallus1k.fa")
