library(seqinr)
library(future)
library(future.apply)
library(glmnet)
library(e1071)
library(Biostrings)

gff_eco = "./GCF_000005845.2_ASM584v2_genomic.gff"
fasta_eco = "./GCF_000005845.2_ASM584v2_genomic.fna"

gff_yeast = "./GCA_000146045.2_R64_genomic.gff"
fasta_yeast = "./GCA_000146045.2_R64_genomic.fna"

get_regions <- function(row,chrom_sizes,is_circular){
  
  chrom_name <- row["seqname"][[1]]
  start <- as.integer(row["start"][[1]]) 
  end <- as.integer(row["end"][[1]])
  strand <- row["strand"][[1]]
  boundary <- chrom_sizes[chrom_name]
  
  if (is_circular) {
    if (strand == "+") {
      region_fw <- rbind(chr = chrom_name, 
                         start_sequence = ifelse(start-300 >= 1,start-300,(boundary-300)+start),
                         end_sequence = start)
      return(region_fw)
    } else {
      region_rw <- rbind(chr = chrom_name, 
                         start_sequence = end,
                         end_sequence = ifelse(end+300 <= boundary,end+300,(end+300)-boundary))
      return(region_rw)
    }
    
  } else {
    if (strand == "+") {
      region_fw <- rbind(chr = chrom_name, 
                         start_sequence = ifelse(start-300 >= 1,start-300+1,1),
                         end_sequence = start)
      return(region_fw)
    } else {
      region_rw <- rbind(chr = chrom_name, 
                         start_sequence = end,
                         end_sequence = ifelse(end+300 <= boundary,end+300-1,boundary))
      return(region_rw)
    }
  }
}

extract_seqs <- function(FASTA,row){
  chr <- row[1]
  start <- row[2]
  end <- row[3]
  gene_name <- rownames(row)
  sequence <- substr(FASTA[chr], start, end)
  return(c(chromosome = chr, gene_name = gene_name, sequence = sequence))
}

get_sequences <- function(gff,fasta,is_circular){
  
  data <- read.table(gff_yeast, 
                     sep = "\t", 
                     header = FALSE, 
                     quote = "",
                     comment.char = "#",
                     stringsAsFactors = FALSE)
  
  colnames(data) <- c("seqname",
                      "source",
                      "feature",
                      "start",
                      "end",
                      "score",
                      "strand",
                      "frame",
                      "group")
  
  chrom_size <- data[data$feature == "region",]$end
  names(chrom_size) <- data$seqname[data$feature == "region"]
  
  genes <- data[data$feature == "gene",]
  
  tmp <- strsplit(x=genes$group, split=";")
  GeneName <- unlist(lapply(X = tmp, function(x)x[1])) 
  GeneName <- sub(pattern="ID=gene-",
                  x=GeneName,
                  replacement = "")
  
  regions_complete <- as.data.frame(t(apply(X = genes,
                                            MARGIN = 1,
                                            FUN = get_regions,
                                            chrom_size,
                                            is_circular)))
  
  colnames(regions_complete) <- c("chromosome", "start", "end")
  rownames(regions_complete) <- c(GeneName)
  
  regions_complete$start <- as.integer(regions_complete$start)
  regions_complete$end <- as.integer(regions_complete$end)
  regions_complete <- regions_complete[regions_complete$start <= regions_complete$end,]
  #names_vector <- paste0("chr",seq_along(chrom_size))
  #names(names_vector) <- names(chrom_size)
  #regions_complete$chromosome <- names_vector[regions_complete$chromosome]
  
  FASTA <- Biostrings::readDNAStringSet(fasta)
  FASTA <- as.character(FASTA)
  
  new_names <- sapply(strsplit(names(FASTA), " "), `[`, 1)
  #`[`, 1 extracts first element of each sublist
  names(FASTA) <- new_names
  
  
  all_seqs <- t(apply(X = regions_complete, MARGIN = 1, FUN = extract_seqs, FASTA = FASTA))
  colnames(all_seqs) <- c('Chr','Sequence')
  return(all_seqs)}
sequences <- (get_sequences(gff_yeast,fasta_yeast,is_circular))
unique(sequences[,1])
unique(nchar(sequences[,2]))
