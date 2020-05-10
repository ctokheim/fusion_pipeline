library(rMAUPS)
library(optparse)
library(stringr)
library(dplyr)

# command line
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="The path to the fusion proteins sequences", metavar="character"),
  make_option(c("-r", "--regex"), type="character", default=NULL,
              help="Regex motif", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="The path to output file", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# search prot sequence
prot_seq <- read.delim(opt$input, sep='\t', stringsAsFactors=F)
regex <- read.delim(opt$regex, sep='\t', stringsAsFactors=F)
#regex <- regex[str_detect(regex$ELMIdentifier, '^DEG_'),]

# search fusion transcripts for every degron motif
output_df <- data.frame()
for (i in 1:nrow(regex)) {
  # regex search
  #deg_id <- regex[i,'ELMIdentifier']
  #tmp_regex <- regex[i, 'Regex']
  deg_id <- regex[i,'Name']
  tmp_regex <- regex[i, 'Degron']
  tmp_result <- searchProtSeq(prot_seq, tmp_regex)

  # merge in gene id
  tmp_result <- left_join(tmp_result, prot_seq[c('ID', 'gene', 'transcript_id')], by=c('UniprotId' = "ID"))

  # concatenate results
  if (nrow(tmp_result)>0) {
    tmp_result['degron'] <- deg_id
    output_df <- rbind(output_df, tmp_result)
  }
}

write.table(output_df, opt$output, quote=F, row.names=F, sep='\t')
