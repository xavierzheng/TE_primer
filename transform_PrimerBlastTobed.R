#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("Usage: transform_PrimerBlastTobed.R <blastn_short_outfmt7.table>

       # AIM: parser the blast outfmt7 table into bed format, only the region with 100% coverage with <2 mismatch and 0 gap were retained for bed output. 
       # Only amplicon size less than 500bp were output to potential_amplicon.bed. 

       # NEED THESE blastn parameters:
        -task blastn-short 
        -evalue 100 
        -outfmt 7 qseqid sseqid qlen slen sstrand pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

 	# Default output: primer_target_site.bed AND potential_amplicon.bed
       ", call.=FALSE)
}


#setwd("/nfs/project1/pxzhe/primer")

FileName <- as.character(x = args[1])

suppressPackageStartupMessages(library(tidyverse, quietly = F))
options(dplyr.summarise.inform=FALSE)

df <- read.table(
  file = FileName, 
  header = F, sep = "\t", comment.char = "#", stringsAsFactors = F
)

colnames(df) <- c("qseqid", "sseqid", "qlen", "slen", "sstrand", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")

df_out <- df %>%
  filter(
    qlen == length
  ) %>%
  filter(
    mismatch < 2
  ) %>%
  filter(
    gapopen == 0
  ) %>%
  select(
    sseqid, sstart, send, sstrand, qseqid, mismatch
  ) %>%
  mutate(
    BED_START = ifelse(
      sstrand == "plus", 
      sstart-1, 
      send-1
    ), 
    BED_END = ifelse(
      sstrand == "plus", 
      send, 
      sstart
    ), 
    BED_STRAND = ifelse(
      sstrand == "plus", 
      "+", 
      "-"
    ), 
    BED_SCORE = mismatch, 
    BED_CHROM = sseqid,
    BED_NAME = qseqid
  ) %>%
  select(
    BED_CHROM, BED_START, BED_END, BED_NAME, BED_SCORE, BED_STRAND
  ) %>%
  arrange(
    BED_CHROM, BED_START
  ) 

df_out %>%
  write.table(
    file = "primer_target_site.bed", 
    col.names = F, 
    row.names = F, 
    sep = "\t", 
    quote = F
  )

#=====================================
# amplicon

PCR_size_threshold <- as.numeric(500)

# cut _F _R and primer set name, set name = GRP
# df_GRP <- df_out %>%
#   mutate(
#     GRP = str_remove(BED_NAME, pattern = "[RF]$"), 
#     GRP_STRAND = str_remove(BED_NAME, pattern = "LTR\\d+"),
#   ) 

df_GRP <- df_out %>%
  mutate(
    GRP = str_split_fixed(BED_NAME, "_",2)[,1],
    GRP_STRAND = str_split_fixed(BED_NAME, "_", 2)[,2]
  )

#unique(df_GRP$GRP)

temp_list <- lapply(unique(df_GRP$GRP), FUN = function(NAM){
  
  # POS STRAND
  df_GRP_POS <- df_GRP %>%
    filter(
      GRP %in% NAM
    ) %>%
    filter(
      BED_STRAND %in% "+" & GRP_STRAND %in% "F" | BED_STRAND %in% "-" & GRP_STRAND %in% "R"
    ) %>%
    mutate(
      CHECK_POINT_POS = ifelse(
        BED_STRAND %in% "+" & GRP_STRAND %in% "F", 
        "TRUE", 
        "FALSE"
      ), 
      CHECK_POINT_NEG = ifelse(
        BED_STRAND %in% "-" & GRP_STRAND %in% "R", 
        "TRUE", 
        "FALSE"
      )
    ) %>%
    arrange(
      BED_CHROM, BED_START
    ) %>%
    mutate(
      NEXT_START = ifelse(
        GRP_STRAND %in% "F" & lead(GRP_STRAND) %in% "R", 
        lead(BED_START), 
        NA
      ), 
      NEXT_END = ifelse(
        GRP_STRAND %in% "F" & lead(GRP_STRAND) %in% "R", 
        lead(BED_END), 
        NA
      ), 
      NEXT_CHROM = ifelse(
        BED_CHROM == lead(BED_CHROM), 
        BED_CHROM, 
        NA
      ),
      AMPLICON_START = BED_START, 
      AMPLICON_END = NEXT_END, 
      AMPLICON_SIZE = NEXT_END - BED_START
    ) %>%
    filter(
      !is.na(NEXT_CHROM)
    ) %>%
    select(
      -starts_with("CHECK_POINT")
    ) %>%
    filter(
      !is.na(NEXT_START)
    ) 
  
  # NEG strand
  df_GRP_NEG <- df_GRP %>%
    filter(
      GRP %in% NAM
    ) %>%
    filter(
      BED_STRAND %in% "-" & GRP_STRAND %in% "F" | BED_STRAND %in% "+" & GRP_STRAND %in% "R"
    ) %>%
    mutate(
      CHECK_POINT_POS = ifelse(
        BED_STRAND %in% "-" & GRP_STRAND %in% "F", 
        "TRUE", 
        "FALSE"
      ), 
      CHECK_POINT_NEG = ifelse(
        BED_STRAND %in% "+" & GRP_STRAND %in% "R", 
        "TRUE", 
        "FALSE"
      )
    ) %>%
    arrange(
      BED_CHROM, BED_START
    ) %>%
    mutate(
      NEXT_START = ifelse(
        GRP_STRAND %in% "R" & lead(GRP_STRAND) %in% "F", 
        lead(BED_START), 
        NA
      ), 
      NEXT_END = ifelse(
        GRP_STRAND %in% "R" & lead(GRP_STRAND) %in% "F", 
        lead(BED_END), 
        NA
      ), 
      NEXT_CHROM = ifelse(
        BED_CHROM == lead(BED_CHROM), 
        BED_CHROM, 
        NA
      ),
      AMPLICON_START = BED_START, 
      AMPLICON_END = NEXT_END, 
      AMPLICON_SIZE = NEXT_END - BED_START
    ) %>%
    filter(
      !is.na(NEXT_CHROM)
    ) %>%
    select(
      -starts_with("CHECK_POINT")
    ) %>%
    filter(
      !is.na(NEXT_START)
    )
  
  df_GRP_bind <- rbind(df_GRP_POS, df_GRP_NEG)
  
  df_amplicon_out <- df_GRP_bind %>%
    filter(
      AMPLICON_SIZE < PCR_size_threshold
    ) %>%
    mutate(
      AMPLICON_STRAND = ifelse(
        GRP_STRAND %in% "F", 
        "+", 
        "-"
      )
    ) %>%
    select(
      BED_CHROM, AMPLICON_START, AMPLICON_END, GRP, AMPLICON_SIZE, AMPLICON_STRAND
    ) %>%
    mutate(
      BED_CHROM = as.character(BED_CHROM), 
      AMPLICON_START = as.numeric(AMPLICON_START), 
      AMPLICON_END = as.numeric(AMPLICON_END), 
      GRP = as.character(GRP), 
      AMPLICON_SIZE = as.numeric(AMPLICON_SIZE), 
      AMPLICON_STRAND = as.character(AMPLICON_STRAND)
    )
  
  return(df_amplicon_out)
  
})

df_amplicon_out <- bind_rows(temp_list)

write.table(
  df_amplicon_out, 
  file = "potential_amplicon.bed", 
  col.names = F, row.names = F, sep = "\t", quote = F
)

print("# potential No. of amplicon around whole genome")
df_amplicon_out %>%
  group_by(GRP) %>%
  summarise(
    COUNT = n()
  ) %>%
  write.table(
    file = "potential_amplicon_summary.txt", 
    col.names = T, row.names = F, sep = "\t", quote = F
  )



#=======================================
# # POS STRAND
# df_GRP_POS <- df_GRP %>%
#   filter(
#     GRP %in% "LTR1"
#   ) %>%
#   filter(
#     BED_STRAND %in% "+" & GRP_STRAND %in% "F" | BED_STRAND %in% "-" & GRP_STRAND %in% "R"
#   ) %>%
#   mutate(
#     CHECK_POINT_POS = ifelse(
#       BED_STRAND %in% "+" & GRP_STRAND %in% "F", 
#       "TRUE", 
#       "FALSE"
#     ), 
#     CHECK_POINT_NEG = ifelse(
#       BED_STRAND %in% "-" & GRP_STRAND %in% "R", 
#       "TRUE", 
#       "FALSE"
#     )
#   ) %>%
#   arrange(
#     BED_CHROM, BED_START
#   ) %>%
#   mutate(
#     NEXT_START = ifelse(
#       GRP_STRAND %in% "F" & lead(GRP_STRAND) %in% "R", 
#       lead(BED_START), 
#       NA
#     ), 
#     NEXT_END = ifelse(
#       GRP_STRAND %in% "F" & lead(GRP_STRAND) %in% "R", 
#       lead(BED_END), 
#       NA
#     ), 
#     NEXT_CHROM = ifelse(
#       BED_CHROM == lead(BED_CHROM), 
#       BED_CHROM, 
#       NA
#     ),
#     AMPLICON_START = BED_START, 
#     AMPLICON_END = NEXT_END, 
#     AMPLICON_SIZE = NEXT_END - BED_START
#   ) %>%
#   filter(
#     !is.na(NEXT_CHROM)
#   ) %>%
#   select(
#     -starts_with("CHECK_POINT")
#   ) %>%
#   filter(
#     !is.na(NEXT_START)
#   ) 
# 
# # NEG strand
# df_GRP_NEG <- df_GRP %>%
#   filter(
#     GRP %in% "LTR1"
#   ) %>%
#   filter(
#     BED_STRAND %in% "-" & GRP_STRAND %in% "F" | BED_STRAND %in% "+" & GRP_STRAND %in% "R"
#   ) %>%
#   mutate(
#     CHECK_POINT_POS = ifelse(
#       BED_STRAND %in% "-" & GRP_STRAND %in% "F", 
#       "TRUE", 
#       "FALSE"
#     ), 
#     CHECK_POINT_NEG = ifelse(
#       BED_STRAND %in% "+" & GRP_STRAND %in% "R", 
#       "TRUE", 
#       "FALSE"
#     )
#   ) %>%
#   arrange(
#     BED_CHROM, BED_START
#   ) %>%
#   mutate(
#     NEXT_START = ifelse(
#       GRP_STRAND %in% "R" & lead(GRP_STRAND) %in% "F", 
#       lead(BED_START), 
#       NA
#     ), 
#     NEXT_END = ifelse(
#       GRP_STRAND %in% "R" & lead(GRP_STRAND) %in% "F", 
#       lead(BED_END), 
#       NA
#     ), 
#     NEXT_CHROM = ifelse(
#       BED_CHROM == lead(BED_CHROM), 
#       BED_CHROM, 
#       NA
#     ),
#     AMPLICON_START = BED_START, 
#     AMPLICON_END = NEXT_END, 
#     AMPLICON_SIZE = NEXT_END - BED_START
#   ) %>%
#   filter(
#     !is.na(NEXT_CHROM)
#   ) %>%
#   select(
#     -starts_with("CHECK_POINT")
#   ) %>%
#   filter(
#     !is.na(NEXT_START)
#   )
# 
# df_GRP_bind <- rbind(df_GRP_POS, df_GRP_NEG)
# 
# df_amplicon_out <- df_GRP_bind %>%
#   filter(
#     AMPLICON_SIZE < PCR_size_threshold
#   ) %>%
#   mutate(
#     AMPLICON_STRAND = ifelse(
#       GRP_STRAND %in% "F", 
#       "+", 
#       "-"
#     )
#   ) %>%
#   select(
#     BED_CHROM, AMPLICON_START, AMPLICON_END, GRP, AMPLICON_SIZE, AMPLICON_STRAND
#   ) 
