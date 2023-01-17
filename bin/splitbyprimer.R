#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
location <- args[3]

library(readr)
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)


locations <- read_table(location)
fastq <- read_table(input, col_names = FALSE)

new_fq <- tibble(id = str_remove(fastq$X1[seq(from = 1, to = nrow(fastq)- 3, by = 4)], "^@"),
                 sequence = fastq$X1[seq(from = 2, to = nrow(fastq)- 2, by = 4)],
                 quality = fastq$X1[seq(from = 4, to = nrow(fastq), by = 4)])


fast1 <- new_fq %>%
  filter(id %in% locations$seqID) %>%
  left_join(locations, by = c("id" = "seqID")) %>%
  mutate(match_seq = str_sub(sequence, start, end)) %>%
  group_by(id, sequence, quality) %>%
  summarise(regex = paste0(match_seq, collapse = "|")) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(splitted = str_split(sequence, regex)) %>%
  ungroup() %>%
  unnest(splitted) %>%
  mutate(loc = str_locate(sequence, splitted),
         splitted_quality = str_sub(quality, loc[,"start"], loc[,"end"])) %>%
  group_by(id) %>%
  mutate(old_id = paste0("@", id),
         x = LETTERS[1:n()]) %>%
  ungroup() %>%
  mutate(sequence = splitted,
         quality = splitted_quality) %>%
  dplyr::select(sequence, quality, old_id) 

fast <- fast1 %>%
  rbind(dplyr::select(mutate(filter(new_fq, !(id %in% fast1$old_id)), old_id = paste0("@",id)), sequence, quality, old_id)) %>%
          left_join(fastq, by = c("old_id" = "X1")) %>%
  mutate(id  = paste0("@",1:nrow(.))) 
        
nfq <- rep("", times = nrow(fast) * 4)
nfq[seq(1, length(nfq)-3, by = 4)] <- fast$id
nfq[seq(2, length(nfq)-2, by = 4)] <- fast$sequence
nfq[seq(3, length(nfq)-1, by = 4)] <- "+"
nfq[seq(4, length(nfq), by = 4)] <- fast$quality


try <- tibble(X1 = nfq) %>%
  left_join(dplyr::select(fast, id, X2, X3, X4, X5, X6), by = c("X1" = "id")) 

write.table(try, output, quote = FALSE, col.names = FALSE, row.names = FALSE, na = "")
