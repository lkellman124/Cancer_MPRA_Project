library(motifbreakR)
library(universalmotif)
library(TFBSTools)
library(here)
library(magrittr)
library(tidyverse)
library(here)
library(furrr)


# By Robin Meyers (with slight adaptations)

out_dir <- file.path(here("output/motifs/"))
# dir.create(out_dir, recursive = T)

### Read in data
mpra_ref <- read_tsv(here("data/new_lib_5.20_jmprimers.txt"))
mpra_ref <- mpra_ref %>% dplyr::select(snp_id, locus, allele, sequence) %>%
  distinct()
mpra_ref <- mpra_ref %>% unite(col = "seq_id", locus, allele, sep = "_", remove=F )

hocomoco_pwms <- motifbreakR::hocomoco %>% convert_motifs("TFBSTools-PWMatrix") %>%
    set_names(map_chr(., ~ name(.))) %>% do.call(PWMatrixList, .)


# function to search all sequences for a single motif and return results in a tidy data frame
search_tf_motif <- function(pwm, seqs, min_score = 0.8) {
    # make sure you understand what's going on in here!
    cat(pwm@name, "\n")
    search_results <- searchSeq(pwm, DNAStringSet(seqs),
                                min.score = min_score)
    search_results_with_at_least_one_hit <- map_lgl(as.list(search_results), ~ length(.) > 0)
    search_results_cleaned <- as.list(search_results)[search_results_with_at_least_one_hit] %>%
        map_dfr(as.data.frame, .id = "sequence")
    return(search_results_cleaned)
}

# Get mpra sequences with their ID as names
mpra_sequences <- mpra_ref$sequence %>% set_names(mpra_ref$seq_id)


# you can manually set the number of cores you want to use
# this will take the maximum
# lower it if your computer starts running too slowly
cores <- availableCores() - 1

plan(multisession, workers = cores)

search_tf_motif_to_file <- function(pwm, seqs, file, min_score = 0.8) {
    # make sure you understand what's going on in here!
    cat(pwm@name, "\n")
    search_results <- searchSeq(pwm, DNAStringSet(seqs),
                                min.score = min_score)
    search_results_with_at_least_one_hit <- map_lgl(as.list(search_results), ~ length(.) > 0)
    search_results_cleaned <- as.list(search_results)[search_results_with_at_least_one_hit] %>%
        map_dfr(as.data.frame, .id = "sequence")
    write_tsv(search_results_cleaned, file)
}

motif_dir <- file.path(out_dir, "motif_matches")
dir.create(motif_dir, recursive = T)

# try on a limited number of motifs and sequences
future_walk(as.list(hocomoco_pwms[1:16]),
            ~ search_tf_motif_to_file(., mpra_sequences[1:100],
                                      file.path(motif_dir, paste0(.@ID, ".tsv")),
                                      min_score = 0.8))


# Now map this function across all TF motifs in HOCOMOCO while searching all MPRA sequences
future_walk(as.list(c(hocomoco_pwms)),
            ~ search_tf_motif_to_file(., mpra_sequences,
                                      file.path(motif_dir, paste0(.@ID, ".tsv")),
                                      min_score = 0.8))




