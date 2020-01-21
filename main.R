library(RIdeogram)
library(tidyverse)

# Load the human_karyotype data from RIdeogram to get start and
# end point for chromosomes
data("human_karyotype")

# Load functions
source("add_segment.R")
source("first_generation.R")
source("meiosis.R")
source("parse_generation.R")

children_per_generation <- 3 # for a static number use a positive number
                             # plan for later: use 0 for dynamic
# mean_offspring <- 3
# sd_offspring <- 2

# Create first generation genome
# first_gen <- first_gen_genome()

# create plot to see how reasonable it is
# ideogram(karyotype = human_karyotype, overlaid = first_gen)
# convertSVG("chromosome.svg", device = "png")
# A couple of things are obvious
# 1. It looks very similar to the real thing. Differences may be
# due to the 1% pr Mbase approximation being of varying accuracy
# from chromosome to chromosome or because I'm comparing with 
# only one grandparent mapping
# 2. It can create crossovers right around the centromere that
# are not realistic.

# Create several first generation genomes and determine how
# many percent of the 0th generation genome is present.
current_generation <- tibble(
  iD = c("1", "2", "3"),
  genome = list(first_gen_genome(), first_gen_genome(), first_gen_genome())
)

analysis <- parse_generation(current_generation)                                                     

# Create second generation genomes and do the same
next_generation <- tibble(
  iD = character(),
  genome = list()
)

for (i in 1:nrow(current_generation)) {
  next_generation <- next_generation %>%
    add_row(iD = c(paste0(current_generation$iD[i], 1),
                   paste0(current_generation$iD[i], 2),
                   paste0(current_generation$iD[i], 3)),
            genome = c(list(meiosis(current_generation$genome[[i]])),
                       list(meiosis(current_generation$genome[[i]])),
                       list(meiosis(current_generation$genome[[i]])))
    )
                   
}

analysis2 <- parse_generation(next_generation)

ideogram(karyotype = human_karyotype, overlaid = analysis2[[3,2]])
convertSVG("chromosome.svg", device = "png")
