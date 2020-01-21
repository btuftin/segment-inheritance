# The first_gen_genome function creates a data fram of a genome
# formatted with columns "Chr", "Start", "End", "Value",
# where value indicates the paternal or maternal chromosome,
# by randomly chosing which chromosome to pick out
# and randomly inserting crossover events under the assumption that
# they occur with a frequency of 1% per megabase

first_gen_genome <- function(){ 
  genome <- NULL
  # We'll ignore X and Y chromosomes because of their peculiar inheritance patterns
  for (i in 1:22) {
    # get endpoint for chromosome
    chr_end <- human_karyotype[i,3]
  
    # Well call the two chromosomes in a pair 1 and -1
    current_chr <- sample(c(-1, 1), 1)
    
    # To simplify the probabilities and dealing with the 
    # potential for multiple crossovers we approximate
    # crossover potential to 1% per megabase
    # and work with chunks of no more than one megabase
    current_chunk_start <- 1
    current_chunk_end <- 1000000
    
    current_segment_start <- 1
    while(current_segment_start < chr_end)
    {
      # Determine if there is a crossover with probability 1%
      # We ignore the small chance of two crossovers within one megabase
      if (sample(c(TRUE, FALSE),1,prob = c(0.01, 0.99))) {
          
          # put the crossover somewhere in the 1 million long block
          crossover <- sample(current_chunk_start:current_chunk_end, 1)
        
          # make sure we don't create a segment beyond the end of the chromosome
          if (crossover >= chr_end - 1) crossover <- chr_end
          # add the segment that ends at this crossover
          genome <- add_segment(genome, i, current_segment_start, 
                                crossover, current_chr)
          # change current_segment information
          current_segment_start <- crossover +1
          current_chr <- -current_chr
      } else {
        # if there wasn't a crossover in this chunk
        # check if we have reached the end of the chromosome
        if (current_chunk_end > chr_end) {
          # if we've reached the end, add the segment we're currently on
          genome <- add_segment(genome, i, current_segment_start, 
                                chr_end, current_chr)
          # set loop termination condition
          current_segment_start <- chr_end
        }
      }
      current_chunk_start <- current_chunk_end + 1
      current_chunk_end <- current_chunk_end + 1000000
    }
  }
  return(genome)
}

