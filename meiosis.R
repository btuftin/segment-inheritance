# The meiosis function takes a data frame with inherited segments,
# formatted with columns "Chr", "Start", "End", "Value",
# representing either the paternal or maternal chromosome,
# runs through it, one chromosome at a time, and constructs the
# resulting chromosome that is passed on, with segments
# passing on in full or in part or not at all due to 
# random choice of maternal or paternal chromosome,
# as well as crossover events.

# Probability of crossover event per 1Mbase
p <- 0.0125

meiosis <- function(genome) {
  # We'll work on one chromosome at a time
  out_genome <- NULL
  chromosomes <- unique(genome$Chr)
  for (i in 1:length(chromosomes)) {
    # randomly determine whether we're inheriting from 
    # this chromosome to begin with
    inheriting <- sample(c(TRUE, FALSE), 1)
    
    # starting with the first segment, determine whether it is inherited
    # in whole or in part
    segments <- genome %>% filter(Chr == chromosomes[i])
    
    # Check if there are any segments on the current chromosome
    if(nrow(segments) == 0) {next}
    
    # if there are more than one segment we will need to determine if there
    # are crossovers inbetween segments in addition to inside, but we can ignore any
    # "lead in"
    current_chunk_start <- segments[[1,2]]
    current_chunk_end <- current_chunk_start + 1000000 - 1
    
    # initialize the crossover point, it needs to persist between
    # iterations to preserve inheritance of close together segments
    old_crossover <- 0

    for (j in 1:nrow(segments)) {
      # Set segment information variables
      upstream_segment_start <- segments[[j,2]]
      upstream_segment_end <- segments[[j,3]]
      downstream_segment_start <- segments[[j,2]]
      
      # determine whether the is "blank space" before this segment
      # and if a crossover occurs there
      while (current_chunk_end < segments[[j,2]]){
        if (sample(c(TRUE, FALSE), 1, p = c(p, 1-p))) {
          # Switch to other chromosome in pair
          inheriting <- !inheriting
        }
        current_chunk_start <- current_chunk_end + 1
        current_chunk_end <- current_chunk_end + 1000000
        old_crossover <- 0
      }
      
      # Our chunk is now guaranteed to include the start and
      # potentially all of a segment
      upstream_segment_start <- segments[[j,2]]
      upstream_segment_end <- segments[[j,3]]
      downstream_segment_start <- segments[[j,2]]
      while(downstream_segment_start < upstream_segment_end) {
        # Determine if there is a crossover in this chunk with probability 1%
        # or if we inherited a crossover from the previous chunk due to
        # close contact segments
        if (sample(c(TRUE, FALSE),1,prob = c(p, 1-p))
            || old_crossover != 0 ) {
          
          # if we didn't inherit a crossover from the previous loop,
          # put the crossover somewhere in the 1 million long chunk
          if (old_crossover == 0) {
            crossover <- sample(current_chunk_start:current_chunk_end, 1)
          }
          
          # if the crossover is beyond the end of the current upstream segment,
          # set the end of the segment to pass on as the current end of the segment
          # and set old_crossover to crossover to pass on to check on the next
          # segment. (This is perhaps unnecessarily exact, but ...)
          if (crossover >= upstream_segment_end - 1) {
            downstream_segment_end <- upstream_segment_end
            old_crossover <- crossover
          } else {
            downstream_segment_end <- crossover
            old_crossover <- 0
          }
          
          # add the segment that ends at this crossover if we're on that part
          # of the pair
          if (inheriting) {
            out_genome <- add_segment(out_genome, 
                                      chromosomes[i], 
                                      downstream_segment_start, 
                                      downstream_segment_end,
                                      segments[[j,4]])
          }
          # change downstream_segment information
          downstream_segment_start <- downstream_segment_end +1
          inheriting <- !inheriting
        } else {
          # if there wasn't a crossover in this chunk
          # check if we have reached the end of this segment
          if (current_chunk_end > upstream_segment_end) {
            # if we've reached the end, add the segment we're currently on
            # if we're on the right chromosome in the pair
            if (inheriting) {
              out_genome <- add_segment(out_genome,
                                    chromosomes[i], 
                                    downstream_segment_start, 
                                    upstream_segment_end, 
                                    segments[[j,4]])
            }
            # set loop termination condition
            downstream_segment_start <- upstream_segment_end
            old_crossover <- 0
          }
        }
        # Unless we're passing on a crossover point
        # move along one chunk
        if (old_crossover == 0) {
          current_chunk_start <- current_chunk_end + 1
          current_chunk_end <- current_chunk_end + 1000000
        }
      }
    }
  }
  return(out_genome)
}