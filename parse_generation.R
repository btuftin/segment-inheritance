# The parse_generation function accepts as input a data frame with the format
# iD, genome, where iD is a string identifying the individual (eg. 3.4 for the
# fourth child of the third child) and genome is a data frame of the format
# described in first_generation.R and meiosis.R  (columns "Chr", 
# "Start", "End", "Value")

# The segment data is compared and used to create statistics and a plottable
# representation of the original genome. Whereas the input genomes have
# up to 22 chromosomes and values -1 and 1 representing the two cromosomes
# of the 0th generation, the output genome has a up to 44 chromosomes
# with all positive values representing the number of individuals
# preserving that segment from that original chromosome

# The longest chromosome has a length in the plot of about 1600 pixels
# and is 248,956,422 bases long, or about 155,000 bases per pixel.
# Each chromosome will be represented in comparison by creating a
# raw vector of length chromosome_length / 78000 which allocates abt two
# values in the vector for each pixel
divisor <- 78000

parse_generation <- function(generation) {
  # initialize return object
  out <- tibble(id = character(),
                object = list())
  
  # Initialize data frame for the results of the first parsing
  vector_frame <- tibble(Chr = character(), origin = character(), vector = list())
  
  # Initialize size tibble, needed for statistics
  chr_sizes_scaled <- integer(length = nrow(human_karyotype)-2)
  for (i in 1:length(chr_sizes_scaled)) {
    chr_sizes_scaled[i] <- round(human_karyotype$End[i]/divisor)
  }
  
  # Initialize unscaled size tibble, needed for final output
  # Could use human_karyotype directly, but doing this for flexibility
  chr_sizes <- integer(length = nrow(human_karyotype)-2)
  for (i in 1:length(chr_sizes_scaled)) {
    chr_sizes[i] <- human_karyotype$End[i]
  }
  
  
  # Iterate through the individuals
  for (i in 1:nrow(generation)) {
    
    # Iterate through the chromosomes
    chromosomes <- unique(generation[[i,2]]$Chr)
    for (j in 1:length(chromosomes)) {
      ChrId <- chromosomes[j]
      
      # iterate through the segments on this chromosome
      segments <- generation[[i,2]] %>% filter(Chr == ChrId)
      
      # If there are no segments for this chromosome skip to the next
      if(nrow(segments) == 0) next
      
      # Create raw vectors to take the values, automatically initialized
      # to zero
      chr_length <- human_karyotype %>% filter(Chr == ChrId) %>% .$End
      vector_length <- round(chr_length / divisor)
      vector_a <- raw(length = vector_length)
      vector_b <- raw(length = vector_length)
      
      for (k in 1:nrow(segments)) {
        segment <- segments[k,] 
        segment_start <- round(segment$Start/divisor)
        segment_end <- round(segment$End/divisor)
        
        # double check for rounding errors
        if (segment_end > vector_length) {segment_end <- vector_length}
        
        # check which chromosome to mark and mark it
        if(segment$Value == 1) {
          vector_a[segment_start:segment_end] <- as.raw(1)
        } else {
          vector_b[segment_start:segment_end] <- as.raw(1)
        }
      }
      # add the vectors to the data frame
      vector_frame <- vector_frame %>%
        add_row(Chr = ChrId, origin = "a", vector = list(vector_a)) %>%
        add_row(Chr = ChrId, origin = "b", vector = list(vector_b))
    }
  }
  
  # Figured out how to initialize empty tibble 
  # The first row in our data frame was a dummy row, so we delete that now
  #  vector_frame <- vector_frame[-1,]
  
  # Now we have a data frame with vectors showing the chromosomes
  # from the 0th generation still present in each invididual.
  # Next step is to parse these to get just one vector per unique Chr and origin value
  # and with weights representing the number of individuals matching each
  # other on that segment
  vector_frame_collapsed <- tibble(Chr = character(), origin = character(), vector = list())
  
  # Iterate through the chromosomes !Note! Id is spread across Chr and origin!
  chromosomes <- unique(vector_frame$Chr)
  for (i in 1:length(chromosomes)) {
    ChrId <- chromosomes[i]
    for (j in c("a", "b")) {
      vf_subset <- vector_frame %>% filter(Chr == ChrId, origin == j)
      if (nrow(vf_subset) == 0) next
      
      chr_vector <- integer(length = length(vf_subset$vector[[1]]))
      # If there is just one such chromosome, add that, else combine them
      if (nrow(vf_subset) == 1) {
        vector_frame_collapsed <- vector_frame_collapsed %>%
          add_row(Chr = ChrId,
                  origin = vf_subset$origin,
                  vector = list(as.integer(vf_subset$vector[[1]])))
      } else {
        for (i in 1:nrow(vf_subset)) {
          chr_vector <- chr_vector + as.integer(vf_subset$vector[[i]])
        }
        vector_frame_collapsed <- vector_frame_collapsed %>%
          add_row(Chr = ChrId,
                  origin = j,
                  vector = list(chr_vector))
      }
    }
  }
  
  
  # Figured out how to initialize empty tibble 
  # The first row in our data frame was a dummy row, so we delete that now
  #  vector_frame_collapsed <- vector_frame_collapsed[-1,]
  
  # We can now do statistics on current coverage
  # start by initializing object and iterating through chromosomes
  coverage_pr_chromosome <- tibble(Chr = character(), coverage = integer())
  
  total_present <- 0
  
  for(i in 1:length(chr_sizes_scaled)) {
    for (o in c("a", "b")) {
      chromosome <- vector_frame_collapsed %>%
        filter(Chr == as.character(i), origin == o)
      
      # If there is data on that chromosome, add that data to the output
      # and add it to the total. If not, add a line with 0 to the output.
      if (nrow(chromosome) != 0) {
        coverage_pr_chromosome <- coverage_pr_chromosome %>%
          add_row(Chr = paste(as.character(i), o),
                  coverage = mean(as.logical(chromosome$vector[[1]])))
        total_present <- total_present + sum(as.logical(chromosome$vector[[1]]))
      } else {
        coverage_pr_chromosome <- coverage_pr_chromosome %>%
          add_row(Chr = paste(as.character(i), o),
                  coverage = 0)
      }
    }
  }
  
  out <- out %>% add_row(id = "Coverage pr chr", object = list(coverage_pr_chromosome))
  
  # All chromosomes complete we take the sum of present blocks and
  # divide by total blocks
  out <- out %>% add_row(id = "Total coverage",
                         object = total_present/(2*sum(chr_sizes_scaled)))
  
  
  # In preparation to create the object that lets us plot the coverage we
  # now need collapse all the pairs down to one, and change the sign on
  # the weight for segments representing two different segments
  vector_frame_final <- tibble(Chr = character(), vector = list())
  
  # Iterate through the chromosomes
  chromosomes <- unique(vector_frame_collapsed$Chr)
  for (i in 1:length(chromosomes)) {
    ChrId <- chromosomes[i]
    vf_subset <- vector_frame_collapsed %>% filter(Chr == ChrId)
    
    # If there is just one such chromosome, add that, else combine them
    if (nrow(vf_subset) == 1) {
      vector_frame_final <- vector_frame_final %>%
        add_row(Chr = ChrId, vector = list(vf_subset$vector[[1]]))
    } else {
      # Adding the two vectors sets the weight
      vector_out <- as.integer(vf_subset$vector[[1]]) + 
        as.integer(vf_subset$vector[[2]])
      
      # Casting them to logical and multiplying creates
      # a logical vector showing where they overlap
      overlap <- as.logical(vf_subset$vector[[1]]) * as.logical(vf_subset$vector[[2]])
      
      # Subtracting the vector * 2 * overlap from itself, flips the sign where
      # there are overlaps
      vector_out <- vector_out - 2*vector_out*overlap
      vector_frame_final <- vector_frame_final %>%
        add_row(Chr = ChrId, vector = list(vector_out))
    }
  }
  
  
  
  # To get back to a representation of the chromosome with a segment pr row
  # and start and endpoints marked as bases we need to parse through our vector
  # frame and note where the values flip
  
  out_genome <- tibble(Chr = character(), Start = integer(),
                       End = integer(), Value = integer())
  
  for(i in 1:length(chr_sizes_scaled)) {
    chromosome <- vector_frame_final %>%
      filter(Chr == as.character(i))
    
    # If there is no data on that chromosome, add a segment across all of it
    # with value 0. If there is, parse the data and add segments for each
    # section with the same weight and add the weight as value
    if (nrow(chromosome) == 0) {
      out_genome <- out_genome %>%
        add_row(Chr = as.character(i),
                Start = 1,
                End = chr_sizes_scaled[i],
                Value = 0)
    } else {
      
      start_point <- 1
      end_point <- 1
      value <- Inf
      current_chr <- chromosome$Chr
      current_vec <- chromosome$vector[[1]]
      
      # For each vector
      for (j in 1:length(current_vec)) {
        if (current_vec[j] != value) {
          # We're at the start of a chromosome or the start of a new
          # segment
          if (value == Inf) {
            value <- current_vec[j]
            end_point <- end_point + divisor
          } else {
            # We've reached a new segment and should add it and update
            # our running values
            out_genome <- out_genome %>%
              add_row(Chr = as.character(i),
                      Start = start_point,
                      End = end_point,
                      Value = value)
            value <- current_vec[j]
            start_point <- end_point + 1
            end_point <- end_point + divisor
          }
        } else {
          # We're on the same segment as the block before so add
          # the appropriate value to the end_point
          end_point <- end_point + divisor
        }
      }
      # We add the last segment
      out_genome <- out_genome %>%
        add_row(Chr = as.character(i),
                Start = start_point,
                End = chr_sizes[i],
                Value = value)
    }
  }
  # We add the plotable coverage object to the output and return from the function
  out <- out %>% add_row(id = "plot_genome", object = list(out_genome))
}