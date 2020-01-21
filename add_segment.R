# add_segment adds a segment to the data.frame genome, or creates it if
# it doesn't exist
# For the purposes of displaying, value is a pseudodensity variable.
# For inheritance, value indicates which parent contributed that
# segment to the root node person

add_segment <- function(genome, chromosome, seg_start, seg_end, value) {
  if(is.null(genome)) {
    genome <- tibble(Chr = as.character(chromosome),
                     Start = as.integer(seg_start),
                     End = as.integer(seg_end),
                     Value = as.integer(value))
  } else {
    genome <- genome %>% add_row(Chr = as.character(chromosome),
                                 Start = as.integer(seg_start),
                                 End = as.integer(seg_end),
                                 Value = as.integer(value))
  }
  return(genome)
}