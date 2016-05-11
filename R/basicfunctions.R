
#' Corner:
#' Display corner of matrix or data.frame
#' @param data.frame or matrix to view
#' @param ncr number of columns/rows to view (default = 5)


corner <- function(data.frame, ncr = 5) data.frame[1:ncr, 1:ncr]

#' Save data.frame with write.table to file test.txt
#' @param df data.frame to save. Formats with quote = F, row.names = F and tab delimited.


save.test <- function(df) write.table(df, "test.txt", quote = F, sep = "\t", row.names = F)

#' Recode multiple values in a vector at once.
#' @param data vector to be recoded
#' @param oldvalue vector of old values
#' @param newvalue vector of new values corresponding to oldvalue



recoderFunc <- function(data, oldvalue, newvalue) {

  # convert any factors to characters

  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)

  # create the return vector

  newvec <- data

  # put recoded values into the correct position in the return vector

  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]

  newvec

}

#' Count instances of individual entries in vector, return object same length as
#' vector
#' @param x vector to count


countIF <- function(x){

  x <- data.frame(Var1 = x)
  x$Order <- 1:nrow(x)

  y <- data.frame(table(x$Var1))

  x <- merge(x, y, all.x = T)

  x <- x[with(x, order(Order)), ]

  x$Freq

}
