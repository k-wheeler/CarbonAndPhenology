##' @name read_Licor
##' @title read_Licor
##' 
##' @author Mike Dietze
##' @export
##' 
##' @param filename  name of the file to read
##' @param sep       file delimiter. defaults to tab
##' @param ...       optional arguements forwarded to read.table
read_Licor <- function(filename, sep = "\t", ...) {
  fbase <- sub(".txt", "", tail(unlist(strsplit(filename, "/")), n = 1))
  print(fbase)
  full <- readLines(filename)
  ## remove meta=data
  start <- grep(pattern = "OPEN", full)
  skip <- grep(pattern = "STARTOFDATA", full)
  for (i in rev(seq_along(start))) {
    full <- full[-(start[i]:(skip[i] + 1 * (i > 1)))]  # +1 is to deal with second header
  }
  full <- full[grep("\t", full)]  ## skip timestamp lines
  dat <- read.table(textConnection(full), header = TRUE, blank.lines.skip = TRUE,
                    sep = sep, ...)
  fname <- rep(fbase, nrow(dat))
  dat <- as.data.frame(cbind(fname, dat))
  return(dat)
} # read_Licor


mat2mcmc.list <- function(w) {
  temp <- list()
  chain.col <- which(colnames(w) == "CHAIN")
  for (i in unique(w[, "CHAIN"])) {
    temp[[i]] <- as.mcmc(w[w[, "CHAIN"] == i, -chain.col])
  }
  return(as.mcmc.list(temp))
}