# Read tsv file and preprocess (1. convert to integer values, 2. process column names, 
# 3. process rownames) it 

read.tsv <- function(filepath) {
  data <- read.csv(filepath, header = TRUE, sep= "", quote = "")
  data <- as.matrix(data.matrix(data))
  rownames(data) <- gsub('"', '', rownames(data))
  colnames(data) <- gsub('X[.]+(.*)...', '\\1', colnames(data))
  colnames(data) <- gsub('[.]','-',colnames(data))
  return(data)  
}
