require("seqinr")

getncbiseq <- function(accession) {
  # first find which ACNUC database the accession is stored in: 
  dbs <- c("genbank","refseq","refseqViruses","bacterial")
  numdbs <- length(dbs)
  for (i in 1:numdbs) {
    db <- dbs[i]
    choosebank(db)
    # check if the sequence is in ACNUC database â€™dbâ€™:
    resquery <- try(query(".tmpquery", paste("AC=", accession)), silent = TRUE) 
    if (!(inherits(resquery, "try-error")))  {
      print(paste("trying: ",db))
      queryname <- "query2"
      thequery <- paste("AC=",accession,sep="") 
      print(thequery)
      # query("queryname","thequery")
      query2 <- query(`queryname`,`thequery`)
      # see if a sequence was retrieved:
      seq <- getSequence(query2$req[[1]]) 
      closebank()
      return(seq) 
    }
    print(paste("accession not in: ",db))
    closebank()
  }
  print(paste("ERROR: accession",accession,"was not found"))
}
retrieveseqs <- function(seqnames,acnucdb) {
  myseqs <- list()   # Make a list to store the sequences
  choosebank(acnucdb)
  for (i in 1:length(seqnames)) {
    seqname <- seqnames[i]
    print(paste("Retrieving sequence",seqname,"..."))
    queryname <- "query2"
    query <- paste("AC=",seqname,sep="")
    query2 <- query(`queryname`,`query`)
    seq <- getSequence(query2$req[[1]]) # Makes a vector "seq" containing the sequence
    myseqs[[i]] <- seq
  }
  closebank()
  return(myseqs)
}
cleanAlignment <- function(alignment, minpcnongap, minpcid)
{
  # make a copy of the alignment to store the new alignment in:
  newalignment <- alignment
  # find the number of sequences in the alignment
  numseqs <- alignment$nb
  # empty the alignment in "newalignment")
  for (j in 1:numseqs) { newalignment$seq[[j]] <- "" }
  # find the length of the alignment
  alignmentlen <- nchar(alignment$seq[[1]])
  # look at each column of the alignment in turn:
  for (i in 1:alignmentlen)
  {
    # see what percent of the letters in this column are non-gaps:
    nongap <- 0
    for (j in 1:numseqs)
    {
      seqj <- alignment$seq[[j]]
      letterij <- substr(seqj,i,i)
      if (letterij != "-") { nongap <- nongap + 1}
    }
    pcnongap <- (nongap*100)/numseqs
    # Only consider this column if at least minpcnongap % of the letters are not gaps:
    if (pcnongap >= minpcnongap)
    {
      # see what percent of the pairs of letters in this column are identical:
      numpairs <- 0; numid <- 0
      # find the letters in all of the sequences in this column:
      for (j in 1:(numseqs-1))
      {
        seqj <- alignment$seq[[j]]
        letterij <- substr(seqj,i,i)
        for (k in (j+1):numseqs)
        {
          seqk <- alignment$seq[[k]]
          letterkj <- substr(seqk,i,i)
          if (letterij != "-" && letterkj != "-")
          {
            numpairs <- numpairs + 1
            if (letterij == letterkj) { numid <- numid + 1}
          }
        }
      }
      pcid <- (numid*100)/(numpairs)
      # Only consider this column if at least %minpcid of the pairs of letters are identical:
      if (pcid >= minpcid)
      {
        for (j in 1:numseqs)
        {
          seqj <- alignment$seq[[j]]
          letterij <- substr(seqj,i,i)
          newalignmentj <- newalignment$seq[[j]]
          newalignmentj <- paste(newalignmentj,letterij,sep="")
          newalignment$seq[[j]] <- newalignmentj
        }
      }
    }
  }
  return(newalignment)
}

## get rid of bootstrap probabilities in a tree output
convert_tree <- function(tree) {
  gsub("\\)[0-9]\\.[0-9]+", ")", tree, perl=TRUE)
}
