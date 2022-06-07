#CodonAlignment.R
#Alignment of fa CDS files. Usage (in command line): `Rscript CodonAlignment.R <myfasta.fa> [<other optional .fa's ... >]
suppressMessages(library(Biostrings))
suppressMessages(library(DECIPHER))
codon_align<-function(inputList){
  inputSeqs<-readDNAStringSet(inputList)#Vector of DNA FASTA files.
  inputSeqs<-RemoveGaps(inputSeqs)#Dealign any alignments for correct sequence length.
  
  inputSeqs<-inputSeqs[]#Filter by some criteria.
  names(inputSeqs)<-lapply(X = names(inputSeqs),function(X){paste(strsplit(X," ")[[1]],collapse= "_")}) #avoids some name parsing issues with other programs.
  
  MAlign<-AlignTranslation(processors = 8,inputSeqs,refinements=10)
  SMAlign<-StaggerAlignment(MAlign)
  return(SMAlign)
}

if (!interactive()){
  inputList<-commandArgs(trailingOnly=TRUE)  
  codon_align(inputList)#CDS fasta files. 
}
