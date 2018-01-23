# This script takes input of a CSV file with the following format:
#
# datafile   M+0    M+1    M+2 ...
# sample1    38     395    NF
# sample2    NF     403    200
# ...
#
# and some arguments about the number of H, C, N, O atoms
# 
# and as output produces corrected isotope peak labels in a new csv file
# with the same labelling as the original.
#
# This program makes some assumptions. Mono-isotopic phosphorus,
# complete isotopic purity of the labelling reagent, and others.
#
# Largely, this paper was used:
# https://www.sciencedirect.com/science/article/pii/S0003269716304250
#
# This paper was also used to help formalise some of the math better,
# and for algorithm implementation tips:
# https://www.sciencedirect.com/science/article/pii/S0003269709001638

if(!require("nnls",character.only=T, quietly=T)) {
  install.packages("nnls")
  library("nnls", character.only=T)
}
if(!require("optparse",character.only=T, quietly=T)) {
  install.packages("optparse")
  library("optparse", character.only=T)
}

option_list <- list (
  make_option(c("--input_file","-i"), type="character",default=NULL,
              help="Name of input csv file (see header for format details)",
              metavar="character"),
  make_option(c("--Label","-L"),type="character",default="H",
              help="Atom type which has been labelled (e.g.: H, C)",
              metavar="character"),
  make_option(c("--output_file","-o"),type="character",default=NULL,
              help="Desired name of output file",metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##

input <- read.table(opt$input_file,header=T,sep=",",
                    stringsAsFactors=F,row.names=1)

# sanitise by making all NF --> 0 then putting to numeric
print("Treating all NF as 0...")
input[which(input == "NF",arr.ind=T)] <- 0

# create output matrix with all -1s so it will be obvious
# if something goes wrong
output <- input
output[,-c(1)] <- -1

isotope_data <- read.table("./isotope_data.csv",header=T,sep=",",
                           stringsAsFactors=F,row.names=1)

atomise <- function(chem_formula) {
  init <- strsplit(chem_formula, "(?=[A-Z])(?<=[0-9])|(?=[A-Z])(?<=[A-Za-z])",
                   fixed=F,perl=T)[[1]]
  # if it did not come with a number to signify how many of that atom there were,
  # there is one of that atom
  init[!grepl("\\d",init)] <- paste0(init[!grepl("\\d",init)], "1")
  split.up <- sapply(init, function(x) strsplit(x, "(?=[0-9])(?<=[A-Za-z])",fixed=F,perl=T)[[1]] )
  amounts <- as.integer(split.up[2,])
  names(amounts) <- split.up[1,]
  return(amounts)
}

for(chem_formula in unique(input$formula)) {
  print("Note that all correction matrices are analytical rather than empirical.")
  print(paste0("Building correction matrices for species of formula: ",chem_formula,"..."))
  
  ## split apart chemical formula into component parts
  atom_list <- atomise(chem_formula)
  
  ## build correction matrices
  ## correction matrices for assumed monoisotopic elements are identity matrices
  ## and therefore do not need to be built
  masterCM <- matrix(rep(0, (atom_list[opt$Label]+1 )**2),nrow=atom_list[opt$Label]+1)
  non_labelled_isos <- atom_list[-which(names(atom_list) == opt$Label)]
  
  ## build rho with a series of convolutions
  rho = 1
  for(atom_type in names(non_labelled_isos)) {
    a_atom <- isotope_data[atom_type,]
    if(any(a_atom == 0)) {
      a_atom <- as.numeric(a_atom[-which(a_atom == 0)])
    } else {
      a_atom <- as.numeric(a_atom)
    }
    for(b in 1:non_labelled_isos[atom_type]) {
      rho = convolve(rho, rev(a_atom),type="o")
    }
  }
  # any rho < 0 is an underflow error
  if(any(rho < 0)) {
    rho[which(rho < 0)] = 0
  }
  
  # the correction matrix just uses elements of rho to fill the columns
  for(col in 1:(atom_list[opt$Label]+1)) {
    masterCM[col:ncol(masterCM),col] <- rho[1:(ncol(masterCM)-col+1)]
  }
  
  # generate the labelled matrix, then multiply by the masterCM
  # to get the final master CM which we can then solve by nnls
  labelled_cm <- matrix(rep(0, (atom_list[opt$Label]+1 )**2),nrow=atom_list[opt$Label]+1)
  labelled_a <- isotope_data[opt$Label,]
  if(any(labelled_a == 0)) {
    labelled_a <- as.numeric(labelled_a[-which(labelled_a == 0)])
  } else {
    labelled_a <- as.numeric(labelled_a)
  }
  start_populating = 1
  for(n in 0:(ncol(labelled_cm)-1)) {
    N = atom_list[opt$Label]
    for(i in 0:(N - n)) {
      labelled_cm[i+start_populating,n+1] <- choose(N-n, i)*(labelled_a[1]**(N-n-i))*(labelled_a[2]**(i))
    }
    start_populating = start_populating + 1
  }
  
  masterCM <- masterCM %*% labelled_cm
  
  
  ### have master thing
  for(expt in rownames(input)) {
    if(input[expt,"formula"] == chem_formula) {
      unpadded_x <- as.numeric(input[expt,2:ncol(input)])
      solution <- nnls(masterCM, c(unpadded_x,rep(0, ncol(masterCM) - length(unpadded_x))))
      output[expt,2:ncol(output)] <- solution$x[1:(ncol(output) - 1)]
    }
  }
  
}

output[,-1] <- round(output[,-1],0)
if(is.null(opt$output_file)) {
  print(paste0("Writing to output file ","corrected_",opt$input_file))
  write.table(output, file=paste0("corrected_",opt$input_file),
              sep=",",quote=F)
} else {
  print(paste0("Writing to output file ", opt$output_file))
  write.table(output, file=opt$output_file,
              sep=",",quote=F)
}


