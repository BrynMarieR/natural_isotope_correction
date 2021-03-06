# This script takes input of a CSV file with the following format:
#
# datafile   chemformula     M+0    M+1    M+2 ...
# sample1    C6H12O6         38     395    NF
# sample2    C6H12O6         NF     403    200
# ...
#
# and the type of atom that was labelled in the experiment (e.g., -L H)
# 
# and as output produces corrected isotope peak labels in a new csv file
# with the same labelling as the original.
#
# This program makes some assumptions. Mono-isotopic phosphorus,
# and others
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
if(!require("gtools",character.only = T,quietly=T)) {
  install.packages("gtools")
  library(gtools,character.only=T)
}
if(!require("plyr",character.only = T,quietly=T)) {
  install.packages("plyr")
  library(plyr,character.only=T)
}

option_list <- list (
  make_option(c("--input_file","-i"), type="character",default=NULL,
              help="Name of input csv file (see header for format details)",
              metavar="character"),
  make_option(c("--Label","-L"),type="character",default="H",
              help="Atom type which has been labelled (e.g.: H, C)",
              metavar="character"),
  make_option(c("--output_file","-o"),type="character",default=NULL,
              help="Desired name of output file",metavar="character"),
  make_option(c("--method","-m"),type="character",default="skewed",
              help="Desired correction method (skewed, classical -- CLASSICAL NOT RECOMMENDED)"),
  make_option(c("--purity","-p"),type="double",default=1.00,
              help="Tracer purity, which should be <= 1.00")
)

opt_parser <- OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##

method = opt$method
purity = as.numeric(opt$purity)
input <- read.table(opt$input_file,header=T,sep=",",
                    stringsAsFactors=F)
#barplot(as.matrix(input[,3:6]),beside=T)


# sanitise by making all NF --> 0 then putting to numeric
print("Treating all NF as 0...")
input[which(input == "NF",arr.ind=T)] <- 0

# create output matrix with all -1s so it will be obvious
# if something goes wrong
output <- input
output[,-c(1,2)] <- -1

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
  masterCM <- diag(atom_list[opt$Label]+1)
  non_labelled_isos <- atom_list[-which(names(atom_list) == opt$Label)]
  
  ## build rho with a series of convolutions
  for(atom_type in names(non_labelled_isos)) {
    rho = 1
    a_atom <- isotope_data[atom_type,]
    if(any(a_atom == 0)) {
      a_atom <- as.numeric(a_atom[-which(a_atom == 0)])
    } else {
      a_atom <- as.numeric(a_atom)
    }
    for(b in 1:non_labelled_isos[atom_type]) {
      rho = convolve(rho, rev(a_atom),type="o")
    }
    # any rho < 0 is an underflow error
    if(any(rho < 0)) {
      rho[which(rho < 0)] = 0
    }
    local_cm <- diag(atom_list[opt$Label] + 1)
    if(length(rho) < ncol(local_cm)) { # pad with zeros
      rho <- c(rho, rep(0, (ncol(local_cm) - length(rho))))  
    }
    for(col in 1:(atom_list[opt$Label] + 1)) {
      local_cm[col:ncol(local_cm), col] <- rho[1:(ncol(local_cm) - col + 1)]
    }
    masterCM <- masterCM %*% local_cm
  }
  
  # generate the labelled matrix, then multiply by the masterCM
  # to get the final master CM which we can then solve by nnls
  labelled_cm <- matrix(rep(0, (atom_list[opt$Label]+1 )**2),nrow=atom_list[opt$Label]+1)
  labelled_a <- isotope_data[opt$Label,]
  
  if(method=="skewed") {

    start_populating = 1
    possible_isotopes = seq(0, ncol(isotope_data) - 1)
    for(n_from_labelled in 0:(ncol(labelled_cm)-1)) {
      N_atoms = atom_list[opt$Label]
      
      for(i in 0:(N_atoms - n_from_labelled)) {
        if(N_atoms != n_from_labelled) {
          result <- combinations(n=ncol(isotope_data), 
                                 r=(N_atoms - n_from_labelled), v=seq(0,N_atoms), repeats.allowed=T)
          result <- result[which(apply(result, 1, sum) == i),]
          if(!is.vector(result)) {
            probability <- sum(apply(result, 1, function(x) {
              tmp <- count(x)
              if(nrow(tmp) < length(possible_isotopes)) {
                missing <- data.frame(cbind(possible_isotopes[-which(possible_isotopes %in% tmp$x)], c(0)))
                colnames(missing) <- c("x","freq")
                tmp <- rbind(tmp, missing)
              }
              sorted.tmp <- tmp[order(tmp$x),] 
              dmultinom(x=sorted.tmp[,"freq"], prob=as.matrix(isotope_data[opt$Label,]))
            }))
          } else {
            tmp <- count(result)
            if(nrow(tmp) < length(possible_isotopes)) {
              missing <- data.frame(cbind(possible_isotopes[-which(possible_isotopes %in% tmp$x)], c(0)))
              colnames(missing) <- c("x","freq")
              tmp <- rbind(tmp, missing)
            }
            sorted.tmp <- tmp[order(tmp$x),] 
            probability <- dmultinom(x=sorted.tmp[,"freq"], prob=as.matrix(isotope_data[opt$Label,]))
          }
        } else {
          probability = 1
        }

        labelled_cm[i+start_populating,n_from_labelled+1] <- probability
      }
      start_populating = start_populating + 1
    }
  } else if (method == "classical") {
    ## build rho with a series of convolutions
    rho = 1

    for(b in 1:atom_list[opt$Label]) {
      rho = convolve(rho, rev(as.matrix(labelled_a)), type="o")
    }

    # any rho < 1e-15 is an underflow error
    if(any(rho < 1e-15)) {
      rho[which(rho < 1e-15)] = 0
    }
    
    # the correction matrix just uses elements of rho to fill the columns
    for(col in 1:(atom_list[opt$Label]+1)) {
      labelled_cm[col:ncol(labelled_cm),col] <- rho[1:(ncol(labelled_cm)-col+1)]
    }    
  }
  
  # correct for purity of the tracer
  for(col in 2:(atom_list[opt$Label] + 1)) {
    tmp_rho <- labelled_cm[col:ncol(labelled_cm),col]
    for (convolve_timer in 1:(col-1)) {
      tmp_rho <- convolve(tmp_rho, rev(c(1-purity, purity)),type="o")
    }
    labelled_cm[1:ncol(labelled_cm),col] <- tmp_rho[1:ncol(labelled_cm)]
  }
  
  masterCM <- masterCM %*% labelled_cm
  
  
  ### have master thing
  for(expt in rownames(input)) {
    if(input[expt,"formula"] == chem_formula) {
      unpadded_x <- as.numeric(input[expt,3:ncol(input)])
      solution <- nnls(masterCM, c(unpadded_x,rep(0, ncol(masterCM) - length(unpadded_x))))
      # set output after normalising to make the M+0 peak 1
      output[expt,3:ncol(output)] <- solution$x[1:(ncol(output) - 2)] / solution$x[1]
    }
  }
  
}

output[,-c(1,2)] <- round(output[,-c(1,2)],3)
if(is.null(opt$output_file)) {
  print(paste0("Writing to output file ","corrected_",opt$input_file))
  write.table(output, file=paste0("corrected_",opt$input_file),
              sep=",",quote=F,row.names=F)
} else {
  print(paste0("Writing to output file ", opt$output_file))
  write.table(output, file=opt$output_file,
              sep=",",quote=F,row.names=F)
}




