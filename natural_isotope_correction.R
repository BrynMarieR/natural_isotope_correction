# This script takes input of a labelled isotope vector (deuterium-labelled)
# and "corrects" it for the natural isotopic abundance of 18 O, 13 C, 15 N.
# based on an algorithm from the following Rabinowitz paper:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909566/
# However, it's my opinion that this algorithm is WRONG
# because it doesn't account for the fact that an additional 18 O
# adds *2* to the expected mass, not 1. The following paper implements
# a more correct approach and I think a future attempt at this problem
# should use this instead:
# https://www.sciencedirect.com/science/article/pii/S0003269716304250

if(!require("nnls",character.only=T, quietly=T)) {
  install.packages("nnls")
  library(nnls, character.only=T)
}

cat("Provide the measured isotopic labeling distribution M, M+1,...,M+n (in format 123,456,...,789) : ")
m.list <- readLines(file("stdin"),n=1,warn=F)
closeAllConnections()
m.list <- unlist(strsplit(m.list, ",",fixed=T))
numeric.m.list <- unlist(lapply(m.list, function(x) {
  tmp <- trimws(x,"both")
  as.integer(tmp)
}))
print(numeric.m.list)

cat("How many hydrogen atoms in the molecule? ")
n_h <- as.integer(readLines(file("stdin"),n=1))
closeAllConnections()

if(n_h >= length(numeric.m.list)) {
  numeric.m.list <- c(numeric.m.list, rep(0, n_h - length(numeric.m.list)))
} else {
  cat("You cannot have fewer hydrogens than isotope masses!")
}

cat("How many carbon atoms in the molecule? ")
n_c <- as.integer(readLines(file("stdin"),n=1))
closeAllConnections()

cat("How many nitrogen atoms in the molecule? ")
n_n <- as.integer(readLines(file("stdin"),n=1))
closeAllConnections()

cat("How many oxygen atoms in the molecule? ")
n_o <- as.integer(readLines(file("stdin"),n=1))
closeAllConnections()

carbon_mat = nitrogen_mat = oxygen_mat = matrix(rep(0, n_h*n_h),nrow=n_h)
for(atom in c("carbon","nitrogen","oxygen")) {
  tmp_mat <- matrix(rep(0, n_h*n_h),nrow=n_h)
  if(atom == "carbon") {
    n_atom = n_c
    nat_abund = 0.0107
  } else if (atom == "nitrogen") {
    n_atom = n_n
    nat_abund = 0.00368
  } else if (atom == "oxygen") {
    n_atom = n_o
    nat_abund = 0.00205
  }
  
  # populate matrix
  i = 0
  for(j in 1:(n_h - 2)) {
    tmp_mat[j,j] <- choose(n_atom, i) * ((1 - nat_abund)**(n_atom - i)) * ((nat_abund)**(i))
    tmp_mat[j+1,j] <- choose(n_atom, i+1) * ((1 - nat_abund)**(n_atom - i+1)) * ((nat_abund)**(i+1))
    tmp_mat[j+2,j] <- choose(n_atom, i+2) * ((1 - nat_abund)**(n_atom - i+2)) * ((nat_abund)**(i+2))
  }
  # fill in last two columns, they are shorter
  tmp_mat[n_h - 1,n_h - 1] <- choose(n_atom, i) * ((1 - nat_abund)**(n_atom - i)) * ((nat_abund)**(i))
  tmp_mat[n_h - 1 + 1,n_h - 1] <- choose(n_atom, i+1) * ((1 - nat_abund)**(n_atom - i+1)) * ((nat_abund)**(i+1))
  tmp_mat[n_h,n_h] <- choose(n_atom, i) * ((1 - nat_abund)**(n_atom - i)) * ((nat_abund)**(i))
  
  if(atom == "carbon") {
    carbon_mat <- tmp_mat
  } else if (atom == "nitrogen") {
    nitrogen_mat <- tmp_mat
  } else if (atom == "oxygen") {
    oxygen_mat <- tmp_mat
  }
}

full_correction_mat <- carbon_mat %*% nitrogen_mat %*% oxygen_mat
solution <- nnls(full_correction_mat, numeric.m.list)
print("Original: ")
print(numeric.m.list)
print("Corrected: ")
print(solution$x)

