# Natural isotope correction

All atoms are made of electrons, protons, and neutrons. An atom's elemental identity is determined by the number of protons, but the number of neutrons can vary. The variance in neutron number causes a slight shift in mass which can be leveraged for scientific purposes.

For example, in nature, about 1% of all carbon is "heavy carbon", with 7 neutrons instead of 6 neutrons. So, we can say there is a certain "natural abundance" of 13-C. When we experimentally modify the number of "heavy" (labelled) atoms in a molecule, we have to take into account the natural abundance of atoms that will LOOK like we labelled them, but which actually are just heavy by chance.

Therefore, I've built this simple natural isotope abundance correction script. The script requires a file called "isotope_data.csv" to be in the same directory. The "isotope_data.csv" file should have information on the natural abundances of M+0, M+1, M+2, &c isotopes of each atom in the molecule under inquiry.

The "samplefile.csv" shows an example input file. Input files should be given as such:

name_of_expt | chem_formula | M+0 | M+1 | M+2 |
-------------|--------------|-----|-----|-----|
myexpt | C10H17N2O14P3 | NF | 305 | 281 |

Peak intensities are given when they exist, and when the peak intensity is not found, the data point should be given as "NF". You can put multiple metabolites in an input file and adjust each one by the relative amounts of each element. The algorithm is "smart" and can read conventionally written chemical formulae. It is also unnecessary to have an "unlabelled" sample -- this simple script takes a purely analytical approach as given in "The importance of accurately correcting for the natural abundance of stable isotopes," (Midani, Wynn, Schnell, Analytical Biochemistry 2017, https://www.sciencedirect.com/science/article/pii/S0003269716304250), with some algorithmic approaches as in "Numerical bias estimation for mass spectrometric mass isotopomer analysis" (Yang,..., Heinzle, Analytical Biochemistry 2009, https://www.sciencedirect.com/science/article/pii/S0003269709001638)

Some notes:

* You will need to install R, an open-source, free programming language, in order to use this script. 

* If there are any elements in your compound which are *not* in the isotope_data file, they will need to be added before the algorithm can handle them.

* You can find out how to use the script (input and output options, as well as the required "-L" flag for which element type you have labelled!) by typing `Rscript natural_isotope_correction_v2.R --help` at the command line. 

* In this script there is an inherent assumption of complete isotopic purity of the labelling reagent. Obviously this is rarely/never true... 
