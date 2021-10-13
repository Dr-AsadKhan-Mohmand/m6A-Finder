1. m6.FASTA contains 1307 sequecnes of m6A sites and 1307 sequnces of non-m6A sities.
2. acf.py extracts auto-corelation features based on physical properties and exports the features in 
CSV format in ACF.CSV file.
3. m6A-Finder.R contains the following segments
i) m6A-Finder.R reads ACF.CSV.
ii) it extracts hexa-nulcleotide composition features.
iii) both types of features are combined.
iv) The mRMR is used to select optimal features.
v)  The SVM RBF is utilized for classification.