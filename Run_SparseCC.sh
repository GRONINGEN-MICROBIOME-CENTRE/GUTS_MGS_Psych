##############################################
#1. Activate conda environment with SparseCC##
############################################## 
ENV=metaSNV
source activate $ENV
#################
#2. Run sparseCC#
#################
#Usage: Compute the correlation between components (e.g. OTUs). By default uses the SparCC algorithm to account for compositional effects. Correlation and covariance (when applies) matrices are written out as txt files. 
#Counts file needs to be a tab delimited text file where columns are samples and rows are components (e.g. OTUS).
#Usage:   SparCC.py counts_file [options]
#Example: SparCC.py example/fake_data.txt -i 20 --cor_file=example/basis_corr/cor_mat_sparcc.out

#Options:
#  -h, --help            show this help message and exit
#  -c COR_FILE, --cor_file=COR_FILE
#                        File to which correlation matrix will be written.
#  -v COV_FILE, --cov_file=COV_FILE
#                        File to which covariance matrix will be written.
#  -a ALGO, --algo=ALGO  Name of algorithm used to compute correlations (SparCC
#                        (default) | pearson | spearman | kendall)
#  -i ITER, --iter=ITER  Number of inference iterations to average over (20
#                        default).
#  -x XITER, --xiter=XITER
#                        Number of exclusion iterations to remove strongly
#                        correlated pairs (10 default).
#  -t TH, --thershold=TH
#                        Correlation strength exclusion threshold (0.1
#                        default).

#To run this script we first need to preprocess the data, see the main R script for that.
#Preprocessing consists of: 1. Removal of low prevalent taxa (0s are considered to be due to undersampling). 2. Formatting of data (transposition, first column taxa, first row Samples). 3. Dividing taxa from Cases and Controls

#Note: I changed core_methods.py so that in the normalize function, it uses np.apply_along_axis independently of the class of the object. Otherwise the output was a (x,) series object that would fail later on because it had no dimensions in its column 

echo "Running SparCC in cases"
#Run in cases
SparCC.py Data/SparseCC_input/Cases.tsv -i 20 -t 0.1 --cor_file=Results/SparseCC/Cases/Corr_matrix_cases.txt
echo "Running SparCC in controls"
#Run in controls
SparCC.py Data/SparseCC_input/Controls.tsv -i 20 -t 0.1 --cor_file=Results/SparseCC/Controls/Corr_matrix_controls.txt

#Run permutations
echo "Permuting cases" 
MakeBootstraps.py Data/SparseCC_input/Cases.tsv -p Data/SparseCC_input/Perm_cases/ -t 'permuted_#'
echo "Permuting controls"
MakeBootstraps.py Data/SparseCC_input/Controls.tsv -p Data/SparseCC_input/Perm_controls/ -t 'permuted_#'

#Get null distribution
echo "Generating NULL cases"
for i in `seq 0 99`; do SparCC.py Data/SparseCC_input/Perm_cases/permuted_$i -c Results/SparseCC/Cases/Perm/simulated_sparcc_$i.txt -i 20 -t 0.1  ; done
echo "Generating NULL controls"
for i in `seq 0 99`; do SparCC.py Data/SparseCC_input/Perm_controls/permuted_$i -c Results/SparseCC/Controls/Perm/simulated_sparcc_$i.txt -i 20 -t 0.1  ; done

#Compute Pvalues from Null
echo "Computing Ps in cases"
PseudoPvals.py Results/SparseCC/Cases/Corr_matrix_cases.txt 'Results/SparseCC/Cases/Perm/simulated_sparcc_#.txt' 100 -o Results/SparseCC/Cases/P_matrix.txt -t one_sided
echo "Computing Ps in controls"
PseudoPvals.py Results/SparseCC/Controls/Corr_matrix_controls.txt 'Results/SparseCC/Controls/Perm/simulated_sparcc_#.txt' 100 -o Results/SparseCC/Controls/P_matrix.txt -t two_sided

