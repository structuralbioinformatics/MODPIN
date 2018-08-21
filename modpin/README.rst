##############################################################################################
#                                                                                            #
#  MODPIN: SET OF PYTHON SCRIPTS TO MODEL AND ANALYZE PROTEIN-PROTEIN INTERACTIONS           #
#                                                                                            #
#  Authorship: LLUIS DOMINGUEZ, JAUME BONET, ORIOL FORNES, BALDO OLIVA                       #
#                                                                                            #
#  STRUCTURAL BIOINFORMATICS GROUP (GRIB-IMIM)                                               #
#  DEPARTMENT OF EXPERIMENTAL AND LIFE SCIENCES                                              #
#  UNIVERSITAT POMPEU FABRA                                                                  #
#                                                                                            #
#  BARCELONA, 1ST OCTOBER 2017                                                               #
#                                                                                            #
##############################################################################################
#
# BEFORE START YOU NEED TO INSTALL SOME PROGRAMS
#
# REQUIREMENTS
# These scripts work with Python 2.7.
# Use for example module load Python/2.7.11
# 
# Install programs 
#     blast-2.2.26  
#     cd-hit
#     clustalw-2.0.12
#     hbplus-3.2
#     reduce-3.23
#     zrank
#     modeller (version 9)
#     Rosetta
#
# We have installed some of them in modppi/src folder with other libraries such as archdb, SBI and functions
# 
# Set aliasses for these programs or load their modules, for example for cd-hit

alias cd-hit = /soft/EB_repo/bio/sequence/programs/goolf/1.7.20/CD-HIT/4.6.4/bin/cd-hit

#
# set-up the configuration file:
#    Addresses for all programs in "Paths"
#    Files with data sources (to be generated)
#    Parameters to run the programs
#    Parameters of the Cluster to run in parallel
#
# PATHS:
#
# Contains information of the paths of the programs required to run the scripts
#
# PARAMETERS:
#
# When generating the data we can change the parameters of the configuration file, for example:
# Use contact by minimum distance (PPI_threshold_type=min)  of 6A (PPI_distance_threshold=6) and store the files as pdb_Min6.ppi pdb_Min6.dat
# then recalculate by CB-CB distance (PPI_threshold_type=cb) of 12 A  (PPI_distance_threshold=12) and store them as pdb_CB12.ppi pdb_CB12.dat
# finally use a more stringent cut-off 4.5A and minimum distance and save it as pdb_Min4.ppi pdb_Min4.dat
# Select either pdb_Min4, pdb_Min6 or pdb_CB12 and use it as pdb.ppi
# Due to some errors in modelling, we suggest the use of pdb_Min4.ppi to generate the database, but then use CB-CB distance of 8.0 for modelling 
#
# To analyze mutants we define a shell (PPI_distance_threshold_shell)  for clustering mutants with the similar interface larger than 10A 
# (ie. from 12 to 20) and a minimum overlap of the same residue-residue interactions of 5% (overlap_interface=5). Also for the analyses a minimum
# number of poses with th same interface is required to cluster them in the same group of study (min_ratio_common_poses = 0.49)
# In order to evaluate the models we also decide wether to relax or not the model (Relax=yes) and if adding hydrogens (full, polar or none)
# Models will be done if they are acceptable according to the Rost curves strength parameter (n_parameter_twilight_zone = 0)
#
# FILES:
#
# Contains information on the files with sequences and template structures
# It will be gnerated during the installation of data
# We use the folder "data" under the installation of modppi, but this can be selected somewhere else
#
# 
# CLUSTER:
#
# Information of the cluster to run in parallel the modelling of a list of sequences
# cluster_name    Name of the cluster
# cluster_queue   Name of the queue (if any, otherwise use None)
# cluster_submit  Method of submission (qsub or sbatch)
# cluster_qstat   Method to list the number of runing queues
# max_jobs_in_queue 
# min_jobs_in_queue
# command_queue   Name of a file with information to submit in the queues (some nodes require to activate the programs necessary for the system)
#
#########################################################################################################################
#
# INSTALLATION OF DATA  (Files will be defined in the configuration file)
#
#
# STEP1:  DATA from PDB
# 
# Import PDB sources and create FastA file
#

python src/functions/PDB2SQL.py -d data/pdbsource -q data/pdbseq -s data/pdbsql -v

\ls `pwd`/data/pdbsource/*/pdb* > data/pdb.list
#new comand that should avoid the memory problem
for i in `pwd`/data/pdbsource/*/*.gz; do echo "$i"; done > output.txt


cat data/pdbseq/PDBseq.fa > data/pdb.fa

#
# Initialize and write the set of PPIs from PDB 
#


python scripts/modppi.py -init_pdb -v


#
# STEP2: DATA from 3DiD
#

python src/functions/3did.py -i data/3did_flat -ppi 3did.ppi -seq 3did.fasta -v --dummy=dummy 
\ls `pwd`/data/3did/*brk > data/3did.list

#
# STEP3: WHOLE  DATASET: Create the universe of sequences
#

cat data/pdbseq/PDBseq.fa >  data/nr.fa
cat data/3did.fasta >> data/nr.fa

#
# Avoid redundancies: Create a non-redundant set of sequences combining PDB and 3DiD 
#

cd-hit -i data/nr.fa -o data/nr90.fa -n 5 -c 0.9

#
# Create a non-redundant set of sequences joining PDB and 3DiD non-redundant
#

cd-hit -i data/pdbseq/PDBseq.fa -o data/pdb90.fa -n 5 -c 0.9
cd-hit -i data/3did.fasta -o data/3did90.fa  -n 5 -c 0.9
cat data/pdb90.fa > data/merged_nr90.fa
cat data/3did90.fa >> data/merged_nr90.fa
cat data/pdb90.fa.clstr > data/merged_nr90.fa.clstr
cat data/3did90.fa.clstr >> data/merged_nr90.fa.clstr

#
# STEP4: Initialize Dataset of sequences
#
#    Modify if necessary the config.ini file to select fasta_list_file  and nr90_list_file accordingly
#

python scripts/modppi.py -init_blast -3did -v

#
#
#  MODELLING TUTORIAL
#
#  EXAMPLE 1: Model building
#
#   Force the modeling even if the files already exist
#   Clean the dummy files
#   Add hydrogen atoms on the models
#   Skip checking the database of known interactions and their sequences
#   
#   NOTE: Sequence names should not use symbol "-" or ":" because both are used to model sub-regions or interactions
#         Mutant forms should be indicated with "_" and the mutation definition (i.e. L53R) at the end of the sequence name (see further in EXAMPLE 3)
#

python scripts/modppi.py -seq example/baxbid.fa -ppi example/baxbid.ppi -o example/BAX_BID -d ./dummy -v --hydrogens -skip -force -3did -n 5 --renumerate

#
#  EXAMPLE 2: Analyze the models generated in example 1
#

python scripts/analysis.py -o example/BAX_BID/ -l BAX_BID -zrank -v -d dummy -boxplot -ppi example/BAX_BID/interactions_done.list -seq example/baxbid.fa

#
#  EXAMPLE 3: RUN IN A CLUSTER SEVERAL SEQUENCES WIYTH MUTANT FORMS (METHOD "MODELIST")
#
#  For the method to work, it requieres that mutant forms are specified in the name of the FastA file
#  For example:
#  a) >NATIVE  and >NATIVE_G84D                                  show a wild type protein named NATIVE and its mutant form G84D
#  b) >sp|P55957_R84W|BID_HUMAN_R84W  and >sp|P55957|BID_HUMAN   show a wild type protein with accession P55957 and entry BID_HUMAN  and its mutant form R84W
#
#  Non-allowed forms are BID_R84W_HUMAN or sp|BID_HUMAN_R84W|P55957_R84W
#  A label can be added to classify your results or repeat a different run. However, the label should not contain symbols "_" , "-", ":" or white spaces, we recommend only alphanumeric symbols
#
#  Step 1: Create models of two large sets of mutations involving rewiring and unrewiring (hydrogens and relaxation won't be added in order to fasten up the obtantion of the models)
#

python scripts/modelist.py -i example/rewiring_mutants.dat  -l REWIRED -seq example/uniprot_Marc_Vidal.fasta -o example/rewired -n 3 -d dummy -3did -v --parallel -opt 
python scripts/modelist.py -i example/unrewiring_mutants.dat  -l UNREWIRED -seq example/uniprot_Marc_Vidal.fasta -o example/unrewired -n 3 -d dummy -3did -v --parallel -opt 

#
# Step 2: Then continue the analysis (add hydrogens and optimize structure)
#         The file with hydrogens will be the same as the PDB input with extension ".h" instead of ".pdb"
#

python scripts/modelist.py -i example/rewiring_mutants.dat  -l REWIRED -seq example/uniprot_Marc_Vidal.fasta -o example/rewired -n 3 -d dummy -3did -v --parallel -opt -a --continue --hydrogens --renumerate 
python scripts/modelist.py -i example/unrewiring_mutants.dat  -l UNREWIRED -seq example/uniprot_Marc_Vidal.fasta -o example/unrewired -n 3 -d dummy -3did -v --parallel -opt -a --continue  --hydrogens --renumerate

# Step 3: To further analyze and compare the whole set of mutations we run the script "select_cluster". 
#         This produces a selection of clustered interfaces for each complex with mutant forms compared with the wild type 
#         according to two statistics (Mann-Whitney and Kolmnorov-Smirnov). 
#         Plots are collected in the output folders  example/rewired and example/unrewired
#         Then we can run Exponential_Averaging_FEP.py that uses the selected folder (i.e. example/rewired and example/unrewired)
#         to calculate the differences of Free Energy between wild-type and mutant forms using Zwanzig formulae
#         The partition finction Z is calculated with the score of the last variable (ie. ddG_mean) and the variable to compare is the
#         first selected (i.e. dGx_mean, from dG_cross of Rosetta)
#

python src/functions/select_cluster.py example/rewired
python src/functions/select_cluster.py example/unrewired
python src/functions/Exponential_Averaging_FEP.py  example/rewired  dGx_mean ddG_mean
python src/functions/Exponential_Averaging_FEP.py  example/unrewired  dGx_mean ddG_mean



# EXAMPLE 4: Extra utilities

# If you forgot to renumearte a file you can fix it when the model is done. 
# The program requires the original sequences in FastA format and an input file to indicate the sequence name of each chain in the complex
# File example/renumerate_P02792_H133P_P02794.dat contains the sequences of each chain P02792_H133P (FRIL_HUMAN_H133P in chain A) and P02794 (FRIH_HUMAN in chain B)
# and the PDB file to rnumber is test4renumber.pdb

python scripts/Renumerate.py -f example/sequences_Marc_Vidal.fasta -p example/test4renumber.pdb -d dummy -l example/renumerate_P02792_H133P_P02794.dat -o example/renumbered.pdb -v -hydro

# EXAMPLE 5: Analyse change of energy in mutations of SKEMPI

# First we need to parse the data of SKEMPI, this can be done with the parser parser_skempi2ModPIN.py
# -Download data from SKEMPI: table skempi_v2.csv and PDB files in folder PDBs
# -Run the parser for a maximum number of mutations in both chains (M=2)
#

python parser_skempi2ModPIN.py -i skempi_v2.csv -p ./PDBs -o modppi_skempi_short -v -m 2

# This generates the files 

 modppi_skempi_short.ddG  
 modppi_skempi_short.fa   
 modppi_skempi_short.ppi

# Then we can run modelist with modppi_skempi_short.ppi and modppi_skempi_short.fa

python scripts/modelist.py -i example/modppi_skempi_short.ppi  -l SKEMPI -seq example/modppi_skempi_short.fa -o example/SKEMPI -n 3 -d dummy -3did -v --parallel -opt  --renumerate 
python scripts/modelist.py -i example/modppi_skempi_short.ppi  -l SKEMPI -seq example/modppi_skempi_short.fa -o example/SKEMPI -n 3 -d dummy -3did -v --parallel -opt -a --continue --hydrogens --renumerate 
 
# For analysing the energies we can run the bash script:


set home="<MODPIN HOME>"

foreach type ( "example/SKEMPI/")

set models=`ls ${type}/models`
echo $models

foreach score ( "ddG_all" "ddG_mean" "dGx_all" "dGx_mean" "zrank" "dSASA_all" "dSASAp_all" "dSASAh_all" "dSASA_mean" "dSASAp_mean" "dSASAh_mean")

python  src/functions/select_cluster.py ${type} ${score}
\rm  $home/${type}/${score}.rank

foreach ene ( "ddG_all" "ddG_mean" "dGx_all" "dGx_mean" "zrank")
python  src/functions/Exponential_Averaging_FEP.py ${type}  ${score} ${ene}
end

foreach i ( $models )

 echo "cat  ${type}/models/${i}/*cluster_*${score}*out|fgrep ".h"  >>  $home/${type}/${score}.rank"
 cat  ${type}/models/${i}/*cluster_*${score}*out|fgrep ".h"  >>  $home/${type}/${score}.rank
 echo "cat  ${type}/models/${i}/*cluster_*${score}*out|fgrep ".h" |sort -g -k 2  > $home/${type}/${score}_${i}.rank"
 cat  ${type}/models/${i}/*cluster_*${score}*out|fgrep ".h" |sort -g -k 2  > $home/${type}/${score}_${i}.rank

end
end
end

# Finally, we can compare the predicted energies with the real ones, for example:

#Using Zwanzig approach on FEP:

python utils/Analyse_ddG_FEP.py -a example/modppi_skempi_short.ddG -b example/SKEMPI/FEP_zrank_Z_with_zrank.out -lp 1.0 -o compare_zrank_ddGskempi -v 

#Using the comparison of energy distributions from clusters of conformations with similar interface
# by selecting the specific conditions of the cluster, such as the percentage of mutations in the interface (PMI)
# the p-value of significance when comparing the distributions of two forms and the ranking of the cluster, whhich depends
# on the number of models or conformations in the cluster

python utils/Analyse_ddG_cluster.py -g example/SKEMPI/modppi_skempi_short.ddG -l "i1p0005r10" -d example/SKEMPI -e zrank -o example/SKEMPI/Selected -i 1.0 -p 0.005 -r 10 -lp 1.0  -v -k 0.01

#But we can also use an automatic search of these conditions with the script:
# where the conditions shown in the input are just the thresholds
# i => the minimum PMI to start, up to 100% in steps of 1%
# p => the minimum p-value to start, up to 1
# r => the maximum rank, starting from 1

python utils/Analyse_ddG_cluster_best.py -g example/SKEMPI/modppi_skempi_short.ddG -l "optim_rank" -d example/SKEMPI -e zrank -o example/SKEMPI/Selected -i 1.0 -p 0.005 -r 10 -lp 1.0  -v -k 1.0


