import sys
import argparse
import os
import re
import ConfigParser
import shutil
import subprocess
import time
import random
from collections import Counter

# Add '.' to sys.path
src_path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(src_path)

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, 'config.ini')
config.read(config_file)

# Add SBI library path to sys.path
src_path = os.path.join(config.get('Paths','modppi_path'),config.get('Paths','sbi_library_path'))
sys.path.append(src_path)

# Add functions path to sys.path
src_path = os.path.join(config.get('Paths','modppi_path'),config.get('Paths','functions_path'))
sys.path.append(src_path)

# Import functions
from SeqIO import *
from functions import *
from SBI.external.blast import blast_parser
from SBI.structure.contacts import Complex
from SBI.sequence import Sequence
from SBI.structure import PDB
from FilterHomologPairs import filter_homologs
from Bio import SeqIO
from Bio import ExPASy
from Bio import AlignIO
from Bio.Align.Applications import *
import BioLib


def main():
    
    #initialize
    options    = parse_user_arguments()
    verbose    = options.show
    zrank      = options.zrank
    foldx      = options.foldx
    ssp        = options.ssp
    rosetta    = options.rosetta
    input_list = os.path.join(options.ppi_list)
    output_dir = options.outdir
    hbplus     = options.hbplus
    label      = options.label
    boxplot    = options.boxplot

    #Check that at least one method to calculate the energy is set
    if not zrank and not foldx and not ssp and not rosetta:
       ssp=True

    if options.ppi_list is None:
      sys.stderr.write("ERROR: Missing argument input list '-i'\n")
      exit(0)
    if not fileExist(input_list):
      sys.stderr.write("ERROR: Missing file %s\n"%input_list)
      exit(0)
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    dummy_dir=os.path.join(options.dummy_dir,"analysis")
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    #Write  distributions for boxplot
    if options.outdir is None:
        output_dir=os.path.abspath("analysis")
    else:
        output_dir=options.outdir
    if not os.path.exists(output_dir):os.makedirs(output_dir)
         

    #Create a parsed list of models (modify their names in case hydrogens need to be added)
    # information from parser_list_of_models is a dictionary:
    # key= 2-tuple (pair1, pair2) with names of the interactors
    # value= list [] with 2-tuples of (name(number),name(file_list)) of the cluster name(pose) and the name of the file with the list of PDB files
    list_of_models=parser_list_of_models(input_list, options)

    #Check the sequences of the pairs and the differences of the interface for each cluster-pose
    # information from compare_interface_sequences is a dictionary:
    # key= (cluster-pose name,(pair1,pair2)) 2-tuple (pair1,pair2) are names of interactors
    # value=  2-tuples of:
    #    2-tuple (ratio1,ratio2) with ratios of structures where a difference in sequence with respect to its wild-type is involved in the interface
    #
    pose_interfaces=compare_interface_sequences(list_of_models,dummy_dir)
    

    binding={}

    #Run zrank on folders if active option zrank
    if zrank:
     if verbose: sys.stdout.write("\t-- Calculate ZRANK...\n")
     try:
       binding_by_zrank= calculate_zrank(list_of_models,dummy_dir,label,verbose)
     except Exception as e:
       sys.stderr.write("FAIL %s\n"%(e))
       exit(0)
     for key,value in binding_by_zrank.iteritems():
         binding.setdefault(key,value)

    #Run foldx on folders if active option foldx
    if foldx:
     if verbose: sys.stdout.write("\t-- Calculate FOLDX...\n")
     try:
       binding_by_foldx= calculate_foldx(list_of_models,dummy_dir,label,verbose)
     except Exception as e:
       sys.stderr.write("FAIL %s\n"%(e))
       exit(0)
     for key,value in binding_by_foldx.iteritems():
         binding.setdefault(key,value)

    #Run SSP on folders if active option ssp
    if ssp:
     if verbose: sys.stdout.write("\t-- Calculate by Split Potentials...\n")
     try:
       binding_by_ssp= calculate_split_potentials(list_of_models,dummy_dir,label,verbose)
     except Exception as e:
       sys.stderr.write("FAIL %s\n"%(e))
       exit(0)
     for key,value in binding_by_ssp.iteritems():
         binding.setdefault(key,value)

    #Run Rosetta on folders (by default is always done)
    if rosetta:
       if verbose: sys.stdout.write("\t-- Calculate ddG with Rosetta...\n")
       try:
         binding_by_rosetta=calculate_rosetta(list_of_models,dummy_dir,label,verbose)
       except Exception as e:
         sys.stderr.write("FAIL %s\n"%(e))
         exit(0)
       for key,value in binding_by_rosetta.iteritems():
           binding.setdefault(key,value)
   
    #Prepare Plots
    energy_type=set()
    pair_list=set()
    pose_list=set()
    plot_list=set()
    for ((a,b),pose,key),outfile  in binding.iteritems():
        energy_type.add(key)
        pair_list.add((a,b))
        pose_list.add(pose)
    for key in energy_type:
        output=os.path.join(output_dir,label+"."+key+"_distribution.dat") 
        name_file=os.path.join(output_dir,label+"."+key+"_distribution")
        fo=open( output,"w")
        for pose in pose_list:
         if pose == 0:
           output_pose=os.path.join(output_dir,label+"."+key+".all_distribution.dat") 
           name_file_pose=os.path.join(output_dir,label+"."+key+".all_distribution")
         else:
           output_pose=os.path.join(output_dir,label+"."+key+".cluster_"+str(pose)+"_distribution.dat") 
           name_file_pose=os.path.join(output_dir,label+"."+key+".cluster_"+str(pose)+"_distribution")
         fp=open( output_pose,"w")
         for a,b in pair_list:
           if binding.has_key(((a,b),pose,key)):
            outfile= binding[((a,b),pose,key)]
            if pose_interfaces.has_key((pose,(a,b))):
             fp.write("%s %s\n"%(a+"-"+b+"#"+str(pose)+interface_affected(pose_interfaces[pose,(a,b)]),outfile))
         fp.close()
         plot_list.add((key,name_file_pose,1,len(pair_list))) 
        ncmax=0
        for a,b in pair_list:
          nc=0
          for pose in pose_list:
            if pose==0:continue
            if binding.has_key(((a,b),pose,key)):
             nc=nc+1
             if nc>ncmax:ncmax=nc
             outfile= binding[((a,b),pose,key)]
             if pose_interfaces.has_key((pose,(a,b))):
              fo.write("%s %s\n"%(a+"-"+b+"#"+str(pose)+interface_affected(pose_interfaces[pose,(a,b)]),outfile))
        fo.close()
        plot_list.add((key,name_file,ncmax,ncmax*len(pair_list)))

    #Automatic BOXPLOTS
    if boxplot:
       boxplot_path = os.path.join(config.get('Paths','modppi_path'),config.get('Paths','functions_path'))
       boxplot_exe  = os.path.join(boxplot_path,"boxplot.py")
       python_path  = config.get('Paths', 'python_path')
       python_exe   = os.path.join(python_path,"python")
       for key,name,nc,nh in plot_list:
         input_list  = name + ".dat"
         boxplot_name= name
         output_name = name + "_statistic.out"
         os.system('%s %s -l %s -o %s -nc %d -nh %d -t "Distribution of %s" >& %s'%(python_exe,boxplot_exe,input_list,boxplot_name,nc,nh,key,output_name))
         if verbose: sys.stdout.write('\t-- boxplot SCORE %s INPUT %s OUTPUT %s STATISTICS %s\n'%(key,input_list,boxplot_name,output_name))

    #Remove dummy folders
    if not verbose: shutil.rmtree(dummy_dir)


def calculate_zrank(list_of_models,dummy_dir,label,verbose):

  #Initialize
  zrank_exe  = config.get('Paths','zrank_path')
  binding={}

  for pair,folders in list_of_models.iteritems():
    for pose,folder in folders:
       #if verbose: sys.stdout.write("\t\t-- Use list %s for %s...\n"%(folder.split("/")[-1],pair))
       if verbose: sys.stdout.write("\t\t-- Use list %s for %s...\n"%(folder,pair))
       output_zrank=folder+".zr.out"
       output_zrank_label=folder+".zrank."+label+".out"
       output_zrank_failed=folder+".zrank."+label+".fail"
       if not fileExist(output_zrank_label):
        try:
         #process=subprocess.Popen([zrank_exe,folder], stderr = subprocess.STDOUT)
         os.system("%s %s"%(zrank_exe,folder))
        except Exception as e:
         sys.stderr.write("ERROR: %s\n"%e)
         dummy_file=os.path.join(dummy_dir,folder.split("/")[-1])
         dummy_list=open(os.path.join(dummy_dir,folder.split("/")[-1]),"w")
         if verbose: sys.stdout.write("\t\t\t-- Create new list for %s (%s)...\n"%(str(pair),dummy_file))
         path="/".join(folder.split("/")[:-1])
         for models_ready in os.listdir(path):
           if models_ready.strip().endswith(".pdb"):
              dummy_list.write("%s\n"%os.path.join(path,models_ready.strip()))
         dummy_list.close()
         if fileExist(output_zrank): shutil.move(output_zrank,output_zrank_failed)
         shutil.move(dummy_file,folder)
        if not fileExist(output_zrank):
         try:
          #process=subprocess.Popen([zrank_exe,folder], stderr = subprocess.STDOUT)
          #process.wait()
          os.system("%s %s"%(zrank_exe,folder))
         except Exception as e:
          sys.stderr.write("ERROR: %s\n"%e)
        if fileExist(output_zrank):
         check_file=open(folder,"r")
         check_zrank=open(output_zrank,"r")
         check_file_set=set()
         check_zrank_set=set()
         for models_ready in check_file:
             check_file_set.add(models_ready.strip())
         for models_ready in check_zrank:
             check_zrank_set.add(models_ready.strip().split()[0])
         check_file.close()
         check_zrank.close()
         if check_file_set!=check_zrank_set:
            sys.stderr.write("ERROR: ZRANK failed for some models\n")
            if verbose: sys.stdout.write("\t\t\t-- Create new list for %s...\n"%(str(pair)))
            dummy_file=os.path.join(dummy_dir,folder.split("/")[-1])
            dummy_list=open(os.path.join(dummy_dir,folder.split("/")[-1]),"w")
            path="/".join(folder.split("/")[:-1])
            for models_ready in os.listdir(path):
              if models_ready.strip().endswith(".pdb"):
                 dummy_list.write("%s\n"%(os.path.join(path,models_ready.strip())))
            dummy_list.close()
            shutil.move(output_zrank,output_zrank_failed)
            shutil.move(dummy_file,folder)
        if not fileExist(output_zrank):
         try:
          #process=subprocess.Popen([zrank_exe,folder], stderr = subprocess.STDOUT)
          #process.wait()
          os.system("%s %s"%(zrank_exe,folder))
         except Exception as e:
          sys.stderr.write("ERROR: %s\n"%e)
        if not fileExist(output_zrank):
         sys.stderr.write("ERROR: %s doesn't exist\n"%output_zrank)
         sys.stderr.write("SKIP:  %s\n"%output_zrank_label)
        else:
         shutil.move(output_zrank,output_zrank_label)
       #When finally there is a result for zrank added
       binding.setdefault((pair,pose,"zrank"),output_zrank_label)

  return binding


def calculate_foldx(list_of_models,dummy_dir,label,verbose):

  #Initialize
  foldx_path  = config.get('Paths','foldx_path')
  foldx_exe   = os.path.join(foldx_path,"foldx")
  rotabase    = os.path.join(foldx_path,"rotabase.txt")
  binding={}

  for pair,folders in list_of_models.iteritems():
    for pose,folder in folders:
       #if verbose: sys.stdout.write("\t\t-- Use list %s for %s...\n"%(folder.split("/")[-1],pair))
       if verbose: sys.stdout.write("\t\t-- Use list %s for %s...\n"%(folder,pair))
       output_foldx=folder+".foldx."+label+".out"
       if verbose: sys.stdout.write("\t\t-- Output %s...\n"%(output_foldx))
       if not fileExist(output_foldx):
         output_dummy=os.path.basename(folder)+".foldx."+label+".out"
         foldx_list=[]
         fo=open(folder,"r")
         if verbose: sys.stdout.write("\t\t-- Open %s ...\n"%(folder))
         for line in fo:
           #print line.strip()
           p=line.strip()
           if verbose: sys.stdout.write("\t\t\t-- Use %s ...\n"%os.path.basename(p))
           try:
             dummy_file=os.path.join(dummy_dir,os.path.basename(p))
             #print "COPY",dummy_file
             shutil.copyfile(p,dummy_file)
             #print "ADD",p
             foldx_list.append(p)
           except:
             sys.stdout.write("\t\t\t\t-- Skip %s copy cpuld not be added\n"%p)
             continue 
         if verbose: sys.stdout.write("\t\t-- Close %s ...\n"%(folder))
         fo.close()
         if verbose:
            print "\t ===  DATA  ==="
            print foldx_list
            print "\t =============="
         cwd = os.getcwd()
         os.chdir(dummy_dir)
         energy={}
         for pdb_check in foldx_list:
           #print "CHECK",pdb_check
           pdb_name=os.path.basename(pdb_check)
           pdb_root="".join(pdb_name.split(".")[0:-1])
           if not fileExist(pdb_name):
              if verbose: sys.stdout.write("\t\t\t-- File no exist %s ...\n"%pdb_name)
              #raise ValueError("Missing file %s"%(pdb_name))
              continue
           # execute foldx optimization
           pdb_repair=pdb_root+"_Repair.pdb"
           pdb_repair_home = os.path.join(os.path.dirname(pdb_check),pdb_repair)
           if verbose: sys.stdout.write("\t\t-- Optimize file %s ...\n"%pdb_name)
           if fileExist(pdb_repair_home):
              shutil.copyfile(pdb_repair_home,pdb_repair)
           else:
              if verbose: sys.stdout.write("%s  --rotabaseLocation %s --repair_Interface ONLY  --command RepairPDB --pdb=%s >& %s"%(foldx_exe,rotabase,pdb_name,pdb_root+"_Repair.log\n"))
              try:
                os.system("%s  --rotabaseLocation %s --repair_Interface ONLY  --command RepairPDB --pdb=%s >& %s"%(foldx_exe,rotabase,pdb_name,pdb_root+"_Repair.log"))
                shutil.copyfile(pdb_repair,pdb_repair_home)
              except:
                print "FAILED FOLDX REPAIR %s"%pdb_repair
           # execute foldx AnalyseComplex
           if not fileExist(pdb_repair):
              if verbose: sys.stdout.write("\t\t\t-- File optimized has failed  %s ...\n"%pdb_repair)
              #raise ValueError("Missing file %s"%(pdb_repair))
           if verbose: sys.stdout.write("%s  --rotabaseLocation %s --command AnalyseComplex --analyseComplexChains=A,B --pdb=%s >& %s"%(foldx_exe,rotabase,pdb_repair,pdb_root+"_AC.log\n"))
           try:
             os.system("%s  --rotabaseLocation %s --command AnalyseComplex --analyseComplexChains=A,B --pdb=%s >& %s"%(foldx_exe,rotabase,pdb_repair,pdb_root+"_AC.log"))
           except:
             print "FAILED FOLDX %s"%pdb_root
           #parse foldx output
           pdb_fxout="Interaction_"+pdb_root+"_Repair_AC.fxout"
           if verbose: sys.stdout.write("\t\t\t-- Open %s ...\n"%(pdb_fxout))
           if fileExist(pdb_fxout):
              fo=open(pdb_fxout,"r")
              for line in fo:
                  word=line.strip().split()
                  if len(word)>0:
                   if word[0].lstrip("./") == pdb_repair:
                     energy.setdefault(pdb_check,float(word[5]))
              fo.close()
           else:
              #raise ValueError("Missing file %s"%(pdb_fxout))
              sys.stdout.write("\t\t\t-- Missing file %s"%(pdb_fxout))
              continue
         if verbose: sys.stdout.write("\t\t\t-- Use FoldX output %s ...\n"%(output_dummy))
         fo=open(output_dummy,"w")
         for pdb,ene in energy.iteritems():
           fo.write("%s\t%s\n"%(pdb,ene))
         fo.close()
         shutil.move(output_dummy,output_foldx)
         os.chdir(cwd)
       else:
         if verbose: sys.stdout.write("\t\t\t-- Use FoldX output %s ...\n"%(os.path.basename(output_foldx)))
       #Add result of foldx
       binding.setdefault((pair,pose,"foldx"),output_foldx)

  return binding

class SPPPI(object):
    '''
    Analyze a protein folds using the split potentials and prints an XML file
    '''
    def __init__(self, receptor, ligand, pot_type):
        '''
        Contructor
        '''
        self.receptor    = receptor
        self.ligand      = ligand
        self.pot_type    = pot_type

 
    def get_energies(self):
        '''
        Compute the split potentials for a ppi
        '''
        receptor=self.receptor
        ligand  =self.ligand

        receptor.set_dssp()
        receptor.normalize_residues()
        ligand.set_dssp()
        ligand.normalize_residues()

        ppi = BioLib.Interaction(receptor, ligand)

        # Compute split potentials
        if self.pot_type == 'CB':
            cutoff = 12
        else: 
            cutoff = 5

        split_potentials = BioLib.SplitPotentialsPPI(c_type=self.pot_type, cutoff=cutoff) # definition of SplitPotentialsPPI class (BioLib.Docking)
        global_energies = split_potentials.calculate_global_energies(ppi, Zscores=True) # calculation of global energies using the method of the SplitPotentialsPPI class

        return global_energies



def calculate_split_potentials(list_of_models,dummy_dir,label,verbose):

  #Initialize
  ssp_list = config.get('Parameters','ssp_score')
  ssp_type = config.get('Parameters','ssp_type')

  ssp_score_list=ssp_list.split(";")

  binding={}

  for pair,folders in list_of_models.iteritems():
   for pose,folder in folders:
    for ssp in ssp_score_list:
       ssp_score=ssp.replace(" ","")
       #if verbose: sys.stdout.write("\t\t-- Use list %s for %s...\n"%(folder.split("/")[-1],pair))
       if verbose: sys.stdout.write("\t\t-- Use list %s for %s...\n"%(folder,pair))
       output_ssp=folder+"."+ssp_score+"."+label+".out"
       if not fileExist(output_ssp):
         output_dummy=os.path.basename(folder)+"."+ssp_score+"."+label+".out"
         ssp_list=[]
         fo=open(folder,"r")
         for line in fo:
           p=line.strip()
           dummy_file=os.path.join(dummy_dir,os.path.basename(p))
           shutil.copyfile(p,dummy_file)
           ssp_list.append(p)
         fo.close()
         cwd = os.getcwd()
         os.chdir(dummy_dir)
         energy={}
         for pdb_check in ssp_list:
           pdb_name=os.path.basename(pdb_check)
           pdb_root="".join(pdb_name.split(".")[0:-1])
           if not fileExist(pdb_name):
              raise ValueError("Missing file %s"%(pdb_name))
           # execute ssp 
           pdb_ppi=BioLib.PDB.read_pdb(pdb_name)
           receptor     = pdb_ppi.pop()
           ligand       = pdb_ppi.pop()
           ssp_ppi      = SPPPI(receptor,ligand,ssp_type)
           ssp_energies = ssp_ppi.get_energies()
           zscore=0
           if ssp_score.lower().startswith("z"):   zscore=1
           if ssp_score.lower().endswith("local"): ssp_ene="D-LOCAL"
           if ssp_score.lower().endswith("pair"):  ssp_ene="D-PAIR"
           if ssp_score.lower().endswith("comb"):  ssp_ene="D-COMB"
           if ssp_score.lower().endswith("s3dc"):  ssp_ene="D-S3DC"
           energy.setdefault(pdb_check,float(ssp_energies[zscore][ssp_ene]))
         if verbose: sys.stdout.write("\t\t\t-- Use Split Potentials output %s ...\n"%(output_dummy))
         fo=open(output_dummy,"w")
         for pdb,ene in energy.iteritems():
           fo.write("%s\t%s\n"%(pdb,ene))
         fo.close()
         shutil.move(output_dummy,output_ssp)
         os.chdir(cwd)
       else:
         if verbose: sys.stdout.write("\t\t\t-- Use Split Potentials output %s ...\n"%(os.path.basename(output_ssp)))
       #Add result of ssp
       binding.setdefault((pair,pose,ssp_score),output_ssp)

  return binding

def calculate_rosetta(list_of_models,dummy_dir,label,verbose):

    binding={}

    if verbose: sys.stdout.write("\t\t-- Cleaning output files...\n")
    for pair,folders in list_of_models.iteritems():
     for pose,folder in folders:
      path_folder="/".join(folder.split("/")[:-1])
      clean_file =folder.split("/")[-1]
      if verbose: sys.stdout.write("\t\t\t Clean directory %s \n"%path_folder)
      for files in os.listdir(path_folder):
        m=re.search("^%s"%clean_file,files)
        if m:
          line=files.lstrip(folder).rstrip().split(".")
          if line[-1]=="out":
            if len(line[0].split("_")) > 1:
             check=line[0].split("_")[1]
             if check == "mean" or check == "all":
                if verbose: sys.stdout.write("\t\t\t\t Remove file %s \n"%files)
                os.remove(os.path.join(path_folder,files))
    for pair,folders in list_of_models.iteritems():
     for pose,folder in folders:
      if verbose: sys.stdout.write("\t\t-- Use list %s for %s...\n"%(folder.split("/")[-1],pair))
      fd=open(folder,"r")
      for line in fd:
        pdb_file=line.strip().split()[0]
        energy_file=pdb_file+"."+"interface_analyzer"+"."+label+".out"
        if verbose: sys.stdout.write("\t\t\t-- Interface Analyzer output: %s\n"%energy_file.split("/")[-1])
        try:
          binding_data=binding_energy(pdb_file,energy_file,dummy_dir)
          energy_type=set()
          for key,score in binding_data.iteritems():
            feature,chain_A,chain_B=key
            feature_all=feature+"_all"
            out_feature_all=folder+"."+feature+"_all."+label+".out"
            fa=open(out_feature_all,"a")
            fa.write("%s\t%f\n"%(pdb_file+"_"+chain_A+chain_B,score))
            fa.close()
            energy_type.add(feature)
          for feature in energy_type:
            feature_mean=feature+"_mean"
            out_feature_mean=folder+"."+feature+"_mean."+label+".out"
            fm=open(out_feature_mean,"a")
            dat=[score for (e,a,b),score in binding_data.iteritems() if feature == e]
            if len(dat)>0:
              mean=sum(dat)/len(dat)
              fm.write("%s\t%f\n"%(pdb_file,mean))
            fm.close()
        except ValueError as e:
          sys.stderr.write("ERROR %s\n"%e)
      fd.close()
      for feature in energy_type:
          binding.setdefault((pair,pose,feature+"_all"),folder+"."+feature+"_all."+label+".out")
          binding.setdefault((pair,pose,feature+"_mean"),folder+"."+feature+"_mean."+label+".out")

    return binding

              
def compare_templates(sa,sb):

   new_sa=set()
   new_sb=set()
   for p in sa:
     a,b = p.rstrip(".pdb").split("::")
     template_a="_".join(a.split("_")[:-1])
     template_b="_".join(b.split("_")[:-1])
     new_sa.add((template_a,template_b))
   for p in sb:
     a,b = p.rstrip(".pdb").split("::")
     template_a="_".join(a.split("_")[:-1])
     template_b="_".join(b.split("_")[:-1])
     new_sb.add((template_a,template_b))


   return (new_sa,new_sb)
          
def binding_energy(pdb_file,energy_file,dummy_dir="/tmp"):
   ''' 
    Binding_energy function
    uses a pdb_file name as input and the name energy_file to write the output of energies
    The function returns a dictionary:
     key= 3-tuple (a,b,c) of
          a=the type of energy score from Rosetta: ddG, dGx, dSASA...
          b=chain name of one partner,A
          c=chain name of the interactor,B
     value=score
   '''
   #Initialize
   rosetta=config.get('Paths','rosetta_path')
   interface_analyzer=config.get('Paths','interface_analyzer')
   database=os.path.join(rosetta,'database')
   binding={}
   if fileExist(pdb_file):
     pdb=PDB(pdb_file) 
   else:
     raise ValueError("Missing file %s"%(pdb_file))
   chains=[x for x in pdb.chain_identifiers]
   for x in chains:
    for y in chains:
     if x==y:continue
     output=energy_file+"_"+x+"_"+y
     commands=[]
     commands.append("%s"%interface_analyzer)
     commands.append("-pack_input")
     commands.append("-pack_separated")
     commands.append("-ignore_zero_occupancy=False")
     commands.append("-in:ignore_unrecognized_res")
     commands.append("-add_regular_scores_to_scorefile")
     commands.append("-run:constant_seed")
     commands.append("-nodelay")
     commands.append("-database %s"%(database))
     commands.append("-s %s"%(pdb_file))
     commands.append("-out:file:score_only %s"%(output))
     commands.append("-fixedchains %s %s"%(x,y))
     if not fileExist(output):
        #print " ".join(commands)
        os.system("%s >& %s"%(" ".join(commands),os.path.join(dummy_dir,"interface_analyzer_"+output.split("/")[-1]+".log")))
        #process=subprocess.Popen(commands, stderr = subprocess.STDOUT, stdout = subprocess.STDOUT)
        #process.wait()
     try:
      f=open(output,"r")
      for line in f:
       word=line.strip().split()
       if word[0]=="SCORE:" and word[1]!="total_score":
        dGx   =float(word[3])
        ddG   =float(word[5])
        dSASAh=float(word[7])
        dSASA =float(word[8])
        dSASAp=float(word[9])
      f.close()
      binding.setdefault(("dGx",x,y),dGx)
      binding.setdefault(("ddG",x,y),ddG)
      binding.setdefault(("dSASAh",x,y),dSASAh)
      binding.setdefault(("dSASA",x,y),dSASA)
      binding.setdefault(("dSASAp",x,y),dSASAp)
     except:
      sys.stderr.write("ERROR: Rosetta FAIL to obtain %s \n"%(output))

   return binding


def parser_list_of_models(input_list,options):

    renumerate=options.renumerate
    fasta_file=options.fasta_file
    dummy_dir=options.dummy_dir
    verbose  =options.show

    if not fileExist(input_list):
      sys.stderr.write("Missing input file %s \n"%input_list)
      exit(0)

    if not fileExist(fasta_file):
      sys.stderr.write("Missing FastA file  %s \n"%fasta_file)
      exit(0)

    #Dictionary with results
    new_models={}
    group_models={}
    #Select the active models
    if verbose: sys.stdout.write('\t\t-- Open list of models %s ...\n'%input_list)
    fd=open(input_list,"r")
    models={}
    for line in fd:
      data=line.strip().split('\t')
      if data[2]=="DONE" or data[2]=="FIXED":
       pair=(data[0],data[1])
       if models.has_key(pair):
          models.pop(pair)
          models.setdefault(pair,data[3])
       else:
          models.setdefault(pair,data[3])
          if verbose: sys.stdout.write('\t\t\t-- Add models for %s from %s ...\n'%(pair,data[3]))
    fd.close()

    #If we need to renumerate PDBs we need to read the FastA file
    proteins= {}
    if verbose: sys.stdout.write('\t\t-- Reading protein sequences to renumerate files...\n')
    for protein in FASTA_iterator(fasta_file):
            name = protein.get_identifier()
            if len(name.split("|"))>1:
               name1 = name.split("|")[1] 
               name2 = name.split("|")[2] 
            if len(name.split("|")) == 2:
               name1 = name.split("|")[1]
               name2 = name.split("|")[1]
            if len(name.split("|"))<=1:
               name1 = name
               name2 = name
            seq  = protein.get_sequence()
            if len(seq) > 0:
                try:
                    proteins.setdefault(name1, ProteinSequence(name1, seq))
                    proteins.setdefault(name2, ProteinSequence(name2, seq))
                except IncorrectSequenceLetter as e:
                    sys.stderr.write('WARNING: %s\n' %e)
                    sys.stdout.write('\t\t-- Skip input sequence: %s Sequence: %s\n' %(name, seq))


    #If we need to add hydrogens / renumerate PDBs
    if options.hbplus:
      for pair,folder in models.iteritems():
        fo=open(folder+".h","w")
        fd=open(folder,"r")
        for line in fd:
         path_array= line.strip().split("/")
         path_file = "/".join([str(x) for x in path_array[:-1]])
         input_file=path_array[-1]
         if not fileExist(os.path.join(path_file,input_file)): continue
         if input_file.split(".")[-1]=="pdb":
          output_file=path_file+"/"+".".join([str(x) for x in input_file.split(".")[:-1]])+".h"
         else:
          output_file=path_file+"/"+input_file+".h"
         if not fileExist(output_file):
          add_hydrogens(config,path_file,input_file,output_file,dummy_dir)
         fo.write("%s\n"%output_file)
         if renumerate:
          if verbose: sys.stdout.write("\t\t\t-- Renumerate residues as original sequence %s\n"%input_file)
          sequences_complex = {}
          if proteins.has_key(pair[0]) and proteins.has_key(pair[1]):
           sequences_complex.setdefault("A",proteins.get(pair[0]))
           sequences_complex.setdefault("B",proteins.get(pair[1]))
           try:
                pdb_renumber=PDB()
                pdb_renumber=renumber_pdb(config,path_file,output_file,sequences_complex,dummy_dir)
                if len(pdb_renumber.chain_identifiers)>1:
                   pdb_renumber.write(output_file,force=True)
                else:
                   if verbose: sys.stdout.write("\t\t\t-- Failed to renumerate %s\n"%output_file)
           except Exception as e:
                sys.stderr.write("WARNING %s\n"%e)
          else:
           sys.stdout.write('\t\t\t-- Skip renumbering of %s (missing FastA sequences %s %s)\n'%(output_file,pair[0],pair[1]))
        fd.close()
        fo.close()
        new_models.setdefault(pair,folder+".h")
        pose_list= os.path.join(path_file,options.label+"_"+pair[0]+"::"+pair[1]+"_cluster_"+str(0)+".list")
        shutil.copy(folder,pose_list)
        group_models.setdefault(pair,[]).append((0,pose_list))
    else:
      new_models=models
      for pair,folder in models.iteritems():
        path_array= folder.strip().split("/")
        path_file = "/".join([str(x) for x in path_array[:-1]])
        pose_list= os.path.join(path_file,options.label+"_"+pair[0]+"::"+pair[1]+"_cluster_"+str(0)+".list")
        shutil.copy(folder,pose_list)
        group_models.setdefault(pair,[]).append((0,pose_list))
        fd=open(folder,"r")
        for line in fd:
          path_array= line.strip().split("/")
          path_file = "/".join([str(x) for x in path_array[:-1]])
          input_file=path_array[-1]
          output_file=path_file+"/"+input_file
          if not fileExist(os.path.join(path_file,input_file)): continue
          if renumerate:
            if verbose: sys.stdout.write("\t\t-- Renumerate residues as original sequences (%s %s) in %s\n"%(pair[0],pair[1],input_file))
            sequences_complex = {}
            if proteins.has_key(pair[0]) and proteins.has_key(pair[1]):
             sequences_complex.setdefault("A",proteins.get(pair[0]))
             sequences_complex.setdefault("B",proteins.get(pair[1]))
             try:
                pdb_renumber=PDB()
                pdb_renumber=renumber_pdb(config,path_file,output_file,sequences_complex,dummy_dir)
                if len(pdb_renumber.chain_identifiers)>1:
                   pdb_renumber.write(output_file,force=True)
                else:
                   if verbose: sys.stdout.write("\t\t\t-- Failed to renumerate %s\n"%output_file)
             except Exception as e:
                sys.stderr.write("WARNING %s\n"%e)
            else:
             sys.stdout.write('\t\t\t-- Skip renumbering of %s (missing FastA sequences %s %s)\n'%(input_file,pair[0],pair[1]))
        fd.close()

    # Informe on the files parsed
    if options.show:
     for pair,folder in new_models.iteritems():
        sys.stdout.write("Parse: %s \n"%(folder))

    # Group files with similar interfaces
    # Defined criteria from config file
    potential=config.get('Parameters','PPI_threshold_type')
    radius   =float(config.get('Parameters','PPI_distance_threshold_shell'))
    overlap  =float(config.get('Parameters','overlap_interface'))
    #Obtain the name of the wild-type
    wt=(None,None)
    wti=(None,None)
    for a,b in new_models.iterkeys():
      wt_A = a.split("_")[0]
      wt_B = b.split("_")[0]
      wt=(wt_A,wt_B)
      wti=(wt_B,wt_A)
    if not new_models.has_key(wt) and new_models.has_key(wti): wt=wti
    #If there is no wild type select the pair with the biggest cluster to use as basic core in substitution of the wild type
    if not new_models.has_key(wt):
     maxsize=0
     for pair,folder in new_models.iteritems():
       try:
         pdbdir,interface,cluster = filter_homologs(folder,potential=potential,radius=radius,overlap=overlap,show=False)
       except IOError as e:
         print e
         continue
       if max([len(s) for c,s in sorted(cluster.items(),key=lambda x: len(x[1]),reverse=True)]) >maxsize: wt=pair
    if wt==(None,None):
      return group_models 
    #Informe on the selection if "verbose"   
    if options.show:
      sys.stdout.write("\t\t -- Clustering similar interfaces for wild-type PPI %s - %s ...\n"%(wt_A,wt_B))

    #Rank the clusters of Wild Type (or select pair) by size and create the lists with grouped structures
    #Keep in the dictionary a list with the name of files (one for each pose) listing the models
    folder_wt = new_models[wt]
    try:
     pdbdir_wt,interface,cluster_wt = filter_homologs(folder_wt,potential=potential,radius=radius,overlap=overlap,show=False)
    except IOError as e:
     print e
    if options.show: sys.stdout.write("\t\t\t Directory %s\n"%pdbdir_wt)

    #Write a list of non redundant interfaces for wt
    fg=open(os.path.join(pdbdir_wt,options.label+"_"+wt[0]+"::"+wt[1]+"_nr_interfaces.list"),"w")
    #Pose 0 is reserved for all structures wityhout clustering, next steps pose will increase by 1
    pose=0
    pose_wt=[]
    for c,s in sorted(cluster_wt.items(),key=lambda x: len(x[1]),reverse=True): 
        if len(s) <=0: continue
        pose     = pose + 1
        pose_list= os.path.join(pdbdir_wt,options.label+"_"+wt[0]+"::"+wt[1]+"_cluster_"+str(pose)+".list")
        if options.show: sys.stdout.write("\t\t\t Seed cluster %s of wild-type is in POSE %d\n"%(c,pose))
        pose_wt.append(c)
        group_models.setdefault(wt,[]).append((pose,pose_list))
        fd=open(pose_list,"w")
        grouped="\t"
        for p in s:
          pdb=pdbdir_wt+"/"+p
          if not fileExist(pdb): continue
          fd.write("%s\n"%(pdb))
          if p != c:
            z=0.0
            if len(interface[c]) >0:
              z=100 * float(len(interface[c].intersection(interface[p])))/len(interface[c])
            pb=p.split(".")
            grouped += ".".join([str(x) for x in pb[1:-1]])
            grouped += " ({0:5.1f})\t".format(z)
        fd.close()
        fg.write("%s\t%s\n"%(pdbdir_wt+"/"+c,grouped))
    fg.close()
    pose_last=pose

    #Compare and group similar poses to the wt for the other pairs
    #the remaining poses are ranked at the end by the size of the clusters
    min_ratio_common_poses=float(config.get('Parameters', 'min_ratio_common_poses'))
    for pair,folder in new_models.iteritems():
       if pair==wt: continue

       if options.show: sys.stdout.write("\t\t -- Clustering similar interfaces for PPI %s - %s ...\n"%(pair[0],pair[1]))
       try:
         pdbdir,interface,cluster = filter_homologs(folder,potential=potential,radius=radius,overlap=overlap,show=False)
       except IOError as e:
         print e
         continue
       if options.show: sys.stdout.write("\t\t\t Directory %s\n"%pdbdir)

       #Select the poses of wt that correspond to each cluster of similar interfaces
       select={}
       for c,s in sorted(cluster.items(),key=lambda x: len(x[1]),reverse=True):
        if len(s) <=0: continue
        pose_pair=None
        maxcommon=0
        ratio_common=0
        for ipose in xrange(1,len(pose_wt)+1):
         if ipose in select.itervalues(): continue
         s_wt=cluster_wt[pose_wt[ipose-1]]
         main_wt,main_s=compare_templates(s_wt,s)
         if len(main_wt.intersection(main_s))>maxcommon: 
            pose_pair=ipose
            maxcommon=len(main_wt.intersection(main_s))
            ratio_common=float(maxcommon)/len(main_wt)
           # if options.show: sys.stdout.write("\t\t\t Ratio of common poses %f (%d/%d) with pose %d \n"%(ratio_common,maxcommon,len(main_wt),ipose))
            if maxcommon > len(main_wt):
               if options.show:sys.stdout.write("Error too large number of common poses %s compared with %s\n"%(repr(main_wt.intersection(main_s)),repr(main_wt)))
        if pose_pair is not None and ratio_common > min_ratio_common_poses:  
           select.setdefault(c,pose_pair)
           if options.show: sys.stdout.write("\t\t\t Seed cluster %s in wild-type POSE %d (ratio %f)\n"%(c,pose_pair,ratio_common))
        
       #To continue other poses after the last one use the iterator pose
       pose=pose_last

       #Write a list of non redundant interfaces for  the pair
       fg=open(os.path.join(pdbdir,options.label+"_"+pair[0]+"::"+pair[1]+"_nr_interfaces.list"),"w")

       for c,s in sorted(cluster.items(),key=lambda x: len(x[1]),reverse=True):
        if len(s) <=0: continue

        #Print the list of the selected pose and add to the list in group_models dictionary
        if select.has_key(c):
         pose_pair=select[c]
         pose_list= os.path.join(pdbdir,options.label+"_"+pair[0]+"::"+pair[1]+"_cluster_"+str(pose_pair)+".list")
         group_models.setdefault(pair,[]).append((pose_pair,pose_list))
         fd=open(pose_list,"w")
         grouped="\t"
         for p in s:
          pdb=pdbdir+"/"+p
          if not fileExist(pdb): continue
          fd.write("%s\n"%(pdb))
          if p != c:
            z=0.0
            if len(interface[c]) >0:
              z=100 * float(len(interface[c].intersection(interface[p])))/len(interface[c])
            pb=p.split(".")
            grouped += ".".join([str(x) for x in pb[1:-1]])
            grouped += " ({0:5.1f})\t".format(z)
         fd.close()
         fg.write("%s\t%s\n"%(pdbdir+"/"+c,grouped))

        #Print the list for a new pose and add to the list in group_models dictionary
        else:
         pose     = pose + 1
         pose_list= os.path.join(pdbdir,options.label+"_"+pair[0]+"::"+pair[1]+"_cluster_"+str(pose)+".list")
         group_models.setdefault(pair,[]).append((pose,pose_list))
         if options.show: sys.stdout.write("\t\t\t Seed cluster %s in new POSE %d \n"%(c,pose))
         fd=open(pose_list,"w")
         grouped="\t"
         for p in s:
          pdb=pdbdir+"/"+p
          if not fileExist(pdb): continue
          fd.write("%s\n"%(pdb))
          if p != c:
            z=0.0
            if len(interface[c]) >0:
              z=100 * float(len(interface[c].intersection(interface[p])))/len(interface[c])
            pb=p.split(".")
            grouped += ".".join([str(x) for x in pb[1:-1]])
            grouped += " ({0:5.1f})\t".format(z)
         fd.close()
         fg.write("%s\t%s\n"%(pdbdir+"/"+c,grouped)) 
       fg.close()
       #Keep the number of last pose for the analysis of next pairs
       pose_last=pose

    return group_models

def interface_affected(pair):

    x,y=pair
# x,y are the ratio of interfaces with different sequence with respect to a wild_type form interaction of proteins a and b, espectively

    ratio_a = int(10*x)
    ratio_b = int(10*y)
# version for 2 ratios
#    return "_i"+str(ratio_a)+"&"+str(ratio_b)
# version single value
    return "_i"+str(max([ratio_a,ratio_b]))

def compare_interface_sequences(list_of_models,dummy_dir):

#Result dicitionary
    pose_interfaces={}
#Define parameters to calculate the interface
    PPI_threshold_type = config.get('Parameters', 'PPI_threshold_type')
    PPI_distance_threshold = float(config.get('Parameters', 'PPI_distance_threshold'))

#Select the name of the wild-type pair of sequences for each pair
#wt_of_pair is a dictionary that assigns for each pair the corresponding wild-type
#indicating if the order of proteins is the same "+1" or inverse "-1", or if there is no wild-type to compare "0"
    wt_of_pair={}
    for pair,cluster_list in list_of_models.iteritems():
      a,b = pair
      wt_A = a.split("_")[0]
      wt_B = b.split("_")[0]
      wt=(wt_A,wt_B)
      wti=(wt_B,wt_A)
      if list_of_models.has_key(wt):
        wt_of_pair.setdefault(pair,(wt,1))
      elif list_of_models.has_key(wti):
        wt_of_pair.setdefault(pair,(wti,-1))
      else:
        wt_of_pair.setdefault(pair,(pair,0))

#Get the list of pair of reference-sequences of each pose and their wild-type references

    main_sequences,main_sequences_wt = extract_main_sequence(list_of_models,wt_of_pair)

# Chek sequences of different forms with respect to wild-type fo each pose
# defined names of use
# a,b is the name of the pair of interacting proteins
# pose is the cluster number to be checked
# seq_a and seq_b are the main sequences of a and b for the corresponding pose
# pa and pb are the ratio of occurrence of the main sequences used for a and b in the  corresponding pose
# wt_a,seq_wt_a,pwa and  wt_b,seq_wt_b,pwb are the names sequecnes and ratio of occurrence of the wild-type
#   protein names corresponding in the pose for a and b (this can change in diferent poses, not the name but the main sequence)
# check_A is the set of 2-tuples aa,residue_number different of sequence A in pose respect to the main wild-type sequence
# check_wt_A is the set of 2-tuples aa,residue_number different of wild type sequence of A respect to A in pose
# idem for check_B and check_wt_B
# Information is duplicated, with check_wt_A and check_wt_B is enough
#
    check_sequence={}
    for ppairs,spairs in main_sequences.iteritems():
        pose,pair      = ppairs
        a,b            = pair
        sp_a,sp_b      = spairs
        seq_a,pa       = sp_a
        seq_b,pb       = sp_b
        wt,order       = wt_of_pair[(a,b)]
        if not main_sequences_wt.has_key((pose,wt)): continue
        s_wt_a,s_wt_b  = main_sequences_wt[(pose,wt)]
        if order > 0:
          wt_a,wt_b    = wt
          seq_wt_a,pwa = s_wt_a
          seq_wt_b,pwb = s_wt_b
        if order < 0:
          wt_b,wt_a    = wt
          seq_wt_a,pwa = s_wt_b
          seq_wt_b,pwb = s_wt_a
        try:
         check_a, check_wt_A = compare_seq_wt(a,seq_a,wt_a,seq_wt_a,dummy_dir)
         check_b, check_wt_B = compare_seq_wt(b,seq_b,wt_b,seq_wt_b,dummy_dir)
         check_sequence.setdefault(ppairs,(check_a, check_wt_A, check_b, check_wt_B))
        except:
         sys.stderr.write("WARNING: Skip pair %s %s\n"%(a,b))

# Get interface aminoacids of Wild-type for each pose
# interface is a dictionary with the pose and the pair as key, 
# the values are 3-tuples: two Counter comprising the data of the lists of Aa in the interface and the number of PDB files in the cluster
    interface_wt={}
    for ppairs,spairs in main_sequences_wt.iteritems():
        pose,pair      = ppairs
        a,b            = pair
        sp_a,sp_b      = spairs
        seq_a,pa       = sp_a
        seq_b,pb       = sp_b
        cluster_list_wt= list_of_models[pair]
        for pose_wt,pose_list_wt in cluster_list_wt:
          if pose!=pose_wt: continue
          if interface_wt.has_key((pair,pose_wt)): continue
          filelist=open(pose_list_wt,"r")
          interface_a=[]
          interface_b=[]
          n=0
          for line in  filelist:
            n=n+1
            pdb_file=line.strip()
            if not fileExist(pdb_file):continue
            pdb=PDB(pdb_file)
            chains=[x for x in pdb.chain_identifiers]
            a_index={}
            b_index={}
            p=pdb.get_chain_by_id(chains[0])
            for k,idx in enumerate(p.protein_idx.split(";")): a_index.setdefault(idx.replace(" ",""),k)
            q=pdb.get_chain_by_id(chains[1])
            for k,idx in enumerate(q.protein_idx.split(";")): b_index.setdefault(idx.replace(" ",""),k)
            #check these are main-sequences for analysis
            if p.protein_sequence != seq_a or q.protein_sequence != seq_b : continue
            # Check contacts
            protein_complex = Complex(pdb, PPI_type = PPI_threshold_type, PPI_distance = PPI_distance_threshold)
            for interaction in  protein_complex.PPInterfaces:
                 order_complex=0
                 if interaction.interactor_id.split("_")[-1] == chains[1] and interaction.protein_id.split("_")[-1] == chains[0]: order_complex=1
                 if interaction.interactor_id.split("_")[-1] == chains[0] and interaction.protein_id.split("_")[-1] == chains[1]: order_complex=-1
                 if order_complex==0: continue
                 if order_complex>0:
                    for aa in interaction.interactor_positions: interface_b.append((aa.single_letter, b_index[aa.identifier.replace(" ","")]))
                    for aa in interaction.protein_positions: interface_a.append((aa.single_letter, a_index[aa.identifier.replace(" ","")]))
                 if order_complex<0:
                    for aa in interaction.interactor_positions: interface_a.append((aa.single_letter, a_index[aa.identifier.replace(" ","")]))
                    for aa in interaction.protein_positions: interface_b.append((aa.single_letter, b_index[aa.identifier.replace(" ","")]))
          #After checking all the structures in the list
          filelist.close()
          #Add the counter dictionary of the list of pairs Aa,residue_number of the interfaces between a and b chains 
          interface_wt.setdefault(ppairs,(Counter(interface_a),Counter(interface_b),n))

#Check if residues involved in the differences with respect to wild-type are also involved in the interface of the wt                   
    
    for pair,cluster_list in list_of_models.iteritems():
      a,b = pair
      wt,order = wt_of_pair[(a,b)]
      #check only pairs with known wild-type
      if order == 0: continue  
      for pose,pose_list in cluster_list:
         ppair     =(pose,pair)
         wt_pair   =(pose,wt)
         #check exist wild-type interface for this pose
         if not interface_wt.has_key(wt_pair): continue
         if not main_sequences.has_key(ppair): continue
         sp_a,sp_b = main_sequences[ppair]
         seq_a,ra  = sp_a
         seq_b,rb  = sp_b
         check_a, check_wt_A, check_b, check_wt_B         = check_sequence[ppair]
         interface_wt_A, interface_wt_B, number_wt_models = interface_wt[wt_pair]
         ratio_A=0
         ratio_B=0
         if number_wt_models > 0:
           affected_residue_A=check_wt_A.intersection(interface_wt_A.keys())
           if len(affected_residue_A) >0:
             largest_number_A = max([interface_wt_A[x] for x in affected_residue_A])
             ratio_A = float(largest_number_A) / number_wt_models
           affected_residue_B=check_wt_B.intersection(interface_wt_B.keys())
           if len(affected_residue_B)>0:
             largest_number_B = max([interface_wt_B[x] for x in affected_residue_B])
             ratio_B = float(largest_number_B) / number_wt_models
           pose_interfaces.setdefault(ppair,(ratio_A,ratio_B))
#return result
    return pose_interfaces

def compare_seq_wt(name,seq,name_wt,seq_wt,dummy_dir):
    clustal_exe  = os.path.join(config.get('Paths','clustal_path'),'clustalw2')
    #create files

    infile=dummy_dir+"/tmp_"+name_wt+"_"+name+".fa"
    outfile=dummy_dir+"/tmp_"+name_wt+"_"+name+".aln"
    dndfile=dummy_dir+"/tmp_"+name_wt+"_"+name+".dnd"
    fd=open(infile,"w")
    fd.write(">WT{0:s}\n{1:s}\n".format(name_wt,seq_wt))
    fd.write(">{0:s}\n{1:s}\n".format(name,seq))
    fd.close()

    try:
      # run clustalw2
      msa_cline=ClustalwCommandline(clustal_exe,infile=infile,outfile=outfile)
      child = subprocess.Popen(str(msa_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell="/bin/bash")
      child.communicate()
      #store alignment in compare
      compare=AlignIO.read(outfile,'clustal')
    except Exception as e:
      sys.stderr.write("ERROR: %s\n"%e)
      return (set(),set())

    #remove temporary fasta and alignment files
    remove_files([infile,outfile,dndfile])

    wt_msa =compare[0].seq
    seq_msa=compare[1].seq
    try:
      len_wt=len(seq_wt)
      len_msa=len(wt_msa)
    except  Exception as e:
      sys.stderr.write("ERROR and SKIP: %s\n"%e)
      return (set(),set())
    sys.stdout.write("\t\t -- Sequence comparison of Wild-Type %s versus the form %s\n"%(name_wt,name))
    difference=""
    for ii in itertools.izip(wt_msa,seq_msa):
      if ii[0]==ii[1]:difference=difference+" "
      else:           difference=difference+"*"
    sys.stdout.write("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n"%(wt_msa,seq_msa,difference))
    #Check sequence affected of wild-type
    residue_number_wt=0
    k=0
    affected_wt_residue=set()
    while (residue_number_wt<len(seq_wt)) and (k<len(wt_msa)):
      if wt_msa[k] == seq_wt[residue_number_wt] and seq_msa[k]==wt_msa[k]:
         residue_number_wt = residue_number_wt + 1
         k = k + 1
      elif  wt_msa[k] == "-":
         affected_wt_residue.add((seq_wt[residue_number_wt],residue_number_wt))
         if residue_number_wt>0: affected_wt_residue.add((seq_wt[residue_number_wt-1],residue_number_wt-1))
         k = k + 1
      elif  seq_msa[k] !=  wt_msa[k] and wt_msa[k] == seq_wt[residue_number_wt]:
         affected_wt_residue.add((seq_wt[residue_number_wt],residue_number_wt))
         residue_number_wt = residue_number_wt + 1
         k = k + 1
      else:
         print("Error: no matching sequences %s [%d] vs %s [%d]\n"%(wt_msa[k],k,seq_wt[residue_number_wt],residue_number_wt))
         residue_number_wt = residue_number_wt + 1
         k = k + 1

    #Check sequence affected of sequence under test
    residue_number=0
    k=0
    affected_residue=set()
    try:
      len_seq=len(seq)
      len_msa=len(seq_msa)
    except  Exception as e:
      sys.stderr.write("ERROR and SKIP: %s\n"%e)
      return (set(),set())

    while (residue_number<len(seq)) and (k<len(seq_msa)):
      if seq_msa[k] == seq[residue_number] and seq_msa[k]==wt_msa[k]:
         residue_number = residue_number + 1
         k = k + 1
      elif  seq_msa[k] == "-":
         affected_residue.add((seq[residue_number],residue_number))
         if residue_number>0: affected_residue.add((seq[residue_number-1],residue_number-1))
         k = k + 1
      elif  seq_msa[k] != wt_msa[k] and  seq_msa[k] == seq[residue_number]:
         affected_residue.add((seq[residue_number],residue_number))
         residue_number = residue_number + 1
         k = k + 1
      else:
         print("Error: no matching sequences %s [%d] vs %s [%d]\n"%(seq_msa[k],k,seq[residue_number],residue_number))
         residue_number_wt = residue_number_wt + 1
         k = k + 1
        
    return (affected_residue,affected_wt_residue)
     



def extract_main_sequence(list_of_models,wt_of_pair):

    main_sequences={}
    main_sequences_wt={}
    for pair,cluster_list in list_of_models.iteritems():
        a,b = pair
        for pose,pose_list in cluster_list:
            sequence={}
            seq_a=None
            seq_b=None
            ratio_a=0.0
            ratio_b=0.0
            filelist=open(pose_list,"r")
            for line in  filelist:
              pdb_file=line.strip()
              if not fileExist(pdb_file):continue
              try:
                pdb=PDB(pdb_file)
                chains=[x for x in pdb.chain_identifiers]
                p=pdb.get_chain_by_id(chains[0])
                q=pdb.get_chain_by_id(chains[1])
                sequence.setdefault(a,[]).append(p.protein_sequence)
                sequence.setdefault(b,[]).append(q.protein_sequence)
              except Exception as e:
                sys.stdout.write("WARNING: file not found, skip %s \n"%pdb_file)
                sys.stderr.write("ERROR: %s\n"%e)
            filelist.close()
            if sequence.has_key(a): 
              seq_a,times_a=Counter(sequence[a]).most_common(1)[0]
              ratio_a=float(times_a)/len(sequence[a])
            if sequence.has_key(b): 
              seq_b,times_b=Counter(sequence[b]).most_common(1)[0]
              ratio_b=float(times_b)/len(sequence[b])
            main_sequences.setdefault((pose,pair),((seq_a,ratio_a),(seq_b,ratio_b)))
            wt,order= wt_of_pair[pair]
            if pair == wt:
                  main_sequences_wt.setdefault((pose,wt),((seq_a,ratio_a),(seq_b,ratio_b)))
    main_sequences_all=( main_sequences , main_sequences_wt)
    return main_sequences_all

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = 'Energy analysis of protein-protein interactions modelled with ModPPI',
        epilog      = '@Oliva\'s lab 2016')
    parser.add_argument('-ppi', '--ppi_list', dest = 'ppi_list', action = 'store',default=None,
                        help = 'File with models done (i.e. "interactions_done.list" )')
    parser.add_argument('-seq', '--fasta', dest = 'fasta_file', action = 'store',
                        help = 'Input file with a list of FASTA sequences of the proteins to test')
    parser.add_argument('-l', '--label_name', dest = 'label', action = 'store', default = "ModPPI",
                        help = 'Label to store the analyses (default is "ModPPI")')
    parser.add_argument('-o', '--output_dir', dest = 'outdir', action = 'store', default = None,
                        help = 'Directory to store the results (default is for global outputs "analysis" and for each model the same as the directory with the 3D-models')
    parser.add_argument('-d', '--dummy_dir', dest = 'dummy_dir', action = 'store', default = '/tmp/modppi_dummy',
                        help = 'Specifies the location of the dummy folder (default is /tmp/modppi_dummy)')
    parser.add_argument('-zrank', '--zrank', dest = 'zrank', action = 'store_true',
                        help = 'Flag to calculate energies with ZRANK (default is False)')
    parser.add_argument('-foldx', '--foldx', dest = 'foldx', action = 'store_true',
                        help = 'Flag to calculate energies with FoldX (default is False)')
    parser.add_argument('-rosetta', '--rosetta', dest = 'rosetta', action = 'store_true',
                        help = 'Flag to calculate energies with ROSETTA (default is False)')
    parser.add_argument('-ssp', '--split-potentials', dest = 'ssp', action = 'store_true',
                        help = 'Flag to calculate energies with Split-Statistic Potentials (default is False unless no other method is used)')
    parser.add_argument('-hydro','--hydrogens', dest = 'hbplus', action = 'store_true',
                        help = 'Flag to include hydrogens (default is False)')
    parser.add_argument('-boxplot','--boxplot', dest = 'boxplot', action = 'store_true',
                        help = 'Flag to make default boxplots (default is False)')
    parser.add_argument('-r','--renumerate', dest = 'renumerate' , action = 'store_true',
                        help = 'Flag to renumber the sequences as in the original FastA (default is False)')
    parser.add_argument('-v', '--verbose', dest = 'show', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
    options = parser.parse_args()
    return options


if __name__ == '__main__':
    main()


