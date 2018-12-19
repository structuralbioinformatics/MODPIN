import ConfigParser
import sys
import os 
from os import listdir
from os.path import isfile, join
import fnmatch

def main():
    
	print ("* ---- Building the config file ---- *") 


	script_dir = os.path.dirname(__file__)
	rel_path = os.path.join("scripts","config.ini")
	abs_file_path = os.path.join(script_dir, rel_path)
	config = ConfigParser.ConfigParser()
	config.read(abs_file_path)
	
	modppi_path =  os.path.dirname(os.path.abspath( __file__ ))
	python_path = sys.executable

	config.set('Paths', 'python_path', python_path)
	config.set('Paths', 'modppi_path', modppi_path)
	config.set('Paths', 'sbi_library_path', os.path.join(".", "src"))
	config.set('Paths', 'functions_path', os.path.join(".", "src", "functions" ))
	config.set('Paths', 'archdb_path',  os.path.join(modppi_path, "src"))
 
	config.set('Parameters','PPI_distance_threshold','8.0')
	config.set('Parameters','PPI_distance_threshold_shell','15.0')
	config.set('Parameters','PPI_threshold_type','cb')
	config.set('Parameters','overlap_interface','35.0')
	config.set('Parameters','min_ratio_common_poses','0.49')
	config.set('Parameters','n_parameter_twilight_zone','0')
	config.set('Parameters','extension_threshold','15')
	config.set('Parameters','hydrogens','full')
	config.set('Parameters','relax','yes')
	config.set('Parameters','ssp_score','ES3DC; ZES3DC; PAIR; ZPAIR; COMB; ZCOMB')
	config.set('Parameters','ssp_type','CB')



	with open(abs_file_path, 'w') as configfile:    
			config.write(configfile)	

	configuration = raw_input ("Press enter to continue with this configuration or enter \"n\" if you prefer to exit and edit config.ini file manually: ")
	if configuration== "n":
		print ("Please edit the config.ini file in the scripts directory of your Modpin installation folder.")
	else:	
                configuration_data = raw_input ("Press enter to configure your own storage of data \"n\" if you prefer to use the default configuration: ")

                if configuration_data =="n" or configuration_data =="N":
       			try:
          			if not os.path.isdir(os.path.join(modppi_path,"data")): os.makedirs(os.path.join(modppi_path,"data"))
          			if not os.path.isdir(os.path.join(modppi_path,"blast")):os.makedirs(os.path.join(modppi_path,"blast"))
        		except:
          			print "ERROR: You don't have permissions in this folder"
          			exit(0)
			config.set('Paths', 'data_path', os.path.join(".", "data"))
			config.set('Paths', 'blast_dir',  os.path.join(".", "blast"))
			config.set('Paths', 'pdb_path',  os.path.join(".", "data", "pdbsource"))
			config.set('Paths', '3did_path',  os.path.join(".", "data", "3did"))
			config.set('Files', 'fasta_list_file', os.path.join(".", "data", "nr.fa"))
			config.set('Files', 'nr90_list_file', os.path.join(".", "data", "merged_nr90.fa"))
			config.set('Files', 'ppi_files', os.path.join(".", "data", "pdb"))
			config.set('Files', 'ddi_files', os.path.join(".", "data", "3did"))
			config.set('Files', 'pdb_list_file', os.path.join(".", "data", "pdb.list"))
			config.set('Files', '3did_list_file', os.path.join(".", "data", "3did.list"))
			config.set('Files', 'database_file', os.path.join(".", "database.fasta"))
		else:
		 	data_path = raw_input("STRUCTURE) Insert the path to structural data files (should have pdbsource and 3did): ")
                	data_path = data_path.lstrip()
		 	while not find_file(data_path,"pdb*") or not  find_file(data_path,"3did"):
                    	  data_path = raw_input("Insert the correct path to structural data files (pdbsource or 3did are not found), or press enter to skip and add it later manually to the config.ini file: ")
                   	  data_path = data_path.lstrip()
                   	  if not data_path:
                   	     break;
                   	  else:
			     continue

                 	if find_file(data_path,"pdb*") and find_file(data_path,"3did"):
                    	  config.set('Paths', 'data_path',data_path)
		 	  config.set('Paths', 'pdb_path',  os.path.join(data_path, "pdbsource"))
		 	  config.set('Paths', '3did_path',  os.path.join(data_path, "3did"))
			  config.set('Files', 'ppi_files', os.path.join(data_path, "pdb"))
			  config.set('Files', 'ddi_files', os.path.join(data_path, "3did"))

		 	blast_dir = raw_input("SEQUENCES) Insert the path to sequence data formatted by BLAST/2.2.26: ")
			blast_dir = blast_dir.lstrip()
		 	while not os.path.isdir(blast_dir):
                    	  blast_dir = raw_input("Insert the correct path to sequence data files (folder not found), or press enter to skip and add it later manually to the config.ini file: ")
                   	  blast_dir = blast_dir.lstrip()
                   	  if not blast_dir:
                   	     break;
                 	  else:
			     continue

                	if os.path.isdir(blast_dir):
                 	    config.set('Paths', 'blast_dir',blast_dir)
			
                        database_fasta = "database.fasta"
                        while not find_file(blast_dir,database_fasta):
                          database_fasta = raw_input("Insert the name of the sequence database ('database.fasta' by default),  or press enter to skip and add it later manually to the config.ini file: ")
                          database_fasta = database_fasta.lstrip()
                          if not database_fasta:
                   	     break;
                 	  else:
			     continue

                 	if find_file(blast_dir,database_fasta):
			  config.set('Files', 'database_file', os.path.join(".", database_fasta))

                        fasta_list_file=os.path.join(data_path,"nr.fa")
		        while not os.path.isfile(fasta_list_file):
				fasta_list_file = raw_input("Insert the correct file of non-redundant sequences 'nr.fa', or press enter if it has to be created: ")
				fasta_list_file = fasta_list_file.lstrip()
				if not fasta_list_file: 
					break
				else:
					continue
                        if os.path.isfile(fasta_list_file):
				config.set('Files', 'fasta_list_file',fasta_list_file)

                        nr90_list_file=os.path.join(data_path,"merged_nr90.fa")
		        while not os.path.isfile(nr90_list_file):
				nr90_list_file = raw_input("Insert the correct file of non-redundant sequences merging PDB and 3DID 'merged_nr90.fa', or press enter if it has to be created: ")
				nr90_list_file = nr90_list_file.lstrip()
				if not nr90_list_file: 
					break
				else:
					continue
                        if os.path.isfile(nr90_list_file):
				config.set('Files', 'nr90_list_file',nr90_list_file)
                            
                        pdb_list_file=os.path.join(data_path,"pdb.list")
		        while not os.path.isfile(pdb_list_file):
				pdb_list_file = raw_input("Insert the correct file of list of PDB 'pdb.list', or press enter if it has to be created: ")
				pdb_list_file = pdb_list_file.lstrip()
				if not pdb_list_file: 
					break
				else:
					continue
                        if os.path.isfile(pdb_list_file):
				config.set('Files', 'pdb_list_file',pdb_list_file)


                        did_list_file=os.path.join(data_path,"3did.list")
		        while not os.path.isfile(did_list_file):
				did_list_file = raw_input("Insert the correct file of list of 3DID '3did.list', or press enter if it has to be created: ")
				did_list_file = did_list_file.lstrip()
				if not did_list_file: 
					break
				else:
					continue
                        if os.path.isfile(did_list_file):
				config.set('Files', '3did_list_file',did_list_file)
                            
                            



		modeller_path = raw_input("1) Insert path to bin folder of Modeller: ")
		modeller_path=modeller_path.lstrip()
		while not find_file(modeller_path, "mod*.*"):
			modeller_path = raw_input("Insert the correct path to bin folder of Modeller, or press enter to skip and add it later manually to the config.ini file: ")
			modeller_path = modeller_path.lstrip()
			if not modeller_path:
				break
			else:
				continue
		if find_file(modeller_path, "mod*.*"):
			config.set('Paths', 'modeller_path', modeller_path)

		blast_path = raw_input("2) Insert path to bin folder of BLAST/2.2.26: ")
		blast_path =  blast_path.lstrip()
		while which(os.path.join(blast_path, "blastpgp"))==False or os.path.isdir(blast_path)==False:
			blast_path = raw_input ("Insert the correct path to executable of BLAST, or press enter to skip and add it later manually to the config.ini file: ")	
			blast_path =  blast_path.lstrip()
			if not blast_path:
				break
			else:
				continue
		if which(os.path.join(blast_path, "blastpgp"))==True:
			config.set('Paths', 'blast_path', blast_path) 
	
		zrank_path = raw_input("3) Insert path to executable of ZRANK: ")
		zrank_path = zrank_path.lstrip()
		while os.path.isfile(zrank_path) == False:
			zrank_path  = raw_input ("Insert the correct path to executable of ZRANK, or press enter to skip and add it later manually to the config.ini file: ")
			zrank_path = zrank_path.lstrip()
			if not zrank_path :
				break
			else:
				continue
		if  os.path.isfile(zrank_path)==True:
			config.set('Paths', 'zrank_path', zrank_path)
	
		foldx_path = raw_input("4) Insert path to executable of FOLDX: ")
		foldx_path = foldx_path.lstrip()
		while os.path.isfile(foldx_path) == False:
			foldx_path  = raw_input ("Insert the correct path to executable of FoldX, or press enter to skip and add it later manually to the config.ini file: ")
			foldx_path = foldx_path.lstrip()
			if not foldx_path :
				break
			else:
				continue
		if  os.path.isfile(os.path.join(foldx_path,"foldx")==True and os.path.isfile(os.path.join(foldx_path,"rotabase.txt")==True:
			config.set('Paths', 'foldx_path', foldx_path)
                else:
                        sys.stdout.write("\t-- Warning!!! Neither 'foldx' nor 'rotabase.txt' were found under FoldX path\n")

		clustal_path = raw_input("5) Insert path to bin folder of clustalw2: ")
		clustal_path = clustal_path.lstrip()
		while which(os.path.join(clustal_path, "clustalw2")) == False or os.path.isdir(clustal_path) == False:
			clustal_path = raw_input ("Insert the correct path to executable of clustalw2, or press enter to skip and add it later manually to the config.ini file: ")
			clustal_path = clustal_path.lstrip()
			if not clustal_path:
				break
			else:
				continue
		if  which(os.path.join(clustal_path, "clustalw2")) == True:
			config.set('Paths', 'clustal_path', clustal_path)

		rosetta_path = raw_input("6) Insert path to the main folder Rosetta: ")
		rosetta_path = rosetta_path.lstrip()
		while os.path.isdir(os.path.join(rosetta_path, "database"))  == False:
			rosetta_path = raw_input ("Insert the correct path to Rosetta main folder which contains the database, or press enter to skip and add it later manually to the config.ini file: ")		
			rosetta_path = rosetta_path.lstrip()
			if not rosetta_path:
				break
			else:
				continue

		if  os.path.isdir(os.path.join(rosetta_path, "database"))  == True:
			config.set('Paths', 'rosetta_path', rosetta_path)

			relax_exe = os.path.join(rosetta_path, "source", "bin", "fixbb.linuxgccrelease")
                        if find_file( os.path.join(rosetta_path, "source", "bin"),"fixbb*"): 
                           sys.stdout.write("\t-- Found fixbb, default to be used for relaxing the structures: %s\n"%found_file(os.path.join(rosetta_path, "source", "bin"),"fixbb*"))
                           relax_exe = os.path.join(rosetta_path, "source", "bin", found_file(os.path.join(rosetta_path, "source", "bin"),"fixbb*"))
                        elif find_file( os.path.join(rosetta_path,"bin"),"fixbb*"):
                           sys.stdout.write("\t-- Found fixbb, default to be used for relaxing the structures: %s\n"%found_file(os.path.join(rosetta_path,"bin"),"fixbb*"))
                           relax_exe = os.path.join(rosetta_path, "bin", found_file(os.path.join(rosetta_path, "bin"),"fixbb*"))
                        elif find_file( os.path.join(rosetta_path,"source","bin"),"relax*"):
                           sys.stdout.write("\t-- Found relax, default to be used for relaxing the structures: %s\n"%found_file(os.path.join(rosetta_path, "source", "bin"),"relax*"))
                           relax_exe = os.path.join(rosetta_path, "source", "bin", found_file(os.path.join(rosetta_path, "source","bin"),"relax*"))
                        elif find_file( os.path.join(rosetta_path,"bin"),"relax*"):
                           sys.stdout.write("\t-- Found relax, default to be used for relaxing the structures: %s\n"%found_file(os.path.join(rosetta_path, "bin"),"relax*"))
                           relax_exe = os.path.join(rosetta_path, "bin", found_file(os.path.join(rosetta_path, "bin"),"relax*"))
			else:
                           sys.stdout.write("\t-- Neither FIXBB nor RELAX executables were found under ROSETTA path\n")
                              
			if os.path.isfile(relax_exe)== True:
				config.set('Paths', 'relax_exe', relax_exe)
			else: 
				while os.path.isfile(relax_exe)  == False:
					relax_exe = raw_input (" -Relaxing backbone-design application fixbb.linuxgccrelease (or relax.linuxgccrelease) not found under Rosetta bin folder. \n  Insert the correct path wih the name of executable, or press enter to skip and add it later manually to the config.ini file: ")
					if not relax_exe:
						break
					else:
						continue 
			if os.path.isfile(relax_exe)== True:
				config.set('Paths', 'relax_exe', relax_exe)


			interface_analyzer = os.path.join(rosetta_path, "source", "bin", "InterfaceAnalyzer.linuxgccrelease")
                        if find_file( os.path.join(rosetta_path, "source", "bin"),"InterfaceAnalyzer*"):
                          sys.stdout.write("\t-- Found InterfaceAnalyzer: %s\n"%found_file(os.path.join(rosetta_path, "source", "bin"),"InterfaceAnalyzer*"))
			  interface_analyzer = os.path.join(rosetta_path, "source", "bin", found_file(os.path.join(rosetta_path, "source", "bin"),"InterfaceAnalyzer*"))
                        elif find_file(os.path.join(rosetta_path, "bin"),"InterfaceAnalyzer*"):
                          sys.stdout.write("\t-- Found InterfaceAnalyzer: %s\n"%found_file(os.path.join(rosetta_path,  "bin"),"InterfaceAnalyzer*"))
			  interface_analyzer = os.path.join(rosetta_path, "bin", found_file(os.path.join(rosetta_path,  "bin"),"InterfaceAnalyzer*"))
			else:
                           sys.stdout.write("\t-- INTERFACEANALYZER executable was not found under ROSETTA path\n")

			if os.path.isfile(interface_analyzer)== True:
				config.set('Paths', 'interface_analyzer', interface_analyzer)
			else: 
				while os.path.isfile(interface_analyzer)  == False:
					interface_analyzer = raw_input(" -Interface-analyzer application InterfaceAnalyzer.linuxgccrelease not found under Rosetta bin folder. \n  Insert the correct path wih the name of executable, or press enter to skip and add it later manually to the config.ini file: ")
					if not interface_analyzer:
						break
					else:
						continue 
			if os.path.isfile(interface_analyzer)== True:
				config.set('Paths', 'interface_analyzer', interface_analyzer)

		

		hbplus_path = raw_input("7) Insert path to executable of hbplus: ")
		hbplus_path = hbplus_path.lstrip()
		while os.path.isfile(hbplus_path) == False:
			hbplus_path = raw_input ("Insert the correct path to executable of hbplus, or press enter to skip and add it later manually to the config.ini file: ")
			hbplus_path = hbplus_path.lstrip()
			if not hbplus_path:
				break
			else:
				continue
		if   os.path.isfile(hbplus_path)  == True:
			config.set('Paths', 'hbplus_path', hbplus_path)

		reduce_path = raw_input("8) Insert path to executable of reduce: ")
		reduce_path = reduce_path.lstrip()
		while os.path.isfile(reduce_path) == False:
			reduce_path = raw_input ("Insert the correct path to executable of reduce, or press enter to skip and add it later manually to the config.ini file: ")
			reduce_path = reduce_path.lstrip()
			if not reduce_path:
				break
			else:
				continue
		if  os.path.isfile(reduce_path)  == True:
			config.set('Paths', 'reduce_path', reduce_path)




		reduce_db_path = raw_input("9) Insert path and name of HET dict. file for reduce: ")
		reduce_db_path = reduce_db_path.lstrip()
		path, filename = os.path.split(reduce_db_path)
		while os.path.isfile(reduce_db_path)== False or filename.startswith("reduce")==False  or filename.endswith(".txt") == False:
			reduce_db_path = raw_input ("Insert the correct path and name of reduce dictionary, or press enter to skip and add it later manually to the config.ini file: ")
		        path, filename = os.path.split(reduce_db_path)
			if not reduce_db_path:
				break
			else:
				continue
		if os.path.isfile(reduce_db_path)==True and filename.startswith("reduce")==True and filename.endswith(".txt") == True:
				config.set('Paths', 'reduce_db_path', reduce_db_path)

                configure_cluster = raw_input ("10) Do you want configure the parallelization? (Y/N): ")
                if configure_cluster.startswith("y") or configure_cluster.startswith("Y"):
                   cluster_name = raw_input ("\t-Cluster name: ")
                   if not cluster_name: cluster_name="None"
                   cluster_queue = raw_input ("\t-Name of queue to be used in the cluster: ")
                   if not cluster_queue: cluster_queue="None"
                   cluster_submit = raw_input ("\t-Submission method; usually this is 'submit' or 'sbatch' (default is 'submit'): ")
                   if not cluster_submit: cluster_submit="submit"
                   cluster_qstat = raw_input ("\t-Command to show queues; usually this is 'qstat' or 'squeue' (default is 'qstat'): ")
                   if not cluster_qstat: cluster_qstat="qstat"
                   command_queue = raw_input ('\t-File with header commands (i.e. shell and environmental variables)\n\t required for running scripts in the cluster: (default is "command_queues.txt" empty file)')
	           if not command_queue: command_queue=os.path.join(".","command_queues.txt")
                   config.set('Cluster','cluster_name',cluster_name)
                   config.set('Cluster','cluster_queue',cluster_queue)
                   config.set('Cluster','cluster_submit',cluster_submit)
                   config.set('Cluster','cluster_qstat',cluster_qstat)
                   config.set('Cluster','max_jobs_in_queue','750')
                   config.set('Cluster','min_jobs_in_queue','5')
                   config.set('Cluster','command_queue', command_queue)


		with open(abs_file_path, 'w') as configfile:    
			config.write(configfile)

                print "\nDone"



def find_file(path, pattern):
			p =pattern
                        found=False
			if os.path.exists(path):
				for file in os.listdir(path):
                                        if found: break
					if fnmatch.fnmatch(file, p):
						found=True
			else:
				return False
                        return found

def found_file(path, pattern):
			p =pattern
                        found=False
			if os.path.exists(path):
				for file in os.listdir(path):
                                        if found: break
					if fnmatch.fnmatch(file, p):
						found=True
                                                found_file=file
			else:
				return False
                        return found_file
	
#check if exe path inserted by user exists
def which(program):
	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return True
		else:
			for path in os.environ["PATH"].split(os.pathsep):
				exe_file = os.path.join(path, program)
				for elem in ext_elem(exe_file):
				  if is_exe(elem):
				    return True
	return False



def is_exe(fpath):
	return os.path.isfile(fpath) and os.access(fpath, os.X_OK) and os.path.isfile(fpath)	
def ext_elem(fpath):
 yield fpath
 for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
	yield fpath + ext 

if __name__ == '__main__':
    main()

