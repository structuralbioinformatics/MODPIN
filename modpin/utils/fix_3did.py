import sys,os
import shutil

path="/home/boliva/PROJECTS/MODPIN/modppi/data/3did"
for file_brk in os.listdir(path):
  pdb=file_brk.strip()
  new=pdb.lower()
  print "mv %s %s"%(os.path.join(path,pdb),os.path.join(path,new))
  shutil.move(os.path.join(path,pdb),os.path.join(path,new))

