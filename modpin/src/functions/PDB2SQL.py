import sys
import argparse
import os
import ConfigParser
import gzip
import shutil
import subprocess

# Add "." to sys.path #
src_path =  os.path.abspath(os.path.dirname(__file__))+"/../../scripts/"
sys.path.append(src_path)

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, 'config.ini')
config.read(config_file)

# Add functions path to sys.path
src_path = config.get('Paths', 'functions_path')
sys.path.append(src_path)

# Add special SBI library path to sys.path
src_path =  os.path.join(config.get('Paths', 'archdb_path'),"archdb")
sys.path.append(src_path)

# Add ArchDB path to sys.path
src_path = config.get('Paths', 'archdb_path')
sys.path.append(src_path)




# Import functions

from SBI               import SBIglobals
from SBI.beans           import Path
from SBI.databases import PDBlink

from archdb.src import Source, PDB
from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'PDB2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database",
                        action = "store",       required = True,  help = "Destination Directory for PDB database.")
    parser.add_argument('-q', '--seqdatabase',  dest     = "seqdatabase",
                        action = "store",       required = True,  help = "Destination Directory for PDBseq database.")
    parser.add_argument('-l', '--listfiles',    dest     = "listfiles",
                        action = "store",       default  = None,  help = "In case one needs to run only a subset.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile",
                        action = "store",       required = True,  help = "Output SQL dir")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download",
                        action = "store_true",  default  = False,  help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose",
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def PDB2SQL(database, seqdatabase, listfiles, sqlfile, skip_download, verbose):
    pdb_connect = PDBlink(local  = database, PDBseq = seqdatabase)
    newsource   = None

    if not skip_download:
        if verbose: sys.stderr.write("Syncronizing PDB database to {0} ...\n".format(database))
        pdb_connect.sync_PDB(log_file=os.path.join(database,'PDB.sync.log'))
        newsource = Source(name = 'PDB', source = pdb_connect.source)
        if verbose: sys.stderr.write("Creating PDBseq in {0} ...\n".format(seqdatabase))
        pdb_connect.make_PDBseq(log_file = os.path.join(seqdatabase,'PDB.seq.log'))
        if verbose: sys.stderr.write("Download Finished.\n")
        outdir = os.path.abspath(os.path.join(sqlfile,'00'))
        Path.mkdir(outdir)
        sql_fd = gzip.open(os.path.join(outdir, '0000.sql.gz'), 'wb')
        sql_fd.write(start_transaction())
        sql_fd.write(newsource.toSQL())
        sql_fd.write(end_transaction())
        sql_fd.close()
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    files2check = set()
    if listfiles is not None:
        fd = open(listfiles)
        for line in fd:
            files2check.add(line.strip())
        fd.close()
        logfd = open(listfiles + ".log","w")
    else:
        logfd = open("PDB2SQL.log","w")
    import traceback
    for pdbfile in pdb_connect.localPDBs:
        try:
            if listfiles is not None and pdbfile not in files2check:
                if len(files2check) == 0: break
                continue
            #else:
               # files2check.add(pdbfile)
               # files2check.remove(pdbfile)
            if verbose: sys.stderr.write("Working file {0}\n".format(pdbfile))
            newPDB = PDB(pdb_file = pdbfile)
            outsqldir  = os.path.join(sqlfile,newPDB.id[1:3].lower())
            Path.mkdir(outsqldir)
            outsqlfile = os.path.join(outsqldir, newPDB.id + '.sql.gz')
            # outsqlfile = os.path.join(os.getcwd(), newPDB.id + '.sql.gz')
            if verbose: sys.stderr.write("\tOutput SQL file is {0}.\n".format(outsqlfile))
            sql_fd = gzip.open(outsqlfile,'wb')
            sql_fd.write(start_transaction())
            sql_fd.write(PDB.preuniprotdeleted())
            sql_fd.write(newPDB.toSQL())
            sql_fd.write(PDB.afteruniprotdeleted())
            sql_fd.write(end_transaction())
            sql_fd.close()
        except KeyboardInterrupt:
            raise
        except:
            if verbose: sys.stderr.write("\tAn error occurred. Check log file\n")
            SBIglobals.alert('error', None, '\tAn error occurred for {0} . Check log file'.format(pdbfile))
            logfd.write("FILE {0}\n".format(pdbfile))
            logfd.write(traceback.format_exc())
            logfd.write("\n");

if __name__ == "__main__":
    SBIglobals.deepdebug = False
    options = options()
    if options.listfiles is not None:
        SBIglobals.output_file(options.listfiles + ".debug")
    PDB2SQL(options.database, options.seqdatabase, options.listfiles, options.sqlfile, options.skip_download, options.verbose)
