import sys
import os
import argparse
import gzip

from SBI.beans           import Path
from SBI.structure.protein import Arch as SSpair

from archdb.src import Source, Arch
from archdb.src import start_transaction, end_transaction


# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'Arch2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database",
                        action = "store",       required = True,  help = "Arch output dir.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile",
                        action = "store",       required = True,  help = "Output SQL dir")
    parser.add_argument('-v', '--verbose',      dest     = "verbose",
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options


def Arch2SQL(database, sqlfile, verbose):

    if verbose:
        sys.stderr.write("Retrieving data from {0} ...\n".format(database))
    newsource = Source(name = 'PDBarch', source = "http://www-pdb.org/")
    outdir = os.path.join(os.path.join(os.path.abspath(sqlfile), '00'))
    Path.mkdir(outdir)
    sql_fd = gzip.open(os.path.join(outdir, '0000.sql.gz'), 'wb')
    sql_fd.write(start_transaction())
    sql_fd.write(newsource.toSQL())
    sql_fd.write(end_transaction())
    sql_fd.close()

    files_list_by_pdb = {}
    subdirs = ['archobj', 'superobj']
    for subdir in subdirs:
        for archobjfile in Path.list_files(os.path.join(database, subdir)):
            if archobjfile.endswith('.archObj'):
                data = tuple(os.path.splitext(os.path.split(archobjfile)[-1])[0].split('_')[2:])
                files_list_by_pdb[data] = archobjfile

    old_pdb = None
    newArchSet = None
    for dofdata in sorted(files_list_by_pdb):
        pdb = dofdata[0] + '_' + dofdata[1]
        if pdb != old_pdb:
            if old_pdb is not None:
                sql_fd.write(newArchSet.toSQL())
                sql_fd.write(end_transaction())
                sql_fd.close()
            outdir = os.path.join(os.path.join(os.path.abspath(sqlfile), dofdata[0][1:3].lower()))
            Path.mkdir(outdir)
            if verbose:
                sys.stderr.write("Retrieving loops from {0} ...\n".format(pdb))
            sql_fd = gzip.open(os.path.join(outdir, pdb + '.sql.gz'), 'wb')
            sql_fd.write(start_transaction())
            if verbose:
                sys.stderr.write("Printing data from {0} ...\n".format(pdb))
            old_pdb    = pdb
            newArchSet = Arch(pdb)
        newArchSet.archs = SSpair.load(files_list_by_pdb[dofdata])

    sql_fd.write(newArchSet.toSQL())
    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose:
        sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    Arch2SQL(options.database, options.sqlfile, options.verbose)
