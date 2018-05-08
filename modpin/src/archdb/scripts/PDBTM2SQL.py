import sys, os, argparse, gzip

from SBI.beans           import Path
from SBI.databases import PDBTMlink

from archdb.src import Source, TM
from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'PDBTM2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database", 
                        action = "store",       required = True,  help = "Destination Directory for PDBTM database.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile", 
                        action = "store",       required = True,  help = "Output SQL file")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download", 
                        action = "store_true",  default  = False,  help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose", 
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def PDBTM2SQL(database, sqlfile, skip_download, verbose):
    pdbtm_connect  = PDBTMlink(local = database)
    newsource      = None
    if not skip_download:
        if verbose: sys.stderr.write("Downloading PDBTM database to {0} ...\n".format(database))
        #pdbtm_connect.download()
        newsource = Source(name = 'enzyme', source = pdbtm_connect.source)
        if verbose: sys.stderr.write("Download Finished.\n")
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    if verbose: sys.stderr.write("Parsing PDBTM.\n")
    if verbose: sys.stderr.write("Writing {0} ....\n".format(sqlfile))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(sqlfile, 'wb')
    sql_fd.write(start_transaction())
    sql_fd.write(TM.prepdbdeleted())
    if newsource is not None:
        sql_fd.write(newsource.toSQL())

    sql_fd.write(TM.regions2SQL())
    for line in pdbtm_connect.localTM:
        tmdata = TM(inline = line)
        sql_fd.write(tmdata.toSQL())

    sql_fd.write(TM.afterpdbdeleted())
    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    PDBTM2SQL(options.database, options.sqlfile, options.skip_download, options.verbose)
