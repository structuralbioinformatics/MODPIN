import sys, os, argparse, gzip

from SBI.beans           import Path
from SBI.databases import SCOPlink

from archdb.src import Source, SCOP
from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'SCOP2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database", 
                        action = "store",       required = True,  help = "Destination Directory for SCOP database.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile", 
                        action = "store",       required = True,  help = "Output SQL file")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download", 
                        action = "store_true",  default  = False,  help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose", 
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def SCOP2SQL(database, sqlfile, skip_download, verbose):
    scop_connect = SCOPlink(local = database)
    newsource      = None
    if not skip_download:
        if verbose: sys.stderr.write("Downloading SCOP database to {0} ...\n".format(database))
        scop_connect.download()
        newsource = Source(name = 'enzyme', source = scop_connect.source)
        if verbose: sys.stderr.write("Download Finished.\n")
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    if verbose: sys.stderr.write("Parsing SCOP.\n")
    if verbose: sys.stderr.write("Writing {0} ....\n".format(sqlfile))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(sqlfile, 'wb')
    sql_fd.write(start_transaction())
    if newsource is not None:
        sql_fd.write(newsource.toSQL())

    transfers = []
    scop_obj = SCOP()
    for line in scop_connect.descriptions:
        scop_obj.add_description(line.strip())
    for line in scop_connect.relations:
        scop_obj.add_relation(line.strip())

    sql_fd.write(SCOP.prepdbdeleted())
    sql_fd.write(scop_obj.toSQL())
    sql_fd.write(SCOP.afterpdbdeleted())
    
    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    SCOP2SQL(options.database, options.sqlfile, options.skip_download, options.verbose)
