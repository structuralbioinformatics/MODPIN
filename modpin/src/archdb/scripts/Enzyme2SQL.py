import sys, os, argparse, gzip

from SBI.beans           import Path
from SBI.databases import Enzymelink

from archdb.src import Source, Enzyme
from archdb.src import start_transaction, end_transaction


#
# deleted and transfered entries should have had their parents set by inheritance. they were not
# THIS MUST BE TAKEN INTO ACCOUNT WHEN REDOING THIS
#
#
# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'Enzyme2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database", 
                        action = "store",       required = True,  help = "Destination Directory for Enzyme database.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile", 
                        action = "store",       required = True,  help = "Output SQL file")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download", 
                        action = "store_true",  default  = False,  help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose", 
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def Enzyme2SQL(database, sqlfile, skip_download, verbose):
    enzyme_connect = Enzymelink(local = database)
    newsource      = None
    if not skip_download:
        if verbose: sys.stderr.write("Downloading Enzyme database to {0} ...\n".format(database))
        enzyme_connect.download()
        newsource = Source(name = 'enzyme', source = enzyme_connect.source)
        if verbose: sys.stderr.write("Download Finished.\n")
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    if verbose: sys.stderr.write("Parsing Enzyme.\n")
    if verbose: sys.stderr.write("Writing {0} ....\n".format(sqlfile))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(sqlfile, 'wb')
    sql_fd.write(start_transaction())
    if newsource is not None:
        sql_fd.write(newsource.toSQL())

    transfers = []
    for enz_line in enzyme_connect.localEnzymes:
        newenz = Enzyme(inline = enz_line)
        sql_fd.write(newenz.toSQL())
        if newenz.has_transfers:
            transfers.append(newenz.transfered2SQL())

    sql_fd.write("".join(transfers))
    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    Enzyme2SQL(options.database, options.sqlfile, options.skip_download, options.verbose)