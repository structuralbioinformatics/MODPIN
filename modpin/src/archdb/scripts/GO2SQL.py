import sys, os, argparse, gzip

from SBI.beans           import Path
from SBI.databases import GOlink

from archdb.src import Source, GOterm
from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'GO2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database", 
                        action = "store",       required = True,  help = "Destination Directory for GO database.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile", 
                        action = "store",       required = True,  help = "Output SQL file")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download", 
                        action = "store_true",  default  = False,  help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose", 
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def GO2SQL(database, sqlfile, skip_download, verbose):
    go_connect = GOlink(local = database)
    newsource  = None
    if not skip_download:
        if verbose: sys.stderr.write("Downloading GO database to {0} ...\n".format(database))
        go_connect.download()
        newsource = Source(name   = 'GO', source = go_connect.source)
        if verbose: sys.stderr.write("Download Finished.\n")
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    with_parents   = []
    with_relations = []

    if verbose: sys.stderr.write("Parsing GO.\n")
    if verbose: sys.stderr.write("Writing {0} ....\n".format(sqlfile))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(sqlfile, 'wb')
    sql_fd.write(start_transaction())
    if newsource is not None:
        sql_fd.write(newsource.toSQL())

    for go_line in go_connect.localGOs:
        newGO = GOterm(inline = go_line)
        sql_fd.write(newGO.toSQL() + "\n")
        if len(newGO.relations) > 0:
            with_relations.append(newGO)
        if len(newGO.parents) > 0:
            with_parents.append(newGO)

    for GO in with_relations:
        sql_fd.write(GO.relations2SQL() + "\n")
    for GO in with_parents:
        sql_fd.write(GO.parents2SQL() + "\n")

    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    GO2SQL(options.database, options.sqlfile, options.skip_download, options.verbose)