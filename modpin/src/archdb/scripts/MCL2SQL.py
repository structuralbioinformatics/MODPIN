import sys
import os
import argparse
import gzip
import re

from SBI.beans           import Path

from archdb.src import Subclass, Cclass, Loop
from archdb.src import start_transaction, end_transaction


# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'MCS2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database",
                        action = "store",       required = True,  help = "MCL output dir.")
    parser.add_argument('-q', '--looplist',     dest     = "looplist",
                        action = "store",       required = True,  help = "List of loops by order in PDB and length of first ss.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile",
                        action = "store",       required = True,  help = "Output SQL dir")
    parser.add_argument('-v', '--verbose',      dest     = "verbose",
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options


def MCL2SQL(database, looplist, sqlfile, verbose):

    Path.mkdir(sqlfile)
    grupfiles = {}
    for mclfile in Path.list_files(database):
        subclasstype = os.path.split(mclfile)[-1].split('.')[1]
        grupfiles.setdefault(subclasstype, []).append(mclfile)

    for subclasstype in grupfiles:
        classification = Cclass(subclasstype)
        loops = readlist(looplist, subclasstype)
        sql_fd = gzip.open(os.path.join(sqlfile, subclasstype + '.sql.gz'), 'wb')
        sql_fd.write(start_transaction())
        for mclfile in sorted(grupfiles[subclasstype]):
            subclassrange = os.path.split(mclfile)[-1].split('.')[2]
            if verbose:
                sys.stderr.write("Retrieving data for subclass {0} {1}...\n".format(subclasstype, subclassrange))
            sql_in = open(mclfile)
            read   = False
            for line in sql_in:
                dataline = line.rstrip('\n')
                #SKIP LINES
                if line.startswith('==') or line.startswith('***') or len(line.strip()) == 0 or line.startswith('---- P R O T E I N    C O D E  ----'):
                    continue
                if line.startswith('CONSENSUS & MULTIPLE ALIGNEMENT IN THE'):
                    data   = line.split(':')[-1].strip().split()
                    classification.subclasses = Subclass(tuple([data[0].strip(), data[3].strip()]), data[4], subclassrange)
                    workscls = classification.lastsubclass
                    read = True
                    continue
                if line.startswith('GLOBAL STATISTICS'):
                    read = False
                    continue
                if read:
                    if line.startswith('        SEQUENCE   ALIGNEMENT                           :'):   parse_mode, counter = 'P', 0
                    elif line.startswith('       ACCESSIBLE SURFACE ALIGNEMENT                    :'): parse_mode, counter = 'E', 0
                    elif line.startswith('           RAMACHANDRAN                                 :'): parse_mode, counter = 'R', 0
                    elif line.startswith('        SECONDARY STRUCTURE                             :'): parse_mode, counter = 'S', 0
                    elif line.startswith('--------- CONSENSUS THORNTON       :'):   workscls.add_consensus(dataline, 'MCL', loops)
                    elif line.startswith('--------- CONSENSUS TOPOLOGY'):           workscls.add_topology(dataline, 'MCL')
                    elif line.startswith('CENTROIDE POLAR COORD.   :'):             workscls.add_coordinates(dataline)
                    elif line.startswith('--------- RAMACHANDRAN PATTERN     :'):   workscls.ram_pat = re.sub('\(X\)','',dataline.split(':')[1].strip().strip('.'))
                    elif line.startswith('--------- SEQUENCE  PATTERN        :'):   workscls.seq_pat = re.sub('\(X\)','',dataline.split(':')[1].strip().strip('.'))
                    elif line.startswith('--------- BURIAL    PATTERN        :'):   workscls.exp_pat = re.sub('\(X\)','',dataline.split(':')[1].strip().strip('.'))

                    elif line.startswith('                             ') and len(dataline) < 400:
                        if parse_mode == 'P':
                            workscls.loops = Loop(info = dataline)
                        if parse_mode == 'E':
                            workscls.loops[counter].add_surface(info = dataline)
                            counter += 1
                        if parse_mode == 'R':
                            workscls.loops[counter].add_ramachandran(info = dataline)
                            counter += 1
                        if parse_mode == 'S':
                            workscls.loops[counter].add_secondary_str(info = dataline)
                            counter += 1

        sql_fd.write(classification.toSQL('MCL'))

        sql_in.close()
        sql_fd.write(end_transaction())
        sql_fd.close()

    if verbose:
        sys.stderr.write("End execution.\n")


def readlist(looplist, cstype):
    loops = {}
    fd = open(looplist)
    for l in fd:
        d = l.strip().split()

        if d[1] == cstype:
            n = d[0].split('_')
            key = tuple(['_'.join(n[0:2]), n[2]])
            loops[key] = int(d[2])
    fd.close()
    return loops

if __name__ == "__main__":
    options = options()
    MCL2SQL(options.database, options.looplist, options.sqlfile, options.verbose)
