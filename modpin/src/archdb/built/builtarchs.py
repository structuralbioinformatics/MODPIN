import sys, os, argparse

from SBI            import SBIglobals
from SBI.beans      import File, Path
from SBI.databases  import PDBlink
from SBI.utilities  import archer

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'builtArchs',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-l', '--list',         dest     = "list", 
                        action = "store",       required = True,  help = "List of selected PDB_chain.")
    parser.add_argument('-p', '--pdbdir',       dest     = "pdbdir", 
                        action = "store",       required = True,  help = "Dir of local PDB database.")
    parser.add_argument('-o', '--output',       dest     = "output", 
                        action = "store",       required = True,  help = "Dir to output the archdata.")
    parser.add_argument('-s', '--sort',         dest     = "sort", 
                        action = "store_true",  default  = False, help = "Skip building loops, only sort them.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose", 
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def PDBlist2dict(inputfile, separator = '_'):
    data = []
    inputf = File(inputfile,'r')
    for line in inputf.descriptor:
        data.append(line.strip().split(separator))
    return data

def builtArchDB(pdblist, pdbdir, outdir):
    pdb_connect = PDBlink(local = pdbdir)
    for pdbinfo in pdblist:
        pdb, chain  = pdbinfo
        subdir      = pdb[1:3].lower()
        pdbfile     = pdb_connect.get_PDB(pdb)
        SBIglobals.alert('verbose', None, 'Processing file: {0}'.format(pdbfile))
        archs       = archer.build_archs(sourcepdb = pdbfile, chain = chain, limit_distance = 25)
        for archkey in archs:
            if len(archs[archkey]) > 0: 
                Path.mkdir(os.path.join(outdir[archkey],subdir))
                Path.mkdir(os.path.join(outdir['STRUC'],subdir))
            for arch in archs[archkey]:
                pyobjName = "_".join([str(arch.aminoacid_distance), arch.type, arch.identifier]) + '.archObj'
                arch.dump(os.path.join(os.path.join(outdir[archkey],subdir),pyobjName))
                arch.format2file(filename  = os.path.join(os.path.join(outdir['STRUC'],subdir),arch.identifier), 
                                 extension = 'pdb', center = True)
                arch.format2file(filename  = os.path.join(os.path.join(outdir['STRUC'],subdir),arch.identifier), 
                                 extension = 'js',  center = True)

if __name__ == "__main__":
    options            = options()
    SBIglobals.verbose = options.verbose
    pdblist            = PDBlist2dict(options.list)
    outdir = {'ARCHS' : os.path.join(options.output,'archobj'), 
              'SUPER' : os.path.join(options.output,'superobj'),
              'STRUC' : os.path.join(options.output,'structures'),
              'DBARC' : os.path.join(options.output,'db'),
              'DBSUP' : os.path.join(options.output,'superdb')}
    if not options.sort:
        builtArchDB(pdblist, options.pdbdir, outdir)
    archer.sortarchs(outdir['ARCHS'],  outdir['DBARC'])
    archer.sortarchs(outdir['SUPER'],  outdir['DBSUP'])
