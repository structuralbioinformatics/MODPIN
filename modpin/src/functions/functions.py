import os, sys, re
import collections
import gzip
import itertools
import json
import numpy
import subprocess
import pwd
import time

# Add "." to sys.path #
sys.path.append(os.path.dirname(__file__))

# Imports jbonet's module #
from SBI.structure       import PDB
from SBI.structure.chain import ChainOfNucleotide
#from SBI.data            import aminoacids3to1, aminoacids_polarity_boolean, aminoacids_surface, nitrogenous_bases, dna_complementary
from SBI.data            import aminoacids3to1, aminoacids_polarity_boolean, aminoacids_surface
from SBI.data            import *

#-------------#
# Classes     #
#-------------#

class DSSP(object):
    '''
    This class defines a DSSP object.
    '''

    def __init__(self, dssp_file):
        self._file = dssp_file
        self._residues = None

    def _check_parsing(self):
        if self._residues == None:
            self._residues = {}
            self._parse_file()

    def _get_file(self):
        return self._file

    def _parse_file(self):
        if os.path.exists(self._get_file()):
            fd = open(self._get_file(), "r")
            for line in fd:
                # Get DSSP info #
                m = re.search("(\d+)\s\S\s[ACDEFGHIKLMNPQRSTVWY]", line)
                if m:
                    chain = line[11:12]
                    residue_num = int(line[5:10])
                    accessible_surface_area = float(line[35:38])
                    secondary_structure = line[16:17]
                    self._residues[(chain, residue_num)] = (accessible_surface_area, secondary_structure)
            fd.close()
        else:
            raise ValueError("Could not open DSSP file %s" % self._get_file())

    def has_residue(self, chain, residue_num):
        self._check_parsing()

        if (chain, residue_num) in self._residues:
            return True

        return False

    def get_accessible_surface_area(self, chain, residue_num):
        self._check_parsing()

        if self.has_residue(chain, residue_num):
            return self._residues[(chain, residue_num)][0]

        return None

    def get_secondary_structure(self, chain, residue_num):
        self._check_parsing()

        secondary_structures = set("EH")

        if self._residues[(chain, residue_num)][1] not in secondary_structures:
            return "C" # residue is in a "coil" region

        return self._residues[(chain, residue_num)][1]

class PirAlignmentMod(object):
    '''
    This class defines a PIR alignment file.
    '''

    def __init__(self, alignment_file):
        self._file = alignment_file
        self._query = None
        self._pdb_name = None
        #self._type = None
        #self._pdb_chains = None
        # self._potential = None
        self._query_alignments = None
        self._hit_alignments = None

        self._parse_file()

    def _get_file(self):
        return self._file

    def _parse_file(self):
        if os.path.exists(self._get_file()):
            self._query_alignments = []
            self._hit_alignments = []
            alignment = None

            fd = open(self._get_file(), "rt")
            for line in fd.readlines():
                # m = re.search("^structure:(.{4}):\d+:.:\d+:.:(dimer|monomer):(\S+):(family|general):", line)
                m = re.search("^structureX:(.{4}):\d+:(\w):.:(\w):", line)
                if m:
                    self._pdb_name = m.group(1)
                    self._pdb_chains = m.group(2), m.group(3)
                    alignment = "structure"
                m = re.search("^sequence:(\S+):", line)
                if m:
                    self._query = m.group(1)
                    alignment = "sequence"
                m = re.search("^(.+)[/*]$", line)
                if m:
                    if alignment == "structure":
                        self._hit_alignments.append(m.group(1))
                    if alignment == "sequence":
                        self._query_alignments.append(m.group(1))
            fd.close()
        else:
            raise ValueError("Could not open statistical potentials file %s" % self._get_file())

    def get_query(self):
        return self._query

    def get_pdb_name(self):
        return self._pdb_name

    def get_pdb_chains(self):
        #return sorted(self._pdb_chains)
        return self._pdb_chains

    def get_type(self):
        return self._type

    # def get_potential(self):
    #     return self._potential

    def get_alignments(self, atype="query"):
        if atype != "query" and atype != "hit":
            atype = "query"

        if atype == "query":
            return self._query_alignments

        return self._hit_alignments

#-------------#
# Accessory   #
#-------------#

def make_subdirs(main, subdirs):
    '''
    This function makes all subdirs listed in "subdirs".
    '''

    for subdir in subdirs:
        if not os.path.exists(os.path.join(main, subdir)):
            os.makedirs(os.path.join(main, subdir))

def remove_files(files):
    '''
    This function removes all files listed in "files".
    '''

    for each_file in files:
        if os.path.exists(each_file):
            os.remove(each_file)

#-------------#
# Parsers     #
#-------------#

def parse_list_file(list_file, gz=False):
    '''
    This function parses any "list" file and returns a generator.
    '''
    #print("Parse list %s\n"%list_file)
    if os.path.exists(list_file):
        fd = None

        if gz: fd = gzip.open(list_file, "rt")
        else:  fd = open(list_file, "rt")
        for line in fd:
            yield line.strip("\n")
        fd.close()
    else:
        raise ValueError("Could not open list file %s" % list_file)

def parse_csv_file(csv_file, gz=False):
    '''
    This function parses any "CSV" file and returns a generator of lists.
    '''

    if os.path.exists(csv_file):
        fd = None

        if gz: fd = gzip.open(csv_file, "rt")
        else:  fd = open(csv_file, "rt")
        for line in fd:
            yield line.strip("\n").split(",")
        fd.close()
    else:
        raise ValueError("Could not open CSV file %s" % csv_file)

def parse_fasta_file(fasta_file, gz=False, clean=True):
    '''
    This function parses any "FASTA" file and iteratively yields sequences as a tuple
    of the form (identifier, sequence).
    '''

    if os.path.exists(fasta_file):
        identifier = ""
        sequence = ""

        for line in parse_list_file(fasta_file, gz):
            line = line.strip()
            if line[0] == ">":
                if sequence != "":
                    yield (identifier, sequence)
                m = re.search("^>(\S+)\|*", line)
                identifier = m.group(1)
                sequence = ""
            else:
                if clean:
                    sequence += re.sub("\W|\d", "X", line.upper())
                else:
                    sequence += line.upper()
        yield (identifier, sequence)
    else:
        raise ValueError("Could not open FASTA file %s" % fasta_file)

#------------#
# SBI module #
#------------#

def get_aminoacid_cb_or_ca(aminoacid):
    """
    Returns the CB {Atom} or the CA {Atom}; otherwise return "None"
    @rtype: {Atom}
    """

    if aminoacid.has_cb:
        return aminoacid.cb
    elif aminoacid.has_ca: # for glycines
        return aminoacid.ca

    return None


#-------------#
# PDB         #
#-------------#

def get_residue_residue_contacts(chain_A, chain_B):
    '''
    This function gets all residue-residue contacts according to the
    definition my Mosca R., Ceol A. & Aloy P., 2013.
    '''
    disulfide_bridge_threshold = 2.56
    hydrogen_bond_threshold = 3.5
    salt_bridge_threshold = 5.5
    van_der_waals_threshold = 5.0
    contacts = {}

    # For each amino acid... #
    for aminoacid_A in chain_A.aminoacids:
        # For each amino acid... #
        for aminoacid_B in chain_B.aminoacids:
            contact = aminoacid_A.distance(aminoacid_B, "min")
            # Skip if amino acids are too apart or distance could not be calcualted #
            if contact[-1] == -1 or contact[-1] > 5.5: continue
            # For each atom... #
            for atom_A in aminoacid_A.atoms:
                atom_type_A = get_atom_type(atom_A.pretty_name)
                # Skip if atom type cannot make contacts #
                if atom_type_A == None: continue
                # For each atom... #
                for atom_B in aminoacid_B.atoms:
                    atom_type_B = get_atom_type(atom_B.pretty_name)
                    # Skip if atom type cannot make contacts #
                    if atom_type_B == None: continue
                    # Get distance #
                    distance = atom_A.distance(atom_B)
                    # Skip if atoms are too apart or distance could not be calcualted #
                    if distance == -1 or distance > 5.5: continue
                    # Get interactions #
                    disulfide_bridge = is_disulfide_bridge(aminoacid_A, atom_type_A, aminoacid_B, atom_type_B, distance)
                    hydrogen_bond = is_hydrogen_bond(aminoacid_A, atom_type_A, aminoacid_B, atom_type_B, distance)
                    salt_bridge = is_salt_bridge(aminoacid_A, atom_type_A, aminoacid_B, atom_type_B, distance)
                    van_der_waals = is_van_der_waals(aminoacid_A, atom_type_A, aminoacid_B, atom_type_B, distance)
                    if disulfide_bridge or hydrogen_bond or salt_bridge or van_der_waals:
                        if atoms_clash(atom_type_A, atom_type_B, distance, disulfide_bridge) == False:
                            contacts.setdefault(((chain_A.chain, aminoacid_A.number), (chain_B.chain, aminoacid_B.number)), [])
                            contacts[(chain_A.chain, aminoacid_A.number), (chain_B.chain, aminoacid_B.number)].append([(aminoacid_A.standard_type, atom_type_A), (aminoacid_B.standard_type, atom_type_B), distance, disulfide_bridge, hydrogen_bond, salt_bridge, van_der_waals])

    return contacts

def get_atom_type(atom_name):
    '''
    This function returns wheter an atom is a "C", "N", "O", "S", or "None" of these.
    '''

    m = re.search("^ C  |\WC\WA|\WC\WB|\WC\WC|\WC\WD|\WC\WE|\WC\WF|\WC\WG|\WC\WH|\WC\WI|\WC\WJ|\WC\WK|\WC\WL|\WC\WM|\WC\WN|\WC\WO|\WC\WP|\WC\WQ|\WC\WR|\WC\WS|\WC\WT|\WC\WU|\WC\WV|\WC\WW|\WC\WX|\WC\WY|\WCA\W|\WCAA|\WCAB|\WCAC|\WCAE|\WCAF|\WCAG|\WCAH|\WCAI|\WCAJ|\WCAK|\WCAL|\WCAM|\WCAN|\WCAO|\WCAP|\WCAQ|\WCAR|\WCAS|\WCAT|\WCAU|\WCAX|\WCAY|\WCAZ|\WCB\W|\WCBA|\WCBB|\WCBC|\WCBD|\WCBE|\WCBF|\WCBG|\WCBH|\WCBI|\WCBX|\WCC\W|\WCCA|\WCCB|\WCCX|\WCD\W|\WCDA|\WCDX|\WCE\W|\WCEA|\WCEB|\WCEC|\WCEX|\WCF\W|\WCFX|\WCG\W|\WCGC|\WCGX|\WCH\W|\WCHA|\WCHB|\WCHC|\WCHD|\WCHX|\WCI\W|\WCIX|\WCJ\W|\WCJX|\WCK\W|\WCKX|\WCL\W|\WCM\W|\WCMA|\WCMB|\WCMC|\WCMD|\WCME|\WCMP|\WCMT|\WCMX|\WCMZ|\WCN\W|\WCNT|\WCNX|\WCO\W|\WCOX|\WCP\W|\WCPX|\WCQ\W|\WCQX|\WCS\W|\WCT\W|\WCW\W|\WCX\W|\WCXD|\WCXE|\WCXF|\WCXG|\WCXN|\WCXO|\WCXP|\WCXQ|\WCXT|\WCXU|\WCXV|\WCYT|\WCZ\W$", atom_name)
    if m:
        return "C"
    m = re.search("^ N  |\WN\WA|\WN\WB|\WN\WC|\WN\WD|\WN\WE|\WN\WF|\WN\WM|\WN\WP|\WN\WS|\WN\WT|\WN\WX|\WN\WY|\WNA\W|\WNAA|\WNAB|\WNAH|\WNAI|\WNAV|\WNAW|\WNAX|\WNB\W|\WNBI|\WNBJ|\WNBK|\WNC\W|\WND\W|\WNE\W|\WNF\W|\WNG\W|\WNH\W|\WNI\W|\WNJ\W|\WNK\W|\WNL\W|\WNLX|\WNM\W|\WNN\W|\WNO\W|\WNP\W|\WNPA|\WNPB|\WNPC|\WNPD|\WNQ\W|\WNR\W|\WNRO|\WNS\W|\WNT\W|\WNW\W|\WNX\W|\WNXT|\WNXU|\WNXV|\WNXW|\WNZ\W$", atom_name)
    if m:
        return "N"
    m = re.search("^ O  |\WO\WA|\WO\WB|\WO\WC|\WO\WD|\WO\WE|\WO\WF|\WO\WG|\WO\WL|\WO\WM|\WO\WP|\WO\WQ|\WO\WR|\WO\WS|\WO\WT|\WO\WV|\WO\WX|\WO\WY|\WOA\W|\WOAB|\WOAD|\WOAP|\WOAX|\WOAY|\WOB\W|\WOBC|\WOC\W|\WOCC|\WOCD|\WOD\W|\WODA|\WODB|\WOE\W|\WOEA|\WOF\W|\WOG\W|\WOGL|\WOH\W|\WOHB|\WOHN|\WOI\W|\WOJ\W|\WOL\W|\WOM\W|\WON\W|\WOO\W|\WOP\W|\WOPP|\WOR\W|\WORA|\WOS\W|\WOT\W|\WOW\W|\WOX\W|\WOXA|\WOXS|\WOXT|\WOXX$", atom_name)
    if m:
        return "O"
    m = re.search("^ S  |\WS\WA|\WS\WG|\WS\WP|\WSBC|\WSC\W|\WSD\W|\WSG\W|\WSP\W$", atom_name)
    if m:
        return "S"

    return None

def is_disulfide_bridge(aminoacid_A, atom_type_A, aminoacid_B, atom_type_B, distance):
    '''
    This function returns whether two atoms form a disulfide bridge or not.
    
    According to definition by Mosca R., Ceol A. & Aloy P., 2013, a disulfide
    bridge is formed by any atom pair S-S from 2 Cys at 2.56A.
    '''

    return aminoacid_A.single_letter == "C" and aminoacid_B.single_letter == "C" and atom_type_A == "S" and atom_type_B == "S" and distance <= 2.56

def is_hydrogen_bond(aminoacid_A, atom_type_A, aminoacid_B, atom_type_B, distance):
    '''
    This function returns whether two atoms form an hydrogen bond or not.
    
    According to definition by Mosca R., Ceol A. & Aloy P., 2013, an hydrogen
    bond is formed by any atom pair N-O and O-N at 3.5A.
    '''


    return ((atom_type_A == "N" and atom_type_B == "O") or (atom_type_A == "O" and atom_type_B == "N")) and distance <= 3.5

def is_salt_bridge(aminoacid_A, atom_type_A, aminoacid_B, atom_type_B, distance):
    '''
    This function returns whether two atoms form a salt bridge or not.
    
    According to definition by Mosca R., Ceol A. & Aloy P., 2013, a salt
    bridge is formed by any atom pair N-O and O-N at 5.5A.
    '''
    negative_residues = ['D', 'E']
    positive_residues = ['R', 'H', 'K', 'S', 'Y']

    if aminoacid_A.single_letter in negative_residues and aminoacid_B.single_letter in positive_residues:
        return atom_type_A == "O" and atom_type_B == "N" and distance <= 5.5

    if aminoacid_A.single_letter in positive_residues and aminoacid_B.single_letter in negative_residues:
        return atom_type_A == "N" and atom_type_B == "O" and distance <= 5.5

    return False

def is_van_der_waals(aminoacid_A, atom_type_A, aminoacid_B, atom_type_B, distance):
    '''
    This function returns whether two atoms form a van der waals interaction
    or not.
    
    According to definition by Mosca R., Ceol A. & Aloy P., 2013, a van der
    waals interaction is formed by any atom pair C-C at 5.0A.
    '''

    return atom_type_A == "C" and atom_type_B == "C" and distance <= 5.0

def atoms_clash(atom_type_A, atom_type_B, distance, disulfide_bridge):
    '''
    This function returns wheter two atoms clash or not.

    According to definition by Mosca R., Ceol A. & Aloy P., 2013, any atom
    pairs at distance less than the sum of the two covalent radii plus 0.5A
    that are not forming a disulfide bridge are considered clashes.
    '''
    covalent_radii = {'C': 0.75, 'N': 0.71, 'O': 0.63}

    if disulfide_bridge == False:
        return distance < (covalent_radii[atom_type_A] + covalent_radii[atom_type_B] + 0.5)

    return False

#-------------#
# Modelling   #
#-------------#

def smith_waterman(A, B, emboss_path):
    '''
    This function aligns a pair of sequences "A" and "B" using water
    from the EMBOSS package.

    @return           = tuple(query alignment, hit alignment, query_start, query_end, hit_start, hit_end)

    type query align. = str
    type hit align.   = str
    type query start  = int
    type query end    = int
    type hit start    = int
    type hit end      = int
    '''
    alignment = []
    positions = []
    

    try:
        process = subprocess.Popen([os.path.join(emboss_path, "water"), "-asequence", A, "-bsequence", B, "-gapopen", "10.0", "-gapextend", "0.5", "-outfile", "stdout"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        for line in process.stdout:
            # Skip commented lines #
            if line.startswith("#"): continue
            # Capture alignment #
            m = re.search("^\S+\s+(\d+)\s+([\w-]+)\s+(\d+)$", line)
            if m:
                alignment.append(m.group(2))
                positions.append([int(m.group(1)), int(m.group(3))])

        query_positions = positions[0::2]
        hit_positions = positions[1::2]

        return("".join(alignment[0::2]), "".join(alignment[1::2]), query_positions[0][0], query_positions[-1][-1], hit_positions[0][0], hit_positions[-1][-1])

    except:
        return(None, None, None, None, None, None)

def get_sequence_to_crystal_correlations(pdb, chain, gapped=True):
    '''
    This function gets the correlation 1 to 1 between the sequence
    and the crystal positions, and returns them.
    
    @return           = correlation, correlation

    type correlation  = dict
    '''
    sequence_to_crystal = {}
    crystal_to_sequence = {}

    # For each protein chain... #
    for protein_chain in pdb.proteins:
        # Skip if not chain of interest #
        if protein_chain.chain == chain:
            # Get protein sequence #
            protein_sequence = protein_chain.gapped_protein_sequence
            if gapped == False:
                protein_sequence = protein_chain.protein_sequence
            # Get protein sequence positions #
            sequence_positions = range(len(protein_sequence))
            # Get crystal positions #
            crystal_positions = [aminoacid.number for aminoacid in protein_chain.aminoacids]
            # For each nucleotide... #
            for aminoacid in protein_sequence:
                if aminoacid == "x":
                    sequence_positions.pop(0)
                    continue
                sequence_to_crystal[sequence_positions[0]] = crystal_positions[0]
                crystal_to_sequence[crystal_positions[0]] = sequence_positions[0]
                # Remove last #
                sequence_positions.pop(0)
                crystal_positions.pop(0)

    return sequence_to_crystal, crystal_to_sequence

def filter_hit(query_alignment, hit_alignment, pdb, contacts, dssp, gaps=True):
    '''
    This function filters a hit if the query does not cover all of the core
    (i.e. region of the hit that contacts the DNA). Additionally, if gaps is
    "False", it filters a hit if either the query or the hit alignments are
    gapped. If gaps is set to "True", it allows insertions/deletions in un-
    structured regions (i.e. "C" secondary structure regions) or up to 1 gap
    in structured regions.

    @return           = boolean
    '''

    # Get sequence/PDB correlation #
    sequence_to_crystal, crystal_to_sequence = get_sequence_to_crystal_correlations(pdb=pdb, chain=pdb.chains[0].chain)
    # Get sequence #
    sequence = pdb.chains[0].gapped_protein_sequence
    # Get core residues #
    core_residues = contacts.get_core_residues()
    # Get secondary structure #
    secondary_structure = ""
    for position in range(len(sequence)):
        if position + 1 in core_residues:
            secondary_structure += "*"
        elif sequence[position] == "x":
            secondary_structure += "C"
        else:
            secondary_structure += dssp.get_secondary_structure(pdb.chains[0].chain, sequence_to_crystal[position])
    # Get start, end hit #
    start = sequence.index(hit_alignment.replace("-", ""))
    end = start + len(hit_alignment.replace("-", "")) - 1
    # Get query secondary structure #
    query_secondary_structure = ""
    positions = range(start, end + 1)
    for i in range(len(hit_alignment)):
        if hit_alignment[i] != "-":
            position = positions.pop(0)
            if query_alignment[i] != "-":
                query_secondary_structure += secondary_structure[position]
            else:
                query_secondary_structure += "-"
        else:
            query_secondary_structure += "-"
    # If contact residues are not conserved, skip #
    if len(re.findall("\*", query_secondary_structure)) != len(re.findall("\*", secondary_structure)):
        return True
    # If insertions/deletions in structured region holding contact residues, skip #
    if re.search("[\*|E]+\-+[\*|E]+", query_secondary_structure) or re.search("[\*|H]+\-+[\*|H]+", query_secondary_structure):
        return True

    return False


def realign(a,b,template):
    """Realign the 'b' gapped sequence which was aligned with 'a' upon a new 'template' of 'a' """
    c=''
    k=0
    j=0
    #print "Original sequence A %s length %d"%(a,len(a))
    #print "Original sequence B %s length %d"%(b,len(b))
    #print "New template for  A %s length %d"%(template,len(template))
    while (k<len(a) or j<len(template)):
        if k<len(a) and j<len(template):
            if (template[j] == a[k] and not template[j] == '-'):
                c+=b[k]
                #print"(%d,%d) %s %s %s"%(j,k,c,a[0:k],template[0:k]) 
                k=k+1
                j=j+1
            elif (template[j] == '-'):
                c+='-'
                #print"(%d,%d) %s %s %s"%(j,k,c,a[0:k],template[0:k]) 
                j=j+1
            elif (a[k] == '-'):
                #print"(%d,%d) %s %s %s"%(j,k,c,a[0:k],template[0:k]) 
                if (k+1<len(a)):
                    print("Skip original sequence known gap %d (.. %s ..)\n"%(k,a[k]))
                else:
                    print("Skip known gap %d (.. %s ..)\n"%(k,a[k]))
                k=k+1
            else:
                #print("Some mistake in (%d,%d) A %s Template %s \n"%(j,k,a[0:k],template[j]))
                c+='-'
                #print"(%d,%d) %s %s %s"%(j,k,c,a[0:k],template[0:k]) 
                k=k+1
                j=j+1
        elif j<len(template):
            c+='-'
            j=j+1 
        elif k<len(a):
            c+=b[k]
            k=k+1
    #print "New alignement of B %s"%c
    return c

def fileExist(file):
    '''
    Check existing files
    '''
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False

def printverbose(f,flag,message):
    '''
    Define verbose to print in file 'f', if 'flag'=True, message given as 'message'
    '''
    if flag: f.write("%s"%(message))

def printfasta (out,name,seq):
    '''
    Print Fasta format
    '''
    out.write(">%s\n"%(name))
    n=0;
    w=""
    for s in seq:
        w += s
        n+=1
        if n == 80:
            out.write("%s\n"%(w))
            n=0
            w=""
    if n>0: out.write("%s\n"%(w))

def PPI_iterator(filename):
    '''
    Iterate on PPIs raw data with alternate names
    '''
    fd=open(filename,"r");
    for line in fd:
        if not line.startswith("#"):
            word=line.strip().split()
            if len(word)>4:
                a=word[0]
                b=word[1]
                alt_a=word[2]
                alt_b=word[3]
                yield (a,b,alt_a,alt_b)
            elif len(word)>1:
                a=word[0]
                b=word[1]
                yield (a,b,"","")
            else:
                continue
    fd.close()

def get_nodes(edges):
    '''
    Gets the list of nodes from a list of edges
    '''
    n=set()
    if len(edges)>0:
        for ppi in edges:
            a,b=ppi
            n.add(a)
            n.add(b)
        if len(n)>0: 
            return n
        else:        
            return None
    else:
        return None



#-------------#
# Cluster     #
#-------------#
def submit_command_to_queue(command, queue=None, max_jobs_in_queue=None, queue_file=None, dummy_dir="/tmp", submit="qsub", qstat="qstat", job_label=""):
    """
    This function submits any {command} to a cluster {queue}.

    @input:
    command {string}
    queue {string} by default it submits to any queue
    max_jobs_in_queue {int} limits the number of jobs in queue
    queue_file is a file with information specific of the cluster for running a queue

    """
    import hashlib

    if max_jobs_in_queue is not None:
        while number_of_jobs_in_queue(qstat) >= max_jobs_in_queue: time.sleep(5)

    cwd = os.path.join(dummy_dir,"sh")
    if not os.path.exists(cwd): os.makedirs(cwd)
    script= os.path.join(cwd,"submit_" + job_label + hashlib.sha224(command).hexdigest() + ".sh")
    if queue_file is not None:
      fd=open(script,"w")
      with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
      fd.close()
      queue_standard.close()
      if queue is not None:
       if  submit=="qsub":
          os.system("%s -q %s %s" % (submit, queue,script))
       elif submit=="sbatch":
          os.system("%s -p %s %s" % (submit, queue,script))
       else:
          os.system("%s %s"% (submit,script))
      else:
        os.system("%s %s"% (submit,script))
    else:
      if queue is not None:
        os.system("echo \"%s\" | %s -q %s" % (command, submit, queue))
      else:
        os.system("echo \"%s\" | %s" % (submit,command))


def number_of_jobs_in_queue(qstat="qstat"):
    """
    This functions returns the number of jobs in queue for a given
    user.

    """

    # Initialize #
    user_name = get_username()
    process = subprocess.check_output([qstat, "-u", user_name])

    return len([line for line in process.split("\n") if user_name in line])

def number_of_relevant_jobs_in_queue(qstat="qstat", job_label=""):
    """
    Returns the number of jobs in the queue that have the current user and the
    current job label.
    """
    user_name = get_username()

    if qstat == "qstat":
        process = subprocess.check_output("{} -f | grep -C 1 {}@".format(qstat, user_name), shell=True)
    else:
        process = subprocess.check_output('{} -u {} --format="%u %j"'.format(qstat, user_name), shell=True)

    job_re = re.compile(r'\ssubmit_{}'.format(re.escape(job_label)))
    return len([line for line in process.split("\n") if user_name in line and job_re.search(line)])

def get_username():
    """
    This functions returns the user name.

    """

    return pwd.getpwuid(os.getuid())[0]

def wait_until_jobs_finish(qstat="qstat", job_label=""):
    """
    Waits until jobs with name formatted as 'submit_{job_label}*' are gone from
    the queue, then returns True. If the process is interrupted, returns False.
    """
    last_num = 1e5
    while last_num > 0:
        try:
            new_num = number_of_relevant_jobs_in_queue(qstat=qstat, job_label=job_label)
            if new_num > 0:
                if new_num != last_num:
                    print("Waiting for {} jobs to complete...".format(new_num))
                time.sleep(5)
            last_num = new_num
        except KeyboardInterrupt:
            return False
    return True

def add_hydrogens(config,path,inp,out):
    #Initialize
    from SBI.structure import PDB
    import shutil
    src_path = config.get('Paths','modppi_path')
    hbplus = config.get('Paths', 'hbplus_path')
    reduce_exe  =  config.get('Paths', 'reduce_path')
    reduce_db  =  config.get('Paths', 'reduce_db_path')
    relax_exe = config.get('Paths', 'relax_exe')
    hydrogen_type = config.get('Parameters','hydrogens')
    relax = config.get('Parameters','relax')
    cwd = os.getcwd()
    os.chdir(path)
    if fileExist(inp):
     if len(inp.split('.')) > 0: output_hbplus = ".".join(inp.split('.')[:-1]) + ".h"
     else: output_hbplus= inp.strip()+".h"
     if hydrogen_type == "full":
       os.system("%s -Quiet %s -DB %s> %s"%(reduce_exe,inp,reduce_db,output_hbplus))
     else:
       os.system("%s -o %s >& hbplus.log"%(hbplus,inp))
     if relax == "yes":
       sys.stdout.write("\t\t\t-- Relaxing the hydrogen-intermediate model %s (see Rosetta output in relax.log and score.sc)...\n"%output_hbplus)
       os.system("%s -s %s -in:file:fullatom -nstruct 1  -packing:repack_only >& relax.log"%(relax_exe,output_hbplus))
       opt_model=".".join(output_hbplus.split('.')[:-1])+"_0001.pdb"
       old_model=".".join(output_hbplus.split('.')[:-1])+"_non_optimized.pdb"
       shutil.move(output_hbplus,old_model)
       if fileExist(opt_model):
          check_pdb=PDB(opt_model)
          if check_pdb.has_protein:
            check_pdb.clean()
            check_pdb.write(output_hbplus)
            try:
              sys.path.remove(opt_model)
            except:
              sys.stdout.write("\t\t\t-- Keeping old file %s ...\n"%opt_model)
          else:
            shutil.copy(old_model,output_hbplus)
       else:
          shutil.copy(old_model,output_hbplus)
     if not fileExist(output_hbplus):
       raise ValueError("Cannot find file with hydrogen atoms")
     else:
      pdb=PDB(output_hbplus)
      pdb.clean()
      pdb.write(out,force=True)
    os.chdir(cwd)

def renumber_pdb(config,path,pdb_name,sequences,dummy_dir):
    ''' 
     Renumber PDB file located in path folder with the real sequences

     path	Folder where PDB file is located
     pdb 	PDB file
     sequences  dictionary of sequences (of ProteinSequence Class from SeqIO) that define the Aa number
                chain identifier is the key of the dictionary
     dummy_dir  Dummy directory to cerate files

    '''

    #Initialize
    from SBI.structure.chain import Chain
    from SBI.sequence import Sequence
    from SBI.structure import PDB
    from Bio import SeqIO
    from Bio import ExPASy
    from Bio import AlignIO
    from Bio.Align import Applications

    clustal_exe  = os.path.join(config.get('Paths','clustal_path'),'clustalw2')
    name_pdb = ".".join(pdb_name.split('/')[-1].split('.')[:-1])
    new_pdb=PDB()
    pdb_file=os.path.join(path,pdb_name)
    pdb=PDB(pdb_file)
    pdb.clean()
    for chain_id,chain_seq in sequences.iteritems():
       name_chain = name_pdb+"_"+chain_id
       name_seq   = chain_seq.get_identifier()
       pdb_chain  = pdb.get_chain_by_id(chain_id)
       new_chain  = Chain(name_pdb,chain_id)
       #define/create files
       infile =dummy_dir+"/tmp_"+name_chain+"_"+name_seq+".fa"
       outfile=dummy_dir+"/tmp_"+name_chain+"_"+name_seq+".aln"
       dndfile=dummy_dir+"/tmp_"+name_chain+"_"+name_seq+".dnd"
       fd=open(infile,"w")
       fd.write(">{0:s}\n{1:s}\n".format(name_chain,pdb_chain.protein_sequence))
       fd.write(">{0:s}\n{1:s}\n".format(name_seq,chain_seq.get_sequence()))
       fd.close()
       try:
         # run clustalw2
         msa_cline=Applications.ClustalwCommandline(clustal_exe,infile=infile,outfile=outfile)
         child = subprocess.Popen(str(msa_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell="/bin/bash")
         child.communicate()
         #store alignment in compare
         alignment=AlignIO.read(outfile,'clustal')
         structure=alignment[0].seq
         reference=alignment[1].seq
         try:
           len_3d =len(structure)
           len_ref=len(reference)
         except Exception as e:
           sys.stderr.write("ERROR: %s\n"%e)
           return e
       except Exception as e:
         sys.stderr.write("ERROR: %s\n"%e)
         return e
       #remove temporary fasta and alignment files
       remove_files([infile,outfile,dndfile])
       #mapping of residues to the original sequence
       mapping=create_mapping(pdb_chain.protein_idx.split(";"),structure,reference)
       #fill the new chain with the correct numbering of residues
       for residue in pdb_chain.aminoacids:
          pair = (str(residue.number),residue.version)
          number,version = mapping.get(pair)
          residue.number=number
          residue.version=version
          new_chain.add_residue(residue)
       #fill the new pdb
       new_pdb.add_chain(new_chain)
      
    return new_pdb  

def create_mapping(idx,structure,reference):
    mapped={}
    n=0
    m=0
    jump=0
    abc=" ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    for i in range(len(reference)):
      idx_reference=m+1
      if structure[i]!="-":
         number =idx[n][0:-1]
         version=idx[n][-1]
         if reference[i] == "-":
           if jump < 53: new_version = abc[jump]
           else: new_version="@"
         else:
           new_version = version
         mapped.setdefault((number,version),(idx_reference,new_version))
         n = n +1 
      if reference[i]!="-":
         jump = 0
         m    = m + 1
      else:
         jump = jump + 1
        
    return mapped 



