import math
import string
import sys
import os

#ALPHABETS


protein_letters = 'ACDEFGHIKLMNPQRSTVWYXUBZ'
rna_letters = 'GAUC'
dna_letters = 'GATC'


#WEIGHTS

protein_weights = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18,
                   'H': 155.16,'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13,
                   'S': 105.09,'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19, 'X': 0.0 ,
                   'U':0.0, 'B':0.0,'Z':0.0}

rna_weights = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}

dna_weights = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}

#COMPLEMENT

dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
rna_complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
dna_transcribe = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
rna_transcribe = {'A': 'T', 'C': 'G', 'G': 'C', 'U': 'A'}


#CODON TABLES

rna_table = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}

rna_table_back = {'A': 'GCU', 'C': 'UGU', None: 'UAA', 'E': 'GAG', 'D': 'GAU', 'G': 'GGU', 'F': 'UUU', 'I': 'AUU', 'H': 'CAU', 'K': 'AAG', 'M': 'AUG', 'L': 'UUG', 'N': 'AAU', 'Q': 'CAG', 'P': 'CCU', 'S': 'UCU', 'R': 'CGU', 'T': 'ACU', 'W': 'UGG', 'V': 'GUU', 'Y': 'UAU'}

rna_stop_codons = ['UAA', 'UAG', 'UGA']
rna_start_codons = ['UUG', 'CUG', 'AUG']



dna_table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}

dna_table_back = {'A': 'GCT', 'C': 'TGT', None: 'TAA', 'E': 'GAG', 'D': 'GAT', 'G': 'GGT', 'F': 'TTT', 'I': 'ATT', 'H': 'CAT', 'K': 'AAG', 'M': 'ATG', 'L': 'TTG', 'N': 'AAT', 'Q': 'CAG', 'P': 'CCT', 'S': 'TCT', 'R': 'CGT', 'T': 'ACT', 'W': 'TGG', 'V': 'GTT', 'Y': 'TAT'}

dna_stop_codons = ['TAA', 'TAG', 'TGA']
dna_start_codons = ['TTG', 'CTG', 'ATG']


# FUNCTIONS

def readframe(init,_sequence,start,stop,table):
  begin=False
  end=False
  proteins=[]
  seqprot=""
  for i in range(init,len(_sequence),3):
    codon=_sequence[i:i+3]
    if len(codon)<3:
      continue
    #print "ReadFrame %d Codon %s Seqprot %s \n"%(init,codon,seqprot)
    if codon in start and not begin:
      begin=True
      seqprot=""
    if begin and codon in stop:
      end=True
      begin=False
      if len(seqprot)>0:
        proteins.append(seqprot)
      continue
    if begin and not end:
      if codon in table:
        seqprot=seqprot+str(table[codon])
      else:
        raise ValueError("Codon %s not found in tables"%(codon))
  if len(seqprot)>0:
    proteins.append(seqprot)
  return proteins


def FASTA_iterator(fasta_filename):
 fd=open(fasta_filename,"r")
 seq=""
 name=""
 for line in fd:
   if line.startswith(">"):
    if len(name)>0 and len(seq)>0 :
        try:
          yield ProteinSequence(name,seq)
        except IncorrectSequenceLetter as e: #for python2.7 and up
          sys.stderr.write("%s"%(e))
          sys.stderr.write("Skip protein: %s\n"%(name))
    name=line.lstrip(">").strip()
    seq=""
   else:
    seq+=line.upper().strip()
 if len(name)>0 and len(seq)>0 :
  try:
    yield ProteinSequence(name,seq)
  except IncorrectSequenceLetter as e: #for python2.7 and up
    sys.stderr.write("%s\n"%e)
    sys.stderr.write("Skip protein: %s\n"%(name))
    
 fd.close()


# CLASSES

class IncorrectSequenceLetter(Exception):
  def __init__(self,letter,class_name):
    super(IncorrectSequenceLetter,self).__init__(
        "The _sequence items: '%s' are not found in the alphabet of class %s\n"%(letter,class_name))

class Sequence(object):
 alphabet=[]
 residue_mw={}
 def __init__(self,name=None,sequence=None):
  self._identifier=str(name)
  self._sequence=str(sequence)
  unknown=[w for w in self._sequence if w not in self.alphabet]
  if len(unknown)>0:
    word=",".join(set(unknown))
    raise IncorrectSequenceLetter(word,self.__class__.__name__)
 def __len__(self):
  return len(self._sequence)
 def __eq__(self,other):
  return self._sequence == other._sequence #This implies that we neglect different identifiers
  #return self._sequence == other._sequence and self.identifier == other.identifier
 def __ne__(self,other):
  return self._sequence != other._sequence
 def __le__(self,other):
  return self._sequence <= other._sequence
 def __ge__(self,other):
  return self._sequence >= other._sequence
 def __lt__(self,other):
  return self._sequence < other._sequence
 def __gt__(self,other):
  return self._sequence > other._sequence
 def __getitem__(self,item):
  return self._sequence[item]
 def __cmp__(self,other):
  if len(self._sequence) == len(other._sequence): return 0
  elif len(self._sequence) < len(other._sequence): return -1
  else: return 1
 def __hash__(self):
  #required to transform lists into set or dict (i.e. hashable elements)
  #this gives a unique id for a _sequence
  return self._sequence.__hash__()
 def __str__(self):
  return "I am a member of the Class %s"%(self.__class__.__name__)
 def __add__(self,other):
  return self.__class__(self.identifier,self._sequence+other._sequence)
 def count_tokens(self):
  return len(self._sequence)
 def get_identifier(self):
  return self._identifier
 def get_sequence(self):
  return self._sequence
 def get_mw(self):
  return sum(self.residue_mw.setdefault(aa,0) for aa in self._sequence)
 @property
 def identifier(self):
  return self._identifier
 @property
 def sequence(self):
  return self._sequence
 @property
 def mw(self):
  return sum(self.residue_mw.setdefault(aa,0) for aa in self._sequence)
 def has_sub_sequence(self,sub_sequence=None):
  if sub_sequence is None:
    raise ValueError("No sub_sequence was given")
  else:
    return sub_sequence.get_sequence() in self._sequence

class NucleotideSequence(Sequence):
 complement={}
 table={}
 stop=[]
 start=[]
 def get_complement(self):
  na= "".join([self.complement[u] for u in self._sequence])
  id=self.identifier
  return self.__class__(id,na)
 def translate(self):
  id=self.identifier
  proteins={}
  #proteins=set()
  for init in range(3):
    proteins.setdefault(init,readframe(init,self._sequence,self.start,self.stop,self.table))
    #proteins.update(readframe(init,self._sequence,self.start,self.stop,self.table))
  return proteins

class DNASequence(NucleotideSequence):
 alphabet=dna_letters
 residue_mw=dna_weights
 complement=dna_complement
 table=dna_table
 start=dna_start_codons
 stop=dna_stop_codons
 transcript=dna_transcribe
 def transcribe(self):
  trans="".join([self.transcript[u] for u in self.get_sequence()])
  id=self.get_identifier()
  return RNASequence(id,trans)

class RNASequence(NucleotideSequence):
 alphabet=rna_letters
 residue_mw=rna_weights
 complement=rna_complement
 table=rna_table
 start=rna_start_codons
 stop=rna_stop_codons
 transcript=rna_transcribe
 def reverse_transcribe(self):
  trans="".join([self.transcript[u] for u in self.get_sequence()])
  id=self.get_identifier()
  return DNASequence(id,trans)

class ProteinSequence(Sequence):
 alphabet=protein_letters
 residue_mw=protein_weights
 def reverse_translate_to_RNA(self):
  rnaseq= "".join([rna_table_back[u] for u in self.get_sequence()])
  id=self.get_identifier()
  return RNASequence(id,rnaseq)
 def reverse_translate_to_DNA(self):
  dnaseq= "".join([dna_table_back[u] for u in self.get_sequence()])
  id=self.get_identifier()
  return DNASequence(id,dnaseq)

