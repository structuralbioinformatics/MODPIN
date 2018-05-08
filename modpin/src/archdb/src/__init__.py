def start_transaction(): return "SET autocommit = 0;\nSTART TRANSACTION;\n"
def end_transaction():   return "COMMIT;\n"

from Source   import Source
from PDBeChem import Element, PDBeChem
from GOterm   import GOterm
from TaxID    import TaxID
from Uniprot  import Uniprot
from Enzyme   import Enzyme
from PDB      import PDB
from Drug     import Drug
from SCOP     import SCOP
from TM       import TM
from CDhit    import CDhit
from Arch     import Arch
from Subclass import Subclass, Cclass, Loop