# from SBI.structure import PDB, ResidueOfNucleotide, ResidueOfAminoAcid
# from SBI.contacts  import PPInterface, PNInterface

from SBI.structure           import PDB
from SBI.structure.contacts  import Complex, InnerContacts
from SBI                         import SBIglobals

tables = {'main'        : 'PDB',
          'old'         : 'oldPDB',
          'chain'       : 'chain',
          'ec'          : 'chain2enzyme',
          'enz'         : 'enzyme',
          'ctaxid'      : 'chain2taxid',
          'taxid'       : 'taxid',
          'taxidold'    : 'taxid_old',
          'uniprot'     : 'uniprot',
          'pdbuniprot'  : 'chain2uniprot',
          'ent2acc'     : 'uniprot_entry2accession',
          'hetatm'      : 'hetero',
          'hetchn'      : 'hetero2PDB',
          'repeidx'     : 'chainrepeatidx',
          'repe'        : 'chainrepeat',
          'site'        : 'site',
          'sitepos'     : 'site2position',
          'h6c'         : 'hetero_contacts6',
          'cont12'      : 'contacts12',
          'chainc12'    : 'chain_contacts12'}

class InnerContacts(InnerContacts):
    def toSQL(self):
        SBIglobals.alert('debug', self, 'Printing Inner Contacts (Protein-Heteroatom) as ArchDB SQL')
        command = ''
        for ht in self.HTcontacts:
            if len(ht) == 0: continue
            SBIglobals.alert('debug', self, '\tInner Contact {0} for chain {1}'.format(ht.identifier, ht.protein.chain))
            command += "SET @chain = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.pdb.id, ht.protein.chain)
            for hetero in ht.heteroatoms_view_innercontact:
                SBIglobals.alert('debug', self, '\tFound {0} contacts for HETATM {1.type} in {1.number}'.format(len(ht.heteroatoms_view_innercontact[hetero]), hetero))
                command += "SET @heteroID = gethetero2PDB(@chain,{0.number},'{0.type}');\n".format(hetero)
                for contact in ht.heteroatoms_view_innercontact[hetero]:
                    command += "INSERT INTO {0}(chain,hetero,position,idxp,type,distance) VALUES (@chain,@heteroID,{2.number},'{2.version}','{2.type}',{1:.3f});\n".format(tables['h6c'], contact.min_distance,contact.aminoacid)
        return command

class Complex(Complex):
    def toSQL(self):
        SBIglobals.alert('debug', self, 'Printing Interfaces as ArchDB SQL')
        command = ''
        for ppi in self.PPInterfaces:
            if len(ppi) == 0: continue
            SBIglobals.alert('debug', self, '\tProtein - Protein Interface {0} for chains {1} - {2}'.format(ppi.identifier, ppi.protein_chain.chain, ppi.protein_interactor.chain))
            command += "SET @chainA = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.pdb.id, ppi.protein_chain.chain)
            command += "SET @chainB = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.pdb.id, ppi.protein_interactor.chain)
            SBIglobals.alert('debug', self, '\tFound a total of {0} contacts'.format(len(ppi.contacts)))
            for contact in ppi.contacts:
                command += "INSERT INTO {0} (alpha_distance, beta_distance, min_distance, type) VALUES ({1:.3f},{2:.3f},{3:.3f},'PPI');\n".format(tables['cont12'], contact.ca_distance, contact.cb_distance, contact.min_distance)
                command += "SET @contact = LAST_INSERT_ID();\n";
                command += "INSERT INTO {0} (contacts12, chain, position, idxp, type) VALUES (@contact,@chainA,{1.number},'{1.version}','{1.type}');\n".format(tables['chainc12'], contact.aminoacid1)
                command += "INSERT INTO {0} (contacts12, chain, position, idxp, type) VALUES (@contact,@chainB,{1.number},'{1.version}','{1.type}');\n".format(tables['chainc12'], contact.aminoacid2)

        for pni in self.PNInterfaces:
            if len(pni) == 0: continue
            SBIglobals.alert('debug', self, '\tProtein - Nucleotide Interface {0} for chains {1} - {2}'.format(pni.identifier, pni.protein.chain, pni.nucleotide.chain))
            command += "SET @chainP = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.pdb.id, pni.protein.chain)
            command += "SET @chainN = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.pdb.id, pni.nucleotide.chain)
            SBIglobals.alert('debug', self, '\tFound a total of {0} contacts'.format(len(pni.contacts)))
            for contact in pni.contacts:
                command += "INSERT INTO {0} (beta_distance, min_distance, type) VALUES ({1:.3f},{2:.3f},'PNI');\n".format(tables['cont12'], contact.cbbackbone_distance, contact.min_distance)
                command += "SET @contact = LAST_INSERT_ID();\n";
                command += "INSERT INTO {0} (contacts12, chain, position, idxp, type) VALUES (@contact,@chainP,{1.number},'{1.version}','{1.type}');\n".format(tables['chainc12'], contact.aminoacid)
                command += "INSERT INTO {0} (contacts12, chain, position, idxp, type) VALUES (@contact,@chainN,{1.number},'{1.version}','{1.type}');\n".format(tables['chainc12'], contact.nucleotide)
        for phi in self.PHInterfaces:
            if len(phi) == 0: continue
            SBIglobals.alert('debug', self, '\tProtein - Heteroatom Interface {0} for chains {1} - {2}'.format(phi.identifier, phi.protein.chain, phi.heteroatoms.chain))
            command += "SET @chainP = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.pdb.id, phi.protein.chain)
            command += "SET @chainH = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.pdb.id, phi.heteroatoms.chain)
            for hetero in phi.heteroatoms_view_interface:
                SBIglobals.alert('debug', self, '\tFound {0} contacts for HETATM {1.type} in {1.number} on {2} with chain {3}'.format(len(phi.heteroatoms_view_interface[hetero]), hetero, phi.heteroatoms.chain, phi.protein.chain))
                command += "SET @heteroID = gethetero2PDB(@chainH,{0.number},'{0.type}');\n".format(hetero)
                for contact in phi.heteroatoms_view_interface[hetero]:
                    command += "INSERT IGNORE INTO {0}(chain,hetero,position,idxp,type,distance) VALUES (@chainP,@heteroID,{2.number},'{2.version}','{2.type}',{1:.3f});\n".format(tables['h6c'], contact.min_distance,contact.aminoacid)
        return command

class PDB(PDB):
    def __init__(self, pdb_file = None):
        SBIglobals.alert('debug', self, 'Loading PDB file {0}'.format(pdb_file))

        super(PDB, self).__init__(pdb_file = pdb_file, dehydrate = True, header = True)

        SBIglobals.alert('debug', self, 'Calculating Inner Contacts for Protein - Heteroatom')
        self.innercontacts = InnerContacts(pdb = self, AA = False,
                                                       NC = False,
                                                       HT = True,  HT_type = "min", HT_distance = 6)

        SBIglobals.alert('debug', self, 'Calculating PP and PN interfaces of the biomolecules in the PDB')
        self.interfaces    = Complex(pdb = self, biomolecule = True, PPI = True, PPI_type = "cb",  PPI_distance = 12,
                                                                     PNI = True, PNI_type = "min", PNI_distance = 8,
                                                                     PHI = True, PHI_type = "min", PHI_distance = 6)

    @staticmethod
    def preuniprotdeleted():
        return "INSERT INTO {0} (entry, source) VALUES ('FAKE_PROTEIN', 'swissprot');\n".format(tables['uniprot'])
    @staticmethod
    def afteruniprotdeleted():
        command  = "DELETE FROM {0} WHERE uniprot = 'FAKE_PROTEIN';\n".format(tables['pdbuniprot'])
        command += "DELETE FROM {0} WHERE entry = 'FAKE_PROTEIN';\n".format(tables['uniprot'])
        return command

    @staticmethod
    def _choose_ec(ec):
        command = []
        data = ec.split('.')
        main = "(SELECT id FROM {0} WHERE id = '{1}')"
        if len(data) == 3: data.append('-')
        command.append(main.format(tables['enz'], '.'.join(data)))
        if data[3] != '-':
            data[3] = '-'
            command.append(main.format(tables['enz'], '.'.join(data)))
        if data[2] != '-':
            data[2] = '-'
            command.append(main.format(tables['enz'], '.'.join(data)))
        if data[1] != '-':
            data[1] = '-'
            command.append(main.format(tables['enz'], '.'.join(data)))
        return 'SET @ec = (COALESCE('+ ",".join(command)+'));\n'

    def toSQL(self):
        pdbheader = self.header

        #if pdbheader.valid_resolution and pdbheader.resolution != 'NULL':
        if pdbheader.experiment.resolution >0 :
            command = "INSERT INTO {0} VALUES ('{1.id}','{2.date}','{2.header}','{2.xpdta}',{2.resolution:.2f},{2.rfactor:.3f},{2.freeR:.3f});\n".format(tables['main'], self, pdbheader)
        else:
            command = "INSERT INTO {0} (pdb,date,header,method) VALUES ('{1.id}','{2.date}','{2.header}','{2.xpdta}');\n".format(tables['main'], self, pdbheader)
        for deprec in pdbheader.deprecated:
            command += "INSERT INTO {0} VALUES ('{1}','{2.id}');\n".format(tables['old'], deprec, self)
        for chain in self.chains:
            m = pdbheader.get_molecule4chain(chain.chain)
            f = chain.first_structure
            l = chain.last_structure
            command += "INSERT INTO {0}(pdb,chain,name,type,start,idxs,end,idxe) VALUES ('{1.pdb}','{1.chain}','{2.name}','{1.chaintype}',{3.number},'{3.version}',{4.number},'{4.version}');\n".format(tables['chain'],chain,m,f,l)
            command += "SET @chain = LAST_INSERT_ID();\n";
            for ec in m.ec:
                command += PDB._choose_ec(ec)
                command += "INSERT INTO {0} VALUES (@chain,@ec);\n".format(tables['ec'])
            for tx in m.taxid:
                command += "SET @taxid = (SELECT COALESCE((SELECT taxid FROM {0} WHERE oldid={2}),(SELECT id FROM {1} WHERE id={2})));\n".format(tables['taxidold'], tables['taxid'], tx)
                command += "INSERT INTO {0} VALUES (@chain,@taxid);\n".format(tables['ctaxid'])
        for dbr in pdbheader.dbrefs:
            if dbr.chain in self._chain_id:
                command += "SET @chain = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.id, dbr.chain)
                command += "SET @uniprot = (COALESCE((SELECT entry FROM {0} WHERE accession='{1}' LIMIT 1), 'FAKE_PROTEIN'));\n".format(tables['ent2acc'], dbr.uniprot)
                command += "INSERT INTO {0} VALUES (@chain,@uniprot,{1.start},'{1.idxs}',{1.end},'{1.idxe}');\n".format(tables['pdbuniprot'], dbr)
        for h in pdbheader.hetero:
            if h.chain in self._chain_id:
                command += "INSERT IGNORE INTO {0} VALUES ('{1.id}','{1._name}','{1.form}', 0, 0);\n".format(tables['hetatm'], h)
                command += "SET @chain = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.id, h.chain)
                dataset = 0
                if int(h.pos) <= self.get_chain_by_id(h.chain).last_structure.number: dataset = 1
                command += "INSERT INTO {0}(chain,position,hetero,inchain) VALUES (@chain, {1.pos}, '{1.id}', {2});\n".format(tables['hetchn'], h, dataset)
        for m in pdbheader.molecules:
            command += "INSERT INTO {0} VALUES ();\n".format(tables['repeidx'])
            for c in pdbheader.molecules[m].chains:
                if c in self._chain_id:  #sometimes a chain is described in the header that has no relation in the coordinates
                    command += "SET @chain = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.id, c)
                    command += "INSERT INTO {0} VALUES (LAST_INSERT_ID(), @chain);\n".format(tables['repe'])
        for s in pdbheader.sites:
            if pdbheader.sites[s].bind is not None:
                bdata = pdbheader.sites[s].bind
                if not isinstance(bdata[1], list):
                    bdata[1] = bdata[1].split()
                    if len(bdata[1]) == 1:
                        bdata[1].append(bdata[1][0][1:].strip())
                        bdata[1][0] = bdata[1][0][0]
                if (len(bdata[1]) == 0): #only one het of this type in the pdb
                    for e in pdbheader.hetero:
                        if bdata[0] == e.id:
                            bdata[1].append(e.chain)
                            bdata[1].append(e.pos)
                import re
                posidx = re.compile('(\-*\d+)(\w*)')
                m      = posidx.match(bdata[1][1])
                pos    = m.group(1)
                idx    = m.group(2)
                command += "SET @chain = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2[1][0]}');\n".format(tables['chain'], self.id, bdata)
                command += "SET @bind = (SELECT nid FROM {0} WHERE chain=@chain AND position={3} AND hetero='{2[0]}');\n".format(tables['hetchn'], self.id, bdata, pos)
                command += "INSERT INTO {0}(pdb,name,description,bind) VALUES ('{1}','{2.id}','{2.desc}', @bind);\n".format(tables['site'],self.id,pdbheader.sites[s])
            else:
                command += "INSERT INTO {0}(pdb,name,description) VALUES ('{1}','{2.id}','{2.desc}');\n".format(tables['site'],self.id,pdbheader.sites[s])
            for p in pdbheader.sites[s].spec:
                command += "SET @chain = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self.id, p[1])
                idx = ' '
                num = p[2]
                try:
                    int(num[-1])
                    idx = ' '
                    num = p[2]
                except:
                    idx = num[-1]
                    num = num[:-1]
                command += "INSERT IGNORE INTO {0} VALUES(LAST_INSERT_ID(),@chain,{1},'{2}','{3}');\n".format(tables['sitepos'],num, idx, p[0])

        command += self.innercontacts.toSQL()
        command += self.interfaces.toSQL()

        return command
