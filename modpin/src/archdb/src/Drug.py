from SBI.databases import Drug

tables = {'drg'      : 'drugBank',
          'drg_cat'  : 'drugBank_category',
          'drg_ext'  : 'drugBank_external_identifiers',
          'drg_grp'  : 'drugBank_group',
          'drg_san'  : 'drugBank_secondary_accession_number',
          'drg_syn'  : 'drugBank_synonym',
          'drg_tar'  : 'drugBank_target', 
          'drg_act'  : 'drugBank_target_action',
          'uni'      : 'uniprot_entry2accession',
          'uniprot'  : 'uniprot'}

class Drug(Drug):
    def __init__(self, inline):
        super(Drug, self).__init__(inline = inline)

    @staticmethod
    def preuniprotdeleted():
        return "INSERT INTO {0} (entry, source) VALUES ('FAKE_PROTEIN', 'swissprot');\n".format(tables['uniprot'])
    @staticmethod
    def afteruniprotdeleted():
        command  = "DELETE FROM {0} WHERE uniprot = 'FAKE_PROTEIN';".format(tables['drg_act'])
        command += "DELETE FROM {0} WHERE uniprot = 'FAKE_PROTEIN';".format(tables['drg_tar'])
        command += "DELETE FROM {0} WHERE entry = 'FAKE_PROTEIN';".format(tables['uniprot'])
        return command

    def toSQL(self):
        h  =     self.dproperty["Hydrophobicity"]        if "Hydrophobicity"    in self.dproperty else 'NULL'
        i  =     self.dproperty["Isoelectric Point"]     if "Isoelectric Point" in self.dproperty else 'NULL'
        mw =     self.dproperty["Molecular Weight"]      if "Molecular Weight"  in self.dproperty else 'NULL'
        mf = "'"+self.dproperty["Molecular Formula"]+"'" if "Molecular Formula" in self.dproperty else 'NULL'
        cas_number          = "'"+self.cas+"'"                     if self.cas  != '' else 'NULL'
        description         = "'"+self.desc.replace("'",r'\'')+"'" if self.desc != '' else 'NULL'
        name                = self.name.replace("'",r'\'')

        command = ""
        command += "INSERT INTO {0}(drugbank_id,name,type,cas_number,description,hydrophobicity,isoelectric_point,molecular_weight,molecular_formula) ".format(tables['drg'])
        command += "VALUES ('{0.dbid}','{7}','{0.type}',{1},{2},{3},{4},{5},{6});\n".format(self, cas_number,description,h,i,mw,mf, name)
        for c in self.categories:
            command += "INSERT INTO {0} (drugbank_id, drugbank_category) VALUES ('{1.dbid}','{2}');\n".format(tables['drg_cat'], self, c)
        for e in self.external:
            command += "INSERT INTO {0} (drugbank_id, identifier, resource) VALUES ('{1.dbid}','{2[1]}','{2[0]}');\n".format(tables['drg_ext'],self,e)
        for g in self.groups:
            command += "INSERT INTO {0} (drugbank_id, drugbank_group) VALUES ('{1.dbid}','{2}');\n".format(tables['drg_grp'], self, g)
        for s in self.secondary_accession_numbers:
            command += "INSERT INTO {0} (drugbank_id, secondary_accession_number) VALUES ('{1.dbid}','{2}');\n".format(tables['drg_san'], self, s)
        for y in self.synonym:
            command += "INSERT INTO {0} (drugbank_id,drugbank_synonym) VALUES ('{1.dbid}','{2}');\n".format(tables['drg_syn'], self, y.replace("'",r'\''))
        for t in self.partners:
            if t[1] != ('',''):
                command += "SET @uni = (SELECT COALESCE((SELECT entry FROM uniprot WHERE entry = '{0[1]}'),(SELECT entry FROM uniprot_entry2accession WHERE accession='{0[0]}' LIMIT 1),'FAKE_PROTEIN'));\n".format(t[1])
                known = 'yes' if t[-1] == 'yes' else 'unknown' if t[-1] == 'unknown' else 'no'
                command += "INSERT INTO {0} (drugbank_id,uniprot,known_action,target_type) VALUES ('{1.dbid}',@uni,'{2}','{3}');\n".format(tables['drg_tar'], self, known, t[0])

                for i in range(2,len(t)):
                    if t[i] != 'yes' and t[i] != 'unknown':
                        command += "INSERT INTO {0} (drugbank_id,uniprot,action) VALUES ('{1.dbid}',@uni, '{2}');\n".format(tables['drg_act'], self, t[i])
        return command

