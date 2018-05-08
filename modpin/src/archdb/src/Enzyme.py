from SBI.databases import Enzyme

tables = {'enz'      : 'enzyme',
          'enz_ext'  : 'enzyme_extended',
          'enz_reac' : 'enzyme_reaction',
          'enz_com'  : 'enzyme_reaction_compounds',
          'enz_cop'  : 'enzyme_compound',
          'enz_cof'  : 'enzyme_cofactor',
          'enz_trs'  : 'enzyme_transfered', 
          'unienz'   : 'uniprot2enzyme'}

class Enzyme(Enzyme):
    def __init__(self, ec = None, description = None, inline = None):
        super(Enzyme, self).__init__(ec, description, inline)

    def toSQL(self):
        command = ""
        for cof in self.cofactors:
            command += "INSERT IGNORE INTO {0} SET name ='{1}', type='cofactor';\n".format(tables['enz_cop'], cof)
        for comp in self.compounds:
            command += "INSERT IGNORE INTO {0} SET name ='{1}', type='compound';\n".format(tables['enz_cop'], comp)

        if self.has_direct_parent and not self.has_transfers and not self.is_deleted:
            command += "INSERT INTO {0} VALUES ('{1.ec}','{1.description}',{1.level},'{1.direct_p}');\n".format(tables['enz'], self)
        else:
            command += "INSERT INTO {0} (id, description, level) VALUES ('{1.ec}','{1.description}',{1.level});\n".format(tables['enz'], self)
        
        if not self.has_transfers and not self.is_deleted:
            for parent in self.parents:
                command += "INSERT INTO {0} VALUES ('{1}','{2.ec}');\n".format(tables['enz_ext'], parent, self)

        for cofactor in self.cofactors:
            command += "SET @cof = (SELECT nid FROM {0} WHERE name='{1}' AND type='cofactor');\n".format(tables['enz_cop'], cofactor)
            command += "INSERT INTO {0} VALUES ('{1.ec}', @cof);\n".format(tables['enz_cof'], self, cofactor)

        for reaction in self.reactions:
            if not reaction.repeated:
                command += "INSERT INTO {0} (ec, description) VALUES ('{1.ec}','{2.description}');\n".format(tables['enz_reac'], self, reaction)
            else:
                command += "SET @old_id = (SELECT nid FROM {0} WHERE description='{1.description}' LIMIT 1);\n".format(tables['enz_reac'], reaction)
                command += "INSERT INTO {0} VALUES (@old_id, '{1.ec}','{2.description}');\n".format(tables['enz_reac'], self, reaction)
            for subst in reaction.substrates:
                command += "SET @compnum = (SELECT nid FROM {0} WHERE name='{1[0]}' AND type='compound');\n".format(tables['enz_cop'], subst)
                command += "INSERT INTO {0} VALUES (LAST_INSERT_ID(),@compnum,'substrate',{1[1]});\n".format(tables['enz_com'], subst)
            for prod in reaction.products:
                command += "SET @compnum = (SELECT nid FROM {0} WHERE name='{1[0]}' AND type='compound');\n".format(tables['enz_cop'], prod)
                command += "INSERT INTO {0} VALUES (LAST_INSERT_ID(),@compnum,'product',{1[1]});\n".format(tables['enz_com'], prod)
        for protein in self.uniprots:
            command += "INSERT INTO {0} VALUES('{1}','{2}');\n".format(tables['unienz'], protein, self.ec)
        return command

    def transfered2SQL(self):
        command = ""
        for trans in self.transfers:
            command += "INSERT INTO {0} VALUES ('{1.ec}','{2}');\n".format(tables['enz_trs'], self, trans)
        return command
