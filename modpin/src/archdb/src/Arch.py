# from SBI.structure.protein import SSpair
import re
tables = {'chain' : 'chain',
          'ld'    : 'loop_description',
          'sld'   : 'superloop_description',
          'l2c'   : 'loop2chain',
          's2c'   : 'superloop2chain',
          'l2s'   : 'loop2superloop'}


class Arch(object):

    def __init__(self, pdb):
        d = pdb.split('_')
        self._pdb   = d[0]
        self._chain = d[1]
        self._archs = []

    @property
    def archs(self):
        return self._archs

    @archs.setter
    def archs(self, value):
        self._archs.append(value)

    def toSQL(self):
        command = "SET @chain = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], self._pdb, self._chain)

        singles = {}
        supers  = {}
        for arch in self._archs:
            # command += str(arch)
            l1  = arch.length_ss1
            l2  = arch.length_ss2
            mi1 = arch._ss1.get_moment_of_inertia_length('f11')
            mi2 = arch._ss2.get_moment_of_inertia_length('f44')
            sec = arch.aminoacid_sequence
            exp = arch.access_surface
            ss  = arch.structure_sequence
            i   = re.search('(-*\d+)(\w*)', arch.initial_position)
            e   = re.search('(-*\d+)(\w*)', arch.end_position)
            ii, ix, ee, ex = i.group(1), ' ' if i.group(2) == '' else i.group(2), e.group(1), ' ' if e.group(2) == '' else e.group(2)
            if not arch.is_superarch:
                table  = tables['ld']
                values = "(loop_id,type,length,ss1L,ss1moiL,ss2L,ss2moiL,distance,theta,rho,delta,sequence,ss,exposition)"
                close  = ");\n"
                archid = arch._source + '_' + str(arch.initial_position)
                singles[int(ii)] = archid.strip()
            else:
                table  = tables['sld']
                values = "(loop_id,type,length,ss1L,ss1moiL,ss2L,ss2moiL,distance,theta,rho,delta,sequence,ss,exposition,internal_ss_count,internal_ss_desc)"
                close  = ",{0._intss},'{0._inttyp}');\n".format(arch)
                archid = arch._source + '_' + arch.initial_position + '_' + str(arch._intss)
                supers[str(ii) + '_' + str(arch._intss)] = archid.strip()
            command += "INSERT INTO {0} {1} VALUES ('{11}','{2.type}',{2._distAA},{4},{5},{6},{7},{2._dist:.3f},{2._angls[0]:.3f},{2._angls[1]:.3f},{2._angls[2]:.3f},'{8}','{9}','{10}'{3}".format(table, values, arch, close, l1, mi1, l2, mi2, sec, ss, exp, archid.strip())
            command += "SET @arch = LAST_INSERT_ID();\n"

            if not arch.is_superarch:
                command += "CALL new_loop2chain(@arch,@chain,{2},'{3}',{4},'{5}','D',{6});\n".format(self._pdb, self._chain, ii, ix, ee, ex, arch._order)
            else:
                command += "CALL new_superloop2chain(@arch,@chain,{2},'{3}',{4},'{5}','D',{6});\n".format(self._pdb, self._chain, ii, ix, ee, ex, arch._order)

        for s in sorted(supers):
            info = s.split('_')
            info[0] = int(re.search('(-*\d+)(\w*)', info[0]).group(1))
            info[1] = int(info[1])
            command += "SET @super = (SELECT nid FROM {0} WHERE loop_id='{1}');\n".format(tables['sld'], supers[s])
            toprint = False
            for a in sorted(singles):
                if info[1] < 0:
                    break
                if a == info[0]:
                    toprint = True
                if toprint:
                    command += "SET @single = (SELECT nid FROM {0} WHERE loop_id='{1}');\n".format(tables['ld'], singles[a])
                    command += "INSERT INTO {0} VALUES (@super,@single);\n".format(tables['l2s'])
                    info[1] -= 1
        return command
