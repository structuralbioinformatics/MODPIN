import re
class Subclass(object):

    def __init__(self, seed, length, lrange = None):
        self._seed   = seed
        if lrange is None or str(length) != '4':
            self._length = length
        else:
            if lrange == '00-04': self._length = 4
            else:                 self._length = 4.5
        self._loops  = []
        self.range   = lrange
        self._identifier = None

        self.dist_range  = None
        self.delta_range = None
        self.theta_range = None
        self.rho_range   = None

        self.dist  = None
        self.delta = None
        self.theta = None
        self.rho   = None

        self._ram_pat = None
        self._seq_pat = None
        self._exp_pat = None

        self.consensus = None

        self.subcode = 0
        self.subclassID = 'NULL'

    @property 
    def length(self): 
        return float(self._length)
        # if str(self._length) not in ['4A','4B']:
        #     return int(self._length)
        # else: return 4
    @property 
    def defining_length(self):
        return self._length
    @property 
    def size(self): return len(self._loops)
    @property 
    def loops(self): return self._loops
    @loops.setter
    def loops(self, loop): 
        self._loops.append(loop)
        self._loops[-1].order = self.length

    @property
    def identifier(self):
        if self._identifier is None:
            self._identifier = "_".join([y for y in sorted(["_".join([x.loopID,x.position]) for x in self.loops])])
        return self._identifier

    @property
    def ram_pat(self):
        return self._ram_pat
    @ram_pat.setter
    def ram_pat(self, value):
        self._ram_pat = self._motif(full = value)

    @property
    def seq_pat(self):
        return self._seq_pat
    @seq_pat.setter
    def seq_pat(self, value):
        self._seq_pat = self._motif(full = value)

    @property
    def exp_pat(self):
        return self._exp_pat
    @exp_pat.setter
    def exp_pat(self, value):
        self._exp_pat = self._motif(full = value)

    def add_topology(self, info, method):

        dist_code = {"A":'0-2', "B":'2-4', "C":'4-6', "D":'6-8', "E":'8-10', "F":'10-12', "G":'12-14', "H":'14-16', "I":'16-18',
                     "J":'18-20', "K":'20-22', "L":'22-24', "M":'24-26', "N":'26-28', "O":'28-30', "P":'30-32', "Q":'32-34', "R":'34-36', "S":'36-38',
                     "T":'38-40', "?":'40-42'}

        rho_code = {"A":'0-45', "B":'45-90', "C":'90-135', "D":'135-180', "E":'180-225', "F":'225-270', "G":'270-315', "H":'315-360', "?":'315-360'}

        #METHOD DEF
        if method == 'DS': #DENSITY SEARCH
            delta_code = {"A":'0-45', "B":'45-90', "C":'90-135', "D":'135-180', "E":'180-225', "F":'225-270', "G":'270-315', "H":'315-360', "?":'135-180',}
            theta_code = {"A":'0-45', "B":'45-90', "C":'90-135', "D":'135-180', "E":'180-225', "F":'225-270', "G":'270-315', "H":'315-360'}
        elif method == 'MCL': #MCL
            delta_code = {"A":'0-30', "B":'30-60', "C":'60-90', "D":'90-120', "E":'120-180', "F":'180-210', "G":'210-240', "H":'240-270',
                          "I":'270-300', "J":'300-330', "K":'330-360'}
            theta_code = {"A":'0-30', "B":'30-60', "C":'60-90', "D":'90-120', "E":'120-180', "F":'180-210', "G":'210-240', "H":'240-270',
                          "I":'270-300', "J":'300-330', "K":'330-360'}


        topology = info.split(':')[1].strip()

        self.dist_range  = dist_code[topology[0]]
        self.delta_range = delta_code[topology[1]]
        self.theta_range = theta_code[topology[2]]
        self.rho_range   = rho_code[topology[3]]

    def add_coordinates(self, info):

        data = info.split(':')[1].strip().split()
        self.dist  = data[2]
        self.delta = data[5]
        self.theta = data[8]
        self.rho   = data[11]

    def _motif(self, full):
        data = [re.sub(r',.\Z','',re.sub(r'[\(|\)]','',x.split('=')[-1])) for x in re.split('\)\(',full)]
        for d in range(len(data)):
            if re.search(r',',data[d]):
                data[d] = "[" + data[d] + "]"
        data = "-".join(data)

        return data

    def add_consensus(self, info, method, loops):
        INFO = info
        info = info.split(':')
        ini = None
        end = None
        if method == 'DS':
            l1 = len(info[0]) + 1 + len(info[1]) - len(info[1].lstrip())
            l2 = l1 + len(info[1].strip())
            z1 = loops[self._seed]
            the_first = True
            for loop in self.loops:
                loop.fix_align(l1 - loop.refer, l2 - loop.refer)
                if the_first:
                    # print "CON", ":".join(info)[loop.refer[0]-1:loop.refer[1]-1]
                    # print info[1]
                    # print INFO
                    # print INFO[loop.prel + z1 -2:loop.prel + z1 +2 + self.length]
                    # self.consensus = ":".join(info)[loop.refer[0]-1:loop.refer[1]-1]
                    # print consensus
                    # self.consensus = self.consensus[z1-2: z1 + 2 + self.length]
                    self.consensus = INFO[loop.prel + z1 -2:loop.prel + z1 +2 + int(self.length)]
                    the_first = False
                # if loop.id in loop_data:
                #     ini = loop.prel + loop.iniseq + loop_data[loop.id] - 2
                #     end = loop.prel + loop.iniseq + loop_data[loop.id] + loop.length + 2
                #     break
        if method == 'MCL':
            l1 = len(info[0]) + 1 + len(info[1]) - len(info[1].lstrip())
            l2 = l1 + len(info[1].strip())
            z0 = []
            the_first = True
            for loop in self.loops:
                z1 = loops[loop.loop_ID_pos()]
                loop.fix_align(l1 - loop.refer, l2 - loop.refer)
                if the_first:
                    z0 = [loop.prel + z1 -4, loop.prel + z1 +4 + int(self.length)]
                    the_first = False
                # else:
                #     if loop.prel + z1 -4 > z0[0]: z0[0] = loop.prel + z1 -4
                #     if loop.prel + z1 +4 + self.length < z0[1]: z0[1] = loop.prel + z1 +4 + self.length

            self.consensus = INFO[z0[0]:z0[1]]
                # the_first = False

    def toSQL(self, classcode, counter_subclass):
        subclassID = ".".join([classcode, str(counter_subclass)])
        command  = "INSERT INTO cluster_subclass (class_nid,subclass,name,size,seq_pat,exp_pat,ram_pat)"
        command += " VALUES (@classID,{2},'{1}',{0.size},'{0.seq_pat}','{0.exp_pat}','{0.ram_pat}');\n".format(self, subclassID, counter_subclass)
        command += "SET @subclassID = LAST_INSERT_ID();\n"
        loop_counter = 1
        for loop in self.loops:
            command += loop.toSQL(loop_counter)
            loop_counter += 1
        command += self.get_loops_params()
        return command


    def get_loops_params(self):
        command  = "SELECT CONCAT(FLOOR(MIN(ld.distance)),'-',CEIL(MAX(ld.distance))), "
        command += "CEIL(MAX(ld.distance)), FLOOR(MIN(ld.distance)),"
        command += "CONCAT(FLOOR(MIN(ld.theta)),'-',CEIL(MAX(ld.theta))), "
        command += "CEIL(MAX(ld.theta)), FLOOR(MIN(ld.theta)), "
        command += "CONCAT(FLOOR(MIN(ld.rho)),'-',CEIL(MAX(ld.rho))), "
        command += "CEIL(MAX(ld.rho)), FLOOR(MIN(ld.rho)),"
        command += "CONCAT(FLOOR(MIN(ld.delta)),'-',CEIL(MAX(ld.delta))), "
        command += "CEIL(MAX(ld.delta)),FLOOR(MIN(ld.delta)) "
        command += "FROM cluster_subclass cs, loop2cluster lc, loop_description ld WHERE "
        command += "cs.nid=@subclassID AND lc.cluster_nid=cs.nid AND ld.nid=lc.loop_nid "
        command += "INTO @Sr, @Sx, @Sm, @Tr, @Tx, @Tm, @Rr, @Rx, @Rm, @Dr, @Dx, @Dm;\n"
        command += "UPDATE cluster_subclass SET dist_range=@Sr, theta_range=@Tr, rho_range=@Rr, delta_range=@Dr, "
        command += "dist_range_min=@Sm, theta_range_min=@Tm, rho_range_min=@Rm, delta_range_min=@Dm, "
        command += "dist_range_max=@Sx, theta_range_max=@Tx, rho_range_max=@Rx, delta_range_max=@Dx "
        command += "WHERE nid=@subclassID;\n"
        return command

    def __cmp__(self, other):
        return self.size.__cmp__(other.size)
        

    # def add_consensus(self, info, loop_data, method):
    #     ini = None
    #     end = None
    #     if method == 'DS':
    #         for loop in self.loops:
    #             if loop.id in loop_data:
    #                 ini = loop.prel + loop.iniseq + loop_data[loop.id] - 2
    #                 end = loop.prel + loop.iniseq + loop_data[loop.id] + loop.length + 2
    #                 break
    #     elif method == 'MCL':
    #         for loop in self.loops:
    #             if loop.id in loop_data:
    #                 new_ini = loop.prel + loop.iniseq + loop_data[loop.id]
    #                 new_end = loop.prel + loop.iniseq + loop_data[loop.id] + loop.length
    #                 if ini is None:
    #                     ini = new_ini
    #                     end = new_end
    #                 else:
    #                     if new_ini < end:
    #                         if new_ini > ini: ini = new_ini
    #                         if new_end < end: end = new_end
    #         if ini is not None:
    #             ini = ini - 4
    #             end = end + 4
    #         else:
    #             ini = 0
    #             end = -1
    #     self.consensus = info[ini:end]
class Consensus(object):
    def __init__(self, length, consensus, gtype, subclasses, method):
        self.consensus = consensus
        self.length = length
        self.gtype = gtype
        self.subclasses = subclasses
        self.method = method
        self._size = -1
    @property
    def size(self):
        if self._size == -1:
            self._size = 0
            for s in self.subclasses:
                self._size += s.size
        return self._size

    def __cmp__(self, other):
        return self.size.__cmp__(other.size)

    def toSQL(self, classcode):
        length = self.length
        if self.method == 'MCL' and length == 4: length = '4S'
        elif self.method == 'MCL' and length == 4.5: length = '4M'
        else: length = int(length)
        # length  = self.length if self.length not in ['4A','4B'] else '4S' if self.length == '4A' else '4M'
        classID = ".".join([self.gtype,str(length),str(classcode)]) 
        command  = "INSERT INTO cluster_class (method,type,length,class,name,size,consensus) "
        command += "VALUES (@databasemethod,'{0.gtype}','{3}',{1},'{2}',{0.size},'{0.consensus}');\n".format(self, classcode, classID, length)
        command += "SET @classID = LAST_INSERT_ID();\n"
        counter_subclass = 1
        for s in self.subclasses:
            command += s.toSQL(classID,counter_subclass)
            counter_subclass += 1
        return command

class Cclass(object):
    def __init__(self, kltype):
        self._kltype       = kltype
        self._subclasses   = {}
        self._lastsubclass = None

    @property 
    def subclasses(self): return self._subclasses
    @subclasses.setter
    def subclasses(self, value): 
        self._subclasses.setdefault(value.length,[]).append(value)
        self._lastsubclass = value

    @property 
    def lastsubclass(self): return self._lastsubclass

    def toSQL(self, method):
        info = {}
        consdata = {}
        donesubclass = set()
        command = "SELECT nid FROM method WHERE name='{0}' INTO @databasemethod;\n".format(method)
        for l in self.subclasses:
            info.setdefault(l,{})
            for s in self.subclasses[l]:
                if method == 'MCL' and s.identifier in donesubclass: continue
                if s.size >= 3:
                    info[l].setdefault(s.consensus,[]).append(s)
                    if method == 'MCL': donesubclass.add(s.identifier)
            for c in info[l]:
                consdata.setdefault(l,[]).append(Consensus(l,c,self._kltype, sorted(info[l][c], reverse=True), method))
                # info[l][c] = sorted(info[l][c], reverse=True)


        for l in sorted(consdata):
            class_count = 1
            # print 'LENGTH ', l
            consdata[l] = sorted(consdata[l], reverse=True)
            for c in consdata[l]:
                # if method == 'MCL' and c.identifier in donesubclass: continue
                command += c.toSQL(class_count)
                class_count += 1
                # if method == 'MCL': donesubclass.add(c.identifier)
                # subclass_count = 0
                # print '\tCONS ', c.consensus, 'SIZE ', c.size
                # for s in c.subclasses:
                #     print '\t\tSIZE ', s.size, ' SEED ', s._seed
       
        return command

    # def toSQL(self):
    #     for cla in sorted(self.subclasses):
    #         print 'LENGTH: ', cla
    #         for sub in self.subclasses[cla]:
    #             print 'SEED: ', sub._seed, ' SIZE: ', sub.size, ' S: ', len(sub.seq_pat.split('-')), ' R: ', len(sub.ram_pat.split('-')), ' E: ', len(sub.exp_pat.split('-'))
    #             print 'CONS: ', sub.consensus
    #             print sub.toSQL()
    #             for loop in sub.loops:
    #                 print loop.toSQL()

class Loop(object):
    def __init__(self, info):
        self.loopID   = info.split(':')[0].strip().split()[0]
        self.position = info.split(':')[0].strip().split()[3]
        self.length   = int(info.split(':')[0].strip().split()[4])
        self.seq      = info.split(':')[1].split('|')[0]
        self.exp      = None
        self.ram      = None
        self.scs      = None

        self.score    = float(info.split(':')[1].split('|')[1].strip())
        self.prel     = None
        self.order    = None
        self.refer    = None

    def loop_ID_pos(self): return tuple([self.loopID,self.position])
    def toSQL(self, order):
        return self.loop_nid_SQL() + self.loop2cluster(order)

    def loop_nid_SQL(self):
        command  = "SELECT lc.loop_id INTO @loopNID FROM loop2chain lc, chain c\n"
        command += " WHERE c.pdb='{0[0]}' and c.chain='{0[1]}'\n".format(self.loopID.split('_'))
        command += "  AND lc.chain=c.nid\n"
        command += "  AND lc.position={0};\n".format(self.position)
        return command

    def loop2cluster(self, order):
        command  = "INSERT INTO loop2cluster (cluster_nid,loop_nid,clust_order,score,seq_ali,ss_ali,exp_ali,ram_ali)\n"
        command += " VALUES (@subclassID,@loopNID,{1},{0.score:.3f},'{0.seq}','{0.scs}','{0.exp}','{0.ram}');\n".format(self, order)
        return command

    def fix_align(self, ini, end):
        self.seq = re.sub(' ','-',self.seq[int(ini):int(end)])
        self.exp = re.sub(' ','-',self.exp[int(ini):int(end)])
        self.ram = re.sub(' ','-',self.ram[int(ini):int(end)])
        self.scs = re.sub(' ','-',self.scs[int(ini):int(end)])

    # @property
    # def iniseq(self):
    #     return len(self.seq) - len(self.seq.lstrip())

    # @property
    # def endseq(self):
    #     return -(len(self.seq) - len(self.seq.rstrip()))

    def add_surface(self, info):

        self.exp   = info.split(':')[1].split('|')[0]
        self.refer = len(info.split(':')[0]) + 1 
        position   = re.search('\[\s*\d+\s*\+\s*(\-*\d+),\s*\d+\s*\+\s*(\-*\d+)\]', info.split(':')[1].split('|')[-1])
        p1,p2      = int(position.group(1)), int(position.group(2))
        # correct    = self.prel + 1 + len(self.exp) - len(self.exp.lstrip())
        self.prel = len(info.split(':')[0]) + 1 + len(info.split(':')[1].split('|')[0]) - len(info.split(':')[1].split('|')[0].lstrip())
        # self.refer.extend([p1+correct, p2+correct])

    def add_ramachandran(self, info):
        self.ram   = info.split(':')[1].split('|')[0]

    def add_secondary_str(self, info):
        self.scs   = info.split(':')[1].split('|')[0]

    

    # def fix_align(self, ini, end):
    #     self.seq = re.sub(' ','-',self.seq[int(ini):int(end)])
    #     self.exp = re.sub(' ','-',self.exp[int(ini):int(end)])
    #     self.ram = re.sub(' ','-',self.ram[int(ini):int(end)])
    #     self.scs = re.sub(' ','-',self.scs[int(ini):int(end)])

    # def __repr__(self):
    #     data = []

    #     data.append("SELECT nid FROM loop_description WHERE loop_id='{0.id}' INTO @loopNID;".format(self))
    #     data.append("INSERT INTO loop2cluster (cluster_nid,loop_nid,clust_order,score,seq_ali,ss_ali,exp_ali,ram_ali)" + \
    #                 " VALUES (LAST_INSERT_ID(),@loopNID,{0.order},{0.score:.3f},'{0.seq}','{0.scs}','{0.exp}','{0.ram}');".format(self))

    #     return "\n".join(data)