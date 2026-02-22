import re
import math
import gzip

def parse_pdb_atom(line, convert_resSeq=True, strip_fields=None,
                   modelSerial=None, keep_line=False):
    if line.startswith('ATOM  '):
        recName = line[:6].strip()
        try:
            #serial = int(line[6:11])
            # there are some files in PDBTM with atoms that have serial numbers
            # >99999 and everything is shifted, check for this
            i = 11
            while line[i].isdigit():
                i += 1
            serial = int(line[6:i])
            line = line[:11]+line[i:]
        except ValueError:
            serial = None
        name = line[12:16]
        altLoc = line[16]
        resName = line[17:20]
        if ((resName[0] == ' ') and (altLoc != ' ') and (name[-1] != ' ') and
            (name[0] == ' ')):
            # probably the line is shifted due to extra space, caused by
            # 4-letter atom names not aligned correctly, e.g. PDBTM 7u8g
            line = line[:12] + line[13:]
            name = line[12:16]
            altLoc = line[16]
            resName = line[17:20]
        chainID = line[21]
        if convert_resSeq:
            try:
                resSeq = int(line[22:26])
            except ValueError:
                resSeq = None
        else:
            resSeq = line[22:26]
        iCode = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        try:
            occupancy = float(line[54:60])
        except ValueError:
            occupancy = None
        try:
            tempFactor = float(line[60:66])
        except ValueError:
            tempFactor = None
        # save this due to pdbqt format
        padding = line[66:76]
        element = line[76:78]
        charge = line[78:80]
        d = dict(
            recName=recName,
            serial=serial,
            name=name,
            altLoc=altLoc,
            resName=resName,
            chainID=chainID,
            resSeq=resSeq,
            iCode=iCode,
            x=x,
            y=y,
            z=z,
            occupancy=occupancy,
            tempFactor=tempFactor,
            element=element,
            charge=charge,
            padding=padding,
        )
        if modelSerial is not None: d['modelSerial'] = modelSerial
        if keep_line: d['line'] = line
        if strip_fields:
            for field in strip_fields:
                d[field] = d[field].strip()
        return d
    else: return None

def parse_pdb_hetatm(line, convert_resSeq=False, strip_fields=None,
                     modelSerial=None, keep_line=False):
    if line.startswith('HETATM'):
        recName = line[:6].strip()
        try:
            #serial = int(line[6:11])
            # there are some files in PDBTM with atoms that have serial numbers
            # >99999 and everything is shifted, check for this
            i = 11
            while line[i].isdigit():
                i += 1
            serial = int(line[6:i])
            line = line[:11]+line[i:]
        except ValueError:
            serial = None
        name = line[12:16]
        altLoc = line[16]
        resName = line[17:20]
        if ((resName[0] == ' ') and (altLoc != ' ') and (name[-1] != ' ') and
            (name[0] == ' ')):
            # probably the line is shifted due to extra space, caused by
            # 4-letter atom names not aligned correctly, e.g. PDBTM 7u8g
            line = line[:12] + line[13:]
            name = line[12:16]
            altLoc = line[16]
            resName = line[17:20]
        chainID = line[21]
        if convert_resSeq:
            try:
                resSeq = int(line[22:26])
            except ValueError:
                resSeq = None
        else:
            resSeq = line[22:26]
        iCode = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        try:
            occupancy = float(line[54:60])
        except ValueError:
            occupancy = None
        try:
            tempFactor = float(line[60:66])
        except ValueError:
            tempFactor = None
        element = line[76:78]
        charge = line[78:80]
        d = dict(
            recName=recName,
            serial=serial,
            name=name,
            altLoc=altLoc,
            resName=resName,
            chainID=chainID,
            resSeq=resSeq,
            iCode=iCode,
            x=x,
            y=y,
            z=z,
            occupancy=occupancy,
            tempFactor=tempFactor,
            element=element,
            charge=charge,
        )
        if modelSerial is not None: d['modelSerial'] = modelSerial
        if keep_line: d['line'] = line
        if strip_fields:
            for field in strip_fields:
                d[field] = d[field].strip()
        return d
    else: return None

def parse_pdb_seqres(line):
    if line.startswith('SEQRES'):
        serNum = int(line[7:10])
        chainID = line[11]
        numRes = int(line[13:17])
        resNames = [
            line[19+x*4:19+3+x*4].strip()
            for x in range(13)
        ]
        resNames = [x for x in resNames if x != '']
        d = dict(
            serNum=serNum,
            chainID=chainID,
            numRes=numRes,
            resNames=resNames,
        )
        return d
    else: return None

def parse_pdb_atoms(fn,
                    first_model=False,
                    # currently this only applies to ATOM
                    allowed_altloc=' A',
                    # ATOM filter and options
                    atom_filter=None,
                    atom_opts=None,
                    # HETATM
                    hetatm=False,           # pass 'seqres' to parse non-std aas
                    hetatm_filter=None,
                    hetatm_opts=None,
                    alt_seqres=None,
                    **kwargs
                   ):
    """
    You can push kwargs to parse_pdb_atom and parse_pdb_hetatm either through
    **kwargs directly, in which case they will be applied to both atom and
    hetatm, or through atom_opts and hetatm_opts, which will be applied
    individually.

    The other **kwargs will be send to both of them.

    FIXME: There is an issue that modified amino acids will be dropped if heteroatoms
    are not wanted. One solution would be to keep heteroatoms if they are part
    of a chain that also has regular atoms.
    """
    atoms = []
    if fn.endswith('gz'):
        f = gzip.open(fn, 'rt', encoding='utf-8')
    else:
        f = open(fn, 'rt', encoding='utf-8')
    if atom_opts is None: atom_opts = dict()
    atom_opts.update(kwargs)
    if hetatm_opts is None: hetatm_opts = dict()
    hetatm_opts.update(kwargs)
    seq_resnames = dict()
    for line in f:
        if line.startswith('SEQRES'):
            # we should parse this for non-standard amino acids, etc.
            seqres = parse_pdb_seqres(line)
            ch = seqres['chainID']
            if ch in seq_resnames:
                seq_resnames[ch].update(seqres['resNames'])
            else:
                seq_resnames[ch] = set(seqres['resNames'])
        if line.startswith('ATOM  '):
            atom = parse_pdb_atom(line, **atom_opts)
            if atom['altLoc'] not in allowed_altloc: continue
            if (atom_filter is not None) and (not atom_filter(atom)): continue
            atoms.append(atom)
        elif line.startswith('HETATM'):
            if hetatm == False: continue
            atom = parse_pdb_hetatm(line, **hetatm_opts)
            if hetatm == 'seqres':
                resnames = seq_resnames.get(atom['chainID'], set())
                if alt_seqres is not None:
                    resnames.update(alt_seqres)
                if atom['resName'] not in resnames: continue
            elif hetatm != True: continue
            if atom['altLoc'] not in allowed_altloc: continue
            if (hetatm_filter is not None) and (not hetatm_filter(atom)): continue
            atoms.append(atom)
        elif line.startswith('MODEL '):
            modelSerial = int(line[10:14])
            atom_opts['modelSerial'] = modelSerial
            hetatm_opts['modelSerial'] = modelSerial
        elif line.startswith('ENDMDL'):
            if 'modelSerial' in atom_opts: del atom_opts['modelSerial']
            if 'modelSerial' in hetatm_opts: del hetatm_opts['modelSerial']
            if first_model: break
    f.close()
    return atoms

def separate_chains(atoms,
                    atom_filter=None,
                    rename_duplicate_chains=False,
                    merge_same_chains=False):
    chains = []
    chain = []
    last_atom = None
    # find segments with same chainID
    for atom in atoms:
        # skip filtered atoms
        if (atom_filter is not None) and (not atom_filter(atom)): continue
        if ((last_atom is not None) and
            ((atom['chainID'] != last_atom['chainID']) or
             (atom.get('modelSerial', None) != 
              last_atom.get('modelSerial', None)))
           ):
            # chain break or next model
            if len(chain) > 0:
                chains.append(chain)
                chain = []
        chain.append(atom)
        last_atom = atom
    # add last one
    if len(chain) > 0:
        chains.append(chain)
        chain = []
    # rename segments with duplicate chainID
    if rename_duplicate_chains:
        seen_chain_ids = set()
        for chain in chains:
            chain_id = chain[0]['chainID']
            i = 2
            while chain_id in seen_chain_ids:
                chain_id = chain[0]['chainID'] + str(i)
                i += 1
            if i > 2:
                # rename chainID in all atoms
                for atom in chain:
                    atom['chainID'] = chain_id
            seen_chain_ids.add(chain_id)
    # merge segments with same chainID
    if merge_same_chains:
        # python >=3.7 keeps directory insertion order
        all_chains = dict()
        for chain in chains:
            chain_id = chain[0]['chainID']
            if chain_id not in all_chains:
                all_chains[chain_id] = chain
            else:
                # merge with existing
                all_chains[chain_id].extend(chain)
        chains = list(all_chains.values())
    return chains

def parse_mmcif_record(line, headers, convert=dict(), rename=dict()):
    l = line.strip().split()
    d = dict()
    for value, key in zip(l, headers):
        d[key] = value
    for key, filt in convert.items():
        if key in d:
            # I think this is a special placeholder for an undefined value
            if d[key] == '.':
                d[key] = None
            else:
                d[key] = filt(d[key])
    for new_name, old_name in rename.items():
        if new_name in headers:
            d[old_name] = d[new_name]
    return d

def parse_mmcif_block(f, block_name, convert=dict(), rename=dict(), record_filter=None):
    lead = block_name + '.'
    headers = []
    for line in f:
        if line.startswith(lead):
            headers.append('.'.join(line.strip().split('.')[1:]))
        elif headers:
            if line.startswith('#'): break
            d = parse_mmcif_record(line, headers, convert=convert,
                                   rename=rename)
            if record_filter and not record_filter(d):
                continue
            yield d

def parse_mmcif_chains(fn, allowed_altloc='.A', hetatm=False, convert=dict(),
                       rename=dict(), chain_field='chainID', only_model=None):
    if fn.endswith('gz'):
        f = gzip.open(fn, 'rt', encoding='utf-8')
    else:
        f = open(fn, 'rt', encoding='utf-8')
    atom_opts = dict(
        # compatibility with pdb
        convert={'Cartn_x': float, 'Cartn_y': float, 'Cartn_z': float,
                 'auth_seq_id': int},
        rename={'auth_asym_id': 'chainID',
                'Cartn_x': 'x',
                'Cartn_y': 'y',
                'Cartn_z': 'z',
                'auth_comp_id': 'resName',
                'auth_seq_id': 'resSeq',
                'auth_atom_id': 'name',
                'group_PDB': 'recName'},
    )
    atom_opts['convert'].update(convert)
    atom_opts['rename'].update(rename)

    seq_resnames = dict()
    for record in parse_mmcif_block(f, '_pdbx_poly_seq_scheme'):
        chain_id = record['asym_id']
        if chain_id not in seq_resnames:
            seq_resnames[chain_id] = set([record['mon_id']])
        else:
            seq_resnames[chain_id].update([record['mon_id']])

    # the blocks can be in any order, sadly... e.g. 1amt.cif.gz in PDB
    f.seek(0)

    pdb1_chains = []
    chain = []
    last_atom = None
    for atom in parse_mmcif_block(
        f, '_atom_site', **atom_opts
    ):
        if atom['recName'] == 'HETATM':
            if hetatm == False: continue
            elif hetatm == 'seqres':
                if ((atom['label_asym_id'] not in seq_resnames) or 
                    (atom['resName'] not in seq_resnames[atom['label_asym_id']])): continue
            elif hetatm != True: continue
        # skip hydrogens for some reason
        if atom['label_atom_id'].startswith('H'): continue
        if (allowed_altloc is not None) and (atom['label_alt_id'] not in allowed_altloc): continue
        if (only_model is not None) and (int(atom['pdbx_PDB_model_num']) != only_model): continue
        chain_id = atom[chain_field]
        if ((last_atom is not None) and 
            (chain_id != last_atom[chain_field])):
            # chain break
            pdb1_chains.append(chain)
            chain = []
        chain.append(atom)
        last_atom = atom
    if len(chain) > 0:
        pdb1_chains.append(chain)
    f.close()
    return pdb1_chains

class Atom:
	def parsepdb(self, line):
		self.pqrline = False

		self.pdbrecname = line[:6]
		# lattam mar olyat, hogy vmd hexa serialt ad az atomoknak, ezert marad string
		self.pdbserial = line[6:11].strip()

		pdbname = line[12:16]
		name = pdbname.strip().upper()
		if (re.match("[A-Z]{2}\d\d", name)):
			name = name[3] + name[:3]
		self.name = name
		self.pdbname = pdbname

		self.pdbaltloc = line[16]

		#self.pdbresname = line[17:20]
		self.pdbresname = line[17:21]
		self.resname = self.pdbresname.strip()
		# hidrogeneket kiszorjuk
		#if (name[1] == "H"):
		#	continue

		self.pdbresid = line[22:27].lstrip()
		self.pdbresidnum = int(self.pdbresid[:-1])
		self.resid = -1

		self.chainid = line[21].strip()
		# Atom.makename miatt inkabb lower()
		self.coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]

		try:
			self.pdboccupancy = line[54:60]
			self.occupancy = float(self.pdboccupancy)
		except IndexError:
			self.pdboccupancy = None
			self.occupancy = 0.0
		except ValueError:
			self.occupancy = 0.0

		try:
			self.pdbtempfactor = line[60:66]
			self.tempfactor = float(self.pdbtempfactor)
		except IndexError:
			self.pdbtempfactor = None
			self.tempfactor = 0.0
		except ValueError:
			self.tempfactor = 0.0

		# save due to pdbqt format
		self.padding = line[66:76]

		self.pdbelement = line[76:78]
		self.element = self.pdbelement.strip()
		if (self.element == ""): self.pdbelement = ""

		self.pdbcharge = line[78:80]
		self.charge = self.pdbcharge.strip()
		if (self.charge == ""): self.pdbcharge = ""
		
		self.pdbline = line

		self.epilogue = line[80:]
		### ez miert???
		#if (not self.epilogue and (line[-1] != '\n')): self.epilogue = '\n'
		if (line[-1] == '\n'):
			if (not self.epilogue): self.epilogue = '\n'
			elif (self.epilogue[-1] != '\n'): self.epilogue += '\n'
		#if (not self.epilogue and (line[-1] == '\n')): self.epilogue = '\n'

	def parsepqr(self, line):
		self.pqrline = True

		self.pdbrecname = line[:6]
		# lattam mar olyat, hogy vmd hexa serialt ad az atomoknak, ezert marad string
		self.pdbserial = line[6:11].strip()

		pdbname = line[12:16]
		name = pdbname.strip().upper()
		if (re.match("[A-Z]{2}\d\d", name)):
			name = name[3] + name[:3]
		self.name = name
		self.pdbname = pdbname

		self.pdbaltloc = line[16]

		self.pdbresname = line[17:20]
		self.resname = self.pdbresname.strip()
		# hidrogeneket kiszorjuk
		#if (name[1] == "H"):
		#	continue

		self.pdbresid = line[22:27].lstrip()
		self.pdbresidnum = int(self.pdbresid[:-1])
		self.resid = -1

		self.chainid = line[21].strip()
		# Atom.makename miatt inkabb lower()
		self.coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]

		self.pdbcharge = line[54:62]
		self.charge = float(self.pdbcharge)

		self.pdbaccessibility = line[62:69]
		self.accessibility = float(self.pdbaccessibility)

		self.pdbline = line

		self.epilogue = line[80:]
		### ez miert???
		#if (not self.epilogue and (line[-1] != '\n')): self.epilogue = '\n'
		if (line[-1] == '\n'):
			if (not self.epilogue): self.epilogue = '\n'
			elif (self.epilogue[-1] != '\n'): self.epilogue += '\n'
		#if (not self.epilogue and (line[-1] == '\n')): self.epilogue = '\n'

	def __setattr__(self, name, value):
		self.__dict__[name] = value
		if (name == "name"):
			pdbname = self.__dict__.get("pdbname", "")
			if (len(value) == 4):
				if re.match("[A-Z]{2}\d\d", pdbname):
					if re.match("[A-Z]{2}\d\d", value): pdbname = value
					else: pdbname = value[1:] + value[1]
				else:
					if re.match("[A-Z]{2}\d\d", value): pdbname = value[3] + value[:3]
					else: pdbname = value
			else:
				if (len(self.__dict__.get("name", "")) == len(value)):
					m = re.match(r"^(\s+)(\S+)(\s+)$", pdbname)
					if (m):
						pdbname = m.group(1) + value + m.group(3)
					else:
						if (re.match(r"^\d[A-Z]", value)):
							pdbname = "%-4s" % value
						else:
							pdbname = "%-3s" % value
				else:
					if (re.match(r"^\d[A-Z]", value)):
						pdbname = "%-4s" % value
					else:
						pdbname = "%-3s" % value
			self.__dict__["pdbname"] = pdbname
		elif (name == "resname"):
			pdbresname = self.__dict__.get("pdbresname", "")
			if (len(self.__dict__.get("resname", "")) == len(value)):
				m = re.match(r"^(\s+)(\S+)(\s+)$", pdbresname)
				if (m):
					pdbresname = m.group(1) + value + m.group(3)
				else:
					pdbresname = "%-3s" % value
			else:
				pdbresname = "%-3s" % value
			self.__dict__["pdbresname"] = pdbresname
	def frompdb(self, name, chainid, origresid, resid, resname, coords):
		self.pqrline = False
		self.pdbname = name
		self.pdborigresid = origresid
		self.pdbresid = resid
		self.resid = -1
		self.resname = resname
		self.chainid = chainid
		self.coords = coords
	def makename(self):
		return "_".join([self.pdbname, self.pdbresid, self.chainid])
	"""
	def __str__(self):
		if (self.chainid): c = self.chainid
		else: c = ""
		if (self.pdbresid): r = self.pdbresid
		else: r = ""
		if (self.pdbname): n = self.pdbname
		else: n = ""
		return "(%s:%s:%s)" % (c, r, n)
	"""
	def __str__(self):
		if (self.pqrline): return self.formatpqr()
		else: return self.formatpdb()
	def formatpdb(self):
		if (type(self.pdbserial) == type(1)): pdbserial = self.pdbserial % 100000
		else: pdbserial = self.pdbserial[-5:]
		try:
			pdbresid = int(self.pdbresid[:-1]) % 10000
			pdbresid = "%d%s" % (pdbresid, self.pdbresid[-1])
		except ValueError:
			pdbresid = self.pdbresid[-5:]
		except TypeError:
			# self.pdbresid is int
			pdbresid = "%d " % (int(self.pdbresid), )
		s = "%s%5s %4s%s%-4s%1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f" % (self.pdbrecname, pdbserial, self.pdbname, self.pdbaltloc, self.pdbresname, self.chainid, pdbresid, self.coords[0], self.coords[1], self.coords[2], self.occupancy, self.tempfactor)
		### %-2s or %2s ??
		if (self.padding): s += self.padding
		else: s += "          "
		if (self.pdbelement): s += "%2s" % self.element
		if (self.pdbcharge): s += "%2s" % self.pdbcharge
		if (self.epilogue): s += self.epilogue
		return s
	def formatpqr(self):
		if (type(self.pdbserial) == type(1)): pdbserial = self.pdbserial % 100000
		else: pdbserial = self.pdbserial
		try:
			pdbresid = int(self.pdbresid[:-1]) % 10000
			pdbresid = "%d%s" % (pdbresid, self.pdbresid[-1])
		except ValueError:
			pdbresid = self.pdbresid[-5:]
		except TypeError:
			# self.pdbresid is int
			pdbresid = "%d " % (int(self.pdbresid), )
		s = "%s%5s %4s%s%3s %1s%5s   %8.3f%8.3f%8.3f%8.4f%7.4f" % (self.pdbrecname, pdbserial, self.pdbname, self.pdbaltloc, self.pdbresname, self.chainid, pdbresid, self.coords[0], self.coords[1], self.coords[2], self.charge, self.accessibility)
		### %-2s or %2s ??
		"""
		if (self.pdbelement): s += "          %2s" % self.element
		if (self.pdbcharge): s += "%2s" % self.pdbcharge
		"""
		if (self.epilogue): s += self.epilogue
		return s
	def dist(self, b):
		return math.sqrt(sum([(self.coords[x]-b.coords[x])**2 for x in range(3)]))

class Residue:
	def __init__(self, name, pdbid, chainid):
		self.name = name
		self.id = -1
		self.pdbid = pdbid
		self.pdbidnum = int(pdbid[:-1])
		self.chainid = chainid
		self.atoms = []
		self.atombypdbname = {}
	def addatom(self, a):
		self.atoms.append(a)
		self.atombypdbname[a.pdbname] = a
		a.res = self
	def delatom(self, a):
		self.atoms.remove(a)
		del self.atombypdbname[a.pdbname]
	def __lt__(self, other):
		try: a = self.sortkey
		except AttributeError: a = None
		try: b = other.sortkey
		except AttributeError: b = None
		return a < b
    # this does not work anymore in python 3.x
	def __cmp__(self, other):
		try: a = self.sortkey
		except AttributeError: a = None
		try: b = other.sortkey
		except AttributeError: b = None
		return cmp(a, b)
	def __str__(self):
		if (self.chainid): c = self.chainid
		else: c = ""
		if (self.pdbid): i = self.pdbid.strip().lower()
		else: i = ""
		if (self.name): n = self.name
		else: n = ""
		return "(%s:%s %s)" % (c, i, n)
	def __setattr__(self, name, value):
		self.__dict__[name] = value
		if (name == "pdbid"):
			self.__dict__["sortkey"] = "%5s" % value
			if (type(value) == type(1)): self.__dict__[name] = "%s " % value
			atoms = self.__dict__.get("atoms")
			if (atoms):
				for a in atoms: a.pdbresid = self.pdbid
		if (name == "chainid"):
			atoms = self.__dict__.get("atoms")
			if (atoms):
				for a in atoms: a.chainid = self.chainid
		if (name == "name"):
			atoms = self.__dict__.get("atoms")
			if (atoms):
				for a in atoms: a.resname = self.name
	def formatpdb(self):
		l = []
		for a in self.atoms:
			l.append(a.formatpdb())
		return "".join(l)
	def formatpqr(self):
		l = []
		for a in self.atoms:
			l.append(a.formatpqr())
		return "".join(l)

class Chain:
	def __init__(self, id):
		self.id = id
		self.residues = []
		self.residuebypdbid = {}
		self.residuebyid = {}
		self.sortedresidues = []
		self.atoms = []
		self.atomsbypdbname = {}
		self.residuesbyname = {}
	def addatom(self, a):
		newres = False
		r = self.residuebypdbid.get(a.pdbresid)
		if (r == None):
			newres = True
			r = Residue(a.resname, a.pdbresid, a.chainid)
			self.residues.append(r)
			self.residuebypdbid[a.pdbresid] = r
			l = self.residuesbyname.get(r.name)
			if (l == None): self.residuesbyname[r.name] = [r]
			else: l.append(r)
			r.chain = self
		r.addatom(a)
		self.atoms.append(a)
		a.chain = self

		l = self.atomsbypdbname.get(a.name)
		if (l == None): self.atomsbypdbname[a.name] = [a]
		else: self.atomsbypdbname[a.name].append(a)
		
		return newres
	def delatom(self, a):
		self.residuebypdbid[a.pdbresid].delatom(a)
		self.atoms.remove(a)
	def genserials(self, reslist):
		l = self.residues[:]
		l.sort()
		for r in l:
			r.id = len(reslist)
			for a in r.atoms:
				a.resid = r.id
			reslist.append(r)
			self.residuebyid[r.id] = r
		self.sortedresidues = l

class PDB:
	def __init__(self):
		self.chains = []
		self.chainbyid = {}
		self.sortedchains = []
		self.residuebyid = []
		self.resnamebyid = []
		self.atoms = []
		self.atombyname = {}
		self.atomsbypdbname = {}
		self.atombypdbserial = {}
		self.residuesbyname = {}
		self.bonds = []
	def addatom(self, a):
		c = self.chainbyid.get(a.chainid)
		if (c == None):
			c = Chain(a.chainid)
			self.chains.append(c)
			self.chainbyid[a.chainid] = c
		if (c.addatom(a)):
			# new residue
			l = self.residuesbyname.get(a.resname)
			if (l == None): self.residuesbyname[a.resname] = [a.res]
			else: l.append(a.res)
		self.atoms.append(a)
		self.atombypdbserial[a.pdbserial] = a
		self.atombyname[a.makename()] = a

		l = self.atomsbypdbname.get(a.name)
		if (l == None): self.atomsbypdbname[a.name] = [a]
		else: self.atomsbypdbname[a.name].append(a)
	def delatom(self, a):
		a.chain.delatom(a)
		#self.chainbyid[a.chainid].delatom(a)
		self.atoms.remove(a)
		del self.atombyname[a.makename()]
		del self.atomsbypdbname[a.name]
		# ellenorizni kellene, hogy meguresedett-e a residue, mert ha igen, azt is torolni kell
	def addbond(self, b):
		for s in b.pdbbound:
			a = self.atombypdbserial[s]
			b.append(a)
		a0 = self.atombypdbserial[b.pdbserial]
		if (not hasattr(a0, "bound")): a0.bound = []
		for a in b:
			a0.bound.append(a)
		b.atom = a0

		self.bonds.append(b)
	def genserials(self):
		l = list(self.chainbyid.values())
		#l.sort(lambda x, y: cmp(x.id, y.id))
		l.sort(key=lambda a: a.id)
		self.residuebyid = []
		for c in l:
			c.genserials(self.residuebyid)
		self.resnamebyid = []
		for r in self.residuebyid:
			self.resnamebyid.append(r.name)
		self.sortedchains = l

class Cryst:
	def parsepdb(self, line):
		self.pdbrecname = line[:6]
		try:
			self.a = float(line[6:15])
			self.b = float(line[15:24])
			self.c = float(line[24:33])
			self.alpha = float(line[33:40])
			self.beta = float(line[40:47])
			self.gamma = float(line[47:54])
			self.sgroup = line[55:66]
			self.z = line[66:70]
			self.epilogue = line[70:]
			if (not self.epilogue): self.epilogue = '\n'
			self.gromacs = False
		except ValueError:
			# gromacs hibas pdb
			self.a = float(line[6:15])
			self.b = float(line[15:24])
			self.c = float(line[24:33])
			self.epilogue = line[33:]
			if (not self.epilogue): self.epilogue = '\n'
			self.gromacs = True
	def formatpdb(self):
		if (self.gromacs):
			return "%s%9.3f%9.3f%9.3f%s" % (self.pdbrecname, self.a, self.b, self.c, self.epilogue)
		else:
			return "%s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s%s%s" % (self.pdbrecname, self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.sgroup, self.z, self.epilogue)
	def __str__(self):
		return self.formatpdb()

class Bond(list):
	def parsepdb(self, line):
		self.pdbserial = line[6:11].strip()
		self.pdbbound = [x.strip() for x in [line[11:16], line[16:21], line[21:26], line[26:31]] if (x.strip() != "")]
		self.epilogue = line[31:]
		if (not self.epilogue): self.epilogue = '\n'
	def formatpdb(self):
		l = [x.pdbserial for x in self]
		l = l + ["" for x in range((4 - (len(l) % 4)) % 4)]
		lines = []
		while (len(l) > 0):
			s = "CONECT%5s%5s%5s%5s%5s" % (self.pdbserial, l[0], l[1], l[2], l[3])
			if (self.epilogue): s += self.epilogue
			lines.append(s)
			del l[0:4]
		return "".join(lines)
	def __str__(self):
		return self.formatpdb()

class PDBFile:
	def __init__(self):
		self.struc = PDB()
		self.cryst = None
		self.pdblines = []
		self.prologue = []
		self.epilogue = []
		self.body = []
	# lines can be file or any iterable
	def parsepdb(self, lines):
		pro = True
		for s in lines:
			if (s.startswith("ATOM  ") or s.startswith("HETATM")):
				pro = False
				a = Atom()
				a.parsepdb(s)
				self.struc.addatom(a)
				self.pdblines.append(a)
				self.body.append(a)
			else:
				if (s.startswith("CRYST1")):
					c = Cryst()
					c.parsepdb(s)
					s = c
					self.cryst = c
				elif (s.startswith("CONECT")):
					b = Bond()
					b.parsepdb(s)
					s = b
					self.struc.addbond(b)
				elif (s.startswith("ENDMDL")):
					break
				if (pro): self.prologue.append(s)
				else: self.epilogue.append(s)
				self.pdblines.append(s)
		self.struc.genserials()

	def formatpdb(self):
		for s in self.pdblines:
			print(str(s), end='')

	def formatpdbstr(self):
		result = []
		for s in self.pdblines:
			result.append(str(s))
		return "".join(result)

	def parsepqr(self, lines):
		pro = True
		for s in lines:
			if (s.startswith("ATOM  ") or s.startswith("HETATM")):
				pro = False
				a = Atom()
				a.parsepqr(s)
				self.struc.addatom(a)
				self.pdblines.append(a)
				self.body.append(a)
			else:
				if (s.startswith("CRYST1")):
					c = Cryst()
					c.parsepdb(s)
					s = c
					self.cryst = c
				elif (s.startswith("CONECT")):
					b = Bond()
					b.parsepdb(s)
					s = b
					self.struc.addbond(b)
				elif (s.startswith("ENDMDL")):
					break
				if (pro): self.prologue.append(s)
				else: self.epilogue.append(s)
				self.pdblines.append(s)
		self.struc.genserials()
	def formatpqr(self):
		for s in self.pdblines:
			print(str(s), end='')
