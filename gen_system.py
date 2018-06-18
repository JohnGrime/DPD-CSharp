import sys, math, random

def GenLipid( head_sites, tail_sites, N_tails ):
	"""
	Generate an example lipid (site names, bonds, angles) given a
	list of head names and tail names, and the number of tails (1 or 2)

	Args:
	  head_sites (list of strings): site names for the head region
	  tail_sites (list of strings): site names for the tail region
	  N_tails (int) : 1 or 2, specify single- or dual-tail lipid

	Returns:
		3-tuple of (sites,bonds,angles), each tuple element is a list
	"""

	b_eq, b_K = 0.7, 100.0
	a_eq, a_K = 3.142, 20.0

	all_sites  = []
	all_bonds  = []
	all_angles = []

	offset = 0

	#
	# Head region
	#

	sites = head_sites
	N = len(sites)
	s,b,a = [],[],[]

	for si in range(0,N):
		s.append( sites[si] )

	for bi in range(0,N-1):
		i = offset + bi + 1
		b.append( [i, i+1, b_eq, b_K] )

	for ai in range(0,N-2):
		i = offset + ai + 1
		a.append( [i ,i+1, i+2, a_eq, a_K] )

	if( len(s) > 0 ): all_sites.append( s )
	if( len(b) > 0 ): all_bonds.append( b )
	if( len(a) > 0 ): all_angles.append( a )

	offset += len(sites)

	#
	# Tail 1
	#

	sites = tail_sites
	N = len(sites)
	s,b,a = [],[],[]

	for si in range(0,N):
		s.append( sites[si] )

	for bi in range(0,N-1):
		i = offset + bi + 1
		b.append( [i, i+1, b_eq, b_K] )

	for ai in range(0,N-2):
		i = offset + ai + 1
		a.append( [i ,i+1, i+2, a_eq, a_K] )

	if( len(s) > 0 ): all_sites.append( s )
	if( len(b) > 0 ): all_bonds.append( b )
	if( len(a) > 0 ): all_angles.append( a )

	offset += len(sites)

	#
	# Tail 2?
	#

	if N_tails == 2:
		sites = tail_sites
		N = len(sites)
		s,b,a = [],[],[]

		for si in range(0,N):
			s.append( sites[si] )

		for bi in range(0,N-1):
			i = offset + bi + 1
			b.append( [i, i+1, b_eq, b_K] )

		for ai in range(0,N-2):
			i = offset + ai + 1
			a.append( [i ,i+1, i+2, a_eq, a_K] )

		if( len(s) > 0 ): all_sites.append( s )
		if( len(b) > 0 ): all_bonds.append( b )
		if( len(a) > 0 ): all_angles.append( a )

	#
	# Interconnect head to tail(s), additional angle if needed
	#

	h_N  = len(head_sites)      # last head site
	t_1 = h_N + 1               # first site of tail 1
	t_2 = t_1 + len(tail_sites) # first site of tail 2
	s,b,a = [],[],[]

	b.append( [h_N, t_1, b_eq, b_K] )

	if( N_tails == 2 ):
		b.append( [h_N, t_2, b_eq, b_K] )
		a.append( [t_1,h_N,t_2, 1.571, 3.0] )

	if( len(s) > 0 ): all_sites.append( s )
	if( len(b) > 0 ): all_bonds.append( b )
	if( len(a) > 0 ): all_angles.append( a )

	return all_sites, all_bonds, all_angles

def PrintSettings( f, settings ):
	"""
	Prints a DPD-style 'settings' section to file.

	Args:
	  f (file): destination for output
	  settings (dictionary: string->string) : settings dictionary
	"""
	print >>f, 'settings'
	for key in settings:
		if key == 'cell':
			x,y,z = settings[key]
			print >>f, '\t', 'cell', str(x), str(y), str(z)
		else:
			print >>f, '\t', key, str(settings[key])
	print >>f, 'end'

def PrintSites( f, site_names ):
	"""
	Prints a DPD-style 'sites' section to file.

	Args:
	  f (file): destination for output
	  site_names (list of strings) : list of site names
	"""
	print >>f, 'sites'
	for name in site_names:
		print >>f, '\t', name
	print >>f, 'end'

def PrintInteractions( f, interactions ):
	"""
	Prints a DPD-style 'intractions' section to file.

	Args:
	  f (file): destination for output
	  interactions (2d list of float) : indexed by site type; interactions[i][j] == pair coefficient for i,j
	"""

	print >>f, 'interactions'
	for a in interactions:
		for b in interactions[a]:
			print >>f, '\t %s %s %.3f' % ( a, b, interactions[a][b] )
	print >>f, 'end'

def PrintMolecule( f, name, count, sites, bonds, angles ):
	"""
	Prints a DPD-style 'molecule' section to file.

	Args:
	  f (file): destination for output
	  count (int) : number of instances of this molecule type
	  sites (list of strings) : member site names for the molecule
	  bonds (list of [i,j,eq,K]) : bonds in the molecule
	  angles (list of [i,j,k,eq,K]) : angles in the molecule
	"""

	bnd_fmt = '\t bond  %5d %5d %.3f %.3f'
	ang_fmt = '\t angle %5d %5d %5d %.3f %.3f'

	print >>f, 'molecule'
	print >>f, ''
	print >>f, '\t name %s' % ( name )
	print >>f, ''
	print >>f, '\t count %s' % ( count )
	print >>f, ''
	idx = 1
	for struct in sites:
		for s in struct:
			print >>f, '\t site %s \t# %d' % ( s, idx )
			idx += 1
		print >>f, ''
	for bond in bonds:
		for b in bond:
			i, j, eq, K = b
			print >>f, bnd_fmt % ( i, j, eq, K )
		print >>f, ''
	for angle in angles:
		for a in angle:
			i, j, k, eq, K = a
			print >>f, ang_fmt % ( i, j, k, eq, K )
		print >>f, ''
	print >>f, 'end'

def PrintCoords( f, coords ):
	"""
	Prints a DPD-style 'coords' section to file.

	Args:
	  f (file): destination for output
	  coords (list of [name,x,y,z] entries ) : coordinate data to write
	"""
	fmt = "\t %s %f %f %f"
	print >>f, 'coords'
	for c in coords:
		name, x, y, z = c
		print >>f, fmt % ( name, x, y, z )
	print >>f, 'end'

def MakePDBAtomLine( serial, name, resName, chainID, resSeq, x, y, z ):
	"""
	Generates a line of text suitable for writing to a PDB file from a PDB atom dictionary.

	Args:
	  atom (dictionary): A PDB atom dictionary

	Returns:
	  A string suitable for writing as an ATOM or HETATM line in a PDB file.
	"""
	PDB_atom_line_format = '%-6.6s%5.5s %4.4s%1.1s%3.3s %1.1s%4d%1.1s   %8.2f%8.2f%8.2f%6.6s%6.6s          %2.2s%2.2s'

	string = PDB_atom_line_format % (
		'ATOM',
		str(serial),
		name,
		'',
		resName,
		chainID,
		resSeq,
		'',
		x, y, z,
		'', '',
		'', '' )
	
	return string

def PrintPDBMolecule( f, atom_upto, mol_upto, resName, atom_sets, coords ):
	"""
	Prints a PDB 'molecule' to file, where each set of atoms in the structure is
	assigned a different chainID.

	Args:
	  f (file): destination for output
	  atom_upto (int) : serial number of first atom in the molecule
	  mol_upto (int) : current molecule number
	  resName (string) : resName to use in output as 'name' of molecule
	  atom_sets (list of string lists) : each set is considered a different chain
	  coords (list of [name,x,y,z] entries ) : coordinate data, indexed using atom_upto + offset

	Returns:
	  Updated atom_upto and mol_upto values after writing the molecule data to file
	"""
	chains = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
	for si in range( 0, len(atom_sets) ):
		for ai in range( 0, len(atom_sets[si]) ):
			serial = atom_upto+ai
			name = atom_sets[si][ai]
			chainID = chains[si]
			resSeq = mol_upto
			x,y,z = coords[(atom_upto+ai)-1][1:4] # -1 as converting unit based to zero-based indices
			line = MakePDBAtomLine( serial, name, resName, chainID, resSeq, x, y, z )
			print >>f, line
		atom_upto += len(atom_sets[si])
	mol_upto += 1
	print >>f, 'TER'

	return atom_upto, mol_upto

def PrintPDBMoleculeBonds( f, offset, bond_sets ):
	"""
	Prints the CONECT entries for a PDB 'molecule' to file.

	Args:
	  f (file): destination for output
	  bonds (list of [i,j,eq,K]) : bonds in the molecule
	"""
	fmt = '%5d%5d'

	all_bonds = []
	for bond_set in bond_sets:
		for b in bond_set:
			i,j,eq,k = b
			bi = min( (offset+i)-1, (offset+j)-1 )
			bj = max( (offset+i)-1, (offset+j)-1 )
			all_bonds.append( [bi,bj] )

	sorted_bonds = sorted( all_bonds )
	base = None
	output = None
	counter = 0
	for bi in range( 0, len(sorted_bonds) ):
		b = sorted_bonds[bi]
		counter += 1

		if( (base!=b[0]) or (counter>4) ):
			if output != None: print >>f, output
			base = b[0]
			counter = 0
			output = 'CONECT%5d' % ( base )

		output += '%5d' % ( b[1] )

	print >>f, output

#
# Main script starts here.
#

#
# Set up default system information
#

settings = {
	'step_no': 1,
	'max_steps': 25000,
	'save_every': 1000,
	'print_every': 1000,
	'delta_t': 0.01,
	'lambda': 0.65,
	'sigma': 3.0,
	'ran1': -1,
	'kBT': 1.0,
	'cell': [10,10,10]
}

site_names = [ "w", "h", "t" ]

interactions = {}
for s in site_names:
	interactions[s] = {}

for i in range(0,len(site_names)):
	for j in range(i,len(site_names)):
		a, b = site_names[i], site_names[j]
		interactions[a][b] = 25.0

rho = 3.0 # target site density
phi = 0.3 # mass fraction of system that is lipid
head_len, tail_len, N_tails = 3, 5, 2

# Interaction weights after Venturoli et al

interactions['w']['w'] = 25.0
interactions['w']['h'] = 15.0
interactions['w']['t'] = 80.0

interactions['h']['h'] = 35.0
interactions['h']['t'] = 80.0

interactions['t']['t'] = 25.0

#
# Get user params from command line, and basic sanity check
#

if len(sys.argv) < 5:
	print ''
	print 'Usage: %s Lx Ly Lz phi' % ( sys.argv[0] )
	print ''
	print 'Where:'
	print '  - Lx, Ly, Lz : simulation dimensions'
	print '  - phi: lipid mass fraction of system'
	print ''
	sys.exit( -1 )

Lx = float(sys.argv[1])
Ly = float(sys.argv[2])
Lz = float(sys.argv[3])
phi = float(sys.argv[4])

if ( (Lx<=0.0) or Lx<=0.0 or Lx<=0.0 ):
	print >>sys.stderr, 'Bad system dimensions', Lx, Ly, Lz
	sys.exit( -1 )

if( (phi<0.0) or (phi>=1.0) ):
	print >>sys.stderr, 'Bad phi', phi
	sys.exit( -1 )

settings['cell'] = [Lx,Ly,Lz]

#
# Generate DPD system & write inut file
#

Lx, Ly, Lz = settings["cell"]

head_sites = [ 'h' for i in range(0,head_len) ]
tail_sites = [ 't' for i in range(0,tail_len) ]

lipid_sites, lipid_bonds, lipid_angles = GenLipid( head_sites, tail_sites, N_tails )

lipid_length = head_len + (N_tails*tail_len)
N_sites_total = int( Lx*Ly*Lz*rho )
N_lipid_molecules = int( (phi*N_sites_total)/lipid_length )
N_water_molecules = N_sites_total - (N_lipid_molecules*lipid_length)

coords = []

# All coords water, to start with
for i in range( 0, N_sites_total ):
	x = random.random() - 0.5
	y = random.random() - 0.5
	z = random.random() - 0.5
	coords.append( ['w', x, y, z] )

# Swap some waters to lipids
flat_names = []
for domain in lipid_sites:
	for s in domain:
		flat_names.append( s );
upto = 0
for li in range( 0, N_lipid_molecules ):
	for i in range( 0, lipid_length ):
		coords[upto][0] = flat_names[i]
		upto += 1


f = sys.stdout

print >>f, '#'
print >>f, '# %s' % ( ' '.join( sys.argv ) )
print >>f, '#'
print >>f, '# %d total sites, %d lipids of length %d (phi = %.3f)' % ( N_sites_total, N_lipid_molecules, lipid_length, float(N_lipid_molecules*lipid_length)/N_sites_total )
print >>f, '#'

print >>f, ''
PrintSettings( f, settings )
print >>f, ''
PrintSites( f, site_names )
print >>f, ''
PrintInteractions( f, interactions )
print >>f, ''
PrintMolecule( f, "lipid", N_lipid_molecules, lipid_sites, lipid_bonds, lipid_angles )
print >>f, ''
PrintMolecule( f, "water", N_sites_total-(N_lipid_molecules*lipid_length), ['w'], [], [] )
print >>f, ''
PrintCoords( f, coords )
print >>f, ''

#
# Write PDB file for more straightforward visualization
#

f = open( 'system.pdb', 'w' )

atom_upto = 1
mol_upto = 1

for mi in range( 0, N_lipid_molecules ):
	atom_upto, mol_upto = PrintPDBMolecule( f, atom_upto, mol_upto, 'LPD', lipid_sites, coords )

for mi in range( 0, N_water_molecules ):
	atom_upto, mol_upto = PrintPDBMolecule( f, atom_upto, mol_upto, 'SOL', [ ['w'] ], coords )

atom_upto = 1

for mi in range( 0, N_lipid_molecules ):
	PrintPDBMoleculeBonds( f, atom_upto, lipid_bonds )
	atom_upto += lipid_length

f.close()
