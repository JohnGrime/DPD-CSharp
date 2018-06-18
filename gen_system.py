import sys, math, random

def GenLipid( head_sites, tail_sites, N_tails ):

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
	print >>f, 'settings'
	for key in settings:
		if key == 'cell':
			x,y,z = settings[key]
			print >>f, '\t', 'cell', str(x), str(y), str(z)
		else:
			print >>f, '\t', key, str(settings[key])
	print >>f, 'end'

def PrintSites( f, site_names ):
	print >>f, 'sites'
	for name in site_names:
		print >>f, '\t', name
	print >>f, 'end'

def PrintInteractions( f, interactions ):
	print >>f, 'interactions'
	for a in interactions:
		for b in interactions[a]:
			print >>f, '\t %s %s %.3f' % ( a, b, interactions[a][b] )
	print >>f, 'end'

def PrintMolecule( f, name, count, sites, bonds, angles ):

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
	fmt = "\t %s %f %f %f"
	print >>f, 'coords'
	for c in coords:
		name, x, y, z = c
		print >>f, fmt % ( name, x, y, z )
	print >>f, 'end'

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

#
# Set up defaults
#

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

# After Venturoli et al

interactions['w']['w'] = 25.0
interactions['w']['h'] = 15.0
interactions['w']['t'] = 80.0

interactions['h']['h'] = 35.0
interactions['h']['t'] = 80.0

interactions['t']['t'] = 25.0

#
# Get user params
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
# Calculate system & write
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
