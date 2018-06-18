using System;
using System.IO;
using System.Collections.Generic;

using static DPD.Misc;

namespace DPD
{

public class DPDIO
{

public static void LoadSim( StreamReader f, DPDSim sim )
{
    parse_dpd_sim( f, sim );

    //
    // How many bonds and angles do we need?
    //
    int N_bonds = 0;
    int N_angles = 0;
    foreach( var mol in sim.molecule_types )
    {
        N_bonds += mol.count * mol.bond_k.Count;
        N_angles += mol.count * mol.angle_k.Count;
    }
    
    Array.Resize( ref sim.bond_site_indices, N_bonds*2 );
    Array.Resize( ref sim.bond_eq, N_bonds );
    Array.Resize( ref sim.bond_k, N_bonds );

    Array.Resize( ref sim.angle_site_indices, N_angles*3 );
    Array.Resize( ref sim.angle_eq, N_angles );
    Array.Resize( ref sim.angle_k, N_angles );

    //
    // Populate system bond and angle arrays using molecule definitions and counts
    //
    int sys_bond_upto = 0;
    int sys_angle_upto = 0;
    int offset = 0; // start index into system sites for the current molecule
    foreach( var mol in sim.molecule_types )
    {
        for( var mol_inst = 0; mol_inst < mol.count; mol_inst++ )
        {
            var N = mol.bond_k.Count;
            for( var bi = 0; bi < N; bi++ )
            {
                var i = offset + mol.bond_site_indices[ (bi*2)+0 ];
                var j = offset + mol.bond_site_indices[ (bi*2)+1 ];

                var b_eq = mol.bond_eq[bi];
                var b_k  = mol.bond_k[bi];

                sim.bond_site_indices[sys_bond_upto*2 +0] = i-1; // -1 : unit based index -> zero based index
                sim.bond_site_indices[sys_bond_upto*2 +1] = j-1;
                sim.bond_eq[sys_bond_upto] = b_eq;
                sim.bond_k[sys_bond_upto] = b_k;

                sys_bond_upto++;
            }

            N = mol.angle_k.Count;
            for( var ai = 0; ai < N; ai++ )
            {
                var i = offset + mol.angle_site_indices[ (ai*3)+0 ];
                var j = offset + mol.angle_site_indices[ (ai*3)+1 ];
                var k = offset + mol.angle_site_indices[ (ai*3)+2 ];

                var a_eq = mol.angle_eq[ai];
                var a_k  = mol.angle_k[ai];

                sim.angle_site_indices[sys_angle_upto*3 +0] = i-1; // -1 : unit based index -> zero based index
                sim.angle_site_indices[sys_angle_upto*3 +1] = j-1;
                sim.angle_site_indices[sys_angle_upto*3 +2] = k-1;
                sim.angle_eq[sys_angle_upto] = a_eq;
                sim.angle_k[sys_angle_upto] = a_k;

                sys_angle_upto++;
            }

            offset += mol.site_internal_ids.Count; // update current site offset for bond indices.
        }
    }

    //
    // Nonbonded exclusion lists; used in nonbonded force calculations.
    // MAX_EXCLUSION_ENTRY suggested to be 4. If a lower value is posible, use it!
    // Only bonds used here; trivial to use angles too, but not needed for DPD.
    //
    var MaxExclusions = DPDSim.MaxExclusionEntries;
    Array.Resize( ref sim.exclude, sim.site_ids.Length*MaxExclusions );
    for( var i=0; i<sim.exclude.Length; i++ ) sim.exclude[i] = -1;
    for( var bi=0; bi<sim.bond_k.Length; bi++ )
    {
        var i = sim.bond_site_indices[(bi*2)+0];
        var j = sim.bond_site_indices[(bi*2)+1];
        
        // j exclusion
        var l = 0;
        while( l<MaxExclusions && sim.exclude[ (i*MaxExclusions)+l ] != -1 ) l++;
        if( l == MaxExclusions ) DPDError( "Too few exclusion entries defined (1)" );
        sim.exclude[ (i*MaxExclusions)+l ] = j;
        
        // k exclusion
        l = 0;
        while( l<MaxExclusions && sim.exclude[ (j*MaxExclusions)+l ] != -1 ) l++;
        if( l == MaxExclusions ) DPDError( "Too few exclusion entries defined (2)" );
        sim.exclude[ (j*MaxExclusions)+l ] = i;
    }

    sim.ZeroMomentum();

    //
    // Check whether the input file specified that the velocities should be initialised by the system...
    //
    if( sim.step_no < 1 )
    {
        Console.WriteLine( "step_no < 1, so assuming system initialisation desired; selecting site velocities from a Gaussian distribution." );
        sim.step_no = 1;
        sim.SetInitialVelocities();
    }

    //
    // Wrap particle positions that have crossed a periodic boundary
    //
    {
        for( var i=0; i<sim.site_ids.Length; i++ )
        {
            var j = i*3;
            var x = sim.r[j+0];
            var y = sim.r[j+1];
            var z = sim.r[j+2];
            VecMinimumImage( x,y,z, ref x, ref y, ref z, ref sim.cell );
            sim.r[j+0] = x;
            sim.r[j+1] = y;
            sim.r[j+2] = z;
        }
    }

    //
    // Link cell data only needs to be set up once if the cell is constant
    //
    sim.SetupCells();
    
    //
    // calculate friction_coefficient from sigma
    // sig^2 = 2*fric*kBT
    //
    sim.fric = (sim.sigma*sim.sigma) / (2.0*sim.target_kBT);
    Console.WriteLine( "Friction coefficient is {0} ( as sigma = {1} )", sim.fric, sim.sigma );
        
    Console.WriteLine( "Bead density is {0} per cubic Rc", ((double)sim.site_ids.Length) / (sim.cell[0]*sim.cell[1]*sim.cell[2]) );
    
    // boot ran1() with the seed provided
    Ran1.ran1( ref sim.ran1_value );
}

public static void SaveSim( StreamWriter f, DPDSim sim )
{
    var N_site_types = sim.site_types.Count;
    var N_mol_types = sim.molecule_types.Count;

    //
    // Settings
    //
    f.WriteLine( "settings" );
    f.WriteLine( "\tstep_no {0}", sim.step_no );
    f.WriteLine( "\tmax_steps {0}", sim.max_steps );
    f.WriteLine( "\tdelta_t {0}", sim.delta_t );
    f.WriteLine( "\tlambda {0}", sim.lambda );
    f.WriteLine( "\tsigma {0}", sim.sigma );
    f.WriteLine( "\tkBT {0}", sim.target_kBT );
    f.WriteLine( "\tran1 {0}", sim.ran1_value );
    f.WriteLine( "\tsave_every {0}", sim.save_every );
    f.WriteLine( "\tprint_every {0}", sim.print_every );
    f.WriteLine( "\tcell {0:G} {1:G} {2:G}\n", sim.cell[0], sim.cell[1], sim.cell[2] );
    f.WriteLine( "end" );

    f.WriteLine( "" );

    //
    // Sites
    //
    f.WriteLine( "sites" );
    for( var i=0; i<N_site_types; i++ )
    {
        f.WriteLine( "\t{0}", sim.site_types[i].name );
    }
    f.WriteLine( "end" );

    f.WriteLine( "" );

    //
    // Interactions
    //
    f.WriteLine( "interactions" );
    for( var i=0; i<N_site_types; i++ )
    {
        for( var j=i; j<N_site_types; j++ )
        {
            f.WriteLine( "\t{0}\t{1}\t{2}\n",
                sim.site_types[i].name, sim.site_types[j].name, sim.interactions[ (i*N_site_types)+j ] );
        }
    }
    f.WriteLine( "end" );

    f.WriteLine( "" );

    //
    // Molecules
    //
    for( var i=0; i<N_mol_types; i++ )
    {
        write_molecule_type( f, sim.molecule_types[i], sim );
        f.WriteLine( "" );
    }
    
    //
    // Coords
    //
    f.WriteLine( "coords" );
    for( var i=0; i<sim.site_ids.Length; i++ )
    {
        var site_id = sim.site_ids[i];
        f.WriteLine( "\t{0}", sim.site_types[site_id].name );
        f.WriteLine( "\t{0,12:F6} {1,12:F6} {2,12:F6}", sim.r[(i*3)+0], sim.r[(i*3)+1], sim.r[(i*3)+2] );
        f.WriteLine( "\t{0,12:F6} {1,12:F6} {2,12:F6}", sim.v[(i*3)+0], sim.v[(i*3)+1], sim.v[(i*3)+2] );
        f.WriteLine( "\t{0,12:F6} {1,12:F6} {2,12:F6}", sim.f[(i*3)+0], sim.f[(i*3)+1], sim.f[(i*3)+2] );
    }
    f.WriteLine( "end" );

    f.WriteLine( "" );
}

public static void SaveTrajectoryFrame( StreamWriter f, DPDSim sim )
{
    f.WriteLine( "ITEM: TIMESTEP" );
    f.WriteLine( "{0}", sim.step_no );

    f.WriteLine( "ITEM: NUMBER OF ATOMS" );
    f.WriteLine( "{0}", sim.site_ids.Length );

    f.WriteLine( "ITEM: BOX BOUNDS pp pp pp" );
    f.WriteLine( "{0:G} {1:G}", -sim.cell[0]/2, +sim.cell[0]/2 );
    f.WriteLine( "{0:G} {1:G}", -sim.cell[1]/2, +sim.cell[1]/2 );
    f.WriteLine( "{0:G} {1:G}", -sim.cell[2]/2, +sim.cell[2]/2 );

    f.WriteLine( "ITEM: ATOMS id type x y z" );
    for( var i=0; i<sim.site_ids.Length; i++ )
    {
        f.WriteLine( "{0} {1} {2:G} {3:G} {4:G}",
            i+1, sim.site_ids[i]+1, sim.r[(i*3)+0], sim.r[(i*3)+1], sim.r[(i*3)+2] );
    }
}

public static void PrintSimInfo( DPDSim sim, double cpu_time )
{
    var com = new double[3];
    var mom = new double[3];

    var N_sites = sim.site_ids.Length;
    var N_bonds = sim.bond_site_indices.Length;
    var N_angles = sim.angle_site_indices.Length;
    
    //print some info
    com[0] = com[1] = com[2] = 0.0;
    mom[0] = mom[1] = mom[2] = 0.0;
    double kBT = 0.0;
    for( var i=0; i<N_sites; i++ )
    {
        com[0] += sim.r[(i*3)+0];
        com[1] += sim.r[(i*3)+1];
        com[2] += sim.r[(i*3)+2];

        var vx = sim.v[(i*3)+0];
        var vy = sim.v[(i*3)+1];
        var vz = sim.v[(i*3)+2];

        mom[0] += vx;
        mom[1] += vy;
        mom[2] += vz;
        
        kBT += 0.5 * ( vx*vx + vy*vy + vz*vz );
    }
    
    kBT = (2.0/3.0) * (kBT/(N_sites-1));
    
    com[0] = com[0] / N_sites;
    com[1] = com[1] / N_sites;
    com[2] = com[2] / N_sites;

    Console.WriteLine( "Step {0}/{1}, sim time {2,12:F2}, CPU time {3:F0}s:",
        sim.step_no, sim.max_steps, sim.step_no * sim.delta_t, cpu_time );

    Console.WriteLine( "\tTotal energy     : {0,12:F6}",
        sim.kinetic_energy + sim.nonbonded_energy + sim.bond_energy + sim.angle_energy );

    Console.WriteLine( "\tKinetic energy   : {0,12:F6} ( Average {1,12:F6}, target kBT {2,12:F6}, current kBT {3,12:F6} )",
        sim.kinetic_energy, sim.kinetic_energy / N_sites, sim.target_kBT, kBT );

    Console.WriteLine( "\tNonbonded energy : {0,12:F6} ( Average {1,12:F6} from {2:F0} collisions )",
        sim.nonbonded_energy, sim.nonbonded_energy / sim.ninteractions, sim.ninteractions );
    
    if( N_bonds > 0 ) Console.WriteLine( "\tBond energy      : {0,12:F6} ( Average {1,12:F6} )",
        sim.bond_energy, sim.bond_energy / N_bonds );
    else Console.WriteLine( "\tBond energy      : 0" );

    if( N_angles > 0 ) Console.WriteLine( "\tAngle energy     : {0,12:F6} ( Average {1,12:F6} )",
        sim.angle_energy, sim.angle_energy / N_angles );
    else Console.WriteLine( "\tangle energy     : 0" );
    
    Console.WriteLine( "\tPressure         : {0,12:F6} {1,12:F6} {2,12:F6}",
        sim.pressure[0], sim.pressure[1], sim.pressure[2] );
    Console.WriteLine( "\t                   {0,12:F6} {1,12:F6} {2,12:F6}",
        sim.pressure[3], sim.pressure[4], sim.pressure[5] );
    Console.WriteLine( "\t                   {0,12:F6} {1,12:F6} {2,12:F6}",
        sim.pressure[6], sim.pressure[7], sim.pressure[8] );

    Console.WriteLine( "\tSystem centre of mass = {0,12:F6} {1,12:F6} {2,12:F6}", com[0], com[1], com[2] );

    Console.WriteLine( "\tNet system momentum   = {0,12:F6} {1,12:F6} {2,12:F6}", mom[0], mom[1], mom[2] );

    Console.WriteLine( "\n" );     
}

//
// Internal helper routines
//

//
// Return the internal id ( ie. index into the site type array ) for a site of type "name"
//
private static int get_site_type_from_name( string name, DPDSim sim )
{
    for( var i=0; i<sim.site_types.Count; i++ )
    {
        if( sim.site_types[i].name == name ) return sim.site_types[i].internal_id;
    }
    return -1;
}

private enum ParseState { None, Settings, Sites, Molecule, Interactions, Coords };

private static void parse_dpd_sim( StreamReader f, DPDSim sim )
{
    var delim = new char[] { ' ', ',', '\t' };
    string linebuffer;
    int cell_found = 0;

    ParseState parse_state = ParseState.None;
    
    int site_upto = 0, line_no = 0;

    var site = new DPDSiteType();
    var mol = new DPDMoleculeType();
    mol.Clear();

    sim.Clear();

    while( (linebuffer = f.ReadLine()) != null )
    {
        line_no++;
        var tokens = linebuffer.Split( delim, StringSplitOptions.RemoveEmptyEntries );

        if( tokens.Length < 1 ) continue;
        if( tokens[0].Length < 1 || tokens[0][0] == '#' ) continue;
        
        // if we hit an end flag
        if( tokens[0] == "end" )
        {
            if( parse_state == ParseState.None ) DPDError( "Unexpected 'end' token" );

            if( parse_state == ParseState.Molecule )
            {
                sim.molecule_types.Add( mol );
                mol = new DPDMoleculeType();
            }

            parse_state = ParseState.None;
        }
        // else if in settings section...
        else if( parse_state == ParseState.Settings )
        {
            var key = tokens[0];

            if( key == "kBT" ) sim.target_kBT = Convert.ToDouble( tokens[1] );
            else if( key == "sigma" ) sim.sigma = Convert.ToDouble( tokens[1] );
            else if( key == "lambda" ) sim.lambda = Convert.ToDouble( tokens[1] );
            else if( key == "delta_t" ) sim.delta_t = Convert.ToDouble( tokens[1] );
            else if( key == "step_no" ) sim.step_no = Convert.ToInt32( tokens[1] );
            else if( key == "max_steps" ) sim.max_steps = Convert.ToInt32( tokens[1] );
            else if( key == "save_every" ) sim.save_every  = Convert.ToInt32( tokens[1] );
            else if( key == "print_every" )sim.print_every = Convert.ToInt32( tokens[1] );
            else if( key == "ran1" ) sim.ran1_value = Convert.ToInt64( tokens[1] );
            else if( key == "dumb" ) sim.i_am_dumb = 1;
            else if( key == "cell" ) 
            {
                cell_found = 1;
                sim.cell[0] = Convert.ToDouble( tokens[1] );
                sim.cell[1] = Convert.ToDouble( tokens[2] );
                sim.cell[2] = Convert.ToDouble( tokens[3] );
            }
            else
            {
                Console.WriteLine( "Unknown entry {0} found on line {1} in settings; ignoring.", key, line_no );
            }
        }
        // if in site section...
        else if( parse_state == ParseState.Sites )
        {
            site.name = tokens[0];
            site.internal_id = sim.site_types.Count;
            sim.site_types.Add( site );
            site = new DPDSiteType();
        }
        else if( parse_state == ParseState.Interactions )
        {
            var N_site_types = sim.site_types.Count;

            var i = get_site_type_from_name( tokens[0], sim );
            var j = get_site_type_from_name( tokens[1], sim );

            // j and l are the internal ids, so use them.
            sim.interactions[ (i*N_site_types)+j ] = Convert.ToDouble( tokens[2] ); // symmetrical!
            sim.interactions[ (j*N_site_types)+i ] = sim.interactions[ (i*N_site_types)+j ];
        }
        // if in molecule section...
        else if( parse_state == ParseState.Molecule )
        {
            var key = tokens[0];

            if( key == "name" ) mol.name = tokens[1];
            else if( key == "count" ) mol.count = Convert.ToInt32( tokens[1] );
            else if( key == "site" )
            {
                var i = get_site_type_from_name( tokens[1], sim );
                if( i == -1 ) 
                {
                    DPDError( "Error on line {0}; unable to find a type id for molecule site of type {1} in {2}. Probably undefined.",
                        line_no, tokens[1], mol.name );
                }
                mol.site_internal_ids.Add( i );
            }
            else if( key == "bond" )
            {
                var N = mol.site_internal_ids.Count;

                // check sites exist!
                for( var i=0; i<2; i++ )
                {
                    var idx = Convert.ToInt32( tokens[1+i] );
                    if( idx < 1 || idx > N )
                    {
                        DPDError( "Error on line {0}; bond index {1} is invalid.", line_no, idx );
                    }
                    mol.bond_site_indices.Add( idx );
                }

                mol.bond_eq.Add( Convert.ToDouble( tokens[3] ) );
                mol.bond_k.Add( Convert.ToDouble( tokens[4] ) );
            }
            else if( key == "angle" )
            {
                var N = mol.site_internal_ids.Count;
                
                // check sites exist!
                for( var i=0; i<3; i++ )
                {
                    var idx = Convert.ToInt32( tokens[1+i] );
                    if( idx < 1 || idx > N )
                    {
                        DPDError( "Error on line {0}; angle index {1} is invalid.", line_no, idx );
                    }
                    mol.angle_site_indices.Add( idx );
                }

                mol.angle_eq.Add( Convert.ToDouble( tokens[4] ) );
                mol.angle_k.Add( Convert.ToDouble( tokens[5] ) );
            }
            else
            {
                Console.WriteLine( "Unknown entry {0} found on line {1} in molecule; ignoring.", key, line_no );
            }
        }
        // if in coords section...
        else if( parse_state == ParseState.Coords )
        {
            if( site_upto > sim.site_ids.Length )
            {
                DPDError( "Error on line {0}; this site should not exist", line_no );
            }

            var i = get_site_type_from_name( tokens[0], sim );
            if( i == -1 )
            {
                DPDError( "Error on line {0}; this site ( {1} ) is not recognised as a defined site type", line_no, tokens[0] );
            }
            sim.site_ids[site_upto] = i;
            
            sim.r[(site_upto*3)+0] = Convert.ToDouble( tokens[1] );
            sim.r[(site_upto*3)+1] = Convert.ToDouble( tokens[2] );
            sim.r[(site_upto*3)+2] = Convert.ToDouble( tokens[3] );
            
            // ignore velocity and force info if start of sim; this can also be used to force reset of info.
            if( sim.step_no > 0 )
            {
                if( tokens.Length > 4 ) // try to get velocities from file
                {
                    if( tokens.Length < 7 )
                    {
                        DPDError( "Error on line {0}; incomplete velocity entry for site! Only {1} value, and expecting 3 ( vx vy vz )", line_no, tokens.Length-4 );
                    }
                    sim.v[(site_upto*3)+0] = Convert.ToDouble( tokens[4] );
                    sim.v[(site_upto*3)+1] = Convert.ToDouble( tokens[5] );
                    sim.v[(site_upto*3)+2] = Convert.ToDouble( tokens[6] );
                }
                if( tokens.Length > 7 ) // try to get forces from file
                {
                    if( tokens.Length < 10 )
                    {
                        DPDError( "Error on line {0}; incomplete force entry for site! Only {1} value, and expecting 3 ( fx fy fz )", line_no, tokens.Length-7 );
                    }
                    sim.f[(site_upto*3)+0] = Convert.ToDouble( tokens[7] );
                    sim.f[(site_upto*3)+1] = Convert.ToDouble( tokens[8] );
                    sim.f[(site_upto*3)+2] = Convert.ToDouble( tokens[9] );
                }
            }
            site_upto++;
        }
        // assume this line defines the parse section
        else
        {
            string key = tokens[0];

            if( key == "settings" ) { parse_state = ParseState.Settings; }
            else if( key == "sites" ) { parse_state = ParseState.Sites; }
            else if( key == "interactions" )
            {
                var N_site_types = sim.site_types.Count;

                if( N_site_types < 1 )
                {
                    DPDError( "Error on line {0}; attempting to define interactions without sites having been specified.", line_no );
                }

                Array.Resize( ref sim.interactions, N_site_types*N_site_types );
                for( var i=0; i<sim.interactions.Length; i++ ) sim.interactions[i] = 0.0;
                
                parse_state = ParseState.Interactions;
            }
            else if( key == "molecule" ) { parse_state = ParseState.Molecule; }
            else if( key == "coords" )
            {               
                var N_sites = 0;
                foreach( var mt in sim.molecule_types )
                {
                    Console.WriteLine( "{0} {1} {2}", mt.name, mt.count, mt.site_internal_ids.Count );
                    N_sites += (mt.site_internal_ids.Count * mt.count);   
                }

                Array.Resize( ref sim.site_ids, N_sites );

                Array.Resize( ref sim.r, N_sites*3 );
                Array.Resize( ref sim.v, N_sites*3 );
                Array.Resize( ref sim.f, N_sites*3 );
                
                Array.Resize( ref sim.v_, N_sites*3 );
                Array.Resize( ref sim.f_, N_sites*3 );

                parse_state = ParseState.Coords;
            }
            else
            {
                DPDError( "Unknown parse section '{0}'' on line {1}", key, line_no );
            }
        }
    }
        
    //
    // Check for some obvious problems with the input.
    //
    if( cell_found == 0 ) DPDError( "cell not defined in DPD sim file" );
    if( sim.cell[0]/2.0 < sim.rcut || sim.cell[1]/2.0 < sim.rcut || sim.cell[2]/2.0 < sim.rcut ) DPDError( "A cell dimension divided by two is smaller than rcut" );
    if( sim.site_types.Count == 0 ) DPDError( "No sites defined in file" );
    if( sim.molecule_types.Count == 0 ) DPDError( "No molecules defined in file" );
    if( sim.interactions.Length == 0 ) DPDError( "No interactions defined in file" );
    if( sim.site_ids.Length == 0 ) DPDError( "No site coordinates defined in file" );
    if( site_upto != sim.site_ids.Length )
    {
        DPDError( "Error in number of sites found; got {0}, wanted {1}\n", site_upto, sim.site_ids.Length );
    }

    //
    // Print system information.
    //
    {       
        Console.WriteLine( "DPD simulation parameters:" );
        Console.WriteLine( "\tstep_no is {0}", sim.step_no );
        Console.WriteLine( "\tmax_steps is {0}", sim.max_steps );
        Console.WriteLine( "\tsave_every is {0}", sim.save_every );
        Console.WriteLine( "\tprint_every is {0}", sim.print_every );

        Console.WriteLine( "\tdelta_t is {0}", sim.delta_t );
        Console.WriteLine( "\tran1 initialised with {0}", sim.ran1_value );
        Console.WriteLine( "\trcut is {0}", sim.rcut );
        Console.WriteLine( "\tlambda is {0}", sim.lambda );
        Console.WriteLine( "\tsigma is {0}", sim.sigma );
        Console.WriteLine( "\tkBT is {0}", sim.target_kBT );
        Console.WriteLine( "\tcell is {0}, {1}, {2}", sim.cell[0], sim.cell[1], sim.cell[2] );
        if( sim.i_am_dumb == 1 ) Console.WriteLine( "\t!!! This simulation is dumb regarding nonbonded interactions !!!" );      
    }

    Console.WriteLine( "{0} site types:", sim.site_types.Count );
    for( var i=0; i<sim.site_types.Count; i++ )
    {
        Console.WriteLine( "\t{0}; (internal {1})", sim.site_types[i].name, sim.site_types[i].internal_id );
    }

    Console.WriteLine( "{0} molecule types:", sim.molecule_types.Count );
    for( int mi=0; mi<sim.molecule_types.Count; mi++ )
    {
        Console.WriteLine( "\t{0}; {1} counts in system:", sim.molecule_types[mi].name, sim.molecule_types[mi].count );

        Console.WriteLine( "\t\t{0} sites", sim.molecule_types[mi].site_internal_ids.Count );
        for( int si=0; si<sim.molecule_types[mi].site_internal_ids.Count; si++ )
        {
            var id = sim.molecule_types[mi].site_internal_ids[si];
            Console.WriteLine( "\t\t\t{0}, (id {1})", sim.site_types[id].name, sim.site_types[id].internal_id );
        }

        var n_bonds = sim.molecule_types[mi].bond_site_indices.Count/2;
        Console.WriteLine( "\t\t{0} bonds", n_bonds );
        for( var bi=0; bi<n_bonds; bi++ )
        {
            Console.WriteLine( "\t\t\t{0} - {1} : eq length {2:F3}, k {3:F3}",
                sim.molecule_types[mi].bond_site_indices[bi*2+0],
                sim.molecule_types[mi].bond_site_indices[bi*2+1],
                sim.molecule_types[mi].bond_eq[bi],
                sim.molecule_types[mi].bond_k[bi] );
        }

        var n_angles = sim.molecule_types[mi].angle_site_indices.Count/3;
        Console.WriteLine( "\t\t{0} angles", n_angles );
        for( var ai=0; ai<n_angles; ai++ )
        {
            Console.WriteLine( "\t\t\t{0} - {1} - {2} : eq length {3:F3}, k {4:F3}",
                sim.molecule_types[mi].angle_site_indices[ai*2+0],
                sim.molecule_types[mi].angle_site_indices[ai*2+1],
                sim.molecule_types[mi].angle_site_indices[ai*2+1],
                sim.molecule_types[mi].angle_eq[ai],
                sim.molecule_types[mi].angle_k[ai] );
        }
    }

    {
        var N_site_types = sim.site_types.Count;
        Console.Write( "{0} interactions (defined in kBT):\n{1,12:S}", N_site_types, "              " );
        for( var i=0; i<N_site_types; i++ )
        {
            Console.Write( "\t{0,12:S}", sim.site_types[i].name );
        }
        Console.WriteLine( "" );
        for( var i=0; i<N_site_types; i++ )
        {
            Console.Write( "\t{0,12:S}", sim.site_types[i].name );
            for( var j=0; j<N_site_types; j++ )
            {
                Console.Write( "\t{0,9:F3}", sim.interactions[(i*N_site_types)+j] );
            }
            Console.WriteLine( "" );
        }
    }

    //
    // check the sites we have tally with the molecular data.
    //
    {
        var l = 0;
        for( var mi=0; mi<sim.molecule_types.Count; mi++ )
        {
            for( var j=0; j<sim.molecule_types[mi].count; j++ )
            {
            for( var k=0; k<sim.molecule_types[mi].site_internal_ids.Count; k++ )
            {
                var expected_type = sim.molecule_types[mi].site_internal_ids[k];
                var actual_type = sim.site_ids[l];
                var expected_name = sim.site_types[expected_type].name;
                var actual_name = sim.site_types[actual_type].name;
                if( actual_type != expected_type )
                {
                    DPDError( "Error in site id; wanted {0} ( \"{1}\" from molecule \"{2}\" number {3}, atom {4} ) but found {5} ( \"{6}\" )",
                        expected_type, expected_name,
                        sim.molecule_types[mi].name,
                        j+1, k+1,
                        actual_type, actual_name );
                }
                l++;
            }
            }
        }
    }

    Console.WriteLine( "Sites ({0}):", sim.site_ids.Length );
    if( sim.site_ids.Length > 10 )
    {
        Console.WriteLine( "\t( printing truncated to first 10 sites )\n" );
        for( var i=0; i<10; i++ )
        {
            Console.WriteLine( "\t{0} (internal {1}), position = {2:F3}, {2:F3}, {2:F3}",
                sim.site_types[ sim.site_ids[i] ].name, sim.site_ids[i], sim.r[i*3], sim.r[(i*3)+1], sim.r[(i*3)+2] );
        }
    }
    else
    {
        for( var i=0; i<sim.site_ids.Length; i++ )
        {
            Console.WriteLine( "\t{0} (internal {1}), position = {2:F3}, {2:F3}, {2:F3}",
                sim.site_types[ sim.site_ids[i] ].name, sim.site_ids[i], sim.r[i*3], sim.r[(i*3)+1], sim.r[(i*3)+2] );
        }
    }
}

private static void write_molecule_type( StreamWriter f, DPDMoleculeType mol, DPDSim sim )
{
    var N_sites = mol.site_internal_ids.Count;
    var N_bonds = mol.bond_site_indices.Count / 2;
    var N_angles = mol.angle_site_indices.Count / 3;

    f.WriteLine( "molecule" );
    f.WriteLine( "\tname {0}", mol.name );
    f.WriteLine( "\tcount {0}", mol.count );
    for( var si=0; si<N_sites; si++ )
    {
        var i = mol.site_internal_ids[si];
        f.WriteLine( "\tsite {0}", sim.site_types[i].name );
    }
    for( var bi=0; bi<N_bonds; bi++ )
    {
        var i = mol.bond_site_indices[bi*2 +0];
        var j = mol.bond_site_indices[bi*2 +1];
        f.WriteLine( "\tbond\t{0}\t{1}\t\t{2}\t{3}\n", i, j, mol.bond_eq[bi], mol.bond_k[bi] );
    }
    for( var ai=0; ai<N_angles; ai++ )
    {
        var i = mol.angle_site_indices[ai*3 +0];
        var j = mol.angle_site_indices[ai*3 +1];
        var k = mol.angle_site_indices[ai*3 +2];
        f.WriteLine( "\tangle\t{0}\t{1}\t{2}\t\t{3}\t{4}\n", i, j, k, mol.angle_eq[ai], mol.angle_k[ai] );
    }
    f.WriteLine( "end\n" );
}

}

}
