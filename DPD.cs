using System;
using System.Collections.Generic;

using static DPD.Misc;

namespace DPD
{
    
public class DPDSiteType
{
    public string name = "";
    public int internal_id = 0;
}

public class DPDMoleculeType
{
    public string name = "";
    public int count = 0;

    public List<int> site_internal_ids = new List<int>();

    public List<int> bond_site_indices = new List<int>();
    public List<double> bond_eq = new List<double>();
    public List<double> bond_k = new List<double>();

    public List<int> angle_site_indices = new List<int>();
    public List<double> angle_eq = new List<double>();
    public List<double> angle_k = new List<double>();

    public void Clear()
    {
        name = "";
        count = 0;

        site_internal_ids.Clear();
        
        bond_site_indices.Clear();
        bond_eq.Clear();
        bond_k.Clear();

        angle_site_indices.Clear();
        angle_eq.Clear();
        angle_k.Clear();
    }
}

public class DPDSim
{
    public readonly static int MaxExclusionEntries = 4;

    //
    // Neighbour cell information
    //
    public readonly static int nneighbours = 14;
    public readonly static int[] neighbour_offsets = new int[] {
         0, 0, 0, // current cell!
         1, 0, 0,
         1, 1, 0,
        -1, 1, 0,
         0, 1, 0,
         0, 0, 1,
        -1, 0, 1,
         1, 0, 1,
        -1,-1, 1,
         0,-1, 1,
         1,-1, 1,
        -1, 1, 1,
         0, 1, 1,
         1, 1, 1
    };

    //
    // System definition (we assume these don't change)
    //
    public double[] cell = new double[3];
    public List<DPDSiteType> site_types = new List<DPDSiteType>();
    public List<DPDMoleculeType> molecule_types = new List<DPDMoleculeType>();

    //
    // Per-particle information
    //
    public int[] site_ids = new int[1];
    public double[] r = new double[1];
    public double[] v = new double[1];
    public double[] f = new double[1];

    //
    // Nonbonded interaction information
    //
    public double rcut, ninteractions;
    public double[] interactions = new double[1];
    public double[] exclude = new double[1];

    //
    // Bond interaction information
    //
    public int[] bond_site_indices = new int[1];
    public double[] bond_eq = new double[1];
    public double[] bond_k  = new double[1];
    
    //
    // Angle interaction information
    //
    public int[] angle_site_indices = new int[1];
    public double[] angle_eq = new double[1];
    public double[] angle_k  = new double[1];

    //
    // Virial-related information
    //
    public double target_kBT;
    public double kinetic_energy, nonbonded_energy, bond_energy, angle_energy;
    public double[] pressure = new double[9];

    //
    // DPD-specific information
    //
    public double lambda, sigma, fric;
    
    //
    // Intervals & integration
    //
    public int step_no, max_steps, save_every, print_every;
    public double delta_t;
    
    //
    // Misc
    //
    public double[] v_ = new double[1];
    public double[] f_ = new double[1];
    public int[] cell_next = new int[1];
    public int[] cell_head = new int[1];
    public int[] cell_neighbours = new int[1];
    public long ran1_value;
    public int i_am_dumb;

    public DPDSim()
    {
        Clear();
    }
        
    public void Clear()
    {
        site_types.Clear();
        molecule_types.Clear();

        Array.Resize( ref site_ids, 0 );
        Array.Resize( ref r, 0 );
        Array.Resize( ref v, 0 );
        Array.Resize( ref f, 0 );

        Array.Resize( ref v_, 0 );
        Array.Resize( ref f_, 0 );

        Array.Resize( ref interactions, 0 );
        Array.Resize( ref exclude, 0 );

        Array.Resize( ref bond_site_indices, 0 );
        Array.Resize( ref bond_eq, 0 );
        Array.Resize( ref bond_k, 0 );

        Array.Resize( ref angle_site_indices, 0 );
        Array.Resize( ref angle_eq, 0 );
        Array.Resize( ref angle_k, 0 );

        for( var i=0; i<3; i++ ) cell[i] = 10.0;
        for( var i=0; i<9; i++ ) pressure[i] = 0.0;

        step_no = 0;
        max_steps = 10000;
        save_every = 1000;
        print_every = 1000;

        delta_t = 0.01;

        lambda = 0.65;
        sigma = 3.0;
        rcut = 1.0;
        target_kBT = 1.0;

        ran1_value = -1;
        
        i_am_dumb = 0;

        ClearEnergyAndPressure();
    }

    public void ClearEnergyAndPressure()
    {
        kinetic_energy = 0.0;
        nonbonded_energy = 0.0;
        bond_energy = 0.0;
        angle_energy = 0.0;

        for( var i=0; i<9; i++ ) pressure[i] = 0.0;

        ninteractions = 0.0;        
    }

    //
    // Use Gaussian distribution to set initial velocities.
    //
    public void SetInitialVelocities()
    {
        // 0.5mv**2 = 3/2kBT.
        // v = sqrt( 3kBT / m )

        double factor = Math.Sqrt( 3.0*target_kBT ) / 3.0; // divide evenly across degs of freedom for v components
        for( var i=0; i<site_ids.Length; i++ )
        {
            var j = i*3;

            v[ j+0 ] = Ran1.gasdev( ref ran1_value ) * factor;
            v[ j+1 ] = Ran1.gasdev( ref ran1_value ) * factor;
            v[ j+2 ] = Ran1.gasdev( ref ran1_value ) * factor;
        }
    }

    //
    // Remove any net momentum from the system (assumes all particles same mass)
    //
    public void ZeroMomentum()
    {
        var N_sites = site_ids.Length;
        var net_m = new double[3];

        net_m[0] = 0.0;
        net_m[1] = 0.0;
        net_m[2] = 0.0;

        for( var i=0; i<N_sites; i++ )
        {
            var j = i*3;
            net_m[0] += v[ j+0 ];
            net_m[1] += v[ j+1 ];
            net_m[2] += v[ j+2 ];
        }

        net_m[0] = net_m[0] / N_sites;
        net_m[1] = net_m[1] / N_sites;
        net_m[2] = net_m[2] / N_sites;

        for( var i=0; i<N_sites; i++ )
        {
            var j = i*3; 
            v[ j+0 ] -= net_m[0];
            v[ j+1 ] -= net_m[1];
            v[ j+2 ] -= net_m[2];
        }
    }

    //
    // Generate the neighbour cells.
    //
    public void SetupCells()
    {   
        //
        // Check for bad counts; minimum number of cells on any dimension is three.
        //
        var ncellx = (int) Math.Floor(cell[0]/rcut);
        var ncelly = (int) Math.Floor(cell[1]/rcut);
        var ncellz = (int) Math.Floor(cell[2]/rcut);
        var ncells = ncellx*ncelly*ncellz;
        
        if( ncellx < 3 || ncelly < 3 || ncellz < 3 ) DPDError( "A cell dimension has fewer than 3 cells; this is a linked list cell error." );         
        
        Array.Resize( ref cell_head, ncells );
        Array.Resize( ref cell_neighbours, ncells*nneighbours );

        for( var cellz=0; cellz<ncellz; cellz++ )
        {
            for( var celly=0; celly<ncelly; celly++ )
            {
                for( var cellx=0; cellx<ncellx; cellx++ )
                {
                    var cell_no = cellx + (celly*ncellx) + (cellz*ncelly*ncellx);
                    for( var l=0; l<nneighbours; l++ )
                    {
                        // get cell coords of current neighbour
                        var i = cellx + neighbour_offsets[ (l*3)+0 ];
                        var j = celly + neighbour_offsets[ (l*3)+1 ];
                        var k = cellz + neighbour_offsets[ (l*3)+2 ];
                        
                        // wrap cell lattice coords across periodic boundaries
                        if( i == -1 ) i = ncellx - 1;
                        else if( i == ncellx ) i = 0;

                        if( j == -1 ) j = ncelly - 1;
                        else if( j == ncelly ) j = 0;

                        if( k == -1 ) k = ncellz - 1;
                        else if( k == ncellz ) k = 0;

                        // record the neighbour index for this cell
                        cell_neighbours[ (cell_no*nneighbours) + l ] = i + (j*ncellx) + (k*ncelly*ncellx);
                    }
                }
            }
        }
    }
}

}
