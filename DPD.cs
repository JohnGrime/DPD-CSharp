using System;
using System.Collections.Generic;

using static DPD.Misc;

namespace DPD
{
    
public class DPDSiteType
{
    public string name;
    public int internal_id;
}

public class DPDMoleculeType
{
    public string name;
    public int count;

    public List<int> site_internal_ids;

    public List<int> bond_site_indices;
    public List<double> bond_eq, bond_k;

    public List<int> angle_site_indices;
    public List<double> angle_eq, angle_k;

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
    public List<int> site_ids = new List<int>();
    public List<double> r, v, f;

    //
    // Nonbonded interaction information
    //
    public double rcut, ninteractions;
    public List<double> interactions = new List<double>();
    public List<int> exclude = new List<int>();

    //
    // Bond interaction information
    //
    public List<int> bond_site_indices = new List<int>();
    public List<double> bond_eq = new List<double>();
    public List<double> bond_k = new List<double>();
    
    //
    // Angle interaction information
    //
    public List<int> angle_site_indices = new List<int>();
    public List<double> angle_eq = new List<double>();
    public List<double> angle_k = new List<double>();

    //
    // Virial-related inforamtion
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
    public List<double> v_ = new List<double>();
    public List<double> f_ = new List<double>();
    public List<int> cell_next = new List<int>();
    public List<int> cell_head = new List<int>();
    public List<int> cell_neighbours = new List<int>();
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

        site_ids.Clear();
        r.Clear();
        v.Clear();
        f.Clear();

        v_.Clear();
        f_.Clear();

        interactions.Clear();
        exclude.Clear();

        bond_site_indices.Clear();
        bond_eq.Clear();
        bond_k.Clear();

        angle_site_indices.Clear();
        angle_eq.Clear();
        angle_k.Clear();

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
        for( var i=0; i<site_ids.Count; i++ )
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
        var N_sites = site_ids.Count;
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
            
        cell_head.Capacity = ncells;
        cell_neighbours.Capacity = ncells*nneighbours;

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
