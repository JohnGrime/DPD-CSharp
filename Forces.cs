using System;
using System.Collections.Generic;

using static DPD.Misc;

namespace DPD
{

public class Forces
{
    //
    // Harmonic bonds:
    // U(r) = 0.5k( r - r_eq )**2
    //
    public static void DoBonds( DPDSim sim )
    {
        double dx, dy, dz;
        double dx_hat, dy_hat, dz_hat;
        double fx, fy, fz;
            
        var N_bonds = sim.bond_site_indices.Length / 2;
        for( var bond_i=0; bond_i<N_bonds; bond_i++ )
        {
            var i = sim.bond_site_indices[ (bond_i*2)+0 ];
            var j = sim.bond_site_indices[ (bond_i*2)+1 ];
            
            dx = sim.r[ (i*3)+0 ] - sim.r[ (j*3)+0 ];
            dy = sim.r[ (i*3)+1 ] - sim.r[ (j*3)+1 ];
            dz = sim.r[ (i*3)+2 ] - sim.r[ (j*3)+2 ];
            
            //
            // Minimum image separation
            //
            VecMinimumImage( dx,dy,dz, ref dx, ref dy, ref dz, ref sim.cell );
            double dr_mag = Math.Sqrt( dx*dx + dy*dy + dz*dz );
            dx_hat = dx / dr_mag;
            dy_hat = dy / dr_mag;
            dz_hat = dz / dr_mag;
            
            //
            // omega = energy, gamma = derivative of energy wrt particle separation
            //
            double K = sim.bond_k[bond_i];
            double r0 = sim.bond_eq[bond_i];
            double omega = 0.5 * K*( (dr_mag-r0)*(dr_mag-r0) );
            double gamma = K*(dr_mag-r0);           
            
            //
            // Accumulate force on particles due to bond
            //
            fx = -gamma*dx_hat;
            fy = -gamma*dy_hat;
            fz = -gamma*dz_hat;
            
            sim.f[ (i*3)+0 ] += fx;
            sim.f[ (i*3)+1 ] += fy;
            sim.f[ (i*3)+2 ] += fz;

            sim.f[ (j*3)+0 ] -= fx;
            sim.f[ (j*3)+1 ] -= fy;
            sim.f[ (j*3)+2 ] -= fz;
            
            //
            // Accumulate bond energy
            //
            sim.bond_energy += omega;

            //
            // Accumulate contribution to virial
            //
            sim.pressure[0] += dx * fx;
            sim.pressure[1] += dx * fy;
            sim.pressure[2] += dx * fz;

            sim.pressure[3] += dy * fx;
            sim.pressure[4] += dy * fy;
            sim.pressure[5] += dy * fz;

            sim.pressure[6] += dz * fx;
            sim.pressure[7] += dz * fy;
            sim.pressure[8] += dz * fz;
        }
    }

    //
    // Harmonic angular potential function.
    // U = 0.5k( theta - theta_eq )^2.
    //
    public static void DoAngles( DPDSim sim )
    {
        const double theta_tol = 0.000001;
        var N_angles = sim.angle_site_indices.Length / 3;

        double dx_ij, dy_ij, dz_ij;
        double dx_jk, dy_jk, dz_jk;
        double fx_i, fy_i, fz_i;
        double fx_k, fy_k, fz_k;
        
        for( var angle_i=0; angle_i<N_angles; angle_i++ )
        {
            //
            // Angle formed by sites i-j-k
            //
            var i = sim.angle_site_indices[ (angle_i*3)+0 ];
            var j = sim.angle_site_indices[ (angle_i*3)+1 ];
            var k = sim.angle_site_indices[ (angle_i*3)+2 ];

            double theta_k = sim.angle_k[angle_i];
            double theta_eq = sim.angle_eq[angle_i];

            //
            // Vector connecting i->j, using minimum image convention
            //
            dx_ij = sim.r[ (j*3)+0 ] - sim.r[ (i*3)+0 ];
            dy_ij = sim.r[ (j*3)+1 ] - sim.r[ (i*3)+1 ];
            dz_ij = sim.r[ (j*3)+2 ] - sim.r[ (i*3)+2 ];
            VecMinimumImage( dx_ij, dy_ij, dz_ij, ref dx_ij, ref dy_ij, ref dz_ij, ref sim.cell );
            double ij_mag = Math.Sqrt( dx_ij*dx_ij + dy_ij*dy_ij + dz_ij*dz_ij );

            //
            // Vector connecting k->j, using minimum image convention
            //
            dx_jk = sim.r[ (j*3)+0 ] - sim.r[ (k*3)+0 ];
            dy_jk = sim.r[ (j*3)+1 ] - sim.r[ (k*3)+1 ];
            dz_jk = sim.r[ (j*3)+2 ] - sim.r[ (k*3)+2 ];
            VecMinimumImage( dx_jk, dy_jk, dz_jk, ref dx_jk, ref dy_jk, ref dz_jk, ref sim.cell );
            double jk_mag = Math.Sqrt( dx_jk*dx_jk + dy_jk*dy_jk + dz_jk*dz_jk );
            
            double ij_dot_jk = dx_ij*dx_jk + dy_ij*dy_jk + dz_ij*dz_jk;
            double ij_mag_jk = ij_mag * jk_mag;
            
            //
            // Following Gromacs; perhaps check cos_theta for -1 < cos_theta < 1?
            //
            double cos_theta = ij_dot_jk / ij_mag_jk;
            double theta = Math.Acos( cos_theta );
            double sin_theta = Math.Sin( theta );
            if( Math.Abs(sin_theta) < theta_tol ) sin_theta = theta_tol;
            
            double st = theta_k * ( theta - theta_eq ) / sin_theta;
            double sth = st * cos_theta;
            double cik = st  / ( ij_mag * jk_mag );
            double cii = sth / ( ij_mag * ij_mag );
            double ckk = sth / ( jk_mag * jk_mag );

            fx_i = -cik*dx_jk + cii*dx_ij;
            fy_i = -cik*dy_jk + cii*dy_ij;
            fz_i = -cik*dz_jk + cii*dz_ij;

            fx_k = -cik*dx_ij + ckk*dx_jk;
            fy_k = -cik*dy_ij + ckk*dy_jk;
            fz_k = -cik*dz_ij + ckk*dz_jk;

            //
            // Force acting on the central site calculated form other forces; we know
            // an angle potential is internal, so net acceleration on the 3 sites should
            // be zero, and hence if two forces known third must cancel out their sum.
            //
            sim.f[ (i*3)+0 ] += fx_i;
            sim.f[ (i*3)+1 ] += fy_i;
            sim.f[ (i*3)+2 ] += fz_i;

            sim.f[ (j*3)+0 ] += -fx_i - fx_k;
            sim.f[ (j*3)+1 ] += -fy_i - fy_k;
            sim.f[ (j*3)+2 ] += -fz_i - fz_k;

            sim.f[ (k*3)+0 ] += fx_k;
            sim.f[ (k*3)+1 ] += fy_k;
            sim.f[ (k*3)+2 ] += fz_k;

            //
            // Pressure tensor contribution
            //
            sim.pressure[0] += dx_ij*fx_i + dx_jk*fx_k;
            sim.pressure[1] += dx_ij*fy_i + dx_jk*fy_k;
            sim.pressure[2] += dx_ij*fz_i + dx_jk*fz_k;

            sim.pressure[3] += dy_ij*fx_i + dy_jk*fx_k;
            sim.pressure[4] += dy_ij*fy_i + dy_jk*fy_k;
            sim.pressure[5] += dy_ij*fz_i + dy_jk*fz_k;

            sim.pressure[6] += dz_ij*fx_i + dz_jk*fx_k;
            sim.pressure[7] += dz_ij*fy_i + dz_jk*fy_k;
            sim.pressure[8] += dz_ij*fz_i + dz_jk*fz_k;

            //
            // Add angle energy to accumulator
            //
            double delta = theta-theta_eq;
            sim.angle_energy += 0.5 * theta_k * ( delta*delta );
        }
    }

    //
    // This is the actual DPD force calculation for two nonbonded sites.
    //
    private static void nonbonded_pair( int i, int j, double cutsq, double sqrt_dt, DPDSim sim )
    {
        double dx, dy, dz;
        double dx_hat, dy_hat, dz_hat;
        double fx, fy, fz;
        double vx, vy, vz;
        
        var i_off = i*3;
        var j_off = j*3;

        //
        // Ignore bound sites. I should really unroll this loop for speed, but this is "nicer" for changing MaxExclusionEntries etc.
        //
        for( var k=0; k<DPDSim.MaxExclusionEntries; k++ )
        {
            if( sim.exclude[ (i*DPDSim.MaxExclusionEntries)+k ] == j ) return;
        }
        
        //
        // Difference in particle positions (minimum image convention is applied) and velocities.
        // Try for the early skip, so we don't waste time reading velocities etc from memory if
        // we don't need to (should also help reduce cache pressure etc).
        //
        dx = sim.r[ i_off+0 ] - sim.r[ j_off+0 ];
        dy = sim.r[ i_off+1 ] - sim.r[ j_off+1 ];
        dz = sim.r[ i_off+2 ] - sim.r[ j_off+2 ];
        VecMinimumImage( dx,dy,dz, ref dx, ref dy, ref dz, ref sim.cell );
        double r_mag = (dx*dx + dy*dy + dz*dz);

        // Avoid sqrt() until we know it's needed, compare squared distances instead
        if( r_mag > cutsq ) return;

        vx = sim.v[i_off+0] - sim.v[j_off+0];
        vy = sim.v[i_off+1] - sim.v[j_off+1];
        vz = sim.v[i_off+2] - sim.v[j_off+2];

        r_mag = Math.Sqrt( r_mag );
        dx_hat = dx / r_mag;
        dy_hat = dy / r_mag;
        dz_hat = dz / r_mag;

        double rhat_dot_v = (dx_hat*vx) + (dy_hat*vy) + (dz_hat*vz);

        //
        // Conservative and random weightings are the same, dissipative
        // weight = random weight^2
        //
        double cons_weight = 1.0 - r_mag/sim.rcut;
        double rand_weight = cons_weight;
        double diss_weight = rand_weight*rand_weight;

        //
        // Conservative interaction parameters (interactions defined in kBT)
        //
        double a = sim.interactions[ (sim.site_ids[i]*sim.site_types.Count)+sim.site_ids[j] ];
        sim.nonbonded_energy += a * ( r_mag*( (r_mag/(2.0*sim.rcut)) - 1 ) + (sim.rcut/2.0) );

        //
        // Accumulate conservative force
        //
        fx = a*cons_weight*dx_hat;
        fy = a*cons_weight*dy_hat;
        fz = a*cons_weight*dz_hat;

        //
        // Accumulate dissipative force
        //
        fx += -sim.fric * diss_weight * rhat_dot_v * dx_hat;
        fy += -sim.fric * diss_weight * rhat_dot_v * dy_hat;
        fz += -sim.fric * diss_weight * rhat_dot_v * dz_hat;

        //
        // Accumulate random force - use gasdev() routine for zero mean and unit variance.
        //
        double random_number = (double) Ran1.gasdev( ref sim.ran1_value );
        fx += sim.sigma * rand_weight * random_number * dx_hat / sqrt_dt;
        fx += sim.sigma * rand_weight * random_number * dy_hat / sqrt_dt;
        fx += sim.sigma * rand_weight * random_number * dz_hat / sqrt_dt;

        //
        // Update sim force arrays
        //
        sim.f[i_off+0] += fx;
        sim.f[i_off+1] += fy;
        sim.f[i_off+2] += fz;

        sim.f[j_off+0] -= fx;
        sim.f[j_off+1] -= fy;
        sim.f[j_off+2] -= fz;
        
        //
        // Pressure tensor contribution
        //
        sim.pressure[0] += dx * fx;
        sim.pressure[1] += dx * fy;
        sim.pressure[2] += dx * fz;

        sim.pressure[3] += dy * fx;
        sim.pressure[4] += dy * fy;
        sim.pressure[5] += dy * fz;

        sim.pressure[6] += dz * fx;
        sim.pressure[7] += dz * fy;
        sim.pressure[8] += dz * fz;

        //
        // Maintain a tally of the number of pair interactions
        //
        sim.ninteractions += 1.0;
    }

    //
    // This is the naive loop over pairs in the system. Simple, but
    // *extremely* inefficient as number of particles or particle density increases.
    //
    public static void DoNonbonded( DPDSim sim )
    {
        var N = sim.site_ids.Length;
        double cutsq = sim.rcut * sim.rcut;
        double sqrt_dt = Math.Sqrt( sim.delta_t );
        
        for( var i=0; i<N; i++ )
        {
            for( var j=i+1; j<N; j++ )
            {
                nonbonded_pair( i, j, cutsq, sqrt_dt, sim );
            }
        }
    }

    //
    // Faster method to evaluate pair interactions.
    // For more info on the linked list cell algorithm see Allen & Tildesley, Computer Simulation of Liquids.
    // Assumes periodic wrap applied before this routine is called.
    //
    public static void DoNonbonded2( DPDSim sim )
    {
        var N = sim.site_ids.Length;

        double cutsq = sim.rcut * sim.rcut;
        double sqrt_dt = Math.Sqrt( sim.delta_t );
        
        var ncellx = (int) Math.Floor(sim.cell[0]/sim.rcut);
        var ncelly = (int) Math.Floor(sim.cell[1]/sim.rcut);
        var ncellz = (int) Math.Floor(sim.cell[2]/sim.rcut);
        if( ncellx < 3 || ncelly < 3 || ncellz < 3 )
        {
            DPDError( "A cell dimension has fewer than 3 cells; this is a linked list cell error." );
        }
        var ncells = ncellx*ncelly*ncellz;

        //
        // If we're using a dynamic simulation cell, uncomment the following line!
        //
        //SetupCells( sim );

        if( sim.cell_head.Length != ncells )
        {
            DPDError( "sim.cell_head.Length ({0}) != ncells ({1}); did you call SetupCells()?",
                sim.cell_head.Length, ncells );
        }

        //
        // Reset cell head indices.
        //
        for( var i=0; i<sim.cell_head.Length; i++ ) sim.cell_head[i] = -1;

        //
        // Assign sites to cells, update cell head indices.
        //
        if( sim.cell_next.Length != N ) Array.Resize( ref sim.cell_next, N );
        for( var i=0; i<N; i++ )
        {
            var j = i*3;

            var cellx = (int) Math.Floor( ((sim.r[j+0]/sim.cell[0]) + 0.5) * ncellx );
            var celly = (int) Math.Floor( ((sim.r[j+1]/sim.cell[1]) + 0.5) * ncelly );
            var cellz = (int) Math.Floor( ((sim.r[j+2]/sim.cell[2]) + 0.5) * ncellz );
            var cell_no = cellx + (celly*ncellx) + (cellz*ncelly*ncellx);

            sim.cell_next[i] = sim.cell_head[cell_no];
            sim.cell_head[cell_no] = i;
        }

        //
        // Faster version of pairwise interactions using neighbour cells.
        //
        for( var cellz=0; cellz<ncellz; cellz++ )
        {
            for( var celly=0; celly<ncelly; celly++ )
            {
                for( var cellx=0; cellx<ncellx; cellx++ )
                {
                    var cell_no = cellx + (celly*ncellx) + (cellz*ncelly*ncellx);
                    //
                    // Loop over contents of current cell
                    //
                    for( var i=sim.cell_head[cell_no]; i!=-1; i=sim.cell_next[i] )
                    {
                        //
                        // Loop over contents of neighbour cells ( includes central cell at offset 0 in cell_neighbours )
                        //
                        for( var k = 0; k < DPDSim.nneighbours; k++ )
                        {
                            var j = k; // dummy initialization, to make sure correct type. Value overwritten below.
                            if( k == 0 ) { j = sim.cell_next[i]; } // same cell: loop over j>i in current cell
                            else { j = sim.cell_head[ sim.cell_neighbours[ (cell_no*DPDSim.nneighbours) + k ] ]; } // different cell: loop over all j in other cell
                            for( ; j!=-1; j = sim.cell_next[j] ) nonbonded_pair( i, j, cutsq, sqrt_dt, sim );
                        }
                    }
                }
            }
        }
    }
}

}
