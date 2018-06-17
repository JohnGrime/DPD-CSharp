using System;
using System.Collections.Generic;

using static DPD.Misc;

namespace DPD
{

public class Forces
{
    public static readonly int MAX_EXCLUSION_ENTRY = 4;

    //
    // Harmonic bonds:
    // U(r) = 0.5k( r - r_eq )**2
    //
    public static void DoBonds( ref DPDSim sim )
    {
        var dr = new double[3];
        var dr_hat = new double[3];
        var f = new double[3];
            
        var N_bonds = sim.bond_site_indices.Count / 2;
        for( var bond_i=0; bond_i<N_bonds; bond_i++ )
        {
            var i = sim.bond_site_indices[ (bond_i*2)+0 ];
            var j = sim.bond_site_indices[ (bond_i*2)+1 ];
            
            dr[0] = sim.r[ (i*3)+0 ] - sim.r[ (j*3)+0 ];
            dr[1] = sim.r[ (i*3)+1 ] - sim.r[ (j*3)+1 ];
            dr[2] = sim.r[ (i*3)+2 ] - sim.r[ (j*3)+2 ];
            
            //
            // Minimum image separation
            //
            VecMinimumImage( ref dr, ref dr, ref sim.cell );
            double dr_mag = Math.Sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );
            dr_hat[0] = dr[0] / dr_mag;
            dr_hat[1] = dr[1] / dr_mag;
            dr_hat[2] = dr[2] / dr_mag;
            
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
            f[0] = -gamma*dr_hat[0];
            f[1] = -gamma*dr_hat[1];
            f[2] = -gamma*dr_hat[2];
            
            sim.f[ (i*3)+0 ] += f[0];
            sim.f[ (i*3)+1 ] += f[1];
            sim.f[ (i*3)+2 ] += f[2];

            sim.f[ (j*3)+0 ] -= f[0];
            sim.f[ (j*3)+1 ] -= f[1];
            sim.f[ (j*3)+2 ] -= f[2];
            
            
            //
            // Accumulate bond energy
            //
            sim.bond_energy += omega;

            //
            // Accumulate contribution to virial
            //
            sim.pressure[0] += dr[0] * f[0];
            sim.pressure[1] += dr[0] * f[1];
            sim.pressure[2] += dr[0] * f[2];

            sim.pressure[3] += dr[1] * f[0];
            sim.pressure[4] += dr[1] * f[1];
            sim.pressure[5] += dr[1] * f[2];

            sim.pressure[6] += dr[2] * f[0];
            sim.pressure[7] += dr[2] * f[1];
            sim.pressure[8] += dr[2] * f[2];
        }
    }

    //
    // Harmonic angular potential function.
    // U = 0.5k( theta - theta_eq )^2.
    //
    public static void DoAngles( ref DPDSim sim )
    {
        const double theta_tol = 0.000001;
        var N_angles = sim.angle_site_indices.Count / 3;

        var r_ij = new double[3];
        var r_jk = new double[3];
        var fi = new double[3];
        var fk = new double[3];
        
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
            r_ij[0] = sim.r[ (j*3)+0 ] - sim.r[ (i*3)+0 ];
            r_ij[1] = sim.r[ (j*3)+1 ] - sim.r[ (i*3)+1 ];
            r_ij[2] = sim.r[ (j*3)+2 ] - sim.r[ (i*3)+2 ];
            VecMinimumImage( ref r_ij, ref r_ij, ref sim.cell );
            double r_ij_sq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
            double r_ij_mag = Math.Sqrt( r_ij_sq );

            //
            // Vector connecting k->j, using minimum image convention
            //
            r_jk[0] = sim.r[ (j*3)+0 ] - sim.r[ (k*3)+0 ];
            r_jk[1] = sim.r[ (j*3)+1 ] - sim.r[ (k*3)+1 ];
            r_jk[2] = sim.r[ (j*3)+2 ] - sim.r[ (k*3)+2 ];
            VecMinimumImage( ref r_jk, ref r_jk, ref sim.cell );
            double r_jk_sq = r_jk[0]*r_jk[0] + r_jk[1]*r_jk[1] + r_jk[2]*r_jk[2];
            double r_jk_mag = Math.Sqrt( r_jk_sq );
            
            double r_ij_dot_r_jk = r_ij[0]*r_jk[0] + r_ij[1]*r_jk[1] + r_ij[2]*r_jk[2];
            double r_ij_mag_r_jk = r_ij_mag * r_jk_mag;
            
            //
            // Following Gromacs; perhaps check cos_theta for -1 < cos_theta < 1?
            //
            double cos_theta = r_ij_dot_r_jk / r_ij_mag_r_jk;
            double theta = Math.Acos( cos_theta );
            double sin_theta = Math.Sin( theta );
            if( Math.Abs(sin_theta) < theta_tol ) sin_theta = theta_tol;
            
            double st = theta_k * ( theta - theta_eq ) / sin_theta;
            double sth = st * cos_theta;
            double cik = st  / ( r_ij_mag * r_jk_mag );
            double cii = sth / ( r_ij_mag * r_ij_mag );
            double ckk = sth / ( r_jk_mag * r_jk_mag );

            fi[0] = -cik * r_jk[0] + cii * r_ij[0];
            fi[1] = -cik * r_jk[1] + cii * r_ij[1];
            fi[2] = -cik * r_jk[2] + cii * r_ij[2];

            fk[0] = -cik * r_ij[0] + ckk * r_jk[0];
            fk[1] = -cik * r_ij[1] + ckk * r_jk[1];
            fk[2] = -cik * r_ij[2] + ckk * r_jk[2];

            //
            // Force acting on the central site calculated form other forces; we know
            // an angle potential is internal, so net acceleration on the 3 sites should
            // be zero, and hence if two forces known third must cancel out their sum.
            //
            sim.f[ (i*3)+0 ] += fi[0];
            sim.f[ (i*3)+1 ] += fi[1];
            sim.f[ (i*3)+2 ] += fi[2];

            sim.f[ (j*3)+0 ] += -fi[0] - fk[0];
            sim.f[ (j*3)+1 ] += -fi[1] - fk[1];
            sim.f[ (j*3)+2 ] += -fi[2] - fk[2];

            sim.f[ (k*3)+0 ] += fk[0];
            sim.f[ (k*3)+1 ] += fk[1];
            sim.f[ (k*3)+2 ] += fk[2];

            //
            // Pressure tensor contribution
            //
            sim.pressure[0] += r_ij[0] * fi[0] + r_jk[0] * fk[0];
            sim.pressure[1] += r_ij[0] * fi[1] + r_jk[0] * fk[1];
            sim.pressure[2] += r_ij[0] * fi[2] + r_jk[0] * fk[2];

            sim.pressure[3] += r_ij[1] * fi[0] + r_jk[1] * fk[0];
            sim.pressure[4] += r_ij[1] * fi[1] + r_jk[1] * fk[1];
            sim.pressure[5] += r_ij[1] * fi[2] + r_jk[1] * fk[2];

            sim.pressure[6] += r_ij[2] * fi[0] + r_jk[2] * fk[0];
            sim.pressure[7] += r_ij[2] * fi[1] + r_jk[2] * fk[1];
            sim.pressure[8] += r_ij[2] * fi[2] + r_jk[2] * fk[2];

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
    private static void nonbonded_pair( int i, int j, double cutsq, double sqrt_dt, ref DPDSim sim )
    {
        var r = new double[3];
        var r_hat = new double[3];
        var f = new double[3];
        var v = new double[3];
        
        // Ignore bound sites. I should really unroll this loop for speed, but this is "nicer" for changing MAX_EXCLUSION_ENTRY etc.
        for( var k=0; k<MAX_EXCLUSION_ENTRY; k++ )
        {
            if( sim.exclude[ (i*MAX_EXCLUSION_ENTRY)+k ] == j ) return;
        }
        
        // offset into the positions arrays for sites i and j. Flat arrays, so each position is 3 consecutive numbers.
        int i_off = i*3;
        int j_off = j*3;

        r[0] = sim.r[ i_off+0 ] - sim.r[ j_off+0 ];
        r[1] = sim.r[ i_off+1 ] - sim.r[ j_off+1 ];
        r[2] = sim.r[ i_off+2 ] - sim.r[ j_off+2 ];

        VecMinimumImage( ref r, ref r, ref sim.cell );

        double r_mag = (r[0]*r[0]) + (r[1]*r[1]) + (r[2]*r[2]);

        // avoids sqrt() where not needed - check the squares instead.
        if( r_mag > cutsq ) return;

        r_mag = Math.Sqrt( r_mag );
        r_hat[0] = r[0] / r_mag;
        r_hat[1] = r[1] / r_mag;
        r_hat[2] = r[2] / r_mag;

        // conservative and random weightings are the same,
        // dissipative weight = random weight^2 
        double cons_weight = 1.0 - r_mag/sim.rcut;
        double rand_weight = cons_weight;
        double diss_weight = rand_weight*rand_weight;

        // conservative contribution (interactions defined in kBT)
        double a = sim.interactions[ (sim.site_ids[i]*sim.site_types.Count)+sim.site_ids[j] ];
        sim.nonbonded_energy += a * ( r_mag*( (r_mag/(2.0*sim.rcut)) - 1 ) + (sim.rcut/2.0) );

        f[0] = a*cons_weight*r_hat[0];
        f[1] = a*cons_weight*r_hat[1];
        f[2] = a*cons_weight*r_hat[2];

        // dissipative & random contribution
        v[0] = sim.v[i_off+0] - sim.v[j_off+0];
        v[1] = sim.v[i_off+1] - sim.v[j_off+1];
        v[2] = sim.v[i_off+2] - sim.v[j_off+2];

        double rhat_dot_v = (r_hat[0]*v[0]) + (r_hat[1]*v[1]) + (r_hat[2]*v[2]);

        // accumulate dissipative force
        f[0] += -sim.fric * diss_weight * rhat_dot_v * r_hat[0];
        f[1] += -sim.fric * diss_weight * rhat_dot_v * r_hat[1];
        f[2] += -sim.fric * diss_weight * rhat_dot_v * r_hat[2];
        // accumulate random force - use gasdev() routine for zero mean and unit variance.
        double random_number = (double) Ran1.gasdev( ref sim.ran1_value );
        f[0] += sim.sigma * rand_weight * random_number * r_hat[0] / sqrt_dt;
        f[1] += sim.sigma * rand_weight * random_number * r_hat[1] / sqrt_dt;
        f[2] += sim.sigma * rand_weight * random_number * r_hat[2] / sqrt_dt;

        //
        // Update sim force array, and add contributions to pressure tensor.
        //
        sim.f[i_off+0] += f[0];
        sim.f[i_off+1] += f[1];
        sim.f[i_off+2] += f[2];

        sim.f[j_off+0] -= f[0];
        sim.f[j_off+1] -= f[1];
        sim.f[j_off+2] -= f[2];
        
        // random pressure tensor contribution
        sim.pressure[0] += r[0] * f[0];
        sim.pressure[1] += r[0] * f[1];
        sim.pressure[2] += r[0] * f[2];

        sim.pressure[3] += r[1] * f[0];
        sim.pressure[4] += r[1] * f[1];
        sim.pressure[5] += r[1] * f[2];

        sim.pressure[6] += r[2] * f[0];
        sim.pressure[7] += r[2] * f[1];
        sim.pressure[8] += r[2] * f[2];

        // keep a tally of number of interactions
        sim.ninteractions += 1.0;
    }

    //
    // This is the naive loop over pairs in the system.
    // Simple, but *extremely* inefficient as number of particles or particle density increases.
    //
    public static void DoNonbonded( ref DPDSim sim )
    {
        var N = sim.site_ids.Count;
        double cutsq = sim.rcut * sim.rcut;   // precalculate
        double sqrt_dt = Math.Sqrt( sim.delta_t ); // precalculate
        
        for( int i=0; i<N; i++ )
        {
            for( int j=i+1; j<N; j++ )
            {
                nonbonded_pair( i, j, cutsq, sqrt_dt, ref sim );
            }
        }
    }

    //
    // Faster method to evaluate pair interactions.
    // For more info on the linked list cell algorithm see Allen & Tildesley, Computer Simulation of Liquids.
    // Assumes periodic wrap applied before this routine is called.
    //
    public static void DoNonbonded2( ref DPDSim sim )
    {
        var N = sim.site_ids.Count;

        double cutsq = sim.rcut * sim.rcut;
        double sqrt_dt = Math.Sqrt( sim.delta_t );
        
        int ncellx = (int) Math.Floor(sim.cell[0]/sim.rcut);
        int ncelly = (int) Math.Floor(sim.cell[1]/sim.rcut);
        int ncellz = (int) Math.Floor(sim.cell[2]/sim.rcut);
        if( ncellx < 3 || ncelly < 3 || ncellz < 3 ) DPDError( "A cell dimension has fewer than 3 cells; this is a linked list cell error." );         

        //
        // If we're using a dynamic simulation cell, uncomment the following line!
        //
        //SetupCells( sim );

        // reset head indices.
        for( var i=0; i<sim.cell_head.Count; i++ ) sim.cell_head[i] = -1;

        // assign sites to cells, and update head lists etc.
        sim.cell_next.Capacity = N;
        for( var i=0; i<N; i++ )
        {
            var j = i*3;
            var cellx = (int) Math.Floor( ((sim.r[j+0]/sim.cell[0]) + 0.5) * ncellx );
            var celly = (int) Math.Floor( ((sim.r[j+1]/sim.cell[1]) + 0.5) * ncelly );
            var cellz = (int) Math.Floor( ((sim.r[j+2]/sim.cell[2]) + 0.5) * ncellz );

            // update linked lists.
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
                    // Loop over contents of current cell
                    for( var i=sim.cell_head[cell_no]; i!=-1; i=sim.cell_next[i] )
                    {
                        // loop over contents of neighbour cells, and current cell ( at offset 0 in cell_neighbours )
                        for( var k = 0; k < DPDSim.nneighbours; k++ )
                        {
                            var j = k; // dummy initialization, to make sure correct type. Value overwritten below.
                            if( k == 0 ) { j = sim.cell_next[i]; } // same cell: loop over j>i in current cell
                            else { j = sim.cell_head[ sim.cell_neighbours[ (cell_no*DPDSim.nneighbours) + k ] ]; } // different cell: loop over all j in other cell
                            for( ; j!=-1; j = sim.cell_next[j] ) nonbonded_pair( i, j, cutsq, sqrt_dt, ref sim );
                        }
                    }
                }
            }
        }
    }
}

}
