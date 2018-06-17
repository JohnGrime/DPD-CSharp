using System;
using System.Collections.Generic;

namespace DPD
{

public class Forces
{
    //
    // Harmonic bonds:
    // U(r) = 0.5k( r - r_eq )**2
    //
    public static void DoBonds( ref DPDSim sim )
    {
        const double theta_tol = 0.000001;
        int N_angles = sim.angle_site_indices.Count / 3;

        double[] r_ij = new double[3];
        double[] r_jk = new double[3];
        double[] fi = new double[3];
        double[] fk = new double[3];
        
        for( int angle_i=0; angle_i<N_angles; angle_i++ )
        {
            //
            // Angle formed by sites i-j-k
            //
            int i = sim.angle_site_indices[ (angle_i*3)+0 ];
            int j = sim.angle_site_indices[ (angle_i*3)+1 ];
            int k = sim.angle_site_indices[ (angle_i*3)+2 ];

            double theta_k = sim.angle_k[angle_i];
            double theta_eq = sim.angle_eq[angle_i];

            //
            // Vector connecting i->j, using minimum image convention
            //
            r_ij[0] = sim.r[ (j*3)+0 ] - sim.r[ (i*3)+0 ];
            r_ij[1] = sim.r[ (j*3)+1 ] - sim.r[ (i*3)+1 ];
            r_ij[2] = sim.r[ (j*3)+2 ] - sim.r[ (i*3)+2 ];
//                VecMinimumImage( r_ij, r_ij, sim.cell );
            double r_ij_sq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
            double r_ij_mag = Math.Sqrt( r_ij_sq );

            //
            // Vector connecting k->j, using minimum image convention
            //
            r_jk[0] = sim.r[ (j*3)+0 ] - sim.r[ (k*3)+0 ];
            r_jk[1] = sim.r[ (j*3)+1 ] - sim.r[ (k*3)+1 ];
            r_jk[2] = sim.r[ (j*3)+2 ] - sim.r[ (k*3)+2 ];
//                VecMinimumImage( r_jk, r_jk, sim.cell );
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
    // Harmonic angular potential function.
    // U = 0.5k( theta - theta_eq )^2.
    //
    public static void DoAngles( ref DPDSim sim )
    {
        const double theta_tol = 0.000001;
        int N_angles = sim.angle_site_indices.Count / 3;

        double[] r_ij = new double[3];
        double[] r_jk = new double[3];
        double[] fi = new double[3];
        double[] fk = new double[3];
        
        for( int angle_i=0; angle_i<N_angles; angle_i++ )
        {
            //
            // Angle formed by sites i-j-k
            //
            int i = sim.angle_site_indices[ (angle_i*3)+0 ];
            int j = sim.angle_site_indices[ (angle_i*3)+1 ];
            int k = sim.angle_site_indices[ (angle_i*3)+2 ];

            double theta_k = sim.angle_k[angle_i];
            double theta_eq = sim.angle_eq[angle_i];

            //
            // Vector connecting i->j, using minimum image convention
            //
            r_ij[0] = sim.r[ (j*3)+0 ] - sim.r[ (i*3)+0 ];
            r_ij[1] = sim.r[ (j*3)+1 ] - sim.r[ (i*3)+1 ];
            r_ij[2] = sim.r[ (j*3)+2 ] - sim.r[ (i*3)+2 ];
//            VecMinimumImage( r_ij, r_ij, sim.cell );
            double r_ij_sq = r_ij[0]*r_ij[0] + r_ij[1]*r_ij[1] + r_ij[2]*r_ij[2];
            double r_ij_mag = Math.Sqrt( r_ij_sq );

            //
            // Vector connecting k->j, using minimum image convention
            //
            r_jk[0] = sim.r[ (j*3)+0 ] - sim.r[ (k*3)+0 ];
            r_jk[1] = sim.r[ (j*3)+1 ] - sim.r[ (k*3)+1 ];
            r_jk[2] = sim.r[ (j*3)+2 ] - sim.r[ (k*3)+2 ];
//            VecMinimumImage( r_jk, r_jk, sim.cell );
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

    public static void DoNonbonded( ref DPDSim sim )
    {
    }

    public static void DoNonbonded2( ref DPDSim sim )
    {
    }
}

}
