using System;
using System.Collections.Generic;

namespace DPD
{

public class Integrator
{
    static public void VelocityVerlet( DPDSim sim )
    {
        var N_sites = sim.site_ids.Length;

        double dt = sim.delta_t;
        double h_dt_sq = 0.5*(dt*dt);
        
        //
        // Advance positions.
        //
        for( var site_i=0; site_i<N_sites; site_i++ )
        {
            var s = site_i*3;
            
            //
            // Calculate new positions
            //
            sim.r[s+0] += dt*sim.v[s+0] + h_dt_sq*sim.f[s+0];
            sim.r[s+1] += dt*sim.v[s+1] + h_dt_sq*sim.f[s+1];
            sim.r[s+2] += dt*sim.v[s+2] + h_dt_sq*sim.f[s+2];
            
            //
            // Apply periodic boundary conditions to new positions
            //
            sim.r[s+0] -= sim.cell[0] * Math.Round( sim.r[s+0]/sim.cell[0], MidpointRounding.AwayFromZero );
            sim.r[s+1] -= sim.cell[1] * Math.Round( sim.r[s+1]/sim.cell[1], MidpointRounding.AwayFromZero );
            sim.r[s+2] -= sim.cell[2] * Math.Round( sim.r[s+2]/sim.cell[2], MidpointRounding.AwayFromZero );

            //
            // Store current velocities and forces
            //
            sim.v_[s+0] = sim.v[s+0];
            sim.v_[s+1] = sim.v[s+1];
            sim.v_[s+2] = sim.v[s+2];

            sim.f_[s+0] = sim.f[s+0];
            sim.f_[s+1] = sim.f[s+1];
            sim.f_[s+2] = sim.f[s+2];
            
            //
            // Calculate temporary velocity prediction
            //
            sim.v[s+0] += sim.lambda*dt*sim.f[s+0];
            sim.v[s+1] += sim.lambda*dt*sim.f[s+1];
            sim.v[s+2] += sim.lambda*dt*sim.f[s+2];

            //
            // Zero force array, as nonbonded, bonded and angle force routines treat it as an accumulator.
            //
            sim.f[s+0] = 0.0;
            sim.f[s+1] = 0.0;
            sim.f[s+2] = 0.0;
        }

        //
        // Get new forces for these positions and temporary velocity prediction
        // Note: DoNonbondedSlow() is the slow method; DoNonbondedFast() is the neighbour cell version.
        //
        if( sim.i_am_dumb == 1 ) Forces.DoNonbondedSlow( sim );
        else Forces.DoNonbondedFast( sim );

        Forces.DoBonds( sim );
        Forces.DoAngles( sim );
        
        //
        // At this point we have:
        // sim->v_ == original velocities for this step, sim->v == velocity prediction.
        // sim->f_ == original forces for this step, sim->f == forces at new positions and predicted velocities.
        //
        
        //
        // Correct the velocities using the two force arrays.
        //
        sim.kinetic_energy = 0.0;
        for( var site_i=0; site_i<N_sites; site_i++ )
        {
            var s = site_i*3;
            
            sim.v[s+0] = sim.v_[s+0] + 0.5*dt*( sim.f_[s+0] + sim.f[s+0] );
            sim.v[s+1] = sim.v_[s+1] + 0.5*dt*( sim.f_[s+1] + sim.f[s+1] );
            sim.v[s+2] = sim.v_[s+2] + 0.5*dt*( sim.f_[s+2] + sim.f[s+2] );

            //
            // Add kinetic contribution to the pressure tensor.
            //
            double pvx = 0.5 * ( sim.v_[s+0] + sim.v[s+0] );
            double pvy = 0.5 * ( sim.v_[s+1] + sim.v[s+1] );
            double pvz = 0.5 * ( sim.v_[s+2] + sim.v[s+2] );

            sim.pressure[0] += pvx*pvx;
            sim.pressure[1] += pvx*pvy;
            sim.pressure[2] += pvx*pvz;

            sim.pressure[3] += pvy*pvx;
            sim.pressure[4] += pvy*pvy;
            sim.pressure[5] += pvy*pvz;

            sim.pressure[6] += pvz*pvx;
            sim.pressure[7] += pvz*pvy;
            sim.pressure[8] += pvz*pvz;

            sim.kinetic_energy += 0.5 * ( (sim.v[s+0]*sim.v[s+0]) + (sim.v[s+1]*sim.v[s+1]) + (sim.v[s+2]*sim.v[s+2]) );
        }
        
        //
        // Divide all pressure tensor elements by sim volume for correct units
        //
        double volume = sim.cell[0]*sim.cell[1]*sim.cell[2];
        for( var i=0; i<9; i++ ) sim.pressure[i] /= volume;            
    }
}

}
