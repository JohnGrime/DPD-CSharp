using System;
using System.IO;

using static DPD.Misc;

namespace DPD
{

class Program
{
    private static void save_checkpoint_and_trajectory( DPDSim sim )
    {

        try
        {
            using( StreamWriter f = File.CreateText("output.checkpoint") )
            {
                DPDIO.SaveSim( f, sim );
            }
        }
        catch( Exception ) { Console.WriteLine( "Unable to write checkpoint." ); }

        try
        {
            using( StreamWriter f = File.AppendText("output.lammpstrj") )
            {
                DPDIO.SaveTrajectoryFrame( f, sim );
            }
        }
        catch( Exception ) { Console.WriteLine( "Unable to write trajectory frame." ); }
    }

    static void Main( string[] args )
    {
    	if( args.Length < 1 )
    	{
    		Console.WriteLine( "Pass in simulation definition file as first parameter!" );
    		System.Environment.Exit( -1 );
    	}

        DPDSim sim = new DPDSim();

        //
        // Load DPD sim data
        //
        try
        {
            using( StreamReader f = File.OpenText(args[0]) )
            {
                DPDIO.LoadSim( f, sim );
            }
        }
        catch( Exception e )
        {
            Console.WriteLine( e );
            DPDError( "Unable to open input file '{0}'", args[0] );
        }

        //
        // Print initial system info
        //
        {
            sim.ClearEnergyAndPressure();

            if( sim.i_am_dumb == 1 ) Forces.DoNonbonded( sim );
            else Forces.DoNonbonded2( sim );

            Forces.DoBonds( sim );
            Forces.DoAngles( sim );

            DPDIO.PrintSimInfo( sim, 0.0 );
        }

        //
        // Main integration loop
        //
        var time_start = DateTime.Now;
        for( ; sim.step_no <= sim.max_steps; sim.step_no++ )
        {
            sim.ClearEnergyAndPressure();
            
            Integrator.VelocityVerlet( sim );
            
            if( sim.step_no % sim.save_every == 0 )
            {
                save_checkpoint_and_trajectory( sim );
            }
            if( sim.step_no % sim.print_every == 0 )
            {
                var time_now = DateTime.Now;
                var elapsed = time_now - time_start;
                DPDIO.PrintSimInfo( sim, elapsed.TotalSeconds );
            }       
        }

        {
            var time_now = DateTime.Now;
            var elapsed = time_now - time_start;
        }

        //
        // Save final information.
        //
        save_checkpoint_and_trajectory( sim );
    }
}

}
