using System;
using System.IO;

using static DPD.Misc;

namespace DPD
{

class Program
{
    private static void save_checkpoint_and_trajectory( ref DPDSim sim )
    {

        try
        {
            using( StreamWriter f = File.CreateText("output.checkpoint") )
            {
//                    SaveSim( ref f, ref sim );
            }
        }
        catch( Exception ) { Console.WriteLine( "Unable to write checkpoint." ); }

        try
        {
            using( StreamWriter f = File.AppendText("output.lammpstrj") )
            {
//                    SaveTrajectoryFrame( ref f, ref sim );
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
//              LoadSim( ref f, ref sim );
            }
        }
        catch( Exception ) { DPDError( "Unable to open input file '{0}'", args[0] ); }

        //
        // Print initial system info
        //
        {
            sim.ClearEnergyAndPressure();

            if( sim.i_am_dumb == 1 ) Forces.DoNonbonded( ref sim );
            else Forces.DoNonbonded2( ref sim );

            Forces.DoBonds( ref sim );
            Forces.DoAngles( ref sim );

//            PrintDPDSimInfo( ref sim, 0.0 );
        }

        //
        // Main integration loop
        //
        var time_start = new DateTime();
        for( ; sim.step_no <= sim.max_steps; sim.step_no++ )
        {
            sim.ClearEnergyAndPressure();
            
            Integrator.VelocityVerlet( ref sim );
            
            if( sim.step_no % sim.save_every == 0 )
            {
                save_checkpoint_and_trajectory( ref sim );
            }
            if( sim.step_no % sim.print_every == 0 )
            {
                var time_now = new DateTime();
                var elapsed = time_start - time_now;
//                PrintSimInfo( ref sim, difftime( time_end, elapsed.TotalSeconds ) );
            }       
        }

        {
            var time_now = new DateTime();
            var elapsed = time_start - time_now;
            Console.WriteLine( "Took ${0} seconds\n", elapsed.TotalSeconds );
        }

        //
        // Save final information.
        //
        save_checkpoint_and_trajectory( ref sim );
    }
}

}
