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
//                    SaveDPDSim( ref f, ref sim );
                }
            }
            catch( Exception ) { Console.WriteLine( "Unable to write checkpoint." ); }

            try
            {
                using( StreamWriter f = File.AppendText("output.lammpstrj") )
                {
//                    SaveDPDSim( ref f, ref sim );
                }
            }
            catch( Exception ) { Console.WriteLine( "Unable to write trajectory frame." ); }
        }

        static void Main( string[] args )
        {
        	if( Environment.GetCommandLineArgs().Length < 2 )
        	{
        		Console.WriteLine( "Pass in simulation definition file as first parameter!" );
        		System.Environment.Exit( -1 );
        	}

            DPDError( "int {0}, float {1}, string {2}", 12, 12.12, "twelve" );
        }
    }
}
