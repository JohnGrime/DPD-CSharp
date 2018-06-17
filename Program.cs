using System;

using static DPD.Misc;

namespace DPD
{
    class Program
    {
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
