using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace DPD
{

public class Misc
{
    //
    // Given the positions of 2 sites, and the simulation cell, this routine calculates
    // the minimum distance between the sites given periodic boundary conditions.
    //
    public static void VecMinimumImage( ref double[] v1, ref double[] dest, ref double[] cell )
    {
        var temp_vec = new double[3];

        temp_vec[0] = v1[0] - cell[0] * Math.Round( v1[0]/cell[0], MidpointRounding.AwayFromZero );
        temp_vec[1] = v1[1] - cell[1] * Math.Round( v1[1]/cell[1], MidpointRounding.AwayFromZero );
        temp_vec[2] = v1[2] - cell[2] * Math.Round( v1[2]/cell[2], MidpointRounding.AwayFromZero );

        dest[0] = temp_vec[0];
        dest[1] = temp_vec[1];
        dest[2] = temp_vec[2];
    }

    public static void DPDError( string fmt, params object[] args )
    {
        var sf = new StackTrace( fNeedFileInfo: true ).GetFrame( 1 );

        Console.WriteLine( "ERROR:" );
        Console.WriteLine( "File: {0}", sf.GetFileName() );
        Console.WriteLine( "Method: {0}", sf.GetMethod() );
        Console.WriteLine( "Line: {0}", sf.GetFileLineNumber() );
        Console.WriteLine( "Message: {0}", string.Format( fmt, args ) );

        Environment.Exit( -1 );
    }
}

}
