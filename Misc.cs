using System;
using System.Diagnostics;
using System.Runtime.CompilerServices;

namespace DPD
{

public class Misc
{
    //
    // Given the positions of 2 sites, and the simulation cell, this routine calculates
    // the minimum distance between the sites given periodic boundary conditions.
    //
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void VecMinimumImage(
            double x,      double y,      double z,
        ref double x_, ref double y_, ref double z_,
        ref double[] cell )
    {
        x_ = x - cell[0] * Math.Round( x/cell[0], MidpointRounding.AwayFromZero );
        y_ = y - cell[1] * Math.Round( y/cell[1], MidpointRounding.AwayFromZero );
        z_ = z - cell[2] * Math.Round( z/cell[2], MidpointRounding.AwayFromZero );
    }

    //
    // On error, print information about where the error occurred as well
    // as the error message.
    //
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
