using System;

namespace DPD
{
    
//
// The ran1() function from Numerical Recipes in C: The Art of Scientific Computing ( ISBN: 0-521-43108-5 ).
// Bad PRNG, but retained for consistency with results from the C++ version of the DPD code.
// Call with idum as negative number to start, then do not alter idum!
//
public class Ran1
{
    // For ran1()
    private static readonly int NTAB = 32;
    private static long iy = 0;
    private static long[] iv = new long[NTAB];

    // For gasdev()
    private static int iset = 0; 
    private static double gset; 

    //
    // Uniform PRNG
    //
    public static double ran1( ref long idum )
    {
        var IA = 16807;
        var IM = 2147483647;
        var AM = (1.0/IM);
        var IQ = 127773;
        var IR = 2836;
        var NDIV = (1+(IM-1)/NTAB);
        var EPS = 1.2e-7;
        var RNMX = (1.0-EPS);

        long j, k;
        double temp;
        
        if( idum <= 0 || iy == 0 )
        {
            //
            // Initialize.
            //
            if (-idum < 1) idum = 1; // Be sure to prevent idum = 0.
            else idum = -idum;

            for( j = NTAB+7; j >= 0; j-- )
            {
                //
                // Load the shuffle table ( after 8 warm-ups ).
                //
                k = idum / IQ; 
                idum = IA*(idum-k*IQ) - IR*k; 
                if( idum < 0 ) idum += IM; 
                if( j < NTAB ) iv[j] = idum; 
            } 
            iy = iv[0]; 
        }
        
        k = idum / IQ; // Start here when not initializing.
        idum = IA*(idum-k*IQ) - IR*k; // Compute idum = (IA*idum) % IM without overflows by Schrage's method
        if( idum < 0 ) idum += IM;
        j = iy / NDIV; // Will be in the range 0..NTAB-1.
        iy = iv[j]; // Output previously stored value and refill the shuffle table
        iv[j] = idum;
        if( (temp = AM*iy) > RNMX ) return RNMX; // Because users don’t expect end point values.
        else return temp;
    }

    //
    // Sample from Gaussian distribution
    //
    public static double gasdev( ref long idum )
    {
        double fac, rsq, v1, v2; 
        
        if( idum < 0 ) iset = 0; // Reinitialize. 

        if( iset == 0 )
        {
            //
            // We don’t have an extra deviate handy, so pick two uniform numbers in the square
            // extending from -1 to +1 in each direction
            //
            do
            {
                v1 = 2.0*ran1( ref idum )-1.0;
                v2 = 2.0*ran1( ref idum )-1.0;
                rsq = v1*v1+v2*v2; // see if they are in the unit circle,
            } while( rsq >= 1.0 || rsq == 0.0 ); // and if they are not, try again.
            fac = Math.Sqrt( -2.0* ( (float)Math.Log((double)rsq)/rsq ) );

            //
            // Box-Muller transformation to get two normal deviates.
            // Return one and save the other for next time.
            //
            gset = v1*fac;
            iset = 1; // Set flag.
            return (float) v2*fac;
        }
        else
        {
            //
            // We have an extra deviate handy, so unset the flag and return it.
            //
            iset = 0;
            return (float) gset;
        } 
    }
}

}
