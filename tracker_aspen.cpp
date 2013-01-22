/*
  Taking some basics from Ben M's tracker_3ds6.C and going to try and explore things in prep for Aspen's 2013 "New Directions in Neutrino Physics" conference....
*/

// Standard
#include <iostream>
#include <stdio.h>
#include <cmath>
// 3rd Party
#include "gsl_errno.h"
// Internal
using namespace std;

/*
Time to declare tons of things at global scope. Yay for != best practices
*/
double c = 0.299792458; // speed of light in meters per ns 
double w0=27.9925;  //nrel cyclotron frequency in GHz
double me = 510998.902; // electron mass in eV
double mesi = 9.10838188e-31; // electron mass in kg
double qesi = 1.60217646e-19; // electron charge in Coulombs; 
double qe = 1; // electron charge has been folded into "B" units, I think
double dept = 0; //common-block e field strength (set in code)
double sweept[500],sweepf[500],sweepd[500]; //sweep parameters
double nstoevoc = 5.34428542119e-28; //momentum: convert ev/c to J*s 
bool eperturb = false;
double eperturb_q = 0 ;
bool bperturb = false;
double bperturb_b = 0 ;
int useconstamp=0;
int debug=0; //debug level flag
int nsweep;
int iswf;
int iswd;
double econst = 1000;

void get_e_field(const double y[], double evec[])
{
    // should be in units "volts per meter" 
    evec[0] = -econst * y[0]; //hyperbolic electrostatic potential = proportional to R
    evec[1] = -econst * y[1]; //repulsive in xy
    evec[2] =  econst * y[2]; // attractive in z  

    if (eperturb) { // add a perturbing potential as due to a point charge 10cm away
        double qor3 = eperturb_q / pow(y[0] * y[0] + y[1] * y[1] + pow(y[2] + 0.1, 2), 1.5);  
        evec[0] = y[0] * qor3;
        evec[1] = y[1] * qor3;
        //static Z component guarantees that oscillation center is in middle
        evec[2] = (y[2] + 0.1) * qor3 - eperturb_q / pow(0.1, 2); 
    }
}

void get_b_field(const double y[], double bvec[])
{
    // magnetic field should be in Tesla 
    bvec[0] = 0;  
    bvec[1] = 0; 
    bvec[2] = w0 / 177; // w0 is in GHz, 177 = frequency when B=1T.  
    //also return total field strength
    bvec[3] = bvec[2];

    if (bperturb) {
        bvec[0] += bperturb_b * y[0] / 2; //x, y components get stronger off axis
        bvec[1] += bperturb_b * y[1] / 2;
        bvec[2] -= bperturb_b * y[2];  // z component gets weaker at +z
        bvec[3] = sqrt(bvec[0] * bvec[0] + bvec[1] * bvec[1] + bvec[2] * bvec[2]);
    }
}

//calculates the drive frequency at the present time 
double wda(double t)
{
    //  cout << "wda " << nsweep << " " << me << endl;
    for (int i=0; i < nsweep-1; i++) {
        if (t >= sweept[i] && t < sweept[i+1]) { //return an interpolated ramp 
            double tmp =  sweepf[i] + (t - sweept[i]) * (sweepf[i+1] - sweepf[i]) / (sweept[i+1] - sweept[i]); 
            return w0 * me / (me + tmp);
        }
    }

    return w0 * me / (me + sweepf[nsweep-1]); 
}

//calculates the electron cyclotron radius at the present time 
double rda(double t)
{
//    cout << "wda " << nsweep << " " << me << endl;
    for (int i=0; i < nsweep-1; i++) {
        if (t >= sweept[i] && t < sweept[i+1]) { //return an interpolated ramp 
            double tmp =  sweepf[i] + (t - sweept[i])*(sweepf[i+1]-sweepf[i])/(sweept[i+1]-sweept[i]); 
        return sqrt(tmp);
        }
    }

    return sqrt(sweepf[nsweep-1]); 
}

//calculates the drive frequency at the present time 
double power(double t)
{
//    cout << nsweep;
    for (int i=0; i<nsweep-1; i++) {
//        cout << "." ;
        if (t >= sweept[i] && t < sweept[i+1]) { //return an interpolated ramp 
//            cout << "power calc: " << i << " " << t << " " << sweepd[i] + (t - sweept[i])*(sweepd[i+1]-sweepd[i])/(sweept[i+1]-sweept[i]) << endl;
            if (useconstamp != 0) { // set ramp amplitude to correct for drive freq.
                return (sweepd[i] + (t - sweept[i]) * (sweepd[i+1] - sweepd[i]) / (sweept[i+1] - sweept[i])) * 31.6227766 / rda(t); //normalized to give unmodified power at E=1000 
            }
            else {
                return sweepd[i] + (t - sweept[i]) * (sweepd[i+1] - sweepd[i]) / (sweept[i+1] - sweept[i]); 
            }
        }
    }

//    cout << " return last " << sweepd[nsweep-1] << endl; 
    return sweepd[nsweep-1]; 
}

int func (double t, const double y[], double f[], void *params)
{  
    //mostly decoupled: synchrotron oscillations 
    //current version is copied exactly from tracker07/tracker_rk4.C 
    //this will need to be coupled in better once we've figured out what ppar means exactly.
    //especially: when we allow w0 to vary spatially (inhomogenous B field) 
   
    double oogamma = me / (me + y[4]);
    double beta = sqrt(1 - oogamma * oogamma);
    f[4] = power(t) * sqrt(2 * y[4] * me) * sin(y[5]) * 0.3 -
         3.201e-9 * w0 * oogamma * w0 * oogamma * beta * beta / (1 - beta * beta); //in eV/nanosecond
    f[5] = wda(t) - w0 * oogamma; //in rad/nanosecond

    // y[0,1,2] = x, y, z
    // y[3] = ppar (momentum parallel to guiding center motion) 
    // y[4] = total particle kinetic energy
    // y[5] = synchrotron orbit phase
    // see Zhang et. al. appendix for possible improvements
    // maybe add parameters for synchrotron phase, energy
    double ex, ey, ez; 
    double bx, by, bz; 
    double ppar, btot, bsquared; 
    
    double bvec[4]; //3-vector; 4th component is total field strength.     
    double evec[3];

    get_e_field(y, evec); // puts LOCAL electric field into evec[3] by reference;
    get_b_field(y, bvec); // puts LOCAL magnetic field into bvec[3] by reference;

    // should be in units "volts per meter" 
    ex = evec[0];
    ey = evec[1];
    ez = evec[2];
    //should be in kg*m/s
    //oogamma provides all of the relativistic correction
    ppar = oogamma * y[3] / mesi;  //and we'll pre-divide by the SI electron mass
    // magnetic field should be in Tesla 
    bx = bvec[0];  
    by = bvec[1]; 
    bz = bvec[2]; // w0 is in GHz, 177 = frequency when B=1T.  
    btot = bvec[3];  //scalar magnetic field strength
    bsquared = bvec[3] * bvec[3];  //scalar magnetic field strength

    f[0] =                  - ey*bz/bsquared  + ez*by/bsquared + ppar*bx/btot; 
    f[1] =  ex*bz/bsquared                    - ez*bx/bsquared + ppar*by/btot; 
    f[2] = -ex*by/bsquared  + ey*bx/bsquared                   + ppar*bz/btot; 
    f[3] = -qesi*ex*bx/btot - qesi*ey*by/btot - qesi*ez*bz/btot; 

    //reduce it all by a factor of 10^9 since the time units are nanoseconds

    f[0] /= 1e9;
    f[1] /= 1e9;
    f[2] /= 1e9;
    f[3] /= 1e9;
    // Doppler shift term: relative phase between particle and microwave changes with Z according to the equation dPhi/dz = c/wda.   So dPhi/dt = dPhi/dz * dz/dt, and f[2] is dz/dt in meters per nanosecond.  Sign is arbitrary; it depends which way the microwaves are propagating.  In principle they could propagate in any direction, and you could have a Dopplar shift depending on dx/dt or whatever.    
    f[5] += wda(t) / c * f[2];  

    if (debug==60) {
        cout << t << " " << f[4] << " " << f[5] << " via " << wda(t) <<
            " - " << w0 / (me + y[4]) << " " <<  y[2] << " " << f[2] << endl;
    }
    if (debug==61) {
        cout << t << " " <<  qesi << " " <<  ppar << " " <<  bz << " " <<
           btot << " " <<  f[2] << " " <<  y[2] << " " <<  endl;
    }

    return GSL_SUCCESS; 
} 

int main() {
    printf("I don't do anything yet!\nThanks for playing.\n");
    return 0;
}
