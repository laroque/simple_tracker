/*
  Taking some basics from Ben M's tracker_3ds6.C and going to try and explore things in prep for Aspen's 2013 "New Directions in Neutrino Physics" conference....
*/

// Standard
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
// 3rd Party
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_matrix.h>
// Internal

// Namespace(s)
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
double wda = w0;//eV,ns,GHz
double power0 = 1e-15;//check this for units

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

/* YOU DON'T HAVE ANY NEED FOR THIS
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
*/

/* YOU ALSO HAVE NO NEED FOR ANY OF THIS
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

*/


int func (double t, const double y[], double f[], void *params)
{  
    //mostly decoupled: synchrotron oscillations 
    //current version is copied exactly from tracker07/tracker_rk4.C 
    //this will need to be coupled in better once we've figured out what ppar means exactly.
    //especially: when we allow w0 to vary spatially (inhomogenous B field) 
   
    double oogamma = me / (me + y[4]);
    double beta = sqrt(1 - oogamma * oogamma);
    f[4] = power0 * sqrt(2 * y[4] * me) * sin(y[5]) * 0.3 -
         3.201e-9 * w0 * oogamma * w0 * oogamma * beta * beta / (1 - beta * beta); //in eV/nanosecond

    f[5] = wda - w0 * oogamma; //in rad/nanosecond
    //f[5] = wda(t) - w0 * oogamma; //in rad/nanosecond

/* AGAIN, DON'T NEED ANYTHING HERE
    // y[0,1,2] = x, y, z
    // y[3] = ppar (momentum parallel to guiding center motion) 
    // y[4] = total particle kinetic energy
    // y[5] = synchrotron orbit phase
    // see Zhang et. al. appendix for possible improvements
//    // maybe add parameters for synchrotron phase, energy
//    double ex, ey, ez; 
//    double bx, by, bz; 
//    double ppar, btot, bsquared; 
//    
//    double bvec[4]; //3-vector; 4th component is total field strength.     
//    double evec[3];
//
//    get_e_field(y, evec); // puts LOCAL electric field into evec[3] by reference;
//    get_b_field(y, bvec); // puts LOCAL magnetic field into bvec[3] by reference;
//
//    // should be in units "volts per meter" 
//    ex = evec[0];
//    ey = evec[1];
//    ez = evec[2];
//    //should be in kg*m/s
//    //oogamma provides all of the relativistic correction
//    ppar = oogamma * y[3] / mesi;  //and we'll pre-divide by the SI electron mass
//    // magnetic field should be in Tesla 
//    bx = bvec[0];  
//    by = bvec[1]; 
//    bz = bvec[2]; // w0 is in GHz, 177 = frequency when B=1T.  
//    btot = bvec[3];  //scalar magnetic field strength
//    bsquared = bvec[3] * bvec[3];  //scalar magnetic field strength
//
//    f[0] =                  - ey*bz/bsquared  + ez*by/bsquared + ppar*bx/btot; 
//    f[1] =  ex*bz/bsquared                    - ez*bx/bsquared + ppar*by/btot; 
//    f[2] = -ex*by/bsquared  + ey*bx/bsquared                   + ppar*bz/btot; 
//    f[3] = -qesi*ex*bx/btot - qesi*ey*by/btot - qesi*ez*bz/btot; 
//
//    //reduce it all by a factor of 10^9 since the time units are nanoseconds
//
//    f[0] /= 1e9;
//    f[1] /= 1e9;
//    f[2] /= 1e9;
//    f[3] /= 1e9;
*/
    f[0]=0;
    f[1]=0;
    f[2]=0;
    f[3]=0;
    /* YOU GUESSED IT.... DON'T WANT ANY OF THIS
    // Doppler shift term: relative phase between particle and microwave changes with Z according to the equation dPhi/dz = c/wda.   So dPhi/dt = dPhi/dz * dz/dt, and f[2] is dz/dt in meters per nanosecond.  Sign is arbitrary; it depends which way the microwaves are propagating.  In principle they could propagate in any direction, and you could have a Dopplar shift depending on dx/dt or whatever.    
//    f[5] += wda(t) / c * f[2];  
//
//    if (debug==60) {
//        cout << t << " " << f[4] << " " << f[5] << " via " << wda(t) <<
//            " - " << w0 / (me + y[4]) << " " <<  y[2] << " " << f[2] << endl;
//    }
//    if (debug==61) {
//        cout << t << " " <<  qesi << " " <<  ppar << " " <<  bz << " " <<
//           btot << " " <<  f[2] << " " <<  y[2] << " " <<  endl;
//    }
    */

    return GSL_SUCCESS; 
} 

int foo (double t1, double ene, ofstream& filename)
{
    // Dammit, I'm just blindly applying Ben's stuff here... it would be better if I knew what was going on
    filename << "something simple" << endl;
    double mu = 10;
    double time = 0;//#Double_t
    double phistep = 0.1;//Double_t
    double txinit = 0.1;//Double_t
    double tyinit = 0.1;//Double_t
    double tzinit = 0.1;//Double_t
    double tpparinit = 0.1;//Double_t
    //double ene = 18001; // is in [eV]
    double dphi = 0;//Double_t
    //double t1 = 4e5;//100.0;//ns
    double h = phistep;
    double y[6] = {txinit, tyinit, tzinit, tpparinit, ene, dphi};

    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step * stepper = gsl_odeiv_step_alloc(T, 6);
    gsl_odeiv_control * controller = gsl_odeiv_control_y_new(1e-8, 0.0);
    gsl_odeiv_evolve * evolver = gsl_odeiv_evolve_alloc(6);
    gsl_odeiv_system sys = {func, NULL, 6, &mu};

    int status = GSL_SUCCESS;
    while (time < t1) {
        status = gsl_odeiv_evolve_apply (evolver, controller, stepper, &sys, &time, t1, &h, y);
        cout << "time is " << time << "     " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] << " " << y[6] << endl;
    }

    return 0;

}

int main(int argc, char* argv[]) {
    //int foop = foo();
   // outputfile << "something 3" << endl;
    if (argc != 3) {
        cout << "usage: $ ./aspen <final_time> <initial_energy>" << endl;
        return 1;
    } else {
        double time_f= atof(argv[1]);
        double energy_i= atof(argv[2]);
        stringstream filename;
        filename << "timeF" << time_f << "energyI" << energy_i << ".txt";
        ifstream dum(filename.str().c_str());
        cout << "filename is: " << filename.str() << endl;
        if (dum) {
            cout << "file already there" << endl;
            dum.close();
            return 1;
        }
        ofstream outputfile;
        outputfile.open(filename.str().c_str());
        //outputfile.open("footext.txt");
        if (outputfile != NULL) {
            int foop = foo(time_f, energy_i, outputfile);
            printf("some stuff in the file\n");
        } else {
            printf("file is NULL\n");
        }
        //int foop = foo(time_f, energy_i, outputfile);
        printf("I don't do anything yet!\nThanks for playing.\n");
        outputfile.close();
        return 0;
    }
}
