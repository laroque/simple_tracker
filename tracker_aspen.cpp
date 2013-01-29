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
       //ROOT in particular
#include "TObject.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
// Internal

// Namespace(s)
using namespace std;

/*
Time to declare tons of things at global scope. Yay for != best practices
*/
//double c = 0.299792458; // speed of light in meters per ns 
double w0=27.9925;  //nrel cyclotron frequency in GHz
double me = 510998.902; // electron mass in eV
//double mesi = 9.10838188e-31; // electron mass in kg
//double qesi = 1.60217646e-19; // electron charge in Coulombs; 
//double qe = 1; // electron charge has been folded into "B" units, I think
//double dept = 0; //common-block e field strength (set in code)
//double sweept[500],sweepf[500],sweepd[500]; //sweep parameters
//double nstoevoc = 5.34428542119e-28; //momentum: convert ev/c to J*s 
bool eperturb = false;
double eperturb_q = 0 ;
bool bperturb = false;
double bperturb_b = 0 ;
//int useconstamp=0;
//int debug=0; //debug level flag
//int nsweep;
//int iswf;
//int iswd;
double econst = 1000;
double wda = w0;//eV,ns,GHz
double power0 = 1e-4;//check this for units

// huzzah for root
TFile* tfout;
TNtuple* ntout;

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

int func (double t, const double y[], double f[], void *params)
{  
    //mostly decoupled: synchrotron oscillations 
    //current version is copied exactly from tracker07/tracker_rk4.C 
    //this will need to be coupled in better once we've figured out what ppar means exactly.
    //especially: when we allow w0 to vary spatially (inhomogenous B field) 
   
    double oogamma = me / (me + y[4]);
    double beta = sqrt(1 - oogamma * oogamma);
    /*f[4] = power0 * sqrt(3 * y[4] * me) * sin(y[5]) * 0.3 -
         3.201e-9 * w0 * oogamma * w0 * oogamma * beta * beta / (1 - beta * beta); //in eV/nanosecond*/
    f[4] = power0 * sqrt(2 * y[4] * me) * sin(y[5]) * 0.3 -
        3.201e-9 * w0 * oogamma * w0 * oogamma * beta * beta / (1 - beta * beta); //is in ev/nanosecond

    f[5] = wda - w0 * oogamma; //in rad/nanosecond
    //f[5] = wda(t) - w0 * oogamma; //in rad/nanosecond

    f[0]=0;
    f[1]=0;
    f[2]=0;
    f[3]=0;

    return GSL_SUCCESS; 
} 

int RunStepper(double t1, double ene, double dphi)
{
    // Dammit, I'm just blindly applying Ben's stuff here... it would be better if I knew what was going on
    double mu = 10;
    double time = 0;//#Double_t
    double phistep = 0.01;//Double_t
    double txinit = 0.1;//Double_t
    double tyinit = 0.1;//Double_t
    double tzinit = 0.1;//Double_t
    double tpparinit = 0.1;//Double_t
    //double ene = 18001; // is in [eV]
    //double dphi = 0;//Double_t
    //double t1 = 4e5;//100.0;//ns
    double h = phistep;
    double y[6] = {txinit, tyinit, tzinit, tpparinit, ene, dphi};

    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step * stepper = gsl_odeiv_step_alloc(T, 6);
    gsl_odeiv_control * controller = gsl_odeiv_control_y_new(1e-8, 0.0);
    gsl_odeiv_evolve * evolver = gsl_odeiv_evolve_alloc(6);
    gsl_odeiv_system sys = {func, NULL, 6, &mu};

    int status = GSL_SUCCESS;
    double wrapphase;
    while (time < t1) {
        status = gsl_odeiv_evolve_apply (evolver, controller, stepper, &sys, &time, t1, &h, y);
        h = phistep; //evolve_apply updates the recommended step size, set it back so that it isn't too big
        wrapphase = y[5] - 2 * TMath::Pi() * int(y[5] / (2 * TMath::Pi()));
        //filename << time << " " << y[4] << " " << y[5] << endl;
        ntout->Fill(0,time,y[4],y[5],wrapphase);
    }
    return 0;

}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cout << "usage: $ ./aspen <final_time> <initial_energy> <initial_phase>" << endl;
        return 1;
    } else {
        double time_f= atof(argv[1]);
        double energy_i = atof(argv[2]);
        double phase_i = atof(argv[3]);
        stringstream filename;
        filename << "timeF" << time_f << "energyI" << energy_i << ".txt";
//        ifstream dum(filename.str().c_str());
//        if (dum) {
//            cout << "Finished this energy previously, not repeating: " << filename.str() << endl;
//            dum.close();
//            return 1;
//        }
//        ofstream outputfile;
//        outputfile.open(filename.str().c_str());
        tfout = new TFile(filename.str().c_str(), "recreate");

        ntout = new TNtuple("nt", "nt", "i:t:ene:ph:rph");
        //if (outputfile != NULL) {
        if (tfout != NULL) {
            wda = w0*(18000.0) / (18000.0 + me);
            RunStepper(time_f, energy_i, phase_i);
        } else {
            printf("file is NULL\n");
            //outputfile.close();
            return 1;
        }
        cout << "Finished writing file: " << filename.str() << " " << ntout->GetEntries() << endl;
        //outputfile.close();
        ntout->Write();
        tfout->Close();
        return 0;
    }
}
