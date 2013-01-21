// 3D simulation of a particle in a Penning trap, including microwave excitation, in the GUIDING CENTER approximation.  This ought to allow us to include various asymmetries and magnetron/axial coupling.

// 3ds4: features added 7/3/08
// 1) Change from non-adaptive to adaptive integrator.   Force the use of a point on either side of each ramp transition.   

// 3ds5: features added 7/10/08
// 1) add "perturbations" to E and B fields.  
// E field perturbation is a point charge displaced along Z.
// B field perturbation is a small field divergence.  

// 3ds6: features added 7/14/08
// try to add relativistic aspects of the Penning trap motion.  In a nutshell: 
// the relativistic effect comes only from a mass shift due to the cyclotron 
// motion.  The axial motion is always nonrelativistic.  

// The equations are from J. Zhang et. al.,, J Phys B 40 1019-1033.
// dr/dt = E x B / (B^2) + Bhat * p/M
// dp/dt = qE dot Bhat
// note: E, B, r are vectors; p is a scalar "p parallel, meaning along the axis"; Bhat is a vector length 1. 


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "TObject.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TArc.h"
#include "TLine.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TApplication.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TMarker.h"
#include "TText.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TVector3.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
using namespace std;

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
    evec[0] = -econst*y[0]; //hyperbolic electrostatic potential = proportional to R
    evec[1] = -econst*y[1]; //repulsive in xy
    evec[2] = econst*y[2]; // attractive in z  

    if (eperturb) // add a perturbing potential as due to a point charge 10cm away
      {        	
	double qor3 = eperturb_q/pow(y[0]*y[0] + y[1]*y[1] + pow(y[2]+0.1,2),1.5);  
	evec[0] = y[0]*qor3;
	evec[1] = y[1]*qor3;
	//static Z component guarantees that oscillation center is in middle
	evec[2] = (y[2]+0.1)*qor3 - eperturb_q/pow(0.1,2); 
      }

}

void get_b_field(const double y[], double bvec[])
{
    // magnetic field should be in Tesla 
    bvec[0] = 0;  
    bvec[1] = 0; 
    bvec[2] = w0/177; // w0 is in GHz, 177 = frequency when B=1T.  
    //also return total field strength
    bvec[3] = bvec[2];

    if (bperturb)
      {
	bvec[0] += bperturb_b*y[0]/2; //x, y components get stronger off axis
	bvec[1] += bperturb_b*y[1]/2;
	bvec[2] -= bperturb_b*y[2];  // z component gets weaker at +z
	bvec[3] = sqrt(bvec[0]*bvec[0]+bvec[1]*bvec[1]+bvec[2]*bvec[2]);
      }
}


//calculates the drive frequency at the present time 
double wda(double t)
{
  //  cout << "wda " << nsweep << " " << me << endl;
  for (int i=0;i<nsweep-1;i++)
    {
      if (t >= sweept[i] && t < sweept[i+1])
	{ //return an interpolated ramp 
	  double tmp =  sweepf[i] + (t - sweept[i])*(sweepf[i+1]-sweepf[i])/(sweept[i+1]-sweept[i]); 
	  return w0*me/(me+tmp);
	}
    }
  return w0*me/(me+sweepf[nsweep-1]); 
}

//calculates the electron cyclotron radius at the present time 
double rda(double t)
{
  //  cout << "wda " << nsweep << " " << me << endl;
  for (int i=0;i<nsweep-1;i++)
    {
      if (t >= sweept[i] && t < sweept[i+1])
	{ //return an interpolated ramp 
	  double tmp =  sweepf[i] + (t - sweept[i])*(sweepf[i+1]-sweepf[i])/(sweept[i+1]-sweept[i]); 
	  return sqrt(tmp);
	}
    }
  return sqrt(sweepf[nsweep-1]); 
}

//calculates the drive frequency at the present time 
double power(double t)
{
  //  cout << nsweep;
  for (int i=0;i<nsweep-1;i++)
    {
      //      cout << "." ;
      if (t >= sweept[i] && t < sweept[i+1])
	{ //return an interpolated ramp 
	  //	  cout << "power calc: " << i << " " << t << " " << sweepd[i] + (t - sweept[i])*(sweepd[i+1]-sweepd[i])/(sweept[i+1]-sweept[i]) << endl;
	  if (useconstamp != 0) // set ramp amplitude to correct for drive freq.
	    {
	      return (sweepd[i] + (t - sweept[i])*(sweepd[i+1]-sweepd[i])/(sweept[i+1]-sweept[i]))*31.6227766/rda(t); //normalized to give unmodified power at E=1000 
	    }
	  else 
	    {
	      return sweepd[i] + (t - sweept[i])*(sweepd[i+1]-sweepd[i])/(sweept[i+1]-sweept[i]); 
	    }
	}
    }
  //  cout << " return last " << sweepd[nsweep-1] << endl; 

  return sweepd[nsweep-1]; 
}


int func (double t, const double y[], double f[],
      void *params)
{  

   //mostly decoupled: synchrotron oscillations 
   //current version is copied exactly from tracker07/tracker_rk4.C 
   //this will need to be coupled in better once we've figured out what ppar means exactly.
   //especially: when we allow w0 to vary spatially (inhomogenous B field) 
   
  double oogamma = me/(me+y[4]);
  double beta=sqrt(1-oogamma*oogamma);
  f[4] = power(t)*sqrt(2*y[4]*me)*sin(y[5])*0.3 - 3.201e-9*w0*oogamma*w0*oogamma*beta*beta/(1-beta*beta); //in eV/nanosecond
  f[5] = wda(t) - w0*oogamma; //in rad/nanosecond

  // y[0,1,2] = x, y, z
  // y[3] = ppar (momentum parallel to guiding center motion) 
  // y[4] = total particle kinetic energy
  // y[5] = synchrotron orbit phase
  // see Zhang et. al. appendix for possible improvements
  // maybe add parameters for synchrotron phase, energy
    double ex , ey, ez ; 
    double bx, by, bz ; 
    double ppar, btot, bsquared; 
    
    double bvec[4]; //3-vector; 4th component is total field strength.     
    double evec[3];

    get_e_field(y,evec); // puts LOCAL electric field into evec[3] by reference;
    get_b_field(y,bvec); // puts LOCAL magnetic field into bvec[3] by reference;

    // should be in units "volts per meter" 
    ex = evec[0];
    ey = evec[1];
    ez = evec[2];
    //should be in kg*m/s
    //oogamma provides all of the relativistic correction
    ppar = oogamma*y[3]/mesi;  //and we'll pre-divide by the SI electron mass
    // magnetic field should be in Tesla 
    bx = bvec[0];  
     by = bvec[1]; 
    bz = bvec[2]; // w0 is in GHz, 177 = frequency when B=1T.  
    btot = bvec[3];  //scalar magnetic field strength
    bsquared = bvec[3]*bvec[3];  //scalar magnetic field strength

   f[0] =                   -ey*bz/bsquared + ez*by/bsquared + ppar*bx/btot; 
   f[1] =  ex*bz/bsquared                   - ez*bx/bsquared + ppar*by/btot; 
   f[2] = -ex*by/bsquared +  ey*bx/bsquared                  + ppar*bz/btot; 
   f[3] = -qesi*ex*bx/btot - qesi*ey*by/btot - qesi*ez*bz/btot; 

   //reduce it all by a factor of 10^9 since the time units are nanoseconds

   f[0] /= 1e9;
   f[1] /= 1e9;
   f[2] /= 1e9;
   f[3] /= 1e9;


  // Doppler shift term: relative phase between particle and microwave changes with Z according to the equation dPhi/dz = c/wda.   So dPhi/dt = dPhi/dz * dz/dt, and f[2] is dz/dt in meters per nanosecond.  Sign is arbitrary; it depends which way the microwaves are propagating.  In principle they could propagate in any direction, and you could have a Dopplar shift depending on dx/dt or whatever.    
  f[5] += wda(t)/c*f[2];  

  if (debug==60) {cout << t << " " << f[4] << " " << f[5] << " via " << wda(t) << " - " << w0/(me+y[4]) << " " <<  y[2] << " " << f[2] << endl;}
  if (debug==61) {cout << t << " " <<  qesi << " " <<  ppar << " " <<  bz << " " <<  btot << " " <<  f[2] << " " <<  y[2] << " " <<  endl; } 
  
   return GSL_SUCCESS; 
} 



int main(int argc, char *argv[])
{
  
//GSL ode-solver from http://www.gnu.org/software/gsl/manual/html_node/ODE-Example-programs.html
//using the ODEs laid out on notebook pg. 27: y[0] = ene, y[1]= phi   
//cyclotron power from Johner 1987: P_tot = (1/(4 pi epsilon0) * (2*q^2*w^2)/(3c) * b_t^2/(1-b^2)

  bool xrandom = false; 
  bool yrandom = false; 
  bool zrandom = false; 
  bool pparrandom = false; 
  int i=0;
  Double_t dept = 0.1; // 0.1 ev per turn
  Double_t wd = w0*me/(me+18000);
  Double_t time=0;
  Double_t endtime=0;
  Double_t ene=18001; // in eV
  Double_t de, w, dt;
  Double_t tcrit[5], tcritb; //time in nanoseconds
  int icrit;
  Double_t dphi=0;
  Double_t beta, dedt,dwdt,dwde;
  int liveupdate=0;
  int speedup=1;
  Double_t giveup=400;
  Double_t edmin=0;
  Double_t edmax=21000;
  Double_t tdmin=0;
  Double_t fixphase=0.75;
  Double_t tdmax=500;
  Double_t eloss=1;
  Double_t xinit = 0.1; 
  Double_t yinit = 0.1; 
  Double_t zinit = 0.1; 
  Double_t pparinit = 0.1; 

  Double_t txinit = 0.1; 
  Double_t tyinit = 0.1; 
  Double_t tzinit = 0.1; 
  Double_t tpparinit = 0.1; 
  
  Double_t tmax = 500;
  //  Double_t sweept[500],sweepf[500],sweepd[500];
  Double_t thisf;
  Double_t phistep=0.1;
  Int_t phasestep=100;
  //  Int_t nsweep=0;
  Int_t isw=0;
  Double_t noiseamp=0;
  Double_t myrand=1;
  Double_t enelist[500]={15000,16000,17000,18000,19000,20000,21000};
  Int_t nphases=1;
  Double_t phaselist[50000]={0.75};
  int nenes=6;

  TString useenefile="";
  bool useenefile_flag=false;
  Int_t fentries=1; 
  Int_t enefile_repeat=1; 
  const int n_fenes_max = 10000;
  Double_t finitlist[n_fenes_max][2];
  int fdraw=5;
  int labels=1;
  int batchmode=0;
  int savetrack=1;
  int saveps=1;
  int saveroot=1;
  int identifier=-1; // reference number of directory; stored with "data"
  int bailedout=0;

  TString dummy, argname,psfilename, psnrgfilename, psdedtfilename, rootfilename;
  Double_t value;

  TApplication theApp("App",&argc,argv);
  TROOT rsession("test", "test");

  //1st command line argument is name of input card, "card.blahblahblah".
  // Notice that the name *must* begin with "card."  

  TString cardname=theApp.Argv(1);
  char cdummy[1000]; 
//  TString cardname="card.test";
  ifstream cardfile(cardname,ios::in);

  cout << "Hello! The cardfile is " << cardname << endl;

  // read input card; each line has a one-word flag and one or more numbers

  while (cardfile >> argname)
    {
      // first, take care of any case which isn't a regular "card command"
      //      cout << argname << endl;
      if (argname.Contains("#")) //skip lines beginning with #.  Note the single quotes (argname[0] is of type char, not type string) 
	{
	  //get the rest of the line and discard it 
	  cardfile.unget(); //fixes special case of a blank line with # on it
	  cardfile.ignore(1000,'\n');
	  continue; 
	}
      else if (argname=="END") {break;} // 
      
      //if the line isn't a #comment or an END, get the first "number", then
      // figure out what the flag word was and decide what to do with it
      
      cardfile >> value;      

      if (!cardfile) { 
	cout << "error in cardfile: " << argname << " not followed by valid value. ( " << value << " instead.)  Exiting." << endl;
	exit(0); } 

      //all valid arguments should appear in this list somewhere.  
      if (argname=="enelist") { //List of initial electron energies.
	nenes=value;                 // first parameter: # of quants in list
	if (nenes > 0)
	  {
	    for (int i=0;i<nenes;i++)    // other parameters:  energies in eV
	      {cardfile >> value;
		enelist[i]=value;
	      }
	  }
	else if (nenes == -1) 	 
	  {
	    Double_t enelow; //where to start
	    Double_t step; // how big the steps are
	    int numenestep; //how many steps to take
	    cardfile >> enelow;
	    cardfile >> step;
	    cardfile >> numenestep;
	    
	    cout << enelow << ", " << step << ", " << numenestep << endl;
	    
	    for (int i=0; i<numenestep; i++){
	      enelist[i]=enelow+i*step;
	      //	      cout << phaselist[i] << endl;
	    }
	    
	  nenes = numenestep;
	  }
	else if (nenes == -2)
	  {
	    int numenestep;
	    Double_t enelow, enehi;
	    
	    cardfile >> numenestep >> enelow >> enehi;	    
	    cout << numenestep << " " << enelow << " " << enehi<<endl;
	    
	    TRandom *enernd = new TRandom();
	    enernd->SetSeed(UInt_t(cardname.Hash()));		

	    for(int i=0; i<numenestep; i++){
	      enelist[i]=enernd->Uniform(enelow,enehi);
	      //	cout << phaselist[i] << endl;
	    }
	    
	    nenes=numenestep;	    
	  }
	else if (nenes == -3) //read a list of energies and phases from a file
	  {
	    cardfile >> enefile_repeat; //how many times to repeat each line
	    cardfile >> useenefile ; //get the filename
	    
	    useenefile_flag=true;
	  }

      }
      else if (argname=="phaselist") { // list of initial orbit phases
	nphases=value;                 //first parameter: # of phases in list
	if (nphases > 0){
	  for (int i=0;i<nphases;i++)    // other parameters: phases (range 0-1)
	    {cardfile >> value;
	      phaselist[i]=value;
	      //	      cout << phaselist[i] << endl;
	    }
	}
	else if (nphases == -1){//allow for a list of constructed phases
	  
	  Double_t phaselow; //where to start
	  Double_t step; // how big the steps are
	  int numphasestep; //how many steps to take
	  
	  cardfile >> phaselow;
	  cardfile >> step;
	  cardfile >> numphasestep;
	  
	  cout << phaselow << ", " << step << ", " << numphasestep << endl;
	  
	  for (int i=0; i<numphasestep; i++){
	    phaselist[i]=phaselow+i*step;
	    //	      cout << phaselist[i] << endl;
	  }
	  
	  nphases = numphasestep;
	}
	else if (nphases==-2){//allow for a list of random phases
	  
	  int numphasestep;
		
		cardfile >> numphasestep;
		cout << numphasestep << endl;
		
		TRandom *phsrnd = new TRandom();		
		phsrnd->SetSeed(UInt_t(cardname.Hash()+1));
		
		for(int i=0; i<numphasestep; i++){
		  phaselist[i]=phsrnd->Uniform(0.0,2*TMath::Pi());
		  //	cout << phaselist[i] << endl;
		}
		
		nphases=numphasestep;
	}
	else{

	   cout << "phases error! check card!" << endl;
	}
      }
      else if (argname=="useconstamp") { useconstamp=value; } //detail of "input power"-"de/dt" convers 
      else if (argname=="speedup"){  speedup=value; }   // if !=0: Abort when e ene falls "giveup" below drive
      else if (argname=="giveup"){  giveup=value; }     // how far below drive to consider "out of bucket"
      else if (argname=="w0"){  w0=value; }              // set zero-energy cyclotron freq in GHz
      else if (argname=="dept") {	  dept=value; }  // ignored
      else if (argname=="edmax"){  edmax=value; }     // max energy to display
      else if (argname=="edmin") {  edmin=value; }  // min energy to display
      else if (argname=="tdmax") {  tdmax=value; }  // max time to display (ns)
      else if (argname=="tdmin") {  tdmin=value; } // min time to display
      else if (argname=="tmax") {  tmax=value; }   // end time of run (ns)
      else if (argname=="eloss") { eloss=value;}    // include energy loss? =1 for correct cyclotron radiation
      else if (argname=="noise") { noiseamp=value;}   // amplitude of random noise (as fraction of drive) 
      else if (argname=="wde") {  wd=w0*me/(me+value); }  //ignored
      else if (argname=="debug") { debug=value; }       //set to -1 for minimal output
      else if (argname=="liveupdate") { liveupdate=value; } // if 1: refresh TCanvas periodically mid-run
      else if (argname=="phistep") {  phistep=value; }      // Tracking step size: =1 for > 1 update per beat.
      else if (argname=="phasestep") {  phasestep=value; }  // Tracking step size: =1 for > 1 update per orbit
      else if (argname=="fdraw") { fdraw=value; }  // what fraction of points to draw/save with TGraph (>5000)
      else if (argname=="saveps") { saveps=value; } // =1 to save a postscript (big!) of the canvas
      else if (argname=="savetrack") { savetrack=value; } // =1 to save a ROOT file of the whole path
      else if (argname=="saveroot") { saveroot=value; } // =1 to save a ROOT summary of init/final states
      else if (argname=="labels") { labels=value;}     // =1 to draw some text at the track endpoint
      else if (argname=="batchmode") {batchmode=value;}      // =1 to exit quietly; =0 to wait for keystroke
      //      else if (argname=="batch") {batchmode=value;}      // =1 to exit quietly; =0 to wait for keystroke
      else if (argname=="identifier") {identifier=value;} // arbitrary number saved in the SAVEROOT tree
      else if (argname=="fixphase"){fixphase=value;}     // ignored
      else if (argname=="xyzinit"){xinit=value; cardfile >> yinit; cardfile >> zinit;}     // init magnetron
      else if (argname=="xrandom"){if ( value!=0)xrandom = true;}  //user may set magnetron values to random; in this case xyzinit values are used as range
      else if (argname=="yrandom"){if (value!=0) yrandom = true;}  
      else if (argname=="zrandom"){if (value!=0)zrandom = true;}   
      else if (argname=="pparinit"){pparinit=value;}  
      else if (argname=="pparrandom"){if (value!=0) pparrandom = true;}     // init magnetron
      else if (argname=="econst") {econst = value; } 
      else if (argname=="bperturb") {bperturb_b = value; bperturb= true; } 
      else if (argname=="eperturb") {eperturb_q = value; eperturb= true;} 
      else if (argname=="econst") {econst = value; } 
      else if (argname=="postscript")	//LEGACY, do not use
      	{ // value is discarded!	  
      	  cardfile >> psfilename;
	}


      // TO tell the cyclotron how to frequency-ramp, you set a handful of reference points;
      // each point specifies a time, a frequency, and a drive power.  The actual freqency 
      // applied to the particle is a linear ramp from one point to the next. Frequencies are specified in
      // terms of the equivalent electron energy.  
      //    Example: four-point ramp (start, raise energy, lower amplitude, hold) 
      //  a) ramp from E=0 to E=1000 eV in 30 microseconds (at 5e-9 "drive power")
      //  b) hold energy at E-1000 for 1000 ns while ramping power from 5e-9 to 1e-10
      //  c) hold energy at 1000, power at 1e-10 for 200us.  
      //
      // sweepf 4 000 1000 1000 1000
      // sweept 4 000 30000 31000 231000
      // sweepd 4 5e-9 5e-9 1e-10 1e-10

      else if (argname=="sweepf")      //sweep frequency specs (in eV-equivalent)
	{nsweep=value;
	  for (int i=0;i<nsweep;i++)
	    {cardfile >> value;
	      sweepf[i]=value;
	    }}
      else if (argname=="sweept") //sweep time specs (in nsec)
	{nsweep=value;
	  for (int i=0;i<nsweep;i++)
	    {cardfile >> value;
	      sweept[i]=value;
	    }}
      else if (argname=="sweepd")  //sweep drive specs (in arbitrary units; usually 1e-8 and below.)
	{nsweep=value;
	  for (int i=0;i<nsweep;i++)
	    {cardfile >> value;
	      sweepd[i]=value;
	      if (debug==3) cout << "s " << sweepd[i] << endl; 
	    }}
      else if (argname=="rootfile")  //LEGACY, do not use	
	{ // value is discarded!	  
	  cardfile >> rootfilename;
	}
    }

  // OK, we've read the card file.   Prepare outputs, if any: 

  TCanvas* c1 = new TCanvas();
//  TPostScript *3myps;
  if (saveps!=0)
    {
      psfilename=cardname;
      psfilename.ReplaceAll("card",4,"ps",2);
      psnrgfilename = psfilename + "nrg.eps";
	  psdedtfilename= psfilename +  "dedt.eps";
      psfilename += ".eps";	
      cout << "psfile " << psfilename << endl;
//      myps = new TPostScript(psfilename);      
    }

  TFile* outfile;
  TNtuple* data;
  Float_t data_array[20];
  TNtuple* track; 
  if (saveroot!=0 || savetrack !=0)
    {
      rootfilename=cardname;
      rootfilename.ReplaceAll("card",4,"root",4);
      rootfilename+=".root";
      rootfilename.ReplaceAll(".root.root",".root");
      outfile = new TFile(rootfilename.Data(),"RECREATE");
    }

  if (saveroot!=0)
    {   data = new TNtuple("cyclo","cyclo","identifier:iene:iphase:drive1:drive2:rampt:rampf:einit:phaseinit:xinit:yinit:zinit:pparinit:efinal:phasefinal:xfinal:yfinal:zfinal:pparfinal");    }

  TH1F* dedtfft = new TH1F("name4","fft of dedt",200, -100,100);

  TMarker* mm = new TMarker();
  TMarker* mwd = new TMarker();
  mm->SetMarkerStyle(1);
  mwd->SetMarkerStyle(1);
  mwd->SetMarkerColor(2);
  
  for (int i=0;i<nsweep;i++)
    { cout << i << " " << sweept[i] << " " << sweepf[i] << " " <<sweepd[i] << endl; } 

  //if we're reading energies and phases from an init file, do it now
  if (useenefile_flag) 
    {
      ifstream tfenefile(useenefile,ios::in);
      fentries = 0 ; 
      while (tfenefile >> finitlist[fentries][0] >> finitlist[fentries][1])
	{
	  fentries++;
	  if (fentries >= n_fenes_max) break;
	}
      if (fentries==0) { cout << "no luck with " << useenefile << ";abort" <<endl; exit(0); } 
      nenes=1;
      nphases=enefile_repeat; 
    }

  //make sure RNG isn't correlated between Z and phase
  TRandom* tr = new TRandom();
  tr->SetSeed(UInt_t(cardname.Hash()+2)); 

  for (int i=0;i<nenes;i++)
    {
      for (int iphase=0;iphase<nphases;iphase++)
	{
	  for (int ifentries=0;ifentries<fentries;ifentries++)
	    {
	  if (savetrack!=0)
	    {
	      TString trackname="track";

	      if (!useenefile_flag)		{
	      trackname+=i;
	      trackname+="_";
	      trackname+=iphase;	      		}
	      else 		{
		  trackname+=ifentries; 
		  trackname+="_";
		  trackname+=iphase; 		}
	      track = new TNtuple(trackname.Data(),trackname.Data(),"ene:phase:x:y:z:ppar:time");
	    }

	  if (!useenefile_flag)
	    {
	      dphi=phaselist[iphase];
	      ene = enelist[i];
	    }
	  else
	    {
	      ene = finitlist[ifentries][0];
	      dphi = finitlist[ifentries][1];
	      cout << "init from file" << ene << " " << dphi << endl;
	    }
	  


       if (xrandom)
	 { txinit = tr->Rndm()*xinit; } 
       else {txinit = xinit; } 

       if (yrandom)
	 { tyinit = tr->Rndm()*yinit; } 
       else {tyinit = yinit; } 

       if (zrandom)
	 { tzinit = tr->Rndm()*zinit; } 
       else {tzinit = zinit; } 

       if (pparrandom)
	 { tpparinit = tr->Rndm()*pparinit; } 
       else {tpparinit = pparinit; } 
       
	  time=0;
	  Double_t lasttime=-100;     
	  if (debug != 0) {  
	    cout << "======================="; 
	    cout << "ene = " << ene << " dphi = " << dphi << endl;
	    cout << "xyzp= " << txinit << " " << tyinit << " " << tzinit << " "<< tpparinit << endl;

	  }

	  // prepare run parameters for Ntuple
	  data_array[0] = identifier; 
	  data_array[1] = i;
	  data_array[2] = iphase;
	  data_array[3] = sweepd[1];
	  data_array[4] = sweepd[2];
	  data_array[5] = sweept[2];
	  data_array[6] = sweepf[2];
	  data_array[7] = ene;
	  data_array[8] = dphi;
	  data_array[9] = txinit;
	  data_array[10] = tyinit;
	  data_array[11] = tzinit;
	  data_array[12] = tpparinit;

	  //	      data->Fill(identifier,i,iphase,sweepd[1],sweepd[2],sweept[1],sweepf[1],enelist[i],ene,phaselist[iphase]);
 const gsl_odeiv_step_type * T 
         = gsl_odeiv_step_rk8pd;
     
       gsl_odeiv_step * stepper 
         = gsl_odeiv_step_alloc (T, 6);
       gsl_odeiv_control * controller 
         = gsl_odeiv_control_y_new (1e-8, 0.0);
       gsl_odeiv_evolve * evolver 
         = gsl_odeiv_evolve_alloc (6);
     
       double mu = 10;
       gsl_odeiv_system sys = {func, NULL, 6, &mu};
       
       //       double time, tmax;
       double t = 0.0, t1 = 100.0;
       double h = phistep; 
       double y[6] = { txinit,tyinit,tzinit,tpparinit,ene,dphi}; // LATER we have to make this initialization card-controllable.  
       double y_err[6];
       double dydt_in[6], dydt_out[6];
     
       /* initialise dydt_in from system parameters */
       //       GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);

       isw=0;
       bailedout=0;
       int j=0;
       
       //perform ODE evolution with adaptive step size control.  
       // we want to force the stepper to use a point right-close-to each 
       // ramp transition, so we let those transition points define tmax.  
       // for each ramp segment.  
       for (int ir=0;ir<2*nsweep;ir++)
	 {
	   if (ir%2==0) //long ramp entirely within each gap 
	     {
	       t1 = sweept[ir/2 + 1] - 0.1; 
	     }
	   else //extra points spanning each gap
	     {
	       t1 = sweept[ir/2 + 1] + 0.1;
	     }
	   
	   while (time < t1)
	     {
	       j++;
	       int status = gsl_odeiv_evolve_apply (evolver, controller, stepper,
                                       &sys,
                                       &time, t1, &h,
                                       y);
		  
		  if (status != GSL_SUCCESS)
		    { cout << "integration error" << endl; break;}
		  
		  //		  time += h;
		  
// 		  dydt_in[0] = dydt_out[0];		  
// 		  dydt_in[1] = dydt_out[1];		  
// 		  dydt_in[2] = dydt_out[2];		  
// 		  dydt_in[3] = dydt_out[3];		  
// 		  dydt_in[4] = dydt_out[4];		  
// 		  dydt_in[5] = dydt_out[5];
		  

		  
		  if (debug == 5) printf ("gsl: %.5e %.5e %.5e\n", time, y[4], y[5]);

		  if (y[4] != y[4]) {if (debug != 0) cout  << "NaN error" << endl; bailedout=1; break;} 
		  if (y[4] <= 0) {if (debug != 0) cout  << "bailout-el!" << endl; bailedout=1; break;} 
		  if (y[4] > 1.0e9 || y[4] < -1.0e9) {if (debug != 0) cout  << "bailout-eh!" << endl; bailedout=1; break;} 

		  if (speedup != 0)
		    {
		      if (y[5] > 6.29 || y[5] < -6.29) {if (debug != 0) cout  << "bailout-p!" << endl; bailedout=1; break;}
		    }
 
		  //		  if (y[5] > 10 || y[5] < -10) {if (debug != 0) cout  << "bailout-p!" << endl; bailedout=1; break;} 

	      if (savetrack != 0)
		{
		  if (j%fdraw==0)
		    {
		      //		      ene = y[4];
		      //		      dphi = y[5];
		      //		      mm->DrawMarker(time,y[4]);
		      //		      mm->DrawMarker(y[5],y[4]);

		      track->Fill(float(y[4]),float(y[5]),float(y[0]),float(y[1]),float(y[2]),float(y[3]),float(time));
		    }
		}
	     }
		 
	 }
	  
	  if (debug==1) cout << "updateing " ;
	  if (batchmode == 0) 
	    {
	      c1->Update();      
	    }
	  endtime=time;
	  time=0;
	  if (debug==1) cout << "done " ;
	  if (bailedout==0)
	    {
	      cout << "final e,phi \t" << y[4] << "\t " << y[5] << endl;
	    }
	  if (savetrack!=0 && bailedout==0) 
	    {
	      track->Write();
	      cout << "good track saved!" << endl;
	    }
	  if (savetrack!=0)
	    {
	      track->Delete(""); //should free up memory
	    }
	  if (saveroot!=0) 
	    {
	      data_array[13] = float(y[4]);
	      data_array[14] = float(y[5]);
	      data_array[15] = float(y[0]);
	      data_array[16] = float(y[1]);
	      data_array[17] = float(y[2]);
	      data_array[18] = float(y[3]);
	      data->Fill(data_array);
	    }

	  //cin.get();

	    }	    
	}
    }

/*  TCanvas* c4 = new TCanvas("dd4","display4",400,600);
  c4->cd();
  dedth->FFT(dedtfft,"RE R2C ES");
  dedtfft->Draw("");
  c4->Update();
*/


  if (saveroot != 0)
  {
    cout << "saving"<<endl;
    data->Write();
    cout << "flushing"<<endl;
    outfile->Flush();
    cout << "done"<<endl;
  }

  
/*
  outfile->Close();
*/
  cout << "done " << endl;
  if (batchmode == 0)
    {  TMarker *mark = (TMarker*)c1->WaitPrimitive("TMarker","Marker"); } 
}
