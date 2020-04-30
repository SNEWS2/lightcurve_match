#include "conversionfactor.hh"
using namespace std;

TGraph *gsigmatab = NULL;

Double_t funcsigma(Double_t Enu){
	if (gsigmatab == NULL) {
		gsigmatab = new TGraph();
		gsigmatab->SetName("gsigmatab");
		double _Enu, sigmatab;
		ifstream ifs("cross_Vissani.dat");
		while (ifs >> _Enu >> sigmatab) {
			//sigmatab*=1e-41; //the data in the table is in the 1e-41 cm^2
			gsigmatab->SetPoint(gsigmatab->GetN(),_Enu,sigmatab);
		}
	}
	if (Enu<gsigmatab->GetX()[0]) return 0;
	return gsigmatab->Eval(Enu, 0, "S");
}

Double_t funcnudistr(Double_t Enu, Double_t norm, Double_t alpha, Double_t avEnu){
	/* 
	 * parametrisation according to 
	 * arXiv:1211.3920v 2High-resolution supernova neutrino spectra represented by a simple fit and
	 * JUNO TDR 1507.05613
	 */
	return norm*pow(Enu,alpha)*exp(-(1+alpha)*Enu/avEnu);
}

//distribution of the neutrino energies in the flux
Double_t funcnudistrwrapper(Double_t *x, Double_t *par){
        double norm = par[0];
        double alpha = par[1];
        double avEnu = par[2];
        double Enu = x[0];
        return funcnudistr(Enu,norm,alpha,avEnu);
}

//distribution of the neutrino interactions in the flux "f x sigma_ibd
Double_t funcnuintdistrwrapper(Double_t *x, Double_t *par){ 
    double norm = par[0];
    double alpha = par[1];
    double avEnu = par[2];
    double Enu = x[0];
    return funcnudistr(Enu,norm,alpha,avEnu)*funcsigma(Enu);
}

/* 
* norm - should come from the condition that Integral from funcnudistr(Enu,norm,alpha,avEnu) is 1 (so it is a probability distribution dNnu/dEnu)
*/
Double_t funcposdistr(Double_t Epos, Double_t norm, Double_t alpha, Double_t avEnu){
	double Enu = Epos + 1.293 + 0.511;  //neglecting the nucleon recoil (and Lorentz factor for c.m); Epos is the positron kinetic energy!
	return funcnudistr(Enu,norm,alpha,avEnu)*funcsigma(Enu); 
}

Double_t funcposdistrwrapper(Double_t *x, Double_t *par){
        double norm = par[0];
        double alpha = par[1];
        double avEnu = par[2];
        double Epos = x[0];
        return funcposdistr(Epos,norm,alpha,avEnu);
}

double conversionfactor(){
	ostringstream oss, oss2;
	//1 joule is equal to 10000000 erg (1e7erg). 1 MeV is 1e6*1.6e-19 J = 1.6e-13 J = 1.6e-6 erg
	double MeVtoerg = 1.6021773e-6;
	double kpctocm = 3.086e+21;
	double Etot = 1e51/MeVtoerg; // we want to estimate number of interactions per 1e51 erg
	double dist = 10*kpctocm; 
	double Ntarget = 1e9/18*6.022e23*2; //number of free protons in the target (1kt) 1 kt = 1e9 gr, H20 - 18 gr/mol, Na = 6.022e23 1/mol, 2 H per molecule H20 
        double alpha = 3; //value 3 suggested in JUNO, alpha 2-5 in "MONTE CARLO STUDY OF SUPERNOVA NEUTRINO SPECTRA FORMATION"
	double avEnu = 14; //MeV, average value used in JUNO's study

        TF1 *fdistr = new TF1("fdistr",funcnudistrwrapper,0,100,3);
        fdistr->SetTitle(";E [MeV]; dN_{#tilde{v_{e}}}/dE");
        fdistr->SetNpx(1000);
        fdistr->SetParameters(1,alpha,avEnu);
        double norm = fdistr->Integral(0,avEnu*5);

        TF1 *fnuintdistr = new TF1("fnuintdistr",funcnuintdistrwrapper,0,100,3);
        fnuintdistr->SetTitle(";E [MeV]; dN_{e^{+}}/dE*#sigma_{IBD}");
        /*
        * 1/norm to normalis prob distr for neutrinos
        * Etot/avEnu -- number of nbe in 100ms at the source (electron antineutrinos)
        * 1/(4.*TMath::Pi()*dist*dist)   -- number of neutrinos at the detector (per surface of the detector)
        * Ntarget - number of targets (free protons)
        * 1e-41 crossection is written in 1e-41 cm^2 units
        */
        double nuintnorm = 1./norm*(Etot/avEnu)/(4.*TMath::Pi()*dist*dist)*Ntarget*1e-41;
        fnuintdistr->SetParameters(nuintnorm,alpha,avEnu);
        fdistr->SetNpx(1000);
        //double totnuint = fnuintdistr->Integral(0,100,1.e-10 ); //ROOT6 
        double totnuint = fnuintdistr->Integral(0,100);
        cout << "<Enu>= " << avEnu << " MeV, total number of interactions in 1 kt per 1e51 erg luminosity at 10 kpc: " << totnuint << endl;
        return totnuint;



}
