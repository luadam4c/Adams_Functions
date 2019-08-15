// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// T-type calcium conductance
// Destexhe, Contreras, Sejnowski, & Steriade 1994
// https://www.ncbi.nlm.nih.gov/pubmed/7527077
// Soplata AE, McCarthy MM, Sherfey J, Lee S, Purdon PL, Brown EN, et al.
// Thalamocortical control of propofol phase-amplitude coupling. PLOS
// Computational Biology. 2017;13: e1005879. doi:10.1371/journal.pcbi.1005879

#ifndef CAT
#define CAT
#include "conductance.hpp"

//inherit conductance class spec
class CaT: public conductance {

public:

    // Specify parameters + Initial conditions
    CaT(double V_, double pcabar_, double shiftm_, double slopem_, double slopeh_)
    {
        // Read in arguments
        V = V_;
        pcabar = p_;
        shiftm = shiftm_;
        shifth = shifth_;
        slopem = slopem_;
        slopeh = slopeh_;

        // TODO: How to get celsius?
        // TODO: Where do we define constants FARADAY and R

        // Set defaults for accessible parameters
        if (isnan(pcabar)) { pcabar = .2e-3; }  // maximum Ca++ permeability
        if (isnan(shiftm)) { shiftm = 1; }      // depolarizing shift of activation curve
        if (isnan(shifth)) { shifth = 1; }      // depolarizing shift of inactivation curve
        if (isnan(slopem)) { slopem = 1; }      // scaling factor for slope of activation curve
        if (isnan(slopeh)) { slopeh = 1; }      // scaling factor for slope of inactivation curve

        // Inaccessible parameters
        qm = 3.6;           // Q10 for activation, based on Coulter et al 1989
        qh = 2.8;           // Q10 for inactivation, based on Coulter et al 1989

        // Parameters computed once
        phim = fast_pow(qm, ( celsius - 23 ) / 10 );
        phih = fast_pow(qh, ( celsius - 23 ) / 10 );
        m = m_inf(V);
        h = h_inf(V);

        // Allow this channel to be approximated
        // TODO: What's this?
        approx_m = 1;
        approx_h = 1;
    }

    // Method declarations
    string getClass(void);
    void integrate(double);
    double efun(double);
    double ghk(double, double, double);
    double m_inf(double);
    double h_inf(double);
    double tau_m(double);
    double tau_h(double);

};

string CaT::getClass(){
    return "CaT";
}

void CaT::integrate(double V) {
    // Dereference the calcium concentrations
    cai = container->Ca_in;
    cao = container->Ca_out;

    // Compute the activation gating variable at this time step
    m = m_inf(V) + (m - m_inf(V)) * exp(-dt / tau_m(V));

    // Compute the inactivation gating variable at this time step
    h = h_inf(V) + (h - h_inf(V)) * exp(-dt / tau_h(V));

    // Compute the new T channel conductance
    ica = pcabar * m * m * h * ghk(v, cai, cao);

    // Add to the calcium current
    container->i_Ca += ica;
}

double efun(double z) {
    if (fabs(z) < 1e-4) {
        // 1st order Taylor approximation of z/(exp(z) - 1)
        return 1 - z/2;
    } else {
        return z / (exp(z) - 1);
    }
}

double ghk(double V, double cai, double cao) {
    // This is ZFV/RT, which is dimensionless 
    //  after applying conversion factor 1e-3 V/mV
    z = (1e-3) * 2.0 * FARADAY * v / (R * (celsius + 273.15 ) );

    // This has units of [mM] = [umol/cm3]
    eco = cao * efun(z);

    // This has units of [mM] = [umol/cm3]
    eci = cai * efun(-z);

    // This has units of [mC/cm3]
    //  after applying conversion factor 1e-3 mmol/umol
    return (1e-3) * 2.0 * FARADAY * (eci - eco);
}

double CaT::m_inf(double V)
{
    return 1.0 / ( 1.0 + exp( (V - shiftm + 57.0) / -6.2 * slopem) );
}

double CaT::h_inf(double V) 
{
    return 1.0 / (1.0 + exp( ( V - shifth + 81.0) / 4.0 * slopeh) ); 
}

double CaT::tau_m(double V) 
{
    return ( 0.612 + 1.0 / 
                ( exp( (V + 132.0 - shiftm) / -16.7) 
                + exp( (V + 16.8 - shiftm) / 18.2) ) ) 
            / (phim * slopem);

}

double CaT::tau_h(double V) 
{
    if ( v < -80 + shifth) {
        return 1.0 * exp((v + 467 - shifth) / 66.6)
               / (phih * slopeh);
    } else {
        return ( 28 + 1.0 * exp( (v + 22 - shifth) / -10.5) ) 
               / (phih * slopeh);
    }
}

#endif
