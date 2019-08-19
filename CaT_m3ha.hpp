// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// T-type calcium conductance
//
// Adapted from:
//      Dendritic Low-Threshold Calcium Currents in Thalamic Relay Cells
//      Alain Destexhe, Mike Neubig, Daniel Ulrich and John Huguenard
//      https://www.jneurosci.org/content/18/10/3574.short

#ifndef CAT_M3HA
#define CAT_M3HA
#include "conductance.hpp"

// Define constants (compatible with NEURON)
//  https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3783
const double FARADAY = 96485.309;
const double R = 8.3134;

// Inherit conductance class specifications
//  Note: Class names must not be all caps
class CaT_m3ha: public conductance {

private:

    // zzFF_by_RT = 4*F*F/(R*T)
    double zzFF_by_RT;

    // K_exp = exp(-ZFV/RT)
    double K_exp;

    // For temperature dependence
    double qm;
    double qh;
    double delta_temp;
    double phim;
    double phih;

public:

    // Note: the following are declared in conductance.hpp
    //      container, gbar, gbar_next, E, m, h, minf, hinf, taum, tauh, p, q, 
    //      dt, temperature, is_calcium

    // Maximum permeability in units of cm/s
    double pbar; 

    // Shifts in activation and inactivation curves
    double shiftm;
    double shifth;
    double slopem;
    double slopeh;

    // Specify parameters + Initial conditions
    CaT_m3ha(double g_, double E_, double m_, double h_, double pbar_, double shiftm_, double shifth_, double slopem_, double slopeh_)
    {
        // Read in arguments
        gbar = g_;
        E = E_;
        m = m_;
        h = h_;
        pbar = pbar_;
        shiftm = shiftm_;
        shifth = shifth_;
        slopem = slopem_;
        slopeh = slopeh_;

        // Note: gbar and E are unused in this model
        //  For safety, E must be set to zero if I = g(V-E) is accidentally used
        E = 0;

        // We can't set as NaN?
        gbar = 0;

        // Set defaults for arguments
        if (isnan(pbar)) { pbar = .2e-3; }  // maximum Ca++ permeability
        if (isnan(shiftm)) { shiftm = 1; }      // depolarizing shift of activation curve
        if (isnan(shifth)) { shifth = 1; }      // depolarizing shift of inactivation curve
        if (isnan(slopem)) { slopem = 1; }      // scaling factor for slope of activation curve
        if (isnan(slopeh)) { slopeh = 1; }      // scaling factor for slope of inactivation curve

        // Inaccessible parameters
        qm = 3.6;      // Q10 for activation, based on Coulter et al 1989
        qh = 2.8;      // Q10 for inactivation, based on Coulter et al 1989

        // specify exponents of m and h
        p = 2;
        q = 1;

        // Whether to allow this channel to be approximated
        approx_m = 0;
        approx_h = 0;

        // This is a calcium channel
        is_calcium = true;
    }

    // Method declarations
    string getClass(void);
    void connect(compartment *);
    void integrate(double, double);
    double GHK_by_Vp(double, double);
    double getCurrent(double);
    double m_inf(double, double);
    double h_inf(double, double);
    double tau_m(double, double);
    double tau_h(double, double);

};

string CaT_m3ha::getClass(){
    return "CaT_m3ha";
}

void CaT_m3ha::connect(compartment * pcomp_) {
    // Set the default gbar, m, h
    conductance::connect(pcomp_);

    // Compute the constant zzFF_by_RT = 4*F*F/(R*T)
    //  Note: units of [mol / umol] * [C / mol)] * [J/(V * mol)] / 
    //                      ([J/(K * mol)] * [K])
    //              = [C/(V * umol)]
    zzFF_by_RT = (1.0e6) * fast_pow(2, 2) * FARADAY * FARADAY / 
                    (R * (temperature + 273.16));

    // Compute temperature dependence
    delta_temp = ( temperature - 23 ) / 10;
    phim = pow(qm, delta_temp);
    phih = pow(qh, delta_temp);
}

double CaT_m3ha::getCurrent(double V) {
    return g * V;
}

void CaT_m3ha::integrate(double V, double Ca) {
    // Compute the new T channel maximal conductance
    // Note: Re-write the GHK equation in the form
    //          I = gV - gE
    //      where E = 0.
    //      units of [cm / s] * [C / (V * L)] * [10 mm / cm]
    //                  = [S * mm / L] * [L / 10^6 mm^3]
    //                  = [uS / mm^2]
    gbar = pbar * GHK_by_Vp(V, Ca) * (10.0);

    // Compute the activation gating variable at this time step
    minf = m_inf(V, Ca);
    taum = tau_m(V, Ca);
    m = minf + (m - minf) * exp(-dt / taum);

    // Compute the inactivation gating variable at this time step
    hinf = h_inf(V, Ca);
    tauh = tau_h(V, Ca);
    h = hinf + (h - hinf) * exp(-dt / tauh);

    // Compute the new T channel conductance
    // Note: this is used to integrate voltage
    g = gbar * m * m * h;

    // Add to the calcium current
    // Note: this updated the calcium concentration
    container->i_Ca += getCurrent(V);
}

double CaT_m3ha::GHK_by_Vp(double V, double Ca) {
    // Compute dimensionless constant K_exp = exp(-ZFV/RT)
    K_exp = exp(-V / (container->RT_by_nF));

    // Note: units of [C/(V * umol)] * [umol / L] = [C/(V * L)]
    return zzFF_by_RT * (Ca - container->Ca_out * K_exp) / (1 - K_exp);
}

double CaT_m3ha::m_inf(double V, double Ca)
{
    return 1.0 / ( 1.0 + exp( (V - shiftm + 57.0) / -6.2 * slopem) );
}

double CaT_m3ha::h_inf(double V, double Ca) 
{
    return 1.0 / (1.0 + exp( ( V - shifth + 81.0) / 4.0 * slopeh) ); 
}

double CaT_m3ha::tau_m(double V, double Ca) 
{
    return ( 0.612 + 1.0 / 
                ( exp( (V + 132.0 - shiftm) / -16.7) 
                + exp( (V + 16.8 - shiftm) / 18.2) ) ) 
            / (phim * slopem);

}

double CaT_m3ha::tau_h(double V, double Ca) 
{
    if ( V < -80 + shifth ) {
        return 1.0 * exp((V + 467 - shifth) / 66.6)
               / (phih * slopeh);
    } else {
        return ( 28 + 1.0 * exp( (V + 22 - shifth) / -10.5) ) 
               / (phih * slopeh);
    }
}

#endif

/*
OLD CODE:

m = m_inf(V);
h = h_inf(V);

K = (4.6469e4) / (temperature + 273.16);

// Dereference the calcium concentrations
cai = container->Ca_in;
cao = container->Ca_out;

ica = pbar * m * m * h * ghk(v, cai, cao);

double efun(double);
double ghk(double, double, double);

double CaT_m3ha::efun(double z) {
    if (fabs(z) < 1e-4) {
        // 1st order Taylor approximation of z/(exp(z) - 1)
        return 1 - z/2;
    } else {
        return z / (exp(z) - 1);
    }
}

double CaT_m3ha::ghk(double V, double cai, double cao) {
    // This is ZFV/RT, which is dimensionless 
    //  after applying conversion factor 1e-3 V/mV
    z = (1e-3) * 2.0 * FARADAY * v / (R * (temperature + 273.15 ) );

    // This has units of [mM] = [umol/cm3]
    eco = cao * efun(z);

    // This has units of [mM] = [umol/cm3]
    eci = cai * efun(-z);

    // This has units of [mC/cm3]
    //  after applying conversion factor 1e-3 mmol/umol
    return (1e-3) * 2.0 * FARADAY * (eci - eco);
}

// What for?
gbar = gbar_next;

*/