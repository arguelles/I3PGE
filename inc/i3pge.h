#ifndef _ICECUBE_H_
#define _ICECUBE_H_

#include <cmath>

#include "gsl/gsl_math.h"
#include <gsl/gsl_rng.h>
#include "gsl/gsl_randist.h"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>

#include <SQuIDS/const.h>

#include "flux.h"
#include "xs.h"
#include "muon_energy_loss.h"
#include "effective_area.h"

#define _IC_AEFF_COSTHBINS_ 11
#define _IC_AEFF_EBINS_ 50

//#define _NAIVE_INTEGRATOR_
//#define _NUM_RANDOM_INT_ 500

#define _VEGAS_INTEGRATOR_
#define _VEGAS_CALLS_ 100

//#define IC_DEBUG_DOIT
//#define IC_DEBUG_AEFF

namespace I3PGE {

double KernelHelper(double *, size_t, void *);

class IceCubeEstimator {
    public:
      enum Sampling { Asimov, Poisson, Gaussian };
    private:
      /// neutrino
      const unsigned int neutrino = 0;
      /// neutrino
      const unsigned int antineutrino = 1;
      /// neutrino type auxiliary variable
      unsigned int NEUTYPE;
      /// SQuIDS units module
      squids::Const units;
      /// Cross section calculator
      std::shared_ptr<I3PGECrossSection> XS;
      /// Effective area calculator
      std::shared_ptr<MuonEffectiveAreaSpline> EA;
      /// Muon Energy Loss
      std::shared_ptr<MuonEnergyLossSpline> MELS;
      /// Flux pointer
      std::shared_ptr<Flux> FLUX = nullptr;
      /// lifetime
      const double T_IC;
      /// bin-bin correlated systematic error
      const double f_sys;
      /// solid angle
      double Omega;
      /// units of the flux
      double flux_units;
      /// energy resolution
      double sigmaE;
      /// target nucleon density
      double NT;
      /// helper gaussian function for resolution integration
      double Normal(double Er, double Et){
          sigmaE = Et;

          double norm = 1.0/(sqrt(2.0*M_PI)*sigmaE);
          double z = (Er-Et)/sigmaE;
          return norm*exp(-z*z/2.0);
      };
      /*
        /// energy - zenith support array
        int cth_size;
        int e_size;
        vector<double> cth_range;
        vector<double> loge_range;
        vector<double> DeltaCosTh;
        vector<double> DeltaE;
      */
    public:
      /// Basic constructor
      IceCubeEstimator(std::shared_ptr<MuonEnergyLossSpline> MELS,std::shared_ptr<MuonEffectiveAreaSpline> EA, std::shared_ptr<I3PGECrossSection> XS,
                       double T_IC,double f_sys = 0.):
      MELS(MELS),EA(EA),XS(XS),T_IC(T_IC),f_sys(f_sys)
      {
          Omega = 2.0*M_PI;
          NT = 1.0;
          flux_units = pow(units.GeV,-1.)*pow(units.cm,-2.)*pow(units.sec,-1.);
      }

      /// Basic Function used to calculate event spectation
      nusquids::marray<double,2>  DoIt(std::pair<std::vector<double>,std::vector<double>> zenith_energy_bins,std::shared_ptr<Flux> FLUX,Sampling mode,
                                       double N0 = double(1.0), double dgamma = double(0.0)) {
          nusquids::marray<double,2> ResultTable{zenith_energy_bins.first.size()-1,zenith_energy_bins.second.size()-1};

          for(unsigned int ci = 0; ci < ResultTable.extent(0)-1; ci++){
              for(unsigned int ei = 0; ei < ResultTable.extent(1)-1; ei++){
                  double costh_edge = zenith_energy_bins.first[ci];
                  double energy_edge = zenith_energy_bins.second[ei];
                  double event_num,event_num_error;

                  double DeltaCT = (zenith_energy_bins.first[ci+1] - zenith_energy_bins.first[ci]);
                  double DeltaE = (zenith_energy_bins.second[ei+1] - zenith_energy_bins.second[ei]);

                  #ifdef _NAIVE_INTEGRATOR_
                  // aqui viene la integracion =MAGIA=

                  const gsl_rng_type * T;
                  gsl_rng * r;
                  T = gsl_rng_default;
                  r = gsl_rng_alloc (T);

                  event_num = 0.0;
                  for(unsigned int ir = 0; ir < _NUM_RANDOM_INT_; ir++){

                      double Er = Er_edge + DeltaE * gsl_rng_uniform(r);
                      double Enu = Er_edge + 3.0 * DeltaE * (gsl_rng_uniform (r) - 0.5);
                      double costh = costh_edge + DeltaCT*gsl_rng_uniform (r);

                      event_num += Kernel(0,costh,enu,
                                    muon_initial_energy,muon_final_energy,
                                    muon_length);
                      event_num += Kernel(1,costh,enu,
                                    muon_initial_energy,muon_final_energy,
                                    muon_length);

                  }

                  gsl_rng_free(r);

                  // normalizamos
                  double DeltaEr = DeltaE;
                  double DeltaEE = DeltaE;
                  event_num = event_num*T_IC*Omega*DeltaEr*DeltaEE*DeltaCT*flux_units/_NUM_RANDOM_INT_;
                  #endif

                  #ifdef _VEGAS_INTEGRATOR_
                  double xl[5] = {
                                   costh_edge,
                                   log10(energy_edge - 3*DeltaE),
                                   log10(energy_edge - 3*DeltaE),
                                   log10(energy_edge - DeltaE),
                                   0.
                                 };

                  double xu[5] = {
                                   costh_edge + DeltaCT,
                                   log10(energy_edge + 3*DeltaE),
                                   log10(energy_edge + 3*DeltaE),
                                   log10(energy_edge + DeltaE),
                                   2.0*units.km
                                 };

                  const gsl_rng_type *T;
                  gsl_rng *r;

                  gsl_monte_function K = { &KernelHelper, 3, this};

                  size_t calls = _VEGAS_CALLS_;

                  gsl_rng_env_setup ();

                  T = gsl_rng_default;
                  r = gsl_rng_alloc (T);

                  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

                    gsl_monte_vegas_integrate (&K, xl, xu, 3, 500, r, s,
                                                   &event_num, &event_num_error);

                  do
                    {
                    gsl_monte_vegas_integrate (&K, xl, xu, 3, calls, r, s,
                                   &event_num, &event_num_error);
                                if (gsl_monte_vegas_chisq(s) == 0)
                                    break;
                    }
                  while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.15);

                  gsl_monte_vegas_free(s);
                  gsl_rng_free(r);
                  #endif

                  double Emean_GeV = 1827.03*units.GeV;
                  ResultTable[ci][ei] = N0*event_num*pow(energy_edge/Emean_GeV,dgamma);
              }
          }

          ///================ now we return what ever the user wants ===================///

          if (mode == Asimov){
            // returns expected event number per bin
            return ResultTable;
          } else if (mode == 1 or mode == 2){
            // returns poisson realizations of the event number per bin
            nusquids::marray<double,2> RealizationTable{ResultTable.extent(0),ResultTable.extent(1)};

            const gsl_rng_type * T;
            gsl_rng * r;
            T = gsl_rng_default;
            r = gsl_rng_alloc (T);

            for(unsigned int i = 0; i < ResultTable.extent(0); i++){
              for(unsigned int j = 0; j < ResultTable.extent(1); j++){
                if( ResultTable[i][j] == 0.0 )
                  RealizationTable[i][j] = 0.0;
                else {
                  if (mode == Poisson) // poisson sampling
                    RealizationTable[i][j] = static_cast<double>(gsl_ran_poisson(r,ResultTable[i][j]));
                  else if (mode == Gaussian){ // gaussian sampling
                    double sigma = sqrt(ResultTable[i][j]) + f_sys*ResultTable[i][j];
                    RealizationTable[i][j] = ResultTable[i][j] + static_cast<double>(gsl_ran_gaussian(r,sigma));
                  }
                }
              }
            }

            gsl_rng_free(r);

            return RealizationTable;
          } else {
              std::cout << "Error : wrong IceCube estimator mode" << std::endl;
              exit(0);
          }
      }

      /// Sets the unit of the input flux
      void Set_FluxUnits(double flux_units_){
        flux_units = flux_units_;
      }
      /// Sets the flux
      void Set_Flux(std::shared_ptr<Flux> FLUX_){
        FLUX = FLUX_;
      }

      void Set_Neutype(unsigned int neutype){
        NEUTYPE = neutype;
      }
    protected:
      double Kernel(unsigned int neutype, double costh, double enu,
                    double muon_initial_energy, double muon_final_energy,
                    double muon_length){
        return T_IC*(*FLUX)(neutype,costh,enu)*(*XS)(neutype,enu,muon_initial_energy)*NT*(*MELS)(muon_initial_energy,muon_final_energy,muon_length)*(*EA)(muon_final_energy);
      }
    public:
      double Kernel(double *x){
          double costh = x[0];
          double enu = pow(10.0,x[1]);
          double muon_initial_energy = pow(10.0,x[2]);
          double muon_final_energy = pow(10.0,x[3]);
          double muon_length = x[4];

          return Kernel(NEUTYPE,costh,enu,muon_initial_energy,muon_final_energy,muon_length);
      }
};

double KernelHelper(double *x, size_t dim, void *params){
    IceCubeEstimator* ICE = (IceCubeEstimator*) params;
    return ICE->Kernel(x);
};

} // close namespace

#endif
