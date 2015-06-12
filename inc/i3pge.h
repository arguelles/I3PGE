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
      /// SQuIDS units module
      squids::Const units;
      /// Cross section calculator
      std::shared_ptr<CrossSection> XS;
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
      IceCubeEstimator(std::shared_ptr<MuonEnergyLossSpline> MELS,std::shared_ptr<MuonEffectiveAreaSpline> EA, std::shared_ptr<CrossSection> XS,
                       double T_IC,double f_sys = 0.):
      MELS(MELS),EA(EA),XS(XS),T_IC(T_IC),f_sys(f_sys)
      {
          Omega = 2.0*M_PI;
          flux_units = pow(units.GeV,-1.)*pow(units.cm,-2.)*pow(units.sec,-1.);
      }

      /// Basic Function used to calculate event spectation
      nusquids::marray<double,2>  DoIt(nusquids::marray<double,2> zenith_energy_bins,std::shared_ptr<Flux> FLUX,Sampling mode,
                                       double N0 = double(1.0), double dgamma = double(0.0)) {
          nusquids::marray<double,2> ResultTable{zenith_energy_bins.extent(0)-1,zenith_energy_bins.extent(1)-1};

          for(unsigned int ci = 0; ci < zenith_energy_bins.extent(0)-1; ci++){
              for(unsigned int ei = 0; ei < zenith_energy_bins.extent(1)-1; ei++){
                  double costh_edge = cth_range[ci];
                  double Er_edge = exp(loge_range[ei]);
                  double event_num,event_num_error;

                  #ifdef _NAIVE_INTEGRATOR_
                  // aqui viene la integracion =MAGIA=

                  const gsl_rng_type * T;
                  gsl_rng * r;
                  T = gsl_rng_default;
                  r = gsl_rng_alloc (T);

                  event_num = 0.0;
                  for(unsigned int ir = 0; ir < _NUM_RANDOM_INT_; ir++){

                      double Er = Er_edge + DeltaE[ei]*gsl_rng_uniform(r);
                      double Enu = Er_edge + 3.0 * DeltaE[ei] * (gsl_rng_uniform (r) - 0.5);
                      double costh = costh_edge + DeltaCosTh[ci]*gsl_rng_uniform (r);

                      event_num += Aeff(costh,Enu,0)*ANS->Evaluate(costh,Enu,1,0)*Gaussian(Er,Enu);
                      event_num += Aeff(costh,Enu,1)*ANS->Evaluate(costh,Enu,1,1)*Gaussian(Er,Enu);

                  }

                  gsl_rng_free(r);

                  // normalizamos
                  double DeltaEr = DeltaE[ei];
                  double DeltaEE = DeltaE[ei];
                  event_num = event_num*T_IC*Omega*DeltaEr*DeltaEE*DeltaCosTh[ci]*flux_units/_NUM_RANDOM_INT_;
                  #endif

                  #ifdef _VEGAS_INTEGRATOR_
                  //double xl[3] = { costh_edge, Er_edge/pc->GeV, 0.0};
                  //double xu[3] = { costh_edge + DeltaCosTh[ci], (Er_edge + DeltaE[ei])/pc->GeV, 1.0};
                  //double xl[3] = { costh_edge, Er_edge/pc->GeV, (Er_edge/pc->GeV)/10.0};
                  //double xu[3] = { costh_edge + DeltaCosTh[ci], (Er_edge + DeltaE[ei])/pc->GeV, ((Er_edge + DeltaE[ei])/pc->GeV)*10.0};
                  double xl[3] = { costh_edge, Er_edge/pc->GeV, ((Er_edge - DeltaE[ei])/pc->GeV)};
                  double xu[3] = { costh_edge + DeltaCosTh[ci], (Er_edge + DeltaE[ei])/pc->GeV, ((Er_edge + DeltaE[ei])/pc->GeV)};

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
                  ResultTable[ci][ei] = N0*event_num*pow(Er_edge/Emean_GeV,dgamma);
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
                    RealizationTable[i][j] = static_cast<double>(gsl_ran_poisson(r,ResultsTable[i][j]));
                  else if (mode == Gaussian){ // gaussian sampling
                    double sigma = sqrt(ResultsTable[i][j]) + f_sys*ResultsTable[i][j];
                    RealizationTable[i][j] = ResultsTable[i][j] + static_cast<double>(gsl_ran_gaussian(r,sigma));
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

    protected:
      double Kernel(double costh, double enu,
                    double muon_initial_energy, double muon_final_energy,
                    double muon_length){
        return T_IC*(*FLUX)(costh,enu)*(*XS)(enu,muon_initial_energy)*NT*(*MELS)(muon_initial_energy,muon_final_energy,muon_length)*(*EA)(muon_final_energy);
      }

      double Kernel(double *x){
          double costh = x[0];
          double Er = x[1]*pc->GeV;
          double Enu = x[2]*pc->GeV;

          if (Enu > exp(loge_range[e_size-1]) or Enu < exp(loge_range[0]))
              return 0.0;
          else	{
              return SQR(units.GeV)*(Aeff(costh,Enu,0)*ANS->Evaluate(costh,Enu,1,0)*Gaussian(Er,Enu)+
              Aeff(costh,Enu,1)*ANS->Evaluate(costh,Enu,1,1)*Gaussian(Er,Enu))*T_IC*Omega*flux_units;
          }
      }
};

double KernelHelper(double *x, size_t dim, void *params){
    IceCubeEstimator* ICE = (IceCubeEstimator*) params;
    return ICE->Kernel(x);
};

} // close namespace

#endif
