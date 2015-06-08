#ifndef _ICECUBE_H_
#define _ICECUBE_H_

#include "gsl/gsl_math.h"
#include <gsl/gsl_rng.h>
#include "gsl/gsl_randist.h"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>

#include <SQuIDS/const.h>

#include "flux.h"
#include "xs.h"

#define _IC_AEFF_PATH_ "/Users/carguelles/Dropbox/NuestrosProyectos/CJDiffuseSterile/c++/data/"
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

enum Sampling { Asimov, Poisson, Gaussian };

class IceCubeEstimator {
    private:
      /// SQuIDS units module
      squids::Const units;
      /// Cross section calculator
      std::shared_ptr<CrossSection> XS;
      /// Flux calculator
      std::shared_ptr<Flux> FLUX;
      /// Effective area calculator
      std::shared_ptr<MuonEffectiveAreaSpline> EA;
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
      double Gaussian(double Er, double Et){
          sigmaE = Et;

          double norm = 1.0/(sqrt(2.0*M_PI)*sigmaE);
          double z = (Er-Et)/sigmaE;
          return norm*exp(-z*z/2.0);
      };
      /// energy - zenith support array
      int cth_size;
      int e_size;
      vector<double> cth_range;
      vector<double> loge_range;
      vector<double> DeltaCosTh;
      vector<double> DeltaE;
    public:
      /// Basic constructor
      IceCubeEstimator(double T_IC,double f_sys = 0.):
      T_IC(T_IC),f_sys(f_sys)
      {
          Omega = 2.0*M_PI;

          string path = _IC_AEFF_PATH_;
          Table AeffTableRaw = quickread(path+"icecube_cw_effarea.dat");
          for(int ci = 0; ci < _IC_AEFF_COSTHBINS_ ; ci++){
              cth_range.push_back(AeffTableRaw[ci*_IC_AEFF_EBINS_][1]);
              Row AeffRow;
              for(int ei = 0; ei < _IC_AEFF_EBINS_ ; ei++) {
                  if(ci == 0)
                      loge_range.push_back(log(AeffTableRaw[ei][0]*pc->GeV));// convert to log(E [eV])
                  AeffRow.push_back(AeffTableRaw[ci*_IC_AEFF_EBINS_+ei][2]*SQR(pc->meter)); // convert to eV^-2
              }
              AeffTable.push_back(AeffRow);
          }

          cth_size = cth_range.size();
          e_size = loge_range.size();

          for(int ci = 0; ci < cth_size-1; ci++)
              DeltaCosTh.push_back(cth_range[ci+1]-cth_range[ci]);
          DeltaCosTh.push_back(DeltaCosTh.back());
          for(int ei = 0; ei < e_size-1; ei++)
              DeltaE.push_back(exp(loge_range[ei+1])-exp(loge_range[ei]));
          DeltaE.push_back(DeltaE.back());

          flux_units = pow(pc->GeV,-1)*pow(pc->cm,-2)*pow(pc->sec,-1);
      }

      /// Basic Function used to calculate event spectation
      nusquids::marray<double,2>  DoIt(AtmNeutrinoState* ANS_in,int mode,double N0 = double(1.0), double dgamma = double(0.0)) {
          ANS = ANS_in;
          Table ResultsTable;
          for(int ci = 0; ci < cth_size; ci++){
              Row ResultsRow;
              for(int ei = 0; ei < e_size; ei++){
                  double costh_edge = cth_range[ci];
                  double Er_edge = exp(loge_range[ei]);
                  double event_num,event_num_error;

                  if( Er_edge < exp(ANS->loge_range[0]) or Er_edge > exp(ANS->loge_range.back())){
                      ResultsRow.push_back(0.0);
                      continue;
                  }

                  #ifdef _NAIVE_INTEGRATOR_
                  // aqui viene la integracion =MAGIA=

                  const gsl_rng_type * T;
                  gsl_rng * r;
                  T = gsl_rng_default;
                  r = gsl_rng_alloc (T);

                  event_num = 0.0;
                  for(int ir = 0; ir < _NUM_RANDOM_INT_; ir++){

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
                    //std::cout << gsl_monte_vegas_chisq(s) << std::endl;
                                if (gsl_monte_vegas_chisq(s) == 0)
                                    break;
                    }
                  while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.15);

                  gsl_monte_vegas_free(s);
                  gsl_rng_free(r);
                  #endif

                  double Emean_GeV = 1827.03*pc->GeV;
                  ResultsRow.push_back(N0*event_num*pow(Er_edge/Emean_GeV,dgamma));
              }
              ResultsTable.push_back(ResultsRow);
          }
          if (mode == 0){
              // returns expected event number per bin
              return ResultsTable;
          } else if (mode == 1 or mode == 2){
              // returns poisson realizations of the event number per bin
              Table RealizationTable;

              const gsl_rng_type * T;
              gsl_rng * r;
              T = gsl_rng_default;
              r = gsl_rng_alloc (T);

              for(int i = 0; i< cth_size; i++){
                  Row RealizationRow;
                  for(int j = 0; j < e_size; j++){
                      if( ResultsTable[i][j] == 0.0 )
                          RealizationRow.push_back(0.0);
                      else {
                          if (mode == 1) // poisson sampling
                            RealizationRow.push_back((double) gsl_ran_poisson(r,ResultsTable[i][j]));
                          if (mode == 2){ // gaussian sampling
                            double sigma = sqrt(ResultsTable[i][j]) + f_sys*ResultsTable[i][j];
                            RealizationRow.push_back(ResultsTable[i][j] + (double) gsl_ran_gaussian(r,sigma));
                          }
                      }
                  }
                  RealizationTable.push_back(RealizationRow);
              }

              gsl_rng_free(r);

              return RealizationTable;
          } else {
              std::cout << "Error : wrong IceCube estimator mode" << std::endl;
              exit(0);
          }
      }
    protected:
      double Kernel(double *x){
          double costh = x[0];
          double Er = x[1]*pc->GeV;
          double Enu = x[2]*pc->GeV;

          if (Enu > exp(loge_range[e_size-1]) or Enu < exp(loge_range[0]))
              return 0.0;
          else	{
              return SQR(pc->GeV)*(Aeff(costh,Enu,0)*ANS->Evaluate(costh,Enu,1,0)*Gaussian(Er,Enu)+
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
