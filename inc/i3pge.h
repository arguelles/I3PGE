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
#include "verosimilitud.h"

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

/*
template<typename DataType>
void ChangeParameters(std::vector<Flux> FLUX,std::vector<DataType> params){
    double N0 = params[0];
    double dgamma = params[1];

    std::vector<NeutrinoState> NST = ANS->GetState();

    // mean energy for atmospheric in 1e2 to 1e6 range
    double Emean_GeV = 1827.03;

    for(int it = 0; it < ANS->cth_size; it++){
        for(int ie = 0 ; ie < ANS->e_size; ie++){
            double E_GeV = exp(ANS->loge_range[ie])/1.0e9;
            NST[it].InitialNeutrinoState[ie] = NST[it].InitialNeutrinoState[ie]*(pow(E_GeV/Emean_GeV,dgamma)*N0);
            NST[it].InitialANeutrinoState[ie] = NST[it].InitialANeutrinoState[ie]*(pow(E_GeV/Emean_GeV,dgamma)*N0);
        }
    }

    ANS->SetState(NST);
}
*/

template<typename DataType>
DataType CalculateVerosimilitud(nusquids::marray<double,2> obs_events,
    nusquids::marray<DataType,2> sim_events,
    std::vector<DataType> params, std::vector<Prior*> priors){
    if (obs_events.size() != sim_events.size()){
        std::cout << "Error : Different obs/sim number" << std::endl;
        exit(0);
    }

    std::vector<double> obs;
    std::vector<DataType> sim;

    for(int it = 0; it < (int) obs_events.size(); it++){
        for(int ei = 0; ei < (int) obs_events[0].size(); ei++){
            if ( sim_events[it][ei] != DataType(0) ){
                obs.push_back(obs_events[it][ei]);
                sim.push_back(sim_events[it][ei]);
            }
        }
    }

    double f = 0.15;
    Verosimilitud<DataType> Verdad(obs,sim,priors,2,f);
    return Verdad.EvaluateVerosimilitud(params);
}

class VerosimilitudMinuitAdapter : public ROOT::Minuit2::FCNBase {
private:
    nusquids::marray<double,2> obs;
    //nusquids::marray<double,2> sim;
    std::shared_ptr<Flux> FLUX;
    std::vector<std::shared_ptr<Prior>> priors;
    std::vector<double> params;
    std::shared_ptr<IceCubeEstimator> ICE;

public:
    VerosimilitudMinuitAdapter(nusquids::marray<double,2> obs,
                               AtmNeutrinoState ANS,
                               std::vector<std::shared_ptr<Prior>> priors,
                               std::shared_ptr<IceCubeEstimator> ICE): obs(obs),ANS(ANS),priors(priors),ICE(ICE){};

    double Evaluate(const std::vector<double>& params) const {
        //cout << "Evaluating coso" << endl;
        AtmNeutrinoState* ANF_TMP = new AtmNeutrinoState(ANS);
        ChangeParameters(ANF_TMP,params);
        //cout << "runeval2" << endl;
        nusquids::marray<double,2> sim = ICE->DoIt(ANF_TMP,0);

        //PrintTable(sim);

        delete ANF_TMP;
        ANF_TMP = NULL;
        //cout << "runeval3" << endl;
        double ver = CalculateVerosimilitud(obs,sim,params,priors);
        //cout << params[0] << " " << params[1] << " " << ver << endl;
        return ver;
        //return CalculateVerosimilitud(obs,sim,params,priors);
    }

    double operator()(const std::vector<double>& params) const{
        return (Evaluate(params));
    }

    double Up() const{
        return 0.5;
    }

};

class VerosimilitudBFGSAdapter : public BFGS_FunctionBase {
private:
    nusquids::marray<double,2> obs;
    AtmNeutrinoState ANS;
    std::vector<std::shared_ptr<Prior>> priors;
    std::vector<double> params;
    std::shared_ptr<IceCubeEstimator> ICE;
    static constexpr int DerivativeDimension = Dynamic;
    template<typename DataType>
    nusquids::marray<DataType,2> ChangeParameters(nusquids::marray<DataType,2> sim,std::vector<DataType> params){
      DataType N0 = params[0];
      DataType dgamma = params[1];

      std::cout << N0 << " " << dgamma << std::endl;

      // mean energy for atmospheric in 1e2 to 1e6 range
      double Emean_GeV = 1827.03;
      std::vector<std::vector<DataType>> ch_sim;
      for(int it = 0; it < ICE->cth_size; it++){
          std::vector<DataType> row_sim;
          for(int ie = 0 ; ie < ICE->e_size; ie++){
              double E_GeV = exp(ICE->loge_range[ie])/1.0e9;
              row_sim.push_back(N0*pow(E_GeV/Emean_GeV,dgamma)*sim[it][ie]);
              //std::cout << N0*pow(E_GeV/Emean_GeV,dgamma)*sim[it][ie] << std::endl;
          }
          ch_sim.push_back(row_sim);
      }
      return ch_sim;
    }
public:
    VerosimilitudBFGSAdapter(nusquids::marray<double,2> obs,
                             std::shared_ptr<Flux> FLUX,
                             std::vector<std::shared_ptr<Prior>> priors,
                             std::shared_ptr<IceCubeEstimator> ICE): obs(obs),ANS(ANS),priors(priors),ICE(ICE){};

    template<typename DataType>
    DataType evalF(std::vector<DataType> params){
        nusquids::marray<DataType,2> sim = ICE->DoIt(&ANS,0);
        //nusquids::marray<DataType,2> ch_sim = ChangeParameters(sim,params);
        DataType ver = CalculateVerosimilitud(obs,sim,params,priors);
        return ver;
    }

    double evalF(std::vector<double> params){
        nusquids::marray<double,2> sim = ICE->DoIt(&ANS,0);
        //vector<vector<double>> ch_sim = ChangeParameters(sim,params);
        double ver = CalculateVerosimilitud(obs,sim,params,priors);
        return ver;
    }

    template<typename DataType>
    DataType operator()(const std::vector<DataType> params){
        return (evalF(params));
    }

    std::pair<double,std::vector<double>> evalFG(std::vector<double> x){
        const size_t size=x.size();
        std::vector<FD<DerivativeDimension>> params(size);
        for(size_t i=0; i<size; i++)
            params[i]=FD<DerivativeDimension>(x[i],i);

        FD<DerivativeDimension> result=evalF<FD<DerivativeDimension>>(params);

        std::vector<double> grad(size);
        for(unsigned int i=0; i<size; i++)
          grad[i]=result.derivative(i);

        return(std::make_pair(result.value(),grad));
    }

};

} // close namespace

#endif
