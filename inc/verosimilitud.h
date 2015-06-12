#ifndef __VEROSIMILITUD_H
#define __VEROSIMILITUD_H

#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <vector>
#include "autodiff.h"

class Prior {
  public:
    std::string name;
    bool fix;
    Prior(){};

    virtual double evaluate(double);
    virtual FD<Dynamic, double> evaluate(FD<Dynamic, double>);

    void FixParameter(void){
        fix = true;
    };

      //~Prior(){};
};

class GaussianPrior : public Prior {
  private:
    double mean;
    double stddev;
    double norm;
  public:
    GaussianPrior(double mean, double stddev):
    mean(mean),stddev(stddev),
    norm(1.0/(sqrt(2.0*M_PI)*stddev))
    {
        name = "GaussianPrior";
    };

    double evaluate(double x){
        double z = (x - mean)/stddev;
        return norm*exp(-z*z/2.0);
    };

    FD<Dynamic, double> evaluate(FD<Dynamic, double> x){
        FD<Dynamic, double> z = (x - mean)/stddev;
        return norm*exp(-z*z/2.0);
    };
};

class UniformPrior : public Prior {
  private:
    double min;
    double max;
  public:
    UniformPrior(double min, double max):
    min(min),max(max)
    {
        name = "UniformPrior";
    };

    double evaluate(double x){
        if (x < min || x > max)
            return double(0.0);
        return double(1.0);
    };

    FD<Dynamic, double> evaluate(FD<Dynamic, double> x){
        if (x < min || x > max)
            return FD<Dynamic, double>(0.0);
        return FD<Dynamic, double>(1.0);
    };
};

template<typename DataType>
class Verosimilitud{
    private:
        std::vector<Prior*> PriorVector;
        std::vector<double> ParameterVector;
        std::vector<double> ObservationVector;
        std::vector<DataType> SimulationVector;

        const int mode;
        int BinNumber;
        const double f;

        DataType Poisson(double obs,DataType sim){
            if (sim == 0)
                return DataType(0.0);
            return exp(-sim)*pow(sim,obs)/gsl_sf_gamma(obs+1.0);
        };
        DataType LogPoisson(double obs,DataType sim){
            //if (sim == 0)
            //    return 0.0;
            return obs*log(sim)-sim-gsl_sf_lngamma(obs+1.0);
        };
        DataType Chi2(double obs,DataType sim){
            return (obs-sim)*(obs-sim)/sim;
        };

        DataType ModChi2(double obs,DataType sim){
            //return (obs-sim)*(obs-sim)/(sim+SQR(f*sim));
            return (obs-sim)*(obs-sim)/(sim);
            //return (obs-sim)*(obs-sim)/(obs+SQR(f*obs));
        };

    public:
        /*
        Verosimilitud(std::vector<double> obs,
                      std::vector<DataType> sim,
                      std::vector<Prior*> priors,
                      int mode):
        ObservationVector(obs),
        SimulationVector(sim),
        PriorVector(priors),
        mode(mode),f(0.)
        {
            if (sim.size() == obs.size()){
                BinNumber = sim.size();
                //std::cout << BinNumber << std::endl;
            }
            else {
                std::cout << "VEROSILIMITUD VECTOR SIZE ERROR" << std::endl;
                exit(0);
            }
        };
        */

        Verosimilitud(std::vector<double> obs,
                      std::vector<DataType> sim,
                      std::vector<Prior*> priors,
                      int mode,double f = 0.0):
        ObservationVector(obs),
        SimulationVector(sim),
        PriorVector(priors),
        mode(mode), f(f)
        {
            if (sim.size() == obs.size()){
                BinNumber = sim.size();
                //std::cout << BinNumber << std::endl;
            }
            else {
                std::cout << "VEROSILIMITUD VECTOR SIZE ERROR" << std::endl;
                exit(0);
            }

          if (mode!=2) std::cout << "VEROSIMILITUD f value suplied but mode not 2." <<  std::endl;
        };

        DataType EvaluateLogLikelihood(std::vector<DataType> params){
            DataType llh = DataType(0.0);
            for(int i = 0; i < BinNumber; i++){
                llh += LogPoisson(ObservationVector[i],SimulationVector[i]);
                //std::cout << i << std::endl;
                //std::cout << ObservationVector[i] << " " << SimulationVector[i] << std::endl;
                //std::cout << LogPoisson(ObservationVector[i],SimulationVector[i]) << std::endl;
            }
            for(unsigned int i = 0; i < PriorVector.size(); i++){
                llh += log(PriorVector[i]->evaluate(params[i]));
            }

            return -llh;
        };

        DataType EvaluateChi2(std::vector<DataType> params){
            DataType llh = DataType(0.0);
            for(int i = 0; i < BinNumber; i++){
                //llh += Chi2(ObservationVector[i],SimulationVector[i]);
                llh += ModChi2(ObservationVector[i],SimulationVector[i]);
            }

            std::cout << "llh1 " << llh << std::endl;
            for(unsigned int i = 0; i < PriorVector.size(); i++){
                llh -= log(PriorVector[i]->evaluate(params[i]));
            }
            std::cout << "llh2 " << llh << std::endl;

            return llh;
        };

        DataType EvaluateVerosimilitud(std::vector<DataType> params){
            if (mode == 0)
                return EvaluateLogLikelihood(params);
            else if (mode == 1 or mode == 2)
                return EvaluateChi2(params);
            else
                return DataType(0.0);
                std::cout << "VEROSILIMITUD MODE ERROR" << std::endl;
                exit(0);
        };

        void UpdateSimulation(std::vector<DataType> sim){
            SimulationVector = sim;
        };

        void UpdateObservation(std::vector<double> obs){
            ObservationVector = obs;
        };

        ~Verosimilitud(){
            ParameterVector.clear();
            ObservationVector.clear();
            SimulationVector.clear();
            //for(int i = 0; i < (int) PriorVector.size(); i++){
            //    delete PriorVector[i];
            //}
            //PriorVector.clear();
        };
};

#endif
