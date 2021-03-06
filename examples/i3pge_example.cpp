#include <algorithm>
#include <deque>
#include <iterator>
#include <set>
#include <string>
#include <iostream>

#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/FCNBase.h>
#include "autodiff.h"
//#include "lbfgsb.h"

#include <boost/lexical_cast.hpp>

#define NO_STD_OUT

#define _LOADPREVIOUSCALCULATION_
#define _LOADPREVIOUSCALCULATION_STD_

#define _SAVE_SIM_
#define _SAVE_SIM_STD_

//#define _WRITE_EFA_

#define _SAVE_DATA_PATH_ "/Users/carguelles/Dropbox/NuestrosProyectos/CJAtmSterile/theory_data/"
#define _SIM_DATA_PATH_ "/Users/carguelles/Dropbox/NuestrosProyectos/CJAtmSterile/sim_data/"

// my libraries
#include "flux.h"
#include "xs.h"
#include "effective_area.h"

#include "i3pge.h"
#include "verosimilitud.h"

int main(int argc, char* argv[]){

    // get mixing angles and square mass difference
    double dm41sq,th24;

    if(argc != 3){
        printf("ERROR:USAGE: IC.exe dmsq41 [eVˆ2] th24 [rad] \n");
        exit(0);
    } else {
        dm41sq = atof(argv[1]);
        th24 = atof(argv[2]);
    }

    //dm41sq = 1.0;
    //th24 = 0.1;
#ifndef NO_STD_OUT
    std::cout << "Creando Caja y configunrandola" << std::endl;
#endif
    // setting up oscillation parameters
    PhysConst* pc = new PhysConst();
    pc->numneu = 4;
    // setup 3-neu-std angles and mass diff.
    pc->th12 = 0.5825;
    pc->th13 = 0.1591;
    pc->th23 = 0.7051;
    pc->delta1 = 0.0;
    pc->dm21sq = 0.758e-4;
    pc->dm31sq = 0.235e-2;
    // setup sterile angle and mass diff.
    pc->th24 = th24;
    pc->dm41sq = dm41sq;
    pc->Refresh();
    // setup oscillation
    SUBox* caja = new SUBox;
    caja->param = pc;
    // setup interactions
    caja->oscillation = true;
    caja->attenuation = true;
    caja->nc_inter = true;
    caja->cc_inter = true;
    caja->tau_regeneration = true;
    caja->neutype = 2;
    caja->numneu = 4;
    caja->body = new EarthAtm;

    vector<double> costh_array = linspace(-1.0,0.2,29);
    vector<double> enu_array_GeV = logspace(1.0e2,1.0e6,149);
    vector<double> enu_array_eV = logspace(1.0e2*pc->GeV,1.0e6*pc->GeV,149);


    //for(int ie = 0; ie < costh_array.size(); ie++)
    //    std::cout << costh_array[ie] << " ";
    //for(int ie = 0; ie < enu_array_GeV.size(); ie++)
    //    std::cout << enu_array_GeV[ie] << " ";
    //exit(0);

    caja->esize = (const int) enu_array_eV.size();
    copy(enu_array_eV.begin(),enu_array_eV.end(),caja->E_range);
#ifndef NO_STD_OUT
    std::cout << "Creando flujo atmosferico" << std::endl;
#endif
    // cargar objetos de new nu flux
    boost::shared_ptr<NewNuFlux::FluxFunction> flux_pion = NewNuFlux::makeFlux("honda2006");
    boost::dynamic_pointer_cast<NewNuFlux::PionKaonAdjustable>(flux_pion)->setRelativeKaonContribution(0);
    boost::dynamic_pointer_cast<NewNuFlux::KneeReweightable>(flux_pion)->setKneeReweightingModel("gaisserH3a_elbert");

    boost::shared_ptr<NewNuFlux::FluxFunction> flux_kaon = NewNuFlux::makeFlux("honda2006");
    boost::dynamic_pointer_cast<NewNuFlux::PionKaonAdjustable>(flux_kaon)->setRelativePionContribution(0);
    boost::dynamic_pointer_cast<NewNuFlux::KneeReweightable>(flux_kaon)->setKneeReweightingModel("gaisserH3a_elbert");

    SUVectorTable INeuStateTable;
    SUVectorTable IANeuStateTable;
    vector<EarthAtm::Track*> Tracks;

    for(unsigned int ci = 0; ci < costh_array.size(); ci++){
        SUVector INeuStateVector;
        SUVector IANeuStateVector;
        for(unsigned int ei = 0; ei < enu_array_GeV.size();ei++){
            double cth = costh_array[ci];
            double enu = enu_array_GeV[ei];
            INeuStateVector.push_back(SU_vector(Projector(pc->electron,pc->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuE,enu,cth) +
                                      SU_vector(Projector(pc->muon,pc->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuMu,enu,cth) +
                                      SU_vector(Projector(pc->tau,pc->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuTau,enu,cth) +
                                      SU_vector(Projector(pc->electron,pc->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuE,enu,cth) +
                                      SU_vector(Projector(pc->muon,pc->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuMu,enu,cth) +
                                      SU_vector(Projector(pc->tau,pc->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuTau,enu,cth));

            IANeuStateVector.push_back(SU_vector(Projector(pc->electron,pc->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuEBar,enu,cth) +
                                       SU_vector(Projector(pc->muon,pc->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuMuBar,enu,cth) +
                                       SU_vector(Projector(pc->tau,pc->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuTauBar,enu,cth) +
                                       SU_vector(Projector(pc->electron,pc->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuEBar,enu,cth) +
                                       SU_vector(Projector(pc->muon,pc->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuMuBar,enu,cth) +
                                       SU_vector(Projector(pc->tau,pc->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuTauBar,enu,cth));
        }
        INeuStateTable.push_back(INeuStateVector);
        IANeuStateTable.push_back(IANeuStateVector);

        Tracks.push_back(new EarthAtm::Track(acos(costh_array[ci])));
    }

    // el flujo no estandar
    AtmNeutrinoState ANF(INeuStateTable,IANeuStateTable,Tracks,caja);

    #ifndef _LOADPREVIOUSCALCULATION_
    // evolucionamos el flujo para cada angulo de zenith
    ANF.PropagateState();
    // guardamos los valores
    string save_path = _SAVE_DATA_PATH_;
    ANF.SaveState(save_path + "diff_"+toString(dm41sq)+"_"+toString(th24)+".dat");
    #else
    // cargamos el calcula anterior
    #ifndef NO_STD_OUT
        std::cout << "Cargando Flujo Precalculado" << std::endl;
    #endif
    string save_path = _SAVE_DATA_PATH_;
    ANF.LoadState(save_path + "diff_"+toString(dm41sq)+"_"+toString(th24)+".dat");
    #endif

    // calculamos el caso STD
    SUBox* caja_std = new SUBox(*caja);

    caja_std->param->th24 = 0.0;
    caja_std->param->dm41sq = 0.0;
    caja_std->param->numneu = 3;
    caja_std->param->Refresh();
    caja_std->numneu = 3;


    SUVectorTable INeuStateTable_STD;
    SUVectorTable IANeuStateTable_STD;
    vector<EarthAtm::Track*> Tracks_STD;

    for(unsigned int ci = 0; ci < costh_array.size(); ci++){
        SUVector INeuStateVector;
        SUVector IANeuStateVector;
        for(unsigned int ei = 0; ei < enu_array_GeV.size();ei++){
            double cth = costh_array[ci];
            double enu = enu_array_GeV[ei];
            INeuStateVector.push_back(SU_vector(Projector(pc->electron,caja_std->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuE,enu,cth) +
                                      SU_vector(Projector(pc->muon,caja_std->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuMu,enu,cth) +
                                      SU_vector(Projector(pc->tau,caja_std->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuTau,enu,cth) +
                                      SU_vector(Projector(pc->electron,caja_std->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuE,enu,cth) +
                                      SU_vector(Projector(pc->muon,caja_std->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuMu,enu,cth) +
                                      SU_vector(Projector(pc->tau,caja_std->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuTau,enu,cth));

            IANeuStateVector.push_back(SU_vector(Projector(pc->electron,caja_std->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuEBar,enu,cth) +
                                       SU_vector(Projector(pc->muon,caja_std->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuMuBar,enu,cth) +
                                       SU_vector(Projector(pc->tau,caja_std->numneu))*flux_pion->getFlux(NewNuFlux::I3Particle::NuTauBar,enu,cth) +
                                       SU_vector(Projector(pc->electron,caja_std->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuEBar,enu,cth) +
                                       SU_vector(Projector(pc->muon,caja_std->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuMuBar,enu,cth) +
                                       SU_vector(Projector(pc->tau,caja_std->numneu))*flux_kaon->getFlux(NewNuFlux::I3Particle::NuTauBar,enu,cth));
        }
        INeuStateTable_STD.push_back(INeuStateVector);
        IANeuStateTable_STD.push_back(IANeuStateVector);

        Tracks_STD.push_back(new EarthAtm::Track(acos(costh_array[ci])));
    }

    // el flujo estandar
    AtmNeutrinoState ANF_STD(INeuStateTable_STD,IANeuStateTable_STD,Tracks_STD,caja_std);

    #ifndef _LOADPREVIOUSCALCULATION_STD_
        ANF_STD.PropagateState();
        string save_path_std = _SAVE_DATA_PATH_;
        ANF_STD.SaveState(save_path_std + "diff_0_0.dat");
    #else
        #ifndef NO_STD_OUT
            std::cout << "Cargando Flujo STD Precalculado" << std::endl;
        #endif
        string save_path_std = _SAVE_DATA_PATH_;
        ANF_STD.LoadState(save_path_std + "diff_0_0.dat");
    #endif

    // creamos un icecube coso
    //IceCubeEstimator ICE(1.0*pc->year);
    IceCubeEstimator* ICE = new IceCubeEstimator(1.0*pc->year,0.15);

	// check effective area

    #ifdef _WRITE_EFA_
	Table effective_area;
	for(std::vector<double>::iterator costh = costh_array.begin(); costh != costh_array.end(); ++costh){
			for(std::vector<double>::iterator enu = enu_array_eV.begin(); enu != enu_array_eV.end(); ++enu){
					Row efa;

					efa.push_back(*enu);
					efa.push_back(*costh);
					efa.push_back(ICE->GetEffectiveArea(*costh,*enu,0));

					effective_area.push_back(efa);
			}
	}
	string efa_path = "./effective_area.dat";
	quickwrite(efa_path,effective_area);
    #endif

    Table events_std = ICE->DoIt(&ANF_STD,2);

    #ifdef _SAVE_SIM_STD_
        string sim_std_path = _SIM_DATA_PATH_;
        quickwrite(sim_std_path + "simdata_0_0.dat",events_std);
    #endif

    Table events_t = ICE->DoIt(&ANF,0);
/*
    string sim_path_t = _SIM_DATA_PATH_;
    quickwrite(sim_path_t + "sim_"+toString(dm41sq)+"_"+toString(th24)+".dat",events_t);
    exit(0);
*/
//    PrintTable(events_std);

    GaussianPrior* SlopePrior = new GaussianPrior(0.0,0.05);
    UniformPrior* NormalizationPrior = new UniformPrior(0.0,100.0);

    std::vector<Prior*> priors;// = {NormalizationPrior,SlopePrior};
    priors.push_back(NormalizationPrior);
    priors.push_back(SlopePrior);

    std::vector<double> params_seed;// = {1.0,0.0};
    params_seed.push_back(1.0);
    params_seed.push_back(0.0);

  // otro coso

    LBFGSB_Driver minimizer;
    minimizer.addParameter(params_seed[0],.1,0.0);
    minimizer.addParameter(params_seed[1],.05);

    minimizer.setChangeTolerance(5e-8);
    minimizer.setHistorySize(20);

    VerosimilitudBFGSAdapter VMA(events_std, ANF, priors, ICE);
    bool succeeded=minimizer.minimize(VerosimilitudBFGSAdapter(VMA));
    double likelihood=minimizer.minimumValue();
    vector<double> finalparams = minimizer.minimumPosition();

    std::cout << likelihood << std::endl;

    /*
    // verosimilitud coso
    VerosimilitudMinuitAdapter VMA(events_std, ANF, priors, ICE);
    ROOT::Minuit2::MnUserParameters minuit_params;
	//copy seed values

    for(unsigned int i = 0; i < params_seed.size(); i++){
        minuit_params.Add(priors[i]->name,params_seed[i],0.1);
    }
    minuit_params.Add("p1",params_seed[0],0.1);
    minuit_params.Add("p2",params_seed[1],0.1);
    // normalization cannot be negative
    minuit_params.SetLowerLimit(0,0.0);
    minuit_params.SetLowerLimit(1,-1.0);
    minuit_params.SetUpperLimit(1,1.0);

	ROOT::Minuit2::VariableMetricMinimizer minimizer;
	ROOT::Minuit2::FunctionMinimum fmin=minimizer.Minimize(VMA, minuit_params, ROOT::Minuit2::MnStrategy(2), 0, .01);

	const ROOT::Minuit2::MnUserParameters& finalParams=fmin.UserParameters();
	for(unsigned int i=0; i<params_seed.size(); i++)
		params_seed[i] = finalParams.Value(i);
	double likelihood=fmin.Fval();
  */

    //std::cout << fmin.NFcn() << " " <<  likelihood << std::endl;
	for(unsigned int i=0; i<finalparams.size(); i++)
        std::cout << finalparams[i] << " ";
    std::cout << likelihood << " ";
    std::cout << dm41sq << " " << th24 << std::endl;

#ifdef _SAVE_SIM_
    ChangeParameters(&ANF,finalparams);
    Table events = ICE->DoIt(&ANF,0);

    string sim_path = _SIM_DATA_PATH_;
    quickwrite(sim_path + "sim_"+toString(dm41sq)+"_"+toString(th24)+".dat",events);
#endif

    return(0);
}
