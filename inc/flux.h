#ifndef _flux_
#define _flux_

#include <cmath>
#include <iostream>
#include <memory>
#include <nuSQuIDS/nuSQUIDS.h>

/*
 * Requerimientos:
 * - Informacion de secciones de choque de neutrinos
 * - Funcion estimadora de flujos
*/

namespace I3PGE {

struct Flux {
  using result_type = double;
  // Contenedor de flujo. En principio no hace nada
  // hasta ser especializado a la nuSQuIDS o
  // siguiendo la prescripcion de NeutrinoFlux
  public:
    virtual result_type EvaluateFlux(unsigned int neutype, double costh, double enu) const = 0;
    result_type operator()(unsigned int neutype,double costh, double enu) const { return EvaluateFlux(neutype,costh,enu);};
};

struct SQUIDSFlux: public Flux {
  private:
      const squids::Const units;
      unsigned int muon = 1;
      unsigned int neutrino = 0; unsigned int antineutrino = 1;
  protected:
      nusquids::nuSQUIDSAtm<> nsqa;
  public:
    using result_type = double;
    result_type EvaluateFlux(unsigned int neutype, double costh, double enu) const override {
      // asssume its always nu_mu or nu_mu bar
      if ( neutype == neutrino )
        return nsqa.EvalFlavor(muon,costh,enu,neutrino);
      else if ( neutype == antineutrino)
        return nsqa.EvalFlavor(muon,costh,enu,antineutrino);
      else {
        throw std::runtime_error("SQUIDSFlux::Error::Unkwown particle type : " + std::to_string(neutype) + ".");
      }
    };
    SQUIDSFlux(std::string nusquids_data_file_path): nsqa(nusquids::nuSQUIDSAtm<>(nusquids_data_file_path)) {};
    SQUIDSFlux(nusquids::nuSQUIDSAtm<>&& nsqa): nsqa(std::move(nsqa)) {};
};

} // end namespace

#endif
