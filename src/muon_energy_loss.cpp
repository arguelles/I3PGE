#include "muon_energy_loss.h"

namespace I3PGE {

double MuonEnergyLossSpline::Evaluate(double costh, double distance, double initial_energy, double final_energy) const {
  int centerbuffer[4];
  double xx[4];

  xx[0] = costh;
  xx[1] = log10(distance);
  xx[2] = log10(initial_energy);
  xx[3] = log10(final_energy);

  if(tablesearchcenters(muon_energy_loss_spline.get(),xx,centerbuffer) == 0)
    return ndsplineeval(muon_energy_loss_spline.get(),xx,centerbuffer,0);
  else
    return 0.0;
}

MuonEnergyLossSpline::MuonEnergyLossSpline(std::string splinepath):
  muon_energy_loss_spline(new splinetable,[](splinetable* t){ splinetable_free(t); delete t;})
{
  // lets photospline :)
  readsplinefitstable(splinepath.c_str(),muon_energy_loss_spline.get());
  assert(muon_energy_loss_spline->ndim>0);
}

} // close namespace
