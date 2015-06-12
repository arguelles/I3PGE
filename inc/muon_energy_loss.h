#ifndef _muon_energy_loss
#define _muon_energy_loss

#include <iostream>
#include <memory>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include <photospline/splinetable.h>
#include <photospline/bspline.h>

namespace I3PGE {

struct MuonEnergyLossSpline {
  private:
    // photospline objects
    std::shared_ptr<splinetable> muon_energy_loss_spline;
    double Evaluate(double costh,double distance,double initial_energy,double final_energy) const;
  public:
    MuonEnergyLossSpline(std::string splinepath);
    template<typename... ArgTypes>
    double operator()(ArgTypes&&... args){ return Evaluate(std::forward<ArgTypes>(args)...);};
};

} // close namespace

#endif
