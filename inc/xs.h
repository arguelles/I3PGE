#ifndef _xs_
#define _xs_

#include <cmath>
#include <memory>
#include <photospline/splinetable.h>
#include <photospline/bspline.h>

namespace I3PGE {

struct CrossSection {
  using result_type = double;
  public:
    virtual result_type DoubleDifferentialCrossSection(int, double,double,double) const = 0;
    template<typename Event>
    double operator()(const Event& e) const {
      return DoubleDifferentialCrossSection((int)e.primaryType,e.injectedEnergy,e.intX,e.intY);
    };
};

struct CrossSectionFromSpline: public CrossSection {
  private:
    // photospline objects
    std::shared_ptr<splinetable> numu_dsdxdy;
    std::shared_ptr<splinetable> numubar_dsdxdy;
    double DoubleDifferentialCrossSection(int neutype, double nuEnergy,double x, double y) const;
  public:
    CrossSectionFromSpline(std::string splinepath, std::string model_name = "");
};

} // close I3PGE namespace

#endif
