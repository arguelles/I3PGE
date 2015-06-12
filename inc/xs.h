#ifndef _xs_
#define _xs_

#include <cmath>
#include <memory>
#include <photospline/splinetable.h>
#include <photospline/bspline.h>

namespace I3PGE {

struct I3PGECrossSection: public nusquids::NeutrinoDISCrossSectionsFromTables {
  private:
    // photospline objects
    std::shared_ptr<splinetable> numu_dsdxdy;
    std::shared_ptr<splinetable> numubar_dsdxdy;
  public:
    I3PGECrossSection(std::string splinepath, std::string model_name = "");
    double DoubleDifferentialCrossSection(double E, double x, double y, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override;
    template<typename... ArgTypes>
    double operator()(ArgTypes&&... args){ return SingleDifferentialCrossSection(std::forward<ArgTypes>(args)...);};
};

} // close I3PGE namespace

#endif
