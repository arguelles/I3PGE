#ifndef _effective_area_
#define _effective_area_

#include <iostream>
#include <memory>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include <photospline/splinetable.h>
#include <photospline/bspline.h>
#include <SQuIDS/Const.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/tools.h>

namespace I3PGE {

struct MuonEffectiveAreaSpline {
  private:
    // photospline objects
    std::shared_ptr<splinetable> effective_area_spline;
    double Evaluate(double costh,double muon_energy) const;
  public:
    template<typename... ArgTypes>
    double operator()(ArgTypes&&... args) const { return Evaluate(std::forward<ArgTypes>(args)...);};
    MuonEffectiveAreaSpline(std::string splinepath);
};

struct NeutrinoEffectiveAreaTable {
  private:
    squids::Const units;
    nusquids::marray<double,2> AeffTable;
    std::vector<double> cth_range;
    std::vector<double> loge_range;
    double Evaluate(double costh, double energy,unsigned int neutype) const;
  protected:
    double LinInter(double x,double xM, double xP, double yM, double yP) const {
      return yM + (yP-yM)*(x-xM)/(xP-xM);
    }
  public:
    template<typename... ArgTypes>
    double operator()(ArgTypes&&... args){ return Evaluate(std::forward<ArgTypes>(args)...);};
    NeutrinoEffectiveAreaTable(std::string tablepath,unsigned int, unsigned int);
};

} // close namespace

#endif
