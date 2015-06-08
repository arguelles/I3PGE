#include "effective_area.h"

namespace I3PGE {

double MuonEffectiveAreaSpline::Evaluate(double costh, double muon_energy) const {
  int centerbuffer[2];
  double xx[2];

  xx[0] = costh;
  xx[1] = log10(muon_energy);

  if(tablesearchcenters(effective_area_spline.get(),xx,centerbuffer) == 0)
    return pow(10.0,ndsplineeval(effective_area_spline.get(),xx,centerbuffer,0));
  else
    return 0.0;
}

MuonEffectiveAreaSpline::MuonEffectiveAreaSpline(std::string splinepath):
  effective_area_spline(new splinetable,[](splinetable* t){ splinetable_free(t); delete t;})
{
  // lets photospline :)
  readsplinefitstable(splinepath.c_str(),effective_area_spline.get());
  assert(effective_area_spline->ndim>0);
}

double NeutrinoEffectiveAreaTable::Evaluate(double costh, double energy,unsigned int neutype) const {
  double logE = log(energy);

  if (logE > loge_range.back() or logE < loge_range.front() or
      costh > cth_range.back() or costh < cth_range.front()){
      return 0.0;
  } else {
      int cth_M = -1;
      int cth_size = cth_range.size();
      for(int i = 0; i < cth_size; i++){
          //std::cout << costh << " " << cth_range[i] << std::endl;
          if ( costh >= cth_range[i] and costh < cth_range[i+1] ){
              cth_M = i;
              break;
          }
      }

      int loge_M = -1;
      int e_size = loge_range.size();
      for(int i = 0; i < e_size; i++){
          //std::cout << logE << " " << loge_range[i] << std::endl;
          if ( logE >= loge_range[i] and logE < loge_range[i+1] ){
              loge_M = i;
              break;
          }
      }

      double AMM,AMP,APM,APP;

      AMM = AeffTable[cth_M][loge_M];
      AMP = AeffTable[cth_M][loge_M+1];
      APM = AeffTable[cth_M+1][loge_M];
      APP = AeffTable[cth_M+1][loge_M+1];

      return 0.5*LinInter(costh,cth_range[cth_M],cth_range[cth_M+1],
                          LinInter(logE,loge_range[loge_M],loge_range[loge_M+1],AMM,AMP),
                          LinInter(logE,loge_range[loge_M],loge_range[loge_M+1],APM,APP));
  }
}

NeutrinoEffectiveAreaTable::NeutrinoEffectiveAreaTable(std::string tablepath,
                                                       unsigned int ic_costh_bin_num,
                                                       unsigned int ic_enu_bin_num){
  nusquids::marray<double,2> AeffTableRaw = nusquids::quickread(tablepath);
  if (AeffTableRaw.size() != ic_costh_bin_num*ic_enu_bin_num)
    throw std::runtime_error("Error at constructing NeutrinoEffectiveAreaTable table size and number of bins\
                              do not match.");

  AeffTable.resize(std::vector<size_t>{ic_costh_bin_num,ic_enu_bin_num});

  for(unsigned int ci = 0; ci < ic_costh_bin_num; ci++){
    cth_range.push_back(AeffTableRaw[ci*ic_enu_bin_num][1]);
    for(unsigned int ei = 0; ei < ic_enu_bin_num; ei++){
      if(ci == 0)
        loge_range.push_back(log(AeffTableRaw[ei][0]*units.GeV));// convert to log(E [eV])
      AeffTable[ci][ei] = AeffTableRaw[ci*ic_enu_bin_num+ei][2]*units.meter*units.meter; // convert to eV^-2
    }
  }
}

} // end namespace

