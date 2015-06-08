#include "lepton_weighter.h"

namespace I3PGE {

double CrossSectionFromSpline::DoubleDifferentialCrossSection(int neutype, double nuEnergy,double x, double y) const {
  int centerbuffer[3];
  double xx[3];

  //int target = 2; // target = isospin
  //int cc = 1; // interaction = charge current

  xx[0] = log10(nuEnergy);
  xx[1] = log10(x);
  xx[2] = log10(y);

  if (neutype == (int) particleType::MuMinus or neutype == (int) particleType::NuMu){
    //int flavor = 3;// muon neutrino
    //return dsdxdy_(&nuEnergy,&x,&y,&flavor,&target,&cc);///sigma;
    if(tablesearchcenters(numu_dsdxdy.get(),xx,centerbuffer) == 0)
      return pow(10.0,ndsplineeval(numu_dsdxdy.get(),xx,centerbuffer,0));
    else
      return 0.0;
  }
  else if (neutype == (int) particleType::MuPlus or neutype == (int) particleType::NuMuBar){
    //int flavor = 4;// muon antineutrino
    //return dsdxdy_(&nuEnergy,&x,&y,&flavor,&target,&cc);///sigma;
    if(tablesearchcenters(numubar_dsdxdy.get(),xx,centerbuffer) == 0)
      return pow(10.0,ndsplineeval(numubar_dsdxdy.get(),xx,centerbuffer,0));
    else
      return 0.0;
  }
  else {
    std::cerr << "CrossSection:CalDDXSPhotoSpline : Bad particle "  << neutype << std::endl;
    exit(1);
  }
}

CrossSectionFromSpline::CrossSectionFromSpline(std::string splinepath, std::string model_name):
  numu_dsdxdy(new splinetable,[](splinetable* t){ splinetable_free(t); delete t;}),
  numubar_dsdxdy(new splinetable,[](splinetable* t){ splinetable_free(t); delete t;})
{
  std::string numu_spl,numubar_spl;

  if( model_name == ""){
    numu_spl = splinepath + "/dsdxdy-numu-N-cc.fits";
    numubar_spl = splinepath + "/dsdxdy-numubar-N-cc.fits";
  } else{
    numu_spl = splinepath + "/dsdxdy-numu-N-cc-"+model_name+".fits";
    numubar_spl = splinepath + "/dsdxdy-numubar-N-cc-"+model_name+".fits";
  }

  // lets photospline :)
  readsplinefitstable(numu_spl.c_str(),numu_dsdxdy.get());
  assert(numu_dsdxdy->ndim>0);

  readsplinefitstable(numubar_spl.c_str(),numubar_dsdxdy.get());
  assert(numubar_dsdxdy->ndim>0);
}

} // end namespace

