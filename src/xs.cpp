#include "xs.h"

namespace I3PGE {

double I3PGECrossSection::DoubleDifferentialCrossSection(double E, double x, double y,
    nusquids::NeutrinoCrossSections::NeutrinoFlavor flavor,
    nusquids::NeutrinoCrossSections::NeutrinoType neutype,
    nusquids::NeutrinoCrossSections::Current current) const {

  if (current != nusquids::NeutrinoCrossSections::Current::CC or
      flavor != nusquids::NeutrinoCrossSections::muon)
    throw std::runtime_error("CrossSection::DoubleDifferentialCrossSection: only nu_mu CC cross section spline exist");

  int centerbuffer[3];
  double xx[3];

  xx[0] = log10(E);
  xx[1] = log10(x);
  xx[2] = log10(y);

  if (neutype == nusquids::NeutrinoCrossSections::NeutrinoType::neutrino){
    //int flavor = 3;// muon neutrino
    //return dsdxdy_(&nuEnergy,&x,&y,&flavor,&target,&cc);///sigma;
    if(tablesearchcenters(numu_dsdxdy.get(),xx,centerbuffer) == 0)
      return pow(10.0,ndsplineeval(numu_dsdxdy.get(),xx,centerbuffer,0));
    else
      return 0.0;
  }
  else if (neutype == nusquids::NeutrinoCrossSections::NeutrinoType::antineutrino){
    //int flavor = 4;// muon antineutrino
    //return dsdxdy_(&nuEnergy,&x,&y,&flavor,&target,&cc);///sigma;
    if(tablesearchcenters(numubar_dsdxdy.get(),xx,centerbuffer) == 0)
      return pow(10.0,ndsplineeval(numubar_dsdxdy.get(),xx,centerbuffer,0));
    else
      return 0.0;
  }
  else {
    throw std::runtime_error("CrossSection::DoubleDifferentialCrossSection: unexpected neutype.");
  }
}

I3PGECrossSection::I3PGECrossSection(std::string splinepath, std::string model_name):
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

