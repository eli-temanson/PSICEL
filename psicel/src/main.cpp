
#ifndef main_cpp
#define main_cpp

#include "catima/catima.h"
#include "catima/gwm_integrators.h"
#include "yaml-cpp/yaml.h"
#include "spdlog/spdlog.h"

#include <iostream>
#include <string>

#include "TROOT.h"
#include "TApplication.h"
#include "THashTable.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"


inline void MyFill(THashTable *rootObj,std::string name,int binsx,double minx,double maxx,double valuex,int binsy,double miny,double maxy,double valuey);
inline void MyFill(THashTable *rootObj,std::string name,int binsx,double minx,double maxx,double valuex);


int main(int argc,char *argv[])
{
    TApplication app("app",&argc,argv);

    spdlog::info("Welcome to spdlog!");

    YAML::Node data;
    try
    {
        data = YAML::LoadFile("test_config.yaml");
        spdlog::info("File opened");
    }
    catch (YAML::ParserException& execp)
    {
        spdlog::info("Unable to open yaml file or parse it");
        return 1;
    }

    spdlog::info("value: {0}",data["num"].as<double>());

    catima::Material target = catima::Material({{2,1,4},{12,6,2}}); // CD2
    target.density(0.7987); // g/cm3
    target.thickness(0.000258); // 516 g/cm2

    catima::Material kapton = catima::get_material(catima::material::Kapton);
    kapton.density(1.42);
    kapton.thickness(0.000994);

    catima::Layers ion_chamber;
    catima::Material ic_dl = catima::get_material(catima::material::Isobutane);
    ic_dl.density(0.00011);
    ic_dl.thickness(0.00036);
    ion_chamber.add(ic_dl);

    catima::Material ic_x = catima::get_material(catima::material::Isobutane);
    ic_x.density(0.00011);
    ic_x.thickness(0.00044);
    ion_chamber.add(ic_x);

    catima::Material ic_y = catima::get_material(catima::material::Isobutane);
    ic_y.density(0.00011);
    ic_y.thickness(0.00044);
    ion_chamber.add(ic_y);

    catima::Material ic_de = catima::get_material(catima::material::Isobutane);
    ic_de.density(0.00011);
    ic_de.thickness(0.00088);
    ion_chamber.add(ic_de);
    
    catima::Material ic_e = catima::get_material(catima::material::Isobutane);
    ic_e.density(0.00011);
    ic_e.thickness(0.0022);
    ion_chamber.add(ic_e);

    catima::Projectile light = catima::Projectile(1,1);
    catima::Projectile heavy = catima::Projectile(25,13,13);


    // =========================================================
    catima::Projectile beam1 = catima::Projectile(25,13,13);
    double ebeam = 102.0; 
    // catima::Projectile beam1 = catima::Projectile(24,12,12);
    // double ebeam = 90.57; 
    // catima::Projectile beam1 = catima::Projectile(24,12,11);
    // double ebeam = 76.12;
    
    spdlog::info("CATIMA SIM");
    spdlog::info("eloss in 1/2 target: {0}",catima::integrate_energyloss(beam1(ebeam/beam1.A),target));
    ebeam = ebeam - catima::integrate_energyloss(beam1(ebeam/beam1.A),target);
    spdlog::info("after 1/2 target: {0}", ebeam);
    
    spdlog::info("eloss in target: " , catima::integrate_energyloss(beam1(ebeam/beam1.A),target));
    ebeam = ebeam - catima::integrate_energyloss(beam1(ebeam/beam1.A),target);
    spdlog::info("after target: {0}",ebeam);
 
    spdlog::info("eloss in kapton: {0}",catima::integrate_energyloss(beam1(ebeam/beam1.A),kapton));
    ebeam = ebeam - catima::integrate_energyloss(beam1(ebeam/beam1.A),kapton);
    spdlog::info("after kapton: {0}",ebeam);

    spdlog::info("eloss in ic_dl: {0}",catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_dl));
    ebeam = ebeam - catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_dl);
    spdlog::info("after ic_dl: {0}",ebeam);

    spdlog::info("eloss in ic_x: {0}",catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_x));
    ebeam = ebeam - catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_x);
    spdlog::info("after ic_x: {0}",ebeam);

    spdlog::info("eloss in ic_y: {0}",catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_y));
    ebeam = ebeam - catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_y);
    spdlog::info("after ic_y: {0}",ebeam);

    double efinal = ebeam; // MeV
    double beam_mean = efinal;
    
    spdlog::info("eloss in ic_de: {0}",catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_de));
    ebeam = ebeam - catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_de);
    spdlog::info("after ic_de: {0}",ebeam);

    spdlog::info("eloss in ic_e: {0}",catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_e));
    ebeam = ebeam - catima::integrate_energyloss(beam1(ebeam/beam1.A),ic_e);
    spdlog::info("after ic_e: {0}",ebeam);
    spdlog::info("CATIMA SIM 2 (END)");

    // =========================================================


    // =========================================================
    spdlog::info("CATIMA SIM REVERSE");
    ebeam = efinal;

    spdlog::info("eloss in ic_y: {0}",catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_y));
    ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_y);
    spdlog::info("before y: {0}",ebeam);

    spdlog::info("eloss in ic_x: {0}",catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_x));
    ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_x);
    spdlog::info("before x: {0}",ebeam);

    spdlog::info("eloss in ic_dl: {0}",catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_dl));
    ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_dl);
    spdlog::info("before dl: {0}",ebeam);

    spdlog::info("eloss in kapton: {0}",catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),kapton));
    ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),kapton);
    spdlog::info("before kapton: {0}",ebeam);

    spdlog::info("eloss in target: {0}",catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),target));
    ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),target);
    spdlog::info("before target: {0}",ebeam);

    spdlog::info("eloss in target: {0}",catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),target));
    ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),target);
    spdlog::info("before target: {0}",ebeam);

    spdlog::info("CATIMA SIM REVERSE (END)");


    // =========================================================
    TFile *out_file = new TFile("data/ic_resolution.root","RECREATE");
    THashTable *hist_table = new THashTable();  

    TRandom3 *Rndm = new TRandom3();
    ebeam = beam_mean;// MeV
    target.thickness(0.000516); // 516 g/cm2
    
    for (int i=0; i<100; i++)
    {
    	ebeam = Rndm->Gaus(beam_mean, 2.0/2.355);
    	// x = dx;
	
        ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_y);
        ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_x);
        ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),ic_dl);
        ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),kapton);
        ebeam = ebeam + catima::reverse_integrate_energyloss(beam1(ebeam/beam1.A),target);

    	// fill the histogram and add some beam energy width, 1 MeV
        MyFill(hist_table,"sim_beam_reco",600,0,120, Rndm->Gaus(ebeam, 0.5/2.355) );
    }

    spdlog::info("Output File: {0}",out_file->GetName());
    out_file->Write(out_file->GetName(), TObject::kOverwrite);
    out_file->Close();

    
    return 0;
}


// =================================================================================================


/*2D histogram fill wrapper*/
inline void MyFill(THashTable *rootObj,std::string name,int binsx,double minx,double maxx,double valuex,int binsy,double miny,double maxy,double valuey)
{
  TH2F *histo = (TH2F*) rootObj->FindObject(name.c_str());
  if(histo != NULL) {
    histo->Fill(valuex, valuey);
  } else {
    TH2F *h = new TH2F(name.c_str(), name.c_str(), binsx, minx, maxx, binsy, miny, maxy);
    h->Fill(valuex, valuey);
    rootObj->Add(h);
  }
}
/*1D histogram fill wrapper*/
inline void MyFill(THashTable *rootObj,std::string name,int binsx,double minx,double maxx,double valuex)
{
  TH1F *histo = (TH1F*) rootObj->FindObject(name.c_str());
  if(histo != NULL) {
    histo->Fill(valuex);
  } else {
    TH1F *h = new TH1F(name.c_str(), name.c_str(), binsx, minx, maxx);
    h->Fill(valuex);
    rootObj->Add(h);
  }
}


#endif
