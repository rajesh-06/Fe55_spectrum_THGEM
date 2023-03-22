//# include <execution>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <vector>
#include<algorithm>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1.h>
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Setup the gas.
  MediumMagboltz gas("ar", 70., "co2", 30.);
  gas.SetTemperature(293.15);//in K
  gas.SetPressure(760.);//in torr
  gas.Initialise(true);  
  // Set the Penning transfer efficiency.
  constexpr double rPenning = 0.51;
  constexpr double lambdaPenning = 0.;
  gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");

  // Load the field map.
  ComponentAnsys123 fm;
  fm.Initialise("ELIST.lis", "NLIST.lis", "MPLIST.lis", "PRNSOL.lis", "mm");
  fm.EnableMirrorPeriodicityX();
  fm.EnableMirrorPeriodicityY();
  fm.PrintRange();

  // Associate the gas with the corresponding field map material.
  fm.SetGas(&gas); 
  fm.PrintMaterials();
  // fm.Check();

  // Dimensions of the GEM [cm]
  constexpr double pitch = 0.08;

  ViewField fieldView;
  ViewFEMesh meshView;
  constexpr bool plotField = false;
  if (plotField) {
    fieldView.SetComponent(&fm);
    // Set the normal vector of the viewing plane (xz plane).
    fieldView.SetPlane(0, -1, 0, 0, 0, 0);
    // Set the plot limits in the current viewing plane.
    fieldView.SetArea(-0.5 * pitch, -0.25, 0.5 * pitch, 0.5);
    fieldView.SetVoltageRange(0., -2460.);
    TCanvas* cf = new TCanvas("cf", "", 600, 600);
    cf->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf);
    fieldView.PlotContour();

    meshView.SetArea(-0.5 * pitch, -0.25, 0.5 * pitch, 0.5); 
    meshView.SetCanvas(cf);
    meshView.SetComponent(&fm);
    meshView.SetPlane(0, -1, 0, 0, 0, 0);
    meshView.SetFillMesh(true);
    meshView.SetColor(2, kGray);
    meshView.Plot(true);
  }

  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(&fm);
  sensor.SetArea(-5 * pitch, -5 * pitch, -0.25,
                  5 * pitch,  5 * pitch,  0.5);

  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);

  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.SetDistanceSteps(2.e-4);

  ViewDrift driftView;
  const double source_energy = 5900;
  const int IE = 26;//ionisation energy of argon in eV
  const int KE = 1;
  double TE_pIon = IE+KE;
  int nEvents = 1000;
  int nPrimaryIonisations = source_energy/ TE_pIon;
    std::ofstream file;
    file.open("../gain.txt");
  //std::vector<int> indexes
  std::vector<double> avg_gain;
  TCanvas C1;
  for(unsigned int event=0; event<=nEvents; ++event){
    std::vector<double> gain_vec;
  	for (unsigned int i = 0; i < nPrimaryIonisations; ++i) { 
  		    std::cout << i+1 << "/" << nPrimaryIonisations;
    		// Randomize the initial position. 
    		const double x0 = -0.5 * pitch + RndmUniform() * pitch;
    		const double y0 = -0.5 * pitch + RndmUniform() * pitch;
    		const double z0 = 0.48; 
    		const double t0 = 0.;
    		const double e0 = 2*RndmUniform();//initial KE of electron (averaged over 1eV)
    		aval.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    		int ne = 0, ni = 0;
    		aval.GetAvalancheSize(ne, ni);

        double gain = aval.GetNumberOfElectronEndpoints();
        gain_vec.push_back(gain);//->Fill(gain);
        std::cout<<" Gain: "<<gain<< "\n";
  	}
    double average_gain = 1+accumulate(gain_vec.begin(), gain_vec.end(), 0)/nPrimaryIonisations;
    avg_gain.push_back(average_gain);
    //delete gain_vec;
    std::cout<<"Average gain of this event: "<<average_gain<<std::endl;
    file<<average_gain<<std::endl;
  }
  
  file.close();
  double max_gain = *std::max_element(avg_gain.begin(), avg_gain.end());
  printf("maxgain=%d",max_gain);
  TH1I *hist = new TH1I("Gain","Gain",40,0,max_gain+1);
  for (int gn:avg_gain) hist->Fill(gn); 
  hist->Draw();


  app.Run();
}

