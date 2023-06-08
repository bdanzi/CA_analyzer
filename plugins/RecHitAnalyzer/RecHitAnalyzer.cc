// -*- C++ -*-
//
// Package:    RecHits/RecHitAnalyzer
// Class:      RecHitAnalyzer
//
//
// Original Author:  Brunella D'Anzi
//         Created:  Mon, 28 May 2023 17:16:39 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class RecHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit RecHitAnalyzer(const edm::ParameterSet&);
  ~RecHitAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void plotClusterShape(const SiStripCluster& cluster, const int sequentialId);
  void plotRecHitPosition(const float x, const float y);
  void plotRecHitGlobalPosition(const float xGlobal, const float yGlobal, const float zGlobal);
  void plotRecHitErrors(const float x, const float y, const float xErr, const float yErr);
  void plotDetIdMapping( std::map<uint32_t, int> detIdToSequentialNumber_);
  void plotProjectionsGlobalPosition();
  void plotProjectionsClusterShape();
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster>> clusterToken_;
  edm::EDGetTokenT<SiStripRecHit2DCollection> recHitToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  std::map<uint32_t, int> detIdToSequentialNumber_;
  int sequentialNumber_;
  unsigned int  maxId;
  unsigned int minId;
  TH2F* recHitPositionHist_;
  TH3F* clusterShapeHist_;
  TH2F* recHitErrorsHist_;
  TH2F* detIdMappingHist_;
  TH3F* recHitPositionGlobalHist_;
};


//
// constructors and destructor
//
RecHitAnalyzer::RecHitAnalyzer(const edm::ParameterSet& iConfig)
  :clusterToken_(consumes<edmNew::DetSetVector<SiStripCluster>>(iConfig.getParameter<edm::InputTag>("recHits"))),
   recHitToken_(consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("recHitCollection"))),
   geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()){
  //now do what ever initialization is needed
  sequentialNumber_ = 0;
  maxId= 0;
  minId= 10000000;
  recHitPositionHist_ = nullptr;
  clusterShapeHist_ = nullptr;
  recHitErrorsHist_ = nullptr;
  detIdMappingHist_ = nullptr;
  recHitPositionGlobalHist_ = nullptr;
}

RecHitAnalyzer::~RecHitAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  edm::Handle<edmNew::DetSetVector<SiStripCluster>> clusterHandle;
  iEvent.getByToken(clusterToken_, clusterHandle);
  edm::ESHandle<TrackerGeometry> theTrackerGeometry = iSetup.getHandle(geomToken_);
  edm::Handle<SiStripRecHit2DCollection> recHitHandle;
  iEvent.getByToken(recHitToken_, recHitHandle);
  // Loop over each individual SiStripRecHit2DCollection pair
      for (const auto& detSetPair : *recHitHandle) {
	const auto& detSet = detSetPair;
    	for (const auto& recHit : detSet) {
	const StripGeomDetUnit* theStripDet = dynamic_cast<const StripGeomDetUnit*>(theTrackerGeometry->idToDet(recHit.geographicalId()));
	// Access information from SiStripRecHit2D                                                                     
	float x = recHit.localPosition().x();
	float y = recHit.localPosition().y();
	float xErr = sqrt(recHit.localPositionError().xx());
	float yErr = sqrt(recHit.localPositionError().yy());
	//auto detId_geographicalId = recHit.geographicalId();
	//int detId = detId_geographicalId.rawId();
	int detId = theStripDet->geographicalId().rawId();
	int sequentialId = 0;                                                                                                                                                               float xGlobal = theStripDet->toGlobal(recHit.localPosition()).x();
	float yGlobal = theStripDet->toGlobal(recHit.localPosition()).y();
	float zGlobal = theStripDet->toGlobal(recHit.localPosition()).z();

	//LocalPoint localPos = recHit.localPosition();
	//LocalError localPosError = recHit.localPositionError();

	// Check if the detId is already mapped, if not, assign a sequential number
        if (detIdToSequentialNumber_.find(detId) == detIdToSequentialNumber_.end()) {                                                                                          
          sequentialId = sequentialNumber_++; // Get the current size as the sequential number                                                                          
	  detIdToSequentialNumber_[detId] = sequentialId;                                                      
        }                                                                                        
	else{
	  sequentialId = detIdToSequentialNumber_[detId];
	}
	std::cout << "DetId of each RecHits: " << detId << " Sequential Number: " << sequentialId  << std::endl;
	std::cout << "Local position (x, y): " << x << ", " << y << std::endl;
	std::cout << "Global position (x, y, z): " << xGlobal << ", " << yGlobal << ", " << zGlobal <<std::endl;
	std::cout << "Local position errors (x, y): " << xErr << ", " << yErr << std::endl;
	
	// Plot the recHit position, cluster shape, recHit errors, and DetId mapping                                                                                                                     
	plotRecHitPosition(x, y);
	plotRecHitGlobalPosition(xGlobal, yGlobal, zGlobal);
	plotClusterShape(*recHit.cluster(),sequentialId);
	plotRecHitErrors(x, y, xErr, yErr);                                                                                                                                                                     
      }
    }
      plotDetIdMapping(detIdToSequentialNumber_);
}


void RecHitAnalyzer::plotProjectionsGlobalPosition() {
  TCanvas* canvas_PlotGlobalPositionProjections = new TCanvas("PlotGlobalPositionProjections", "Projections", 1200, 400);
  canvas_PlotGlobalPositionProjections->Divide(3, 1); // Divide the canvas into 3 sub-pads
  // Plot the projection on the XY plane
  canvas_PlotGlobalPositionProjections->cd(1);
  TH2D* projXY = dynamic_cast<TH2D*>(recHitPositionGlobalHist_->Project3D("xy"));
  projXY->Draw("COLZ");

  // Plot the projection on the XZ plane
  canvas_PlotGlobalPositionProjections->cd(2);
  TH2D* projXZ = dynamic_cast<TH2D*>(recHitPositionGlobalHist_->Project3D("xz"));
  projXZ->Draw("COLZ");

  // Plot the projection on the YZ plane
  canvas_PlotGlobalPositionProjections->cd(3);
  TH2D* projYZ = dynamic_cast<TH2D*>(recHitPositionGlobalHist_->Project3D("yz"));
  projYZ->Draw("COLZ");

  // Save the canvas as an image file
  canvas_PlotGlobalPositionProjections->SaveAs("projectionsGlobalPosition.png");

  // Clean up memory
  delete canvas_PlotGlobalPositionProjections;
  delete projXY;
  delete projXZ;
  delete projYZ;
}

void RecHitAnalyzer::plotProjectionsClusterShape() {
  TCanvas* canvas = new TCanvas("PlotClusterShapeProjections", "Cluster Shape Projections", 1200, 400);
  canvas->Divide(3, 1); // Divide the canvas into 3 sub-pads                                                                                                                         
  // Plot the projection on the XY plane                                                                                                                                             
  canvas->cd(1);
  TH2D* projXY = dynamic_cast<TH2D*>(clusterShapeHist_->Project3D("xy"));
  projXY->Draw("COLZ");

  // Plot the projection on the XZ plane                                                                                                                                             
  canvas->cd(2);
  TH2D* projXZ = dynamic_cast<TH2D*>(clusterShapeHist_->Project3D("xz"));
  projXZ->Draw("COLZ");

  // Plot the projection on the YZ plane                                                                                                                                             
  canvas->cd(3);
  TH2D* projYZ = dynamic_cast<TH2D*>(clusterShapeHist_->Project3D("yz"));
  projYZ->Draw("COLZ");

  // Save the canvas as an image file                                                                                                                                                
  canvas->SaveAs("projectionsClusterShape.png");

  // Clean up memory                                                                                                                                                                 
  delete canvas;
  delete projXY;
  delete projXZ;
  delete projYZ;
}
void RecHitAnalyzer::plotClusterShape(const SiStripCluster& cluster, const int sequentialId) {
  // Plot the cluster shape (3D plot)
      auto  amplitudes = cluster.amplitudes();
      const int clusterSize = amplitudes.size();
      const int clusterCharge = cluster.charge();
      for (size_t i = 0; i < amplitudes.size(); ++i) {
	//std::cout << "Strip Number " << i << std::endl;
	clusterShapeHist_->Fill(clusterCharge, amplitudes[i], clusterSize);
      }
    }
    
    void RecHitAnalyzer::plotRecHitPosition(const float x, const float y) {
      
      // Plot the recHit position (2D plot)
      recHitPositionHist_->Fill(x, y);
      
    }

void RecHitAnalyzer::plotRecHitGlobalPosition(const float xGlobal, const float yGlobal, const float zGlobal) {

  // Plot the recHit position (3D plot)                                                                                                                                              
  recHitPositionGlobalHist_->Fill(xGlobal, yGlobal, zGlobal);

}


    void RecHitAnalyzer::plotRecHitErrors(const float x, const float y, const float xErr, const float yErr) {
      
      // Plot the recHit errors (2D plot)
      recHitErrorsHist_->Fill(xErr, yErr);
    
    }

    void RecHitAnalyzer::plotDetIdMapping(std::map<uint32_t, int> detIdToSequentialNumber_) {
      // Plot the DetId mapping (2D plot)
      
      for (const auto& pair : detIdToSequentialNumber_) {
	detIdMappingHist_->Fill(pair.second, pair.first);
	if(pair.first > maxId) maxId = pair.first;
	if(pair.first < minId) minId = pair.first;
      }
      
    }
    // ------------ method called once each job just before starting event loop  ------------
    void RecHitAnalyzer::beginJob() {
      
      detIdMappingHist_ = new TH2F("detIdMappingHist", "DetId Mapping", sequentialNumber_, 0, sequentialNumber_, maxId - minId, minId, maxId);
      detIdMappingHist_->GetXaxis()->SetTitle("Sequential Number");
      detIdMappingHist_->GetYaxis()->SetTitle("DetId");

      recHitErrorsHist_ = new TH2F("recHitErrorsHist", "RecHit Errors", 100, -0.001, 0.2, 500, -0.5, 30.0);
      recHitErrorsHist_->GetXaxis()->SetTitle("x (mm)");
      recHitErrorsHist_->GetYaxis()->SetTitle("y (mm)");
      
      recHitPositionHist_ = new TH2F("recHitPositionHist", "RecHit Position", 100, -5.5, 5.5, 100, -1.0, 1.0);
      recHitPositionHist_->GetXaxis()->SetTitle("x (mm)");
      recHitPositionHist_->GetYaxis()->SetTitle("y (mm)");
      
      recHitPositionGlobalHist_ = new TH3F("recHitPositionGlobalHist", "Global Position Rec Strip Hits", 400, -200.0, 200.0, 400, -200.0, 200.0, 600, -700.0, 700.0);
      recHitPositionGlobalHist_->GetXaxis()->SetTitle("x (mm)");
      recHitPositionGlobalHist_->GetYaxis()->SetTitle("y (mm)");
      recHitPositionGlobalHist_->GetZaxis()->SetTitle("z (mm)");

      clusterShapeHist_ = new TH3F("clusterShapeHist", "Cluster Shape", 1500, 0, 1500, 256, 0, 256, 10, 0, 10);
      clusterShapeHist_->GetXaxis()->SetTitle("Cluster Charge  (ehp)");
      clusterShapeHist_->GetYaxis()->SetTitle("ADC Counts");
      clusterShapeHist_->GetZaxis()->SetTitle("Cluster Size");
      clusterShapeHist_->GetXaxis()->SetTitleOffset(1.5); // Adjust the offset for the X-axis label
      clusterShapeHist_->GetYaxis()->SetTitleOffset(1.5); // Adjust the offset for the Y-axis label
      clusterShapeHist_->GetZaxis()->SetTitleOffset(1.5); // Adjust the offset for the Z-axis label 

    }
    
    // ------------ method called once each job just after ending the event loop  ------------
    void RecHitAnalyzer::endJob() {
      TCanvas* canvas_recHitPosition = new TCanvas("canvas_recHitPosition", "RecHit Position", 800, 600);
      TCanvas* canvas_detIdMapping = new TCanvas("canvas_detIdMapping", "DetId Mapping", 1000, 600);
      TCanvas* canvas_recHitErrors = new TCanvas("canvas_recHitErrors", "RecHit Errors", 800, 600);
      TCanvas* canvas_clusterShapeHist = new TCanvas("canvas_clusterShapeHist", "Cluster Shape", 1000, 600);
      TCanvas* canvas_recHitPositionGlobalHist = new TCanvas("canvas_recHitPositionGlobalHist", "Global Position Strip Rec Hits", 1000, 600);
      plotProjectionsGlobalPosition();
      plotProjectionsClusterShape();

      canvas_recHitPosition->cd();
      // Draw the histogram on the canvas                                                                                                                                          
      recHitPositionHist_->Draw("hist");
      // Save the canvas as a PNG image                                                                                                                                        
      canvas_recHitPosition->SaveAs("recHitPosition.png");                                                                                                                                delete canvas_recHitPosition;

      // Draw the histogram on the canvas                                                    
      canvas_recHitPositionGlobalHist->cd();                                                                                           
      recHitPositionGlobalHist_->Draw();
      canvas_recHitPositionGlobalHist->SaveAs("recHitGlobalPosition.png");
      delete canvas_recHitPositionGlobalHist;

      canvas_detIdMapping->cd();
      detIdMappingHist_->Draw("hist");
      // Save the canvas as a PNG image                                                                                                                                              
      canvas_detIdMapping->SaveAs("detIdMapping.png");
      delete canvas_detIdMapping;                                                                                                                                                                 
      // Save the canvas as a PNG image                                                                                                                                                                 
      canvas_recHitErrors->cd();
      recHitErrorsHist_->Draw("hist");
      canvas_recHitErrors->SaveAs("recHitErrors.png");
      delete canvas_recHitErrors;
      // Save the canvas as a PNG image                                                                                                                                                                 
      canvas_clusterShapeHist->cd();
      TPaveStats* stats = dynamic_cast<TPaveStats*>(clusterShapeHist_->GetListOfFunctions()->FindObject("stats"));
      // Relocate the stats box
      if (stats) {
	stats->SetX1NDC(0.25); // Set the new X-coordinate of the lower edge
	stats->SetX2NDC(0.55); // Set the new X-coordinate of the upper edge
	stats->SetY1NDC(0.75); // Set the new Y-coordinate of the lower edge
	stats->SetY2NDC(0.95); // Set the new Y-coordinate of the upper edge
	stats->Draw();
	canvas_clusterShapeHist->Update();
	}
      clusterShapeHist_->Draw();
      canvas_clusterShapeHist->SaveAs("clusterShapeHist.png"); 
      delete canvas_clusterShapeHist;

      TFile output("RecHitAnalyzer.root", "RECREATE");
      recHitPositionHist_->Write();
      recHitPositionGlobalHist_->Write();
      clusterShapeHist_->Write();
      recHitErrorsHist_->Write();
      detIdMappingHist_->Write();
      output.Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void RecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalyzer);
