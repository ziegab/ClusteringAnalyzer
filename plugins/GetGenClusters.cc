// -*- C++ -*-
//
// Package:    test/ClusteringAnalyzer
// Class:      GetGenClusters
//
/**\class GetGenClusters GetGenClusters.cc test/ClusteringAnalyzer/plugins/GetGenClusters.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabija Ziemyte
//         Created:  Mon, 19 Aug 2024
//
//

// system include files
#include <memory>
#include <array>
#include <iostream>
#include <cmath>
#include <map>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

/*
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
*/

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "boost/format.hpp"
// #include "ClusteringAnalyzer.h"

//root files
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

//using reco::TrackCollection;

class GetGenClusters : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit GetGenClusters(const edm::ParameterSet&);
  ~GetGenClusters() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  struct hit {
    float eta;
    float phi;
    float x;
    float y;
    int layer;
    float weight;
    int clusterId;
  };

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  TFile* root_file;
  TTree* tree;

  std::map<TString, TH2*> m_Histos2D;

  int nev = 0 ;

  void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight);

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<pat::Photon>> photonToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  std::vector<int> EventsToScan_;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;

#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GetGenClusters::GetGenClusters(const edm::ParameterSet& iConfig):
        photonToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag > ("photons"))),
        genPartToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag > ("genParticles"))),

        caloGeometryToken_(esConsumes()),
        EventsToScan_(iConfig.getParameter<std::vector<int>>("EventsToScan"))
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
  //usesResource("TFileService");
  edm::Service<TFileService> fs;

  TString hname;
  
}

GetGenClusters::~GetGenClusters() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void GetGenClusters::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  TString hname;

  std::vector<hit> rechitsEB, rechitsEEp, rechitsEEm;

  std::vector<float> v_eta;
  std::vector<float> v_phi;
  std::vector<int> v_layer;
  std::vector<float> v_weight;

  std::vector<float> v_x_pEE;
  std::vector<float> v_y_pEE;
  std::vector<float> v_eta_pEE;
  std::vector<float> v_phi_pEE;
  std::vector<int> v_layer_pEE;
  std::vector<float> v_weight_pEE;

  std::vector<float> v_x_mEE;
  std::vector<float> v_y_mEE;
  std::vector<float> v_eta_mEE;
  std::vector<float> v_phi_mEE;
  std::vector<int> v_layer_mEE;
  std::vector<float> v_weight_mEE;

  std::vector<float> v_x_pES1;
  std::vector<float> v_y_pES1;
  std::vector<float> v_eta_pES1;
  std::vector<float> v_phi_pES1;
  std::vector<int> v_layer_pES1;
  std::vector<float> v_weight_pES1;

  std::vector<float> v_x_mES1;
  std::vector<float> v_y_mES1;
  std::vector<float> v_eta_mES1;
  std::vector<float> v_phi_mES1;
  std::vector<int> v_layer_mES1;
  std::vector<float> v_weight_mES1;

  std::vector<float> v_x_pES2;
  std::vector<float> v_y_pES2;
  std::vector<float> v_eta_pES2;
  std::vector<float> v_phi_pES2;
  std::vector<int> v_layer_pES2;
  std::vector<float> v_weight_pES2;

  std::vector<float> v_x_mES2;
  std::vector<float> v_y_mES2;
  std::vector<float> v_eta_mES2;
  std::vector<float> v_phi_mES2;
  std::vector<int> v_layer_mES2;
  std::vector<float> v_weight_mES2;

  nev++;

  //======= PHOTONS =======
  edm::Handle<std::vector<pat::Photon> > gamma;
  iEvent.getByToken(photonToken_, gamma);
  std::cout << "PAT photons detected: " << gamma->size() << std::endl;

  for(size_t i = 0; i < gamma->size(); ++ i) {
    const pat::Photon & p = (*gamma)[i];
    std::cout << "Found a photon with status, pt, eta, phi, energy : " << p.status() << ", " << p.pt() << ", " << p.eta() << ", " << p.phi() << ", " << p.energy() << std::endl;
    // hname = Form("PATPhoton_event%i", nev);
    // FillHist2D(hname, p.eta()*(85/1.4), p.phi()*(180/3.14159265), p.energy()) ;
  }

  //======= GEN PARTICLES =======
  edm::Handle<std::vector<pat::PackedGenParticle> > genPart;
  iEvent.getByToken(genPartToken_, genPart);
  std::cout << "Gen particles collection size: " << genPart->size() << std::endl;

 // PdgId for photon is 22, final state particles have status == 1
  
  for(size_t i = 0; i < genPart->size(); ++ i) {
    const pat::PackedGenParticle & p = (*genPart)[i];
    if ((p.pdgId() == 22) & (p.status()==1) & (p.pt()>10)){
      std::cout << "Found a photon with status, pt, eta, phi, energy : " << p.status() << ", " << p.pt() << ", " << p.eta() << ", " << p.phi() << ", " << p.energy() << std::endl;
      // hname = Form("genPhoton_event%i", nev);
      // FillHist2D(hname, p.eta(), p.phi(), p.energy()) ;
    }
  }


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void GetGenClusters::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void GetGenClusters::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GetGenClusters::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GetGenClusters);
