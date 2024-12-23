// -*- C++ -*-
//
// Package:    test/ClusteringAnalyzer
// Class:      ClusteringAnalyzer
//
/**\class ClusteringAnalyzer ClusteringAnalyzer.cc test/ClusteringAnalyzer/plugins/ClusteringAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Garvita Agarwal
//         Created:  Sun, 24 Sep 2023 20:45:30 GMT
//
//

// system include files
#include <memory>
#include <array>
#include <iostream>
#include <cmath>
#include <map>

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
// #include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
#include "ClusteringAnalyzerAOD.h"

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

class ClusteringAnalyzerAOD : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ClusteringAnalyzerAOD(const edm::ParameterSet&);
  ~ClusteringAnalyzerAOD() override;

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
  //edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> ecalRecHitsEBToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> ecalRecHitsEEToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> ecalRecHitsESToken_;
  // edm::EDGetTokenT<vector<pat::Photon>> photonToken_;
  edm::EDGetTokenT<vector<reco::GenParticle>> genPartToken_;

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
ClusteringAnalyzerAOD::ClusteringAnalyzerAOD(const edm::ParameterSet& iConfig):
	//tracksToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
	ecalRecHitsEBToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag > ("ecalRechitsEB"))),
	ecalRecHitsEEToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag > ("ecalRechitsEE"))),
	ecalRecHitsESToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag > ("ecalRechitsES"))),
        // photonToken_(consumes<vector<pat::Photon>>(iConfig.getParameter<edm::InputTag > ("slimmedPhotons"))),
        genPartToken_(consumes<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag > ("genParticles"))),

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
  // const int PHI_BINS = 363;
  // const int ETA_BINS = 192;
  // const int ETA_BINS_EB = 170;

  // double philim = 3.1494;
  // double etabins[ETA_BINS+1] = 
  //  {-3.0, -2.65, -2.5, -2.32, -2.17, -2.04, -1.93, -1.83, -1.74, -1.65, -1.57,
  //   -1.479, -1.4616, -1.4442, -1.4268, -1.4094, -1.392, -1.3746, -1.3572, -1.3398, -1.3224, -1.305, -1.2876, -1.2702, -1.2528,
  //   -1.2354, -1.218, -1.2006, -1.1832, -1.1658, -1.1484, -1.131, -1.1136, -1.0962, -1.0788, -1.0614, -1.044, -1.0266, -1.0092,
  //   -0.9918, -0.9744, -0.957, -0.9396, -0.9222, -0.9048, -0.8874, -0.87, -0.8526, -0.8352, -0.8178, -0.8004, -0.783, -0.7656,
  //   -0.7482, -0.7308, -0.7134, -0.696, -0.6786, -0.6612, -0.6438, -0.6264, -0.609, -0.5916, -0.5742, -0.5568, -0.5394, -0.522,
  //   -0.5046, -0.4872, -0.4698, -0.4524, -0.435, -0.4176, -0.4002, -0.3828, -0.3654, -0.348, -0.3306, -0.3132, -0.2958, -0.2784,
  //   -0.261, -0.2436, -0.2262, -0.2088, -0.1914, -0.174, -0.1566, -0.1392, -0.1218, -0.1044, -0.087, -0.0696, -0.0522, -0.0348,
  //   -0.0174, 0., 0.0174, 0.0348, 0.0522, 0.0696, 0.087, 0.1044, 0.1218, 0.1392, 0.1566, 0.174, 0.1914, 0.2088, 0.2262, 0.2436,
  //   0.261, 0.2784, 0.2958, 0.3132, 0.3306, 0.348, 0.3654, 0.3828, 0.4002, 0.4176, 0.435, 0.4524, 0.4698, 0.4872, 0.5046, 0.522,
  //   0.5394, 0.5568, 0.5742, 0.5916, 0.609, 0.6264, 0.6438, 0.6612, 0.6786, 0.696, 0.7134, 0.7308, 0.7482, 0.7656, 0.783, 0.8004,
  //   0.8178, 0.8352, 0.8526, 0.87, 0.8874, 0.9048, 0.9222, 0.9396, 0.957, 0.9744, 0.9918, 1.0092, 1.0266, 1.044, 1.0614, 1.0788,
  //   1.0962, 1.1136, 1.131, 1.1484, 1.1658, 1.1832, 1.2006, 1.218, 1.2354, 1.2528, 1.2702, 1.2876, 1.305, 1.3224, 1.3398, 1.3572,
  //   1.3746, 1.392, 1.4094, 1.4268, 1.4442, 1.4616, 1.479, 1.57, 1.65, 1.74, 1.83, 1.93, 2.04, 2.17, 2.32, 2.50, 2.65, 3.0};

  // double etabins_EB[ETA_BINS+1] = 
  //  {-1.479, -1.4616, -1.4442, -1.4268, -1.4094, -1.392, -1.3746, -1.3572, -1.3398, -1.3224, -1.305, -1.2876, -1.2702, -1.2528,
  //   -1.2354, -1.218, -1.2006, -1.1832, -1.1658, -1.1484, -1.131, -1.1136, -1.0962, -1.0788, -1.0614, -1.044, -1.0266, -1.0092,
  //   -0.9918, -0.9744, -0.957, -0.9396, -0.9222, -0.9048, -0.8874, -0.87, -0.8526, -0.8352, -0.8178, -0.8004, -0.783, -0.7656,
  //   -0.7482, -0.7308, -0.7134, -0.696, -0.6786, -0.6612, -0.6438, -0.6264, -0.609, -0.5916, -0.5742, -0.5568, -0.5394, -0.522,
  //   -0.5046, -0.4872, -0.4698, -0.4524, -0.435, -0.4176, -0.4002, -0.3828, -0.3654, -0.348, -0.3306, -0.3132, -0.2958, -0.2784,
  //   -0.261, -0.2436, -0.2262, -0.2088, -0.1914, -0.174, -0.1566, -0.1392, -0.1218, -0.1044, -0.087, -0.0696, -0.0522, -0.0348,
  //   -0.0174, 0., 0.0174, 0.0348, 0.0522, 0.0696, 0.087, 0.1044, 0.1218, 0.1392, 0.1566, 0.174, 0.1914, 0.2088, 0.2262, 0.2436,
  //   0.261, 0.2784, 0.2958, 0.3132, 0.3306, 0.348, 0.3654, 0.3828, 0.4002, 0.4176, 0.435, 0.4524, 0.4698, 0.4872, 0.5046, 0.522,
  //   0.5394, 0.5568, 0.5742, 0.5916, 0.609, 0.6264, 0.6438, 0.6612, 0.6786, 0.696, 0.7134, 0.7308, 0.7482, 0.7656, 0.783, 0.8004,
  //   0.8178, 0.8352, 0.8526, 0.87, 0.8874, 0.9048, 0.9222, 0.9396, 0.957, 0.9744, 0.9918, 1.0092, 1.0266, 1.044, 1.0614, 1.0788,
  //   1.0962, 1.1136, 1.131, 1.1484, 1.1658, 1.1832, 1.2006, 1.218, 1.2354, 1.2528, 1.2702, 1.2876, 1.305, 1.3224, 1.3398, 1.3572,
  //   1.3746, 1.392, 1.4094, 1.4268, 1.4442, 1.4616, 1.479};

  for (unsigned int i=0; i < EventsToScan_.size(); i++){
    // ===== event maps ====== //
    // hname = Form("genPhoton_event%i", EventsToScan_[i]) ;
    // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , ETA_BINS, etabins, PHI_BINS, -1*philim, philim);
    //  hname = Form("rechitMap_event%i", EventsToScan_[i]) ;
    // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , ETA_BINS, etabins, PHI_BINS, -1*philim, philim);
    // hname = Form("clusteredMap_event%i", EventsToScan_[i]) ;
    // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , ETA_BINS, etabins, PHI_BINS, -1*philim, philim);
    // // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 171, -1.479, 1.479, 361, -3.15, 3.15);

    // ===== rechit maps for detector geometry ====== //
    hname = Form("rechitMapEB_event%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 170, -85, 85, 360, 0, 360);
    
    hname = Form("rechitMapEEp_event%i", EventsToScan_[i]) ;
    // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 101, -50, 50, 101, -50, 50);
     m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 100, 0, 100, 100, 0, 100);
    hname = Form("rechitMapEEm_event%i", EventsToScan_[i]) ;
    // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 101, -50, 50, 101, -50, 50);
     m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 100, 0, 100, 100, 0, 100);

    hname = Form("rechitMapES1p_event%i", EventsToScan_[i]) ;
    //m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 2520, -20, 20, 2520, -20, 20);
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 40, 0, 40, 1280, 0, 1280);
    hname = Form("rechitMapES1m_event%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 40, 0, 40, 1280, 0, 1280);

    hname = Form("rechitMapES2p_event%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1280, 0, 1280, 40, 0, 40);
    hname = Form("rechitMapES2m_event%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1280, 0, 1280, 40, 0, 40);


    hname = Form("clusteredMapEB_event%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 170, -85, 85, 360, 0, 360);

    hname = Form("clusteredMapEEp_event%i", EventsToScan_[i]) ;
    // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 101, -50, 50, 101, -50, 50);
     m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 100, 0, 100, 100, 0, 100);
    hname = Form("clusteredMapEEm_event%i", EventsToScan_[i]) ;
    // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 101, -50, 50, 101, -50, 50);
     m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 100, 0, 100, 100, 0, 100);
    
    hname = Form("clusteredMapES1p_event%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 40, 0, 40, 1280, 0, 1280);
    hname = Form("clusteredMapES1m_event%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 40, 0, 40, 1280, 0, 1280);

    hname = Form("clusteredMapES2p_event%i", EventsToScan_[i]) ;
    //m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 2520, -20, 20, 2520, -20, 20);
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1280, 0, 1280, 40, 0, 40);
    hname = Form("clusteredMapES2m_event%i", EventsToScan_[i]) ;
    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1280, 0, 1280, 40, 0, 40);
  }
  
}

ClusteringAnalyzerAOD::~ClusteringAnalyzerAOD() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void ClusteringAnalyzerAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
  // edm::Handle<std::vector<pat::Photon> > gamma;
  // iEvent.getByToken(photonToken_, gamma);


  //======= GEN PARTICLES =======
  edm::Handle<std::vector<reco::GenParticle> > genPart;
  iEvent.getByToken(genPartToken_, genPart);
  std::cout << "EB Gen particles collection size: " << genPart->size() << std::endl;

 // PdgId for photon is 22, final state particles have status == 1
  
  for(size_t i = 0; i < genPart->size(); ++ i) {
    const reco::GenParticle & p = (*genPart)[i];
    if ((p.pdgId() == 22) & (p.status()==1) & (p.pt()>10)){
      std::cout << "Found a photon with status, pt, eta, phi, energy : " << p.status() << ", " << p.pt() << ", " << p.eta() << ", " << p.phi() << ", " << p.energy() << std::endl;
      // hname = Form("genPhoton_event%i", nev);
      // FillHist2D(hname, p.eta(), p.phi(), p.energy()) ;
    }
  }

  //======= ECAL BARREL ==========
  edm::ESHandle<CaloGeometry> caloGeometry = iSetup.getHandle(caloGeometryToken_);
  const CaloSubdetectorGeometry* EBgeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);

  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> handle_RecHitsEB;
  iEvent.getByToken(ecalRecHitsEBToken_, handle_RecHitsEB);
  const EcalRecHitCollection* ecalRecHitsEB = handle_RecHitsEB.product();
  std::cout << "rechit collection size: " << ecalRecHitsEB->size() << std::endl;

  if (ecalRecHitsEB->size() >0){
    for(const EcalRecHit &iEB : *ecalRecHitsEB) {
      EBDetId detID(iEB.id());
      std::cout << "DetID: " << detID << std::endl;
      auto cell = EBgeom->getGeometry(detID);
      float eta, phi, E;
      int layer = 0;
      eta = cell->getPosition().eta(); phi = cell->getPosition().phi(); E = iEB.energy();
      std::cout << "eta: "<< eta << " phi: "<< phi <<" Energy: " << E << std::endl;
      std::cout << "ieta: "<< detID.ieta() << " iphi: "<< detID.iphi() <<" Energy: " << E << std::endl;

      v_eta.push_back(detID.ieta());
      v_phi.push_back(detID.iphi());
      v_layer.push_back(layer);
      v_weight.push_back(E);

      hname = Form("rechitMapEB_event%i", nev);
      FillHist2D(hname, detID.ieta(), detID.iphi(), E);

      // hname = Form("rechitMap_event%i", nev);
      // FillHist2D(hname, eta, phi, E);
    }
  }

  //======= ECAL ENDCAP ==========
  const CaloSubdetectorGeometry* EEgeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap);

  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> handle_RecHitsEE;
  iEvent.getByToken(ecalRecHitsEEToken_, handle_RecHitsEE);
  const EcalRecHitCollection* ecalRecHitsEE = handle_RecHitsEE.product();
  std::cout << "EE rechit collection size: " << ecalRecHitsEE->size() << std::endl;

  if (ecalRecHitsEE->size() >0){
    for(const EcalRecHit &iEE : *ecalRecHitsEE) {
      EEDetId detID(iEE.id());
      std::cout << "DetID: " << detID << std::endl;
      auto cell = EEgeom->getGeometry(detID);
      float eta, phi, E;
      int layer = 0;
      eta = cell->getPosition().eta(); phi = cell->getPosition().phi(); E = iEE.energy();
      std::cout << "ix: "<< detID.ix() << " iy: "<< detID.iy() << std::endl;
      std::cout << "x: "<< cell->getPosition().x() << " y: "<< cell->getPosition().y() <<" z: " << cell->getPosition().z() << std::endl;
      std::cout << "eta: "<< cell->getPosition().eta() << " phi: "<< cell->getPosition().phi() <<" Energy: " << E << std::endl;
      if ((eta>=0)&&(E>1.0e-10)){
        v_x_pEE.push_back(detID.ix());
        v_y_pEE.push_back(detID.iy());
        v_eta_pEE.push_back(eta);
        v_phi_pEE.push_back(phi);
        v_layer_pEE.push_back(layer);
        v_weight_pEE.push_back(E);

        hname = Form("rechitMapEEp_event%i", nev);
        FillHist2D(hname, detID.ix(), detID.iy(), E);
      }
      else if ((eta<0)&&(E>1.0e-10)){
        v_x_mEE.push_back(detID.ix());
        v_y_mEE.push_back(detID.iy());
        v_eta_mEE.push_back(eta);
        v_phi_mEE.push_back(phi);
        v_layer_mEE.push_back(layer);
        v_weight_mEE.push_back(E);

        hname = Form("rechitMapEEm_event%i", nev);
        FillHist2D(hname, detID.ix(), detID.iy(), E);
      }
      // hname = Form("rechitMap_event%i", nev);
      // FillHist2D(hname, eta, phi, E);
    }
  }

  //======= ECAL PRESHOWER ==========
  const CaloSubdetectorGeometry* ESgeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalPreshower);

  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> handle_RecHitsES;
  iEvent.getByToken(ecalRecHitsESToken_, handle_RecHitsES);
  const EcalRecHitCollection* ecalRecHitsES = handle_RecHitsES.product();
  std::cout << "Pre-shower rechit collection size: " << ecalRecHitsES->size() << std::endl;

  if (ecalRecHitsES->size() >0){
    for(const EcalRecHit &iES : *ecalRecHitsES) {
      ESDetId detID(iES.id());
      std::cout << "DetID: " << detID << std::endl;
      auto cell = ESgeom->getGeometry(detID);
      float eta, phi, E;
      int layer = 0;
      eta = cell->getPosition().eta(); phi = cell->getPosition().phi(); E = iES.energy();
      std::cout << "ix: "<< detID.six() << " iy: "<< detID.siy() << std::endl;
      std::cout << "x: "<< cell->getPosition().x() << " y: "<< cell->getPosition().y() <<" z: " << cell->getPosition().z() << std::endl;
      std::cout << "eta: "<< cell->getPosition().eta() << " phi: "<< cell->getPosition().phi() <<" Energy: " << E << std::endl;

      if (detID.plane()==1){//shower profiles in the y (vertical) direction
        if ((eta>=0)&&(E>1.0e-10)){
          v_x_pES1.push_back((detID.six()-1)*32+ detID.strip());
          v_y_pES1.push_back((detID.siy()-1)*32+detID.strip());
          v_eta_pES1.push_back(eta);
          v_phi_pES1.push_back(phi);
          v_layer_pES1.push_back(layer);
          v_weight_pES1.push_back(E);

          hname = Form("rechitMapES1p_event%i", nev);
          FillHist2D(hname, detID.six(), (detID.siy()-1)*32+detID.strip(), E);
        }
        else if ((eta<0)&&(E>1.0e-10)){
          v_x_mES1.push_back((detID.six()-1)*32+ detID.strip());
          v_y_mES1.push_back((detID.siy()-1)*32+detID.strip());
          v_eta_mES1.push_back(eta);
          v_phi_mES1.push_back(phi);
          v_layer_mES1.push_back(layer);
          v_weight_mES1.push_back(E);

          hname = Form("rechitMapES1m_event%i", nev);
          FillHist2D(hname, detID.six(), (detID.siy()-1)*32+detID.strip(), E);
        }
      }
      else if (detID.plane()==2){//shower profiles in the x (horizontal) direction
        if ((eta>=0)&&(E>1.0e-10)){
          v_x_pES2.push_back((detID.six()-1)*32+ detID.strip());
          v_y_pES2.push_back((detID.siy()-1)*32+detID.strip());
          v_eta_pES2.push_back(eta);
          v_phi_pES2.push_back(phi);
          v_layer_pES2.push_back(layer);
          v_weight_pES2.push_back(E);

          hname = Form("rechitMapES2p_event%i", nev);
          FillHist2D(hname, (detID.six()-1)*32+ detID.strip(), detID.siy(), E);
        }
        else if ((eta<0)&&(E>1.0e-10)){
          v_x_mES2.push_back((detID.six()-1)*32+ detID.strip());
          v_y_mES2.push_back((detID.siy()-1)*32+detID.strip());
          v_eta_mES2.push_back(eta);
          v_phi_mES2.push_back(phi);
          v_layer_mES2.push_back(layer);
          v_weight_mES2.push_back(E);

          hname = Form("rechitMapES2m_event%i", nev);
          FillHist2D(hname, (detID.six()-1)*32+ detID.strip(), detID.siy(), E);
        }
      }
      // hname = Form("rechitMap_event%i", nev);
      // FillHist2D(hname, eta, phi, E);
    }
  }
  
  //======= Now apply clustering to each event =======
  const float dc_EB = 0.0174*5; //(for ECAL Barrel) chosen based on the shower size and the lateral granularity of detectors
  const float dc_EE = 10; //(for ECAL Endcap) chosen based on the shower size and the lateral granularity of detectors 
  const float dc_ES = 64; //(for ECAL preshower) chosen based on the shower size and the lateral granularity of detectors
  const float rhoc = 5; //chosen to exclude noise
  const float outlierDeltaFactor = 2; //chosen based on the shower sizes and separations
  std::cout << "\n###########################################################" <<std::endl;
  std::cout << "Setting cut-off distance 'dc' for EB and EE to: "<< dc_EB << dc_EE << std::endl;
  std::cout << "###########################################################\n" <<std::endl;
  std::vector<int> v_clusterID_EB, v_clusterID_EEp, v_clusterID_EEm, v_clusterID_ES1p, v_clusterID_ES1m, v_clusterID_ES2p, v_clusterID_ES2m;
  if (v_eta.size()>0)   {v_clusterID_EB  = mainRunEBAOD(v_eta,   v_phi,   v_layer,     v_weight,     dc_EB, rhoc, outlierDeltaFactor, 0, 10, 0);}
  if (v_x_pEE.size()>0) {v_clusterID_EEp = mainRunEEAOD(v_x_pEE, v_y_pEE, v_layer_pEE, v_weight_pEE, dc_EE, rhoc, outlierDeltaFactor, 0, 10, 0);}
  if (v_x_mEE.size()>0) {v_clusterID_EEm = mainRunEEAOD(v_x_mEE, v_y_mEE, v_layer_mEE, v_weight_mEE, dc_EE, rhoc, outlierDeltaFactor, 0, 10, 0);}
  if (v_x_pES1.size()>0) {v_clusterID_ES1p = mainRunESAOD(v_x_pES1, v_y_pES1, v_layer_pES1, v_weight_pES1, dc_ES, rhoc, outlierDeltaFactor, 0, 10, 0);}
  if (v_x_mES1.size()>0) {v_clusterID_ES1m = mainRunESAOD(v_x_mES1, v_y_mES1, v_layer_mES1, v_weight_mES1, dc_ES, rhoc, outlierDeltaFactor, 0, 10, 0);}
  if (v_x_pES2.size()>0) {v_clusterID_ES2p = mainRunESAOD(v_x_pES2, v_y_pES2, v_layer_pES2, v_weight_pES2, dc_ES, rhoc, outlierDeltaFactor, 0, 10, 0);}
  if (v_x_mES2.size()>0) {v_clusterID_ES2m = mainRunESAOD(v_x_mES2, v_y_mES2, v_layer_mES2, v_weight_mES2, dc_ES, rhoc, outlierDeltaFactor, 0, 10, 0);}

  hname = Form("clusteredMapEB_event%i", nev);
  for (int i =0; i<int(v_clusterID_EB.size()); i++){
    if (v_clusterID_EB[i] == -1){v_clusterID_EB[i] = 999;}
    FillHist2D(hname, v_eta[i], v_phi[i], v_clusterID_EB[i]+1);
  }
  hname = Form("clusteredMapEEp_event%i", nev);
  for (int i =0; i<int(v_clusterID_EEp.size()); i++){
    if (v_clusterID_EEp[i] == -1){v_clusterID_EEp[i] = 999;}
    FillHist2D(hname, v_x_pEE[i], v_y_pEE[i], v_clusterID_EEp[i]+1);
  }
   hname = Form("clusteredMapEEm_event%i", nev);
  for (int i =0; i<int(v_clusterID_EEm.size()); i++){
    if (v_clusterID_EEm[i] == -1){v_clusterID_EEm[i] = 999;}
    FillHist2D(hname, v_x_mEE[i], v_y_mEE[i], v_clusterID_EEm[i]+1);
  }
  hname = Form("clusteredMapES1p_event%i", nev);
  for (int i =0; i<int(v_clusterID_ES1p.size()); i++){
    if (v_clusterID_ES1p[i] == -1){v_clusterID_ES1p[i] = 999;}
    FillHist2D(hname, v_x_pES1[i], v_y_pES1[i], v_clusterID_ES1p[i]+1);
  }
   hname = Form("clusteredMapES1m_event%i", nev);
  for (int i =0; i<int(v_clusterID_ES1m.size()); i++){
    if (v_clusterID_ES1m[i] == -1){v_clusterID_ES1m[i] = 999;}
    FillHist2D(hname, v_x_mES1[i], v_y_mES1[i], v_clusterID_ES1m[i]+1);
  }
  hname = Form("clusteredMapES2p_event%i", nev);
  for (int i =0; i<int(v_clusterID_ES2p.size()); i++){
    if (v_clusterID_ES2p[i] == -1){v_clusterID_ES2p[i] = 999;}
    FillHist2D(hname, v_x_pES2[i], v_y_pES2[i], v_clusterID_ES2p[i]+1);
  }
   hname = Form("clusteredMapES2m_event%i", nev);
  for (int i =0; i<int(v_clusterID_ES2m.size()); i++){
    if (v_clusterID_ES2m[i] == -1){v_clusterID_ES2m[i] = 999;}
    FillHist2D(hname, v_x_mES2[i], v_y_mES2[i], v_clusterID_ES2m[i]+1);
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void ClusteringAnalyzerAOD::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void ClusteringAnalyzerAOD::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ClusteringAnalyzerAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

void ClusteringAnalyzerAOD::FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight) 
{
  std::map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    std::cout << "%FillHist2D -- Could not find histogram with name: " << histName << std::endl;
  else
    hid->second->Fill(value1, value2, weight);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ClusteringAnalyzerAOD);
