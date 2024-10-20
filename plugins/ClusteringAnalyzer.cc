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
#include <fstream>
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
#include "ClusteringAnalyzer.h"

//root files
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TFormula.h>
#include <TVector3.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

//using reco::TrackCollection;

class ClusteringAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ClusteringAnalyzer(const edm::ParameterSet&);
  ~ClusteringAnalyzer() override;

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
  std::map<TString, TH1*> m_Histos1D;

  int nev = 0 ;

  void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight);
  // void FillHist1D(const TString& histName, const Double_t& value, const double& weight);

  // ----------member data ---------------------------
  //edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> ecalRecHitsEBToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> ecalRecHitsEEToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> ecalRecHitsESToken_;
  edm::EDGetTokenT<vector<pat::Photon>> photonToken_;
  edm::EDGetTokenT<vector<pat::PackedGenParticle>> genPartToken_;
  edm::EDGetTokenT<vector<reco::GenParticle>> genPartRecoToken_;
  // edm::EDGetTokenT<vector<pat::Electron>> electronToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  std::vector<int> EventsToScan_;
  int nEvent = 0;
  int nEventPAT = 0;

  TTree *tree2;
  std::vector<double> mHiggs;
  std::vector<double> betaval;
  std::vector<double> gammaval;
  std::vector<int> nPatPho;
  std::vector<bool> passesEvent;
  std::vector<int> nClusters;
  std::vector<double> cluster_E;
  // std::vector<double> cluster_eta;
  // std::vector<double> cluster_phi;

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
ClusteringAnalyzer::ClusteringAnalyzer(const edm::ParameterSet& iConfig):
	//tracksToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
	ecalRecHitsEBToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag > ("ecalRechitsEB"))),
	ecalRecHitsEEToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag > ("ecalRechitsEE"))),
	ecalRecHitsESToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag > ("ecalRechitsES"))),
        photonToken_(consumes<vector<pat::Photon>>(iConfig.getParameter<edm::InputTag > ("photons"))),
        genPartToken_(consumes<vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag > ("genParticles"))),
        genPartRecoToken_(consumes<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag > ("genRecoParticles"))),
        // electronToken_(consumes<vector<pat::Electron>>(iConfig.getParameter<edm::InputTag > ("electrons"))),

        caloGeometryToken_(esConsumes()),
        EventsToScan_(iConfig.getParameter<std::vector<int>>("EventsToScan"))
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
  //usesResource("TFileService");
  edm::Service<TFileService> fs2;
  tree2 = fs2->make<TTree>("Events", "Events");
  // TString hname2; // initializes variable for histogram name, to be reused for each histogram

  // hname2 = Form("mH");
  // m_Histos1D[hname2] = fs2->make<TH1F>(hname2, "Mass of Higgs", 100, 0., 100.);
  // hname2 = Form("betaval");
  // m_Histos1D[hname2] = fs2->make<TH1F>(hname2, "Beta Value", 100, 0.9, 1.);
  // hname2 = Form("gammaval");
  // m_Histos1D[hname2] = fs2->make<TH1F>(hname2, "Lorentz Factor Value", 100, 0., 250.);
  // hname2 = Form("nphoton");
  // m_Histos1D[hname2] = fs2->make<TH1I>(hname2, "Number of PAT Photons", 6, 0, 6);
  // hname2 = Form("passEvent");
  // m_Histos1D[hname2] = fs2->make<TH1I>(hname2, "PAT Photons Exist", 2, 0, 2);
  // hname2 = Form("nClsuters");
  // m_Histos1D[hname2] = fs2->make<TH1I>(hname2, "Number of CLUE Clusters", 6, 0, 6);
  // hname2 = Form("cluster_E");
  // m_Histos1D[hname2] = fs2->make<TH1F>(hname2, "Sum of Clustered Crystal Energy");


  tree2->Branch("mH", &mHiggs);
  tree2->GetBranch("mH")->SetTitle("Mass of Higgs");
  tree2->Branch("betaval", &betaval);
  tree2->GetBranch("betaval")->SetTitle("Beta Value");
  tree2->Branch("gammaval", &gammaval);
  tree2->GetBranch("gammaval")->SetTitle("Lorentz Factor Value");
  tree2->Branch("nphoton", &nPatPho);
  tree2->GetBranch("nphoton")->SetTitle("Number of PAT Photons");
  tree2->Branch("passEvent", &passesEvent);
  tree2->GetBranch("passEvent")->SetTitle("PAT Photons Exist");
  tree2->Branch("nClusters", &nClusters);
  tree2->GetBranch("nClusters")->SetTitle("Total Number of Clusters per Event");
  tree2->Branch("cluster_E", &cluster_E);
  tree2->GetBranch("cluster_E")->SetTitle("Sum of Clustered Crystal Energy");
  // tree2->Branch("cluster_eta", &cluster_eta);
  // tree2->GetBranch("cluster_eta")->SetTitle("Eta of the Cluster Seed");
  // tree2->Branch("cluster_phi", &cluster_phi);
  // tree2->GetBranch("cluster_phi")->SetTitle("Phi of the Cluster Seed");
 
  edm::Service<TFileService> fs;

  TString hname; // initializes variable for histogram name, to be reused for each histogram

  //  // ~~~~~~~~~ Include L183-L239 to make CLUE plots ~~~~~~~~~
  // for (unsigned int i=0; i < EventsToScan_.size(); i++){
  //   // ===== event maps ====== //
  //   // hname = Form("PATPhoton_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 170, -85, 85, 360, 0, 360);
  //   // m_Histos2D[hname]->SetMarkerStyle(29);
  //   // hname = Form("genPhoton_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , ETA_BINS, etabins, PHI_BINS, -1*philim, philim);
  //   //  hname = Form("rechitMap_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , ETA_BINS, etabins, PHI_BINS, -1*philim, philim);
  //   // hname = Form("clusteredMap_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , ETA_BINS, etabins, PHI_BINS, -1*philim, philim);
  //   // // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 171, -1.479, 1.479, 361, -3.15, 3.15);

  //   // ===== rechit maps for detector geometry ====== //
  //   hname = Form("rechitMapEB_event%i", EventsToScan_[i]) ;
  //   m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 170, -85, 85, 360, 0, 360);
    
  //   hname = Form("rechitMapEEp_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 101, -50, 50, 101, -50, 50);
  //    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 100, 0, 100, 100, 0, 100);
  //   hname = Form("rechitMapEEm_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 101, -50, 50, 101, -50, 50);
  //    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 100, 0, 100, 100, 0, 100);

  //   // hname = Form("rechitMapES1p_event%i", EventsToScan_[i]) ;
  //   // //m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 2520, -20, 20, 2520, -20, 20);
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 40, 0, 40, 1280, 0, 1280);
  //   // hname = Form("rechitMapES1m_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 40, 0, 40, 1280, 0, 1280);

  //   // hname = Form("rechitMapES2p_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1280, 0, 1280, 40, 0, 40);
  //   // hname = Form("rechitMapES2m_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1280, 0, 1280, 40, 0, 40);


  //   hname = Form("clusteredMapEB_event%i", EventsToScan_[i]) ;
  //   m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 170, -85, 85, 360, 0, 360);

  //   hname = Form("clusteredMapEEp_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 101, -50, 50, 101, -50, 50);
  //    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 100, 0, 100, 100, 0, 100);
  //   hname = Form("clusteredMapEEm_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 101, -50, 50, 101, -50, 50);
  //    m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 100, 0, 100, 100, 0, 100);
    
  //   // hname = Form("clusteredMapES1p_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 40, 0, 40, 1280, 0, 1280);
  //   // hname = Form("clusteredMapES1m_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 40, 0, 40, 1280, 0, 1280);

  //   // hname = Form("clusteredMapES2p_event%i", EventsToScan_[i]) ;
  //   // //m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 2520, -20, 20, 2520, -20, 20);
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1280, 0, 1280, 40, 0, 40);
  //   // hname = Form("clusteredMapES2m_event%i", EventsToScan_[i]) ;
  //   // m_Histos2D[hname] = fs->make<TH2F>(hname, hname , 1280, 0, 1280, 40, 0, 40);
  // }
  
}

ClusteringAnalyzer::~ClusteringAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void ClusteringAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  nEvent++; // increment the total event count
  mHiggs.clear();
  betaval.clear();
  gammaval.clear();
  nPatPho.clear();
  passesEvent.clear();
  nClusters.clear();
  cluster_E.clear();
  // cluster_eta.clear();
  // cluster_phi.clear();

  TString hname; // initializes variable for histogram name, to be reused for each histogram

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


  // double higgspatpt = 0;
  TLorentzVector higgspatvec;
  //======= PAT PHOTONS =======
  edm::Handle<std::vector<pat::Photon> > patpho;
  iEvent.getByToken(photonToken_, patpho);
  int nPhotons = 0;
  bool passEvent = false;
  std::cout << "\n PAT photons detected: " << patpho->size() << std::endl;


  for(size_t i = 0; i < patpho->size(); ++ i) {
    const pat::Photon & p = (*patpho)[i];
    std::cout << "Found a PAT photon with pt, eta, phi, energy : " << p.pt() << ", " << p.eta() << ", " << p.phi() << ", " << p.energy() << std::endl;
    if ((p.pt()>10)){ // & p.photonID("mvaPhoID-RunIIFall17-v2-wp90")){ // pt and ID cut on PAT photons
      nPhotons++; // increment the photon count in each step
      passEvent = true;
      TLorentzVector tempphovec;
      tempphovec.SetPtEtaPhiE(p.pt(), p.eta(), p.phi(), p.energy());
      higgspatvec += tempphovec;
    }
    if (passEvent) nEventPAT++;
    // hname = Form("PATPhoton_event%i", nev);
    // FillHist2D(hname, p.eta()*(85/1.4), p.phi()*(180/3.14159265), p.energy()) ;
  }
  // higgspatpt = higgspatvec.Pt();

//   //======= ELECTRONS =======
//   edm::Handle<std::vector<pat::Electron> > elec;
//   iEvent.getByToken(electronToken_, elec);
//   std::cout << "PAT electrons detected: " << elec->size() << std::endl;

  //======= GEN PARTICLES =======
  edm::Handle<std::vector<pat::PackedGenParticle> > genPart;
  edm::Handle<std::vector<reco::GenParticle> > genRecoPart;
  iEvent.getByToken(genPartToken_, genPart);
  iEvent.getByToken(genPartRecoToken_, genRecoPart);

  std::cout << "Gen particles collection size: " << genRecoPart->size() << std::endl;

  double mH = 0; // Higgs mass from gen higgs
  double mH_photon = 0; // Higgs mass from reconstructed gen photon
  double beta = -99; // v/c, max value is 1
  double gamma = -99; // 1/sqrt(1-beta^2), always positive

  TFormula gamma_L("lorentzFactor", "1/sqrt(1 - x^2)");

 // PdgId for photon is 22, for higgs is 25
 // Final state particles have status == 1
 // Particles produced in hard scatter have a status between 21-29
 // We use a combination of this info to get the generator level info we want
  int numgenpho = 0;
  // double higgsmass = 0;
  // TLorentzVector genhiggsvec;
  TLorentzVector four_vec;
  TVector3 boost_vec;
  TLorentzVector pho_vec;
  for(size_t i = 0; i < genRecoPart->size(); ++ i) {
    // const pat::PackedGenParticle & gen = (*genPart)[i];
    const reco::GenParticle & gen = (*genRecoPart)[i];
    if ((gen.pdgId() == 25)){ //} & (gen.status()<30) & (gen.status()>20)){
      std::cout << "Found a Higgs with status, pt, eta, phi, energy : " << gen.status() << ", " << gen.pt() << ", " << gen.eta() << ", " << gen.phi() << ", " << gen.energy() << std::endl;
      mH = gen.mass();
      four_vec.SetPtEtaPhiE(gen.pt(), gen.eta(), gen.phi(), gen.energy());
      boost_vec = four_vec.BoostVector();
      beta = boost_vec.Mag(); // Magnitude of the boost vector
      gamma = gamma_L.Eval(beta); // Lorentz factor
      std::cout << "Boost: " << beta << ", gamma_L: " << gamma << std::endl;
    }

    // if ((gen.pdgId() == 22) & (gen.status()==1) & (gen.pt()>10)){
    if ((gen.pdgId() == 22) & (gen.status()==1)){ // Final state photons
      std::cout << "Found a photon with status, pt, eta, phi, energy : " << gen.status() << ", " << gen.pt() << ", " << gen.eta() << ", " << gen.phi() << ", " << gen.energy() << std::endl;
      numgenpho += 1;
      TLorentzVector tempgenvec;
      tempgenvec.SetPtEtaPhiE(gen.pt(), gen.eta(), gen.phi(), gen.energy());
      pho_vec += tempgenvec;
      // hname = Form("genPhoton_event%i", nev);
      // FillHist2D(hname, gen.eta(), gen.phi(), gen.energy()) ;
    }
  }
  mH_photon = pho_vec.M();
  if ((mH - mH_photon)>0.0001){
    std::cout << "!!!! Disparity in gen part mass calc !!!!" << std::endl;
    std::cout << mH << "!=" << mH_photon << std::endl;
  }

  // ~~~~~~~ nTupling for GEN/PAT info ~~~~~~~
  mHiggs.push_back(mH);
  betaval.push_back(beta);
  gammaval.push_back(gamma);
  nPatPho.push_back(nPhotons);
  passesEvent.push_back(passEvent);

  //======= ECAL BARREL ==========
  edm::ESHandle<CaloGeometry> caloGeometry = iSetup.getHandle(caloGeometryToken_);
  const CaloSubdetectorGeometry* EBgeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);

  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> handle_RecHitsEB;
  iEvent.getByToken(ecalRecHitsEBToken_, handle_RecHitsEB);
  const EcalRecHitCollection* ecalRecHitsEB = handle_RecHitsEB.product();
  std::cout << "rechit collection size: " << ecalRecHitsEB->size() << std::endl;

  double EBtotE = 0;
  // TLorentzVector EBtotpT;
  if (ecalRecHitsEB->size() >0){
    for(const EcalRecHit &iEB : *ecalRecHitsEB) {
      EBDetId detID(iEB.id());
      // std::cout << "DetID: " << detID << std::endl;
      auto cell = EBgeom->getGeometry(detID);
      // float eta, phi, E;
      float E;
      // float pT;
      // pT = cell->getPosition().rho();
      int layer = 0;
      // eta = cell->getPosition().eta(); phi = cell->getPosition().phi(); E = iEB.energy();
      E = iEB.energy();
      // std::cout << "eta: "<< eta << " phi: "<< phi <<" Energy: " << E << std::endl;
      // std::cout << "ieta: "<< detID.ieta() << " iphi: "<< detID.iphi() <<" Energy: " << E << std::endl;
      // std::cout << "pT:" << detID.ipt() << std::endl;

      EBtotE += E;
      // TLorentzVector tempEB;
      // // auto temppT = 0;
      // // temppT += E*cosh(detID.ieta());
      // tempEB.SetPtEtaPhiE(detID.irho(), detID.ieta(), detID.iphi(), E);
      // EBtotpT += tempEB;
      v_eta.push_back(detID.ieta());
      v_phi.push_back(detID.iphi());
      v_layer.push_back(layer);
      v_weight.push_back(E);

      // hname = Form("rechitMapEB_event%i", nev);
      // FillHist2D(hname, detID.ieta(), detID.iphi(), E);

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

  double EEtotE = 0;
  double EEtotpT = 0;
  if (ecalRecHitsEE->size() >0){
    for(const EcalRecHit &iEE : *ecalRecHitsEE) {
      EEDetId detID(iEE.id());
      // std::cout << "DetID: " << detID << std::endl;
      auto cell = EEgeom->getGeometry(detID);
      float eta, phi, E;
      int layer = 0;
      eta = cell->getPosition().eta(); phi = cell->getPosition().phi(); E = iEE.energy();
      // std::cout << "ix: "<< detID.ix() << " iy: "<< detID.iy() << std::endl;
      // std::cout << "x: "<< cell->getPosition().x() << " y: "<< cell->getPosition().y() <<" z: " << cell->getPosition().z() << std::endl;
      // std::cout << "eta: "<< cell->getPosition().eta() << " phi: "<< cell->getPosition().phi() <<" Energy: " << E << std::endl;
      EEtotE += E;
      EEtotpT += E*cosh(eta);
      if ((eta>=0)&&(E>1.0e-10)){
        // EEtotE += E;
        v_x_pEE.push_back(detID.ix());
        v_y_pEE.push_back(detID.iy());
        v_eta_pEE.push_back(eta);
        v_phi_pEE.push_back(phi);
        v_layer_pEE.push_back(layer);
        v_weight_pEE.push_back(E);

        // hname = Form("rechitMapEEp_event%i", nev);
        // FillHist2D(hname, detID.ix(), detID.iy(), E);
      }
      else if ((eta<0)&&(E>1.0e-10)){
        // EEtotE += E;
        v_x_mEE.push_back(detID.ix());
        v_y_mEE.push_back(detID.iy());
        v_eta_mEE.push_back(eta);
        v_phi_mEE.push_back(phi);
        v_layer_mEE.push_back(layer);
        v_weight_mEE.push_back(E);

        // hname = Form("rechitMapEEm_event%i", nev);
        // FillHist2D(hname, detID.ix(), detID.iy(), E);
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

  // double EStotE = 0;
  // double EStotpT = 0;
  if (ecalRecHitsES->size() >0){
    for(const EcalRecHit &iES : *ecalRecHitsES) {
      ESDetId detID(iES.id());
      // std::cout << "DetID: " << detID << std::endl;
      auto cell = ESgeom->getGeometry(detID);
      float eta, phi, E;
      int layer = 0;
      eta = cell->getPosition().eta(); phi = cell->getPosition().phi(); E = iES.energy();
      // std::cout << "ix: "<< detID.six() << " iy: "<< detID.siy() << std::endl;
      // std::cout << "x: "<< cell->getPosition().x() << " y: "<< cell->getPosition().y() <<" z: " << cell->getPosition().z() << std::endl;
      // std::cout << "eta: "<< cell->getPosition().eta() << " phi: "<< cell->getPosition().phi() <<" Energy: " << E << std::endl;

      // EStotE += E;
      // EStotpT += E*cosh(eta);
      if (detID.plane()==1){//shower profiles in the y (vertical) direction
        if ((eta>=0)&&(E>1.0e-10)){
          v_x_pES1.push_back((detID.six()-1)*32+ detID.strip());
          v_y_pES1.push_back((detID.siy()-1)*32+detID.strip());
          v_eta_pES1.push_back(eta);
          v_phi_pES1.push_back(phi);
          v_layer_pES1.push_back(layer);
          v_weight_pES1.push_back(E);

          // hname = Form("rechitMapES1p_event%i", nev);
          // FillHist2D(hname, detID.six(), (detID.siy()-1)*32+detID.strip(), E);
        }
        else if ((eta<0)&&(E>1.0e-10)){
          v_x_mES1.push_back((detID.six()-1)*32+ detID.strip());
          v_y_mES1.push_back((detID.siy()-1)*32+detID.strip());
          v_eta_mES1.push_back(eta);
          v_phi_mES1.push_back(phi);
          v_layer_mES1.push_back(layer);
          v_weight_mES1.push_back(E);

          // hname = Form("rechitMapES1m_event%i", nev);
          // FillHist2D(hname, detID.six(), (detID.siy()-1)*32+detID.strip(), E);
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

          // hname = Form("rechitMapES2p_event%i", nev);
          // FillHist2D(hname, (detID.six()-1)*32+ detID.strip(), detID.siy(), E);
        }
        else if ((eta<0)&&(E>1.0e-10)){
          v_x_mES2.push_back((detID.six()-1)*32+ detID.strip());
          v_y_mES2.push_back((detID.siy()-1)*32+detID.strip());
          v_eta_mES2.push_back(eta);
          v_phi_mES2.push_back(phi);
          v_layer_mES2.push_back(layer);
          v_weight_mES2.push_back(E);

          // hname = Form("rechitMapES2m_event%i", nev);
          // FillHist2D(hname, (detID.six()-1)*32+ detID.strip(), detID.siy(), E);
        }
      }
      // hname = Form("rechitMap_event%i", nev);
      // FillHist2D(hname, eta, phi, E);
    }
  }
  
  //======= Now apply clustering to each event =======
  // const float dc_EB = 0.0174*5; //(for ECAL Barrel) chosen based on the shower size and the lateral granularity of detectors
  const float dc_EB = 2; //(for ECAL Barrel) chosen based on the shower size and the lateral granularity of detectors
  const float dc_EE = 1; //(for ECAL Endcap) chosen based on the shower size and the lateral granularity of detectors 
  const float dc_ES = 2; //(for ECAL preshower) chosen based on the shower size and the lateral granularity of detectors
  const float rhoc = 15; //chosen to exclude noise (Originally 5)
  const float outlierDeltaFactor = 3; //chosen based on the shower sizes and separations
  const int Nevents = EventsToScan_.size();
  std::cout << "\n###########################################################" <<std::endl;
  std::cout << "Setting cut-off distance 'dc' for EB and EE to: "<< dc_EB << dc_EE << std::endl;
  std::cout << "###########################################################\n" <<std::endl;
  std::vector<int> v_clusterID_EB, v_clusterID_EEp, v_clusterID_EEm, v_clusterID_ES1p, v_clusterID_ES1m, v_clusterID_ES2p, v_clusterID_ES2m;
  if (v_eta.size()>0)   {v_clusterID_EB  = mainRunEB(v_eta,   v_phi,   v_layer,     v_weight,     dc_EB, rhoc, outlierDeltaFactor, 0, Nevents, 0);} 
  if (v_x_pEE.size()>0) {v_clusterID_EEp = mainRunEE(v_x_pEE, v_y_pEE, v_layer_pEE, v_weight_pEE, dc_EE, rhoc, outlierDeltaFactor, 0, Nevents, 0);}
  if (v_x_mEE.size()>0) {v_clusterID_EEm = mainRunEE(v_x_mEE, v_y_mEE, v_layer_mEE, v_weight_mEE, dc_EE, rhoc, outlierDeltaFactor, 0, Nevents, 0);}
  if (v_x_pES1.size()>0) {v_clusterID_ES1p = mainRunES(v_x_pES1, v_y_pES1, v_layer_pES1, v_weight_pES1, dc_ES, rhoc, outlierDeltaFactor, 0, Nevents, 0);}
  if (v_x_mES1.size()>0) {v_clusterID_ES1m = mainRunES(v_x_mES1, v_y_mES1, v_layer_mES1, v_weight_mES1, dc_ES, rhoc, outlierDeltaFactor, 0, Nevents, 0);}
  if (v_x_pES2.size()>0) {v_clusterID_ES2p = mainRunES(v_x_pES2, v_y_pES2, v_layer_pES2, v_weight_pES2, dc_ES, rhoc, outlierDeltaFactor, 0, Nevents, 0);}
  if (v_x_mES2.size()>0) {v_clusterID_ES2m = mainRunES(v_x_mES2, v_y_mES2, v_layer_mES2, v_weight_mES2, dc_ES, rhoc, outlierDeltaFactor, 0, Nevents, 0);}

  int totClusters = 0;
  auto maxElementEB = std::max_element(v_clusterID_EB.begin(), v_clusterID_EB.end());
  if (v_eta.size()>0) {
      std::cout << "Number of EB clusters found: " << *maxElementEB+1 << std::endl;
      totClusters += *maxElementEB+1;
    }
  auto maxElementEEp = std::max_element(v_clusterID_EEp.begin(), v_clusterID_EEp.end());
  if (v_x_pEE.size()>0) {
      std::cout << "Number of EEp clusters found: " << *maxElementEEp+1 << std::endl;
      totClusters += *maxElementEEp+1;
    }
  auto maxElementEEm = std::max_element(v_clusterID_EEm.begin(), v_clusterID_EEm.end());
  if (v_x_mEE.size()>0) {
      std::cout << "Number of EEm clusters found: " << *maxElementEEm+1 << std::endl;
      totClusters += *maxElementEEm+1;
    }

  // std::cout << "Events" << EventsToScan_ << std::endl;

  double higgsE = EBtotE + EEtotE;
  // ~~~~~~~ nTupling for GEN/PAT info ~~~~~~~
  nClusters.push_back(totClusters);
  cluster_E.push_back(higgsE);
  // cluster_eta.push_back();
  // cluster_phi.push_back();

  // // ~~~~~~~~~~~~~~ Making CSV File ~~~~~~~~~~~~~~~~~
  // double higgsE = EBtotE + EEtotE;
  // std::string filename = "clustering10.0Ma100.3000vpTME1.csv";
  // if (numgenpho==2){
  //   std::ofstream outFile(filename, std::ios::app);
  //   // outFile.open();
  //   // if (!outFile.is_open()) {
  //   //   std::cerr << "Error opening file: " << filename << std::endl;
  //   //   return 1;
  //   // }
  //   outFile << totClusters << "," << higgspt << "," << mH << "," << higgsE << std::endl;
  //   outFile.close();
  // }


  // std::cout << "" << std::endl;
  // hname = Form("clusteredMapEB_event%i", nev);
  // for (int i =0; i<int(v_clusterID_EB.size()); i++){
  //   if (v_clusterID_EB[i] == -1){v_clusterID_EB[i] = 999;}
  //   FillHist2D(hname, v_eta[i], v_phi[i], v_clusterID_EB[i]+1);
  // }
  // hname = Form("clusteredMapEEp_event%i", nev);
  // for (int i =0; i<int(v_clusterID_EEp.size()); i++){
  //   if (v_clusterID_EEp[i] == -1){v_clusterID_EEp[i] = 999;}
  //   FillHist2D(hname, v_x_pEE[i], v_y_pEE[i], v_clusterID_EEp[i]+1);
  // }
  //  hname = Form("clusteredMapEEm_event%i", nev);
  // for (int i =0; i<int(v_clusterID_EEm.size()); i++){
  //   if (v_clusterID_EEm[i] == -1){v_clusterID_EEm[i] = 999;}
  //   FillHist2D(hname, v_x_mEE[i], v_y_mEE[i], v_clusterID_EEm[i]+1);
  // }
  // // hname = Form("clusteredMapES1p_event%i", nev);
  // // for (int i =0; i<int(v_clusterID_ES1p.size()); i++){
  // //   if (v_clusterID_ES1p[i] == -1){v_clusterID_ES1p[i] = 999;}
  // //   FillHist2D(hname, v_x_pES1[i], v_y_pES1[i], v_clusterID_ES1p[i]+1);
  // // }
  // //  hname = Form("clusteredMapES1m_event%i", nev);
  // // for (int i =0; i<int(v_clusterID_ES1m.size()); i++){
  // //   if (v_clusterID_ES1m[i] == -1){v_clusterID_ES1m[i] = 999;}
  // //   FillHist2D(hname, v_x_mES1[i], v_y_mES1[i], v_clusterID_ES1m[i]+1);
  // // }
  // // hname = Form("clusteredMapES2p_event%i", nev);
  // // for (int i =0; i<int(v_clusterID_ES2p.size()); i++){
  // //   if (v_clusterID_ES2p[i] == -1){v_clusterID_ES2p[i] = 999;}
  // //   FillHist2D(hname, v_x_pES2[i], v_y_pES2[i], v_clusterID_ES2p[i]+1);
  // // }
  // //  hname = Form("clusteredMapES2m_event%i", nev);
  // // for (int i =0; i<int(v_clusterID_ES2m.size()); i++){
  // //   if (v_clusterID_ES2m[i] == -1){v_clusterID_ES2m[i] = 999;}
  // //   FillHist2D(hname, v_x_mES2[i], v_y_mES2[i], v_clusterID_ES2m[i]+1);
  // // }

  tree2->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void ClusteringAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void ClusteringAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ClusteringAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

void ClusteringAnalyzer::FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight) 
{
  std::map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    std::cout << "%FillHist2D -- Could not find histogram with name: " << histName << std::endl;
  else
    hid->second->Fill(value1, value2, weight);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ClusteringAnalyzer);
