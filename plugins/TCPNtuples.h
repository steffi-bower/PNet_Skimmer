#ifndef PNet_Skimmer_MiniAODCleaner_TCPNtuples_h
#define PNet_Skimmer_MiniAODCleaner_TCPNtuples_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "PNet_Skimmer/interface/JetInfoDS.h"
#include "PNet_Skimmer/interface/MuonInfoDS.h"
#include "PNet_Skimmer/interface/ElectronInfoDS.h"
#include "PNet_Skimmer/interface/TauInfoDS.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

using namespace edm;
using namespace std;

class TCPNtuples : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit TCPNtuples(const edm::ParameterSet&);
  ~TCPNtuples() override {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void Reset() {
    event_ = -9999;
    run_ = -9999;
    lumiblock_ = -9999;
    rJetInfoData->clear();
    stdJetsInfoData->clear();
    muonInfoData->clear();
    electronInfoData->clear();
    //lowPtElectronInfoData->clear();
    tauInfoDataUnCleaned->clear();
    tauInfoDataECleaned->clear();
    //tauInfoDataLowPtECleaned->clear();
    tauInfoDataMCleaned->clear();
    tauInfoDataBoosted->clear();
    metInfo_.pt = -9999.;
    metInfo_.phi = -9999.;
    metInfo_.ptUncor = -9999.;
    metInfo_.phiUncor = -9999.;
    metInfo_.ptJECUp = -9999.;
    metInfo_.phiJECUp = -9999.;
    metInfo_.ptJERUp = -9999.;
    metInfo_.phiJERUp = -9999.;
    metInfo_.ptUncUp = -9999.;
    metInfo_.phiUncUp = -9999.;
    metInfo_.ptJECDown = -9999.;
    metInfo_.phiJECDown = -9999.;
    metInfo_.ptJERDown = -9999.;
    metInfo_.phiJERDown = -9999.;
    metInfo_.ptUncDown = -9999.;
    metInfo_.phiUncDown = -9999.;
  }

  JetInfoDS* rJetInfoData;
  JetInfoDS* stdJetsInfoData;
  MuonInfoDS* muonInfoData;
  ElectronInfoDS* electronInfoData;
  //ElectronInfoDS* lowPtElectronInfoData;
  TauInfoDS* tauInfoDataUnCleaned;
  TauInfoDS* tauInfoDataECleaned;
  //TauInfoDS* tauInfoDataLowPtECleaned;
  TauInfoDS* tauInfoDataMCleaned;
  TauInfoDS* tauInfoDataBoosted;

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

  TTree *tree;
  edm::EDGetTokenT< std::vector<pat::MET> > MET_;
  edm::EDGetTokenT< std::vector<pat::Jet> > stdJets_;
  edm::EDGetTokenT< std::vector<pat::Jet> > reclusteredJets_;
  edm::EDGetTokenT< std::vector<pat::Muon> > Muons_;
  edm::EDGetTokenT< std::vector<pat::Electron> > Electrons_;
  //edm::EDGetTokenT< std::vector<pat::Electron> > LowPtElectrons_;
  string idScoreCut_;
  edm::EDGetTokenT< std::vector<reco::Vertex> > Vertices_;
  edm::EDGetTokenT<double> rhoTag_;
  EffectiveAreas effectiveAreas_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausUnCleaned_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausECleaned_;
  //edm::EDGetTokenT< std::vector<pat::Tau> > TausLowPtECleaned_;
  edm::EDGetTokenT< std::vector<pat::Tau> > TausMCleaned_;
  edm::EDGetTokenT< reco::VertexCompositePtrCandidateCollection > secondaryVertices_;
  //edm::EDGetTokenT< edm::TriggerResults > TriggerResults_;

  //edm::EDGetTokenT< std::vector<pat::Tau> > TausBoosted_;

  int event_;
  int run_;
  int lumiblock_;

  struct MetInfo {
    MetInfo() {
      pt = phi = 0.;
      ptUncor = phiUncor = 0.;
      ptJECUp = phiJECUp = 0.;
      ptJERUp = phiJERUp = 0.;
      ptUncUp = phiUncUp = 0.;
      ptJECDown = phiJECDown = 0.;
      ptJERDown = phiJERDown = 0.;
      ptUncDown = phiUncDown = 0.;
    }
    float pt, phi;
    float ptUncor, phiUncor;
    float ptJECUp, phiJECUp;
    float ptJERUp, phiJERUp;
    float ptUncUp, phiUncUp;
    float ptJECDown, phiJECDown;
    float ptJERDown, phiJERDown;
    float ptUncDown, phiUncDown;
  };
  
  MetInfo metInfo_;

  void fillTauInfoDS(const std::vector<pat::Tau>& TauCollection, int whichColl);
  float deltaR(float phi1, float phi2, float eta1, float eta2);
  bool applyPileupJetID(const pat::Jet & jet, const std::string & level);
};

#endif
