// system include files
#include <memory>
#include <iostream>
#include <regex>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TTree.h"

using namespace edm;
using namespace std;

class TCPTrigNtuples : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {
public:

  explicit TCPTrigNtuples(const edm::ParameterSet&);
  ~TCPTrigNtuples() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override;
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;


  edm::EDGetTokenT< edm::TriggerResults > TriggerResults_;
  
  TTree *tree;
  int event_;
  bool isIsoMu_;
  bool isIsoEle_;
  bool isDoubleMu_;
  bool isDoubleMuMass8_;
  bool isDoubleMuSS_;
  bool isDoubleMuMass3p8_;
  bool isDoubleEG_;
  bool isMuonEG_;
  bool isBPH_;
};

TCPTrigNtuples::TCPTrigNtuples(const edm::ParameterSet& iConfig) :
  TriggerResults_(consumes< edm::TriggerResults >(iConfig.getParameter<edm::InputTag>("TriggerResults"))) {
  usesResource(TFileService::kSharedResource);
  //std::cout << "debug0" << "\n";
}

void TCPTrigNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TCPTrigNtuples::beginJob() {
  
  edm::Service<TFileService> fs;
  fs->mkdir( "trigger" );
  
  tree = fs->make<TTree>("triggerTree", "");
  
  tree->Branch("event", &event_, "event/I");
  
  tree->Branch("isIsoMu", &isIsoMu_, "isIsoMu/O");
  tree->Branch("isIsoEle", &isIsoEle_, "isIsoEle/O");
  tree->Branch("isDoubleMu", &isDoubleMu_, "isDoubleMu/O");
  tree->Branch("isDoubleMuSS", &isDoubleMuSS_, "isDoubleMuSS/O");
  tree->Branch("isDoubleMuMass8", &isDoubleMuMass8_, "isDoubleMuMass8/O");
  tree->Branch("isDoubleMuMass3p8", &isDoubleMuMass3p8_, "isDoubleMuMass3p8/O");
  tree->Branch("isDoubleEG", &isDoubleEG_, "isDoubleEG/O");
  tree->Branch("isMuonEG", &isMuonEG_, "isMuonEG/O");
  tree->Branch("isBPH", &isBPH_, "isBPH/O");
}

void TCPTrigNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  int Event = iEvent.id().event();
  event_ = Event;

  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByToken(TriggerResults_, triggerResultsHandle);
  TriggerResults triggerResults = *triggerResultsHandle;
  auto & names = iEvent.triggerNames(*triggerResultsHandle);

  //std::cout << "debug3" << "\n";


  isIsoMu_ = 0;
  isIsoEle_ = 0;
  isDoubleMu_ = 0;
  isDoubleMuMass8_ = 0;
  isDoubleMuMass3p8_=0;
  isDoubleMuSS_=0;
  isDoubleEG_ = 0;
  isMuonEG_ = 0;
  isBPH_=0;
  //std::cout << "debug4" << "\n";

  for (unsigned i = 0, n = triggerResults.size(); i < n; ++i){
    std::string Trigger;
    Trigger = names.triggerName(i);
    std::regex HLT_Mu9_IP6("(HLT_Mu9_IP6_part)(.*)");
    if ( std::regex_match(Trigger, HLT_Mu9_IP6) && (triggerResults.accept(i) == 1) ) isBPH_ = 1;


  //if ( triggerResults.accept(names.triggerIndex("HLT_IsoMu24_v13")) ) isIsoMu_ = 1;
  //if ( triggerResults.accept(names.triggerIndex("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v19")) ) isDoubleEG_ = 1;
  //if ( triggerResults.accept(names.triggerIndex("HLT_Mu18_Mu9_SameSign_DZ_v4")) ) isDoubleMuSS_ = 1;

  //if ( triggerResults.accept(names.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v5")) ) isDoubleMuMass8_ = 1;
  //if ( triggerResults.accept(names.triggerIndex("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v5")) ) isDoubleMuMass3p8_ = 1;
  //if ( triggerResults.accept(names.triggerIndex("HLT_Mu37_TkMu27_v5")) ) isDoubleMu_ = 1;

  //if ( triggerResults.accept(names.triggerIndex("HLT_Ele32_WPTight_Gsf_v15"))) isIsoEle_ = 1;

  //if (triggerResults.accept(names.triggerIndex("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15")) ) isMuonEG_ = 1;

  }
  tree->Fill();
  //std::cout<<"PostTriggerbuild\n";
}

DEFINE_FWK_MODULE(TCPTrigNtuples);
