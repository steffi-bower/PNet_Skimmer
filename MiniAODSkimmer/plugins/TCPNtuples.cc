#include "BoostedDiTau/MiniAODSkimmer/plugins/TCPNtuples.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

TCPNtuples::TCPNtuples(const edm::ParameterSet& iConfig) :
  MET_(consumes< vector<pat::MET> > (iConfig.getParameter<edm::InputTag>("METCollection"))),
  stdJets_(consumes< vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("STDJetCollection"))),
  reclusteredJets_(consumes< vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("ReclusteredJetCollection"))),
  Muons_(consumes< vector<pat::Muon> > (iConfig.getParameter<edm::InputTag>("MuonCollection"))),
  Electrons_(consumes< vector<pat::Electron> > (iConfig.getParameter<edm::InputTag>("ElectronCollection"))),
  //LowPtElectrons_(consumes< vector<pat::Electron> > (iConfig.getParameter<edm::InputTag>("LowPtElectronCollection"))),
  //idScoreCut_(iConfig.getParameter<string>("LowPtEIdScoreCut")),
  Vertices_(consumes< vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("VertexCollection"))),
  rhoTag_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoTag"))),
  effectiveAreas_((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
  TausUnCleaned_(consumes< vector<pat::Tau> > (iConfig.getParameter<edm::InputTag>("UnCleanedTauCollection"))),
  TausECleaned_(consumes< vector<pat::Tau> > (iConfig.getParameter<edm::InputTag>("ECleanedTauCollection"))),
  //TausLowPtECleaned_(consumes< vector<pat::Tau> > (iConfig.getParameter<edm::InputTag>("LowPtECleanedTauCollection"))),
  TausMCleaned_(consumes< vector<pat::Tau> > (iConfig.getParameter<edm::InputTag>("MCleanedTauCollection"))),  
  secondaryVertices_(consumes<reco::VertexCompositePtrCandidateCollection>  (iConfig.getParameter<edm::InputTag>("sVertices"))){
  //TriggerResults_(consumes< edm::TriggerResults >(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  //TausBoosted_(consumes< vector<pat::Tau> > (iConfig.getParameter<edm::InputTag>("BoostedTauCollection"))) {
  usesResource(TFileService::kSharedResource);

  }

void TCPNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TCPNtuples::beginJob() {
  
  edm::Service<TFileService> fs;
  fs->mkdir( "analysis" );
  
  tree = fs->make<TTree>("analysisTree", "");
  
  tree->Branch("run", &run_, "run/I");
  tree->Branch("lumiblock", &lumiblock_, "lumiblock/I");
  tree->Branch("event", &event_, "event/I");
  
  rJetInfoData = new JetInfoDS();
  stdJetsInfoData = new JetInfoDS();
  muonInfoData = new MuonInfoDS();
  electronInfoData = new ElectronInfoDS();
  //lowPtElectronInfoData = new ElectronInfoDS();
  tauInfoDataUnCleaned = new TauInfoDS();
  tauInfoDataECleaned = new TauInfoDS();
  //tauInfoDataLowPtECleaned = new TauInfoDS();
  tauInfoDataMCleaned = new TauInfoDS();
  tauInfoDataBoosted = new TauInfoDS();
  
  tree->Branch("ReclusteredJets", "JetInfoDS", &rJetInfoData);
  tree->Branch("StdJets", "JetInfoDS", &stdJetsInfoData);
  tree->Branch("Muons", "MuonInfoDS", &muonInfoData);
  tree->Branch("Electrons", "ElectronInfoDS", &electronInfoData);
  //tree->Branch("LowPtElectrons", "ElectronInfoDS", &lowPtElectronInfoData);
  tree->Branch("TausUnCleaned", "TauInfoDS", &tauInfoDataUnCleaned);
  tree->Branch("TausECleaned", "TauInfoDS", &tauInfoDataECleaned);
  //tree->Branch("TausLowPtECleaned", "TauInfoDS", &tauInfoDataLowPtECleaned);
  tree->Branch("TausMCleaned", "TauInfoDS", &tauInfoDataMCleaned);
  //tree->Branch("TausBoosted", "TauInfoDS", &tauInfoDataBoosted);
  tree->Branch("Mets", &metInfo_, "pt/F:phi/F:ptUncor/F:phiUncor/F:ptJECUp/F:phiJECUp/F:ptJERUp/F:phiJERUp/F:ptUncUp/F:phiUncUp/F:ptJECDown/F:phiJECDown/F:ptJERDown/F:phiJERDown/F:ptUncDown/F:phiUncDown/F");
}

void TCPNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  Reset();

  int Event = iEvent.id().event();
  int Run = iEvent.id().run();
  int LumiBlock = iEvent.id().luminosityBlock();

  run_ = Run;
  lumiblock_ = LumiBlock;
  event_ = Event;

  //edm::Handle<edm::TriggerResults> triggerResultsHandle;
  //iEvent.getByToken(TriggerResults_, triggerResultsHandle);
  //TriggerResults triggerResults = *triggerResultsHandle;
  //auto & names = iEvent.triggerNames(*triggerResultsHandle);
  //if (
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part0_v3"))||
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part1_v3"))||
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part2_v3"))||
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part3_v3"))||
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part4_v3"))
  //){
    edm::Handle< std::vector<pat::Jet> > rJetsHandle;
    iEvent.getByToken(reclusteredJets_, rJetsHandle);
    auto rJets = *rJetsHandle;
    edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVerticesH;
    iEvent.getByToken(secondaryVertices_, secondaryVerticesH);
    reco::VertexCompositePtrCandidateCollection svColl = *secondaryVerticesH;
    edm::Handle< std::vector<reco::Vertex> > VerticesHandle;
    iEvent.getByToken(Vertices_, VerticesHandle);
    auto Vertices = *VerticesHandle;
    auto PrimaryVertex = Vertices[0];
    if (rJets.size() > 0) {
      for (unsigned int i = 0; i < rJets.size(); ++i) {
        auto jet = rJets[i];
        if (jet.pt() < 10 || abs(jet.eta()) > 2.4) continue;
        float NHF  = jet.neutralHadronEnergyFraction();
        float NEMF = jet.neutralEmEnergyFraction();
        float CHF  = jet.chargedHadronEnergyFraction();
        float MUF  = jet.muonEnergyFraction();
        float CEMF = jet.chargedEmEnergyFraction();
        auto NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
        //auto NumNeutralParticles = jet.neutralMultiplicity();
        auto CHM = jet.chargedMultiplicity();
        
        bool jetID = CHM>0 && 
        CHF>0 && 
        NumConst>1 && 
        NEMF<0.9 && 
        NHF < 0.9;

        bool jetIDLepVeto = CEMF<0.8 && 
        CHM>0 && 
        CHF>0 && 
        NumConst>1 && 
        NEMF<0.9 && 
        MUF <0.8 && 
        NHF < 0.9;
        if (jetID) {
          JetInfo j;
          //if (not jet.genJet()) j.genJet = 0;
          //else j.genJet = 1;

          j.pt = jet.pt();
          j.eta = jet.eta();
          j.phi = jet.phi();
          j.mass = jet.mass();
          j.ptuncor = jet.correctedP4(0).Pt();
          j.flavorprobb = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
          j.flavorprobbb = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
          j.flavorproblepb = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
          j.pnet_score = 0.0;
          j.ditau2017v1 = 

          j.deepcsv = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
          if (jetIDLepVeto) j.id = 2;
          else j.id = 1;
          int jetpuid = 0;
          if(applyPileupJetID(jet,"loose"))  jetpuid += 1;
          if(applyPileupJetID(jet,"medium")) jetpuid += 2;
          if(applyPileupJetID(jet,"tight"))  jetpuid += 4;
          j.puid=jetpuid;
          if (jetpuid<1) continue;

          rJetInfoData->push_back(j);
        }
      }
    }
    edm::Handle< std::vector<pat::Jet> > stdJetsHandle;
    iEvent.getByToken(stdJets_, stdJetsHandle);
    auto stdJets = *stdJetsHandle;
    
    if (stdJets.size() > 0) {
      for (unsigned int i = 0; i < stdJets.size(); ++i) {
        auto jet = stdJets[i];
        if (jet.pt() < 10 || abs(jet.eta()) > 2.4) continue;
        float NHF  = jet.neutralHadronEnergyFraction();
        float NEMF = jet.neutralEmEnergyFraction();
        float CHF  = jet.chargedHadronEnergyFraction();
        float MUF  = jet.muonEnergyFraction();
        float CEMF = jet.chargedEmEnergyFraction();
        auto NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
        //auto NumNeutralParticles = jet.neutralMultiplicity();
        auto CHM = jet.chargedMultiplicity();
        
        bool jetID = CHM>0 && 
        CHF>0 && 
        NumConst>1 && 
        NEMF<0.9 && 
        NHF < 0.9;

        bool jetIDLepVeto = CEMF<0.8 && 
        CHM>0 && 
        CHF>0 && 
        NumConst>1 && 
        NEMF<0.9 && 
        MUF <0.8 && 
        NHF < 0.9;
        if (jetID) {

          JetInfo j;
          j.pt = jet.pt();
          j.eta = jet.eta();
          j.phi = jet.phi();
          j.mass = jet.mass();
          //if (not jet.genJet()) j.genJet = 0;
          //else j.genJet = 1;


          j.ptuncor = jet.correctedP4(0).Pt();
          j.flavorprobb = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
          j.flavorprobbb = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
          j.flavorproblepb = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
          j.pnet_score = jet.userFloat("PNetScore");
          j.ditau2017v1 = jet.userFloat("ditau2017v1");
          j.ditau2017MDv1 = jet.userFloat("ditau2017MDv1");
          j.ditau2017v2 = jet.userFloat("ditau2017v2");
          j.ditau2017MDv2 = jet.userFloat("ditau2017MDv2");
          j.deepcsv = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
          if (jetIDLepVeto) j.id = 2;
          else j.id = 1;
          int jetpuid = 0;
          if(applyPileupJetID(jet,"loose"))  jetpuid += 1;
          if(applyPileupJetID(jet,"medium")) jetpuid += 2;
          if(applyPileupJetID(jet,"tight"))  jetpuid += 4;
          j.puid=jetpuid;
          if (jetpuid<1) continue;
          stdJetsInfoData->push_back(j);
        }
      }
    }
    edm::Handle< std::vector<pat::Muon> > MuonsHandle;
    iEvent.getByToken(Muons_, MuonsHandle);
    auto Muons = *MuonsHandle;



    if (Muons.size() > 0) {
      for (unsigned int i = 0; i < Muons.size(); ++i) {
        auto muon = Muons[i];
        if (muon.pt() < 3 || abs(muon.eta()) > 2.4) continue;
        MuonInfo m;
        m.pt = muon.pt();
        m.eta = muon.eta();
        m.phi = muon.phi();
        m.mass = muon.mass();
        m.charge = muon.charge();
        m.dB = muon.dB();
        m.edB = muon.edB();
        m.isBPHMuon = muon.triggered("HLT_Mu9_IP6_*");
        m.isMuIsoTrigMuon = muon.triggered("HLT_IsoMu24_v*");
        m.isDoubleMu = muon.triggered("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*");
        m.iso = (muon.pfIsolationR04().sumChargedHadronPt+max(0.,muon.pfIsolationR04().sumPhotonEt+muon.pfIsolationR04().sumNeutralHadronEt-0.5*muon.pfIsolationR04().sumPUPt))/muon.pt();
        if (muon.isTightMuon(PrimaryVertex)) m.id = 3;
        else if (muon.isMediumMuon()) m.id = 2;
        else if (muon.isLooseMuon()) m.id = 1;
        else m.id=-1;
        m.dxy = muon.muonBestTrack()->dxy();
        m.dz = muon.muonBestTrack()->dz();
        muonInfoData->push_back(m);
      }
    }

    edm::Handle< std::vector<pat::Electron> > ElectronsHandle;
    iEvent.getByToken(Electrons_, ElectronsHandle);
    auto Electrons = *ElectronsHandle;

    edm::Handle<double> pRho;
    iEvent.getByToken(rhoTag_, pRho);

    if (Electrons.size() > 0) {
      for (unsigned int i = 0; i < Electrons.size(); ++i) {
        auto electron = Electrons[i];
        
        bool isLoose = false;
        bool isMedium = false;
        bool isTight = false;
        bool isLooseRelIso = false;
        bool isMediumRelIso = false;
        bool isTightRelIso = false;

        // ====== implement impact parameter cuts =======
        // reference: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria_for_V

        // ---  full5x5_sigmaIetaIeta ---
        double sigmaIetaIeta = electron.full5x5_sigmaIetaIeta();

        // --- fabs(dEtaSeed) ---
        double dEtaSeed = fabs(electron.superCluster().isNonnull() && electron.superCluster()->seed().isNonnull() ? electron.deltaEtaSuperClusterTrackAtVtx() - electron.superCluster()->eta() + electron.superCluster()->seed()->eta() : std::numeric_limits<float>::max()); 
        
        // --- fabs(dPhiIn) ---
        double dPhiIn = fabs(electron.deltaPhiSuperClusterTrackAtVtx());
        
        // --- variables for H/E cuts ---
        double HoE = electron.hadronicOverEm();
        double rho = pRho.isValid() ? (*pRho) : 0; 
        double energy = electron.superCluster()->energy();

        // --- variables for relIsoWithEffectiveArea ---
        double chad = electron.pfIsolationVariables().sumChargedHadronPt;
        double nhad = electron.pfIsolationVariables().sumNeutralHadronEt;
        double pho = electron.pfIsolationVariables().sumPhotonEt;
        double elePt = electron.pt();
        double eleEta = electron.superCluster()->eta();
        double eArea = effectiveAreas_.getEffectiveArea(fabs(eleEta));
        double relIsoWithEffectiveArea = (chad + std::max(0.0, nhad + pho - rho*eArea)) / elePt;

        // --- variables for fabs(1/E-1/p) ---
        double eInverseMinusPInverse = fabs(1.0 - electron.eSuperClusterOverP())*(1.0/electron.ecalEnergy());

        // --- expected missing inner hits ---
        int mHits = electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);

        // --- pass conversion veto ---
        bool isPassConVeto = electron.passConversionVeto();

        // ========= select electrons in different cut-based ID accordingly ==========
        if (fabs(eleEta) <= 1.479) {
    isLoose = (sigmaIetaIeta < 0.0112) &&
      (dEtaSeed < 0.00377) &&
      (dPhiIn < 0.0884) &&
      (HoE < 0.05 + 1.16/energy + 0.0324*rho/energy) &&
      (eInverseMinusPInverse < 0.193) &&
      (mHits <= 1) &&
      (isPassConVeto == true);
    isLooseRelIso = (relIsoWithEffectiveArea < 0.112 + 0.506/elePt);
    
    isMedium = (sigmaIetaIeta < 0.0106) &&
      (dEtaSeed < 0.0032) &&
      (dPhiIn < 0.0547) &&
      (HoE < 0.046 + 1.16/energy + 0.0324*rho/energy) &&
      (eInverseMinusPInverse < 0.184) &&
      (mHits <= 1) &&
      (isPassConVeto == true);
    isMediumRelIso = (relIsoWithEffectiveArea < 0.0478 + 0.506/elePt);

    isTight = (sigmaIetaIeta < 0.0104) &&
      (dEtaSeed < 0.00255) &&
      (dPhiIn < 0.022) &&
      (HoE < 0.026 + 1.15/energy + 0.0324*rho/energy) &&
      (eInverseMinusPInverse < 0.159) &&
      (mHits <= 1) &&
      (isPassConVeto == true);
    isTightRelIso = (relIsoWithEffectiveArea < 0.0287 + 0.506/elePt);
        }// endif (fabs(eleEta) <= 1.479)

        else {
    isLoose = (sigmaIetaIeta < 0.0425) &&
      (dEtaSeed < 0.00674) &&
      (dPhiIn < 0.169) &&
      (HoE < 0.0441 + 2.54/energy + 0.183*rho/energy) &&
      (eInverseMinusPInverse < 0.111) &&
      (mHits <= 1) &&
      (isPassConVeto == true);
    isLooseRelIso = (relIsoWithEffectiveArea < 0.108 + 0.963/elePt);

    isMedium = (sigmaIetaIeta < 0.0387) &&
      (dEtaSeed < 0.00632) &&
      (dPhiIn < 0.0394) &&
      (HoE < 0.0275 + 2.52/energy + 0.183*rho/energy) &&
      (eInverseMinusPInverse < 0.0721) &&
      (mHits <= 1) &&
      (isPassConVeto == true);
    isMediumRelIso = (relIsoWithEffectiveArea < 0.0658 + 0.963/elePt);
    
    isTight = (sigmaIetaIeta < 0.0353) &&
      (dEtaSeed < 0.00501) &&
      (dPhiIn < 0.0236) &&
      (HoE < 0.0188 + 2.06/energy + 0.183*rho/energy) &&
      (eInverseMinusPInverse < 0.0197) &&
      (mHits <= 1) &&
      (isPassConVeto == true);
    isTightRelIso = (relIsoWithEffectiveArea < 0.0445 + 0.963/elePt);
        } // end else (fabs(eleEta) > 1.479)

        if (electron.pt() < 7 || electron.eta() > 2.5 || !isLoose) continue;
        ElectronInfo e;
        e.pt = electron.pt();
        e.eta = electron.eta();
        e.phi = electron.phi();
        e.mass = electron.mass();
        e.charge = electron.charge();
        if (isTight) {
    e.id = 3;
        } else if (isMedium) {
    e.id = 2;
        } else {
    e.id = 1;
        }
        if (isTightRelIso) {
    e.iso = 3;
        } else if (isMediumRelIso) {
    e.iso = 2;
        } else if (isLooseRelIso) {
    e.iso = 1;
        } else {
    e.iso = 0;
        }
        e.lowptid = -9999;
        //std::cout << e.iso << '\n';
        math::XYZPointF p1 = electron.trackPositionAtVtx();
        math::XYZPoint p2 = PrimaryVertex.position();
        float dz = abs(p1.z()-p2.z());
        float dxy = sqrt((p1.x()-p2.x())*(p1.x()-p2.x()) + (p1.y()-p2.y())*(p1.y()-p2.y()));
        e.dxy = dxy;
        e.dz = dz;
        e.dB=electron.dB();
        e.edB = electron.edB();
        electronInfoData->push_back(e);
      }
    }
  /*
    edm::Handle< std::vector<pat::Electron> > LowPtElectronsHandle;
    iEvent.getByToken(LowPtElectrons_, LowPtElectronsHandle);
    auto LowPtElectrons = *LowPtElectronsHandle;

    if (LowPtElectrons.size() > 0) {
      for (unsigned int i = 0; i < LowPtElectrons.size(); ++i) {
        auto electron = LowPtElectrons[i];
        
        if ( electron.pt() < 1 || electron.eta() > 2.5 || electron.electronID("ID") < std::stof(idScoreCut_) ) continue;
        ElectronInfo e;
        e.pt = electron.pt();
        e.eta = electron.eta();
        e.phi = electron.phi();
        e.mass = electron.mass();
        e.charge = electron.charge();
        e.id = -9999;
        e.lowptid = electron.electronID("ID");
        //std::cout << e.iso << '\n';
        math::XYZPointF p1 = electron.trackPositionAtVtx();
        math::XYZPoint p2 = PrimaryVertex.position();
        float dz = abs(p1.z()-p2.z());
        float dxy = sqrt((p1.x()-p2.x())*(p1.x()-p2.x()) + (p1.y()-p2.y())*(p1.y()-p2.y()));
        e.dxy = dxy;
        e.dz = dz;
        lowPtElectronInfoData->push_back(e);
      }
    }
  */
    edm::Handle< std::vector<pat::Tau> > TausUnCleanedHandle;
    iEvent.getByToken(TausUnCleaned_, TausUnCleanedHandle);
    auto TausUnCleaned = *TausUnCleanedHandle;
    fillTauInfoDS(TausUnCleaned, 1);

    edm::Handle< std::vector<pat::Tau> > TausECleanedHandle;
    iEvent.getByToken(TausECleaned_, TausECleanedHandle);
    auto TausECleaned = *TausECleanedHandle;
    fillTauInfoDS(TausECleaned, 2);

    edm::Handle< std::vector<pat::Tau> > TausMCleanedHandle;
    iEvent.getByToken(TausMCleaned_, TausMCleanedHandle);
    auto TausMCleaned = *TausMCleanedHandle;
    fillTauInfoDS(TausMCleaned, 3);

    /*edm::Handle< std::vector<pat::Tau> > TausBoostedHandle;
    iEvent.getByToken(TausBoosted_, TausBoostedHandle);
    auto TausBoosted = *TausBoostedHandle;
    fillTauInfoDS(TausBoosted, 4);
    */
    /*edm::Handle< std::vector<pat::Tau> > TausLowPtECleanedHandle;
    iEvent.getByToken(TausLowPtECleaned_, TausLowPtECleanedHandle);
    auto TausLowPtECleaned = *TausLowPtECleanedHandle;
    fillTauInfoDS(TausLowPtECleaned, 5);
  */
    edm::Handle< std::vector<pat::MET> > METHandle;
    iEvent.getByToken(MET_, METHandle);
    auto met = METHandle->front();

    metInfo_.pt = met.pt();
    metInfo_.phi = met.phi();
    metInfo_.ptUncor = met.uncorPt();
    metInfo_.phiUncor = met.uncorPhi();
    metInfo_.ptJECUp = met.shiftedPt(pat::MET::JetEnUp);
    metInfo_.phiJECUp = met.shiftedPhi(pat::MET::JetEnUp);
    metInfo_.ptJERUp = met.shiftedPt(pat::MET::JetResUp);
    metInfo_.phiJERUp = met.shiftedPhi(pat::MET::JetResUp);
    metInfo_.ptUncUp = met.shiftedPt(pat::MET::UnclusteredEnUp);
    metInfo_.phiUncUp = met.shiftedPhi(pat::MET::UnclusteredEnUp);
    metInfo_.ptJECDown = met.shiftedPt(pat::MET::JetEnDown);
    metInfo_.phiJECDown = met.shiftedPhi(pat::MET::JetEnDown);
    metInfo_.ptJERDown = met.shiftedPt(pat::MET::JetResDown);
    metInfo_.phiJERDown = met.shiftedPhi(pat::MET::JetResDown);
    metInfo_.ptUncDown = met.shiftedPt(pat::MET::UnclusteredEnDown);
    metInfo_.phiUncDown = met.shiftedPhi(pat::MET::UnclusteredEnDown);

    tree->Fill();
  //}
}

void TCPNtuples::fillTauInfoDS(const std::vector<pat::Tau>& Taus, int whichColl) {
  if (Taus.size() > 0) {
    for (unsigned int i = 0; i < Taus.size(); ++i) {
      auto tau = Taus[i];
      if (tau.pt() < 10 || abs(tau.eta()) > 2.4) continue; //decayModeFinding > 0.5 by default
      if (!tau.tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") and !tau.tauID("byVVVLooseDeepTau2017v2p1VSjet")) continue;
      TauInfo t;
      t.pt = tau.pt();
      t.eta = tau.eta();
      t.phi = tau.phi();
      t.mass = tau.mass();
      t.charge=tau.charge();
      t.decaymode = tau.decayMode();
      t.mvaidraw = tau.tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017");
      t.deepidraw = tau.tauID("byDeepTau2017v2p1VSeraw");
      //t.dxy = abs(tau.leadChargedHadrCand().get()->dxy(PV));
      //t.dz = abs(tau.leadChargedHadrCand().get()->dz(PV));
      if (tau.tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017")) t.mvaid = 7;
      else if (tau.tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017")) t.mvaid = 6;
      else if (tau.tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017")) t.mvaid = 5;
      else if (tau.tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017")) t.mvaid = 4;
      else if (tau.tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017")) t.mvaid = 3;
      else if (tau.tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017")) t.mvaid = 2;
      else if (tau.tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017")) t.mvaid = 1;
      else t.mvaid = -1;
      if (tau.tauID("byVVTightDeepTau2017v2p1VSjet")) t.deepid = 7;
      else if (tau.tauID("byVTightDeepTau2017v2p1VSjet")) t.deepid = 6;
      else if (tau.tauID("byTightDeepTau2017v2p1VSjet")) t.deepid = 5;
      else if (tau.tauID("byMediumDeepTau2017v2p1VSjet")) t.deepid = 4;
      else if (tau.tauID("byLooseDeepTau2017v2p1VSjet")) t.deepid = 3;
      else if (tau.tauID("byVLooseDeepTau2017v2p1VSjet")) t.deepid = 2;
      else if (tau.tauID("byVVLooseDeepTau2017v2p1VSjet")) t.deepid = 1;
      else if (tau.tauID("byVVVLooseDeepTau2017v2p1VSjet")) t.deepid = 0;
      else t.deepid = -1;
      if (whichColl == 1) tauInfoDataUnCleaned->push_back(t);
      if (whichColl == 2) tauInfoDataECleaned->push_back(t);
      if (whichColl == 3) tauInfoDataMCleaned->push_back(t);
      if (whichColl == 4) tauInfoDataBoosted->push_back(t);
      //if (whichColl == 5) tauInfoDataLowPtECleaned->push_back(t);
    }
  }
}

float TCPNtuples::deltaR(float phi1, float phi2, float eta1, float eta2) {
  
  const float dphi = reco::deltaPhi(phi1, phi2);
  const float deta = eta1 - eta2;
  const float dr = std::sqrt(deta*deta + dphi*dphi);
  return dr;
}
bool TCPNtuples::applyPileupJetID(const pat::Jet & jet, const std::string & level){
  bool passpuid = false;  
  if(jet.hasUserInt("pileupJetIdUpdated:fullId")){
    if(level == "loose"  and (bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 2)) or jet.pt() > 50)) passpuid = true;
    if(level == "medium" and (bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 1)) or jet.pt() > 50)) passpuid = true;
    if(level == "tight"  and (bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 0)) or jet.pt() > 50)) passpuid = true;
  }
  else if (jet.userInt("pileupJetId:fullId")){
    if(level == "loose"  and (bool(jet.userInt("pileupJetId:fullId") & (1 << 2)) or jet.pt() > 50)) passpuid = true;
    if(level == "medium" and (bool(jet.userInt("pileupJetId:fullId") & (1 << 1)) or jet.pt() > 50)) passpuid = true;
    if(level == "tight"  and (bool(jet.userInt("pileupJetId:fullId") & (1 << 0)) or jet.pt() > 50)) passpuid = true;
  }
  return passpuid;
}	
//define this as a plug-in
DEFINE_FWK_MODULE(TCPNtuples);
  
