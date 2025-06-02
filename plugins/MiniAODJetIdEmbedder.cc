/*
 * Embed PF Jet IDs (see https://twiki.cern.ch/twiki/bin/view/CMS/JetID)
 * into pat::Jets
 *
 * Author: Evan K. Friis, UW Madison
 */


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "FWCore/Utilities/interface/StreamID.h"
// CMSSW data formats
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include "PNet_Skimmer/interface/json.hpp"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace cms::Ort;
using namespace btagbtvdeep;

class MiniAODJetIdEmbedder : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
  public:
    MiniAODJetIdEmbedder(const edm::ParameterSet& pset, const ONNXRuntime *);
    virtual ~MiniAODJetIdEmbedder(){}
    void produce(edm::Event& evt, const edm::EventSetup& es);
    static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &);
    static void globalEndJob(const ONNXRuntime *);

  private:
    edm::EDGetTokenT<edm::View<pat::Jet> > srcToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> svToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> ditau2017v1Token_;
    edm::EDGetTokenT<edm::ValueMap<float>> ditau2017MDv1Token_;
    edm::EDGetTokenT<edm::ValueMap<float>> ditau2017v2Token_;
    edm::EDGetTokenT<edm::ValueMap<float>> ditau2017MDv2Token_;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;

    std::vector<std::string> input_names_; // names of each input group - the ordering is important!
    std::vector<std::vector<int64_t>> input_shapes_; // shapes of each input group (-1 for dynamic axis)
    std::vector<unsigned> input_sizes_; // total length of each input vector
    std::unordered_map<std::string, PreprocessParams> prep_info_map_; // preprocessing info for each input group
    FloatArrays data_; // each stream hosts its own data

    const reco::Vertex *pv_ = nullptr;
    std::vector<float> jet_pfcand_pt;
    std::vector<float> jet_pfcand_pt_log;
    std::vector<float> jet_pfcand_mass;
    std::vector<float> jet_pfcand_energy_log;
    std::vector<float> jet_pfcand_energy;
    std::vector<float> jet_pfcand_calofraction;
    std::vector<float> jet_pfcand_hcalfraction;
    std::vector<float> jet_pfcand_dz;
    std::vector<float> jet_pfcand_dxy;
    std::vector<float> jet_pfcand_dzsig;
    std::vector<float> jet_pfcand_dxysig;
    std::vector<float> jet_pfcand_frompv;
    std::vector<float> jet_pfcand_nstriphits;
    std::vector<float> jet_pfcand_npixhits;
    std::vector<unsigned int> jet_pfcand_id;
    std::vector<float> jet_pfcand_nlostinnerhits;
    std::vector<float> jet_pfcand_charge;
    std::vector<float> jet_pfcand_puppiw;
    std::vector<float> jet_pfcand_pperp_ratio;
    std::vector<float> jet_pfcand_ppara_ratio;
    std::vector<float> jet_pfcand_deta;
    std::vector<float> jet_pfcand_eta;
    std::vector<float> jet_pfcand_dphi;
    std::vector<float> jet_pfcand_phi;
    std::vector<float> jet_pfcand_etarel;


    std::vector<float> jet_pfcand_track_chi2;
    std::vector<float> jet_pfcand_track_qual;

    std::vector<float> jet_pfcand_trackjet_dxy;
    std::vector<float> jet_pfcand_trackjet_dxysig;
    std::vector<float> jet_pfcand_trackjet_d3d;
    std::vector<float> jet_pfcand_trackjet_d3dsig;
    std::vector<float> jet_pfcand_trackjet_dist;
    std::vector<float> jet_pfcand_trackjet_decayL;


    // jet-to-secondary vertex
    std::vector<float> jet_sv_pt;
    std::vector<float> jet_sv_pt_log;
    std::vector<float> jet_sv_eta;
    std::vector<float> jet_sv_phi;
    std::vector<float> jet_sv_deta;
    std::vector<float> jet_sv_dphi;
    std::vector<float> jet_sv_mass;
    std::vector<float> jet_sv_energy;
    std::vector<float> jet_sv_energy_log;
    std::vector<float> jet_sv_chi2;
    std::vector<float> jet_sv_dxy;
    std::vector<float> jet_sv_dxysig;
    std::vector<float> jet_sv_d3d;
    std::vector<float> jet_sv_d3dsig;
    std::vector<float> jet_sv_ntrack;
    std::vector<unsigned int> jet_sv_ijet;

    const std::vector<std::string>   pf_points_{"jet_pfcand_deta",
                                                "jet_pfcand_dphi"};
    const std::vector<std::string> pf_features_{"jet_pfcand_pt_log",
                                                "jet_pfcand_mass",
                                                "jet_pfcand_energy_log",
                                                "jet_pfcand_track_qual",
                                                "jet_pfcand_track_chi2",
                                                "jet_pfcand_trackjet_d3d",
                                                "jet_pfcand_trackjet_d3dsig",
                                                "jet_pfcand_trackjet_dist",
                                                "jet_pfcand_trackjet_decayL",
                                                "jet_pfcand_dz",
                                                "jet_pfcand_dxy",
                                                "jet_pfcand_dxysig",
                                                "jet_pfcand_dzsig",
                                                "jet_pfcand_nlostinnerhits",
                                                "jet_pfcand_nstriphits",
                                                "jet_pfcand_npixhits",
                                                "jet_pfcand_frompv",
                                                "jet_pfcand_ppara_ratio",
                                                "jet_pfcand_pperp_ratio",
                                                "jet_pfcand_highpurity",
                                                "jet_pfcand_puppiw",
                                                "jet_pfcand_charge"
                                                };
    const std::vector<std::string>     pf_mask_{"pfcand_mask"};
    const std::vector<std::string>   sv_points_{"jet_sv_dphi",
                                                "jet_sv_deta"};
    const std::vector<std::string> sv_features_{
                                                "jet_sv_pt_log",
                                                "jet_sv_mass",
                                                "jet_sv_energy_log",
                                                "jet_sv_d3d",
                                                "jet_sv_d3dsig",
                                                "jet_sv_ntrack",
                                                "jet_sv_chi2",
                                                "jet_sv_dxy",
                                                "jet_sv_dxysig"
                                                };
    const std::vector<std::string>     sv_mask_{"sv_mask"};
    // Sorters to order object collections in decreasing order of pT
    template<typename T> 
    class PatPtSorter {
    public:
      bool operator()(const T& i, const T& j) const {
        return (i.pt() > j.pt());
      }
    };

    PatPtSorter<pat::Muon>     muonSorter;
    PatPtSorter<pat::Electron> electronSorter;
    PatPtSorter<pat::Tau>      tauSorter;
    PatPtSorter<pat::PackedCandidate>  packedPFCandidateSorter;
};

int my_center_norm_pad(const std::vector<float> &input,
                    float center,
                    float norm_factor,
                    unsigned min_length,
                    unsigned max_length,
                    std::vector<float> &datavec,
                    int startval,
                    float pad_value,
                    float replace_inf_value,
                    float min,
                    float max) {
  assert(min <= pad_value && pad_value <= max);
  assert(min_length <= max_length);

  unsigned target_length = std::clamp((unsigned)input.size(), min_length, max_length);
  for (unsigned i = 0; i < target_length; ++i) {
    if (i < input.size()) {
      datavec[i + startval] = std::clamp((btagbtvdeep::catch_infs(input[i], replace_inf_value) - center) * norm_factor, min, max);
    } else {
      datavec[i + startval] = pad_value;
    }
  }
  return target_length;
}

MiniAODJetIdEmbedder::MiniAODJetIdEmbedder(const edm::ParameterSet& pset, const ONNXRuntime *cache):
track_builder_token_(
          esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder")))
 {
  srcToken_ = consumes<edm::View<pat::Jet>>(pset.getParameter<edm::InputTag>("src"));
  //genParticlesToken_ = consumes<edm::View<reco::GenParticle>>(pset.getParameter<edm::InputTag>("genParticles"));
  vtxToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertices"));
  svToken_ = consumes<edm::View<reco::VertexCompositePtrCandidate>>(pset.getParameter<edm::InputTag>("secondaryVertices"));

  auto json_path = pset.getParameter<edm::FileInPath>("preprocess_json");

  input_names_.clear(); input_shapes_.clear();
  std::ifstream ifs(edm::FileInPath(json_path).fullPath());
  nlohmann::json js = nlohmann::json::parse(ifs);
  js.at("input_names").get_to(input_names_);
  for (const auto &group_name : input_names_) {
    const auto &group_pset = js.at(group_name);
    auto &prep_params = prep_info_map_[group_name];
    group_pset.at("var_names").get_to(prep_params.var_names);
    if (group_pset.contains("var_length")) {
      prep_params.min_length = group_pset.at("var_length");
      prep_params.max_length = prep_params.min_length;
      input_shapes_.push_back({1, (int64_t)prep_params.var_names.size(), prep_params.min_length});
    } else {
      prep_params.min_length = group_pset.at("min_length");
      prep_params.max_length = group_pset.at("max_length");
      input_shapes_.push_back({1, (int64_t)prep_params.var_names.size(), -1});
    }
    const auto &var_info_pset = group_pset.at("var_infos");
    for (const auto &var_name : prep_params.var_names) {
      const auto &var_pset = var_info_pset.at(var_name);
      double median = var_pset.at("median");
      double norm_factor = var_pset.at("norm_factor");
      double replace_inf_value = var_pset.at("replace_inf_value");
      double lower_bound = var_pset.at("lower_bound");
      double upper_bound = var_pset.at("upper_bound");
      double pad = var_pset.contains("pad") ? double(var_pset.at("pad")) : 0;
      prep_params.var_info_map[var_name] =
          PreprocessParams::VarInfo(median, norm_factor, replace_inf_value, lower_bound, upper_bound, pad);
    }
    if (&data_ != nullptr) {
      const auto &len = input_sizes_.emplace_back(prep_params.max_length * prep_params.var_names.size());
      data_.emplace_back(len, 0);
    }
  }
  if (pset.exists("ditau2017v1")) { ditau2017v1Token_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ditau2017v1")); }
  if (pset.exists("ditau2017MDv1")) { ditau2017MDv1Token_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ditau2017MDv1")); }
  if (pset.exists("ditau2017v2")) { ditau2017v2Token_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ditau2017v2")); }
  if (pset.exists("ditau2017MDv2")) { ditau2017MDv2Token_ = consumes<edm::ValueMap<float> >(pset.getParameter<edm::InputTag>("ditau2017MDv2")); }

  produces<pat::JetCollection>("slimmedJetsWIDs").setBranchAlias( "slimmedJetsWIDs");
}

std::unique_ptr<ONNXRuntime> MiniAODJetIdEmbedder::initializeGlobalCache(const edm::ParameterSet &iConfig) {
  return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void MiniAODJetIdEmbedder::globalEndJob(const ONNXRuntime *cache) {}

void MiniAODJetIdEmbedder::produce(edm::Event& evt, const edm::EventSetup& es) {

  //edm::Handle<edm::TriggerResults> triggerResultsHandle;
  //evt.getByToken(TriggerResults_, triggerResultsHandle);
  //TriggerResults triggerResults = *triggerResultsHandle;
  //auto & names = evt.triggerNames(*triggerResultsHandle);
  //if (
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part0_v3"))||
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part1_v3"))||
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part2_v3"))||
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part3_v3"))||
  //  triggerResults.accept(names.triggerIndex("HLT_Mu9_IP6_part4_v3"))
  //) {
    std::unique_ptr<pat::JetCollection> output(new pat::JetCollection);

    edm::Handle<edm::View<pat::Jet> > input;
    evt.getByToken(srcToken_, input);

    output->reserve(input->size());

    edm::Handle<reco::VertexCollection> vtxHandle;
    evt.getByToken(vtxToken_, vtxHandle);
    pv_ = &vtxHandle->at(0);

    edm::Handle<edm::View<reco::VertexCompositePtrCandidate> > svHandle;
    evt.getByToken(svToken_, svHandle);


    edm::ESHandle<TransientTrackBuilder> track_builder_ = es.getHandle(track_builder_token_);

    TrackInfoBuilder trackinfo(track_builder_);


    edm::Handle<edm::ValueMap<float> > ditau2017v1;
    bool ditau2017v1Valid = evt.getByToken(ditau2017v1Token_, ditau2017v1);

    edm::Handle<edm::ValueMap<float> > ditau2017MDv1;
    bool ditau2017MDv1Valid = evt.getByToken(ditau2017MDv1Token_, ditau2017MDv1);

    edm::Handle<edm::ValueMap<float> > ditau2017v2;
    bool ditau2017v2Valid = evt.getByToken(ditau2017v2Token_, ditau2017v2);

    edm::Handle<edm::ValueMap<float> > ditau2017MDv2;
    bool ditau2017MDv2Valid = evt.getByToken(ditau2017MDv2Token_, ditau2017MDv2);




    for (size_t i = 0; i < input->size(); ++i) {
      pat::Jet jet = input->at(i);
      TVector3 jet_direction      (jet.momentum().Unit().x(),jet.momentum().Unit().y(),jet.momentum().Unit().z());
      GlobalVector jet_global_vec (jet.px(),jet.py(),jet.pz());	     

      // https://twiki.cern.ch/twiki/bin/view/CMS/JetID
      bool loose = true;
      bool tight = true;
      bool tightLepVeto = true;
      if (std::abs(jet.eta()) <= 2.6) {
        if (jet.neutralHadronEnergyFraction() >= 0.99) {
          loose = false;
        }
        if (jet.neutralHadronEnergyFraction() >= 0.90) {
          tight = false;
          tightLepVeto = false;
        }
        if (jet.neutralEmEnergyFraction() >= 0.99) {
          loose = false;
        }
        if (jet.neutralEmEnergyFraction() >= 0.90) {
          tight = false;
          tightLepVeto = false;
        }
        if (jet.chargedMultiplicity()+jet.neutralMultiplicity() <= 1){
          loose = false;
          tight = false;
          tightLepVeto = false;
        }
        if (jet.muonEnergyFraction() >= 0.8)
          {
            tightLepVeto = false;
          }
        if (jet.chargedHadronEnergyFraction() == 0) {
          loose = false;
          tight = false;
          tightLepVeto = false;
        }
        if (jet.chargedHadronMultiplicity() == 0) {
          loose = false;
          tight = false;
          tightLepVeto = false;
        }
        if (jet.chargedEmEnergyFraction() >= 0.80) {
          tightLepVeto = false;
        }
      }

      if (std::abs(jet.eta()) >2.6 && std::abs(jet.eta()) <= 2.7) {
        if (jet.neutralHadronEnergyFraction() >= 0.99) {
          loose = false;
        }
        if (jet.neutralHadronEnergyFraction() >= 0.90) {
          tight = false;
          tightLepVeto = false;
        }
        if (jet.neutralEmEnergyFraction() >= 0.99) {
          loose = false;
          tight = false;
          tightLepVeto = false;
        }
        if (jet.muonEnergyFraction() >= 0.8){
          tightLepVeto = false;
        }
        if (jet.chargedHadronMultiplicity() == 0) {
          loose = false;
          tight = false;
          tightLepVeto = false;
        }
        if (jet.chargedEmEnergyFraction() >= 0.80) {
          tightLepVeto = false;
        }
      }

      if (std::abs(jet.eta()) > 2.7 && std::abs(jet.eta()) <= 3.0) {
        if (jet.neutralEmEnergyFraction() >= 0.99 or jet.neutralEmEnergyFraction()<=0.02) {
          loose = false;
          tight = false;
        }
        if (jet.neutralMultiplicity()<=2) {
          loose = false;
          tight = false;
        }
      }

      if (std::abs(jet.eta()) > 3.0) {
        if (jet.neutralEmEnergyFraction() >= 0.90) {
          loose = false;
          tight = false;
        }
        if (jet.neutralHadronEnergyFraction() <= 0.2) {
          tight = false;
          tightLepVeto = false;
        }
        if (jet.neutralMultiplicity()<=10) {
          loose = false;
          tight = false;
        }
      }
      jet.addUserFloat("idLoose", loose);
      jet.addUserFloat("idTight", tight);
      jet.addUserFloat("idTightLepVeto", tightLepVeto);

      // Pileup discriminant
      bool passPU = true;
      float jpumva = jet.userFloat("pileupJetId:fullDiscriminant");
      if(jet.pt() > 20)
        {
          if(fabs(jet.eta()) > 3.)
            {
              if(jpumva <= -0.45) passPU = false;
            }
          else if(fabs(jet.eta()) > 2.75)
            {
              if(jpumva <= -0.55) passPU = false;
            }
          else if(fabs(jet.eta()) > 2.5)
            {
              if(jpumva <= -0.6) passPU = false;
            }
          else if(jpumva <= -0.63) passPU = false;
        }
      else
        {
          if(fabs(jet.eta()) > 3.)
            {
              if(jpumva <= -0.95) passPU = false;
            }
          else if(fabs(jet.eta()) > 2.75)
            {
              if(jpumva <= -0.94) passPU = false;
            }
          else if(fabs(jet.eta()) > 2.5)
            {
              if(jpumva <= -0.96) passPU = false;
            }
          else if(jpumva <= -0.95) passPU = false;
        }

      jet.addUserFloat("puID", float(passPU));

      // Add matching to merged gen bs

      // Add ParticleNet score
      float PNetScore = -1.;
      if(jet.pt() > 10. && abs(jet.eta()) < 2.4) {

        // create jet features
        DeepBoostedJetFeatures features;
        for (const auto &name : pf_points_)
          features.add(name);
        for (const auto &name : pf_features_)
          features.add(name);
        for (const auto &name : pf_mask_)
          features.add(name);
        for (const auto &name : sv_points_)
          features.add(name);
        for (const auto &name : sv_features_)
          features.add(name);
        for (const auto &name : sv_mask_)
          features.add(name);

        // build trackinfo
        GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());

        // reserve space

        int nConstituents = jet.numberOfSourceCandidatePtrs();
        for (const auto &name : pf_points_)
          features.reserve(name, nConstituents);
        for (const auto &name : pf_features_)
          features.reserve(name, nConstituents);
        for (const auto &name : pf_mask_)
          features.reserve(name, nConstituents);


        // fill pfcand features

        for(int k = 0; k < nConstituents; k++){
          reco::CandidatePtr pfcand_ptr = jet.sourceCandidatePtr(k);

          const auto *pfcand = dynamic_cast<const pat::PackedCandidate *>(&(*pfcand_ptr));
          if(pfcand->pt() < .1) continue;
          if(pfcand->puppiWeight()<0) continue;
                  //if (pfcand->hasTrackDetails()==false){////std::cout<<"HASTrack\n";
                  //continue;
                  //}

          features.fill("pfcand_mask", 1);
          TVector3 pfcand_momentum (pfcand->momentum().x(),pfcand->momentum().y(),pfcand->momentum().z());
          ////std::cout<<"pfcand Momentum\n";
          features.fill("jet_pfcand_pt_log", btagbtvdeep::catch_infs(std::isnan(std::log(pfcand->pt())) ? 0 : std::log(pfcand->pt())));
          //features.fill("jet_pfcand_pt",pfcand->pt());
          features.fill("jet_pfcand_energy_log", btagbtvdeep::catch_infs(std::isnan(std::log(pfcand->energy())) ? 0 : std::log(pfcand->energy())));
          //features.fill("jet_pfcand_energy", pfcand->energy());
          //features.fill("jet_pfcand_calofraction", pfcand->caloFraction());
          //features.fill("jet_pfcand_hcalfraction", pfcand->hcalFraction());
          //features.fill("jet_pfcand_id", pfcand->pdgId());
          //features.fill("jet_pfcand_phi", pfcand->phi());
          features.fill("jet_pfcand_highpurity", pfcand->trackHighPurity());
          features.fill("jet_pfcand_frompv", pfcand->fromPV());
          features.fill("jet_pfcand_nstriphits", pfcand->stripLayersWithMeasurement());
          features.fill("jet_pfcand_npixhits", pfcand->numberOfPixelHits());
          features.fill("jet_pfcand_dphi", jet_direction.DeltaPhi(pfcand_momentum));
          features.fill("jet_pfcand_mass", pfcand->mass());
          features.fill("jet_pfcand_deta", (jet_direction.Eta()-pfcand_momentum.Eta()));
          features.fill("jet_pfcand_puppiw", pfcand->puppiWeight());
          features.fill("jet_pfcand_charge", pfcand->charge());
          features.fill("jet_pfcand_nlostinnerhits", pfcand->lostInnerHits());
          //features.fill("jet_pfcand_etarel",reco::btau::etaRel(jetDir,pfcand->momentum()));

          // impact parameters
          //features.fill("jet_pfcand_dz", pfcand->dz(vtxHandle->front().position()));
          //features.fill("jet_pfcand_dzsig", fabs(pfcand->dz(vtxHandle->front().position()))/pfcand->dzError());
          //features.fill("jet_pfcand_dxy", btagbtvdeep::catch_infs(pfcand->dxy(vtxHandle->front().position())));
          //features.fill("jet_pfcand_dxysig", fabs(pfcand->dxy(vtxHandle->front().position()))/pfcand->dxyError());
                  // impact parameters
          features.fill("jet_pfcand_dz", btagbtvdeep::catch_infs(pfcand->dz()));
          features.fill("jet_pfcand_dzsig", fabs(pfcand->bestTrack() ? btagbtvdeep::catch_infs(pfcand->dz() / pfcand->dzError()) : 0));
          features.fill("jet_pfcand_dxy", btagbtvdeep::catch_infs(pfcand->dxy()));
          features.fill("jet_pfcand_dxysig", fabs(pfcand->bestTrack() ? btagbtvdeep::catch_infs(pfcand->dxy() / pfcand->dxyError()) : 0));
          
          features.fill("jet_pfcand_pperp_ratio", jet_direction.Perp(pfcand_momentum)/pfcand_momentum.Mag());
          features.fill("jet_pfcand_ppara_ratio", jet_direction.Dot(pfcand_momentum)/pfcand_momentum.Mag());
          if (pfcand->bestTrack()) {
            const auto* trk = pfcand->bestTrack();
            features.fill("jet_pfcand_track_chi2", btagbtvdeep::catch_infs(trk->normalizedChi2()));
            features.fill("jet_pfcand_track_qual", trk->qualityMask());

            reco::TransientTrack transientTrack = track_builder_->build(*trk);
            Measurement1D meas_ip2d    = IPTools::signedTransverseImpactParameter(transientTrack,jet_global_vec,*pv_).second;
            Measurement1D meas_ip3d    = IPTools::signedImpactParameter3D(transientTrack,jet_global_vec,*pv_).second;
            Measurement1D meas_jetdist = IPTools::jetTrackDistance(transientTrack,jet_global_vec,*pv_).second;
            Measurement1D meas_decayl  = IPTools::signedDecayLength3D(transientTrack,jet_global_vec,*pv_).second;
            features.fill("jet_pfcand_trackjet_d3d", meas_ip3d.value());
            features.fill("jet_pfcand_trackjet_d3dsig", fabs(meas_ip3d.significance()));
            features.fill("jet_pfcand_trackjet_dist", -meas_jetdist.value());
            features.fill("jet_pfcand_trackjet_decayL", meas_decayl.value());
                    ////std::cout<<"pfcand trackjet\n";

          }
          else {
            features.fill("jet_pfcand_track_chi2", 0.);
            features.fill("jet_pfcand_track_qual", 0.);
            features.fill("jet_pfcand_trackjet_d3d", 0.);
            features.fill("jet_pfcand_trackjet_d3dsig", 0.);
            features.fill("jet_pfcand_trackjet_dist", 0.);
            features.fill("jet_pfcand_trackjet_decayL", 0.);
          }
        }
              ////std::cout<<"Fill cands\n";
              ////std::cout<<"Fill cands\n";

        // find SVs associated with the jet
        std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
        for(const auto &sv : *svHandle){
          if(reco::deltaR(jet, sv) < 0.4){
            jetSVs.push_back(&sv);
          }
        }

        // sort by dxy significance
        std::sort(jetSVs.begin(), jetSVs.end(), [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) { return btagbtvdeep::sv_vertex_comparator(*sva, *svb, *pv_); });

        // reserve space
        int nsv = jetSVs.size();
        for (const auto &name : sv_points_)
          features.reserve(name, nsv);
        for (const auto &name : sv_features_)
          features.reserve(name, nsv);
        for (const auto &name : sv_mask_)
          features.reserve(name, nsv);

        // fill associated secondary vertices information
        for (const auto *jetsv : jetSVs) {
                  ////std::cout<<"SVBegindxdydzetc\n";

          features.fill("sv_mask", 1);
          features.fill("jet_sv_pt_log", std::isnan(std::log(jetsv->pt())) ? 0 : std::log(jetsv->pt()));
          features.fill("jet_sv_energy_log", std::isnan(std::log(jetsv->energy())) ? 0 : std::log(jetsv->energy()));
          //features.fill("jet_sv_eta", jetsv->eta());
          //features.fill("jet_sv_phi", jetsv->phi());        
          features.fill("jet_sv_deta", jetsv->eta()-jet.eta());
          features.fill("jet_sv_dphi", jetsv->phi()-jet.phi());
          features.fill("jet_sv_mass", jetsv->mass());
          features.fill("jet_sv_chi2", jetsv->vertexNormalizedChi2());
          reco::Vertex::CovarianceMatrix csv;
          jetsv->fillVertexCovariance(csv);
          reco::Vertex svtx(jetsv->vertex(), csv);
          VertexDistanceXY dxy;

          auto valxy = dxy.signedDistance(svtx, vtxHandle->at(0), jet_global_vec);
          features.fill("jet_sv_dxy", valxy.value());
          features.fill("jet_sv_dxysig", fabs(valxy.significance()));
          VertexDistance3D d3d;
          auto val3d = d3d.signedDistance(svtx, vtxHandle->at(0), jet_global_vec);
          jet_sv_d3d.push_back(val3d.value());
          features.fill("jet_sv_d3d", val3d.value());
          features.fill("jet_sv_d3dsig", fabs(val3d.significance()));
          features.fill("jet_sv_ntrack", jetsv->numberOfDaughters());
                  ////std::cout<<"SVs done\n";

        }

        features.check_consistency(pf_points_);
        features.check_consistency(pf_features_);
        features.check_consistency(pf_mask_);
        features.check_consistency(sv_points_);
        features.check_consistency(sv_features_);
        features.check_consistency(sv_mask_);
        //////std::cout<<"fill SVs\n";

        // make input
        // following: https://github.com/cms-sw/cmssw/blob/4e6da5e/RecoBTag/ONNXRuntime/plugins/BoostedJetONNXJetTagsProducer.cc#L196
        for (unsigned igroup = 0; igroup < input_names_.size(); ++igroup) {
          const auto &group_name = input_names_[igroup];
          const auto &prep_params = prep_info_map_.at(group_name);
          auto &group_values = data_[igroup];
          group_values.resize(input_sizes_[igroup]);
          // first reset group_values to 0
          std::fill(group_values.begin(), group_values.end(), 0);
          unsigned curr_pos = 0;
          // transform/pad
          for (unsigned i = 0; i < prep_params.var_names.size(); ++i) {
            const auto &varname = prep_params.var_names[i];
            const auto &raw_value = features.get(varname);
            const auto &info = prep_params.info(varname);
            int insize = my_center_norm_pad(raw_value,
                                        info.center,
                                        info.norm_factor,
                                        prep_params.min_length,
                                        prep_params.max_length,
                                        group_values,
                                        curr_pos,
                                        info.pad,
                                        info.replace_inf_value,
                                        info.lower_bound,
                                        info.upper_bound);
            curr_pos += insize;
            if (i == 0 && (!input_shapes_.empty())) {
              input_shapes_[igroup][2] = insize;
            }
          }
          group_values.resize(curr_pos);
        }

        // get output
        std::vector<float> outputs = globalCache()->run(input_names_, data_, input_shapes_)[0];
        PNetScore = outputs[0];
      }

      jet.addUserFloat("PNetScore", float(PNetScore));
      float ditau2017v1Value = 0;
      float ditau2017MDv1Value = 0;
      float ditau2017v2Value = 0;
      float ditau2017MDv2Value = 0;
      edm::Ref<edm::View<pat::Jet> > jRef(input,i);
      if (ditau2017v1Valid) {
        ditau2017v1Value = (*ditau2017v1)[jRef];

      }
      if (ditau2017MDv1Valid) ditau2017MDv1Value = (*ditau2017MDv1)[jRef];
      if (ditau2017v2Valid) ditau2017v2Value = (*ditau2017v2)[jRef];
      if (ditau2017MDv2Valid) ditau2017MDv2Value = (*ditau2017MDv2)[jRef];
      jet.addUserFloat("ditau2017v1",ditau2017v1Value);
      jet.addUserFloat("ditau2017MDv1",ditau2017MDv1Value);
      jet.addUserFloat("ditau2017v2",ditau2017v2Value);
      jet.addUserFloat("ditau2017MDv2",ditau2017MDv2Value);
      output->push_back(jet);
    }

    evt.put(std::move(output),"slimmedJetsWIDs");
  //}
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MiniAODJetIdEmbedder);
