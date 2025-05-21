// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "PhysicsTools/SelectorUtils/interface/CutApplicatorBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"


#include "TMath.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "CommonTools/Egamma/interface/ConversionTools.h"

//#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/RefToBaseVector.h"

using namespace edm;
using namespace std;
//
// class declaration
//

class LowPtElectronFilter : public edm::stream::EDFilter<> {
public:
      explicit LowPtElectronFilter(const edm::ParameterSet&);
  ~LowPtElectronFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  
  
private:
  edm::EDGetTokenT<pat::ElectronCollection> electronSrc_;
  virtual void beginStream(edm::StreamID) override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;

  string idScoreCut_;
  
};

LowPtElectronFilter::LowPtElectronFilter(const edm::ParameterSet& iConfig):
  electronSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  idScoreCut_(iConfig.getParameter<string>("LowPtEIdScoreCut"))

  //trk_(consumes<reco::GsfTrackCollection>(iConfig.getParameter<edm::InputTag>("Tracks")))
{
  //now do what ever initialization is needed
  produces<pat::ElectronCollection>( "MiniLowPtElectron" );
  produces<pat::ElectronRefVector>("LowPtElectronRef");
}


LowPtElectronFilter::~LowPtElectronFilter()
{
  
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  
}

// ------------ method called on each new Event  ------------
bool LowPtElectronFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;
 
  Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronSrc_,electrons);
  unique_ptr<pat::ElectronCollection> passedelectrons(new pat::ElectronCollection);
  unique_ptr<pat::ElectronRefVector> passedelectronRef(new pat::ElectronRefVector);
  Handle<reco::VertexCollection> Vertex;


  for(pat::ElectronCollection::const_iterator iele = electrons->begin() ; iele !=electrons->end(); ++iele)
    { 
      pat::ElectronRef ERef(electrons,iele-electrons->begin()); 

      if((iele->pt()>1) && (iele->electronID("ID")>std::stof(idScoreCut_)))

	    {


	      passedelectrons->push_back(*iele);
	      passedelectronRef->push_back(ERef);
	      //cout<< "EEPushback" <<endl;

	    }

	}


  iEvent.put(move(passedelectrons), "MiniLowPtElectron");
  iEvent.put(move(passedelectronRef),"LowPtElectronRef");
  return true;
  //else return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
LowPtElectronFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
LowPtElectronFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  ElectronFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ElectronFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ElectronFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ElectronFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LowPtElectronFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(LowPtElectronFilter);
