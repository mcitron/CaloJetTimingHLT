// -*- C++ -*-
//
// Package:    CaloJetTiming/CaloJetTimingFilter
// Class:      CaloJetTimingFilter
//
/**\class CaloJetTimingFilter CaloJetTimingFilter.cc CaloJetTiming/CaloJetTimingFilter/plugins/CaloJetTimingFilter.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Matthew Daniel Citron
//         Created:  Thu, 15 Jul 2021 03:19:13 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
//
// class declaration
//
namespace edm {
  class ConfigurationDescriptions;
}

class CaloJetTimingFilter : public HLTFilter {
    public:
	explicit CaloJetTimingFilter(const edm::ParameterSet& iConfig);
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs& filterproduct) const;
    private:


	//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	edm::InputTag jetLabel_;
	edm::InputTag jetTimeLabel_;
        unsigned int minJets_;
        double timeThresh_;
        double minPt_;
	edm::EDGetTokenT<reco::CaloJetCollection> jetInputToken;
	edm::EDGetTokenT<edm::ValueMap<float>> jetTimesInputToken;

	// ----------member data ---------------------------
#ifdef THIS_IS_AN_EVENT_EXAMPLE
	edm::EDGetTokenT<ExampleData> exampleToken_;
#endif
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
CaloJetTimingFilter::CaloJetTimingFilter(const edm::ParameterSet& iConfig) : HLTFilter(iConfig){
    jetLabel_= iConfig.getParameter<edm::InputTag>("jets");
    jetTimeLabel_= iConfig.getParameter<edm::InputTag>("jetTimes");
    minJets_= iConfig.getParameter<unsigned int>("minJets");
    timeThresh_ = iConfig.getParameter<double>("timeThresh");
    minPt_ = iConfig.getParameter<double>("minJetPt");
    jetInputToken = consumes<std::vector<reco::CaloJet>>(jetLabel_);
    jetTimesInputToken = consumes<edm::ValueMap<float>>(jetTimeLabel_);
    //now do what ever initialization is needed
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool CaloJetTimingFilter::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup,trigger::TriggerFilterObjectWithRefs& filterproduct) const {
    bool accept = false;
    int ijet = 0;
    edm::Handle<reco::CaloJetCollection> jets;
    iEvent.getByToken(jetInputToken, jets);
    edm::Handle<edm::ValueMap<float>> jetTimes;
    iEvent.getByToken(jetTimesInputToken, jetTimes);
    unsigned int njets = 0;
    for (auto const& c : *jets) {
	reco::CaloJetRef calojetref(jets, ijet);
	if((*jetTimes)[calojetref] > timeThresh_ && c.pt() > minPt_ && fabs(c.eta()) < 1.48) njets++;
	ijet ++;
    }
    accept = njets >= minJets_;
    return accept;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CaloJetTimingFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    makeHLTFilterDescription(desc);
    desc.add<edm::InputTag>("jets", edm::InputTag("hltDisplacedHLTCaloJetCollectionProducerMidPt"));
    desc.add<edm::InputTag>("jetTimes", edm::InputTag("hltDisplacedHLTCaloJetCollectionProducerMidPtTiming"));
    desc.add<unsigned int>("minJets", 1);
    desc.add<double>("timeThresh", 1.);
    desc.add<double>("minJetPt", 40.);
    descriptions.add("caloJetTimingFilter", desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CaloJetTimingFilter);
