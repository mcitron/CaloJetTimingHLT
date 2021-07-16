// -*- C++ -*-
//
// Package:    CaloJetTiming/CaloJetTimingProducer
// Class:      CaloJetTimingProducer
//
/**\class CaloJetTimingProducer CaloJetTimingProducer.cc CaloJetTiming/CaloJetTimingProducer/plugins/CaloJetTimingProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthew Daniel Citron
//         Created:  Wed, 14 Jul 2021 19:38:20 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"


//
// class declaration
//

class CaloJetTimingProducer : public edm::stream::EDProducer<> {
public:
  explicit CaloJetTimingProducer(const edm::ParameterSet&);
  ~CaloJetTimingProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::CaloJetCollection> jetInputToken;
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
CaloJetTimingProducer::CaloJetTimingProducer(const edm::ParameterSet& iConfig)
{
  //register your products
/* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
*/
    produces<edm::ValueMap<float>>("hltDisplacedHLTCaloJetCollectionProducerMidPtTiming");
  //now do what ever other initialization is needed
    jetInputToken = consumes<std::vector<reco::CaloJet>>(edm::InputTag("hltDisplacedHLTCaloJetCollectionProducerMidPt","","HLTX"));
}

CaloJetTimingProducer::~CaloJetTimingProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void CaloJetTimingProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    std::vector<float> jetTimings;
    Handle<reco::CaloJetCollection> jets;
    iEvent.getByToken(jetInputToken, jets); 
    for (auto const& c : *jets) {
	// for (EcalRecHitCollection::const_iterator i=ecalRecHits->begin(); i!=ecalRecHits->end(); i++) {
	// }
	jetTimings.push_back(c.pt());
    }
    std::unique_ptr<edm::ValueMap<float> > jetTimings_out(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler jetTimings_filler(*jetTimings_out);
    jetTimings_filler.insert(jets, jetTimings.begin(), jetTimings.end());
    jetTimings_filler.fill();
    int ijet = 0;
    for (auto const& c : *jets) {
	reco::CaloJetRef calojetref(jets, ijet);
	std::cout << (*jetTimings_out)[calojetref] <<std::endl;
	ijet ++;
    }
    iEvent.put(std::move(jetTimings_out), "hltDisplacedHLTCaloJetCollectionProducerMidPtTiming");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void CaloJetTimingProducer::beginStream(edm::StreamID) {
    // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void CaloJetTimingProducer::endStream() {
    // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   CaloJetTimingProducer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   CaloJetTimingProducer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   CaloJetTimingProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   CaloJetTimingProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CaloJetTimingProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloJetTimingProducer);
