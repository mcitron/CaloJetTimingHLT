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

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

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
  edm::InputTag jetLabel_;
  edm::InputTag ecalEBLabel_;
  edm::InputTag ecalEELabel_;
  bool barrelOnly_;
  edm::EDGetTokenT<reco::CaloJetCollection> jetInputToken;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ecalRecHitsEBToken;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ecalRecHitsEEToken;
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
    produces<edm::ValueMap<float>>("");
  //now do what ever other initialization is needed
    jetLabel_= iConfig.getParameter<edm::InputTag>("jets");
    ecalEBLabel_= iConfig.getParameter<edm::InputTag>("ebRecHitsColl");
    ecalEELabel_= iConfig.getParameter<edm::InputTag>("eeRecHitsColl");
    barrelOnly_ = iConfig.getParameter<bool>("barrelOnly");
    jetInputToken = consumes<std::vector<reco::CaloJet>>(jetLabel_);
    ecalRecHitsEBToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(ecalEBLabel_);
    ecalRecHitsEEToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(ecalEELabel_);
    // jetInputToken = consumes<std::vector<reco::CaloJet>>(edm::InputTag("hltDisplacedHLTCaloJetCollectionProducerMidPt","","HLTX"));
    // ecalRecHitsEBToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(edm::InputTag("hltEcalRecHit","EcalRecHitsEB"));
    // ecalRecHitsEEToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(edm::InputTag("hltEcalRecHit","EcalRecHitsEE"));
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
    Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ecalRecHitsEB;
    iEvent.getByToken(ecalRecHitsEBToken,ecalRecHitsEB);
    Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ecalRecHitsEE;
    iEvent.getByToken(ecalRecHitsEEToken,ecalRecHitsEE);
    std::vector<bool> skipCellEB(ecalRecHitsEB->size(),false); 
    std::vector<bool> skipCellEE(ecalRecHitsEB->size(),false); 
    edm::ESHandle<CaloGeometry> pG; 
    iSetup.get<CaloGeometryRecord>().get(pG); 
    for (auto const& c : *jets) {
	int iCell = -1;
	float weightedTimeCell = 0;
	float totalEmEnergyCell = 0;
	unsigned int nCells = 0;
	for (EcalRecHitCollection::const_iterator i=ecalRecHitsEB->begin(); i!=ecalRecHitsEB->end(); i++) {
	    iCell++;
	    if (skipCellEB[iCell]) continue;
	    if ((i->checkFlag(EcalRecHit::kSaturated) || i->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || i->checkFlag(EcalRecHit::kPoorReco) || i->checkFlag(EcalRecHit::kWeird) || i->checkFlag(EcalRecHit::kDiWeird))) continue;
            if (i->energy() < 0.5) continue;
            if (i->timeError() < 0. || i->timeError() > 100) continue;
            if (i->time() < -12.5 || i->time() > 12.5) continue;
            GlobalPoint p=pG->getPosition(i->detid());
            if (reco::deltaR(c,p) > 0.4) continue;
	    weightedTimeCell += i->time()*i->energy()*sin(p.theta());
	    totalEmEnergyCell += i->energy()*sin(p.theta());
            nCells++;
	}
	iCell = -1;
	if (!barrelOnly_){
	    for (EcalRecHitCollection::const_iterator i=ecalRecHitsEE->begin(); i!=ecalRecHitsEE->end(); i++) {
		iCell++;
		if (skipCellEE[iCell]) continue;
		if ((i->checkFlag(EcalRecHit::kSaturated) || i->checkFlag(EcalRecHit::kLeadingEdgeRecovered) || i->checkFlag(EcalRecHit::kPoorReco) || i->checkFlag(EcalRecHit::kWeird) || i->checkFlag(EcalRecHit::kDiWeird))) continue;
		if (i->energy() < 0.5) continue;
		if (i->timeError() <= 0. || i->timeError() > 100) continue;
		if (i->time() < -12.5 || i->time() > 12.5) continue;
		GlobalPoint p=pG->getPosition(i->detid());
		if (reco::deltaR(c,p) > 0.4) continue;
		weightedTimeCell += i->time()*i->energy()*sin(p.theta());
		totalEmEnergyCell += i->energy()*sin(p.theta());
		nCells++;
	    }
	}
	if (totalEmEnergyCell > 10 && nCells > 5){
	    jetTimings.push_back(weightedTimeCell/totalEmEnergyCell);
	} 
	else{
	    jetTimings.push_back(-50);
	}
    }
    std::unique_ptr<edm::ValueMap<float> > jetTimings_out(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler jetTimings_filler(*jetTimings_out);
    jetTimings_filler.insert(jets, jetTimings.begin(), jetTimings.end());
    jetTimings_filler.fill();
    iEvent.put(std::move(jetTimings_out), "");
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
    desc.add<edm::InputTag>("jets", edm::InputTag(""));
    desc.add<bool>("barrelOnly", false);
    desc.add<edm::InputTag>("ebRecHitsColl", edm::InputTag("hltEcalRecHit","EcalRecHitsEB"));
    desc.add<edm::InputTag>("eeRecHitsColl", edm::InputTag("hltEcalRecHit","EcalRecHitsEE"));
    descriptions.add("caloJetTimingProducer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloJetTimingProducer);
