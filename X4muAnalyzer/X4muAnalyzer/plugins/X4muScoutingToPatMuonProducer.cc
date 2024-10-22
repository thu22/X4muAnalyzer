/*
 * Developed by Yiyang Zhao for Run-3 X-->J/psiJ/psi-->4mu Scouting Analysis
 * 2024-10
 * X4muScoutingToPatMuonProducer unpacks Run3ScoutingMuon formats to pat format
 * and cross link Run3ScoutingMuon - Run3ScoutingTrack
*/

#include <memory>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDPutToken.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

//#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
//#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
//#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

//#include "DataFormats/JetReco/interface/PFJet.h"
//#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
//#include "fastjet/contrib/SoftKiller.hh"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

/* class X4muScoutingToPatMuonProducer : public edm::stream::EDProducer<> {
public:
  explicit X4muScoutingToPatMuonProducer(edm::ParameterSet const &iConfig);
  ~X4muScoutingToPatMuonProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &iEvent, edm::EventSetup const &setup) override;
  void endStream() override {}

  reco::Track createTrack(Run3ScoutingMuon const& scoutingMuon);
  reco::Track createTrack(Run3ScoutingTrack const& scoutingTrack);
  reco::Vertex createVertex(Run3ScoutingVertex const& scoutingVertex);
  reco::Muon createMuon(Run3ScoutingMuon const& scoutingMuon);

  void buildHitPattern(Run3ScoutingMuon const& scoutingMuon, Run3ScoutingTrack const& scoutingTrack, reco::Track & recoTrack);

  void createMuons(edm::Handle<std::vector<Run3ScoutingMuon>> scoutingmuonHandle,
                          std::unique_ptr<reco::MuonCollection> &scoutingmuons);
  
  void clearVars();

private:
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> input_scoutingmuon_token_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingTrack>> input_scoutingtrack_token_;

  //std::vector<reco::Muon> RecoMuon_;
  std::vector<int> MuonCharge_;
  std::vector<float> MuonPt_;
  std::vector<float> MuonEta_;
  std::vector<float> MuonPhi_;
  //std::vector<int8_t> vertexIndex_;
  std::vector<float> normchi2_;
  std::vector<float> trkPt_;
  std::vector<float> trkEta_;
  std::vector<float> trkPhi_;
};

//
// constructors and destructor
//
X4muScoutingToPatMuonProducer::X4muScoutingToPatMuonProducer( edm::ParameterSet const &iConfig)
    : input_scoutingmuon_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingMuon"))),
    input_scoutingtrack_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingTrack")))
  {
  //register products
  produces<pat::MuonCollection>();
  produces<edm::ValueMap<int>>("charge");
  produces<edm::ValueMap<float>>("pt");
  produces<edm::ValueMap<float>>("eta");
  produces<edm::ValueMap<float>>("phi");
  //produces<edm::ValueMap<std::vector<int>>>("vertexIndex");
  produces<edm::ValueMap<float>>("normchi2");
  produces<edm::ValueMap<float>>("trkPt");
  produces<edm::ValueMap<float>>("trkEta");
  produces<edm::ValueMap<float>>("trkPhi");
}

X4muScoutingToPatMuonProducer::~X4muScoutingToPatMuonProducer() = default;

reco::Muon X4muScoutingToPatMuonProducer::createMuon(Run3ScoutingMuon const& scoutingMuon) {
  auto m = 0.1056583755;
  auto q = scoutingMuon.charge();

  float px = scoutingMuon.pt() * cos(scoutingMuon.phi());
  float py = scoutingMuon.pt() * sin(scoutingMuon.phi());
  float pz = scoutingMuon.pt() * sinh(scoutingMuon.eta());
  float p = scoutingMuon.pt() * cosh(scoutingMuon.eta());
  float energy = std::sqrt(p * p + m * m);
  reco::Particle::LorentzVector p4(px, py, pz, energy);

  static const reco::Muon dummy;
  auto recomuon = reco::Muon(q, p4, createTrack(scoutingMuon).vertex());
  auto patmuon = pat::Muon(recomuon);

  //RecoMuon_.push_back(recomuon);
  MuonCharge_.push_back(q);
  MuonPt_.push_back(scoutingMuon.pt());
  MuonEta_.push_back(scoutingMuon.eta());
  MuonPhi_.push_back(scoutingMuon.phi());
  //vertexIndex_.push_back(scoutingMuon.vtxIndx());
  normchi2_.push_back(scoutingMuon.normalizedChi2());
  trkPt_.push_back(scoutingMuon.trk_pt());
  trkEta_.push_back(scoutingMuon.trk_eta());
  trkPhi_.push_back(scoutingMuon.trk_phi());

  //std::cout << "Muon px: " << px << std::endl;

  return patmuon;
}

// ------------ method called to produce the data  ------------
void X4muScoutingToPatMuonProducer::produce(edm::Event &iEvent, edm::EventSetup const &setup) {
  using namespace edm;

  Handle<std::vector<Run3ScoutingMuon>> scoutingmuonHandle;
  iEvent.getByToken(input_scoutingmuon_token_, scoutingmuonHandle);

  if(!scoutingmuonHandle.isValid()) {
    return;
  }
  //else std::cout << "ScoutingMuon is valid" << std::endl;

  auto pfcands = std::make_unique<pat::MuonCollection>();

  createMuons(scoutingmuonHandle, pfcands);

  edm::OrphanHandle<pat::MuonCollection> oh = iEvent.put(std::move(pfcands));

  std::unique_ptr<edm::ValueMap<int>> charge_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_charge(*charge_VM);
  filler_charge.insert(oh, MuonCharge_.begin(), MuonCharge_.end());
  filler_charge.fill();
  iEvent.put(std::move(charge_VM), "charge");

  std::unique_ptr<edm::ValueMap<float>> pt_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_pt(*pt_VM);
  filler_pt.insert(oh, MuonPt_.begin(), MuonPt_.end());
  filler_pt.fill();
  iEvent.put(std::move(pt_VM), "pt");

  std::unique_ptr<edm::ValueMap<float>> eta_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_eta(*eta_VM);
  filler_eta.insert(oh, MuonEta_.begin(), MuonEta_.end());
  filler_eta.fill();
  iEvent.put(std::move(eta_VM), "eta");

  std::unique_ptr<edm::ValueMap<float>> phi_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_phi(*phi_VM);
  filler_phi.insert(oh, MuonPhi_.begin(), MuonPhi_.end());
  filler_phi.fill();
  iEvent.put(std::move(phi_VM), "phi"); */

  /* std::unique_ptr<edm::ValueMap<int>> vertexIndex_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_vertexIndex(*vertexIndex_VM);
  filler_vertexIndex.insert(oh, vertexIndex_.begin(), vertexIndex_.end());
  filler_vertexIndex.fill();
  iEvent.put(std::move(vertexIndex_VM), "vertexIndex"); */

/*   std::unique_ptr<edm::ValueMap<float>> normchi2_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_normchi2(*normchi2_VM);
  filler_normchi2.insert(oh, normchi2_.begin(), normchi2_.end());
  filler_normchi2.fill();
  iEvent.put(std::move(normchi2_VM), "normchi2");

  std::unique_ptr<edm::ValueMap<float>> trkPt_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkPt(*trkPt_VM);
  filler_trkPt.insert(oh, trkPt_.begin(), trkPt_.end());
  filler_trkPt.fill();
  iEvent.put(std::move(trkPt_VM), "trkPt");

  std::unique_ptr<edm::ValueMap<float>> trkEta_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkEta(*trkEta_VM);
  filler_trkEta.insert(oh, trkEta_.begin(), trkEta_.end());
  filler_trkEta.fill();
  iEvent.put(std::move(trkEta_VM), "trkEta");

  std::unique_ptr<edm::ValueMap<float>> trkPhi_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkPhi(*trkPhi_VM);
  filler_trkPhi.insert(oh, trkPhi_.begin(), trkPhi_.end());
  filler_trkPhi.fill();
  iEvent.put(std::move(trkPhi_VM), "trkPhi");

  clearVars();
}

void X4muScoutingToPatMuonProducer::clearVars() {
  //RecoMuon_.clear();
  MuonCharge_.clear();
  MuonPt_.clear();
  MuonEta_.clear();
  MuonPhi_.clear();
  //vertexIndex_.clear();
  normchi2_.clear();
  trkPt_.clear();
  trkEta_.clear();
  trkPhi_.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void X4muScoutingToPatMuonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingMuon", edm::InputTag("hltScoutingMuonPackerVtx"));
  desc.add<edm::InputTag>("scoutingTrack", edm::InputTag("hltScoutingTrackPacker"));
  descriptions.add("X4muScoutingToPatMuonProducer", desc);
}


reco::Track X4muScoutingToPatMuonProducer::createTrack(Run3ScoutingTrack const& scoutingTrack){
    float chi2 = scoutingTrack.tk_chi2();
    float ndof = scoutingTrack.tk_ndof();
    reco::TrackBase::Point referencePoint(scoutingTrack.tk_vx(), scoutingTrack.tk_vy(), scoutingTrack.tk_vz());
    
    float px = scoutingTrack.tk_pt() * cos(scoutingTrack.tk_phi());
    float py = scoutingTrack.tk_pt() * sin(scoutingTrack.tk_phi());
    float pz = scoutingTrack.tk_pt() * sinh(scoutingTrack.tk_eta());
    reco::TrackBase::Vector momentum(px, py, pz);
    
    int charge = scoutingTrack.tk_charge();
    
    std::vector<float> cov_vec(15); // 5*(5+1)/2 = 15
    cov_vec[0] = scoutingTrack.tk_qoverp_Error() * scoutingTrack.tk_qoverp_Error(); // cov(0, 0)
    cov_vec[1] = scoutingTrack.tk_qoverp_lambda_cov(); // cov(0, 1)
    cov_vec[3] = scoutingTrack.tk_qoverp_phi_cov(); // cov(0, 2)
    cov_vec[6] = scoutingTrack.tk_qoverp_dxy_cov(); // cov(0, 3)
    cov_vec[10] = scoutingTrack.tk_qoverp_dsz_cov(); // cov(0, 4)
    cov_vec[2] = scoutingTrack.tk_lambda_Error() * scoutingTrack.tk_lambda_Error(); // cov(1, 1)
    cov_vec[4] = scoutingTrack.tk_lambda_phi_cov(); // cov(1, 2)
    cov_vec[7] = scoutingTrack.tk_lambda_dxy_cov(); // cov(1, 3)
    cov_vec[11] = scoutingTrack.tk_lambda_dsz_cov(); // cov(1, 4)
    cov_vec[5] = scoutingTrack.tk_phi_Error() * scoutingTrack.tk_phi_Error(); // cov(2, 2) 
    cov_vec[8] = scoutingTrack.tk_phi_dxy_cov(); // cov(2, 3)
    cov_vec[12] = scoutingTrack.tk_phi_dsz_cov(); // cov(2, 4)
    cov_vec[9] = scoutingTrack.tk_dxy_Error() * scoutingTrack.tk_dxy_Error(); // cov(3, 3)
    cov_vec[13] = scoutingTrack.tk_dxy_dsz_cov(); // cov(3, 4)
    cov_vec[14] = scoutingTrack.tk_dsz_Error() * scoutingTrack.tk_dsz_Error(); // cov(4, 4)
    reco::TrackBase::CovarianceMatrix cov(cov_vec.begin(), cov_vec.end());

    reco::TrackBase::TrackAlgorithm algo(reco::TrackBase::undefAlgorithm); // undefined
    reco::TrackBase::TrackQuality quality(reco::TrackBase::confirmed); // confirmed
    
    // the rests are default: t0 = 0, beta = 0, covt0t0 = -1, covbetabeta = -1 
  
    reco::Track recoTrack(chi2, ndof, referencePoint, momentum, charge, cov, algo, quality);

    return recoTrack;
}

reco::Track X4muScoutingToPatMuonProducer::createTrack(Run3ScoutingMuon const& scoutingMuon){
    float chi2 = scoutingMuon.trk_chi2();
    float ndof = scoutingMuon.trk_ndof();
    reco::TrackBase::Point referencePoint(scoutingMuon.trk_vx(), scoutingMuon.trk_vy(), scoutingMuon.trk_vz());
    
    float px = scoutingMuon.trk_pt() * cos(scoutingMuon.trk_phi());
    float py = scoutingMuon.trk_pt() * sin(scoutingMuon.trk_phi());
    float pz = scoutingMuon.trk_pt() * sinh(scoutingMuon.trk_eta());
    reco::TrackBase::Vector momentum(px, py, pz);
    
    int charge = scoutingMuon.charge();
    
    std::vector<float> cov_vec(15); // 5*(5+1)/2 = 15
    cov_vec[0] = scoutingMuon.trk_qoverpError() * scoutingMuon.trk_qoverpError(); // cov(0, 0)
    cov_vec[1] = scoutingMuon.trk_qoverp_lambda_cov(); // cov(0, 1)
    cov_vec[3] = scoutingMuon.trk_qoverp_phi_cov(); // cov(0, 2)
    cov_vec[6] = scoutingMuon.trk_qoverp_dxy_cov(); // cov(0, 3)
    cov_vec[10] = scoutingMuon.trk_qoverp_dsz_cov(); // cov(0, 4)
    cov_vec[2] = scoutingMuon.trk_lambdaError() * scoutingMuon.trk_lambdaError(); // cov(1, 1)
    cov_vec[4] = scoutingMuon.trk_lambda_phi_cov(); // cov(1, 2)
    cov_vec[7] = scoutingMuon.trk_lambda_dxy_cov(); // cov(1, 3)
    cov_vec[11] = scoutingMuon.trk_lambda_dsz_cov(); // cov(1, 4)
    cov_vec[5] = scoutingMuon.trk_phiError() * scoutingMuon.trk_phiError(); // cov(2, 2) 
    cov_vec[8] = scoutingMuon.trk_phi_dxy_cov(); // cov(2, 3)
    cov_vec[12] = scoutingMuon.trk_phi_dsz_cov(); // cov(2, 4)
    cov_vec[9] = scoutingMuon.trk_dxyError() * scoutingMuon.trk_dxyError(); // cov(3, 3)
    cov_vec[13] = scoutingMuon.trk_dxy_dsz_cov(); // cov(3, 4)
    cov_vec[14] = scoutingMuon.trk_dszError() * scoutingMuon.trk_dszError(); // cov(4, 4)
    reco::TrackBase::CovarianceMatrix cov(cov_vec.begin(), cov_vec.end());

    reco::TrackBase::TrackAlgorithm algo(reco::TrackBase::undefAlgorithm); // undefined
    reco::TrackBase::TrackQuality quality(reco::TrackBase::confirmed); // confirmed
    
    // the rests are default: t0 = 0, beta = 0, covt0t0 = -1, covbetabeta = -1 
  
    reco::Track recoTrack(chi2, ndof, referencePoint, momentum, charge, cov, algo, quality);

    return recoTrack;
}

void X4muScoutingToPatMuonProducer::createMuons(
    edm::Handle<std::vector<Run3ScoutingMuon>> scoutingmuonHandle,
    std::unique_ptr<pat::MuonCollection> &scoutingmuons) {
  for (unsigned int icand = 0; icand < scoutingmuonHandle->size(); ++icand) {
    auto &scoutingmuon = (*scoutingmuonHandle)[icand];

    auto pfcand = createMuon(scoutingmuon);
    if (pfcand.energy() != 0)
      scoutingmuons->push_back(pfcand);
  }
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(X4muScoutingToPatMuonProducer); */