/*
 * Developed by Yiyang Zhao for Run-3 X-->J/psiEta_c-->2mu 6pai Analysis
 * 2024-11
 * X4muPatSecondaryVertexProducer loop reco::Muon + reco::Track
 * and fit secondary vertex for dimuon + six pais
 */

#include <memory>

#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

// #include "DataFormats/JetReco/interface/PFJet.h"
// #include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/Common/interface/View.h>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
// #include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/DeepCopyPointer.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"
#include "DataFormats/CLHEP/interface/Migration.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace muon;
using namespace trigger;

struct ParticleInfo
{
    // 粒子电荷
    int Pai1Charge;
    int Pai2Charge;
    int Pai3Charge;
    int Pai4Charge;
    int Pai5Charge;
    int Pai6Charge;
    int Muon1Charge;
    int Muon2Charge;

    // 粒子的四动量 (pT, Eta, Phi, Mass)
    TLorentzVector Pai1P4, Pai2P4, Pai3P4, Pai4P4, Pai5P4, Pai6P4;
    TLorentzVector Mu1P4, Mu2P4;
    TLorentzVector EtaP4, JpsiP4;
    TLorentzVector EtaNoMassP4, JpsiNoMassP4;
    TLorentzVector XP4;

    // Muon相关信息
    bool Muon1isSoft;
    bool Muon2isSoft;
    bool Muon1isTight;
    bool Muon2isTight;
    bool Muon1isMedium;
    bool Muon2isMedium;
    bool Muon1isLoose;
    bool Muon2isLoose;
    bool Muon1isTracker;
    bool Muon2isTracker;
    bool Muon1isHighPt;
    bool Muon2isHighPt;
    bool Muon1isGlobal;
    bool Muon2isGlobal;

    int Pai1Hit;
    int Pai2Hit;
    int Pai3Hit;
    int Pai4Hit;
    int Pai5Hit;
    int Pai6Hit;
    int Pai1fromPV;
    int Pai2fromPV;
    int Pai3fromPV;
    int Pai4fromPV;
    int Pai5fromPV;
    int Pai6fromPV;
    double Pai1Trknormchi2;
    double Pai2Trknormchi2;
    double Pai3Trknormchi2;
    double Pai4Trknormchi2;
    double Pai5Trknormchi2;
    double Pai6Trknormchi2;

    // Jpsi拟合信息
    double Etanormchi2;
    double Jpsinormchi2;
    double EtaNoMassnormchi2;
    double JpsiNoMassnormchi2;
    float EtaNoMassMassE;
    float JpsiNoMassMassE;
    float EtaMassDiff;

    // X粒子信息
    double Xnormchi2;

    double minDR;

    // 构造函数，初始化所有变量
    ParticleInfo()
    {
        // 初始化电荷
        Pai1Charge = Pai2Charge = Pai3Charge = Pai4Charge = Pai5Charge = Pai6Charge = 0;
        Muon1Charge = Muon2Charge = 0;

        // 初始化四动量（设为零四动量）
        Pai1P4.SetXYZM(0, 0, 0, 0);
        Pai2P4.SetXYZM(0, 0, 0, 0);
        Pai3P4.SetXYZM(0, 0, 0, 0);
        Pai4P4.SetXYZM(0, 0, 0, 0);
        Pai5P4.SetXYZM(0, 0, 0, 0);
        Pai6P4.SetXYZM(0, 0, 0, 0);
        Mu1P4.SetXYZM(0, 0, 0, 0);
        Mu2P4.SetXYZM(0, 0, 0, 0);
        EtaP4.SetXYZM(0, 0, 0, 0);
        JpsiP4.SetXYZM(0, 0, 0, 0);
        EtaNoMassP4.SetXYZM(0, 0, 0, 0);
        JpsiNoMassP4.SetXYZM(0, 0, 0, 0);
        XP4.SetXYZM(0, 0, 0, 0);

        // 初始化Muon信息
        Muon1isSoft = Muon2isSoft = false;
        Muon1isTight = Muon2isTight = false;
        Muon1isMedium = Muon2isMedium = false;
        Muon1isLoose = Muon2isLoose = false;
        Muon1isTracker = Muon2isTracker = false;
        Muon1isHighPt = Muon2isHighPt = false;
        Muon1isGlobal = Muon2isGlobal = false;

        Pai1Hit = Pai2Hit = Pai3Hit = Pai4Hit = Pai5Hit = Pai6Hit = 0;
        Pai1fromPV = Pai2fromPV = Pai3fromPV = Pai4fromPV = Pai5fromPV = Pai6fromPV = 0;
        Pai1Trknormchi2 = Pai2Trknormchi2 = Pai3Trknormchi2 = Pai4Trknormchi2 = Pai5Trknormchi2 = Pai6Trknormchi2 = 0.0;

        // 初始化Jpsi拟合信息
        Etanormchi2 = Jpsinormchi2 = 0.0;
        EtaNoMassnormchi2 = JpsiNoMassnormchi2 = 0.0;
        EtaNoMassMassE = JpsiNoMassMassE = 0.0;

        Xnormchi2 = 0.0;

        minDR = -1;
    }
};

class X4muPatSecondaryVertexProducer : public edm::stream::EDProducer<>
{
public:
    explicit X4muPatSecondaryVertexProducer(edm::ParameterSet const &iConfig);
    ~X4muPatSecondaryVertexProducer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
    void beginStream(edm::StreamID) override {}
    void produce(edm::Event &iEvent, edm::EventSetup const &iSetup) override;

    void Loop(std::vector<int> selectedPai1, std::vector<int> selectedPai2, std::vector<int> selectedPai3, std::vector<int> selectedPai4, std::vector<int> selectedPai5, std::vector<int> selectedPai6,
              std::vector<pat::PackedCandidate> goodTracks,
              std::vector<pat::Muon> dimuonsSort1, std::vector<pat::Muon> dimuonsSort2,
              double myEtamass, double myEtamasserr, double myJpsimass, double myJpsimasserr, double myExmass, double JvPorbcut,
              double myMumass, double myMumasserr, double myPaimass, double myPaimasserr, double myMesonWindow,
              double MassMinCut, const reco::Vertex *pv, const MagneticField &bFieldHandle,
              bool doIso, double MuIso, double PaiIso, unsigned int max_loop);
    void endStream() override {}

    void clearVars();

private:
    UInt_t getTriggerBits(const edm::Event &);
    const edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tracks_;
    const edm::EDGetTokenT<std::vector<pat::Muon>> muons_;
    const edm::EDGetTokenT<std::vector<reco::Vertex>> vtxToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerresults_;
    std::vector<std::string> FilterNames_;
    const double input_MesonPaiMass_c;
    const double input_MesonPaiMassErr_c;
    const double input_MesonMuMass_c;
    const double input_MesonMuMassErr_c;
    const double input_MesonMassWindow_c;
    const double input_ExMesonMass_c;
    const bool doIso_c;
    const double MuIso_c;
    const double PaiIso_c;
    const bool doTrigger_c;
    const double vProb_c;
    const int selectionType_c;
    const unsigned int maxLoop_c;
    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_; // 声明 magneticFieldToken_

    TTree *X4muTree;
    std::string file_name;

    UInt_t run;
    ULong64_t event;
    UInt_t lumiblock;

    UInt_t trigger;

    vector<ParticleInfo> possibleInfos;

    long npairs = 0;

    int Pai1Charge_[36];
    int Pai2Charge_[36];
    int Pai3Charge_[36];
    int Pai4Charge_[36];
    int Pai5Charge_[36];
    int Pai6Charge_[36];
    int Mu1Charge_[36];
    int Mu2Charge_[36];
    float Pai1Pt_[36];
    float Pai2Pt_[36];
    float Pai3Pt_[36];
    float Pai4Pt_[36];
    float Pai5Pt_[36];
    float Pai6Pt_[36];
    float Mu1Pt_[36];
    float Mu2Pt_[36];
    float Pai1Eta_[36];
    float Pai2Eta_[36];
    float Pai3Eta_[36];
    float Pai4Eta_[36];
    float Pai5Eta_[36];
    float Pai6Eta_[36];
    float Mu1Eta_[36];
    float Mu2Eta_[36];
    float Pai1Phi_[36];
    float Pai2Phi_[36];
    float Pai3Phi_[36];
    float Pai4Phi_[36];
    float Pai5Phi_[36];
    float Pai6Phi_[36];
    float Mu1Phi_[36];
    float Mu2Phi_[36];
    float Pai1Mass_[36];
    float Pai2Mass_[36];
    float Pai3Mass_[36];
    float Pai4Mass_[36];
    float Pai5Mass_[36];
    float Pai6Mass_[36];
    float Mu1Mass_[36];
    float Mu2Mass_[36];
    int Pai1Hit_[36];
    int Pai2Hit_[36];
    int Pai3Hit_[36];
    int Pai4Hit_[36];
    int Pai5Hit_[36];
    int Pai6Hit_[36];
    int Pai1fromPV_[36];
    int Pai2fromPV_[36];
    int Pai3fromPV_[36];
    int Pai4fromPV_[36];
    int Pai5fromPV_[36];
    int Pai6fromPV_[36];
    double Pai1Trknormchi2_[36];
    double Pai2Trknormchi2_[36];
    double Pai3Trknormchi2_[36];
    double Pai4Trknormchi2_[36];
    double Pai5Trknormchi2_[36];
    double Pai6Trknormchi2_[36];
    bool Muon1isSoft_[36];
    bool Muon2isSoft_[36];
    bool Muon1isTight_[36];
    bool Muon2isTight_[36];
    bool Muon1isMedium_[36];
    bool Muon2isMedium_[36];
    bool Muon1isLoose_[36];
    bool Muon2isLoose_[36];
    bool Muon1isTracker_[36];
    bool Muon2isTracker_[36];
    bool Muon1isHighPt_[36];
    bool Muon2isHighPt_[36];
    bool Muon1isGlobal_[36];
    bool Muon2isGlobal_[36];
    double J1normchi2_[36];
    double J2normchi2_[36];
    double J1NoMassnormchi2_[36];
    double J2NoMassnormchi2_[36];
    float J1Pt_[36];
    float J2Pt_[36];
    float J1NoMassPt_[36];
    float J2NoMassPt_[36];
    float J1Eta_[36];
    float J2Eta_[36];
    float J1NoMassEta_[36];
    float J2NoMassEta_[36];
    float J1Phi_[36];
    float J2Phi_[36];
    float J1NoMassPhi_[36];
    float J2NoMassPhi_[36];
    float J1Mass_[36];
    float J2Mass_[36];
    float J1NoMassMass_[36];
    float J1NoMassMassE_[36];
    float J2NoMassMass_[36];
    float J2NoMassMassE_[36];
    float EtaMassDiff_[36];
    float XpT_[36];
    float Xeta_[36];
    float Xphi_[36];
    float Xmass_[36];
    double Xnormchi2_[36];
    double minDR_[36];
    int type_[36];
};

//
// constructors and destructor
//
X4muPatSecondaryVertexProducer::X4muPatSecondaryVertexProducer(edm::ParameterSet const &iConfig)
    : tracks_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      muons_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("patMuon"))),
      vtxToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
      triggerresults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
      FilterNames_(iConfig.getParameter<std::vector<std::string>>("FilterNames")),
      input_MesonPaiMass_c(iConfig.getParameter<double>("MesonPaiMass")),
      input_MesonPaiMassErr_c(iConfig.getParameter<double>("MesonPaiMassErr")),
      input_MesonMuMass_c(iConfig.getParameter<double>("MesonMuMass")),
      input_MesonMuMassErr_c(iConfig.getParameter<double>("MesonMuMassErr")),
      input_MesonMassWindow_c(iConfig.getParameter<double>("MesonMassWindow")),
      input_ExMesonMass_c(iConfig.getParameter<double>("ExMesonMass")),
      doIso_c(iConfig.getParameter<bool>("doIso")),
      MuIso_c(iConfig.getParameter<double>("MuIso")),
      PaiIso_c(iConfig.getParameter<double>("PaiIso")),
      doTrigger_c(iConfig.getParameter<bool>("doTrigger")),
      vProb_c(iConfig.getParameter<double>("vProb")),
      selectionType_c(iConfig.getParameter<int>("selectionType")),
      maxLoop_c(iConfig.getParameter<unsigned int>("maxLoop")),
      magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
    edm::Service<TFileService> fs;
    X4muTree = fs->make<TTree>("X4muTree", "Tree of X4muTree");

    X4muTree->Branch("run", &run, "run/i");
    X4muTree->Branch("event", &event, "event/i");
    X4muTree->Branch("lumiblock", &lumiblock, "lumiblock/i");
    X4muTree->Branch("npairs", &npairs, "npairs/l");
    X4muTree->Branch("trigger", &trigger, "trigger/i");
    X4muTree->Branch("Pai1Charge", Pai1Charge_, "Pai1Charge[36]/I");
    X4muTree->Branch("Pai2Charge", Pai2Charge_, "Pai2Charge[36]/I");
    X4muTree->Branch("Pai3Charge", Pai3Charge_, "Pai3Charge[36]/I");
    X4muTree->Branch("Pai4Charge", Pai4Charge_, "Pai4Charge[36]/I");
    X4muTree->Branch("Pai5Charge", Pai5Charge_, "Pai5Charge[36]/I");
    X4muTree->Branch("Pai6Charge", Pai6Charge_, "Pai6Charge[36]/I");
    X4muTree->Branch("Muon1Charge", Mu1Charge_, "Muon1Charge[36]/I");
    X4muTree->Branch("Muon2Charge", Mu2Charge_, "Muon2Charge[36]/I");
    X4muTree->Branch("Pai1Pt", Pai1Pt_, "Pai1Pt[36]/F");
    X4muTree->Branch("Pai2Pt", Pai2Pt_, "Pai2Pt[36]/F");
    X4muTree->Branch("Pai3Pt", Pai3Pt_, "Pai3Pt[36]/F");
    X4muTree->Branch("Pai4Pt", Pai4Pt_, "Pai4Pt[36]/F");
    X4muTree->Branch("Pai5Pt", Pai5Pt_, "Pai5Pt[36]/F");
    X4muTree->Branch("Pai6Pt", Pai6Pt_, "Pai6Pt[36]/F");
    X4muTree->Branch("Muon1Pt", Mu1Pt_, "Muon1Pt[36]/F");
    X4muTree->Branch("Muon2Pt", Mu2Pt_, "Muon2Pt[36]/F");
    X4muTree->Branch("Pai1Eta", Pai1Eta_, "Pai1Eta[36]/F");
    X4muTree->Branch("Pai2Eta", Pai2Eta_, "Pai2Eta[36]/F");
    X4muTree->Branch("Pai3Eta", Pai3Eta_, "Pai3Eta[36]/F");
    X4muTree->Branch("Pai4Eta", Pai4Eta_, "Pai4Eta[36]/F");
    X4muTree->Branch("Pai5Eta", Pai5Eta_, "Pai5Eta[36]/F");
    X4muTree->Branch("Pai6Eta", Pai6Eta_, "Pai6Eta[36]/F");
    X4muTree->Branch("Muon1Eta", Mu1Eta_, "Muon1Eta[36]/F");
    X4muTree->Branch("Muon2Eta", Mu2Eta_, "Muon2Eta[36]/F");
    X4muTree->Branch("Pai1Phi", Pai1Phi_, "Pai1Phi[36]/F");
    X4muTree->Branch("Pai2Phi", Pai2Phi_, "Pai2Phi[36]/F");
    X4muTree->Branch("Pai3Phi", Pai3Phi_, "Pai3Phi[36]/F");
    X4muTree->Branch("Pai4Phi", Pai4Phi_, "Pai4Phi[36]/F");
    X4muTree->Branch("Pai5Phi", Pai5Phi_, "Pai5Phi[36]/F");
    X4muTree->Branch("Pai6Phi", Pai6Phi_, "Pai6Phi[36]/F");
    X4muTree->Branch("Muon1Phi", Mu1Phi_, "Muon1Phi[36]/F");
    X4muTree->Branch("Muon2Phi", Mu2Phi_, "Muon2Phi[36]/F");
    X4muTree->Branch("Pai1Mass", Pai1Mass_, "Pai1Mass[36]/F");
    X4muTree->Branch("Pai2Mass", Pai2Mass_, "Pai2Mass[36]/F");
    X4muTree->Branch("Pai3Mass", Pai3Mass_, "Pai3Mass[36]/F");
    X4muTree->Branch("Pai4Mass", Pai4Mass_, "Pai4Mass[36]/F");
    X4muTree->Branch("Pai5Mass", Pai5Mass_, "Pai5Mass[36]/F");
    X4muTree->Branch("Pai6Mass", Pai6Mass_, "Pai6Mass[36]/F");
    X4muTree->Branch("Pai1Hit", Pai1Hit_, "Pai1Hit[36]/I");
    X4muTree->Branch("Pai2Hit", Pai2Hit_, "Pai2Hit[36]/I");
    X4muTree->Branch("Pai3Hit", Pai3Hit_, "Pai3Hit[36]/I");
    X4muTree->Branch("Pai4Hit", Pai4Hit_, "Pai4Hit[36]/I");
    X4muTree->Branch("Pai5Hit", Pai5Hit_, "Pai5Hit[36]/I");
    X4muTree->Branch("Pai6Hit", Pai6Hit_, "Pai6Hit[36]/I");
    X4muTree->Branch("Pai1fromPV", Pai1fromPV_, "Pai1fromPV[36]/I");
    X4muTree->Branch("Pai2fromPV", Pai2fromPV_, "Pai2fromPV[36]/I");
    X4muTree->Branch("Pai3fromPV", Pai3fromPV_, "Pai3fromPV[36]/I");
    X4muTree->Branch("Pai4fromPV", Pai4fromPV_, "Pai4fromPV[36]/I");
    X4muTree->Branch("Pai5fromPV", Pai5fromPV_, "Pai5fromPV[36]/I");
    X4muTree->Branch("Pai6fromPV", Pai6fromPV_, "Pai6fromPV[36]/I");
    X4muTree->Branch("Pai1Trknormchi2", Pai1Trknormchi2_, "Pai1Trknormchi2[36]/D");
    X4muTree->Branch("Pai2Trknormchi2", Pai2Trknormchi2_, "Pai2Trknormchi2[36]/D");
    X4muTree->Branch("Pai3Trknormchi2", Pai3Trknormchi2_, "Pai3Trknormchi2[36]/D");
    X4muTree->Branch("Pai4Trknormchi2", Pai4Trknormchi2_, "Pai4Trknormchi2[36]/D");
    X4muTree->Branch("Pai5Trknormchi2", Pai5Trknormchi2_, "Pai5Trknormchi2[36]/D");
    X4muTree->Branch("Pai6Trknormchi2", Pai6Trknormchi2_, "Pai6Trknormchi2[36]/D");
    X4muTree->Branch("Muon1Mass", Mu1Mass_, "Muon1Mass[36]/F");
    X4muTree->Branch("Muon2Mass", Mu2Mass_, "Muon2Mass[36]/F");
    X4muTree->Branch("Muon1isSoft", Muon1isSoft_, "Muon1isSoft[36]/O");
    X4muTree->Branch("Muon2isSoft", Muon2isSoft_, "Muon2isSoft[36]/O");
    X4muTree->Branch("Muon1isTight", Muon1isTight_, "Muon1isTight[36]/O");
    X4muTree->Branch("Muon2isTight", Muon2isTight_, "Muon2isTight[36]/O");
    X4muTree->Branch("Muon1isMedium", Muon1isMedium_, "Muon1isMedium[36]/O");
    X4muTree->Branch("Muon2isMedium", Muon2isMedium_, "Muon2isMedium[36]/O");
    X4muTree->Branch("Muon1isLoose", Muon1isLoose_, "Muon1isLoose[36]/O");
    X4muTree->Branch("Muon2isLoose", Muon2isLoose_, "Muon2isLoose[36]/O");
    X4muTree->Branch("Muon1isTracker", Muon1isTracker_, "Muon1isTracker[36]/O");
    X4muTree->Branch("Muon2isTracker", Muon2isTracker_, "Muon2isTracker[36]/O");
    X4muTree->Branch("Muon1isHighPt", Muon1isHighPt_, "Muon1isHighPt[36]/O");
    X4muTree->Branch("Muon2isHighPt", Muon2isHighPt_, "Muon2isHighPt[36]/O");
    X4muTree->Branch("J1normchi2", J1normchi2_, "J1normchi2[36]/D");
    X4muTree->Branch("J2normchi2", J2normchi2_, "J2normchi2[36]/D");
    X4muTree->Branch("J1NoMassnormchi2", J1NoMassnormchi2_, "J1NoMassnormchi2[36]/D");
    X4muTree->Branch("J2NoMassnormchi2", J2NoMassnormchi2_, "J2NoMassnormchi2[36]/D");
    X4muTree->Branch("J1Pt", J1Pt_, "J1Pt[36]/F");
    X4muTree->Branch("J2Pt", J2Pt_, "J2Pt[36]/F");
    X4muTree->Branch("J1NoMassPt", J1NoMassPt_, "J1NoMassPt[36]/F");
    X4muTree->Branch("J2NoMassPt", J2NoMassPt_, "J2NoMassPt[36]/F");
    X4muTree->Branch("J1Eta", J1Eta_, "J1Eta[36]/F");
    X4muTree->Branch("J2Eta", J2Eta_, "J2Eta[36]/F");
    X4muTree->Branch("J1NoMassEta", J1NoMassEta_, "J1NoMassEta[36]/F");
    X4muTree->Branch("J2NoMassEta", J2NoMassEta_, "J2NoMassEta[36]/F");
    X4muTree->Branch("J1Phi", J1Phi_, "J1Phi[36]/F");
    X4muTree->Branch("J2Phi", J2Phi_, "J2Phi[36]/F");
    X4muTree->Branch("J1NoMassPhi", J1NoMassPhi_, "J1NoMassPhi[36]/F");
    X4muTree->Branch("J2NoMassPhi", J2NoMassPhi_, "J2NoMassPhi[36]/F");
    X4muTree->Branch("J1Mass", J1Mass_, "J1Mass[36]/F");
    X4muTree->Branch("J2Mass", J2Mass_, "J2Mass[36]/F");
    X4muTree->Branch("J1NoMassMass", J1NoMassMass_, "J1NoMassMass[36]/F");
    X4muTree->Branch("J1NoMassMassE", J1NoMassMassE_, "J1NoMassMassE[36]/F");
    X4muTree->Branch("J2NoMassMass", J2NoMassMass_, "J2NoMassMass[36]/F");
    X4muTree->Branch("J2NoMassMassE", J2NoMassMassE_, "J2NoMassMassE[36]/F");
    X4muTree->Branch("EtaMassDiff", EtaMassDiff_, "EtaMassDiff[36]/F");
    X4muTree->Branch("XpT", XpT_, "XpT[36]/F");
    X4muTree->Branch("Xeta", Xeta_, "Xeta[36]/F");
    X4muTree->Branch("Xphi", Xphi_, "Xphi[36]/F");
    X4muTree->Branch("Xmass", Xmass_, "Xmass[36]/F");
    X4muTree->Branch("Xnormchi2", Xnormchi2_, "Xnormchi2[36]/D");
    X4muTree->Branch("minDR", minDR_, "minDR[36]/D");
    X4muTree->Branch("type", type_, "type[36]/I");
}

X4muPatSecondaryVertexProducer::~X4muPatSecondaryVertexProducer() = default;

// ------------ method called to produce the data  ------------
void X4muPatSecondaryVertexProducer::produce(edm::Event &iEvent, edm::EventSetup const &iSetup)
{
    using namespace edm;
    using namespace std;
    using namespace reco;

    const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_);
    const auto myEtamass = input_MesonPaiMass_c;
    const auto myEtamasserr = input_MesonPaiMassErr_c;
    const auto myJpsimass = input_MesonMuMass_c;
    const auto myJpsimasserr = input_MesonMuMassErr_c;
    const auto myMesonWindow = input_MesonMassWindow_c;
    const auto myExmass = input_ExMesonMass_c;
    const auto JvPorbcut = vProb_c;

    const auto doIso = doIso_c;
    const auto MuIso = MuIso_c;
    const auto PaiIso = PaiIso_c;
    const auto doTrigger = doTrigger_c;
    const auto selectingType = selectionType_c;
    const unsigned int max_loop = maxLoop_c;

    edm::Handle<std::vector<pat::PackedCandidate>> tracks;
    iEvent.getByToken(tracks_, tracks);
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muons_, muons);

    run = iEvent.id().run();
    event = iEvent.id().event();
    lumiblock = iEvent.id().luminosityBlock();
    trigger = getTriggerBits(iEvent);

    double myMumass = 0.1056583755;
    double myMumasserr = myMumass * 1e-6;
    double myPaimass = 0.13957039;
    double myPaimasserr = 0.00000018;
    double ProbExcut = 0.03;
    double MassMinCut = 0.01;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    bool goodVtx = false;
    const reco::Vertex *pv;
    for (const reco::Vertex &vtx : *vertices)
    {
        if (vtx.isFake() || !vtx.isValid())
            continue;
        goodVtx = true;
        pv = &vtx;
        break;
    }
    if (!goodVtx)
        return;

    if (doTrigger && !trigger)
        return;

    if (!tracks.isValid() || !muons.isValid())
    {
        std::cout << "reco::Track or pat::Muon collection is not valid" << std::endl;
        return;
    }

    if (tracks->size() < 6 || muons->size() < 2)
    {
        std::cout << "reco::Track or pat::Muon size is too small: " << tracks->size() << " " << muons->size() << std::endl;
        return;
    }

    std::vector<pat::PackedCandidate> goodTracks;
    for (std::vector<pat::PackedCandidate>::const_iterator iTrack = tracks->begin(); iTrack != tracks->end(); ++iTrack)
    {
        if (abs(iTrack->pdgId()) != 211)  //Check
            continue;
        if (iTrack->charge() == 0 || iTrack->pt() < 0.5 || abs(iTrack->eta()) > 2.5)
            continue;
        if (!(iTrack->hasTrackDetails()) || iTrack->numberOfHits() < 5 || iTrack->bestTrack()->normalizedChi2() > 8 || !(iTrack->bestTrack()->quality(reco::Track::highPurity)))
            continue;
        if (iTrack->fromPV() < 1)
            continue;
        //Check not muon
        goodTracks.push_back(*iTrack);
    }

    if (goodTracks.size() < 6)
    {
        // std::cout << "goodTracks size too small: " << goodTracks.size() << std::endl;
        return;
    }

    // Step 1: Dimuon vertex fitting
    std::vector<pat::Muon> dimuons1;
    std::vector<pat::Muon> dimuons2;
    long loop1 = 0;
    for (auto muon1 = muons->begin(); muon1 != muons->end(); ++muon1)
    {
        for (auto muon2 = muon1 + 1; muon2 != muons->end(); ++muon2)
        {
            loop1++;
            if (loop1 > max_loop)
            {
                // std::cout << "loop1 > max_loop" << std::endl;
                goto step1;
            }
            // Ensure opposite charge
            if (muon1->charge() + muon2->charge() != 0)
                continue;
            // Check delta R between muons
            if (doIso && deltaR(*muon1, *muon2) < MuIso)
                continue;

            // Perform dimuon vertex fit without mass constraint
            KinematicParticleFactoryFromTransientTrack pmumuFactory;
            TrackRef muTrack1 = muon1->track();
            TrackRef muTrack2 = muon2->track();

            if (muTrack1.isNull() || muTrack2.isNull())
                continue;

            TransientTrack muonTT1(muTrack1, &(bFieldHandle));
            TransientTrack muonTT2(muTrack2, &(bFieldHandle));

            // initial chi2 and ndf before kinematic fits.
            float chi = 0.;
            float ndf = 0.;
            ParticleMass mu_mass = myMumass;
            float mu_sigma = myMumasserr;

            vector<RefCountedKinematicParticle> dimuonParticles;
            dimuonParticles.push_back(pmumuFactory.particle(muonTT1, mu_mass, chi, ndf, mu_sigma));
            dimuonParticles.push_back(pmumuFactory.particle(muonTT2, mu_mass, chi, ndf, mu_sigma));

            KinematicParticleVertexFitter dimuonFitter;
            RefCountedKinematicTree dimuonVertexFitTree;
            try
            {
                dimuonVertexFitTree = dimuonFitter.fit(dimuonParticles);
            }
            catch (const std::exception &e)
            {
                continue;
            }

            if (!dimuonVertexFitTree->isValid())
                continue;

            dimuonVertexFitTree->movePointerToTheTop();
            RefCountedKinematicParticle dimuon = dimuonVertexFitTree->currentParticle();
            RefCountedKinematicVertex dimuonVertex = dimuonVertexFitTree->currentDecayVertex();

            if (!dimuon->currentState().isValid() || !dimuonVertex->vertexIsValid())
                continue;

            dimuonVertexFitTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitMuon1 = dimuonVertexFitTree->currentParticle();
            dimuonVertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitMuon2 = dimuonVertexFitTree->currentParticle();

            if (!fitMuon1->currentState().isValid() || !fitMuon2->currentState().isValid())
                continue;

            if ((dimuonVertex->chiSquared()) < 0 || (dimuonVertex->degreesOfFreedom()) <= 0 || (dimuonVertex->chiSquared()) > 9999.9)
                continue;
            double vProb1 = ChiSquaredProbability((double)(dimuonVertex->chiSquared()), (double)(dimuonVertex->degreesOfFreedom()));
            if (vProb1 < JvPorbcut)
                continue;
            if (fitMuon1->currentState().mass() <= MassMinCut || fitMuon2->currentState().mass() <= MassMinCut)
                continue;
            if (dimuon->currentState().kinematicParametersError().matrix()(6, 6) < 0)
                continue;
            float Jpsimasserr = sqrt(dimuon->currentState().kinematicParametersError().matrix()(6, 6));
            //Store dimuon p4
            if (abs(dimuon->currentState().mass() - myExmass) < 2.0 * Jpsimasserr)
                continue;
            if (abs(dimuon->currentState().mass() - myJpsimass) > 3.0 * Jpsimasserr || abs(dimuon->currentState().mass() - myJpsimass) > myMesonWindow)
                continue;

            // Mass Constraint Fit for Jpsi
            KinematicParticleVertexFitter kpvFitter;
            KinematicParticleFitter csFitter;
            ParticleMass JpsiM = myJpsimass;
            float jpsi_m_sigma = myJpsimasserr;
            KinematicConstraint *jpsi_c = new MassKinematicConstraint(JpsiM, jpsi_m_sigma);
            vector<RefCountedKinematicParticle> dimuonRef;
            RefCountedKinematicParticle Jpsi_part;

            dimuonRef.push_back(pmumuFactory.particle(muonTT1, mu_mass, chi, ndf, mu_sigma));
            dimuonRef.push_back(pmumuFactory.particle(muonTT2, mu_mass, chi, ndf, mu_sigma));

            RefCountedKinematicTree Jpsi;
            RefCountedKinematicTree JpsiNoMC;
            RefCountedKinematicParticle MyChi1_part;
            try
            {
                JpsiNoMC = kpvFitter.fit(dimuonRef);
            }
            catch (const std::exception &e)
            {
                continue;
            }
            if (JpsiNoMC->isEmpty())
                continue;
            try
            {
                Jpsi = csFitter.fit(jpsi_c, JpsiNoMC);
            }
            catch (const std::exception &e)
            {
                continue;
            }
            if (Jpsi->isEmpty())
                continue;

            Jpsi->movePointerToTheTop();
            RefCountedKinematicVertex myJpsiVertex = Jpsi->currentDecayVertex();
            if (!myJpsiVertex->vertexIsValid())
                continue;
            if ((myJpsiVertex->chiSquared()) < 0 || (myJpsiVertex->degreesOfFreedom()) <= 0 || (myJpsiVertex->chiSquared()) > 9999.9)
                continue;
            double vProb2 = ChiSquaredProbability((double)(myJpsiVertex->chiSquared()), (double)(myJpsiVertex->degreesOfFreedom()));
            if (vProb2 < JvPorbcut)
                continue;

            // Store valid dimuons
            dimuons1.push_back(*muon1);
            dimuons2.push_back(*muon2);
        }
    }
step1:

    if (dimuons1.size() < 1)
        return;

    // Sort dimuons1 and dimuons2 according to dimuon's pT
    std::vector<std::pair<pat::Muon, pat::Muon>> dimuons;
    for (int i = 0; i < (int)dimuons1.size(); i++)
    {
        dimuons.push_back(std::make_pair(dimuons1[i], dimuons2[i]));
    }
    std::sort(dimuons.begin(), dimuons.end(), [](const std::pair<pat::Muon, pat::Muon> &a, const std::pair<pat::Muon, pat::Muon> &b)
              {
        float dimuon_pT_a = (a.first.p4() + a.second.p4()).pt();
        float dimuon_pT_b = (b.first.p4() + b.second.p4()).pt();
        return dimuon_pT_a > dimuon_pT_b; });
    std::vector<pat::Muon> dimuonsSort1;
    std::vector<pat::Muon> dimuonsSort2;
    for (int i = 0; i < (int)dimuons.size(); i++)
    {
        dimuonsSort1.push_back(dimuons[i].first);
        dimuonsSort2.push_back(dimuons[i].second);
    }

    // std::cout << "dimuons size: " << dimuons1.size() << std::endl;

    // Step 2: Loop 3 pais
    long loop2 = 0;

    // Loop over combinations of 3 Pai candidates
    // Store 3 pais in 7 vector pair according to there total charge
    std::vector<int> selectedPai1Charge3;
    std::vector<int> selectedPai2Charge3;
    std::vector<int> selectedPai3Charge3;
    std::vector<int> selectedPai1ChargeM3;
    std::vector<int> selectedPai2ChargeM3;
    std::vector<int> selectedPai3ChargeM3;

    for (std::vector<pat::PackedCandidate>::iterator iPai1 = goodTracks.begin(); iPai1 != goodTracks.end(); ++iPai1)
    {
        for (std::vector<pat::PackedCandidate>::iterator iPai2 = iPai1 + 1; iPai2 != goodTracks.end(); ++iPai2)
        {
            for (std::vector<pat::PackedCandidate>::iterator iPai3 = iPai2 + 1; iPai3 != goodTracks.end(); ++iPai3)
            {
                loop2++;
                if (loop2 > max_loop)
                {
                    // std::cout << "loop2 > max_loop" << std::endl;
                    goto step2;
                }
                int Tcharge = iPai1->charge() + iPai2->charge() + iPai3->charge();
                if (abs(Tcharge) != 3)
                    continue;
                reco::Track PaiTrack1 = *(iPai1->bestTrack());
                reco::Track PaiTrack2 = *(iPai2->bestTrack());
                reco::Track PaiTrack3 = *(iPai3->bestTrack());

                TransientTrack paiTT1(PaiTrack1, &(bFieldHandle));
                TransientTrack paiTT2(PaiTrack2, &(bFieldHandle));
                TransientTrack paiTT3(PaiTrack3, &(bFieldHandle));

                if (!paiTT1.isValid() || !paiTT2.isValid() || !paiTT3.isValid())
                    continue;

                // check deltaR between tracks
                if (doIso && (abs(deltaR(*iPai1, *iPai2)) < PaiIso || abs(deltaR(*iPai1, *iPai3)) < PaiIso || abs(deltaR(*iPai2, *iPai3)) < PaiIso))
                    continue;

                ParticleMass pai_mass = myPaimass; // pdg mass
                float pai_sigma = myPaimasserr;

                KinematicParticleFactoryFromTransientTrack pmumuFactory;
                float chi = 0.;
                float ndf = 0.;
                vector<RefCountedKinematicParticle> paiParticles;
                paiParticles.push_back(pmumuFactory.particle(paiTT1, pai_mass, chi, ndf, pai_sigma));
                paiParticles.push_back(pmumuFactory.particle(paiTT2, pai_mass, chi, ndf, pai_sigma));
                paiParticles.push_back(pmumuFactory.particle(paiTT3, pai_mass, chi, ndf, pai_sigma));

                KinematicParticleVertexFitter paifitter;
                RefCountedKinematicTree etaVertexFitTree;
                try
                {
                    etaVertexFitTree = paifitter.fit(paiParticles);
                }
                catch (const std::exception &e)
                {
                    continue;
                }

                if (etaVertexFitTree->isValid())
                {
                    etaVertexFitTree->movePointerToTheTop();
                    RefCountedKinematicParticle Eta_vFit = etaVertexFitTree->currentParticle();
                    RefCountedKinematicVertex Eta_vFit_vertex = etaVertexFitTree->currentDecayVertex();

                    if (!Eta_vFit->currentState().isValid() || !Eta_vFit_vertex->vertexIsValid())
                        continue;
                    if ((Eta_vFit_vertex->chiSquared()) < 0 || (Eta_vFit_vertex->degreesOfFreedom()) <= 0 || (Eta_vFit_vertex->chiSquared()) > 9999.9)
                        continue;

                    double vProb1 = ChiSquaredProbability((double)(Eta_vFit_vertex->chiSquared()), (double)(Eta_vFit_vertex->degreesOfFreedom()));
                    if (vProb1 < JvPorbcut)
                        continue;

                    if (Tcharge == 3)
                    {
                        selectedPai1Charge3.push_back(iPai1 - goodTracks.begin());
                        selectedPai2Charge3.push_back(iPai2 - goodTracks.begin());
                        selectedPai3Charge3.push_back(iPai3 - goodTracks.begin());
                    }
                    else if (Tcharge == -3)
                    {
                        selectedPai1ChargeM3.push_back(iPai1 - goodTracks.begin());
                        selectedPai2ChargeM3.push_back(iPai2 - goodTracks.begin());
                        selectedPai3ChargeM3.push_back(iPai3 - goodTracks.begin());
                    }
                }
            }
        }
    }
step2:

    if (selectedPai1Charge3.size() < 1 || selectedPai1ChargeM3.size() < 1)
        return;

    // Step3: Loop 6 pai according to charge pair zero
    Loop(selectedPai1Charge3, selectedPai2Charge3, selectedPai3Charge3, selectedPai1ChargeM3, selectedPai2ChargeM3, selectedPai3ChargeM3,
         goodTracks,
         dimuonsSort1, dimuonsSort2,
         myEtamass, myEtamasserr, myJpsimass, myJpsimasserr, myExmass, JvPorbcut,
         myMumass, myMumasserr, myPaimass, myPaimasserr, myMesonWindow,
         MassMinCut, pv, bFieldHandle,
         doIso, MuIso, PaiIso, max_loop);

    npairs = possibleInfos.size();
    if (npairs > 0)
    {
        // Sort possibleInfos vector according to X's PT
        std::sort(possibleInfos.begin(), possibleInfos.end(), [](const ParticleInfo &a, const ParticleInfo &b)
                  { return a.XP4.Pt() > b.XP4.Pt(); });

        // Fill first 36 X candidates
        for (int i = 0; (i < npairs) && (i < 36); i++)
        {
            Pai1Charge_[i] = possibleInfos[i].Pai1Charge;
            Pai2Charge_[i] = possibleInfos[i].Pai2Charge;
            Pai3Charge_[i] = possibleInfos[i].Pai3Charge;
            Pai4Charge_[i] = possibleInfos[i].Pai4Charge;
            Pai5Charge_[i] = possibleInfos[i].Pai5Charge;
            Pai6Charge_[i] = possibleInfos[i].Pai6Charge;
            Mu1Charge_[i] = possibleInfos[i].Muon1Charge;
            Mu2Charge_[i] = possibleInfos[i].Muon2Charge;
            Pai1Pt_[i] = possibleInfos[i].Pai1P4.Pt();
            Pai2Pt_[i] = possibleInfos[i].Pai2P4.Pt();
            Pai3Pt_[i] = possibleInfos[i].Pai3P4.Pt();
            Pai4Pt_[i] = possibleInfos[i].Pai4P4.Pt();
            Pai5Pt_[i] = possibleInfos[i].Pai5P4.Pt();
            Pai6Pt_[i] = possibleInfos[i].Pai6P4.Pt();
            Mu1Pt_[i] = possibleInfos[i].Mu1P4.Pt();
            Mu2Pt_[i] = possibleInfos[i].Mu2P4.Pt();
            Pai1Eta_[i] = possibleInfos[i].Pai1P4.Eta();
            Pai2Eta_[i] = possibleInfos[i].Pai2P4.Eta();
            Pai3Eta_[i] = possibleInfos[i].Pai3P4.Eta();
            Pai4Eta_[i] = possibleInfos[i].Pai4P4.Eta();
            Pai5Eta_[i] = possibleInfos[i].Pai5P4.Eta();
            Pai6Eta_[i] = possibleInfos[i].Pai6P4.Eta();
            Mu1Eta_[i] = possibleInfos[i].Mu1P4.Eta();
            Mu2Eta_[i] = possibleInfos[i].Mu2P4.Eta();
            Pai1Phi_[i] = possibleInfos[i].Pai1P4.Phi();
            Pai2Phi_[i] = possibleInfos[i].Pai2P4.Phi();
            Pai3Phi_[i] = possibleInfos[i].Pai3P4.Phi();
            Pai4Phi_[i] = possibleInfos[i].Pai4P4.Phi();
            Pai5Phi_[i] = possibleInfos[i].Pai5P4.Phi();
            Pai6Phi_[i] = possibleInfos[i].Pai6P4.Phi();
            Mu1Phi_[i] = possibleInfos[i].Mu1P4.Phi();
            Mu2Phi_[i] = possibleInfos[i].Mu2P4.Phi();
            Pai1Mass_[i] = possibleInfos[i].Pai1P4.M();
            Pai2Mass_[i] = possibleInfos[i].Pai2P4.M();
            Pai3Mass_[i] = possibleInfos[i].Pai3P4.M();
            Pai4Mass_[i] = possibleInfos[i].Pai4P4.M();
            Pai5Mass_[i] = possibleInfos[i].Pai5P4.M();
            Pai6Mass_[i] = possibleInfos[i].Pai6P4.M();
            Mu1Mass_[i] = possibleInfos[i].Mu1P4.M();
            Mu2Mass_[i] = possibleInfos[i].Mu2P4.M();
            Muon1isSoft_[i] = possibleInfos[i].Muon1isSoft;
            Muon2isSoft_[i] = possibleInfos[i].Muon2isSoft;
            Muon1isTight_[i] = possibleInfos[i].Muon1isTight;
            Muon2isTight_[i] = possibleInfos[i].Muon2isTight;
            Muon1isMedium_[i] = possibleInfos[i].Muon1isMedium;
            Muon2isMedium_[i] = possibleInfos[i].Muon2isMedium;
            Muon1isLoose_[i] = possibleInfos[i].Muon1isLoose;
            Muon2isLoose_[i] = possibleInfos[i].Muon2isLoose;
            Muon1isTracker_[i] = possibleInfos[i].Muon1isTracker;
            Muon2isTracker_[i] = possibleInfos[i].Muon2isTracker;
            Muon1isHighPt_[i] = possibleInfos[i].Muon1isHighPt;
            Muon2isHighPt_[i] = possibleInfos[i].Muon2isHighPt;
            Muon1isGlobal_[i] = possibleInfos[i].Muon1isGlobal;
            Muon2isGlobal_[i] = possibleInfos[i].Muon2isGlobal;
            Pai1Hit_[i] = possibleInfos[i].Pai1Hit;
            Pai2Hit_[i] = possibleInfos[i].Pai2Hit;
            Pai3Hit_[i] = possibleInfos[i].Pai3Hit;
            Pai4Hit_[i] = possibleInfos[i].Pai4Hit;
            Pai5Hit_[i] = possibleInfos[i].Pai5Hit;
            Pai6Hit_[i] = possibleInfos[i].Pai6Hit;
            Pai1fromPV_[i] = possibleInfos[i].Pai1fromPV;
            Pai2fromPV_[i] = possibleInfos[i].Pai2fromPV;
            Pai3fromPV_[i] = possibleInfos[i].Pai3fromPV;
            Pai4fromPV_[i] = possibleInfos[i].Pai4fromPV;
            Pai5fromPV_[i] = possibleInfos[i].Pai5fromPV;
            Pai6fromPV_[i] = possibleInfos[i].Pai6fromPV;
            Pai1Trknormchi2_[i] = possibleInfos[i].Pai1Trknormchi2;
            Pai2Trknormchi2_[i] = possibleInfos[i].Pai2Trknormchi2;
            Pai3Trknormchi2_[i] = possibleInfos[i].Pai3Trknormchi2;
            Pai4Trknormchi2_[i] = possibleInfos[i].Pai4Trknormchi2;
            Pai5Trknormchi2_[i] = possibleInfos[i].Pai5Trknormchi2;
            Pai6Trknormchi2_[i] = possibleInfos[i].Pai6Trknormchi2;
            J1normchi2_[i] = possibleInfos[i].Jpsinormchi2;
            J2normchi2_[i] = possibleInfos[i].Etanormchi2;
            J1NoMassnormchi2_[i] = possibleInfos[i].JpsiNoMassnormchi2;
            J2NoMassnormchi2_[i] = possibleInfos[i].EtaNoMassnormchi2;
            J1Pt_[i] = possibleInfos[i].JpsiP4.Pt();
            J2Pt_[i] = possibleInfos[i].EtaP4.Pt();
            J1NoMassPt_[i] = possibleInfos[i].JpsiNoMassP4.Pt();
            J2NoMassPt_[i] = possibleInfos[i].EtaNoMassP4.Pt();
            J1Eta_[i] = possibleInfos[i].JpsiP4.Eta();
            J2Eta_[i] = possibleInfos[i].EtaP4.Eta();
            J1NoMassEta_[i] = possibleInfos[i].JpsiNoMassP4.Eta();
            J2NoMassEta_[i] = possibleInfos[i].EtaNoMassP4.Eta();
            J1Phi_[i] = possibleInfos[i].JpsiP4.Phi();
            J2Phi_[i] = possibleInfos[i].EtaP4.Phi();
            J1NoMassPhi_[i] = possibleInfos[i].JpsiNoMassP4.Phi();
            J2NoMassPhi_[i] = possibleInfos[i].EtaNoMassP4.Phi();
            J1Mass_[i] = possibleInfos[i].JpsiP4.M();
            J2Mass_[i] = possibleInfos[i].EtaP4.M();
            J1NoMassMass_[i] = possibleInfos[i].JpsiNoMassP4.M();
            J1NoMassMassE_[i] = possibleInfos[i].JpsiNoMassMassE;
            J2NoMassMass_[i] = possibleInfos[i].EtaNoMassP4.M();
            J2NoMassMassE_[i] = possibleInfos[i].EtaNoMassMassE;
            EtaMassDiff_[i] = possibleInfos[i].EtaMassDiff;
            XpT_[i] = possibleInfos[i].XP4.Pt();
            Xeta_[i] = possibleInfos[i].XP4.Eta();
            Xphi_[i] = possibleInfos[i].XP4.Phi();
            Xmass_[i] = possibleInfos[i].XP4.M();
            Xnormchi2_[i] = possibleInfos[i].Xnormchi2;
            minDR_[i] = possibleInfos[i].minDR;
            type_[i] = 1;
        }
        X4muTree->Fill();
        //std::cout << "npairs: " << npairs << " XMass: " << Xmass_[0] << " J1Mass: " << J1Mass_[0] << " J2Mass: " << J2Mass_[0] << " Pai1Mass: " << Pai1Mass_[0] << " Pai2Mass: " << Pai2Mass_[0] << " Pai3Mass: " << Pai3Mass_[0] << " Mu2Mass: " << Mu2Mass_[0] << std::endl;
    }
    npairs = 0;
    clearVars();
}

UInt_t X4muPatSecondaryVertexProducer::getTriggerBits(const edm::Event &iEvent)
{
    UInt_t trigger = 0;
    edm::Handle<edm::TriggerResults> triggerresults;
    iEvent.getByToken(triggerresults_, triggerresults);
    if (triggerresults.isValid())
    {
        const edm::TriggerNames &TheTriggerNames = iEvent.triggerNames(*triggerresults);
        for (unsigned int i = 0; i < FilterNames_.size(); i++)
        {
            bool matched = false;
            for (int version = 1; (version < 99 && (!matched)); version++)
            {
                std::stringstream ss;
                ss << FilterNames_[i] << "_v" << version;
                unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
                if (bit < triggerresults->size() && triggerresults->accept(bit) && !triggerresults->error(bit))
                    matched = true;
            }
            if (matched)
                trigger += (1 << i);
        }
    }
    // else std::cout << "MMrootupler::getTriggerBits: *** NO triggerResults found *** " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
    return trigger;
}

void X4muPatSecondaryVertexProducer::Loop(std::vector<int> selectedPai1, std::vector<int> selectedPai2, std::vector<int> selectedPai3, std::vector<int> selectedPai4, std::vector<int> selectedPai5, std::vector<int> selectedPai6,
                                          std::vector<pat::PackedCandidate> goodTracks,
                                          std::vector<pat::Muon> dimuonsSort1, std::vector<pat::Muon> dimuonsSort2,
                                          double myEtamass, double myEtamasserr, double myJpsimass, double myJpsimasserr, double myExmass, double JvPorbcut,
                                          double myMumass, double myMumasserr, double myPaimass, double myPaimasserr, double myMesonWindow,
                                          double MassMinCut, const reco::Vertex *pv, const MagneticField &bFieldHandle,
                                          bool doIso, double MuIso, double PaiIso, unsigned int max_loop)
{
    long loop4 = 0;
    for (int iPai1N = 0; iPai1N < (int)selectedPai1.size(); iPai1N++)
    {
        for (int iPai4N = 0; iPai4N < (int)selectedPai4.size(); iPai4N++)
        {
            // Check 6 pai overlap
            if (selectedPai1[iPai1N] == selectedPai4[iPai4N] || selectedPai1[iPai1N] == selectedPai5[iPai4N] || selectedPai1[iPai1N] == selectedPai6[iPai4N] || selectedPai2[iPai1N] == selectedPai4[iPai4N] || selectedPai2[iPai1N] == selectedPai5[iPai4N] || selectedPai2[iPai1N] == selectedPai6[iPai4N] || selectedPai3[iPai1N] == selectedPai4[iPai4N] || selectedPai3[iPai1N] == selectedPai5[iPai4N] || selectedPai3[iPai1N] == selectedPai6[iPai4N])
                continue;

            std::vector<pat::PackedCandidate>::iterator iPai1 = goodTracks.begin() + selectedPai1[iPai1N];
            std::vector<pat::PackedCandidate>::iterator iPai2 = goodTracks.begin() + selectedPai2[iPai1N];
            std::vector<pat::PackedCandidate>::iterator iPai3 = goodTracks.begin() + selectedPai3[iPai1N];
            std::vector<pat::PackedCandidate>::iterator iPai4 = goodTracks.begin() + selectedPai4[iPai4N];
            std::vector<pat::PackedCandidate>::iterator iPai5 = goodTracks.begin() + selectedPai5[iPai4N];
            std::vector<pat::PackedCandidate>::iterator iPai6 = goodTracks.begin() + selectedPai6[iPai4N];

            if (iPai1->charge() + iPai2->charge() + iPai3->charge() + iPai4->charge() + iPai5->charge() + iPai6->charge() != 0)
                continue;
            if (iPai1->pt() < 0.5 || iPai2->pt() < 0.5 || iPai3->pt() < 0.5 || iPai4->pt() < 0.5 || iPai5->pt() < 0.5 || iPai6->pt() < 0.5)
                continue;
            if (doIso && (abs(deltaR(*iPai1, *iPai2)) < PaiIso || abs(deltaR(*iPai1, *iPai3)) < PaiIso || abs(deltaR(*iPai1, *iPai4)) < PaiIso || abs(deltaR(*iPai1, *iPai5)) < PaiIso || abs(deltaR(*iPai1, *iPai6)) < PaiIso || abs(deltaR(*iPai2, *iPai3)) < PaiIso || abs(deltaR(*iPai2, *iPai4)) < PaiIso || abs(deltaR(*iPai2, *iPai5)) < PaiIso || abs(deltaR(*iPai2, *iPai6)) < PaiIso || abs(deltaR(*iPai3, *iPai4)) < PaiIso || abs(deltaR(*iPai3, *iPai5)) < PaiIso || abs(deltaR(*iPai3, *iPai6)) < PaiIso || abs(deltaR(*iPai4, *iPai5)) < PaiIso || abs(deltaR(*iPai4, *iPai6)) < PaiIso || abs(deltaR(*iPai5, *iPai6)) < PaiIso))
                continue;
            // Calculate the invariant mass of 6 pais
            TLorentzVector pai1, pai2, pai3, pai4, pai5, pai6;
            pai1.SetPtEtaPhiM(iPai1->pt(), iPai1->eta(), iPai1->phi(), myPaimass);
            pai2.SetPtEtaPhiM(iPai2->pt(), iPai2->eta(), iPai2->phi(), myPaimass);
            pai3.SetPtEtaPhiM(iPai3->pt(), iPai3->eta(), iPai3->phi(), myPaimass);
            pai4.SetPtEtaPhiM(iPai4->pt(), iPai4->eta(), iPai4->phi(), myPaimass);
            pai5.SetPtEtaPhiM(iPai5->pt(), iPai5->eta(), iPai5->phi(), myPaimass);
            pai6.SetPtEtaPhiM(iPai6->pt(), iPai6->eta(), iPai6->phi(), myPaimass);
            TLorentzVector X = pai1 + pai2 + pai3 + pai4 + pai5 + pai6;
            if (abs(X.M() - myEtamass) > myMesonWindow) // eta mass
                continue;

            //6pai fit
            //Store 6 pai mass

            for (std::vector<pat::Muon>::iterator iMuon1 = dimuonsSort1.begin(); iMuon1 != dimuonsSort1.end(); ++iMuon1)
            {
                std::vector<pat::Muon>::iterator iMuon2 = dimuonsSort2.begin() + (iMuon1 - dimuonsSort1.begin());
                loop4++;
                if (loop4 > max_loop)
                {
                    // std::cout << "loop3 > max_loop" << std::endl;
                    return;
                }

                if (iMuon1->charge() + iMuon2->charge() != 0)
                    continue;
                if (iMuon1->pt() < 0.5 || iMuon2->pt() < 0.5)
                    continue;
                if (doIso && (abs(deltaR(*iMuon1, *iPai1)) < MuIso || abs(deltaR(*iMuon1, *iPai2)) < MuIso || abs(deltaR(*iMuon1, *iPai3)) < MuIso || abs(deltaR(*iMuon1, *iPai4)) < MuIso || abs(deltaR(*iMuon1, *iPai5)) < MuIso || abs(deltaR(*iMuon1, *iPai6)) < MuIso || abs(deltaR(*iMuon2, *iPai1)) < MuIso || abs(deltaR(*iMuon2, *iPai2)) < MuIso || abs(deltaR(*iMuon2, *iPai3)) < MuIso || abs(deltaR(*iMuon2, *iPai4)) < MuIso || abs(deltaR(*iMuon2, *iPai5)) < MuIso || abs(deltaR(*iMuon2, *iPai6)) < MuIso || abs(deltaR(*iMuon1, *iMuon2)) < MuIso))
                    continue;

                reco::Track PaiTrack1 = *(iPai1->bestTrack());
                reco::Track PaiTrack2 = *(iPai2->bestTrack());
                reco::Track PaiTrack3 = *(iPai3->bestTrack());
                reco::Track PaiTrack4 = *(iPai4->bestTrack());
                reco::Track PaiTrack5 = *(iPai5->bestTrack());
                reco::Track PaiTrack6 = *(iPai6->bestTrack());
                TrackRef MuonTrack1 = iMuon1->track();
                TrackRef MuonTrack2 = iMuon2->track();

                TransientTrack paiTT1(PaiTrack1, &(bFieldHandle));
                TransientTrack paiTT2(PaiTrack2, &(bFieldHandle));
                TransientTrack paiTT3(PaiTrack3, &(bFieldHandle));
                TransientTrack paiTT4(PaiTrack4, &(bFieldHandle));
                TransientTrack paiTT5(PaiTrack5, &(bFieldHandle));
                TransientTrack paiTT6(PaiTrack6, &(bFieldHandle));
                TransientTrack muonTT1(MuonTrack1, &(bFieldHandle));
                TransientTrack muonTT2(MuonTrack2, &(bFieldHandle));

                if (!paiTT1.isValid() || !paiTT2.isValid() || !paiTT3.isValid() || !paiTT4.isValid() || !paiTT5.isValid() || !paiTT6.isValid() || !muonTT1.isValid() || !muonTT2.isValid())
                    continue;

                ParticleMass pai_mass = myPaimass; // pdg mass
                float pai_sigma = myPaimasserr;
                ParticleMass mu_mass = myMumass;
                float mu_sigma = myMumasserr;

                KinematicParticleFactoryFromTransientTrack pmumuFactory;

                // initial chi2 and ndf before kinematic fits.
                float chi1 = 0.;
                float ndf1 = 0.;
                float chi2 = 0.;
                float ndf2 = 0.;
                float chi3 = 0.;
                float ndf3 = 0.;

                vector<RefCountedKinematicParticle> PaiParticles, MuParticles, AllParticles;
                PaiParticles.push_back(pmumuFactory.particle(paiTT1, pai_mass, chi1, ndf1, pai_sigma));
                PaiParticles.push_back(pmumuFactory.particle(paiTT2, pai_mass, chi1, ndf1, pai_sigma));
                PaiParticles.push_back(pmumuFactory.particle(paiTT3, pai_mass, chi1, ndf1, pai_sigma));
                PaiParticles.push_back(pmumuFactory.particle(paiTT4, pai_mass, chi1, ndf1, pai_sigma));
                PaiParticles.push_back(pmumuFactory.particle(paiTT5, pai_mass, chi1, ndf1, pai_sigma));
                PaiParticles.push_back(pmumuFactory.particle(paiTT6, pai_mass, chi1, ndf1, pai_sigma));
                MuParticles.push_back(pmumuFactory.particle(muonTT1, mu_mass, chi2, ndf2, mu_sigma));
                MuParticles.push_back(pmumuFactory.particle(muonTT2, mu_mass, chi2, ndf2, mu_sigma));
                AllParticles.push_back(pmumuFactory.particle(paiTT1, pai_mass, chi3, ndf3, pai_sigma));
                AllParticles.push_back(pmumuFactory.particle(paiTT2, pai_mass, chi3, ndf3, pai_sigma));
                AllParticles.push_back(pmumuFactory.particle(paiTT3, pai_mass, chi3, ndf3, pai_sigma));
                AllParticles.push_back(pmumuFactory.particle(paiTT4, pai_mass, chi3, ndf3, pai_sigma));
                AllParticles.push_back(pmumuFactory.particle(paiTT5, pai_mass, chi3, ndf3, pai_sigma));
                AllParticles.push_back(pmumuFactory.particle(paiTT6, pai_mass, chi3, ndf3, pai_sigma));
                AllParticles.push_back(pmumuFactory.particle(muonTT1, mu_mass, chi3, ndf3, mu_sigma));
                AllParticles.push_back(pmumuFactory.particle(muonTT2, mu_mass, chi3, ndf3, mu_sigma));

                if (PaiParticles.size() < 6 || MuParticles.size() < 2 || AllParticles.size() < 8)
                    continue;

                KinematicParticleVertexFitter fitter1, fitter2, Allfitter;
                RefCountedKinematicTree etaVertexFitTree, jpsiVertexFitTree, XVertexFitTree;
                try
                {
                    etaVertexFitTree = fitter1.fit(PaiParticles);
                    jpsiVertexFitTree = fitter2.fit(MuParticles);
                    XVertexFitTree = Allfitter.fit(AllParticles);
                }
                catch (const std::exception &e)
                {
                    continue;
                    // std::cout << "Exception thrown in fitting" << std::endl;
                }

                if (XVertexFitTree->isValid() && etaVertexFitTree->isValid() && jpsiVertexFitTree->isValid())
                {
                    XVertexFitTree->movePointerToTheTop();
                    etaVertexFitTree->movePointerToTheTop();
                    jpsiVertexFitTree->movePointerToTheTop();

                    RefCountedKinematicParticle eta_vFit = etaVertexFitTree->currentParticle();
                    RefCountedKinematicParticle jpsi_vFit = jpsiVertexFitTree->currentParticle();
                    RefCountedKinematicParticle X_vFit = XVertexFitTree->currentParticle();

                    XVertexFitTree->movePointerToTheTop();
                    etaVertexFitTree->movePointerToTheTop();
                    jpsiVertexFitTree->movePointerToTheTop();

                    RefCountedKinematicVertex eta_vFit_vertex = etaVertexFitTree->currentDecayVertex();
                    RefCountedKinematicVertex jpsi_vFit_vertex = jpsiVertexFitTree->currentDecayVertex();
                    RefCountedKinematicVertex X_vFit_vertex = XVertexFitTree->currentDecayVertex();

                    if (!eta_vFit->currentState().isValid() || !jpsi_vFit->currentState().isValid() || !X_vFit->currentState().isValid())
                        continue;
                    if (!eta_vFit_vertex->vertexIsValid() || !jpsi_vFit_vertex->vertexIsValid() || !X_vFit_vertex->vertexIsValid())
                        continue;

                    XVertexFitTree->movePointerToTheFirstChild();
                    RefCountedKinematicParticle Xpai1 = XVertexFitTree->currentParticle();
                    XVertexFitTree->movePointerToTheNextChild();
                    RefCountedKinematicParticle Xpai2 = XVertexFitTree->currentParticle();
                    XVertexFitTree->movePointerToTheNextChild();
                    RefCountedKinematicParticle Xpai3 = XVertexFitTree->currentParticle();
                    XVertexFitTree->movePointerToTheNextChild();
                    RefCountedKinematicParticle Xpai4 = XVertexFitTree->currentParticle();
                    XVertexFitTree->movePointerToTheNextChild();
                    RefCountedKinematicParticle Xpai5 = XVertexFitTree->currentParticle();
                    XVertexFitTree->movePointerToTheNextChild();
                    RefCountedKinematicParticle Xpai6 = XVertexFitTree->currentParticle();
                    XVertexFitTree->movePointerToTheNextChild();
                    RefCountedKinematicParticle Xmu1 = XVertexFitTree->currentParticle();
                    XVertexFitTree->movePointerToTheNextChild();
                    RefCountedKinematicParticle Xmu2 = XVertexFitTree->currentParticle();

                    if (!Xpai1->currentState().isValid() || !Xpai2->currentState().isValid() || !Xpai3->currentState().isValid() || !Xpai4->currentState().isValid() || !Xpai5->currentState().isValid() || !Xpai6->currentState().isValid() || !Xmu1->currentState().isValid() || !Xmu2->currentState().isValid())
                        continue;
                    if ((X_vFit_vertex->chiSquared()) < 0 || (X_vFit_vertex->degreesOfFreedom()) <= 0 || (X_vFit_vertex->chiSquared()) > 9999.9)
                        continue;
                    if ((eta_vFit_vertex->chiSquared()) < 0 || (eta_vFit_vertex->degreesOfFreedom()) <= 0 || (eta_vFit_vertex->chiSquared()) > 9999.9)
                        continue;
                    if ((jpsi_vFit_vertex->chiSquared()) < 0 || (jpsi_vFit_vertex->degreesOfFreedom()) <= 0 || (jpsi_vFit_vertex->chiSquared()) > 9999.9)
                        continue;

                    // std::cout << "X mass: " << X_vFit->currentState().mass() << std::endl;

                    double vProbE = ChiSquaredProbability((double)(eta_vFit_vertex->chiSquared()), (double)(eta_vFit_vertex->degreesOfFreedom()));
                    double vProbJ = ChiSquaredProbability((double)(jpsi_vFit_vertex->chiSquared()), (double)(jpsi_vFit_vertex->degreesOfFreedom()));
                    double vProbX = ChiSquaredProbability((double)(X_vFit_vertex->chiSquared()), (double)(X_vFit_vertex->degreesOfFreedom()));
                    if (vProbJ < JvPorbcut || vProbE < JvPorbcut || vProbX < JvPorbcut)
                        continue;
                    if (Xpai1->currentState().mass() <= MassMinCut || Xpai2->currentState().mass() <= MassMinCut || Xpai3->currentState().mass() <= MassMinCut || Xpai4->currentState().mass() <= MassMinCut || Xpai5->currentState().mass() <= MassMinCut || Xpai6->currentState().mass() <= MassMinCut || Xmu1->currentState().mass() <= MassMinCut || Xmu2->currentState().mass() <= MassMinCut)
                        continue;
                    if (eta_vFit->currentState().kinematicParametersError().matrix()(6, 6) < 0 || jpsi_vFit->currentState().kinematicParametersError().matrix()(6, 6) < 0 || X_vFit->currentState().kinematicParametersError().matrix()(6, 6) < 0)
                        continue;
                    float Etamasserr = sqrt(eta_vFit->currentState().kinematicParametersError().matrix()(6, 6));
                    float Jpsimasserr = sqrt(jpsi_vFit->currentState().kinematicParametersError().matrix()(6, 6));

                    // std::cout << "Peobe mass" << eta_vFit->currentState().mass() << " +/- " << Etamasserr << std::endl;

                    if (abs(eta_vFit->currentState().mass() - myExmass) < 2.0 * Etamasserr || abs(jpsi_vFit->currentState().mass() - myExmass) < 2.0 * Jpsimasserr)
                        continue;
                    if (abs(eta_vFit->currentState().mass() - myEtamass) > 3.0 * Etamasserr || abs(jpsi_vFit->currentState().mass() - myJpsimass) > 3.0 * Jpsimasserr || abs(eta_vFit->currentState().mass() - myEtamass) > myMesonWindow || abs(jpsi_vFit->currentState().mass() - myJpsimass) > myMesonWindow)
                        continue;

                    // std::cout << "Eta mass: " << eta_vFit->currentState().mass() << " +/- " << Etamasserr << std::endl;

                    KinematicParameters EtaNoMassKtmp = eta_vFit->currentState().kinematicParameters();
                    KinematicParameters JpsiNoMassKtmp = jpsi_vFit->currentState().kinematicParameters();
                    TLorentzVector EtaNoMassP4tmp;
                    EtaNoMassP4tmp.SetPxPyPzE(EtaNoMassKtmp.momentum().x(), EtaNoMassKtmp.momentum().y(), EtaNoMassKtmp.momentum().z(), sqrt(EtaNoMassKtmp.momentum().mag2() + eta_vFit->currentState().mass() * eta_vFit->currentState().mass()));
                    TLorentzVector JpsiNoMassP4tmp;
                    JpsiNoMassP4tmp.SetPxPyPzE(JpsiNoMassKtmp.momentum().x(), JpsiNoMassKtmp.momentum().y(), JpsiNoMassKtmp.momentum().z(), sqrt(JpsiNoMassKtmp.momentum().mag2() + jpsi_vFit->currentState().mass() * jpsi_vFit->currentState().mass()));

                    // Mass Constraint Fit for J/psi and Eta_c
                    KinematicParticleVertexFitter kpvFitter;
                    KinematicParticleFitter csFitter;
                    ParticleMass jpsiM = myJpsimass;
                    float jpsi_m_sigma = myJpsimasserr;
                    ParticleMass etaM = myEtamass;
                    float eta_m_sigma = myEtamasserr;
                    KinematicConstraint *jpsi_c = new MassKinematicConstraint(jpsiM, jpsi_m_sigma);
                    KinematicConstraint *eta_c = new MassKinematicConstraint(etaM, eta_m_sigma);
                    vector<RefCountedKinematicParticle> muonP, paiP;
                    RefCountedKinematicParticle Jpsi_part, Eta_part;
                    vector<RefCountedKinematicParticle> Chi_1;
                    RefCountedKinematicTree Chi1_bTree;
                    RefCountedKinematicParticle MyChi1_part;
                    muonP.push_back(pmumuFactory.particle(muonTT1, mu_mass, chi2, ndf2, mu_sigma));
                    muonP.push_back(pmumuFactory.particle(muonTT2, mu_mass, chi2, ndf2, mu_sigma));
                    paiP.push_back(pmumuFactory.particle(paiTT1, pai_mass, chi1, ndf1, pai_sigma));
                    paiP.push_back(pmumuFactory.particle(paiTT2, pai_mass, chi1, ndf1, pai_sigma));
                    paiP.push_back(pmumuFactory.particle(paiTT3, pai_mass, chi1, ndf1, pai_sigma));
                    paiP.push_back(pmumuFactory.particle(paiTT4, pai_mass, chi1, ndf1, pai_sigma));
                    paiP.push_back(pmumuFactory.particle(paiTT5, pai_mass, chi1, ndf1, pai_sigma));
                    paiP.push_back(pmumuFactory.particle(paiTT6, pai_mass, chi1, ndf1, pai_sigma));
                    RefCountedKinematicTree JpsiTree;
                    RefCountedKinematicTree EtaTree;
                    RefCountedKinematicTree JpsinoMCJJ;
                    RefCountedKinematicTree EtanoMCJJ;
                    try
                    {
                        JpsinoMCJJ = kpvFitter.fit(muonP);
                        EtanoMCJJ = kpvFitter.fit(paiP);
                    }
                    catch (const std::exception &e)
                    {
                        continue;
                    }
                    if (JpsinoMCJJ->isEmpty() || EtanoMCJJ->isEmpty())
                        continue;
                    try
                    {
                        JpsiTree = csFitter.fit(jpsi_c, JpsinoMCJJ);
                        EtaTree = csFitter.fit(eta_c, EtanoMCJJ);
                    }
                    catch (const std::exception &e)
                    {
                        std::cout << "dimuon and six pais vertex exception with dimuon and six pais constrained to Jpsi + Eta" << std::endl;
                        continue;
                    }
                    if (JpsiTree->isValid() && EtaTree->isValid())
                    {
                        JpsiTree->movePointerToTheTop();
                        EtaTree->movePointerToTheTop();
                        Jpsi_part = JpsiTree->currentParticle();
                        Eta_part = EtaTree->currentParticle();
                        Chi_1.push_back(Eta_part);
                        Chi_1.push_back(Jpsi_part);
                        bool isagoodfit = true;
                        try
                        {
                            Chi1_bTree = kpvFitter.fit(Chi_1);
                        }
                        catch (const std::exception &e)
                        {
                            isagoodfit = false;
                            std::cout << "exception with Jpsi + Eta vertex fit" << std::endl;
                            continue;
                        }
                        if (Chi1_bTree->isValid() && isagoodfit)
                        {
                            Chi1_bTree->movePointerToTheTop();
                            RefCountedKinematicVertex myEtaJpsiVertex = Chi1_bTree->currentDecayVertex();
                            Chi1_bTree->movePointerToTheFirstChild();
                            RefCountedKinematicParticle etaInXVertex = Chi1_bTree->currentParticle();
                            Chi1_bTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle jpsiInXVertex = Chi1_bTree->currentParticle();

                            JpsiTree->movePointerToTheFirstChild();
                            RefCountedKinematicParticle XMassmu1 = JpsiTree->currentParticle();
                            JpsiTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle XMassmu2 = JpsiTree->currentParticle();
                            EtaTree->movePointerToTheFirstChild();
                            RefCountedKinematicParticle XMasspai1 = EtaTree->currentParticle();
                            EtaTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle XMasspai2 = EtaTree->currentParticle();
                            EtaTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle XMasspai3 = EtaTree->currentParticle();
                            EtaTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle XMasspai4 = EtaTree->currentParticle();
                            EtaTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle XMasspai5 = EtaTree->currentParticle();
                            EtaTree->movePointerToTheNextChild();
                            RefCountedKinematicParticle XMasspai6 = EtaTree->currentParticle();

                            if (!XMassmu1->currentState().isValid() || !XMassmu2->currentState().isValid() || !XMasspai1->currentState().isValid() || !XMasspai2->currentState().isValid() || !XMasspai3->currentState().isValid() || !XMasspai4->currentState().isValid() || !XMasspai5->currentState().isValid() || !XMasspai6->currentState().isValid())
                                continue;

                            JpsiTree->movePointerToTheTop();
                            EtaTree->movePointerToTheTop();
                            RefCountedKinematicVertex myEtaDecayVtx = EtaTree->currentDecayVertex();
                            RefCountedKinematicVertex myJpsiDecayVtx = JpsiTree->currentDecayVertex();

                            // std::cout << "Eta Chi2 " << myEtaDecayVtx->chiSquared() << " Eta degreesOfFreedom " << myEtaDecayVtx->degreesOfFreedom() << std::endl;

                            if ((myEtaDecayVtx->chiSquared()) < 0 || (myEtaDecayVtx->degreesOfFreedom()) <= 0 || (myEtaDecayVtx->chiSquared()) > 9999.9)
                                continue;
                            if ((myJpsiDecayVtx->chiSquared()) < 0 || (myJpsiDecayVtx->degreesOfFreedom()) <= 0 || (myJpsiDecayVtx->chiSquared()) > 9999.9)
                                continue;
                            if ((myEtaJpsiVertex->chiSquared()) < 0 || (myEtaJpsiVertex->degreesOfFreedom()) <= 0 || (myEtaJpsiVertex->chiSquared()) > 9999.9)
                                continue;

                            double vProbEMC = ChiSquaredProbability((double)(myEtaDecayVtx->chiSquared()), (double)(myEtaDecayVtx->degreesOfFreedom()));
                            double vProbJMC = ChiSquaredProbability((double)(myJpsiDecayVtx->chiSquared()), (double)(myJpsiDecayVtx->degreesOfFreedom()));
                            double vProbXMC = ChiSquaredProbability((double)(myEtaJpsiVertex->chiSquared()), (double)(myEtaJpsiVertex->degreesOfFreedom()));
                            if (vProbJMC < JvPorbcut || vProbEMC < JvPorbcut || vProbXMC < JvPorbcut)
                                continue;
                            if (XMassmu1->currentState().mass() <= MassMinCut || XMassmu2->currentState().mass() <= MassMinCut || XMasspai1->currentState().mass() <= MassMinCut || XMasspai2->currentState().mass() <= MassMinCut || XMasspai3->currentState().mass() <= MassMinCut || XMasspai4->currentState().mass() <= MassMinCut || XMasspai5->currentState().mass() <= MassMinCut || XMasspai6->currentState().mass() <= MassMinCut)
                                continue;

                            KinematicParameters mu1Ktmp, mu2Ktmp, pai1Ktmp, pai2Ktmp, pai3Ktmp, pai4Ktmp, pai5Ktmp, pai6Ktmp, JpsiKtmp, EtaKtmp, XKtmp;
                            TLorentzVector mu1P4tmp, mu2P4tmp, pai1P4tmp, pai2P4tmp, pai3P4tmp, pai4P4tmp, pai5P4tmp, pai6P4tmp, JpsiP4tmp, EtaP4tmp, XP4tmp;
                            mu1Ktmp = XMassmu1->currentState().kinematicParameters();
                            mu2Ktmp = XMassmu2->currentState().kinematicParameters();
                            pai1Ktmp = XMasspai1->currentState().kinematicParameters();
                            pai2Ktmp = XMasspai2->currentState().kinematicParameters();
                            pai3Ktmp = XMasspai3->currentState().kinematicParameters();
                            pai4Ktmp = XMasspai4->currentState().kinematicParameters();
                            pai5Ktmp = XMasspai5->currentState().kinematicParameters();
                            pai6Ktmp = XMasspai6->currentState().kinematicParameters();
                            JpsiKtmp = jpsiInXVertex->currentState().kinematicParameters();
                            EtaKtmp = etaInXVertex->currentState().kinematicParameters();
                            Chi1_bTree->movePointerToTheTop();
                            MyChi1_part = Chi1_bTree->currentParticle();
                            XKtmp = MyChi1_part->currentState().kinematicParameters();
                            mu1P4tmp.SetPxPyPzE(mu1Ktmp.momentum().x(), mu1Ktmp.momentum().y(), mu1Ktmp.momentum().z(), sqrt(mu1Ktmp.momentum().mag2() + XMassmu1->currentState().mass() * XMassmu1->currentState().mass()));
                            mu2P4tmp.SetPxPyPzE(mu2Ktmp.momentum().x(), mu2Ktmp.momentum().y(), mu2Ktmp.momentum().z(), sqrt(mu2Ktmp.momentum().mag2() + XMassmu2->currentState().mass() * XMassmu2->currentState().mass()));
                            pai1P4tmp.SetPxPyPzE(pai1Ktmp.momentum().x(), pai1Ktmp.momentum().y(), pai1Ktmp.momentum().z(), sqrt(pai1Ktmp.momentum().mag2() + XMasspai1->currentState().mass() * XMasspai1->currentState().mass()));
                            pai2P4tmp.SetPxPyPzE(pai2Ktmp.momentum().x(), pai2Ktmp.momentum().y(), pai2Ktmp.momentum().z(), sqrt(pai2Ktmp.momentum().mag2() + XMasspai2->currentState().mass() * XMasspai2->currentState().mass()));
                            pai3P4tmp.SetPxPyPzE(pai3Ktmp.momentum().x(), pai3Ktmp.momentum().y(), pai3Ktmp.momentum().z(), sqrt(pai3Ktmp.momentum().mag2() + XMasspai3->currentState().mass() * XMasspai3->currentState().mass()));
                            pai4P4tmp.SetPxPyPzE(pai4Ktmp.momentum().x(), pai4Ktmp.momentum().y(), pai4Ktmp.momentum().z(), sqrt(pai4Ktmp.momentum().mag2() + XMasspai4->currentState().mass() * XMasspai4->currentState().mass()));
                            pai5P4tmp.SetPxPyPzE(pai5Ktmp.momentum().x(), pai5Ktmp.momentum().y(), pai5Ktmp.momentum().z(), sqrt(pai5Ktmp.momentum().mag2() + XMasspai5->currentState().mass() * XMasspai5->currentState().mass()));
                            pai6P4tmp.SetPxPyPzE(pai6Ktmp.momentum().x(), pai6Ktmp.momentum().y(), pai6Ktmp.momentum().z(), sqrt(pai6Ktmp.momentum().mag2() + XMasspai6->currentState().mass() * XMasspai6->currentState().mass()));
                            JpsiP4tmp.SetPxPyPzE(JpsiKtmp.momentum().x(), JpsiKtmp.momentum().y(), JpsiKtmp.momentum().z(), sqrt(JpsiKtmp.momentum().mag2() + jpsiInXVertex->currentState().mass() * jpsiInXVertex->currentState().mass()));
                            EtaP4tmp.SetPxPyPzE(EtaKtmp.momentum().x(), EtaKtmp.momentum().y(), EtaKtmp.momentum().z(), sqrt(EtaKtmp.momentum().mag2() + etaInXVertex->currentState().mass() * etaInXVertex->currentState().mass()));
                            XP4tmp.SetPxPyPzE(XKtmp.momentum().x(), XKtmp.momentum().y(), XKtmp.momentum().z(), sqrt(XKtmp.momentum().mag2() + MyChi1_part->currentState().mass() * MyChi1_part->currentState().mass()));
                            double Etanormchi2 = (double)myEtaDecayVtx->chiSquared() / (double)myEtaDecayVtx->degreesOfFreedom();
                            double Jpsinormchi2 = (double)myJpsiDecayVtx->chiSquared() / (double)myJpsiDecayVtx->degreesOfFreedom();
                            double Xnormchi2 = (double)myEtaJpsiVertex->chiSquared() / (double)myEtaJpsiVertex->degreesOfFreedom();
                            // Store Candidate in vector
                            ParticleInfo Xtmp;

                            Xtmp.Pai1Charge = iPai1->charge();
                            Xtmp.Pai2Charge = iPai2->charge();
                            Xtmp.Pai3Charge = iPai3->charge();
                            Xtmp.Pai4Charge = iPai4->charge();
                            Xtmp.Pai5Charge = iPai5->charge();
                            Xtmp.Pai6Charge = iPai6->charge();
                            Xtmp.Muon1Charge = iMuon1->charge();
                            Xtmp.Muon2Charge = iMuon2->charge();
                            Xtmp.Pai1P4 = pai1P4tmp;
                            Xtmp.Pai2P4 = pai2P4tmp;
                            Xtmp.Pai3P4 = pai3P4tmp;
                            Xtmp.Pai4P4 = pai4P4tmp;
                            Xtmp.Pai5P4 = pai5P4tmp;
                            Xtmp.Pai6P4 = pai6P4tmp;
                            Xtmp.Mu1P4 = mu1P4tmp;
                            Xtmp.Mu2P4 = mu2P4tmp;
                            Xtmp.JpsiP4 = JpsiP4tmp;
                            Xtmp.EtaP4 = EtaP4tmp;
                            Xtmp.XP4 = XP4tmp;
                            Xtmp.EtaNoMassP4 = EtaNoMassP4tmp;
                            Xtmp.JpsiNoMassP4 = JpsiNoMassP4tmp;

                            const reco::Muon &muon1 = *iMuon1;
                            const reco::Muon &muon2 = *iMuon2;
                            const reco::Vertex &vertex = *pv;
                            Xtmp.Muon1isSoft = (bool)muon::isSoftMuon(muon1, vertex, false);
                            Xtmp.Muon2isSoft = (bool)muon::isSoftMuon(muon2, vertex, false);
                            Xtmp.Muon1isLoose = (bool)muon::isLooseMuon(muon1);
                            Xtmp.Muon2isLoose = (bool)muon::isLooseMuon(muon2);
                            Xtmp.Muon1isMedium = (bool)muon::isMediumMuon(muon1);
                            Xtmp.Muon2isMedium = (bool)muon::isMediumMuon(muon2);
                            Xtmp.Muon1isTight = (bool)muon::isTightMuon(muon1, vertex);
                            Xtmp.Muon2isTight = (bool)muon::isTightMuon(muon2, vertex);
                            Xtmp.Muon1isHighPt = (bool)muon::isHighPtMuon(muon1, vertex);
                            Xtmp.Muon2isHighPt = (bool)muon::isHighPtMuon(muon2, vertex);
                            Xtmp.Muon1isGlobal = (bool)muon1.isGlobalMuon();
                            Xtmp.Muon2isGlobal = (bool)muon2.isGlobalMuon();
                            Xtmp.Muon1isTracker = (bool)muon1.isTrackerMuon();
                            Xtmp.Muon2isTracker = (bool)muon2.isTrackerMuon();
                            Xtmp.Pai1Hit = iPai1->numberOfHits();
                            Xtmp.Pai2Hit = iPai2->numberOfHits();
                            Xtmp.Pai3Hit = iPai3->numberOfHits();
                            Xtmp.Pai4Hit = iPai4->numberOfHits();
                            Xtmp.Pai5Hit = iPai5->numberOfHits();
                            Xtmp.Pai6Hit = iPai6->numberOfHits();
                            Xtmp.Pai1fromPV = iPai1->fromPV();
                            Xtmp.Pai2fromPV = iPai2->fromPV();
                            Xtmp.Pai3fromPV = iPai3->fromPV();
                            Xtmp.Pai4fromPV = iPai4->fromPV();
                            Xtmp.Pai5fromPV = iPai5->fromPV();
                            Xtmp.Pai6fromPV = iPai6->fromPV();
                            Xtmp.Pai1Trknormchi2 = iPai1->bestTrack()->normalizedChi2();
                            Xtmp.Pai2Trknormchi2 = iPai2->bestTrack()->normalizedChi2();
                            Xtmp.Pai3Trknormchi2 = iPai3->bestTrack()->normalizedChi2();
                            Xtmp.Pai4Trknormchi2 = iPai4->bestTrack()->normalizedChi2();
                            Xtmp.Pai5Trknormchi2 = iPai5->bestTrack()->normalizedChi2();
                            Xtmp.Pai6Trknormchi2 = iPai6->bestTrack()->normalizedChi2();

                            Xtmp.Etanormchi2 = Etanormchi2;
                            Xtmp.Jpsinormchi2 = Jpsinormchi2;
                            Xtmp.Xnormchi2 = Xnormchi2;
                            Xtmp.EtaNoMassnormchi2 = (double)(eta_vFit_vertex->chiSquared()) / (double)(eta_vFit_vertex->degreesOfFreedom());
                            Xtmp.JpsiNoMassnormchi2 = (double)(jpsi_vFit_vertex->chiSquared()) / (double)(jpsi_vFit_vertex->degreesOfFreedom());
                            Xtmp.EtaNoMassMassE = Etamasserr;
                            Xtmp.JpsiNoMassMassE = Jpsimasserr;
                            Xtmp.EtaMassDiff = X.M() - eta_vFit->currentState().mass();

                            double minDR_tmp = 999;
                            if (abs(deltaR(*iPai1, *iPai2)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai1, *iPai2));
                            if (abs(deltaR(*iPai1, *iPai3)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai1, *iPai3));
                            if (abs(deltaR(*iPai1, *iPai4)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai1, *iPai4));
                            if (abs(deltaR(*iPai1, *iPai5)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai1, *iPai5));
                            if (abs(deltaR(*iPai1, *iPai6)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai1, *iPai6));
                            if (abs(deltaR(*iPai2, *iPai3)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai2, *iPai3));
                            if (abs(deltaR(*iPai2, *iPai4)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai2, *iPai4));
                            if (abs(deltaR(*iPai2, *iPai5)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai2, *iPai5));
                            if (abs(deltaR(*iPai2, *iPai6)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai2, *iPai6));
                            if (abs(deltaR(*iPai3, *iPai4)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai3, *iPai4));
                            if (abs(deltaR(*iPai3, *iPai5)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai3, *iPai5));
                            if (abs(deltaR(*iPai3, *iPai6)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai3, *iPai6));
                            if (abs(deltaR(*iPai4, *iPai5)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai4, *iPai5));
                            if (abs(deltaR(*iPai4, *iPai6)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai4, *iPai6));
                            if (abs(deltaR(*iPai5, *iPai6)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai5, *iPai6));
                            if (abs(deltaR(*iPai1, *iMuon1)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai1, *iMuon1));
                            if (abs(deltaR(*iPai1, *iMuon2)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai1, *iMuon2));
                            if (abs(deltaR(*iPai2, *iMuon1)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai2, *iMuon1));
                            if (abs(deltaR(*iPai2, *iMuon2)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai2, *iMuon2));
                            if (abs(deltaR(*iPai3, *iMuon1)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai3, *iMuon1));
                            if (abs(deltaR(*iPai3, *iMuon2)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai3, *iMuon2));
                            if (abs(deltaR(*iPai4, *iMuon1)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai4, *iMuon1));
                            if (abs(deltaR(*iPai4, *iMuon2)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai4, *iMuon2));
                            if (abs(deltaR(*iPai5, *iMuon1)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai5, *iMuon1));
                            if (abs(deltaR(*iPai5, *iMuon2)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai5, *iMuon2));
                            if (abs(deltaR(*iPai6, *iMuon1)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai6, *iMuon1));
                            if (abs(deltaR(*iPai6, *iMuon2)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iPai6, *iMuon2));
                            if (abs(deltaR(*iMuon1, *iMuon2)) < minDR_tmp)
                                minDR_tmp = abs(deltaR(*iMuon1, *iMuon2));
                            Xtmp.minDR = minDR_tmp;

                            // Fill
                            possibleInfos.push_back(Xtmp);
                        }
                    }
                }
            }
        }
    }
}

void X4muPatSecondaryVertexProducer::clearVars()
{
    run = 0;
    event = 0;
    lumiblock = 0;
    trigger = 0;
    possibleInfos.clear();

    for (size_t i = 0; i < 36; i++)
    {
        Pai1Charge_[i] = 0;
        Pai2Charge_[i] = 0;
        Pai3Charge_[i] = 0;
        Pai4Charge_[i] = 0;
        Pai5Charge_[i] = 0;
        Pai6Charge_[i] = 0;
        Mu1Charge_[i] = 0;
        Mu2Charge_[i] = 0;
        Pai1Pt_[i] = 0;
        Pai2Pt_[i] = 0;
        Pai3Pt_[i] = 0;
        Pai4Pt_[i] = 0;
        Pai5Pt_[i] = 0;
        Pai6Pt_[i] = 0;
        Mu1Pt_[i] = 0;
        Mu2Pt_[i] = 0;
        Pai1Eta_[i] = 0;
        Pai2Eta_[i] = 0;
        Pai3Eta_[i] = 0;
        Pai4Eta_[i] = 0;
        Pai5Eta_[i] = 0;
        Pai6Eta_[i] = 0;
        Mu1Eta_[i] = 0;
        Mu2Eta_[i] = 0;
        Pai1Phi_[i] = 0;
        Pai2Phi_[i] = 0;
        Pai3Phi_[i] = 0;
        Pai4Phi_[i] = 0;
        Pai5Phi_[i] = 0;
        Pai6Phi_[i] = 0;
        Mu1Phi_[i] = 0;
        Mu2Phi_[i] = 0;
        Pai1Mass_[i] = 0;
        Pai2Mass_[i] = 0;
        Pai3Mass_[i] = 0;
        Pai4Mass_[i] = 0;
        Pai5Mass_[i] = 0;
        Pai6Mass_[i] = 0;
        Mu1Mass_[i] = 0;
        Mu2Mass_[i] = 0;
        Muon1isSoft_[i] = 0;
        Muon2isSoft_[i] = 0;
        Muon1isTight_[i] = 0;
        Muon2isTight_[i] = 0;
        Muon1isMedium_[i] = 0;
        Muon2isMedium_[i] = 0;
        Muon1isLoose_[i] = 0;
        Muon2isLoose_[i] = 0;
        Muon1isTracker_[i] = 0;
        Muon2isTracker_[i] = 0;
        Muon1isHighPt_[i] = 0;
        Muon2isHighPt_[i] = 0;
        Muon1isGlobal_[i] = 0;
        Muon2isGlobal_[i] = 0;
        Pai1Hit_[i] = 0;
        Pai2Hit_[i] = 0;
        Pai3Hit_[i] = 0;
        Pai4Hit_[i] = 0;
        Pai5Hit_[i] = 0;
        Pai6Hit_[i] = 0;
        Pai1fromPV_[i] = 0;
        Pai2fromPV_[i] = 0;
        Pai3fromPV_[i] = 0;
        Pai4fromPV_[i] = 0;
        Pai5fromPV_[i] = 0;
        Pai6fromPV_[i] = 0;
        Pai1Trknormchi2_[i] = 0;
        Pai2Trknormchi2_[i] = 0;
        Pai3Trknormchi2_[i] = 0;
        Pai4Trknormchi2_[i] = 0;
        Pai5Trknormchi2_[i] = 0;
        Pai6Trknormchi2_[i] = 0;
        J1normchi2_[i] = 0;
        J2normchi2_[i] = 0;
        J1Pt_[i] = 0;
        J2Pt_[i] = 0;
        J1Eta_[i] = 0;
        J2Eta_[i] = 0;
        J1Phi_[i] = 0;
        J2Phi_[i] = 0;
        J1Mass_[i] = 0;
        J2Mass_[i] = 0;
        XpT_[i] = 0;
        Xeta_[i] = 0;
        Xphi_[i] = 0;
        Xmass_[i] = 0;
        Xnormchi2_[i] = 0;

        J1NoMassnormchi2_[i] = 0;
        J2NoMassnormchi2_[i] = 0;
        J1NoMassPt_[i] = 0;
        J2NoMassPt_[i] = 0;
        J1NoMassEta_[i] = 0;
        J2NoMassEta_[i] = 0;
        J1NoMassPhi_[i] = 0;
        J2NoMassPhi_[i] = 0;
        J1NoMassMass_[i] = 0;
        J2NoMassMass_[i] = 0;
        J1NoMassMassE_[i] = 0;
        J2NoMassMassE_[i] = 0;
        EtaMassDiff_[i] = 0;

        minDR_[i] = 0;

        type_[i] = 0;
    }
}

void X4muPatSecondaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("tracks", edm::InputTag("tracks"));
    desc.add<edm::InputTag>("patMuon", edm::InputTag("patMuon"));
    desc.add<edm::InputTag>("vertices", edm::InputTag("vertices"));
    desc.add<edm::InputTag>("TriggerResults");
    desc.add<std::vector<std::string>>("FilterNames");
    desc.add<double>("MesonPaiMass");
    desc.add<double>("MesonPaiMassErr");
    desc.add<double>("MesonMuMass");
    desc.add<double>("MesonMuMassErr");
    desc.add<double>("MesonMassWindow");
    desc.add<double>("ExMesonMass");
    desc.add<bool>("doIso");
    desc.add<double>("MuIso");
    desc.add<double>("PaiIso");
    desc.add<bool>("doTrigger");
    desc.add<double>("vProb");
    desc.add<int>("selectionType");
    desc.add<unsigned int>("maxLoop");

    descriptions.add("X4muPatSecondaryVertexProducer", desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(X4muPatSecondaryVertexProducer);