/*
 * Developed by Yiyang Zhao for Run-3 X-->J/psiJ/psi-->4mu Scouting Analysis
 * 2024-10
 * X4muSecondaryVertexProducer loop reco::Muon
 * and fit secondary vertex for dimuon + diJ/psi
 */

#include <memory>

#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDPutToken.h"
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
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "MagneticField/Engine/interface/MagneticField.h"
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

class X4muSecondaryVertexProducer : public edm::stream::EDProducer<>
{
public:
   explicit X4muSecondaryVertexProducer(edm::ParameterSet const &iConfig);
   ~X4muSecondaryVertexProducer() override;

   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
   void beginStream(edm::StreamID) override {}
   void produce(edm::Event &iEvent, edm::EventSetup const &iSetup) override;
   void endStream() override {}

   void clearVars();

private:
   const edm::EDGetTokenT<std::vector<reco::Muon>> input_recomuon_token_;
   const edm::EDGetTokenT<std::vector<reco::Track>> input_recoTrack_token_;
   const double input_MesonMassBig_c;
   const double input_MesonMassBigErr_c;
   const double input_MesonMassSmall_c;
   const double input_MesonMassSmallErr_c;
   const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_; // 声明 magneticFieldToken_

   TTree *X4muTree;
   std::string file_name;

   /* UInt_t    run;
   ULong64_t event;
   UInt_t    lumiblock; */

   int Muon1Charge_[36];
   int Muon2Charge_[36];
   int Muon3Charge_[36];
   int Muon4Charge_[36];
   float Muon1Pt_[36];
   float Muon2Pt_[36];
   float Muon3Pt_[36];
   float Muon4Pt_[36];
   float Muon1Eta_[36];
   float Muon2Eta_[36];
   float Muon3Eta_[36];
   float Muon4Eta_[36];
   float Muon1Phi_[36];
   float Muon2Phi_[36];
   float Muon3Phi_[36];
   float Muon4Phi_[36];
   float Muon1Mass_[36];
   float Muon2Mass_[36];
   float Muon3Mass_[36];
   float Muon4Mass_[36];
   float J1normchi2_[36];
   float J2normchi2_[36];
   float J1NoMassnormchi2_[36];
   float J2NoMassnormchi2_[36];
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
   float XpT_[36];
   float Xeta_[36];
   float Xphi_[36];
   float Xmass_[36];
   float Xnormchi2_[36];
   int type_[36];
};

//
// constructors and destructor
//
X4muSecondaryVertexProducer::X4muSecondaryVertexProducer(edm::ParameterSet const &iConfig)
    : input_recomuon_token_(consumes(iConfig.getParameter<edm::InputTag>("recoMuon"))),
      input_recoTrack_token_(consumes(iConfig.getParameter<edm::InputTag>("recoTrack"))),
      input_MesonMassBig_c(iConfig.getParameter<double>("MesonMassBig")),
      input_MesonMassBigErr_c(iConfig.getParameter<double>("MesonMassBigErr")),
      input_MesonMassSmall_c(iConfig.getParameter<double>("MesonMassSmall")),
      input_MesonMassSmallErr_c(iConfig.getParameter<double>("MesonMassSmallErr")),
      magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
   edm::Service<TFileService> fs;
   X4muTree = fs->make<TTree>("X4muTree", "Tree of X4muTree");

   X4muTree->Branch("Muon1Charge", Muon1Charge_);
   X4muTree->Branch("Muon2Charge", Muon2Charge_);
   X4muTree->Branch("Muon3Charge", Muon3Charge_);
   X4muTree->Branch("Muon4Charge", Muon4Charge_);
   X4muTree->Branch("Muon1Pt", Muon1Pt_);
   X4muTree->Branch("Muon2Pt", Muon2Pt_);
   X4muTree->Branch("Muon3Pt", Muon3Pt_);
   X4muTree->Branch("Muon4Pt", Muon4Pt_);
   X4muTree->Branch("Muon1Eta", Muon1Eta_);
   X4muTree->Branch("Muon2Eta", Muon2Eta_);
   X4muTree->Branch("Muon3Eta", Muon3Eta_);
   X4muTree->Branch("Muon4Eta", Muon4Eta_);
   X4muTree->Branch("Muon1Phi", Muon1Phi_);
   X4muTree->Branch("Muon2Phi", Muon2Phi_);
   X4muTree->Branch("Muon3Phi", Muon3Phi_);
   X4muTree->Branch("Muon4Phi", Muon4Phi_);
   X4muTree->Branch("Muon1Mass", Muon1Mass_);
   X4muTree->Branch("Muon2Mass", Muon2Mass_);
   X4muTree->Branch("Muon3Mass", Muon3Mass_);
   X4muTree->Branch("Muon4Mass", Muon4Mass_);
   X4muTree->Branch("J1normchi2", J1normchi2_);
   X4muTree->Branch("J2normchi2", J2normchi2_);
   X4muTree->Branch("J1NoMassnormchi2", J1NoMassnormchi2_);
   X4muTree->Branch("J2NoMassnormchi2", J2NoMassnormchi2_);
   X4muTree->Branch("J1Pt", J1Pt_);
   X4muTree->Branch("J2Pt", J2Pt_);
   X4muTree->Branch("J1NoMassPt", J1NoMassPt_);
   X4muTree->Branch("J2NoMassPt", J2NoMassPt_);
   X4muTree->Branch("J1Eta", J1Eta_);
   X4muTree->Branch("J2Eta", J2Eta_);
   X4muTree->Branch("J1NoMassEta", J1NoMassEta_);
   X4muTree->Branch("J2NoMassEta", J2NoMassEta_);
   X4muTree->Branch("J1Phi", J1Phi_);
   X4muTree->Branch("J2Phi", J2Phi_);
   X4muTree->Branch("J1NoMassPhi", J1NoMassPhi_);
   X4muTree->Branch("J2NoMassPhi", J2NoMassPhi_);
   X4muTree->Branch("J1Mass", J1Mass_);
   X4muTree->Branch("J2Mass", J2Mass_);
   X4muTree->Branch("J1NoMassMass", J1NoMassMass_);
   X4muTree->Branch("J1NoMassMassE", J1NoMassMassE_);
   X4muTree->Branch("J2NoMassMass", J2NoMassMass_);
   X4muTree->Branch("J2NoMassMassE", J2NoMassMassE_);
   X4muTree->Branch("XpT", XpT_);
   X4muTree->Branch("Xeta", Xeta_);
   X4muTree->Branch("Xphi", Xphi_);
   X4muTree->Branch("Xmass", Xmass_);
   X4muTree->Branch("Xnormchi2", Xnormchi2_);
   X4muTree->Branch("type", type_);
}

X4muSecondaryVertexProducer::~X4muSecondaryVertexProducer() = default;

// ------------ method called to produce the data  ------------
void X4muSecondaryVertexProducer::produce(edm::Event &iEvent, edm::EventSetup const &iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_);
   const auto myYmass = input_MesonMassBig_c;
   const auto myJmass = input_MesonMassSmall_c;
   const auto myYmasserr = input_MesonMassBigErr_c;
   const auto myJmasserr = input_MesonMassSmallErr_c;

   Handle<std::vector<reco::Muon>> muonHandle;
   Handle<std::vector<reco::Track>> trackHandle;
   iEvent.getByToken(input_recomuon_token_, muonHandle);
   iEvent.getByToken(input_recoTrack_token_, trackHandle);

   double myMumass = 0.1056583755;
   double myMumasserr = myMumass * 1e-6;
   double JvPorbcut = 0.0001;
   double MassMinCut = 0.001;
   long npairs = 0;

   if (!trackHandle.isValid())
   {
      std::cout << "reco::Track collection is not valid" << std::endl;
      return;
   }

   if (trackHandle->size() < 4)
   {
      std::cout << "reco::Track size is too small: " << trackHandle->size() << std::endl;
      return;
   }

   // Loop 4 Muons
   for (std::vector<reco::Track>::const_iterator iMuon1 = trackHandle->begin(); iMuon1 != trackHandle->end(); ++iMuon1)
   {
      for (std::vector<reco::Track>::const_iterator iMuon2 = iMuon1 + 1; iMuon2 != trackHandle->end(); ++iMuon2)
      {
         for (std::vector<reco::Track>::const_iterator iMuon3 = iMuon2 + 1; iMuon3 != trackHandle->end(); ++iMuon3)
         {
            for (std::vector<reco::Track>::const_iterator iMuon4 = iMuon3 + 1; iMuon4 != trackHandle->end(); ++iMuon4)
            {
               // float mu4pairChi[3] = {}; //[0] for 12+34, [1] for 13+24, [2] for 14+23

               reco::Track muTrack1 = *iMuon1;
               reco::Track muTrack2 = *iMuon2;
               reco::Track muTrack3 = *iMuon3;
               reco::Track muTrack4 = *iMuon4;

               TransientTrack muonTT1(muTrack1, &(bFieldHandle));
               TransientTrack muonTT2(muTrack2, &(bFieldHandle));
               TransientTrack muonTT3(muTrack3, &(bFieldHandle));
               TransientTrack muonTT4(muTrack4, &(bFieldHandle));

               ParticleMass muon_mass = myMumass; // pdg mass
               float muon_sigma = myMumasserr;

               if (!muonTT1.isValid() || !muonTT2.isValid() || !muonTT3.isValid() || !muonTT4.isValid())
                  continue;

               TLorentzVector mu11P4tmp, mu12P4tmp, mu13P4tmp, mu14P4tmp, J11P4tmp, J12P4tmp, X1P4tmp;
               TLorentzVector mu21P4tmp, mu22P4tmp, mu23P4tmp, mu24P4tmp, J21P4tmp, J22P4tmp, X2P4tmp;
               TLorentzVector mu31P4tmp, mu32P4tmp, mu33P4tmp, mu34P4tmp, J31P4tmp, J32P4tmp, X3P4tmp;
               TLorentzVector J11NoMassP4tmp, J12NoMassP4tmp, J21NoMassP4tmp, J22NoMassP4tmp, J31NoMassP4tmp, J32NoMassP4tmp;
               KinematicParameters J11NoMassKtmp, J12NoMassKtmp, J21NoMassKtmp, J22NoMassKtmp, J31NoMassKtmp, J32NoMassKtmp;
               KinematicParameters mu11Ktmp, mu12Ktmp, mu13Ktmp, mu14Ktmp, J11Ktmp, J12Ktmp, X1Ktmp;
               KinematicParameters mu21Ktmp, mu22Ktmp, mu23Ktmp, mu24Ktmp, J21Ktmp, J22Ktmp, X2Ktmp;
               KinematicParameters mu31Ktmp, mu32Ktmp, mu33Ktmp, mu34Ktmp, J31Ktmp, J32Ktmp, X3Ktmp;
               float Muon11Mass, Muon12Mass, Muon13Mass, Muon14Mass;
               float Muon21Mass, Muon22Mass, Muon23Mass, Muon24Mass;
               float Muon31Mass, Muon32Mass, Muon33Mass, Muon34Mass;
               float J11Mass, J12Mass, J21Mass, J22Mass, J31Mass, J32Mass;
               float J11NoMassMass, J12NoMassMass, J21NoMassMass, J22NoMassMass, J31NoMassMass, J32NoMassMass;
               float J11NoMassMassE, J12NoMassMassE, J21NoMassMassE, J22NoMassMassE, J31NoMassMassE, J32NoMassMassE;
               float X1Mass, X2Mass, X3Mass;
               float J11normchi2 = 999.0;
               float J11NoMassnormchi2 = 999.0;
               float J12normchi2 = 999.0;
               float J12NoMassnormchi2 = 999.0;
               float J21normchi2 = 999.0;
               float J21NoMassnormchi2 = 999.0;
               float J22normchi2 = 999.0;
               float J22NoMassnormchi2 = 999.0;
               float J31normchi2 = 999.0;
               float J31NoMassnormchi2 = 999.0;
               float J32normchi2 = 999.0;
               float J32NoMassnormchi2 = 999.0;
               float X1normchi2 = 999.0;
               float X2normchi2 = 999.0;
               float X3normchi2 = 999.0;
               int Y1pos = -1;
               int Y2pos = -1;
               int Y3pos = -1;

               // 12+34
               if ((iMuon1->charge() + iMuon2->charge()) == 0 && (iMuon3->charge() + iMuon4->charge()) == 0)
               {
                  KinematicParticleFactoryFromTransientTrack pmumuFactory;

                  // initial chi2 and ndf before kinematic fits.
                  float chi1 = 0.;
                  float chi2 = 0.;
                  float ndf1 = 0.;
                  float ndf2 = 0.;
                  float chi = 0.;
                  float ndf = 0.;

                  vector<RefCountedKinematicParticle> dimuon1Particles, dimuon2Particles, mu4Particles;
                  dimuon1Particles.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi1, ndf1, muon_sigma));
                  dimuon1Particles.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi1, ndf1, muon_sigma));
                  dimuon2Particles.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi2, ndf2, muon_sigma));
                  dimuon2Particles.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi2, ndf2, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi, ndf, muon_sigma));

                  if (dimuon1Particles.size() < 2 || dimuon2Particles.size() < 2 || mu4Particles.size() < 4)
                     goto match1;

                  KinematicParticleVertexFitter fitter1, fitter2, mu4fitter;
                  RefCountedKinematicTree psiVertexFitTree1, psiVertexFitTree2, XVertexFitTree;
                  psiVertexFitTree1 = fitter1.fit(dimuon1Particles);
                  psiVertexFitTree2 = fitter2.fit(dimuon2Particles);
                  XVertexFitTree = mu4fitter.fit(mu4Particles);

                  if (psiVertexFitTree1->isValid() && psiVertexFitTree2->isValid() && XVertexFitTree->isValid())
                  {
                     psiVertexFitTree1->movePointerToTheTop();
                     psiVertexFitTree2->movePointerToTheTop();
                     XVertexFitTree->movePointerToTheTop();

                     RefCountedKinematicParticle psi_vFit1 = psiVertexFitTree1->currentParticle();
                     RefCountedKinematicParticle psi_vFit2 = psiVertexFitTree2->currentParticle();
                     RefCountedKinematicParticle X_vFit = XVertexFitTree->currentParticle();

                     RefCountedKinematicVertex psi_vFit_vertex1 = psiVertexFitTree1->currentDecayVertex();
                     RefCountedKinematicVertex psi_vFit_vertex2 = psiVertexFitTree2->currentDecayVertex();
                     RefCountedKinematicVertex X_vFit_vertex = XVertexFitTree->currentDecayVertex();

                     if (!psi_vFit1->currentState().isValid() || !psi_vFit2->currentState().isValid() || !X_vFit->currentState().isValid())
                        goto match1;
                     if (!psi_vFit_vertex1->vertexIsValid() || !psi_vFit_vertex2->vertexIsValid() || !X_vFit_vertex->vertexIsValid())
                        goto match1;

                     // KinematicParameters Jpara1 = psi_vFit1->currentState().kinematicParameters();
                     // KinematicParameters Jpara2 = psi_vFit2->currentState().kinematicParameters();
                     // KinematicParameters Xpara = X_vFit->currentState().kinematicParameters();

                     XVertexFitTree->movePointerToTheFirstChild();
                     RefCountedKinematicParticle mu1 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu2 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu3 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu4 = XVertexFitTree->currentParticle();

                     double vProb1 = ChiSquaredProbability((double)(psi_vFit_vertex1->chiSquared()), (double)(psi_vFit_vertex1->degreesOfFreedom()));
                     double vProb2 = ChiSquaredProbability((double)(psi_vFit_vertex2->chiSquared()), (double)(psi_vFit_vertex2->degreesOfFreedom()));
                     // double XvProb = ChiSquaredProbability((double)(X_vFit_vertex->chiSquared()), (double)(X_vFit_vertex->degreesOfFreedom()));
                     if (vProb1 < JvPorbcut || vProb2 < JvPorbcut)
                        goto match1;
                     if (mu1->currentState().mass() <= MassMinCut || mu2->currentState().mass() <= MassMinCut || mu3->currentState().mass() <= MassMinCut || mu4->currentState().mass() <= MassMinCut || psi_vFit1->currentState().mass() <= MassMinCut || psi_vFit2->currentState().mass() <= MassMinCut || X_vFit->currentState().mass() <= MassMinCut)
                        goto match1;
                     if (psi_vFit1->currentState().kinematicParametersError().matrix()(6, 6) < 0 || psi_vFit2->currentState().kinematicParametersError().matrix()(6, 6) < 0 || X_vFit->currentState().kinematicParametersError().matrix()(6, 6) < 0)
                        goto match1;
                     float Jpsi1masserr = sqrt(psi_vFit1->currentState().kinematicParametersError().matrix()(6, 6));
                     float Jpsi2masserr = sqrt(psi_vFit2->currentState().kinematicParametersError().matrix()(6, 6));
                     if (psi_vFit1->currentState().mass() > (myYmass - 3.0 * Jpsi1masserr) && psi_vFit1->currentState().mass() < (myYmass + 3.0 * Jpsi1masserr) && psi_vFit2->currentState().mass() > (myJmass - 3.0 * Jpsi2masserr) && psi_vFit2->currentState().mass() < (myJmass + 3.0 * Jpsi2masserr))
                     {
                        Y1pos = 1;
                        J11NoMassKtmp = psi_vFit1->currentState().kinematicParameters();
                        J12NoMassKtmp = psi_vFit2->currentState().kinematicParameters();
                        J11NoMassP4tmp.SetPxPyPzE(J11NoMassKtmp.momentum().x(), J11NoMassKtmp.momentum().y(), J11NoMassKtmp.momentum().z(), psi_vFit1->currentState().mass());
                        J12NoMassP4tmp.SetPxPyPzE(J12NoMassKtmp.momentum().x(), J12NoMassKtmp.momentum().y(), J12NoMassKtmp.momentum().z(), psi_vFit2->currentState().mass());
                        J11NoMassMass = psi_vFit1->currentState().mass();
                        J12NoMassMass = psi_vFit2->currentState().mass();
                        J11NoMassMassE = Jpsi1masserr;
                        J12NoMassMassE = Jpsi2masserr;
                     }
                     if (psi_vFit2->currentState().mass() > (myYmass - 3.0 * Jpsi2masserr) && psi_vFit2->currentState().mass() < (myYmass + 3.0 * Jpsi2masserr) && psi_vFit1->currentState().mass() > (myJmass - 3.0 * Jpsi1masserr) && psi_vFit1->currentState().mass() < (myJmass + 3.0 * Jpsi1masserr) && vProb2 > vProb1)
                     {
                        Y1pos = 2;
                        J11NoMassKtmp = psi_vFit2->currentState().kinematicParameters();
                        J12NoMassKtmp = psi_vFit1->currentState().kinematicParameters();
                        J11NoMassP4tmp.SetPxPyPzE(J11NoMassKtmp.momentum().x(), J11NoMassKtmp.momentum().y(), J11NoMassKtmp.momentum().z(), psi_vFit2->currentState().mass());
                        J12NoMassP4tmp.SetPxPyPzE(J12NoMassKtmp.momentum().x(), J12NoMassKtmp.momentum().y(), J12NoMassKtmp.momentum().z(), psi_vFit1->currentState().mass());
                        J11NoMassMass = psi_vFit2->currentState().mass();
                        J12NoMassMass = psi_vFit1->currentState().mass();
                        J11NoMassMassE = Jpsi2masserr;
                        J12NoMassMassE = Jpsi1masserr;
                     }
                     else
                        goto match1;

                     // Mass Constraint Fit for J/psi
                     KinematicParticleVertexFitter kpvFitter;
                     KinematicParticleFitter csFitter;
                     ParticleMass jp1, jp2;
                     float jp_m_sigma1, jp_m_sigma2;
                     KinematicConstraint *jpsi_c1;
                     KinematicConstraint *jpsi_c2;
                     vector<RefCountedKinematicParticle> muonP12, muonP34;
                     RefCountedKinematicParticle Jpsi1_part, Jpsi2_part;
                     vector<RefCountedKinematicParticle> Chi_1;
                     muonP12.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi1, ndf1, muon_sigma));
                     muonP12.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi1, ndf1, muon_sigma));
                     muonP34.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi2, ndf2, muon_sigma));
                     muonP34.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi2, ndf2, muon_sigma));
                     RefCountedKinematicTree Jpsi1 = kpvFitter.fit(muonP12);
                     RefCountedKinematicTree Jpsi2 = kpvFitter.fit(muonP34);
                     RefCountedKinematicTree Jpsi1noMCJJ = kpvFitter.fit(muonP12);
                     RefCountedKinematicTree Jpsi2noMCJJ = kpvFitter.fit(muonP34);
                     RefCountedKinematicTree Chi1_bTree;
                     RefCountedKinematicParticle MyChi1_part;
                     if (Y1pos == 1)
                     {
                        jp1 = myYmass;
                        jp_m_sigma1 = myYmasserr;
                        jp2 = myJmass;
                        jp_m_sigma2 = myJmasserr;
                        jpsi_c1 = new MassKinematicConstraint(jp1, jp_m_sigma1);
                        jpsi_c2 = new MassKinematicConstraint(jp2, jp_m_sigma2);
                        try
                        {
                           Jpsi1 = csFitter.fit(jpsi_c1, Jpsi1noMCJJ);
                           Jpsi2 = csFitter.fit(jpsi_c2, Jpsi2noMCJJ);
                        }
                        catch (VertexException const &x)
                        {
                           std::cout << "mu12 vertex exception with mass constrainted to J!" << std::endl;
                        }
                        if (Jpsi1->isEmpty() != true && Jpsi2->isValid() == true)
                        {
                           Jpsi1->movePointerToTheTop();
                           Jpsi2->movePointerToTheTop();
                           Jpsi1_part = Jpsi1->currentParticle();
                           Jpsi2_part = Jpsi2->currentParticle();
                           Chi_1.push_back(Jpsi1_part);
                           Chi_1.push_back(Jpsi2_part);
                           bool isagoodfit = true;
                           try
                           {
                              Chi1_bTree = kpvFitter.fit(Chi_1);
                           }
                           catch (VertexException const &x)
                           {
                              isagoodfit = false;
                              cout << "mu12 and mu34 vertex exception with mu12 and mu34 constrained to JJ" << endl;
                           }
                           if (Chi1_bTree->isValid() && isagoodfit)
                           {
                              Chi1_bTree->movePointerToTheTop();
                              RefCountedKinematicVertex myJpsi1Jpsi2Vertex = Chi1_bTree->currentDecayVertex();
                              Chi1_bTree->movePointerToTheFirstChild();
                              RefCountedKinematicParticle jpsi1In2JpsiVertex = Chi1_bTree->currentParticle();
                              Chi1_bTree->movePointerToTheNextChild();
                              RefCountedKinematicParticle jpsi2In2JpsiVertex = Chi1_bTree->currentParticle();

                              RefCountedKinematicVertex myJpsi1DecayVtx = Jpsi1->currentDecayVertex();
                              RefCountedKinematicVertex myJpsi2DecayVtx = Jpsi2->currentDecayVertex();
                              mu11Ktmp = mu1->currentState().kinematicParameters();
                              mu12Ktmp = mu2->currentState().kinematicParameters();
                              mu13Ktmp = mu3->currentState().kinematicParameters();
                              mu14Ktmp = mu4->currentState().kinematicParameters();
                              J11Ktmp = jpsi1In2JpsiVertex->currentState().kinematicParameters();
                              J12Ktmp = jpsi2In2JpsiVertex->currentState().kinematicParameters();
                              Chi1_bTree->movePointerToTheTop();
                              MyChi1_part = Chi1_bTree->currentParticle();
                              X1Ktmp = MyChi1_part->currentState().kinematicParameters();
                              mu11P4tmp.SetPxPyPzE(mu11Ktmp.momentum().x(), mu11Ktmp.momentum().y(), mu11Ktmp.momentum().z(), mu1->currentState().mass());
                              mu12P4tmp.SetPxPyPzE(mu12Ktmp.momentum().x(), mu12Ktmp.momentum().y(), mu12Ktmp.momentum().z(), mu2->currentState().mass());
                              mu13P4tmp.SetPxPyPzE(mu13Ktmp.momentum().x(), mu13Ktmp.momentum().y(), mu13Ktmp.momentum().z(), mu3->currentState().mass());
                              mu14P4tmp.SetPxPyPzE(mu14Ktmp.momentum().x(), mu14Ktmp.momentum().y(), mu14Ktmp.momentum().z(), mu4->currentState().mass());
                              J11P4tmp.SetPxPyPzE(J11Ktmp.momentum().x(), J11Ktmp.momentum().y(), J11Ktmp.momentum().z(), jpsi1In2JpsiVertex->currentState().mass());
                              J12P4tmp.SetPxPyPzE(J12Ktmp.momentum().x(), J12Ktmp.momentum().y(), J12Ktmp.momentum().z(), jpsi2In2JpsiVertex->currentState().mass());
                              X1P4tmp.SetPxPyPzE(X1Ktmp.momentum().x(), X1Ktmp.momentum().y(), X1Ktmp.momentum().z(), MyChi1_part->currentState().mass());
                              Muon11Mass = mu1->currentState().mass();
                              Muon12Mass = mu2->currentState().mass();
                              Muon13Mass = mu3->currentState().mass();
                              Muon14Mass = mu4->currentState().mass();
                              J11normchi2 = myJpsi1DecayVtx->chiSquared() / myJpsi1DecayVtx->degreesOfFreedom();
                              J12normchi2 = myJpsi2DecayVtx->chiSquared() / myJpsi2DecayVtx->degreesOfFreedom();
                              J11Mass = jpsi1In2JpsiVertex->currentState().mass();
                              J12Mass = jpsi2In2JpsiVertex->currentState().mass();
                              X1Mass = MyChi1_part->currentState().mass();
                              X1normchi2 = (double)myJpsi1Jpsi2Vertex->chiSquared() / (double)myJpsi1Jpsi2Vertex->degreesOfFreedom();
                           }
                        }
                     }
                     else if (Y1pos == 2)
                     {
                        jp2 = myYmass;
                        jp_m_sigma2 = myYmasserr;
                        jp1 = myJmass;
                        jp_m_sigma1 = myJmasserr;
                        jpsi_c2 = new MassKinematicConstraint(jp2, jp_m_sigma2);
                        jpsi_c1 = new MassKinematicConstraint(jp1, jp_m_sigma1);
                        try
                        {
                           Jpsi2 = csFitter.fit(jpsi_c2, Jpsi2noMCJJ);
                           Jpsi1 = csFitter.fit(jpsi_c1, Jpsi1noMCJJ);
                        }
                        catch (VertexException const &x)
                        {
                           std::cout << "mu34 vertex exception with mass constrainted to J!" << std::endl;
                        }
                        if (Jpsi2->isEmpty() != true && Jpsi1->isValid() == true)
                        {
                           Jpsi2->movePointerToTheTop();
                           Jpsi1->movePointerToTheTop();
                           Jpsi1_part = Jpsi2->currentParticle();
                           Jpsi2_part = Jpsi1->currentParticle();
                           Chi_1.push_back(Jpsi1_part);
                           Chi_1.push_back(Jpsi2_part);
                           bool isagoodfit = true;
                           try
                           {
                              Chi1_bTree = kpvFitter.fit(Chi_1);
                           }
                           catch (VertexException const &x)
                           {
                              isagoodfit = false;
                              cout << "mu12 and mu34 vertex exception with mu12 and mu34 constrained to JJ" << endl;
                           }
                           if (Chi1_bTree->isValid() && isagoodfit)
                           {
                              Chi1_bTree->movePointerToTheTop();
                              RefCountedKinematicVertex myJpsi1Jpsi2Vertex = Chi1_bTree->currentDecayVertex();
                              Chi1_bTree->movePointerToTheFirstChild();
                              RefCountedKinematicParticle jpsi1In2JpsiVertex = Chi1_bTree->currentParticle();
                              Chi1_bTree->movePointerToTheNextChild();
                              RefCountedKinematicParticle jpsi2In2JpsiVertex = Chi1_bTree->currentParticle();

                              RefCountedKinematicVertex myJpsi1DecayVtx = Jpsi2->currentDecayVertex();
                              RefCountedKinematicVertex myJpsi2DecayVtx = Jpsi1->currentDecayVertex();
                              mu11Ktmp = mu1->currentState().kinematicParameters();
                              mu12Ktmp = mu2->currentState().kinematicParameters();
                              mu13Ktmp = mu3->currentState().kinematicParameters();
                              mu14Ktmp = mu4->currentState().kinematicParameters();
                              J11Ktmp = jpsi1In2JpsiVertex->currentState().kinematicParameters();
                              J12Ktmp = jpsi2In2JpsiVertex->currentState().kinematicParameters();
                              Chi1_bTree->movePointerToTheTop();
                              MyChi1_part = Chi1_bTree->currentParticle();
                              X1Ktmp = MyChi1_part->currentState().kinematicParameters();
                              mu11P4tmp.SetPxPyPzE(mu11Ktmp.momentum().x(), mu11Ktmp.momentum().y(), mu11Ktmp.momentum().z(), mu1->currentState().mass());
                              mu12P4tmp.SetPxPyPzE(mu12Ktmp.momentum().x(), mu12Ktmp.momentum().y(), mu12Ktmp.momentum().z(), mu2->currentState().mass());
                              mu13P4tmp.SetPxPyPzE(mu13Ktmp.momentum().x(), mu13Ktmp.momentum().y(), mu13Ktmp.momentum().z(), mu3->currentState().mass());
                              mu14P4tmp.SetPxPyPzE(mu14Ktmp.momentum().x(), mu14Ktmp.momentum().y(), mu14Ktmp.momentum().z(), mu4->currentState().mass());
                              J11P4tmp.SetPxPyPzE(J11Ktmp.momentum().x(), J11Ktmp.momentum().y(), J11Ktmp.momentum().z(), jpsi1In2JpsiVertex->currentState().mass());
                              J12P4tmp.SetPxPyPzE(J12Ktmp.momentum().x(), J12Ktmp.momentum().y(), J12Ktmp.momentum().z(), jpsi2In2JpsiVertex->currentState().mass());
                              X1P4tmp.SetPxPyPzE(X1Ktmp.momentum().x(), X1Ktmp.momentum().y(), X1Ktmp.momentum().z(), MyChi1_part->currentState().mass());
                              Muon11Mass = mu1->currentState().mass();
                              Muon12Mass = mu2->currentState().mass();
                              Muon13Mass = mu3->currentState().mass();
                              Muon14Mass = mu4->currentState().mass();
                              J11normchi2 = myJpsi1DecayVtx->chiSquared() / myJpsi1DecayVtx->degreesOfFreedom();
                              J12normchi2 = myJpsi2DecayVtx->chiSquared() / myJpsi2DecayVtx->degreesOfFreedom();
                              J11Mass = jpsi1In2JpsiVertex->currentState().mass();
                              J12Mass = jpsi2In2JpsiVertex->currentState().mass();
                              X1Mass = MyChi1_part->currentState().mass();
                              X1normchi2 = (double)myJpsi1Jpsi2Vertex->chiSquared() / (double)myJpsi1Jpsi2Vertex->degreesOfFreedom();
                           }
                        }
                     }
                  }
               }
match1:
               // 13+24
               if ((iMuon1->charge() + iMuon3->charge()) == 0 && (iMuon2->charge() + iMuon4->charge()) == 0)
               {
                  KinematicParticleFactoryFromTransientTrack pmumuFactory;

                  // initial chi2 and ndf before kinematic fits.
                  float chi1 = 0.;
                  float chi2 = 0.;
                  float ndf1 = 0.;
                  float ndf2 = 0.;
                  float chi = 0.;
                  float ndf = 0.;

                  vector<RefCountedKinematicParticle> dimuon1Particles, dimuon2Particles, mu4Particles;
                  dimuon1Particles.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi1, ndf1, muon_sigma));
                  dimuon1Particles.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi1, ndf1, muon_sigma));
                  dimuon2Particles.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi2, ndf2, muon_sigma));
                  dimuon2Particles.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi2, ndf2, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi, ndf, muon_sigma));

                  if (dimuon1Particles.size() < 2 || dimuon2Particles.size() < 2 || mu4Particles.size() < 4)
                     goto match2;

                  KinematicParticleVertexFitter fitter1, fitter2, mu4fitter;
                  RefCountedKinematicTree psiVertexFitTree1, psiVertexFitTree2, XVertexFitTree;
                  psiVertexFitTree1 = fitter1.fit(dimuon1Particles);
                  psiVertexFitTree2 = fitter2.fit(dimuon2Particles);
                  XVertexFitTree = mu4fitter.fit(mu4Particles);

                  if (psiVertexFitTree1->isValid() && psiVertexFitTree2->isValid() && XVertexFitTree->isValid())
                  {
                     psiVertexFitTree1->movePointerToTheTop();
                     psiVertexFitTree2->movePointerToTheTop();
                     XVertexFitTree->movePointerToTheTop();

                     RefCountedKinematicParticle psi_vFit1 = psiVertexFitTree1->currentParticle();
                     RefCountedKinematicParticle psi_vFit2 = psiVertexFitTree2->currentParticle();
                     RefCountedKinematicParticle X_vFit = XVertexFitTree->currentParticle();

                     RefCountedKinematicVertex psi_vFit_vertex1 = psiVertexFitTree1->currentDecayVertex();
                     RefCountedKinematicVertex psi_vFit_vertex2 = psiVertexFitTree2->currentDecayVertex();
                     RefCountedKinematicVertex X_vFit_vertex = XVertexFitTree->currentDecayVertex();

                     if (!psi_vFit1->currentState().isValid() || !psi_vFit2->currentState().isValid() || !X_vFit->currentState().isValid())
                        goto match2;
                     if (!psi_vFit_vertex1->vertexIsValid() || !psi_vFit_vertex2->vertexIsValid() || !X_vFit_vertex->vertexIsValid())
                        goto match2;

                     // KinematicParameters Jpara1 = psi_vFit1->currentState().kinematicParameters();
                     // KinematicParameters Jpara2 = psi_vFit2->currentState().kinematicParameters();
                     // KinematicParameters Xpara = X_vFit->currentState().kinematicParameters();

                     XVertexFitTree->movePointerToTheFirstChild();
                     RefCountedKinematicParticle mu1 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu3 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu2 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu4 = XVertexFitTree->currentParticle();

                     double vProb1 = ChiSquaredProbability((double)(psi_vFit_vertex1->chiSquared()), (double)(psi_vFit_vertex1->degreesOfFreedom()));
                     double vProb2 = ChiSquaredProbability((double)(psi_vFit_vertex2->chiSquared()), (double)(psi_vFit_vertex2->degreesOfFreedom()));
                     // double XvProb = ChiSquaredProbability((double)(X_vFit_vertex->chiSquared()), (double)(X_vFit_vertex->degreesOfFreedom()));
                     if (vProb1 < JvPorbcut || vProb2 < JvPorbcut)
                        goto match2;
                     if (mu1->currentState().mass() <= MassMinCut || mu2->currentState().mass() <= MassMinCut || mu3->currentState().mass() <= MassMinCut || mu4->currentState().mass() <= MassMinCut || psi_vFit1->currentState().mass() <= MassMinCut || psi_vFit2->currentState().mass() <= MassMinCut || X_vFit->currentState().mass() <= MassMinCut)
                        goto match2;
                     if (psi_vFit1->currentState().kinematicParametersError().matrix()(6, 6) < 0 || psi_vFit2->currentState().kinematicParametersError().matrix()(6, 6) < 0 || X_vFit->currentState().kinematicParametersError().matrix()(6, 6) < 0)
                        goto match2;
                     float Jpsi1masserr = sqrt(psi_vFit1->currentState().kinematicParametersError().matrix()(6, 6));
                     float Jpsi2masserr = sqrt(psi_vFit2->currentState().kinematicParametersError().matrix()(6, 6));
                     if (psi_vFit1->currentState().mass() > (myYmass - 3.0 * Jpsi1masserr) && psi_vFit1->currentState().mass() < (myYmass + 3.0 * Jpsi1masserr) && psi_vFit2->currentState().mass() > (myJmass - 3.0 * Jpsi2masserr) && psi_vFit2->currentState().mass() < (myJmass + 3.0 * Jpsi2masserr))
                     {
                        Y2pos = 1;
                        J21NoMassKtmp = psi_vFit1->currentState().kinematicParameters();
                        J22NoMassKtmp = psi_vFit2->currentState().kinematicParameters();
                        J21NoMassP4tmp.SetPxPyPzE(J21NoMassKtmp.momentum().x(), J21NoMassKtmp.momentum().y(), J21NoMassKtmp.momentum().z(), psi_vFit1->currentState().mass());
                        J22NoMassP4tmp.SetPxPyPzE(J22NoMassKtmp.momentum().x(), J22NoMassKtmp.momentum().y(), J22NoMassKtmp.momentum().z(), psi_vFit2->currentState().mass());
                        J21NoMassMass = psi_vFit1->currentState().mass();
                        J22NoMassMass = psi_vFit2->currentState().mass();
                        J21NoMassMassE = Jpsi1masserr;
                        J22NoMassMassE = Jpsi2masserr;
                     }
                     if (psi_vFit2->currentState().mass() > (myYmass - 3.0 * Jpsi2masserr) && psi_vFit2->currentState().mass() < (myYmass + 3.0 * Jpsi2masserr) && psi_vFit1->currentState().mass() > (myJmass - 3.0 * Jpsi1masserr) && psi_vFit1->currentState().mass() < (myJmass + 3.0 * Jpsi1masserr) && vProb2 > vProb1)
                     {
                        Y2pos = 2;
                        J21NoMassKtmp = psi_vFit2->currentState().kinematicParameters();
                        J22NoMassKtmp = psi_vFit1->currentState().kinematicParameters();
                        J21NoMassP4tmp.SetPxPyPzE(J21NoMassKtmp.momentum().x(), J21NoMassKtmp.momentum().y(), J21NoMassKtmp.momentum().z(), psi_vFit2->currentState().mass());
                        J22NoMassP4tmp.SetPxPyPzE(J22NoMassKtmp.momentum().x(), J22NoMassKtmp.momentum().y(), J22NoMassKtmp.momentum().z(), psi_vFit1->currentState().mass());
                        J21NoMassMass = psi_vFit2->currentState().mass();
                        J22NoMassMass = psi_vFit1->currentState().mass();
                        J21NoMassMassE = Jpsi2masserr;
                        J22NoMassMassE = Jpsi1masserr;
                     }
                     else
                        goto match2;

                     // Mass Constraint Fit for J/psi
                     KinematicParticleVertexFitter kpvFitter;
                     KinematicParticleFitter csFitter;
                     ParticleMass jp1, jp2;
                     float jp_m_sigma1, jp_m_sigma2;
                     KinematicConstraint *jpsi_c1;
                     KinematicConstraint *jpsi_c2;
                     vector<RefCountedKinematicParticle> muonP12, muonP34;
                     RefCountedKinematicParticle Jpsi1_part, Jpsi2_part;
                     vector<RefCountedKinematicParticle> Chi_1;
                     muonP12.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi1, ndf1, muon_sigma));
                     muonP12.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi1, ndf1, muon_sigma));
                     muonP34.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi2, ndf2, muon_sigma));
                     muonP34.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi2, ndf2, muon_sigma));
                     RefCountedKinematicTree Jpsi1 = kpvFitter.fit(muonP12);
                     RefCountedKinematicTree Jpsi2 = kpvFitter.fit(muonP34);
                     RefCountedKinematicTree Jpsi1noMCJJ = kpvFitter.fit(muonP12);
                     RefCountedKinematicTree Jpsi2noMCJJ = kpvFitter.fit(muonP34);
                     RefCountedKinematicTree Chi1_bTree;
                     RefCountedKinematicParticle MyChi1_part;
                     if (Y2pos == 1)
                     {
                        jp1 = myYmass;
                        jp_m_sigma1 = myYmasserr;
                        jp2 = myJmass;
                        jp_m_sigma2 = myJmasserr;
                        jpsi_c1 = new MassKinematicConstraint(jp1, jp_m_sigma1);
                        jpsi_c2 = new MassKinematicConstraint(jp2, jp_m_sigma2);
                        try
                        {
                           Jpsi1 = csFitter.fit(jpsi_c1, Jpsi1noMCJJ);
                           Jpsi2 = csFitter.fit(jpsi_c2, Jpsi2noMCJJ);
                        }
                        catch (VertexException const &x)
                        {
                           std::cout << "mu12 vertex exception with mass constrainted to J!" << std::endl;
                        }
                        if (Jpsi1->isEmpty() != true && Jpsi2->isValid() == true)
                        {
                           Jpsi1->movePointerToTheTop();
                           Jpsi2->movePointerToTheTop();
                           Jpsi1_part = Jpsi1->currentParticle();
                           Jpsi2_part = Jpsi2->currentParticle();
                           Chi_1.push_back(Jpsi1_part);
                           Chi_1.push_back(Jpsi2_part);
                           bool isagoodfit = true;
                           try
                           {
                              Chi1_bTree = kpvFitter.fit(Chi_1);
                           }
                           catch (VertexException const &x)
                           {
                              isagoodfit = false;
                              cout << "mu12 and mu34 vertex exception with mu12 and mu34 constrained to JJ" << endl;
                           }
                           if (Chi1_bTree->isValid() && isagoodfit)
                           {
                              Chi1_bTree->movePointerToTheTop();
                              RefCountedKinematicVertex myJpsi1Jpsi2Vertex = Chi1_bTree->currentDecayVertex();
                              Chi1_bTree->movePointerToTheFirstChild();
                              RefCountedKinematicParticle jpsi1In2JpsiVertex = Chi1_bTree->currentParticle();
                              Chi1_bTree->movePointerToTheNextChild();
                              RefCountedKinematicParticle jpsi2In2JpsiVertex = Chi1_bTree->currentParticle();

                              RefCountedKinematicVertex myJpsi1DecayVtx = Jpsi1->currentDecayVertex();
                              RefCountedKinematicVertex myJpsi2DecayVtx = Jpsi2->currentDecayVertex();
                              mu21Ktmp = mu1->currentState().kinematicParameters();
                              mu22Ktmp = mu2->currentState().kinematicParameters();
                              mu23Ktmp = mu3->currentState().kinematicParameters();
                              mu24Ktmp = mu4->currentState().kinematicParameters();
                              J21Ktmp = jpsi1In2JpsiVertex->currentState().kinematicParameters();
                              J22Ktmp = jpsi2In2JpsiVertex->currentState().kinematicParameters();
                              Chi1_bTree->movePointerToTheTop();
                              MyChi1_part = Chi1_bTree->currentParticle();
                              X2Ktmp = MyChi1_part->currentState().kinematicParameters();
                              mu21P4tmp.SetPxPyPzE(mu21Ktmp.momentum().x(), mu21Ktmp.momentum().y(), mu21Ktmp.momentum().z(), mu1->currentState().mass());
                              mu22P4tmp.SetPxPyPzE(mu22Ktmp.momentum().x(), mu22Ktmp.momentum().y(), mu22Ktmp.momentum().z(), mu2->currentState().mass());
                              mu23P4tmp.SetPxPyPzE(mu23Ktmp.momentum().x(), mu23Ktmp.momentum().y(), mu23Ktmp.momentum().z(), mu3->currentState().mass());
                              mu24P4tmp.SetPxPyPzE(mu24Ktmp.momentum().x(), mu24Ktmp.momentum().y(), mu24Ktmp.momentum().z(), mu4->currentState().mass());
                              J21P4tmp.SetPxPyPzE(J21Ktmp.momentum().x(), J21Ktmp.momentum().y(), J21Ktmp.momentum().z(), jpsi1In2JpsiVertex->currentState().mass());
                              J22P4tmp.SetPxPyPzE(J22Ktmp.momentum().x(), J22Ktmp.momentum().y(), J22Ktmp.momentum().z(), jpsi2In2JpsiVertex->currentState().mass());
                              X2P4tmp.SetPxPyPzE(X2Ktmp.momentum().x(), X2Ktmp.momentum().y(), X2Ktmp.momentum().z(), MyChi1_part->currentState().mass());
                              Muon21Mass = mu1->currentState().mass();
                              Muon22Mass = mu2->currentState().mass();
                              Muon23Mass = mu3->currentState().mass();
                              Muon24Mass = mu4->currentState().mass();
                              J21normchi2 = myJpsi1DecayVtx->chiSquared() / myJpsi1DecayVtx->degreesOfFreedom();
                              J22normchi2 = myJpsi2DecayVtx->chiSquared() / myJpsi2DecayVtx->degreesOfFreedom();
                              J21Mass = jpsi1In2JpsiVertex->currentState().mass();
                              J22Mass = jpsi2In2JpsiVertex->currentState().mass();
                              X2Mass = MyChi1_part->currentState().mass();
                              X2normchi2 = (double)myJpsi1Jpsi2Vertex->chiSquared() / (double)myJpsi1Jpsi2Vertex->degreesOfFreedom();
                           }
                        }
                     }
                     else if (Y2pos == 2)
                     {
                        jp2 = myYmass;
                        jp_m_sigma2 = myYmasserr;
                        jp1 = myJmass;
                        jp_m_sigma1 = myJmasserr;
                        jpsi_c2 = new MassKinematicConstraint(jp2, jp_m_sigma2);
                        jpsi_c1 = new MassKinematicConstraint(jp1, jp_m_sigma1);
                        try
                        {
                           Jpsi2 = csFitter.fit(jpsi_c2, Jpsi2noMCJJ);
                           Jpsi1 = csFitter.fit(jpsi_c1, Jpsi1noMCJJ);
                        }
                        catch (VertexException const &x)
                        {
                           std::cout << "mu34 vertex exception with mass constrainted to J!" << std::endl;
                        }
                        if (Jpsi2->isEmpty() != true && Jpsi1->isValid() == true)
                        {
                           Jpsi2->movePointerToTheTop();
                           Jpsi1->movePointerToTheTop();
                           Jpsi1_part = Jpsi2->currentParticle();
                           Jpsi2_part = Jpsi1->currentParticle();
                           Chi_1.push_back(Jpsi1_part);
                           Chi_1.push_back(Jpsi2_part);
                           bool isagoodfit = true;
                           try
                           {
                              Chi1_bTree = kpvFitter.fit(Chi_1);
                           }
                           catch (VertexException const &x)
                           {
                              isagoodfit = false;
                              cout << "mu12 and mu34 vertex exception with mu12 and mu34 constrained to JJ" << endl;
                           }
                           if (Chi1_bTree->isValid() && isagoodfit)
                           {
                              Chi1_bTree->movePointerToTheTop();
                              RefCountedKinematicVertex myJpsi1Jpsi2Vertex = Chi1_bTree->currentDecayVertex();
                              Chi1_bTree->movePointerToTheFirstChild();
                              RefCountedKinematicParticle jpsi1In2JpsiVertex = Chi1_bTree->currentParticle();
                              Chi1_bTree->movePointerToTheNextChild();
                              RefCountedKinematicParticle jpsi2In2JpsiVertex = Chi1_bTree->currentParticle();

                              RefCountedKinematicVertex myJpsi1DecayVtx = Jpsi2->currentDecayVertex();
                              RefCountedKinematicVertex myJpsi2DecayVtx = Jpsi1->currentDecayVertex();
                              mu21Ktmp = mu1->currentState().kinematicParameters();
                              mu22Ktmp = mu2->currentState().kinematicParameters();
                              mu23Ktmp = mu3->currentState().kinematicParameters();
                              mu24Ktmp = mu4->currentState().kinematicParameters();
                              J21Ktmp = jpsi1In2JpsiVertex->currentState().kinematicParameters();
                              J22Ktmp = jpsi2In2JpsiVertex->currentState().kinematicParameters();
                              Chi1_bTree->movePointerToTheTop();
                              MyChi1_part = Chi1_bTree->currentParticle();
                              X2Ktmp = MyChi1_part->currentState().kinematicParameters();
                              mu21P4tmp.SetPxPyPzE(mu21Ktmp.momentum().x(), mu21Ktmp.momentum().y(), mu21Ktmp.momentum().z(), mu1->currentState().mass());
                              mu22P4tmp.SetPxPyPzE(mu22Ktmp.momentum().x(), mu22Ktmp.momentum().y(), mu22Ktmp.momentum().z(), mu2->currentState().mass());
                              mu23P4tmp.SetPxPyPzE(mu23Ktmp.momentum().x(), mu23Ktmp.momentum().y(), mu23Ktmp.momentum().z(), mu3->currentState().mass());
                              mu24P4tmp.SetPxPyPzE(mu24Ktmp.momentum().x(), mu24Ktmp.momentum().y(), mu24Ktmp.momentum().z(), mu4->currentState().mass());
                              J21P4tmp.SetPxPyPzE(J21Ktmp.momentum().x(), J21Ktmp.momentum().y(), J21Ktmp.momentum().z(), jpsi1In2JpsiVertex->currentState().mass());
                              J22P4tmp.SetPxPyPzE(J22Ktmp.momentum().x(), J22Ktmp.momentum().y(), J22Ktmp.momentum().z(), jpsi2In2JpsiVertex->currentState().mass());
                              X2P4tmp.SetPxPyPzE(X2Ktmp.momentum().x(), X2Ktmp.momentum().y(), X2Ktmp.momentum().z(), MyChi1_part->currentState().mass());
                              Muon21Mass = mu1->currentState().mass();
                              Muon22Mass = mu2->currentState().mass();
                              Muon23Mass = mu3->currentState().mass();
                              Muon24Mass = mu4->currentState().mass();
                              J21normchi2 = myJpsi1DecayVtx->chiSquared() / myJpsi1DecayVtx->degreesOfFreedom();
                              J22normchi2 = myJpsi2DecayVtx->chiSquared() / myJpsi2DecayVtx->degreesOfFreedom();
                              J21Mass = jpsi1In2JpsiVertex->currentState().mass();
                              J22Mass = jpsi2In2JpsiVertex->currentState().mass();
                              X2Mass = MyChi1_part->currentState().mass();
                              X2normchi2 = (double)myJpsi1Jpsi2Vertex->chiSquared() / (double)myJpsi1Jpsi2Vertex->degreesOfFreedom();
                           }
                        }
                     }
                  }
               }
match2:
               // 14+23
               if ((iMuon1->charge() + iMuon4->charge()) == 0 && (iMuon2->charge() + iMuon3->charge()) == 0)
               {
                  KinematicParticleFactoryFromTransientTrack pmumuFactory;

                  // initial chi2 and ndf before kinematic fits.
                  float chi1 = 0.;
                  float chi2 = 0.;
                  float ndf1 = 0.;
                  float ndf2 = 0.;
                  float chi = 0.;
                  float ndf = 0.;

                  vector<RefCountedKinematicParticle> dimuon1Particles, dimuon2Particles, mu4Particles;
                  dimuon1Particles.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi1, ndf1, muon_sigma));
                  dimuon1Particles.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi1, ndf1, muon_sigma));
                  dimuon2Particles.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi2, ndf2, muon_sigma));
                  dimuon2Particles.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi2, ndf2, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi, ndf, muon_sigma));
                  mu4Particles.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi, ndf, muon_sigma));

                  if (dimuon1Particles.size() < 2 || dimuon2Particles.size() < 2 || mu4Particles.size() < 4)
                     goto match3;

                  KinematicParticleVertexFitter fitter1, fitter2, mu4fitter;
                  RefCountedKinematicTree psiVertexFitTree1, psiVertexFitTree2, XVertexFitTree;
                  psiVertexFitTree1 = fitter1.fit(dimuon1Particles);
                  psiVertexFitTree2 = fitter2.fit(dimuon2Particles);
                  XVertexFitTree = mu4fitter.fit(mu4Particles);

                  if (psiVertexFitTree1->isValid() && psiVertexFitTree2->isValid() && XVertexFitTree->isValid())
                  {
                     psiVertexFitTree1->movePointerToTheTop();
                     psiVertexFitTree2->movePointerToTheTop();
                     XVertexFitTree->movePointerToTheTop();

                     RefCountedKinematicParticle psi_vFit1 = psiVertexFitTree1->currentParticle();
                     RefCountedKinematicParticle psi_vFit2 = psiVertexFitTree2->currentParticle();
                     RefCountedKinematicParticle X_vFit = XVertexFitTree->currentParticle();

                     RefCountedKinematicVertex psi_vFit_vertex1 = psiVertexFitTree1->currentDecayVertex();
                     RefCountedKinematicVertex psi_vFit_vertex2 = psiVertexFitTree2->currentDecayVertex();
                     RefCountedKinematicVertex X_vFit_vertex = XVertexFitTree->currentDecayVertex();

                     if (!psi_vFit1->currentState().isValid() || !psi_vFit2->currentState().isValid() || !X_vFit->currentState().isValid())
                        goto match3;
                     if (!psi_vFit_vertex1->vertexIsValid() || !psi_vFit_vertex2->vertexIsValid() || !X_vFit_vertex->vertexIsValid())
                        goto match3;

                     // KinematicParameters Jpara1 = psi_vFit1->currentState().kinematicParameters();
                     // KinematicParameters Jpara2 = psi_vFit2->currentState().kinematicParameters();
                     // KinematicParameters Xpara = X_vFit->currentState().kinematicParameters();

                     XVertexFitTree->movePointerToTheFirstChild();
                     RefCountedKinematicParticle mu1 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu3 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu2 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu4 = XVertexFitTree->currentParticle();

                     double vProb1 = ChiSquaredProbability((double)(psi_vFit_vertex1->chiSquared()), (double)(psi_vFit_vertex1->degreesOfFreedom()));
                     double vProb2 = ChiSquaredProbability((double)(psi_vFit_vertex2->chiSquared()), (double)(psi_vFit_vertex2->degreesOfFreedom()));
                     // double XvProb = ChiSquaredProbability((double)(X_vFit_vertex->chiSquared()), (double)(X_vFit_vertex->degreesOfFreedom()));
                     if (vProb1 < JvPorbcut || vProb2 < JvPorbcut)
                        goto match3;
                     if (mu1->currentState().mass() <= MassMinCut || mu2->currentState().mass() <= MassMinCut || mu3->currentState().mass() <= MassMinCut || mu4->currentState().mass() <= MassMinCut || psi_vFit1->currentState().mass() <= MassMinCut || psi_vFit2->currentState().mass() <= MassMinCut || X_vFit->currentState().mass() <= MassMinCut)
                        goto match3;
                     if (psi_vFit1->currentState().kinematicParametersError().matrix()(6, 6) < 0 || psi_vFit2->currentState().kinematicParametersError().matrix()(6, 6) < 0 || X_vFit->currentState().kinematicParametersError().matrix()(6, 6) < 0)
                        goto match3;
                     float Jpsi1masserr = sqrt(psi_vFit1->currentState().kinematicParametersError().matrix()(6, 6));
                     float Jpsi2masserr = sqrt(psi_vFit2->currentState().kinematicParametersError().matrix()(6, 6));
                     if (psi_vFit1->currentState().mass() > (myYmass - 3.0 * Jpsi1masserr) && psi_vFit1->currentState().mass() < (myYmass + 3.0 * Jpsi1masserr) && psi_vFit2->currentState().mass() > (myJmass - 3.0 * Jpsi2masserr) && psi_vFit2->currentState().mass() < (myJmass + 3.0 * Jpsi2masserr))
                     {
                        Y3pos = 1;
                        J31NoMassKtmp = psi_vFit1->currentState().kinematicParameters();
                        J32NoMassKtmp = psi_vFit2->currentState().kinematicParameters();
                        J31NoMassP4tmp.SetPxPyPzE(J31NoMassKtmp.momentum().x(), J31NoMassKtmp.momentum().y(), J31NoMassKtmp.momentum().z(), psi_vFit1->currentState().mass());
                        J32NoMassP4tmp.SetPxPyPzE(J32NoMassKtmp.momentum().x(), J32NoMassKtmp.momentum().y(), J32NoMassKtmp.momentum().z(), psi_vFit2->currentState().mass());
                        J31NoMassMass = psi_vFit1->currentState().mass();
                        J32NoMassMass = psi_vFit2->currentState().mass();
                        J31NoMassMassE = Jpsi1masserr;
                        J32NoMassMassE = Jpsi2masserr;
                     }
                     if (psi_vFit2->currentState().mass() > (myYmass - 3.0 * Jpsi2masserr) && psi_vFit2->currentState().mass() < (myYmass + 3.0 * Jpsi2masserr) && psi_vFit1->currentState().mass() > (myJmass - 3.0 * Jpsi1masserr) && psi_vFit1->currentState().mass() < (myJmass + 3.0 * Jpsi1masserr) && vProb2 > vProb1)
                     {
                        Y3pos = 2;
                        J31NoMassKtmp = psi_vFit2->currentState().kinematicParameters();
                        J32NoMassKtmp = psi_vFit1->currentState().kinematicParameters();
                        J31NoMassP4tmp.SetPxPyPzE(J31NoMassKtmp.momentum().x(), J31NoMassKtmp.momentum().y(), J31NoMassKtmp.momentum().z(), psi_vFit2->currentState().mass());
                        J32NoMassP4tmp.SetPxPyPzE(J32NoMassKtmp.momentum().x(), J32NoMassKtmp.momentum().y(), J32NoMassKtmp.momentum().z(), psi_vFit1->currentState().mass());
                        J31NoMassMass = psi_vFit2->currentState().mass();
                        J32NoMassMass = psi_vFit1->currentState().mass();
                        J31NoMassMassE = Jpsi2masserr;
                        J32NoMassMassE = Jpsi1masserr;
                     }
                     else
                        goto match3;

                     // Mass Constraint Fit for J/psi
                     KinematicParticleVertexFitter kpvFitter;
                     KinematicParticleFitter csFitter;
                     ParticleMass jp1, jp2;
                     float jp_m_sigma1, jp_m_sigma2;
                     KinematicConstraint *jpsi_c1;
                     KinematicConstraint *jpsi_c2;
                     vector<RefCountedKinematicParticle> muonP12, muonP34;
                     RefCountedKinematicParticle Jpsi1_part, Jpsi2_part;
                     vector<RefCountedKinematicParticle> Chi_1;
                     muonP12.push_back(pmumuFactory.particle(muonTT1, muon_mass, chi1, ndf1, muon_sigma));
                     muonP12.push_back(pmumuFactory.particle(muonTT4, muon_mass, chi1, ndf1, muon_sigma));
                     muonP34.push_back(pmumuFactory.particle(muonTT2, muon_mass, chi2, ndf2, muon_sigma));
                     muonP34.push_back(pmumuFactory.particle(muonTT3, muon_mass, chi2, ndf2, muon_sigma));
                     RefCountedKinematicTree Jpsi1 = kpvFitter.fit(muonP12);
                     RefCountedKinematicTree Jpsi2 = kpvFitter.fit(muonP34);
                     RefCountedKinematicTree Jpsi1noMCJJ = kpvFitter.fit(muonP12);
                     RefCountedKinematicTree Jpsi2noMCJJ = kpvFitter.fit(muonP34);
                     RefCountedKinematicTree Chi1_bTree;
                     RefCountedKinematicParticle MyChi1_part;
                     if (Y3pos == 1)
                     {
                        jp1 = myYmass;
                        jp_m_sigma1 = myYmasserr;
                        jp2 = myJmass;
                        jp_m_sigma2 = myJmasserr;
                        jpsi_c1 = new MassKinematicConstraint(jp1, jp_m_sigma1);
                        jpsi_c2 = new MassKinematicConstraint(jp2, jp_m_sigma2);
                        try
                        {
                           Jpsi1 = csFitter.fit(jpsi_c1, Jpsi1noMCJJ);
                           Jpsi2 = csFitter.fit(jpsi_c2, Jpsi2noMCJJ);
                        }
                        catch (VertexException const &x)
                        {
                           std::cout << "mu12 vertex exception with mass constrainted to J!" << std::endl;
                        }
                        if (Jpsi1->isEmpty() != true && Jpsi2->isValid() == true)
                        {
                           Jpsi1->movePointerToTheTop();
                           Jpsi2->movePointerToTheTop();
                           Jpsi1_part = Jpsi1->currentParticle();
                           Jpsi2_part = Jpsi2->currentParticle();
                           Chi_1.push_back(Jpsi1_part);
                           Chi_1.push_back(Jpsi2_part);
                           bool isagoodfit = true;
                           try
                           {
                              Chi1_bTree = kpvFitter.fit(Chi_1);
                           }
                           catch (VertexException const &x)
                           {
                              isagoodfit = false;
                              cout << "mu12 and mu34 vertex exception with mu12 and mu34 constrained to JJ" << endl;
                           }
                           if (Chi1_bTree->isValid() && isagoodfit)
                           {
                              Chi1_bTree->movePointerToTheTop();
                              RefCountedKinematicVertex myJpsi1Jpsi2Vertex = Chi1_bTree->currentDecayVertex();
                              Chi1_bTree->movePointerToTheFirstChild();
                              RefCountedKinematicParticle jpsi1In2JpsiVertex = Chi1_bTree->currentParticle();
                              Chi1_bTree->movePointerToTheNextChild();
                              RefCountedKinematicParticle jpsi2In2JpsiVertex = Chi1_bTree->currentParticle();

                              RefCountedKinematicVertex myJpsi1DecayVtx = Jpsi1->currentDecayVertex();
                              RefCountedKinematicVertex myJpsi2DecayVtx = Jpsi2->currentDecayVertex();
                              mu31Ktmp = mu1->currentState().kinematicParameters();
                              mu32Ktmp = mu2->currentState().kinematicParameters();
                              mu33Ktmp = mu3->currentState().kinematicParameters();
                              mu34Ktmp = mu4->currentState().kinematicParameters();
                              J31Ktmp = jpsi1In2JpsiVertex->currentState().kinematicParameters();
                              J32Ktmp = jpsi2In2JpsiVertex->currentState().kinematicParameters();
                              Chi1_bTree->movePointerToTheTop();
                              MyChi1_part = Chi1_bTree->currentParticle();
                              X3Ktmp = MyChi1_part->currentState().kinematicParameters();
                              mu31P4tmp.SetPxPyPzE(mu31Ktmp.momentum().x(), mu31Ktmp.momentum().y(), mu31Ktmp.momentum().z(), mu1->currentState().mass());
                              mu32P4tmp.SetPxPyPzE(mu32Ktmp.momentum().x(), mu32Ktmp.momentum().y(), mu32Ktmp.momentum().z(), mu2->currentState().mass());
                              mu33P4tmp.SetPxPyPzE(mu33Ktmp.momentum().x(), mu33Ktmp.momentum().y(), mu33Ktmp.momentum().z(), mu3->currentState().mass());
                              mu34P4tmp.SetPxPyPzE(mu34Ktmp.momentum().x(), mu34Ktmp.momentum().y(), mu34Ktmp.momentum().z(), mu4->currentState().mass());
                              J31P4tmp.SetPxPyPzE(J31Ktmp.momentum().x(), J31Ktmp.momentum().y(), J31Ktmp.momentum().z(), jpsi1In2JpsiVertex->currentState().mass());
                              J32P4tmp.SetPxPyPzE(J32Ktmp.momentum().x(), J32Ktmp.momentum().y(), J32Ktmp.momentum().z(), jpsi2In2JpsiVertex->currentState().mass());
                              X3P4tmp.SetPxPyPzE(X3Ktmp.momentum().x(), X3Ktmp.momentum().y(), X3Ktmp.momentum().z(), MyChi1_part->currentState().mass());
                              Muon31Mass = mu1->currentState().mass();
                              Muon32Mass = mu2->currentState().mass();
                              Muon33Mass = mu3->currentState().mass();
                              Muon34Mass = mu4->currentState().mass();
                              J31normchi2 = myJpsi1DecayVtx->chiSquared() / myJpsi1DecayVtx->degreesOfFreedom();
                              J32normchi2 = myJpsi2DecayVtx->chiSquared() / myJpsi2DecayVtx->degreesOfFreedom();
                              J31Mass = jpsi1In2JpsiVertex->currentState().mass();
                              J32Mass = jpsi2In2JpsiVertex->currentState().mass();
                              X3Mass = MyChi1_part->currentState().mass();
                              X3normchi2 = (double)myJpsi1Jpsi2Vertex->chiSquared() / (double)myJpsi1Jpsi2Vertex->degreesOfFreedom();
                           }
                        }
                     }
                     else if (Y3pos == 2)
                     {
                        jp2 = myYmass;
                        jp_m_sigma2 = myYmasserr;
                        jp1 = myJmass;
                        jp_m_sigma1 = myJmasserr;
                        jpsi_c2 = new MassKinematicConstraint(jp2, jp_m_sigma2);
                        jpsi_c1 = new MassKinematicConstraint(jp1, jp_m_sigma1);
                        try
                        {
                           Jpsi2 = csFitter.fit(jpsi_c2, Jpsi2noMCJJ);
                           Jpsi1 = csFitter.fit(jpsi_c1, Jpsi1noMCJJ);
                        }
                        catch (VertexException const &x)
                        {
                           std::cout << "mu34 vertex exception with mass constrainted to J!" << std::endl;
                        }
                        if (Jpsi2->isEmpty() != true && Jpsi1->isValid() == true)
                        {
                           Jpsi2->movePointerToTheTop();
                           Jpsi1->movePointerToTheTop();
                           Jpsi1_part = Jpsi2->currentParticle();
                           Jpsi2_part = Jpsi1->currentParticle();
                           Chi_1.push_back(Jpsi1_part);
                           Chi_1.push_back(Jpsi2_part);
                           bool isagoodfit = true;
                           try
                           {
                              Chi1_bTree = kpvFitter.fit(Chi_1);
                           }
                           catch (VertexException const &x)
                           {
                              isagoodfit = false;
                              cout << "mu12 and mu34 vertex exception with mu12 and mu34 constrained to JJ" << endl;
                           }
                           if (Chi1_bTree->isValid() && isagoodfit)
                           {
                              Chi1_bTree->movePointerToTheTop();
                              RefCountedKinematicVertex myJpsi1Jpsi2Vertex = Chi1_bTree->currentDecayVertex();
                              Chi1_bTree->movePointerToTheFirstChild();
                              RefCountedKinematicParticle jpsi1In2JpsiVertex = Chi1_bTree->currentParticle();
                              Chi1_bTree->movePointerToTheNextChild();
                              RefCountedKinematicParticle jpsi2In2JpsiVertex = Chi1_bTree->currentParticle();

                              RefCountedKinematicVertex myJpsi1DecayVtx = Jpsi2->currentDecayVertex();
                              RefCountedKinematicVertex myJpsi2DecayVtx = Jpsi1->currentDecayVertex();
                              mu31Ktmp = mu1->currentState().kinematicParameters();
                              mu32Ktmp = mu2->currentState().kinematicParameters();
                              mu33Ktmp = mu3->currentState().kinematicParameters();
                              mu34Ktmp = mu4->currentState().kinematicParameters();
                              J31Ktmp = jpsi1In2JpsiVertex->currentState().kinematicParameters();
                              J32Ktmp = jpsi2In2JpsiVertex->currentState().kinematicParameters();
                              Chi1_bTree->movePointerToTheTop();
                              MyChi1_part = Chi1_bTree->currentParticle();
                              X3Ktmp = MyChi1_part->currentState().kinematicParameters();
                              mu31P4tmp.SetPxPyPzE(mu31Ktmp.momentum().x(), mu31Ktmp.momentum().y(), mu31Ktmp.momentum().z(), mu1->currentState().mass());
                              mu32P4tmp.SetPxPyPzE(mu32Ktmp.momentum().x(), mu32Ktmp.momentum().y(), mu32Ktmp.momentum().z(), mu2->currentState().mass());
                              mu33P4tmp.SetPxPyPzE(mu33Ktmp.momentum().x(), mu33Ktmp.momentum().y(), mu33Ktmp.momentum().z(), mu3->currentState().mass());
                              mu34P4tmp.SetPxPyPzE(mu34Ktmp.momentum().x(), mu34Ktmp.momentum().y(), mu34Ktmp.momentum().z(), mu4->currentState().mass());
                              J31P4tmp.SetPxPyPzE(J31Ktmp.momentum().x(), J31Ktmp.momentum().y(), J31Ktmp.momentum().z(), jpsi1In2JpsiVertex->currentState().mass());
                              J32P4tmp.SetPxPyPzE(J32Ktmp.momentum().x(), J32Ktmp.momentum().y(), J32Ktmp.momentum().z(), jpsi2In2JpsiVertex->currentState().mass());
                              X3P4tmp.SetPxPyPzE(X3Ktmp.momentum().x(), X3Ktmp.momentum().y(), X3Ktmp.momentum().z(), MyChi1_part->currentState().mass());
                              Muon31Mass = mu1->currentState().mass();
                              Muon32Mass = mu2->currentState().mass();
                              Muon33Mass = mu3->currentState().mass();
                              Muon34Mass = mu4->currentState().mass();
                              J31normchi2 = myJpsi1DecayVtx->chiSquared() / myJpsi1DecayVtx->degreesOfFreedom();
                              J32normchi2 = myJpsi2DecayVtx->chiSquared() / myJpsi2DecayVtx->degreesOfFreedom();
                              J31Mass = jpsi1In2JpsiVertex->currentState().mass();
                              J32Mass = jpsi2In2JpsiVertex->currentState().mass();
                              X3Mass = MyChi1_part->currentState().mass();
                              X3normchi2 = (double)myJpsi1Jpsi2Vertex->chiSquared() / (double)myJpsi1Jpsi2Vertex->degreesOfFreedom();
                           }
                        }
                     }
                  }
               }
match3:
               // Choose Best Combination to Save
               if (X1normchi2 <= X2normchi2 && X1normchi2 <= X3normchi2 && X1normchi2 >= 0 && X1normchi2 < 900 && npairs < 36)
               {
                  Muon1Charge_[npairs] = iMuon1->charge();
                  Muon2Charge_[npairs] = iMuon2->charge();
                  Muon3Charge_[npairs] = iMuon3->charge();
                  Muon4Charge_[npairs] = iMuon4->charge();
                  Muon1Pt_[npairs] = mu11P4tmp.Pt();
                  Muon2Pt_[npairs] = mu12P4tmp.Pt();
                  Muon3Pt_[npairs] = mu13P4tmp.Pt();
                  Muon4Pt_[npairs] = mu14P4tmp.Pt();
                  Muon1Eta_[npairs] = mu11P4tmp.Eta();
                  Muon2Eta_[npairs] = mu12P4tmp.Eta();
                  Muon3Eta_[npairs] = mu13P4tmp.Eta();
                  Muon4Eta_[npairs] = mu14P4tmp.Eta();
                  Muon1Phi_[npairs] = mu11P4tmp.Phi();
                  Muon2Phi_[npairs] = mu12P4tmp.Phi();
                  Muon3Phi_[npairs] = mu13P4tmp.Phi();
                  Muon4Phi_[npairs] = mu14P4tmp.Phi();
                  Muon1Mass_[npairs] = Muon11Mass;
                  Muon2Mass_[npairs] = Muon12Mass;
                  Muon3Mass_[npairs] = Muon13Mass;
                  Muon4Mass_[npairs] = Muon14Mass;
                  J1normchi2_[npairs] = J11normchi2;
                  J2normchi2_[npairs] = J12normchi2;
                  J1Pt_[npairs] = J11P4tmp.Pt();
                  J2Pt_[npairs] = J12P4tmp.Pt();
                  J1Eta_[npairs] = J11P4tmp.Eta();
                  J2Eta_[npairs] = J12P4tmp.Eta();
                  J1Phi_[npairs] = J11P4tmp.Phi();
                  J2Phi_[npairs] = J12P4tmp.Phi();
                  J1Mass_[npairs] = J11Mass;
                  J2Mass_[npairs] = J12Mass;
                  XpT_[npairs] = X1P4tmp.Pt();
                  Xeta_[npairs] = X1P4tmp.Eta();
                  Xphi_[npairs] = X1P4tmp.Phi();
                  Xmass_[npairs] = X1Mass;
                  Xnormchi2_[npairs] = X1normchi2;

                  J1NoMassnormchi2_[npairs] = J11NoMassnormchi2;
                  J2NoMassnormchi2_[npairs] = J12NoMassnormchi2;
                  J1NoMassPt_[npairs] = J11NoMassP4tmp.Pt();
                  J2NoMassPt_[npairs] = J12NoMassP4tmp.Pt();
                  J1NoMassEta_[npairs] = J11NoMassP4tmp.Eta();
                  J2NoMassEta_[npairs] = J12NoMassP4tmp.Eta();
                  J1NoMassPhi_[npairs] = J11NoMassP4tmp.Phi();
                  J2NoMassPhi_[npairs] = J12NoMassP4tmp.Phi();
                  J1NoMassMass_[npairs] = J11NoMassMass;
                  J2NoMassMass_[npairs] = J12NoMassMass;
                  J1NoMassMassE_[npairs] = J11NoMassMassE;
                  J2NoMassMassE_[npairs] = J12NoMassMassE;
                  type_[npairs] = 1;

                  npairs++;
               }
               else if (X2normchi2 <= X1normchi2 && X2normchi2 <= X3normchi2 && X2normchi2 >= 0 && X2normchi2 < 900 && npairs < 36)
               {
                  Muon1Charge_[npairs] = iMuon1->charge();
                  Muon2Charge_[npairs] = iMuon2->charge();
                  Muon3Charge_[npairs] = iMuon3->charge();
                  Muon4Charge_[npairs] = iMuon4->charge();
                  Muon1Pt_[npairs] = mu21P4tmp.Pt();
                  Muon2Pt_[npairs] = mu22P4tmp.Pt();
                  Muon3Pt_[npairs] = mu23P4tmp.Pt();
                  Muon4Pt_[npairs] = mu24P4tmp.Pt();
                  Muon1Eta_[npairs] = mu21P4tmp.Eta();
                  Muon2Eta_[npairs] = mu22P4tmp.Eta();
                  Muon3Eta_[npairs] = mu23P4tmp.Eta();
                  Muon4Eta_[npairs] = mu24P4tmp.Eta();
                  Muon1Phi_[npairs] = mu21P4tmp.Phi();
                  Muon2Phi_[npairs] = mu22P4tmp.Phi();
                  Muon3Phi_[npairs] = mu23P4tmp.Phi();
                  Muon4Phi_[npairs] = mu24P4tmp.Phi();
                  Muon1Mass_[npairs] = Muon21Mass;
                  Muon2Mass_[npairs] = Muon22Mass;
                  Muon3Mass_[npairs] = Muon23Mass;
                  Muon4Mass_[npairs] = Muon24Mass;
                  J1normchi2_[npairs] = J21normchi2;
                  J2normchi2_[npairs] = J22normchi2;
                  J1Pt_[npairs] = J21P4tmp.Pt();
                  J2Pt_[npairs] = J22P4tmp.Pt();
                  J1Eta_[npairs] = J21P4tmp.Eta();
                  J2Eta_[npairs] = J22P4tmp.Eta();
                  J1Phi_[npairs] = J21P4tmp.Phi();
                  J2Phi_[npairs] = J22P4tmp.Phi();
                  J1Mass_[npairs] = J21Mass;
                  J2Mass_[npairs] = J22Mass;
                  XpT_[npairs] = X2P4tmp.Pt();
                  Xeta_[npairs] = X2P4tmp.Eta();
                  Xphi_[npairs] = X2P4tmp.Phi();
                  Xmass_[npairs] = X2Mass;
                  Xnormchi2_[npairs] = X2normchi2;

                  J1NoMassnormchi2_[npairs] = J21NoMassnormchi2;
                  J2NoMassnormchi2_[npairs] = J22NoMassnormchi2;
                  J1NoMassPt_[npairs] = J21NoMassP4tmp.Pt();
                  J2NoMassPt_[npairs] = J22NoMassP4tmp.Pt();
                  J1NoMassEta_[npairs] = J21NoMassP4tmp.Eta();
                  J2NoMassEta_[npairs] = J22NoMassP4tmp.Eta();
                  J1NoMassPhi_[npairs] = J21NoMassP4tmp.Phi();
                  J2NoMassPhi_[npairs] = J22NoMassP4tmp.Phi();
                  J1NoMassMass_[npairs] = J21NoMassMass;
                  J2NoMassMass_[npairs] = J22NoMassMass;
                  J1NoMassMassE_[npairs] = J21NoMassMassE;
                  J2NoMassMassE_[npairs] = J22NoMassMassE;
                  type_[npairs] = 2;

                  npairs++;
               }
               else if (X3normchi2 <= X1normchi2 && X3normchi2 <= X2normchi2 && X3normchi2 >= 0 && X3normchi2 < 900 && npairs < 36)
               {
                  Muon1Charge_[npairs] = iMuon1->charge();
                  Muon2Charge_[npairs] = iMuon2->charge();
                  Muon3Charge_[npairs] = iMuon3->charge();
                  Muon4Charge_[npairs] = iMuon4->charge();
                  Muon1Pt_[npairs] = mu31P4tmp.Pt();
                  Muon2Pt_[npairs] = mu32P4tmp.Pt();
                  Muon3Pt_[npairs] = mu33P4tmp.Pt();
                  Muon4Pt_[npairs] = mu34P4tmp.Pt();
                  Muon1Eta_[npairs] = mu31P4tmp.Eta();
                  Muon2Eta_[npairs] = mu32P4tmp.Eta();
                  Muon3Eta_[npairs] = mu33P4tmp.Eta();
                  Muon4Eta_[npairs] = mu34P4tmp.Eta();
                  Muon1Phi_[npairs] = mu31P4tmp.Phi();
                  Muon2Phi_[npairs] = mu32P4tmp.Phi();
                  Muon3Phi_[npairs] = mu33P4tmp.Phi();
                  Muon4Phi_[npairs] = mu34P4tmp.Phi();
                  Muon1Mass_[npairs] = Muon31Mass;
                  Muon2Mass_[npairs] = Muon32Mass;
                  Muon3Mass_[npairs] = Muon33Mass;
                  Muon4Mass_[npairs] = Muon34Mass;
                  J1normchi2_[npairs] = J31normchi2;
                  J2normchi2_[npairs] = J32normchi2;
                  J1Pt_[npairs] = J31P4tmp.Pt();
                  J2Pt_[npairs] = J32P4tmp.Pt();
                  J1Eta_[npairs] = J31P4tmp.Eta();
                  J2Eta_[npairs] = J32P4tmp.Eta();
                  J1Phi_[npairs] = J31P4tmp.Phi();
                  J2Phi_[npairs] = J32P4tmp.Phi();
                  J1Mass_[npairs] = J31Mass;
                  J2Mass_[npairs] = J32Mass;
                  XpT_[npairs] = X3P4tmp.Pt();
                  Xeta_[npairs] = X3P4tmp.Eta();
                  Xphi_[npairs] = X3P4tmp.Phi();
                  Xmass_[npairs] = X3Mass;
                  Xnormchi2_[npairs] = X3normchi2;

                  J1NoMassnormchi2_[npairs] = J31NoMassnormchi2;
                  J2NoMassnormchi2_[npairs] = J32NoMassnormchi2;
                  J1NoMassPt_[npairs] = J31NoMassP4tmp.Pt();
                  J2NoMassPt_[npairs] = J32NoMassP4tmp.Pt();
                  J1NoMassEta_[npairs] = J31NoMassP4tmp.Eta();
                  J2NoMassEta_[npairs] = J32NoMassP4tmp.Eta();
                  J1NoMassPhi_[npairs] = J31NoMassP4tmp.Phi();
                  J2NoMassPhi_[npairs] = J32NoMassP4tmp.Phi();
                  J1NoMassMass_[npairs] = J31NoMassMass;
                  J2NoMassMass_[npairs] = J32NoMassMass;
                  J1NoMassMassE_[npairs] = J31NoMassMassE;
                  J2NoMassMassE_[npairs] = J32NoMassMassE;
                  type_[npairs] = 3;

                  npairs++;
               }
            }
         }
      }
   }

   if (npairs > 0)
   {
      X4muTree->Fill();
      // std::cout << "npairs: " << npairs << " XMass: " << Xmass_[0] << " J1Mass: " << J1Mass_[0] << " J2Mass: " << J2Mass_[0] << " Muon1Mass: " << Muon1Mass_[0] << " Muon2Mass: " << Muon2Mass_[0] << " Muon3Mass: " << Muon3Mass_[0] << " Muon4Mass: " << Muon4Mass_[0] << std::endl;
   }
   npairs = 0;
   clearVars();
}

void X4muSecondaryVertexProducer::clearVars()
{
   for (size_t i = 0; i < 36; i++)
   {
      Muon1Charge_[i] = 0;
      Muon2Charge_[i] = 0;
      Muon3Charge_[i] = 0;
      Muon4Charge_[i] = 0;
      Muon1Pt_[i] = 0;
      Muon2Pt_[i] = 0;
      Muon3Pt_[i] = 0;
      Muon4Pt_[i] = 0;
      Muon1Eta_[i] = 0;
      Muon2Eta_[i] = 0;
      Muon3Eta_[i] = 0;
      Muon4Eta_[i] = 0;
      Muon1Phi_[i] = 0;
      Muon2Phi_[i] = 0;
      Muon3Phi_[i] = 0;
      Muon4Phi_[i] = 0;
      Muon1Mass_[i] = 0;
      Muon2Mass_[i] = 0;
      Muon3Mass_[i] = 0;
      Muon4Mass_[i] = 0;
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
      type_[i] = 0;
   }
}

void X4muSecondaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   edm::ParameterSetDescription desc;
   desc.add<edm::InputTag>("recoMuon", edm::InputTag("recoMuon"));
   desc.add<edm::InputTag>("recoTrack", edm::InputTag("recoTrack"));
   desc.add<double>("MesonMassBig");
   desc.add<double>("MesonMassBigErr");
   desc.add<double>("MesonMassSmall");
   desc.add<double>("MesonMassSmallErr");
   descriptions.add("X4muSecondaryVertexProducer", desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(X4muSecondaryVertexProducer);