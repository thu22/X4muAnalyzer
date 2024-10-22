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
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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
   float J1Pt_[36];
   float J2Pt_[36];
   float J1Eta_[36];
   float J2Eta_[36];
   float J1Phi_[36];
   float J2Phi_[36];
   float J1Mass_[36];
   float J2Mass_[36];
   float XpT_[36];
   float Xeta_[36];
   float Xphi_[36];
   float Xmass_[36];
   float Xnormchi2_[36];
   /* std::vector<int> Muon3Charge_;
   std::vector<int> Muon4Charge_;
   std::vector<float> Muon1Pt_;
   std::vector<float> Muon2Pt_;
   std::vector<float> Muon3Pt_;
   std::vector<float> Muon4Pt_;
   std::vector<float> Muon1Eta_;
   std::vector<float> Muon2Eta_;
   std::vector<float> Muon3Eta_;
   std::vector<float> Muon4Eta_;
   std::vector<float> Muon1Phi_;
   std::vector<float> Muon2Phi_;
   std::vector<float> Muon3Phi_;
   std::vector<float> Muon4Phi_;
   std::vector<float> Muon1Mass_;
   std::vector<float> Muon2Mass_;
   std::vector<float> Muon3Mass_;
   std::vector<float> Muon4Mass_;
   std::vector<float> J1normchi2_;
   std::vector<float> J2normchi2_;
   std::vector<float> J1Pt_;
   std::vector<float> J2Pt_;
   std::vector<float> J1Eta_;
   std::vector<float> J2Eta_;
   std::vector<float> J1Phi_;
   std::vector<float> J2Phi_;
   std::vector<float> J1Mass_;
   std::vector<float> J2Mass_;
   std::vector<float> XpT_;
   std::vector<float> Xeta_;
   std::vector<float> Xphi_;
   std::vector<float> Xmass_;
   std::vector<float> Xnormchi2_; */
};

//
// constructors and destructor
//
X4muSecondaryVertexProducer::X4muSecondaryVertexProducer(edm::ParameterSet const &iConfig)
    : input_recomuon_token_(consumes(iConfig.getParameter<edm::InputTag>("recoMuon"))),
      input_recoTrack_token_(consumes(iConfig.getParameter<edm::InputTag>("recoTrack"))),
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
   X4muTree->Branch("J1Pt", J1Pt_);
   X4muTree->Branch("J2Pt", J2Pt_);
   X4muTree->Branch("J1Eta", J1Eta_);
   X4muTree->Branch("J2Eta", J2Eta_);
   X4muTree->Branch("J1Phi", J1Phi_);
   X4muTree->Branch("J2Phi", J2Phi_);
   X4muTree->Branch("J1Mass", J1Mass_);
   X4muTree->Branch("J2Mass", J2Mass_);
   X4muTree->Branch("XpT", XpT_);
   X4muTree->Branch("Xeta", Xeta_);
   X4muTree->Branch("Xphi", Xphi_);
   X4muTree->Branch("Xmass", Xmass_);
   X4muTree->Branch("Xnormchi2", Xnormchi2_);

   // register products
   /* produces<edm::ValueMap<int>>("Muon1Charge");
   produces<edm::ValueMap<int>>("Muon2Charge");
   produces<edm::ValueMap<int>>("Muon3Charge");
   produces<edm::ValueMap<int>>("Muon4Charge");
   produces<edm::ValueMap<float>>("Muon1Pt");
   produces<edm::ValueMap<float>>("Muon2Pt");
   produces<edm::ValueMap<float>>("Muon3Pt");
   produces<edm::ValueMap<float>>("Muon4Pt");
   produces<edm::ValueMap<float>>("Muon1Eta");
   produces<edm::ValueMap<float>>("Muon2Eta");
   produces<edm::ValueMap<float>>("Muon3Eta");
   produces<edm::ValueMap<float>>("Muon4Eta");
   produces<edm::ValueMap<float>>("Muon1Phi");
   produces<edm::ValueMap<float>>("Muon2Phi");
   produces<edm::ValueMap<float>>("Muon3Phi");
   produces<edm::ValueMap<float>>("Muon4Phi");
   produces<edm::ValueMap<float>>("Muon1Mass");
   produces<edm::ValueMap<float>>("Muon2Mass");
   produces<edm::ValueMap<float>>("Muon3Mass");
   produces<edm::ValueMap<float>>("Muon4Mass");
   produces<edm::ValueMap<float>>("J1normchi2");
   produces<edm::ValueMap<float>>("J2normchi2");
   produces<edm::ValueMap<float>>("J1Pt");
   produces<edm::ValueMap<float>>("J2Pt");
   produces<edm::ValueMap<float>>("J1Eta");
   produces<edm::ValueMap<float>>("J2Eta");
   produces<edm::ValueMap<float>>("J1Phi");
   produces<edm::ValueMap<float>>("J2Phi");
   produces<edm::ValueMap<float>>("J1Mass");
   produces<edm::ValueMap<float>>("J2Mass");
   produces<edm::ValueMap<float>>("XpT");
   produces<edm::ValueMap<float>>("Xeta");
   produces<edm::ValueMap<float>>("Xphi");
   produces<edm::ValueMap<float>>("Xmass");
   produces<edm::ValueMap<int>>("Xnormchi2"); */
}

X4muSecondaryVertexProducer::~X4muSecondaryVertexProducer() = default;

// ------------ method called to produce the data  ------------
void X4muSecondaryVertexProducer::produce(edm::Event &iEvent, edm::EventSetup const &iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   const MagneticField &bFieldHandle = iSetup.getData(magneticFieldToken_); // Todo BField//

   Handle<std::vector<reco::Muon>> muonHandle;
   Handle<std::vector<reco::Track>> trackHandle;
   iEvent.getByToken(input_recomuon_token_, muonHandle);
   iEvent.getByToken(input_recoTrack_token_, trackHandle);

   float myMumass = 0.1056583755;
   double myMumasserr = myMumass * 1e-6;
   float JvPorbcut = 0.0001;
   long npairs = 0;

   // Store all possible 4 muons with 12 pair + 34 pair
   /* vector<TLorentzVector> *mu1P4 = new vector<TLorentzVector>;
   vector<TLorentzVector> *mu2P4 = new vector<TLorentzVector>;
   vector<TLorentzVector> *mu3P4 = new vector<TLorentzVector>;
   vector<TLorentzVector> *mu4P4 = new vector<TLorentzVector>;
   vector<TLorentzVector> *J1P4 = new vector<TLorentzVector>;
   vector<TLorentzVector> *J2P4 = new vector<TLorentzVector>;
   vector<TLorentzVector> *XP4 = new vector<TLorentzVector>;
   vector<float> *Dimuon1vProb = new vector<float>;
   vector<float> *Dimuon2vProb = new vector<float>;
   vector<float> *JJvProb = new vector<float>; */

   if (!trackHandle.isValid())
   {
      std::cout << "reco::Track collection is not valid" << std::endl;
      return;
   }

   if (trackHandle->size() < 4) // 4
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
               // bool mu4pair[3] = {false, false, false}; //[0] for 12+34, [1] for 13+24, [2] for 14+23

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
                     continue;

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
                        continue;
                     if (!psi_vFit_vertex1->vertexIsValid() || !psi_vFit_vertex2->vertexIsValid() || !X_vFit_vertex->vertexIsValid())
                        continue;

                     //KinematicParameters Jpara1 = psi_vFit1->currentState().kinematicParameters();
                     //KinematicParameters Jpara2 = psi_vFit2->currentState().kinematicParameters();
                     //KinematicParameters Xpara = X_vFit->currentState().kinematicParameters();

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
                     //double XvProb = ChiSquaredProbability((double)(X_vFit_vertex->chiSquared()), (double)(X_vFit_vertex->degreesOfFreedom()));
                     if (vProb1 < JvPorbcut || vProb2 < JvPorbcut)
                        continue;
                     if (mu1->currentState().mass() <= 0 || mu2->currentState().mass() <= 0 || mu3->currentState().mass() <= 0 || mu4->currentState().mass() <= 0 || psi_vFit1->currentState().mass() <= 0 || psi_vFit2->currentState().mass() <= 0 || X_vFit->currentState().mass() <= 0)
                        continue;

                     if (npairs < 36)
                     {
                        TLorentzVector mu1P4tmp, mu2P4tmp, mu3P4tmp, mu4P4tmp, J1P4tmp, J2P4tmp, XP4tmp;
                        KinematicParameters mu1Ktmp, mu2Ktmp, mu3Ktmp, mu4Ktmp, J1Ktmp, J2Ktmp, XKtmp;
                        mu1Ktmp = mu1->currentState().kinematicParameters();
                        mu2Ktmp = mu2->currentState().kinematicParameters();
                        mu3Ktmp = mu3->currentState().kinematicParameters();
                        mu4Ktmp = mu4->currentState().kinematicParameters();
                        J1Ktmp = psi_vFit1->currentState().kinematicParameters();
                        J2Ktmp = psi_vFit2->currentState().kinematicParameters();
                        XKtmp = X_vFit->currentState().kinematicParameters();
                        Muon1Charge_[npairs] = iMuon1->charge();
                        Muon2Charge_[npairs] = iMuon2->charge();
                        Muon3Charge_[npairs] = iMuon3->charge();
                        Muon4Charge_[npairs] = iMuon4->charge();
                        mu1P4tmp.SetPxPyPzE(mu1Ktmp.momentum().x(), mu1Ktmp.momentum().y(), mu1Ktmp.momentum().z(), mu1->currentState().mass());
                        mu2P4tmp.SetPxPyPzE(mu2Ktmp.momentum().x(), mu2Ktmp.momentum().y(), mu2Ktmp.momentum().z(), mu2->currentState().mass());
                        mu3P4tmp.SetPxPyPzE(mu3Ktmp.momentum().x(), mu3Ktmp.momentum().y(), mu3Ktmp.momentum().z(), mu3->currentState().mass());
                        mu4P4tmp.SetPxPyPzE(mu4Ktmp.momentum().x(), mu4Ktmp.momentum().y(), mu4Ktmp.momentum().z(), mu4->currentState().mass());
                        J1P4tmp.SetPxPyPzE(J1Ktmp.momentum().x(), J1Ktmp.momentum().y(), J1Ktmp.momentum().z(), psi_vFit1->currentState().mass());
                        J2P4tmp.SetPxPyPzE(J2Ktmp.momentum().x(), J2Ktmp.momentum().y(), J2Ktmp.momentum().z(), psi_vFit2->currentState().mass());
                        XP4tmp.SetPxPyPzE(XKtmp.momentum().x(), XKtmp.momentum().y(), XKtmp.momentum().z(), X_vFit->currentState().mass());
                        Muon1Pt_[npairs] = mu1P4tmp.Pt();
                        Muon2Pt_[npairs] = mu2P4tmp.Pt();
                        Muon3Pt_[npairs] = mu3P4tmp.Pt();
                        Muon4Pt_[npairs] = mu4P4tmp.Pt();
                        Muon1Eta_[npairs] = mu1P4tmp.Eta();
                        Muon2Eta_[npairs] = mu2P4tmp.Eta();
                        Muon3Eta_[npairs] = mu3P4tmp.Eta();
                        Muon4Eta_[npairs] = mu4P4tmp.Eta();
                        Muon1Phi_[npairs] = mu1P4tmp.Phi();
                        Muon2Phi_[npairs] = mu2P4tmp.Phi();
                        Muon3Phi_[npairs] = mu3P4tmp.Phi();
                        Muon4Phi_[npairs] = mu4P4tmp.Phi();
                        Muon1Mass_[npairs] = mu1->currentState().mass();
                        Muon2Mass_[npairs] = mu2->currentState().mass();
                        Muon3Mass_[npairs] = mu3->currentState().mass();
                        Muon4Mass_[npairs] = mu4->currentState().mass();
                        J1normchi2_[npairs] = psi_vFit1->chiSquared() / psi_vFit1->degreesOfFreedom();
                        J2normchi2_[npairs] = psi_vFit2->chiSquared() / psi_vFit2->degreesOfFreedom();
                        J1Pt_[npairs] = J1P4tmp.Pt();
                        J2Pt_[npairs] = J2P4tmp.Pt();
                        J1Eta_[npairs] = J1P4tmp.Eta();
                        J2Eta_[npairs] = J2P4tmp.Eta();
                        J1Phi_[npairs] = J1P4tmp.Phi();
                        J2Phi_[npairs] = J2P4tmp.Phi();
                        J1Mass_[npairs] = psi_vFit1->currentState().mass();
                        J2Mass_[npairs] = psi_vFit2->currentState().mass();
                        XpT_[npairs] = XP4tmp.Pt();
                        Xeta_[npairs] = XP4tmp.Eta();
                        Xphi_[npairs] = XP4tmp.Phi();
                        Xmass_[npairs] = X_vFit->currentState().mass();
                        Xnormchi2_[npairs] = (double) X_vFit->chiSquared() / (double) X_vFit->degreesOfFreedom();
                     }
                     npairs++;
                     // mu4pair[0] = true;
                  }
               }

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
                     continue;

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
                        continue;
                     if (!psi_vFit_vertex1->vertexIsValid() || !psi_vFit_vertex2->vertexIsValid() || !X_vFit_vertex->vertexIsValid())
                        continue;

                     //KinematicParameters Jpara1 = psi_vFit1->currentState().kinematicParameters();
                     //KinematicParameters Jpara2 = psi_vFit2->currentState().kinematicParameters();
                     //KinematicParameters Xpara = X_vFit->currentState().kinematicParameters();

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
                     //double XvProb = ChiSquaredProbability((double)(X_vFit_vertex->chiSquared()), (double)(X_vFit_vertex->degreesOfFreedom()));
                     if (vProb1 < JvPorbcut || vProb2 < JvPorbcut)
                        continue;
                     if (mu1->currentState().mass() <= 0 || mu2->currentState().mass() <= 0 || mu3->currentState().mass() <= 0 || mu4->currentState().mass() <= 0 || psi_vFit1->currentState().mass() <= 0 || psi_vFit2->currentState().mass() <= 0 || X_vFit->currentState().mass() <= 0)
                        continue;

                     if (npairs < 36)
                     {
                        TLorentzVector mu1P4tmp, mu2P4tmp, mu3P4tmp, mu4P4tmp, J1P4tmp, J2P4tmp, XP4tmp;
                        KinematicParameters mu1Ktmp, mu2Ktmp, mu3Ktmp, mu4Ktmp, J1Ktmp, J2Ktmp, XKtmp;
                        mu1Ktmp = mu1->currentState().kinematicParameters();
                        mu2Ktmp = mu2->currentState().kinematicParameters();
                        mu3Ktmp = mu3->currentState().kinematicParameters();
                        mu4Ktmp = mu4->currentState().kinematicParameters();
                        J1Ktmp = psi_vFit1->currentState().kinematicParameters();
                        J2Ktmp = psi_vFit2->currentState().kinematicParameters();
                        XKtmp = X_vFit->currentState().kinematicParameters();
                        Muon1Charge_[npairs] = iMuon1->charge();
                        Muon2Charge_[npairs] = iMuon3->charge();
                        Muon3Charge_[npairs] = iMuon2->charge();
                        Muon4Charge_[npairs] = iMuon4->charge();
                        mu1P4tmp.SetPxPyPzE(mu1Ktmp.momentum().x(), mu1Ktmp.momentum().y(), mu1Ktmp.momentum().z(), mu1->currentState().mass());
                        mu2P4tmp.SetPxPyPzE(mu3Ktmp.momentum().x(), mu3Ktmp.momentum().y(), mu3Ktmp.momentum().z(), mu3->currentState().mass());
                        mu3P4tmp.SetPxPyPzE(mu2Ktmp.momentum().x(), mu2Ktmp.momentum().y(), mu2Ktmp.momentum().z(), mu2->currentState().mass());
                        mu4P4tmp.SetPxPyPzE(mu4Ktmp.momentum().x(), mu4Ktmp.momentum().y(), mu4Ktmp.momentum().z(), mu4->currentState().mass());
                        J1P4tmp.SetPxPyPzE(J1Ktmp.momentum().x(), J1Ktmp.momentum().y(), J1Ktmp.momentum().z(), psi_vFit1->currentState().mass());
                        J2P4tmp.SetPxPyPzE(J2Ktmp.momentum().x(), J2Ktmp.momentum().y(), J2Ktmp.momentum().z(), psi_vFit2->currentState().mass());
                        XP4tmp.SetPxPyPzE(XKtmp.momentum().x(), XKtmp.momentum().y(), XKtmp.momentum().z(), X_vFit->currentState().mass());
                        Muon1Pt_[npairs] = mu1P4tmp.Pt();
                        Muon2Pt_[npairs] = mu2P4tmp.Pt();
                        Muon3Pt_[npairs] = mu3P4tmp.Pt();
                        Muon4Pt_[npairs] = mu4P4tmp.Pt();
                        Muon1Eta_[npairs] = mu1P4tmp.Eta();
                        Muon2Eta_[npairs] = mu2P4tmp.Eta();
                        Muon3Eta_[npairs] = mu3P4tmp.Eta();
                        Muon4Eta_[npairs] = mu4P4tmp.Eta();
                        Muon1Phi_[npairs] = mu1P4tmp.Phi();
                        Muon2Phi_[npairs] = mu2P4tmp.Phi();
                        Muon3Phi_[npairs] = mu3P4tmp.Phi();
                        Muon4Phi_[npairs] = mu4P4tmp.Phi();
                        Muon1Mass_[npairs] = mu1->currentState().mass();
                        Muon2Mass_[npairs] = mu3->currentState().mass();
                        Muon3Mass_[npairs] = mu2->currentState().mass();
                        Muon4Mass_[npairs] = mu4->currentState().mass();
                        J1normchi2_[npairs] = psi_vFit1->chiSquared() / psi_vFit1->degreesOfFreedom();
                        J2normchi2_[npairs] = psi_vFit2->chiSquared() / psi_vFit2->degreesOfFreedom();
                        J1Pt_[npairs] = J1P4tmp.Pt();
                        J2Pt_[npairs] = J2P4tmp.Pt();
                        J1Eta_[npairs] = J1P4tmp.Eta();
                        J2Eta_[npairs] = J2P4tmp.Eta();
                        J1Phi_[npairs] = J1P4tmp.Phi();
                        J2Phi_[npairs] = J2P4tmp.Phi();
                        J1Mass_[npairs] = psi_vFit1->currentState().mass();
                        J2Mass_[npairs] = psi_vFit2->currentState().mass();
                        XpT_[npairs] = XP4tmp.Pt();
                        Xeta_[npairs] = XP4tmp.Eta();
                        Xphi_[npairs] = XP4tmp.Phi();
                        Xmass_[npairs] = X_vFit->currentState().mass();
                        Xnormchi2_[npairs] = (double) X_vFit->chiSquared() / (double) X_vFit->degreesOfFreedom();
                     }
                     npairs++;
                     // mu4pair[1] = true;
                  }
               }

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
                     continue;

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
                        continue;
                     if (!psi_vFit_vertex1->vertexIsValid() || !psi_vFit_vertex2->vertexIsValid() || !X_vFit_vertex->vertexIsValid())
                        continue;

                     //KinematicParameters Jpara1 = psi_vFit1->currentState().kinematicParameters();
                     //KinematicParameters Jpara2 = psi_vFit2->currentState().kinematicParameters();
                     //KinematicParameters Xpara = X_vFit->currentState().kinematicParameters();

                     XVertexFitTree->movePointerToTheFirstChild();
                     RefCountedKinematicParticle mu1 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu4 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu2 = XVertexFitTree->currentParticle();
                     XVertexFitTree->movePointerToTheNextChild();
                     RefCountedKinematicParticle mu3 = XVertexFitTree->currentParticle();

                     double vProb1 = ChiSquaredProbability((double)(psi_vFit_vertex1->chiSquared()), (double)(psi_vFit_vertex1->degreesOfFreedom()));
                     double vProb2 = ChiSquaredProbability((double)(psi_vFit_vertex2->chiSquared()), (double)(psi_vFit_vertex2->degreesOfFreedom()));
                     //double XvProb = ChiSquaredProbability((double)(X_vFit_vertex->chiSquared()), (double)(X_vFit_vertex->degreesOfFreedom()));
                     if (vProb1 < JvPorbcut || vProb2 < JvPorbcut)
                        continue;
                     if (mu1->currentState().mass() <= 0 || mu2->currentState().mass() <= 0 || mu3->currentState().mass() <= 0 || mu4->currentState().mass() <= 0 || psi_vFit1->currentState().mass() <= 0 || psi_vFit2->currentState().mass() <= 0 || X_vFit->currentState().mass() <= 0)
                        continue;

                     if (npairs < 36)
                     {
                        TLorentzVector mu1P4tmp, mu2P4tmp, mu3P4tmp, mu4P4tmp, J1P4tmp, J2P4tmp, XP4tmp;
                        KinematicParameters mu1Ktmp, mu2Ktmp, mu3Ktmp, mu4Ktmp, J1Ktmp, J2Ktmp, XKtmp;
                        mu1Ktmp = mu1->currentState().kinematicParameters();
                        mu2Ktmp = mu2->currentState().kinematicParameters();
                        mu3Ktmp = mu3->currentState().kinematicParameters();
                        mu4Ktmp = mu4->currentState().kinematicParameters();
                        J1Ktmp = psi_vFit1->currentState().kinematicParameters();
                        J2Ktmp = psi_vFit2->currentState().kinematicParameters();
                        XKtmp = X_vFit->currentState().kinematicParameters();
                        Muon1Charge_[npairs] = iMuon1->charge();
                        Muon2Charge_[npairs] = iMuon4->charge();
                        Muon3Charge_[npairs] = iMuon2->charge();
                        Muon4Charge_[npairs] = iMuon3->charge();
                        mu1P4tmp.SetPxPyPzE(mu1Ktmp.momentum().x(), mu1Ktmp.momentum().y(), mu1Ktmp.momentum().z(), mu1->currentState().mass());
                        mu2P4tmp.SetPxPyPzE(mu4Ktmp.momentum().x(), mu4Ktmp.momentum().y(), mu4Ktmp.momentum().z(), mu4->currentState().mass());
                        mu3P4tmp.SetPxPyPzE(mu2Ktmp.momentum().x(), mu2Ktmp.momentum().y(), mu2Ktmp.momentum().z(), mu2->currentState().mass());
                        mu4P4tmp.SetPxPyPzE(mu3Ktmp.momentum().x(), mu3Ktmp.momentum().y(), mu3Ktmp.momentum().z(), mu3->currentState().mass());
                        J1P4tmp.SetPxPyPzE(J1Ktmp.momentum().x(), J1Ktmp.momentum().y(), J1Ktmp.momentum().z(), psi_vFit1->currentState().mass());
                        J2P4tmp.SetPxPyPzE(J2Ktmp.momentum().x(), J2Ktmp.momentum().y(), J2Ktmp.momentum().z(), psi_vFit2->currentState().mass());
                        XP4tmp.SetPxPyPzE(XKtmp.momentum().x(), XKtmp.momentum().y(), XKtmp.momentum().z(), X_vFit->currentState().mass());
                        Muon1Pt_[npairs] = mu1P4tmp.Pt();
                        Muon2Pt_[npairs] = mu2P4tmp.Pt();
                        Muon3Pt_[npairs] = mu3P4tmp.Pt();
                        Muon4Pt_[npairs] = mu4P4tmp.Pt();
                        Muon1Eta_[npairs] = mu1P4tmp.Eta();
                        Muon2Eta_[npairs] = mu2P4tmp.Eta();
                        Muon3Eta_[npairs] = mu3P4tmp.Eta();
                        Muon4Eta_[npairs] = mu4P4tmp.Eta();
                        Muon1Phi_[npairs] = mu1P4tmp.Phi();
                        Muon2Phi_[npairs] = mu2P4tmp.Phi();
                        Muon3Phi_[npairs] = mu3P4tmp.Phi();
                        Muon4Phi_[npairs] = mu4P4tmp.Phi();
                        Muon1Mass_[npairs] = mu1->currentState().mass();
                        Muon2Mass_[npairs] = mu4->currentState().mass();
                        Muon3Mass_[npairs] = mu2->currentState().mass();
                        Muon4Mass_[npairs] = mu3->currentState().mass();
                        J1normchi2_[npairs] = psi_vFit1->chiSquared() / psi_vFit1->degreesOfFreedom();
                        J2normchi2_[npairs] = psi_vFit2->chiSquared() / psi_vFit2->degreesOfFreedom();
                        J1Pt_[npairs] = J1P4tmp.Pt();
                        J2Pt_[npairs] = J2P4tmp.Pt();
                        J1Eta_[npairs] = J1P4tmp.Eta();
                        J2Eta_[npairs] = J2P4tmp.Eta();
                        J1Phi_[npairs] = J1P4tmp.Phi();
                        J2Phi_[npairs] = J2P4tmp.Phi();
                        J1Mass_[npairs] = psi_vFit1->currentState().mass();
                        J2Mass_[npairs] = psi_vFit2->currentState().mass();
                        XpT_[npairs] = XP4tmp.Pt();
                        Xeta_[npairs] = XP4tmp.Eta();
                        Xphi_[npairs] = XP4tmp.Phi();
                        Xmass_[npairs] = X_vFit->currentState().mass();
                        Xnormchi2_[npairs] = (double) X_vFit->chiSquared() / (double) X_vFit->degreesOfFreedom();
                     }
                     npairs++;
                     // mu4pair[2] = true;
                  }
               }
            }
         }
      }
   }

   if (npairs > 0) 
   {
      X4muTree->Fill();
      //std::cout << "npairs: " << npairs << " XMass: " << Xmass_[0] << " J1Mass: " << J1Mass_[0] << " J2Mass: " << J2Mass_[0] << " Muon1Mass: " << Muon1Mass_[0] << " Muon2Mass: " << Muon2Mass_[0] << " Muon3Mass: " << Muon3Mass_[0] << " Muon4Mass: " << Muon4Mass_[0] << std::endl;
   }

   /*    std::unique_ptr<edm::ValueMap<int>> Muon1Charge(new edm::ValueMap<int>());
      std::unique_ptr<edm::ValueMap<int>> Muon2Charge(new edm::ValueMap<int>());
      std::unique_ptr<edm::ValueMap<int>> Muon3Charge(new edm::ValueMap<int>());
      std::unique_ptr<edm::ValueMap<int>> Muon4Charge(new edm::ValueMap<int>());
      std::unique_ptr<edm::ValueMap<float>> Muon1Pt(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon2Pt(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon3Pt(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon4Pt(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon1Eta(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon2Eta(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon3Eta(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon4Eta(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon1Phi(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon2Phi(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon3Phi(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Muon4Phi(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J1normchi2(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J2normchi2(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J1Pt(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J2Pt(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J1Eta(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J2Eta(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J1Phi(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J2Phi(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J1Mass(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> J2Mass(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<float>> Xmass(new edm::ValueMap<float>());
      std::unique_ptr<edm::ValueMap<int>> Xnormchi2(new edm::ValueMap<int>());

      edm::ValueMap<int>::Filler filler_Muon1Charge(*Muon1Charge);
      edm::ValueMap<int>::Filler filler_Muon2Charge(*Muon2Charge);
      edm::ValueMap<int>::Filler filler_Muon3Charge(*Muon3Charge);
      edm::ValueMap<int>::Filler filler_Muon4Charge(*Muon4Charge);
      edm::ValueMap<float>::Filler filler_Muon1Pt(*Muon1Pt);
      edm::ValueMap<float>::Filler filler_Muon2Pt(*Muon2Pt);
      edm::ValueMap<float>::Filler filler_Muon3Pt(*Muon3Pt);
      edm::ValueMap<float>::Filler filler_Muon4Pt(*Muon4Pt);
      edm::ValueMap<float>::Filler filler_Muon1Eta(*Muon1Eta);
      edm::ValueMap<float>::Filler filler_Muon2Eta(*Muon2Eta);
      edm::ValueMap<float>::Filler filler_Muon3Eta(*Muon3Eta);
      edm::ValueMap<float>::Filler filler_Muon4Eta(*Muon4Eta);
      edm::ValueMap<float>::Filler filler_Muon1Phi(*Muon1Phi);
      edm::ValueMap<float>::Filler filler_Muon2Phi(*Muon2Phi);
      edm::ValueMap<float>::Filler filler_Muon3Phi(*Muon3Phi);
      edm::ValueMap<float>::Filler filler_Muon4Phi(*Muon4Phi);
      edm::ValueMap<float>::Filler filler_J1normchi2(*J1normchi2);
      edm::ValueMap<float>::Filler filler_J2normchi2(*J2normchi2);
      edm::ValueMap<float>::Filler filler_J1Pt(*J1Pt);
      edm::ValueMap<float>::Filler filler_J2Pt(*J2Pt);
      edm::ValueMap<float>::Filler filler_J1Eta(*J1Eta);
      edm::ValueMap<float>::Filler filler_J2Eta(*J2Eta);
      edm::ValueMap<float>::Filler filler_J1Phi(*J1Phi);
      edm::ValueMap<float>::Filler filler_J2Phi(*J2Phi);
      edm::ValueMap<float>::Filler filler_J1Mass(*J1Mass);
      edm::ValueMap<float>::Filler filler_J2Mass(*J2Mass);
      edm::ValueMap<float>::Filler filler_Xmass(*Xmass);
      edm::ValueMap<int>::Filler filler_Xnormchi2(*Xnormchi2);

      filler_Muon1Charge.insert(oh, Muon1Charge_.begin(), Muon1Charge_.end());
      filler_Muon2Charge.insert(oh, Muon2Charge_.begin(), Muon2Charge_.end());
      filler_Muon3Charge.insert(oh, Muon3Charge_.begin(), Muon3Charge_.end());
      filler_Muon4Charge.insert(oh, Muon4Charge_.begin(), Muon4Charge_.end());
      filler_Muon1Pt.insert(oh, Muon1Pt_.begin(), Muon1Pt_.end());
      filler_Muon2Pt.insert(oh, Muon2Pt_.begin(), Muon2Pt_.end());
      filler_Muon3Pt.insert(oh, Muon3Pt_.begin(), Muon3Pt_.end());
      filler_Muon4Pt.insert(oh, Muon4Pt_.begin(), Muon4Pt_.end());
      filler_Muon1Eta.insert(oh, Muon1Eta_.begin(), Muon1Eta_.end());
      filler_Muon2Eta.insert(oh, Muon2Eta_.begin(), Muon2Eta_.end());
      filler_Muon3Eta.insert(oh, Muon3Eta_.begin(), Muon3Eta_.end());
      filler_Muon4Eta.insert(oh, Muon4Eta_.begin(), Muon4Eta_.end());
      filler_Muon1Phi.insert(oh, Muon1Phi_.begin(), Muon1Phi_.end());
      filler_Muon2Phi.insert(oh, Muon2Phi_.begin(), Muon2Phi_.end());
      filler_Muon3Phi.insert(oh, Muon3Phi_.begin(), Muon3Phi_.end());
      filler_Muon4Phi.insert(oh, Muon4Phi_.begin(), Muon4Phi_.end());
      filler_J1normchi2.insert(oh, J1normchi2_.begin(), J1normchi2_.end());
      filler_J2normchi2.insert(oh, J2normchi2_.begin(), J2normchi2_.end());
      filler_J1Pt.insert(oh, J1Pt_.begin(), J1Pt_.end());
      filler_J2Pt.insert(oh, J2Pt_.begin(), J2Pt_.end());
      filler_J1Eta.insert(oh, J1Eta_.begin(), J1Eta_.end());
      filler_J2Eta.insert(oh, J2Eta_.begin(), J2Eta_.end());
      filler_J1Phi.insert(oh, J1Phi_.begin(), J1Phi_.end());
      filler_J2Phi.insert(oh, J2Phi_.begin(), J2Phi_.end());
      filler_J1Mass.insert(oh, J1Mass_.begin(), J1Mass_.end());
      filler_J2Mass.insert(oh, J2Mass_.begin(), J2Mass_.end());
      filler_Xmass.insert(oh, Xmass_.begin(), Xmass_.end());
      filler_Xnormchi2.insert(oh, Xnormchi2_.begin(), Xnormchi2_.end());

      filler_Muon1Charge.fill();
      filler_Muon2Charge.fill();
      filler_Muon3Charge.fill();
      filler_Muon4Charge.fill();
      filler_Muon1Pt.fill();
      filler_Muon2Pt.fill();
      filler_Muon3Pt.fill();
      filler_Muon4Pt.fill();
      filler_Muon1Eta.fill();
      filler_Muon2Eta.fill();
      filler_Muon3Eta.fill();
      filler_Muon4Eta.fill();
      filler_Muon1Phi.fill();
      filler_Muon2Phi.fill();
      filler_Muon3Phi.fill();
      filler_Muon4Phi.fill();
      filler_J1normchi2.fill();
      filler_J2normchi2.fill();
      filler_J1Pt.fill();
      filler_J2Pt.fill();
      filler_J1Eta.fill();
      filler_J2Eta.fill();
      filler_J1Phi.fill();
      filler_J2Phi.fill();
      filler_J1Mass.fill();
      filler_J2Mass.fill();
      filler_Xmass.fill();
      filler_Xnormchi2.fill();

      iEvent.put(std::move(Muon1Charge), "Muon1Charge");
      iEvent.put(std::move(Muon2Charge), "Muon2Charge");
      iEvent.put(std::move(Muon3Charge), "Muon3Charge");
      iEvent.put(std::move(Muon4Charge), "Muon4Charge");
      iEvent.put(std::move(Muon1Pt), "Muon1Pt");
      iEvent.put(std::move(Muon2Pt), "Muon2Pt");
      iEvent.put(std::move(Muon3Pt), "Muon3Pt");
      iEvent.put(std::move(Muon4Pt), "Muon4Pt");
      iEvent.put(std::move(Muon1Eta), "Muon1Eta");
      iEvent.put(std::move(Muon2Eta), "Muon2Eta");
      iEvent.put(std::move(Muon3Eta), "Muon3Eta");
      iEvent.put(std::move(Muon4Eta), "Muon4Eta");
      iEvent.put(std::move(Muon1Phi), "Muon1Phi");
      iEvent.put(std::move(Muon2Phi), "Muon2Phi");
      iEvent.put(std::move(Muon3Phi), "Muon3Phi");
      iEvent.put(std::move(Muon4Phi), "Muon4Phi");
      iEvent.put(std::move(J1normchi2), "J1normchi2");
      iEvent.put(std::move(J2normchi2), "J2normchi2");
      iEvent.put(std::move(J1Pt), "J1Pt");
      iEvent.put(std::move(J2Pt), "J2Pt");
      iEvent.put(std::move(J1Eta), "J1Eta");
      iEvent.put(std::move(J2Eta), "J2Eta");
      iEvent.put(std::move(J1Phi), "J1Phi");
      iEvent.put(std::move(J2Phi), "J2Phi");
      iEvent.put(std::move(J1Mass), "J1Mass");
      iEvent.put(std::move(J2Mass), "J2Mass");
      iEvent.put(std::move(Xmass), "Xmass");
      iEvent.put(std::move(Xnormchi2), "Xnormchi2"); */
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
   }
}

void X4muSecondaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
   edm::ParameterSetDescription desc;
   desc.add<edm::InputTag>("recoMuon", edm::InputTag("recoMuon"));
   desc.add<edm::InputTag>("recoTrack", edm::InputTag("recoTrack"));

   descriptions.add("X4muSecondaryVertexProducer", desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(X4muSecondaryVertexProducer);