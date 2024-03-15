#define EDM_ML_DEBUG
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/Math/interface/angle_units.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDSectorForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetSector.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "HepMC/GenRanges.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "MTDHit.h"

class MtdTracksValidation : public DQMEDAnalyzer {
public:
  explicit MtdTracksValidation(const edm::ParameterSet&);
  ~MtdTracksValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const std::pair<bool, bool> checkAcceptance(
      const reco::Track&, const edm::Event&, const edm::EventSetup&, size_t&, float&, float&, float&, float&);

  const bool mvaGenSel(const HepMC::GenParticle&, const float&);
  const bool mvaTPSel(const TrackingParticle&);
  const bool mvaRecSel(const reco::TrackBase&, const reco::Vertex&, const double&, const double&);
  const bool mvaGenRecMatch(const HepMC::GenParticle&, const double&, const reco::TrackBase&, const bool&);
  const edm::Ref<std::vector<TrackingParticle>>* getMatchedTP(const reco::TrackBaseRef&);

  const unsigned long int uniqueId(const uint32_t x, const EncodedEventId& y) {
    const uint64_t a = static_cast<uint64_t>(x);
    const uint64_t b = static_cast<uint64_t>(y.rawId());

    if (x < y.rawId())
      return (b << 32) | a;
    else
      return (a << 32) | b;
  }

  bool isETL(const double eta) const { return (std::abs(eta) > trackMinEtlEta_) && (std::abs(eta) < trackMaxEtlEta_); }

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinPt_;
  const float trackMaxBtlEta_;
  const float trackMinEtlEta_;
  const float trackMaxEtlEta_;

  static constexpr double etacutGEN_ = 4.;               // |eta| < 4;
  static constexpr double etacutREC_ = 3.;               // |eta| < 3;
  static constexpr double pTcut_ = 0.7;                  // PT > 0.7 GeV
  static constexpr double deltaZcut_ = 0.1;              // dz separation 1 mm
  static constexpr double deltaPTcut_ = 0.05;            // dPT < 5%
  static constexpr double deltaDRcut_ = 0.03;            // DeltaR separation
  static constexpr double depositBTLthreshold_ = 1;      // threshold for energy deposit in BTL cell [MeV]
  static constexpr double depositETLthreshold_ = 0.001;  // threshold for energy deposit in ETL cell [MeV]
  static constexpr double rBTL_ = 110.0;
  static constexpr double zETL_ = 290.0;
  static constexpr double etaMatchCut_ = 0.05;
  static constexpr double cluDRradius_ = 0.05;  // to cluster rechits around extrapolated track

  const reco::RecoToSimCollection* r2s_;
  const reco::SimToRecoCollection* s2r_;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<edm::HepMCProduct> HepMCProductToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<CrossingFrame<PSimHit>> btlSimHitsToken_;
  edm::EDGetTokenT<CrossingFrame<PSimHit>> etlSimHitsToken_;
  edm::EDGetTokenT<FTLRecHitCollection> btlRecHitsToken_;
  edm::EDGetTokenT<FTLRecHitCollection> etlRecHitsToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> TofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> TofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> TofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> BetaToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> DtSignPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> DtSignKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> DtSignPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ProbPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ProbKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ProbPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> TimeChisqPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> TimeChisqKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> TimeChisqPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> BtlMatchTimeChi2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> EtlMatchTimeChi2Token_;

  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;

  edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeoToken_;
  edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;
  edm::ESGetToken<MTDDetLayerGeometry, MTDRecoGeometryRecord> mtdlayerToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> builderToken_;
  edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particleTableToken_;

  MonitorElement* meBTLTrackRPTime_;
  MonitorElement* meBTLTrackEffEtaTot_;
  MonitorElement* meBTLTrackEffPhiTot_;
  MonitorElement* meBTLTrackEffPtTot_;
  MonitorElement* meBTLTrackEffEtaMtd_;
  MonitorElement* meBTLTrackEffPhiMtd_;
  MonitorElement* meBTLTrackEffPtMtd_;
  MonitorElement* meBTLTrackPtRes_;

  MonitorElement* meETLTrackRPTime_;
  MonitorElement* meETLTrackEffEtaTot_[2];
  MonitorElement* meETLTrackEffPhiTot_[2];
  MonitorElement* meETLTrackEffPtTot_[2];
  MonitorElement* meETLTrackEffEtaMtd_[2];
  MonitorElement* meETLTrackEffPhiMtd_[2];
  MonitorElement* meETLTrackEffPtMtd_[2];
  MonitorElement* meETLTrackEffEta2Mtd_[2];
  MonitorElement* meETLTrackEffPhi2Mtd_[2];
  MonitorElement* meETLTrackEffPt2Mtd_[2];
  MonitorElement* meETLTrackPtRes_;

  MonitorElement* meTracktmtd_;
  MonitorElement* meTrackt0Src_;
  MonitorElement* meTrackSigmat0Src_;
  MonitorElement* meTrackt0Pid_;
  MonitorElement* meTrackSigmat0Pid_;
  MonitorElement* meTrackt0SafePid_;
  MonitorElement* meTrackSigmat0SafePid_;

  MonitorElement* meTrackDtSign_[3]; 
  MonitorElement* meTrackProb_[3];
  MonitorElement* meTrackTimeChi2_PionTrack_[3]; 
  MonitorElement* meTrackTimeChi2_KaonTrack_[3];
  MonitorElement* meTrackTimeChi2_ProtonTrack_[3];
  MonitorElement* meTrackBtlMatchChi2_[3];
  MonitorElement* meTrackEtlMatchChi2_[3];
  MonitorElement* meRatioMatchedTracksVsP_BTL_num_;
  MonitorElement* meRatioMatchedTracksVsP_BTL_den_;
  MonitorElement* meRatioMatchedTracksVsP_ETL_num_;
  MonitorElement* meRatioMatchedTracksVsP_ETL_den_;

  MonitorElement* meTruePi_;
  MonitorElement* meTrueK_;
  MonitorElement* meTrueP_;

  MonitorElement* mePIDPi_;
  MonitorElement* mePIDK_;
  MonitorElement* mePIDP_;

  MonitorElement* meTruePiNoPID_;
  MonitorElement* meTruePiAsPi_;
  MonitorElement* meTruePiAsK_;
  MonitorElement* meTruePiAsP_;

  MonitorElement* meTrueKNoPID_;
  MonitorElement* meTrueKAsPi_;
  MonitorElement* meTrueKAsK_;
  MonitorElement* meTrueKAsP_;

  MonitorElement* meTruePNoPID_;
  MonitorElement* meTruePAsPi_;
  MonitorElement* meTruePAsK_;
  MonitorElement* meTruePAsP_;  

  MonitorElement* meTOF_;
  MonitorElement* meTmtd_;
  MonitorElement* meInvBeta_[3];
  MonitorElement* meInvBetaVsP_;

  MonitorElement* meDeltaChisq_truePion_vsKaon_;
  MonitorElement* meDeltaChisq_truePion_vsProton_;
  MonitorElement* meDeltaChisq_trueKaon_vsPion_;
  MonitorElement* meDeltaChisq_trueKaon_vsProton_;
  MonitorElement* meDeltaChisq_trueProton_vsPion_;
  MonitorElement* meDeltaChisq_trueProton_vsKaon_;

  MonitorElement* meTrackNumHits_;
  MonitorElement* meTrackNumHitsNT_;
  MonitorElement* meTrackMVAQual_;
  MonitorElement* meTrackPathLenghtvsEta_;

  MonitorElement* meTrackPtTot_;
  MonitorElement* meMVATrackEffPtTot_;
  MonitorElement* meMVATrackMatchedEffPtTot_;
  MonitorElement* meMVATrackMatchedEffPtMtd_;
  MonitorElement* meExtraPtMtd_;
  MonitorElement* meExtraPtEtl2Mtd_;

  MonitorElement* meBTLTrackMatchedTPPtResMtd_;
  MonitorElement* meETLTrackMatchedTPPtResMtd_;
  MonitorElement* meETLTrackMatchedTP2PtResMtd_;
  MonitorElement* meBTLTrackMatchedTPPtRatioGen_;
  MonitorElement* meETLTrackMatchedTPPtRatioGen_;
  MonitorElement* meETLTrackMatchedTP2PtRatioGen_;
  MonitorElement* meBTLTrackMatchedTPPtRatioMtd_;
  MonitorElement* meETLTrackMatchedTPPtRatioMtd_;
  MonitorElement* meETLTrackMatchedTP2PtRatioMtd_;
  MonitorElement* meBTLTrackMatchedTPPtResvsPtMtd_;
  MonitorElement* meETLTrackMatchedTPPtResvsPtMtd_;
  MonitorElement* meETLTrackMatchedTP2PtResvsPtMtd_;
  MonitorElement* meBTLTrackMatchedTPDPtvsPtGen_;
  MonitorElement* meETLTrackMatchedTPDPtvsPtGen_;
  MonitorElement* meETLTrackMatchedTP2DPtvsPtGen_;
  MonitorElement* meBTLTrackMatchedTPDPtvsPtMtd_;
  MonitorElement* meETLTrackMatchedTPDPtvsPtMtd_;
  MonitorElement* meETLTrackMatchedTP2DPtvsPtMtd_;

  MonitorElement* meTrackMatchedTPEffPtTot_;
  MonitorElement* meTrackMatchedTPEffPtMtd_;
  MonitorElement* meTrackMatchedTPEffPtEtl2Mtd_;
  MonitorElement* meTrackMatchedTPmtdEffPtTot_;
  MonitorElement* meTrackMatchedTPmtdEffPtMtd_;
  MonitorElement* meTrackEtaTot_;
  MonitorElement* meMVATrackEffEtaTot_;
  MonitorElement* meMVATrackMatchedEffEtaTot_;
  MonitorElement* meMVATrackMatchedEffEtaMtd_;
  MonitorElement* meExtraEtaMtd_;
  MonitorElement* meExtraEtaEtl2Mtd_;
  MonitorElement* meTrackMatchedTPEffEtaTot_;
  MonitorElement* meTrackMatchedTPEffEtaMtd_;
  MonitorElement* meTrackMatchedTPEffEtaEtl2Mtd_;
  MonitorElement* meTrackMatchedTPmtdEffEtaTot_;
  MonitorElement* meTrackMatchedTPmtdEffEtaMtd_;
  MonitorElement* meMVATrackResTot_;
  MonitorElement* meMVATrackPullTot_;
  MonitorElement* meMVATrackZposResTot_;

  MonitorElement* meExtraPhiAtBTL_;
  MonitorElement* meExtraPhiAtBTLmatched_;
  MonitorElement* meExtraBTLeneInCone_;
  MonitorElement* meExtraMTDfailExtenderEta_;
  MonitorElement* meExtraMTDfailExtenderPt_;
};

// ------------ constructor and destructor --------------
MtdTracksValidation::MtdTracksValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinPt_(iConfig.getParameter<double>("trackMinimumPt")),
      trackMaxBtlEta_(iConfig.getParameter<double>("trackMaximumBtlEta")),
      trackMinEtlEta_(iConfig.getParameter<double>("trackMinimumEtlEta")),
      trackMaxEtlEta_(iConfig.getParameter<double>("trackMaximumEtlEta")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTagV"));
  HepMCProductToken_ = consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("inputTagH"));
  trackingParticleCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  simToRecoAssociationToken_ =
      consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  recoToSimAssociationToken_ =
      consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  btlSimHitsToken_ = consumes<CrossingFrame<PSimHit>>(iConfig.getParameter<edm::InputTag>("btlSimHits"));
  etlSimHitsToken_ = consumes<CrossingFrame<PSimHit>>(iConfig.getParameter<edm::InputTag>("etlSimHits"));
  btlRecHitsToken_ = consumes<FTLRecHitCollection>(iConfig.getParameter<edm::InputTag>("btlRecHits"));
  etlRecHitsToken_ = consumes<FTLRecHitCollection>(iConfig.getParameter<edm::InputTag>("etlRecHits"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));

  TofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofPi"));
  TofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofK"));
  TofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofP"));
  BetaToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("beta"));
  DtSignPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dtSignPi"));
  DtSignKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dtSignK"));
  DtSignPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dtSignP"));
  ProbPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPi"));
  ProbKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probK"));
  ProbPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probP"));
  TimeChisqPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("timeChisqPi"));
  TimeChisqKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("timeChisqK"));
  TimeChisqPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("timeChisqP"));
  BtlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2"));
  EtlMatchTimeChi2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2"));

  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  mtdgeoToken_ = esConsumes<MTDGeometry, MTDDigiGeometryRecord>();
  mtdtopoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>();
  mtdlayerToken_ = esConsumes<MTDDetLayerGeometry, MTDRecoGeometryRecord>();
  magfieldToken_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  builderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  particleTableToken_ = esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>();
}

MtdTracksValidation::~MtdTracksValidation() {}

// ------------ method called for each event  ------------
void MtdTracksValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));
  auto RecVertexHandle = makeValid(iEvent.getHandle(RecVertexToken_));

  std::unordered_map<uint32_t, MTDHit> m_btlHits;
  std::unordered_map<uint32_t, MTDHit> m_etlHits;
  std::unordered_map<uint32_t, std::set<unsigned long int>> m_btlTrkPerCell;
  std::unordered_map<uint32_t, std::set<unsigned long int>> m_etlTrkPerCell;
  std::map<TrackingParticleRef, std::vector<uint32_t>> m_tp2detid;

  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  const auto& t0Src = iEvent.get(t0SrcToken_);
  const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& t0Safe = iEvent.get(t0SafePidToken_);
  const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);

  const auto& tofPi = iEvent.get(TofPiToken_);
  const auto& tofK = iEvent.get(TofKToken_);
  const auto& tofP = iEvent.get(TofPToken_);
  const auto& beta = iEvent.get(BetaToken_);

  const auto& dtSignPi = iEvent.get(DtSignPiToken_);
  const auto& dtSignK = iEvent.get(DtSignKToken_);
  const auto& dtSignP = iEvent.get(DtSignPToken_);

  const auto& probPi  = iEvent.get(ProbPiToken_);
  const auto& probK   = iEvent.get(ProbKToken_);
  const auto& probP   = iEvent.get(ProbPToken_);

  const auto& timeChisqPi = iEvent.get(TimeChisqPiToken_);
  const auto& timeChisqK = iEvent.get(TimeChisqKToken_);
  const auto& timeChisqP = iEvent.get(TimeChisqPToken_);

  const auto& btlMatchTimeChi2 = iEvent.get(BtlMatchTimeChi2Token_);
  const auto& etlMatchTimeChi2 = iEvent.get(EtlMatchTimeChi2Token_);

  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& trackAssoc = iEvent.get(trackAssocToken_);
  const auto& pathLength = iEvent.get(pathLengthToken_);

  const auto& primRecoVtx = *(RecVertexHandle.product()->begin());

  // generator level information (HepMC format)
  auto GenEventHandle = makeValid(iEvent.getHandle(HepMCProductToken_));
  const HepMC::GenEvent* mc = GenEventHandle->GetEvent();
  double zsim = convertMmToCm((*(mc->vertices_begin()))->position().z());
  double tsim = (*(mc->vertices_begin()))->position().t() * CLHEP::mm / CLHEP::c_light;

  auto pdt = iSetup.getHandle(particleTableToken_);
  const HepPDT::ParticleDataTable* pdTable = pdt.product();

  auto simToRecoH = makeValid(iEvent.getHandle(simToRecoAssociationToken_));
  s2r_ = simToRecoH.product();

  auto recoToSimH = makeValid(iEvent.getHandle(recoToSimAssociationToken_));
  r2s_ = recoToSimH.product();

  //Fill maps with simhits accumulated per DetId

  auto btlSimHitsHandle = makeValid(iEvent.getHandle(btlSimHitsToken_));
  MixCollection<PSimHit> btlSimHits(btlSimHitsHandle.product());
  for (auto const& simHit : btlSimHits) {
    if (simHit.tof() < 0 || simHit.tof() > 25.)
      continue;
    DetId id = simHit.detUnitId();
    auto const thisHId = uniqueId(simHit.trackId(), simHit.eventId());
    m_btlTrkPerCell[id.rawId()].insert(thisHId);
    auto simHitIt = m_btlHits.emplace(id.rawId(), MTDHit()).first;
    // --- Accumulate the energy (in MeV) of SIM hits in the same detector cell
    (simHitIt->second).energy += convertUnitsTo(0.001_MeV, simHit.energyLoss());
  }

  auto etlSimHitsHandle = makeValid(iEvent.getHandle(etlSimHitsToken_));
  MixCollection<PSimHit> etlSimHits(etlSimHitsHandle.product());
  for (auto const& simHit : etlSimHits) {
    if (simHit.tof() < 0 || simHit.tof() > 25.) {
      continue;
    }
    DetId id = simHit.detUnitId();
    auto const thisHId = uniqueId(simHit.trackId(), simHit.eventId());
    m_etlTrkPerCell[id.rawId()].insert(thisHId);
    auto simHitIt = m_etlHits.emplace(id.rawId(), MTDHit()).first;
    // --- Accumulate the energy (in MeV) of SIM hits in the same detector cell
    (simHitIt->second).energy += convertUnitsTo(0.001_MeV, simHit.energyLoss());
  }

  //Fill map of DetId per ref to TP

  auto tpHandle = makeValid(iEvent.getHandle(trackingParticleCollectionToken_));
  TrackingParticleCollection tpColl = *(tpHandle.product());
  size_t tpindex(0);
  for (auto tp = tpColl.begin(); tp != tpColl.end(); tp++, ++tpindex) {
    TrackingParticleRef tpref(iEvent.getHandle(trackingParticleCollectionToken_), tpindex);
    if (tp->eventId().bunchCrossing() == 0 && tp->eventId().event() == 0) {
      if (!mvaTPSel(*tp))
        continue;
      for (const auto& simTrk : tp->g4Tracks()) {
        auto const thisTId = uniqueId(simTrk.trackId(), simTrk.eventId());
        for (auto const& cell : m_btlTrkPerCell) {
          if (m_btlHits[cell.first].energy < depositBTLthreshold_) {
            continue;
          }
          for (auto const& simtrack : cell.second) {
            if (thisTId == simtrack) {
              m_tp2detid[tpref].emplace_back(cell.first);
              break;
            }
          }
        }
        for (auto const& cell : m_etlTrkPerCell) {
          if (m_etlHits[cell.first].energy < depositETLthreshold_) {
            continue;
          }
          for (auto const& simtrack : cell.second) {
            if (thisTId == simtrack) {
              m_tp2detid[tpref].emplace_back(cell.first);
              break;
            }
          }
        }
      }
    }
  }

  unsigned int index = 0;

  // flag to select events with reco vertex close to true simulated primary vertex, or PV fake (particle guns)
  const bool isGoodVtx = std::abs(primRecoVtx.z() - zsim) < deltaZcut_ || primRecoVtx.isFake();

  // --- Loop over all RECO tracks ---
  for (const auto& trackGen : *GenRecTrackHandle) {
    const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
    index++;

    if (trackAssoc[trackref] == -1) {
      LogInfo("mtdTracks") << "Extended track not associated";
      continue;
    }

    const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(RecTrackToken_), trackAssoc[trackref]);
    const reco::Track& track = *mtdTrackref;

    bool isBTL = false;
    bool isETL = false;
    bool twoETLdiscs = false;
    bool noCrack = std::abs(trackGen.eta()) < trackMaxBtlEta_ || std::abs(trackGen.eta()) > trackMinEtlEta_;

    if (track.pt() >= trackMinPt_ && std::abs(track.eta()) <= trackMaxEtlEta_) {
      meTracktmtd_->Fill(tMtd[trackref]);
      if (std::round(SigmatMtd[trackref] - Sigmat0Pid[trackref]) != 0) {
        LogWarning("mtdTracks")
            << "TimeError associated to refitted track is different from TimeError stored in tofPID "
               "sigmat0 ValueMap: this should not happen";
      }

      meTrackt0Src_->Fill(t0Src[trackref]);
      meTrackSigmat0Src_->Fill(Sigmat0Src[trackref]);

      meTrackt0Pid_->Fill(t0Pid[trackref]);
      meTrackSigmat0Pid_->Fill(Sigmat0Pid[trackref]);
      meTrackt0SafePid_->Fill(t0Safe[trackref]);
      meTrackSigmat0SafePid_->Fill(Sigmat0Safe[trackref]);
      meTrackMVAQual_->Fill(mtdQualMVA[trackref]);

      // ---------------------------------------------
      // ----------- SIGMA(TOF) VALIDATION -----------
      // ---------------------------------------------
      
      double tof[3] = {tofPi[trackref], tofK[trackref], tofP[trackref]};

      double dtSigns[3] = {dtSignPi[trackref], dtSignK[trackref], dtSignP[trackref]};
      for (int i = 0; i < 3; i++) meTrackDtSign_[i]->Fill(dtSigns[i]);

      double probs[3] = {probPi[trackref], probK[trackref], probP[trackref]};
      for(int i = 0; i < 3; i++) meTrackProb_[i]->Fill(probs[i]);

      double timeChisqs[3] = {timeChisqPi[trackref], timeChisqK[trackref], timeChisqP[trackref]};

      int imax = std::max_element(probs, probs + 3) - probs; //find most likely hypothesis
      string pNames[3] = {"Pi", "K", "P"};

      //--- defining PID status
      unsigned int noPIDtype = 0;
      if (probPi[trackref] == -1) {
        noPIDtype = 1;
      } else if (isnan(probPi[trackref])) {
        noPIDtype = 2;
      } else if (probPi[trackref] == 1 && probK[trackref] == 0 && probP[trackref] == 0) {
        noPIDtype = 3;
      }
      bool noPID = noPIDtype > 0;
      bool isPi = !noPID && 1. - probPi[trackref] < 0.75;
      bool isK = !noPID && !isPi && probK[trackref] > probP[trackref];
      bool isP = !noPID && !isPi && !isK;
      //--- end PID status      

      const reco::TrackBaseRef mytrkrefb(trackref);
      auto mytp_info = getMatchedTP(mytrkrefb);

      // --- PID PURITY/EFFICIENCY STUDIES

      if(imax == 0){ //if track identified as pion
        mePIDPi_->Fill(track.p());
      } else if(imax == 1) { //track ID as kaon
        mePIDK_->Fill(track.p());
      } else { // track ID as proton
        mePIDP_->Fill(track.p());
      }

      if(mytp_info && Sigmat0Safe[trackref] > 0.){
        if (std::abs((*mytp_info)->pdgId()) == 211) {
          meTruePi_->Fill(track.p());

          meDeltaChisq_truePion_vsKaon_->Fill(timeChisqs[0] - timeChisqs[1]);
          meDeltaChisq_truePion_vsProton_->Fill(timeChisqs[0] - timeChisqs[2]);

          for (int i = 0; i < 3; i++) {
            meTrackTimeChi2_PionTrack_[i]->Fill(timeChisqs[i]);
          }

          if (noPID) {
            meTruePiNoPID_->Fill(track.p());
          } else if (isPi) {
            meTruePiAsPi_->Fill(track.p());
          } else if (isK) {
            meTruePiAsK_->Fill(track.p());
          } else if (isP) {
            meTruePiAsP_->Fill(track.p());
          }

        } else if (std::abs((*mytp_info)->pdgId()) == 321) {
          meTrueK_->Fill(track.p());

          meDeltaChisq_trueKaon_vsPion_->Fill(timeChisqs[1] - timeChisqs[0]);
          meDeltaChisq_trueKaon_vsProton_->Fill(timeChisqs[1] - timeChisqs[2]);

          for (int i = 0; i < 3; i++) {
            meTrackTimeChi2_KaonTrack_[i]->Fill(timeChisqs[i]);
          }

          if (noPID) {
            meTrueKNoPID_->Fill(track.p());
          } else if (isPi) {
            meTrueKAsPi_->Fill(track.p());
          } else if (isK) {
            meTrueKAsK_->Fill(track.p());
          } else if (isP) {
            meTrueKAsP_->Fill(track.p());
          }

        } else if (std::abs((*mytp_info)->pdgId()) == 2212) {
          meTrueP_->Fill(track.p());

          meDeltaChisq_trueProton_vsPion_->Fill(timeChisqs[2] - timeChisqs[0]);
          meDeltaChisq_trueProton_vsKaon_->Fill(timeChisqs[2] - timeChisqs[1]);

          for (int i = 0; i < 3; i++) {
            meTrackTimeChi2_ProtonTrack_[i]->Fill(timeChisqs[i]);
          }

          if (noPID) {
            meTruePNoPID_->Fill(track.p());
          } else if (isPi) {
            meTruePAsPi_->Fill(track.p());
          } else if (isK) {
            meTruePAsK_->Fill(track.p());
          } else if (isP) {
            meTruePAsP_->Fill(track.p());
          }

        }
      }

      // MTD hit match
      meTrackBtlMatchChi2_[imax]->Fill(btlMatchTimeChi2[trackref]);
      meTrackEtlMatchChi2_[imax]->Fill(etlMatchTimeChi2[trackref]);
      if(std::abs(track.eta()) < trackMaxBtlEta_){
        meRatioMatchedTracksVsP_BTL_den_->Fill(track.p());
        if(Sigmat0Safe[trackref] > 0.) meRatioMatchedTracksVsP_BTL_num_->Fill(track.p());
      } else if(std::abs(track.eta()) > trackMinEtlEta_ && std::abs(track.eta()) < trackMaxEtlEta_){
        meRatioMatchedTracksVsP_ETL_den_->Fill(track.p());
        if(Sigmat0Safe[trackref] > 0.) meRatioMatchedTracksVsP_ETL_num_->Fill(track.p());
      }

      // HSCP studies
      if(Sigmat0Safe[trackref] > 0.) { // if time info attached to track
        meTOF_->Fill(tof[imax]);
        meTmtd_->Fill(tMtd[trackref]);
        meInvBeta_[imax]->Fill(1/beta[trackref]);
        meInvBetaVsP_->Fill(track.p(), 1/beta[trackref]);
      }

      // --------------------------------------------------------------
      // --------------------------------------------------------------


      meTrackPathLenghtvsEta_->Fill(std::abs(track.eta()), pathLength[trackref]);

      if (std::abs(track.eta()) < trackMaxBtlEta_) {
        // --- all BTL tracks (with and without hit in MTD) ---
        meBTLTrackEffEtaTot_->Fill(track.eta());
        meBTLTrackEffPhiTot_->Fill(track.phi());
        meBTLTrackEffPtTot_->Fill(track.pt());

        bool MTDBtl = false;
        int numMTDBtlvalidhits = 0;
        for (const auto hit : track.recHits()) {
          if (hit->isValid() == false)
            continue;
          MTDDetId Hit = hit->geographicalId();
          if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 1)) {
            MTDBtl = true;
            numMTDBtlvalidhits++;
          }
        }
        meTrackNumHits_->Fill(numMTDBtlvalidhits);

        // --- keeping only tracks with last hit in MTD ---
        if (MTDBtl == true) {
          isBTL = true;
          meBTLTrackEffEtaMtd_->Fill(track.eta());
          meBTLTrackEffPhiMtd_->Fill(track.phi());
          meBTLTrackEffPtMtd_->Fill(track.pt());
          meBTLTrackRPTime_->Fill(track.t0());
          meBTLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
        }
        if (isBTL && Sigmat0Safe[trackref] < 0.) {
          meTrackNumHitsNT_->Fill(numMTDBtlvalidhits);
        }
      }  //loop over (geometrical) BTL tracks

      else {
        // --- all ETL tracks (with and without hit in MTD) ---
        if ((track.eta() < -trackMinEtlEta_) && (track.eta() > -trackMaxEtlEta_)) {
          meETLTrackEffEtaTot_[0]->Fill(track.eta());
          meETLTrackEffPhiTot_[0]->Fill(track.phi());
          meETLTrackEffPtTot_[0]->Fill(track.pt());
        }

        if ((track.eta() > trackMinEtlEta_) && (track.eta() < trackMaxEtlEta_)) {
          meETLTrackEffEtaTot_[1]->Fill(track.eta());
          meETLTrackEffPhiTot_[1]->Fill(track.phi());
          meETLTrackEffPtTot_[1]->Fill(track.pt());
        }

        bool MTDEtlZnegD1 = false;
        bool MTDEtlZnegD2 = false;
        bool MTDEtlZposD1 = false;
        bool MTDEtlZposD2 = false;
        int numMTDEtlvalidhits = 0;
        for (const auto hit : track.recHits()) {
          if (hit->isValid() == false)
            continue;
          MTDDetId Hit = hit->geographicalId();
          if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 2)) {
            isETL = true;
            ETLDetId ETLHit = hit->geographicalId();

            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZnegD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZnegD2 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZposD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZposD2 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
          }
        }
        meTrackNumHits_->Fill(-numMTDEtlvalidhits);
        if (isETL && Sigmat0Safe[trackref] < 0.) {
          meTrackNumHitsNT_->Fill(-numMTDEtlvalidhits);
        }

        // --- keeping only tracks with last hit in MTD ---
        if ((track.eta() < -trackMinEtlEta_) && (track.eta() > -trackMaxEtlEta_)) {
          twoETLdiscs = (MTDEtlZnegD1 == true) && (MTDEtlZnegD2 == true);
          if ((MTDEtlZnegD1 == true) || (MTDEtlZnegD2 == true)) {
            meETLTrackEffEtaMtd_[0]->Fill(track.eta());
            meETLTrackEffPhiMtd_[0]->Fill(track.phi());
            meETLTrackEffPtMtd_[0]->Fill(track.pt());
            if (twoETLdiscs) {
              meETLTrackEffEta2Mtd_[0]->Fill(track.eta());
              meETLTrackEffPhi2Mtd_[0]->Fill(track.phi());
              meETLTrackEffPt2Mtd_[0]->Fill(track.pt());
            }
          }
        }
        if ((track.eta() > trackMinEtlEta_) && (track.eta() < trackMaxEtlEta_)) {
          twoETLdiscs = (MTDEtlZposD1 == true) && (MTDEtlZposD2 == true);
          if ((MTDEtlZposD1 == true) || (MTDEtlZposD2 == true)) {
            meETLTrackEffEtaMtd_[1]->Fill(track.eta());
            meETLTrackEffPhiMtd_[1]->Fill(track.phi());
            meETLTrackEffPtMtd_[1]->Fill(track.pt());
            if (twoETLdiscs) {
              meETLTrackEffEta2Mtd_[1]->Fill(track.eta());
              meETLTrackEffPhi2Mtd_[1]->Fill(track.phi());
              meETLTrackEffPt2Mtd_[1]->Fill(track.pt());
            }
          }
        }
      }

      LogDebug("MtdTracksValidation") << "Track p/pt = " << track.p() << " " << track.pt() << " eta " << track.eta()
                                      << " BTL " << isBTL << " ETL " << isETL << " 2disks " << twoETLdiscs;

      // TrackingParticle based matching

      const reco::TrackBaseRef trkrefb(trackref);
      auto tp_info = getMatchedTP(trkrefb);

      meTrackPtTot_->Fill(trackGen.pt());
      meTrackEtaTot_->Fill(std::abs(trackGen.eta()));
      if (tp_info != nullptr && mvaTPSel(**tp_info)) {
        if (track.pt() < 12.) {
          if (isBTL) {
            meBTLTrackMatchedTPPtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                               std::abs(trackGen.pt() - (*tp_info)->pt()));
            meBTLTrackMatchedTPPtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
            meBTLTrackMatchedTPPtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
            meBTLTrackMatchedTPPtResvsPtMtd_->Fill(
                (*tp_info)->pt(), std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
            meBTLTrackMatchedTPDPtvsPtGen_->Fill((*tp_info)->pt(),
                                                 (trackGen.pt() - (*tp_info)->pt()) / (*tp_info)->pt());
            meBTLTrackMatchedTPDPtvsPtMtd_->Fill((*tp_info)->pt(), (track.pt() - (*tp_info)->pt()) / (*tp_info)->pt());
          }
          if (isETL && !twoETLdiscs) {
            meETLTrackMatchedTPPtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                               std::abs(trackGen.pt() - (*tp_info)->pt()));
            meETLTrackMatchedTPPtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
            meETLTrackMatchedTPPtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
            meETLTrackMatchedTPPtResvsPtMtd_->Fill(
                (*tp_info)->pt(), std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
            meETLTrackMatchedTPDPtvsPtGen_->Fill((*tp_info)->pt(),
                                                 (trackGen.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
            meETLTrackMatchedTPDPtvsPtMtd_->Fill((*tp_info)->pt(),
                                                 (track.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
          }
          if (isETL && twoETLdiscs) {
            meETLTrackMatchedTP2PtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                                std::abs(trackGen.pt() - (*tp_info)->pt()));
            meETLTrackMatchedTP2PtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
            meETLTrackMatchedTP2PtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
            meETLTrackMatchedTP2PtResvsPtMtd_->Fill(
                (*tp_info)->pt(), std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
            meETLTrackMatchedTP2DPtvsPtGen_->Fill((*tp_info)->pt(),
                                                  (trackGen.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
            meETLTrackMatchedTP2DPtvsPtMtd_->Fill((*tp_info)->pt(),
                                                  (track.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
          }
        }
        const bool withMTD = (m_tp2detid.find(*tp_info) != m_tp2detid.end());
        LogDebug("MtdTracksValidation") << "Matched with selected TP, MTD sim hits association: " << withMTD;
        if (noCrack) {
          meTrackMatchedTPEffPtTot_->Fill(trackGen.pt());
          if (withMTD) {
            meTrackMatchedTPmtdEffPtTot_->Fill(trackGen.pt());
          }
        }
        meTrackMatchedTPEffEtaTot_->Fill(std::abs(trackGen.eta()));
        if (withMTD) {
          meTrackMatchedTPmtdEffEtaTot_->Fill(std::abs(trackGen.eta()));
        }
        if (isBTL || isETL) {
          if (noCrack) {
            meTrackMatchedTPEffPtMtd_->Fill(trackGen.pt());
            if (isBTL || twoETLdiscs) {
              meTrackMatchedTPEffPtEtl2Mtd_->Fill(trackGen.pt());
            }
            if (withMTD) {
              meTrackMatchedTPmtdEffPtMtd_->Fill(trackGen.pt());
            }
          }
          meTrackMatchedTPEffEtaMtd_->Fill(std::abs(trackGen.eta()));
          if (isBTL || twoETLdiscs) {
            meTrackMatchedTPEffEtaEtl2Mtd_->Fill(std::abs(trackGen.eta()));
          }
          if (withMTD) {
            meTrackMatchedTPmtdEffEtaMtd_->Fill(std::abs(trackGen.eta()));
          }
        }

        size_t nlayers(0);
        float extrho(0.);
        float exteta(0.);
        float extphi(0.);
        float selvar(0.);
        auto accept = checkAcceptance(trackGen, iEvent, iSetup, nlayers, extrho, exteta, extphi, selvar);
        if (accept.first && std::abs(exteta) < trackMaxBtlEta_) {
          meExtraPhiAtBTL_->Fill(angle_units::operators::convertRadToDeg(extphi));
          meExtraBTLeneInCone_->Fill(selvar);
        }
        if (accept.second) {
          if (std::abs(exteta) < trackMaxBtlEta_) {
            meExtraPhiAtBTLmatched_->Fill(angle_units::operators::convertRadToDeg(extphi));
          }
          if (noCrack) {
            meExtraPtMtd_->Fill(trackGen.pt());
            if (nlayers == 2) {
              meExtraPtEtl2Mtd_->Fill(trackGen.pt());
            }
          }
          meExtraEtaMtd_->Fill(std::abs(trackGen.eta()));
          if (nlayers == 2) {
            meExtraEtaEtl2Mtd_->Fill(std::abs(trackGen.eta()));
          }
          if (accept.first && accept.second && !(isBTL || isETL)) {
            edm::LogInfo("MtdTracksValidation")
                << "MtdTracksValidation: extender fail in " << iEvent.id().run() << " " << iEvent.id().event()
                << " pt= " << trackGen.pt() << " eta= " << trackGen.eta();
            meExtraMTDfailExtenderEta_->Fill(std::abs(trackGen.eta()));
            if (noCrack) {
              meExtraMTDfailExtenderPt_->Fill(trackGen.pt());
            }
          }
        }

      }  // TP matching
    }

    if (isGoodVtx) {
      const bool vtxFake = primRecoVtx.isFake();

      if (mvaRecSel(trackGen, primRecoVtx, t0Safe[trackref], Sigmat0Safe[trackref])) {
        // reco-gen matching used for MVA quality flag

        if (noCrack) {
          meMVATrackEffPtTot_->Fill(trackGen.pt());
        }
        meMVATrackEffEtaTot_->Fill(std::abs(trackGen.eta()));

        double dZ = trackGen.vz() - zsim;
        double dT(-9999.);
        double pullT(-9999.);
        if (Sigmat0Safe[trackref] != -1.) {
          dT = t0Safe[trackref] - tsim;
          pullT = dT / Sigmat0Safe[trackref];
        }
        for (const auto& genP : mc->particle_range()) {
          // select status 1 genParticles and match them to the reconstructed track

          float charge = pdTable->particle(HepPDT::ParticleID(genP->pdg_id())) != nullptr
                             ? pdTable->particle(HepPDT::ParticleID(genP->pdg_id()))->charge()
                             : 0.f;
          if (mvaGenSel(*genP, charge)) {
            if (mvaGenRecMatch(*genP, zsim, trackGen, vtxFake)) {
              meMVATrackZposResTot_->Fill(dZ);
              if (noCrack) {
                meMVATrackMatchedEffPtTot_->Fill(trackGen.pt());
              }
              meMVATrackMatchedEffEtaTot_->Fill(std::abs(trackGen.eta()));
              if (isBTL || isETL) {
                meMVATrackResTot_->Fill(dT);
                meMVATrackPullTot_->Fill(pullT);
                if (noCrack) {
                  meMVATrackMatchedEffPtMtd_->Fill(trackGen.pt());
                }
                meMVATrackMatchedEffEtaMtd_->Fill(std::abs(trackGen.eta()));
              }
              break;
            }
          }
        }
      }
    }  // MC truth matich analysis for good PV
  }    //RECO tracks loop
}

const std::pair<bool, bool> MtdTracksValidation::checkAcceptance(const reco::Track& track,
                                                                 const edm::Event& iEvent,
                                                                 edm::EventSetup const& iSetup,
                                                                 size_t& nlayers,
                                                                 float& extrho,
                                                                 float& exteta,
                                                                 float& extphi,
                                                                 float& selvar) {
  bool isMatched(false);
  nlayers = 0;
  extrho = 0.;
  exteta = -999.;
  extphi = -999.;
  selvar = 0.;

  auto geometryHandle = iSetup.getTransientHandle(mtdgeoToken_);
  const MTDGeometry* geom = geometryHandle.product();
  auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
  const MTDTopology* topology = topologyHandle.product();

  auto layerHandle = iSetup.getTransientHandle(mtdlayerToken_);
  const MTDDetLayerGeometry* layerGeo = layerHandle.product();

  auto magfieldHandle = iSetup.getTransientHandle(magfieldToken_);
  const MagneticField* mfield = magfieldHandle.product();

  auto ttrackBuilder = iSetup.getTransientHandle(builderToken_);

  auto tTrack = ttrackBuilder->build(track);
  TrajectoryStateOnSurface tsos = tTrack.outermostMeasurementState();
  float theMaxChi2 = 500.;
  float theNSigma = 10.;
  std::unique_ptr<MeasurementEstimator> theEstimator =
      std::make_unique<Chi2MeasurementEstimator>(theMaxChi2, theNSigma);
  SteppingHelixPropagator prop(mfield, anyDirection);

  auto btlRecHitsHandle = makeValid(iEvent.getHandle(btlRecHitsToken_));
  auto etlRecHitsHandle = makeValid(iEvent.getHandle(etlRecHitsToken_));

  edm::LogVerbatim("MtdTracksValidation")
      << "MtdTracksValidation: extrapolating track, pt= " << track.pt() << " eta= " << track.eta();

  //try BTL
  bool inBTL = false;
  float eneSum(0.);
  const std::vector<const DetLayer*>& layersBTL = layerGeo->allBTLLayers();
  for (const DetLayer* ilay : layersBTL) {
    std::pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos, prop, *theEstimator);
    if (!comp.first)
      continue;
    if (!inBTL) {
      inBTL = true;
      extrho = comp.second.globalPosition().perp();
      exteta = comp.second.globalPosition().eta();
      extphi = comp.second.globalPosition().phi();
      edm::LogVerbatim("MtdTracksValidation") << "MtdTracksValidation: extrapolation at BTL surface, rho= " << extrho
                                              << " eta= " << exteta << " phi= " << extphi;
    }
    std::vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos, prop, *theEstimator);
    for (const auto& detWithState : compDets) {
      const auto& det = detWithState.first;

      // loop on compatible rechits and check energy in a fixed size cone around the extrapolation point

      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: DetId= " << det->geographicalId().rawId()
          << " gp= " << detWithState.second.globalPosition().x() << " " << detWithState.second.globalPosition().y()
          << " " << detWithState.second.globalPosition().z() << " rho= " << detWithState.second.globalPosition().perp()
          << " eta= " << detWithState.second.globalPosition().eta()
          << " phi= " << detWithState.second.globalPosition().phi();

      for (const auto& recHit : *btlRecHitsHandle) {
        BTLDetId detId = recHit.id();
        DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(topology->getMTDTopologyMode()));
        const MTDGeomDet* thedet = geom->idToDet(geoId);
        if (thedet == nullptr)
          throw cms::Exception("MtdTracksValidation") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                                      << detId.rawId() << ") is invalid!" << std::dec << std::endl;
        if (geoId == det->geographicalId()) {
          const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
          const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

          Local3DPoint local_point(0., 0., 0.);
          local_point = topo.pixelToModuleLocalPoint(local_point, detId.row(topo.nrows()), detId.column(topo.nrows()));
          const auto& global_point = thedet->toGlobal(local_point);
          edm::LogVerbatim("MtdTracksValidation")
              << "MtdTracksValidation: Hit id= " << detId.rawId() << " ene= " << recHit.energy()
              << " dr= " << reco::deltaR(global_point, detWithState.second.globalPosition());
          if (reco::deltaR(global_point, detWithState.second.globalPosition()) < cluDRradius_) {
            eneSum += recHit.energy();
            //extrho = detWithState.second.globalPosition().perp();
            //exteta = detWithState.second.globalPosition().eta();
            //extphi = detWithState.second.globalPosition().phi();
          }
        }
      }
    }
    if (eneSum > depositBTLthreshold_) {
      nlayers++;
      selvar = eneSum;
      isMatched = true;
      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: BTL matched, energy= " << eneSum << " #layers= " << nlayers;
    }
  }
  if (inBTL) {
    return std::make_pair(inBTL, isMatched);
  }

  //try ETL
  bool inETL = false;
  const std::vector<const DetLayer*>& layersETL = layerGeo->allETLLayers();
  for (const DetLayer* ilay : layersETL) {
    size_t hcount(0);
    const BoundDisk& disk = static_cast<const MTDSectorForwardDoubleLayer*>(ilay)->specificSurface();
    const double diskZ = disk.position().z();
    if (tsos.globalPosition().z() * diskZ < 0)
      continue;  // only propagate to the disk that's on the same side
    std::pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos, prop, *theEstimator);
    if (!comp.first)
      continue;
    if (!inETL) {
      inETL = true;
      extrho = comp.second.globalPosition().perp();
      exteta = comp.second.globalPosition().eta();
      extphi = comp.second.globalPosition().phi();
    }
    edm::LogVerbatim("MtdTracksValidation") << "MtdTracksValidation: extrapolation at ETL surface, rho= " << extrho
                                            << " eta= " << exteta << " phi= " << extphi;
    std::vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos, prop, *theEstimator);
    for (const auto& detWithState : compDets) {
      const auto& det = detWithState.first;

      // loop on compatible rechits and check hits in a fixed size cone around the extrapolation point

      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: DetId= " << det->geographicalId().rawId()
          << " gp= " << detWithState.second.globalPosition().x() << " " << detWithState.second.globalPosition().y()
          << " " << detWithState.second.globalPosition().z() << " rho= " << detWithState.second.globalPosition().perp()
          << " eta= " << detWithState.second.globalPosition().eta()
          << " phi= " << detWithState.second.globalPosition().phi();

      for (const auto& recHit : *etlRecHitsHandle) {
        ETLDetId detId = recHit.id();
        DetId geoId = detId.geographicalId();
        const MTDGeomDet* thedet = geom->idToDet(geoId);
        if (thedet == nullptr)
          throw cms::Exception("MtdTracksValidation") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                                      << detId.rawId() << ") is invalid!" << std::dec << std::endl;
        if (geoId == det->geographicalId()) {
          const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
          const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

          Local3DPoint local_point(topo.localX(recHit.row()), topo.localY(recHit.column()), 0.);
          const auto& global_point = thedet->toGlobal(local_point);
          edm::LogVerbatim("MtdTracksValidation")
              << "MtdTracksValidation: Hit id= " << detId.rawId() << " time= " << recHit.time()
              << " dr= " << reco::deltaR(global_point, detWithState.second.globalPosition());
          if (reco::deltaR(global_point, detWithState.second.globalPosition()) < cluDRradius_) {
            hcount++;
            if (hcount == 1) {
              //extrho = detWithState.second.globalPosition().perp();
              //exteta = detWithState.second.globalPosition().eta();
              //extphi = detWithState.second.globalPosition().phi();
            }
          }
        }
      }
    }
    if (hcount > 0) {
      nlayers++;
      selvar = (float)hcount;
      isMatched = true;
      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: ETL matched, counts= " << hcount << " #layers= " << nlayers;
    }
  }

  if (!inBTL && !inETL) {
    edm::LogVerbatim("MtdTracksValidation")
        << "MtdTracksValidation: track not extrapolating to MTD: pt= " << track.pt() << " eta= " << track.eta()
        << " phi= " << track.phi() << " vz= " << track.vz()
        << " vxy= " << std::sqrt(track.vx() * track.vx() + track.vy() * track.vy());
  }
  return std::make_pair(inETL, isMatched);
}

// ------------ method for histogram booking ------------
void MtdTracksValidation::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking
  meBTLTrackRPTime_ = ibook.book1D("TrackBTLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meBTLTrackEffEtaTot_ = ibook.book1D("TrackBTLEffEtaTot", "Track efficiency vs eta (Tot);#eta_{RECO}", 100, -1.6, 1.6);
  meBTLTrackEffPhiTot_ =
      ibook.book1D("TrackBTLEffPhiTot", "Track efficiency vs phi (Tot);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackEffPtTot_ = ibook.book1D("TrackBTLEffPtTot", "Track efficiency vs pt (Tot);pt_{RECO} [GeV]", 50, 0, 10);
  meBTLTrackEffEtaMtd_ = ibook.book1D("TrackBTLEffEtaMtd", "Track efficiency vs eta (Mtd);#eta_{RECO}", 100, -1.6, 1.6);
  meBTLTrackEffPhiMtd_ =
      ibook.book1D("TrackBTLEffPhiMtd", "Track efficiency vs phi (Mtd);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackEffPtMtd_ = ibook.book1D("TrackBTLEffPtMtd", "Track efficiency vs pt (Mtd);pt_{RECO} [GeV]", 50, 0, 10);
  meBTLTrackPtRes_ =
      ibook.book1D("TrackBTLPtRes", "Track pT resolution  ;pT_{Gentrack}-pT_{MTDtrack}/pT_{Gentrack} ", 100, -0.1, 0.1);
  meETLTrackRPTime_ = ibook.book1D("TrackETLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meETLTrackEffEtaTot_[0] =
      ibook.book1D("TrackETLEffEtaTotZneg", "Track efficiency vs eta (Tot) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEtaTot_[1] =
      ibook.book1D("TrackETLEffEtaTotZpos", "Track efficiency vs eta (Tot) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhiTot_[0] =
      ibook.book1D("TrackETLEffPhiTotZneg", "Track efficiency vs phi (Tot) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhiTot_[1] =
      ibook.book1D("TrackETLEffPhiTotZpos", "Track efficiency vs phi (Tot) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPtTot_[0] =
      ibook.book1D("TrackETLEffPtTotZneg", "Track efficiency vs pt (Tot) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPtTot_[1] =
      ibook.book1D("TrackETLEffPtTotZpos", "Track efficiency vs pt (Tot) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffEtaMtd_[0] =
      ibook.book1D("TrackETLEffEtaMtdZneg", "Track efficiency vs eta (Mtd) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEtaMtd_[1] =
      ibook.book1D("TrackETLEffEtaMtdZpos", "Track efficiency vs eta (Mtd) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhiMtd_[0] =
      ibook.book1D("TrackETLEffPhiMtdZneg", "Track efficiency vs phi (Mtd) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhiMtd_[1] =
      ibook.book1D("TrackETLEffPhiMtdZpos", "Track efficiency vs phi (Mtd) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPtMtd_[0] =
      ibook.book1D("TrackETLEffPtMtdZneg", "Track efficiency vs pt (Mtd) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPtMtd_[1] =
      ibook.book1D("TrackETLEffPtMtdZpos", "Track efficiency vs pt (Mtd) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffEta2Mtd_[0] =
      ibook.book1D("TrackETLEffEta2MtdZneg", "Track efficiency vs eta (Mtd 2 hit) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEta2Mtd_[1] =
      ibook.book1D("TrackETLEffEta2MtdZpos", "Track efficiency vs eta (Mtd 2 hit) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhi2Mtd_[0] = ibook.book1D(
      "TrackETLEffPhi2MtdZneg", "Track efficiency vs phi (Mtd 2 hit) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhi2Mtd_[1] = ibook.book1D(
      "TrackETLEffPhi2MtdZpos", "Track efficiency vs phi (Mtd 2 hit) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPt2Mtd_[0] =
      ibook.book1D("TrackETLEffPt2MtdZneg", "Track efficiency vs pt (Mtd 2 hit) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPt2Mtd_[1] =
      ibook.book1D("TrackETLEffPt2MtdZpos", "Track efficiency vs pt (Mtd 2 hit) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackPtRes_ =
      ibook.book1D("TrackETLPtRes", "Track pT resolution;pT_{Gentrack}-pT_{MTDtrack}/pT_{Gentrack} ", 100, -0.1, 0.1);

  meTracktmtd_ = ibook.book1D("Tracktmtd", "Track time from TrackExtenderWithMTD;tmtd [ns]", 150, 1, 16);
  meTrackt0Src_ = ibook.book1D("Trackt0Src", "Track time from TrackExtenderWithMTD;t0Src [ns]", 100, -1.5, 1.5);
  meTrackSigmat0Src_ =
      ibook.book1D("TrackSigmat0Src", "Time Error from TrackExtenderWithMTD; #sigma_{t0Src} [ns]", 100, 0, 0.1);

  meTrackt0Pid_ = ibook.book1D("Trackt0Pid", "Track t0 as stored in TofPid;t0 [ns]", 100, -1, 1);
  meTrackSigmat0Pid_ = ibook.book1D("TrackSigmat0Pid", "Sigmat0 as stored in TofPid; #sigma_{t0} [ns]", 100, 0, 0.1);
  meTrackt0SafePid_ = ibook.book1D("Trackt0SafePID", "Track t0 Safe as stored in TofPid;t0 [ns]", 100, -1, 1);
  meTrackSigmat0SafePid_ =
      ibook.book1D("TrackSigmat0SafePID", "Sigmat0 Safe as stored in TofPid; #sigma_{t0} [ns]", 100, 0, 0.1);
  
  meTrackDtSign_[0] = ibook.book1D("TrackDtSign_Pion", "DeltaT(trk, vtx) significance under pion hypothesis;sign(t_0) [ns]", 40, 0, 20);
  meTrackDtSign_[1] = ibook.book1D("TrackDtSign_Kaon", "DeltaT(trk, vtx) significance under kaon hypothesis;sign(t_0) [ns]", 40, 0, 20);
  meTrackDtSign_[2] = ibook.book1D("TrackDtSign_Proton", "DeltaT(trk, vtx) significance under proton hypothesis;sign(t_0) [ns]", 40, 0, 20);

  meTrackProb_[0] = ibook.book1D("TrackProb_Pion", "Probability of pion hp from vertex time compatibility;P", 100, 0, 1);
  meTrackProb_[1] = ibook.book1D("TrackProb_Kaon", "Probability of kaon hp from vertex time compatibility;P", 100, 0, 1);
  meTrackProb_[2] = ibook.book1D("TrackProb_Proton", "Probability of proton hp from vertex time compatibility;P", 100, 0, 1);

  meTrackTimeChi2_PionTrack_[0] = ibook.book1D("TrackTimeChi2_PionTrack_PionHp", "Time Chi2 under pion hypothesis for pion-ID track;#chi^{2}_t", 100, 0, 100);
  meTrackTimeChi2_PionTrack_[1] = ibook.book1D("TrackTimeChi2_PionTrack_KaonHp", "Time Chi2 under kaon hypothesis for pion-ID track;#chi^{2}_t", 100, 0, 100);
  meTrackTimeChi2_PionTrack_[2] = ibook.book1D("TrackTimeChi2_PionTrack_ProtonHp", "Time Chi2 under proton hypothesis for pion-ID track;#chi^{2}_t", 100, 0, 100);
  meTrackTimeChi2_KaonTrack_[0] = ibook.book1D("TrackTimeChi2_KaonTrack_PionHp", "Time Chi2 under pion hypothesis for kaon-ID track;#chi^{2}_t", 100, 0, 100);
  meTrackTimeChi2_KaonTrack_[1] = ibook.book1D("TrackTimeChi2_KaonTrack_KaonHp", "Time Chi2 under kaon hypothesis for kaon-ID track;#chi^{2}_t", 100, 0, 100);
  meTrackTimeChi2_KaonTrack_[2] = ibook.book1D("TrackTimeChi2_KaonTrack_ProtonHp", "Time Chi2 under proton hypothesis for kaon-ID track;#chi^{2}_t", 100, 0, 100);
  meTrackTimeChi2_ProtonTrack_[0] = ibook.book1D("TrackTimeChi2_ProtonTrack_PionHp", "Time Chi2 under pion hypothesis for proton-ID track;#chi^{2}_t", 100, 0, 100);
  meTrackTimeChi2_ProtonTrack_[1] = ibook.book1D("TrackTimeChi2_ProtonTrack_KaonHp", "Time Chi2 under kaon hypothesis for proton-ID track;#chi^{2}_t", 100, 0, 100);
  meTrackTimeChi2_ProtonTrack_[2] = ibook.book1D("TrackTimeChi2_ProtonTrack_ProtonHp", "Time Chi2 under proton hypothesis for proton-ID track;#chi^{2}_t", 100, 0, 100);

  meDeltaChisq_truePion_vsKaon_   = ibook.book1D("DeltaChisq_truePion_vsKaon", "Delta Chi2 between pion and kaon hypothesis for pion track;#chi^{2}_{#pi} - #chi^{2}_{K}", 80, -20, 20);
  meDeltaChisq_truePion_vsProton_ = ibook.book1D("DeltaChisq_truePion_vsProton", "Delta Chi2 between pion and proton hypothesis for pion track;#chi^{2}_{#pi} - #chi^{2}_{p}", 80, -20, 20);

  meDeltaChisq_trueKaon_vsPion_   = ibook.book1D("DeltaChisq_trueKaon_vsPion", "Delta Chi2 between kaon and pion hypothesis for kaon track;#chi^{2}_{K} - #chi^{2}_{#pi}", 80, -20, 20);
  meDeltaChisq_trueKaon_vsProton_ = ibook.book1D("DeltaChisq_trueKaon_vsProton", "Delta Chi2 between kaon and proton hypothesis for kaon track;#chi^{2}_{K} - #chi^{2}_{p}", 80, -20, 20);

  meDeltaChisq_trueProton_vsPion_ = ibook.book1D("DeltaChisq_trueProton_vsPion", "Delta Chi2 between proton and pion hypothesis for proton track;#chi^{2}_{p} - #chi^{2}_{#pi}", 80, -20, 20);
  meDeltaChisq_trueProton_vsKaon_ = ibook.book1D("DeltaChisq_trueProton_vsKaon", "Delta Chi2 between proton and kaon hypothesis for proton track;#chi^{2}_{p} - #chi^{2}_{K}", 80, -20, 20);

  meTrackBtlMatchChi2_[0] = ibook.book1D("TrackBtlMatchChi2_Pion", "Time chi2 with BTL hit match under pion hp.;#chi^{2}_t", 30, 0 , 15);
  meTrackBtlMatchChi2_[1] = ibook.book1D("TrackBtlMatchChi2_Kaon", "Time chi2 with BTL hit match under kaon hp.;#chi^{2}_t", 30, 0 , 15);
  meTrackBtlMatchChi2_[2] = ibook.book1D("TrackBtlMatchChi2_Proton", "Time chi2 with BTL hit match under proton hp.;#chi^{2}_t", 30, 0 , 15);
  meTrackEtlMatchChi2_[0] = ibook.book1D("TrackEtlMatchChi2_Pion", "Time chi2 with ETL hit match under pion hp.;#chi^{2}_t", 30, 0 , 15);  
  meTrackEtlMatchChi2_[1] = ibook.book1D("TrackEtlMatchChi2_Kaon", "Time chi2 with ETL hit match under kaon hp.;#chi^{2}_t", 30, 0 , 15);
  meTrackEtlMatchChi2_[2] = ibook.book1D("TrackEtlMatchChi2_Proton", "Time chi2 with ETL hit match under proton hp.;#chi^{2}_t", 30, 0 , 15);

  meRatioMatchedTracksVsP_BTL_num_ = ibook.book1D("RatioMatchedTracksVsP_BTL_num", "Number of tracks matched to MTD BTL hit vs p; p [GeV]", 10, 0., 10.);
  meRatioMatchedTracksVsP_BTL_den_ = ibook.book1D("RatioMatchedTracksVsP_BTL_den", "Number of all BTL tracks vs p; p [GeV]", 10, 0., 10.);
  meRatioMatchedTracksVsP_ETL_num_ = ibook.book1D("RatioMatchedTracksVsP_ETL_num", "Number of tracks matched to MTD ETL hit vs p; p [GeV]", 10, 0., 10.);
  meRatioMatchedTracksVsP_ETL_den_ = ibook.book1D("RatioMatchedTracksVsP_ETL_den", "Number of all ETL tracks vs p; p [GeV]", 10, 0., 10.);

  meTruePi_ = ibook.book1D("pTruePi", "True pi momentum spectrum;p [GeV]", 25, 0., 10.);
  meTrueK_ = ibook.book1D("pTrueK", "True K momentum spectrum;p [GeV]", 25, 0., 10.);
  meTrueP_ = ibook.book1D("pTrueP", "True P momentum spectrum;p [GeV]", 25, 0., 10.);

  mePIDPi_ = ibook.book1D("pPIDPi", "PID-pi momentum spectrum;p [GeV]", 25, 0., 10.);
  mePIDK_ = ibook.book1D("pPIDK", "PID-K momentum spectrum;p [GeV]", 25, 0., 10.);
  mePIDP_ = ibook.book1D("pPIDP", "PID-P momentum spectrum;p [GeV]", 25, 0., 10.);

  meTruePiNoPID_ =
      ibook.book1D("pTruePiNoPID", "True pi NoPID momentum spectrum;p [GeV]", 25, 0., 10.);
  meTruePiAsPi_ =
      ibook.book1D("pTruePiAsPi", "True pi as pi momentum spectrum;p [GeV]", 25, 0., 10.);
  meTruePiAsK_ =
      ibook.book1D("pTruePiAsK", "True pi as k momentum spectrum;p [GeV]", 25, 0., 10.);
  meTruePiAsP_ =
      ibook.book1D("pTruePiAsP", "True pi as p momentum spectrum;p [GeV]", 25, 0., 10.);

  meTrueKNoPID_ =
      ibook.book1D("pTrueKNoPID", "True k NoPID momentum spectrum;p [GeV]", 25, 0., 10.);
  meTrueKAsPi_ =
      ibook.book1D("pTrueKAsPi", "True k as pi momentum spectrum;p [GeV]", 25, 0., 10.);
  meTrueKAsK_ =
      ibook.book1D("pTrueKAsK", "True k as k momentum spectrum;p [GeV]", 25, 0., 10.);
  meTrueKAsP_ =
      ibook.book1D("pTrueKAsP", "True k as p momentum spectrum;p [GeV]", 25, 0., 10.);

  meTruePNoPID_ =
      ibook.book1D("pTruePNoPID", "True p NoPID momentum spectrum;p [GeV]", 25, 0., 10.);
  meTruePAsPi_ =
      ibook.book1D("pTruePAsPi", "True p as pi momentum spectrum;p [GeV]", 25, 0., 10.);
  meTruePAsK_ =
      ibook.book1D("pTruePAsK", "True p as k momentum spectrum;p [GeV]", 25, 0., 10.);
  meTruePAsP_ =
      ibook.book1D("pTruePAsP", "True p as p momentum spectrum;p [GeV]", 25, 0., 10.);

  // HSCP studies
  
  meTOF_ = ibook.book1D("TOF", "Time of flight under reco PID;TOF [ns]", 100, 0, 20);
  meTmtd_ = ibook.book1D("Tmtd", "Time at MTD;t(MTD) [ns]", 100, 0, 20);
  meInvBeta_[0] = ibook.book1D("InvBeta_Pion", "1/#beta for PID-pion tracks;1/#beta", 100, 0.9, 1.1); //xmax = 1.7 in TDR
  meInvBeta_[1] = ibook.book1D("InvBeta_Kaon", "1/#beta for PID-kaon tracks;1/#beta", 100, 0.9, 1.1); //xmax = 1.7 in TDR
  meInvBeta_[2] = ibook.book1D("InvBeta_Proton", "1/#beta for PID-proton tracks;1/#beta", 100, 0.9, 1.1); //xmax = 1.7 in TDR
  meInvBetaVsP_ = ibook.book2D("BetaVsP", "1/#beta vs p; p [GeV]; 1/#beta", 25, 0, 5, 100, 0.9, 1.1);

  meTrackNumHits_ = ibook.book1D("TrackNumHits", "Number of valid MTD hits per track ; Number of hits", 10, -5, 5);
  meTrackNumHitsNT_ = ibook.book1D(
      "TrackNumHitsNT", "Number of valid MTD hits per track no time associated; Number of hits", 10, -5, 5);
  meTrackMVAQual_ = ibook.book1D("TrackMVAQual", "Track MVA Quality as stored in Value Map ; MVAQual", 100, 0, 1);
  meTrackPathLenghtvsEta_ = ibook.bookProfile(
      "TrackPathLenghtvsEta", "MTD Track pathlength vs MTD track Eta;|#eta|;Pathlength", 100, 0, 3.2, 100.0, 400.0, "S");

  meMVATrackEffPtTot_ = ibook.book1D("MVAEffPtTot", "Pt of tracks associated to LV; track pt [GeV] ", 110, 0., 11.);
  meMVATrackMatchedEffPtTot_ =
      ibook.book1D("MVAMatchedEffPtTot", "Pt of tracks associated to LV matched to GEN; track pt [GeV] ", 110, 0., 11.);
  meMVATrackMatchedEffPtMtd_ = ibook.book1D(
      "MVAMatchedEffPtMtd", "Pt of tracks associated to LV matched to GEN with time; track pt [GeV] ", 110, 0., 11.);

  meExtraPtMtd_ = ibook.book1D("ExtraPtMtd", "Pt of tracks extrapolated to hits; track pt [GeV] ", 110, 0., 11.);
  meExtraPtEtl2Mtd_ =
      ibook.book1D("ExtraPtEtl2Mtd", "Pt of tracks extrapolated to hits, 2 ETL layers; track pt [GeV] ", 110, 0., 11.);

  meTrackPtTot_ = ibook.book1D("TrackPtTot", "Pt of tracks ; track pt [GeV] ", 110, 0., 11.);
  meTrackMatchedTPEffPtTot_ =
      ibook.book1D("MatchedTPEffPtTot", "Pt of tracks  matched to TP; track pt [GeV] ", 110, 0., 11.);
  meTrackMatchedTPEffPtMtd_ =
      ibook.book1D("MatchedTPEffPtMtd", "Pt of tracks  matched to TP with time; track pt [GeV] ", 110, 0., 11.);
  meTrackMatchedTPEffPtEtl2Mtd_ = ibook.book1D(
      "MatchedTPEffPtEtl2Mtd", "Pt of tracks  matched to TP with time, 2 ETL hits; track pt [GeV] ", 110, 0., 11.);

  meBTLTrackMatchedTPPtResMtd_ = ibook.book1D(
      "TrackMatchedTPBTLPtResMtd",
      "Pt resolution of tracks matched to TP-BTL hit  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
      100,
      0.,
      4.);
  meETLTrackMatchedTPPtResMtd_ = ibook.book1D(
      "TrackMatchedTPETLPtResMtd",
      "Pt resolution of tracks matched to TP-ETL hit  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
      100,
      0.,
      4.);
  meETLTrackMatchedTP2PtResMtd_ = ibook.book1D(
      "TrackMatchedTPETL2PtResMtd",
      "Pt resolution of tracks matched to TP-ETL 2hits  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
      100,
      0.,
      4.);
  meBTLTrackMatchedTPPtRatioGen_ = ibook.book1D(
      "TrackMatchedTPBTLPtRatioGen", "Pt ratio of Gentracks (BTL)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
  meETLTrackMatchedTPPtRatioGen_ = ibook.book1D(
      "TrackMatchedTPETLPtRatioGen", "Pt ratio of Gentracks (ETL 1hit)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
  meETLTrackMatchedTP2PtRatioGen_ = ibook.book1D(
      "TrackMatchedTPETL2PtRatioGen", "Pt ratio of Gentracks (ETL 2hits)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
  meBTLTrackMatchedTPPtRatioMtd_ = ibook.book1D("TrackMatchedTPBTLPtRatioMtd",
                                                "Pt ratio of tracks matched to TP-BTL hits  ;pT_{MTDtrack}/pT_{truth} ",
                                                100,
                                                0.9,
                                                1.1);
  meETLTrackMatchedTPPtRatioMtd_ = ibook.book1D("TrackMatchedTPETLPtRatioMtd",
                                                "Pt ratio of tracks matched to TP-ETL hits  ;pT_{MTDtrack}/pT_{truth} ",
                                                100,
                                                0.9,
                                                1.1);
  meETLTrackMatchedTP2PtRatioMtd_ =
      ibook.book1D("TrackMatchedTPETL2PtRatioMtd",
                   "Pt ratio of tracks matched to TP-ETL 2hits  ;pT_{MTDtrack}/pT_{truth} ",
                   100,
                   0.9,
                   1.1);
  meBTLTrackMatchedTPPtResvsPtMtd_ = ibook.bookProfile("TrackMatchedTPBTLPtResvsPtMtd",
                                                       "Pt resolution of tracks matched to TP-BTL hit vs Pt;pT_{truth} "
                                                       "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                                                       20,
                                                       0.7,
                                                       10.,
                                                       0.,
                                                       4.,
                                                       "s");
  meETLTrackMatchedTPPtResvsPtMtd_ = ibook.bookProfile("TrackMatchedTPETLPtResvsPtMtd",
                                                       "Pt resolution of tracks matched to TP-ETL hit vs Pt;pT_{truth} "
                                                       "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                                                       20,
                                                       0.7,
                                                       10.,
                                                       0.,
                                                       4.,
                                                       "s");
  meETLTrackMatchedTP2PtResvsPtMtd_ =
      ibook.bookProfile("TrackMatchedTPETL2PtResvsPtMtd",
                        "Pt resolution of tracks matched to TP-ETL 2hits Pt pT;pT_{truth} "
                        "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                        20,
                        0.7,
                        10.,
                        0.,
                        4.,
                        "s");
  meBTLTrackMatchedTPDPtvsPtGen_ = ibook.bookProfile(
      "TrackMatchedTPBTLDPtvsPtGen",
      "Pt relative difference of Gentracks (BTL) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
      20,
      0.7,
      10.,
      -0.1,
      0.1,
      "s");
  meETLTrackMatchedTPDPtvsPtGen_ = ibook.bookProfile(
      "TrackMatchedTPETLDPtvsPtGen",
      "Pt relative difference of Gentracks (ETL 1hit) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
      20,
      0.7,
      10.,
      -0.1,
      0.1,
      "s");
  meETLTrackMatchedTP2DPtvsPtGen_ = ibook.bookProfile(
      "TrackMatchedTPETL2DPtvsPtGen",
      "Pt relative difference  of Gentracks (ETL 2hits) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
      20,
      0.7,
      10.,
      -0.1,
      0.1,
      "s");
  meBTLTrackMatchedTPDPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPBTLDPtvsPtMtd",
                                                     "Pt relative difference of tracks matched to TP-BTL hits vs "
                                                     "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                     20,
                                                     0.7,
                                                     10.,
                                                     -0.1,
                                                     0.1,
                                                     "s");
  meETLTrackMatchedTPDPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPETLDPtvsPtMtd",
                                                     "Pt relative difference of tracks matched to TP-ETL hits vs "
                                                     "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                     20,
                                                     0.7,
                                                     10.,
                                                     -0.1,
                                                     0.1,
                                                     "s");
  meETLTrackMatchedTP2DPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPETL2DPtvsPtMtd",
                                                      "Pt relative difference of tracks matched to TP-ETL 2hits vs "
                                                      "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                      20,
                                                      0.7,
                                                      10.,
                                                      -0.1,
                                                      0.1,
                                                      "s");

  meTrackMatchedTPmtdEffPtTot_ =
      ibook.book1D("MatchedTPmtdEffPtTot", "Pt of tracks  matched to TP-mtd hit; track pt [GeV] ", 110, 0., 11.);
  meTrackMatchedTPmtdEffPtMtd_ = ibook.book1D(
      "MatchedTPmtdEffPtMtd", "Pt of tracks  matched to TP-mtd hit with time; track pt [GeV] ", 110, 0., 11.);

  meMVATrackEffEtaTot_ = ibook.book1D("MVAEffEtaTot", "Eta of tracks associated to LV; track eta ", 66, 0., 3.3);
  meMVATrackMatchedEffEtaTot_ =
      ibook.book1D("MVAMatchedEffEtaTot", "Eta of tracks associated to LV matched to GEN; track eta ", 66, 0., 3.3);
  meMVATrackMatchedEffEtaMtd_ = ibook.book1D(
      "MVAMatchedEffEtaMtd", "Eta of tracks associated to LV matched to GEN with time; track eta ", 66, 0., 3.3);

  meExtraEtaMtd_ = ibook.book1D("ExtraEtaMtd", "Eta of tracks extrapolated to hits; track eta ", 66, 0., 3.3);
  meExtraEtaEtl2Mtd_ =
      ibook.book1D("ExtraEtaEtl2Mtd", "Eta of tracks extrapolated to hits, 2 ETL layers; track eta ", 66, 0., 3.3);

  meTrackEtaTot_ = ibook.book1D("TrackEtaTot", "Eta of tracks ; track eta ", 66, 0., 3.3);
  meTrackMatchedTPEffEtaTot_ =
      ibook.book1D("MatchedTPEffEtaTot", "Eta of tracks  matched to TP; track eta ", 66, 0., 3.3);
  meMVATrackEffEtaTot_ = ibook.book1D("MVAEffEtaTot", "Eta of tracks ; track eta ", 66, 0., 3.3);
  meTrackMatchedTPEffEtaMtd_ =
      ibook.book1D("MatchedTPEffEtaMtd", "Eta of tracks  matched to TP with time; track eta ", 66, 0., 3.3);
  meTrackMatchedTPEffEtaEtl2Mtd_ = ibook.book1D(
      "MatchedTPEffEtaEtl2Mtd", "Eta of tracks  matched to TP with time, 2 ETL hits; track eta ", 66, 0., 3.3);

  meTrackMatchedTPmtdEffEtaTot_ =
      ibook.book1D("MatchedTPmtdEffEtaTot", "Eta of tracks  matched to TP-mtd hit; track eta ", 66, 0., 3.3);
  meTrackMatchedTPmtdEffEtaMtd_ =
      ibook.book1D("MatchedTPmtdEffEtaMtd", "Eta of tracks  matched to TP-mtd hit with time; track eta ", 66, 0., 3.3);

  meMVATrackResTot_ = ibook.book1D(
      "MVATrackRes", "t_{rec} - t_{sim} for LV associated tracks; t_{rec} - t_{sim} [ns] ", 120, -0.15, 0.15);
  meMVATrackPullTot_ =
      ibook.book1D("MVATrackPull", "Pull for associated tracks; (t_{rec}-t_{sim})/#sigma_{t}", 50, -5., 5.);
  meMVATrackZposResTot_ = ibook.book1D(
      "MVATrackZposResTot", "Z_{PCA} - Z_{sim} for associated tracks;Z_{PCA} - Z_{sim} [cm] ", 100, -0.1, 0.1);

  meExtraPhiAtBTL_ =
      ibook.book1D("ExtraPhiAtBTL", "Phi at BTL surface of extrapolated tracks; phi [deg]", 720, -180., 180.);
  meExtraPhiAtBTLmatched_ = ibook.book1D("ExtraPhiAtBTLmatched",
                                         "Phi at BTL surface of extrapolated tracksi matched with BTL hits; phi [deg]",
                                         720,
                                         -180.,
                                         180.);
  meExtraBTLeneInCone_ = ibook.book1D(
      "ExtraBTLeneInCone", "BTL reconstructed energy in cone arounnd extrapolated track; E [MeV]", 100, 0., 50.);
  meExtraMTDfailExtenderEta_ =
      ibook.book1D("ExtraMTDfailExtenderEta",
                   "Eta of tracks extrapolated to MTD with no track extender match to hits; track eta",
                   66,
                   0.,
                   3.3);
  ;
  meExtraMTDfailExtenderPt_ =
      ibook.book1D("ExtraMTDfailExtenderPt",
                   "Pt of tracks extrapolated to MTD with no track extender match to hits; track pt [GeV] ",
                   110,
                   0.,
                   11.);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdTracksValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("btlSimHits", edm::InputTag("mix", "g4SimHitsFastTimerHitsBarrel"));
  desc.add<edm::InputTag>("etlSimHits", edm::InputTag("mix", "g4SimHitsFastTimerHitsEndcap"));
  desc.add<edm::InputTag>("btlRecHits", edm::InputTag("mtdRecHits", "FTLBarrel"));
  desc.add<edm::InputTag>("etlRecHits", edm::InputTag("mtdRecHits", "FTLEndcap"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));

  desc.add<edm::InputTag>("tofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
  desc.add<edm::InputTag>("tofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
  desc.add<edm::InputTag>("tofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
  desc.add<edm::InputTag>("beta", edm::InputTag("trackExtenderWithMTD:generalTrackBeta"));
  desc.add<edm::InputTag>("dtSignPi", edm::InputTag("tofPID:dtSignPi"));
  desc.add<edm::InputTag>("dtSignK", edm::InputTag("tofPID:dtSignK"));
  desc.add<edm::InputTag>("dtSignP", edm::InputTag("tofPID:dtSignP"));
  desc.add<edm::InputTag>("probPi", edm::InputTag("tofPID:probPi"));
  desc.add<edm::InputTag>("probK", edm::InputTag("tofPID:probK"));
  desc.add<edm::InputTag>("probP", edm::InputTag("tofPID:probP"));
  desc.add<edm::InputTag>("timeChisqPi", edm::InputTag("tofPID:timeChisqPi"));
  desc.add<edm::InputTag>("timeChisqK", edm::InputTag("tofPID:timeChisqK"));
  desc.add<edm::InputTag>("timeChisqP", edm::InputTag("tofPID:timeChisqP"));
  desc.add<edm::InputTag>("btlMatchTimeChi2", edm::InputTag("trackExtenderWithMTD:btlMatchTimeChi2"));
  desc.add<edm::InputTag>("etlMatchTimeChi2", edm::InputTag("trackExtenderWithMTD:etlMatchTimeChi2"));

  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  // desc.add<double>("trackMinimumPt", 0.7);  // [GeV]
  desc.add<double>("trackMinimumPt", 0.);  // [GeV]
  desc.add<double>("trackMaximumBtlEta", 1.5);
  desc.add<double>("trackMinimumEtlEta", 1.6);
  desc.add<double>("trackMaximumEtlEta", 3.);
  desc.addUntracked<bool>("optionalPlots", true);

  descriptions.add("mtdTracksValid", desc);
}

const bool MtdTracksValidation::mvaGenSel(const HepMC::GenParticle& gp, const float& charge) {
  bool match = false;
  if (gp.status() != 1) {
    return match;
  }
  match = charge != 0.f && gp.momentum().perp() > pTcut_ && std::abs(gp.momentum().eta()) < etacutGEN_;
  return match;
}

const bool MtdTracksValidation::mvaTPSel(const TrackingParticle& tp) {
  bool match = false;
  if (tp.status() != 1) {
    return match;
  }
  auto x_pv = tp.parentVertex()->position().x();
  auto y_pv = tp.parentVertex()->position().y();
  auto z_pv = tp.parentVertex()->position().z();

  auto r_pv = std::sqrt(x_pv * x_pv + y_pv * y_pv);

  match = tp.charge() != 0 && tp.pt() > pTcut_ && std::abs(tp.eta()) < etacutGEN_ && r_pv < rBTL_ && z_pv < zETL_;
  return match;
}

const bool MtdTracksValidation::mvaRecSel(const reco::TrackBase& trk,
                                          const reco::Vertex& vtx,
                                          const double& t0,
                                          const double& st0) {
  bool match = false;
  match = trk.pt() > pTcut_ && std::abs(trk.eta()) < etacutREC_ &&
          (std::abs(trk.vz() - vtx.z()) <= deltaZcut_ || vtx.isFake());
  if (st0 > 0.) {
    match = match && std::abs(t0 - vtx.t()) < 3. * st0;
  }
  return match;
}

const bool MtdTracksValidation::mvaGenRecMatch(const HepMC::GenParticle& genP,
                                               const double& zsim,
                                               const reco::TrackBase& trk,
                                               const bool& vtxFake) {
  bool match = false;
  double dR = reco::deltaR(genP.momentum(), trk.momentum());
  double genPT = genP.momentum().perp();
  match = std::abs(genPT - trk.pt()) < trk.pt() * deltaPTcut_ && dR < deltaDRcut_ &&
          (std::abs(trk.vz() - zsim) < deltaZcut_ || vtxFake);
  return match;
}

const edm::Ref<std::vector<TrackingParticle>>* MtdTracksValidation::getMatchedTP(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return nullptr;

  //matched TP equal to any TP associated to in time events
  for (const auto& tp : found->val) {
    if (tp.first->eventId().bunchCrossing() == 0)
      return &tp.first;
  }

  // reco track not matched to any TP from vertex
  return nullptr;
}

DEFINE_FWK_MODULE(MtdTracksValidation);
