#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

// sim, reco hits and clusters
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"

#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToRecoClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociatorBaseImpl.h"

// reco track and vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// TrackingParticle
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

// Fastjet
#include <fastjet/internal/base.hh>
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

// HepPDTRecord
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// pile-up
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// associator
#include "SimTracker/VertexAssociation/interface/calculateVertexSharedTracks.h"

// vertexing
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"

// simulated vertex
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"

// DQM
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Ntuplizer utilities
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TFile.h"

// set namespaces
using namespace edm;



// class declaration
class Primary4DVertexValidation : public DQMEDAnalyzer {
  typedef math::XYZTLorentzVector LorentzVector;

  // auxiliary class holding simulated vertices
  struct simPrimaryVertex {
    simPrimaryVertex(double x1, double y1, double z1, double t1)
        : x(x1),
          y(y1),
          z(z1),
          t(t1),
          pt(0),
          ptsq(0),
          closest_vertex_distance_z(-1.),
          nGenTrk(0),
          num_matched_reco_tracks(0),
          average_match_quality(0.0) {
      ptot.setPx(0);
      ptot.setPy(0);
      ptot.setPz(0);
      ptot.setE(0);
      p4 = LorentzVector(0, 0, 0, 0);
      r = sqrt(x * x + y * y);
    };
    double x, y, z, r, t;
    HepMC::FourVector ptot;
    LorentzVector p4;
    double pt;
    double ptsq;
    double closest_vertex_distance_z;
    int nGenTrk;
    int num_matched_reco_tracks;
    float average_match_quality;
    EncodedEventId eventId;
    TrackingVertexRef sim_vertex;
    int OriginalIndex = -1;

    unsigned int nwosmatch = 0;                    // number of recvertices dominated by this simevt (by wos)
    unsigned int nwntmatch = 0;                    // number of recvertices dominated by this simevt  (by nt)
    std::vector<unsigned int> wos_dominated_recv;  // list of dominated recv (by wos, size==nwosmatch)

    std::map<unsigned int, double> wnt;  // weighted number of tracks in recvtx (by index)
    std::map<unsigned int, double> wos;  // sum of wos in recvtx (by index) // oops -> this was int before 04-22
    double sumwos = 0;                   // sum of wos in any recvtx
    double sumwnt = 0;                   // sum of weighted tracks
    unsigned int rec = NOT_MATCHED;      // best match (NO_MATCH if not matched)
    unsigned int matchQuality = 0;       // quality flag

    void addTrack(unsigned int irecv, double twos, double twt) {
      sumwnt += twt;
      if (wnt.find(irecv) == wnt.end()) {
        wnt[irecv] = twt;
      } else {
        wnt[irecv] += twt;
      }

      sumwos += twos;
      if (wos.find(irecv) == wos.end()) {
        wos[irecv] = twos;
      } else {
        wos[irecv] += twos;
      }
    }
  };

  // auxiliary class holding reconstructed vertices
  struct recoPrimaryVertex {
    recoPrimaryVertex(double x1, double y1, double z1)
        : x(x1),
          y(y1),
          z(z1),
          pt(0),
          ptsq(0),
          closest_vertex_distance_z(-1.),
          nRecoTrk(0),
          num_matched_sim_tracks(0),
          ndof(0.),
          recVtx(nullptr) {
      r = sqrt(x * x + y * y);
    };
    double x, y, z, r;
    double pt;
    double ptsq;
    double closest_vertex_distance_z;
    int nRecoTrk;
    int num_matched_sim_tracks;
    double ndof;
    const reco::Vertex* recVtx;
    reco::VertexBaseRef recVtxRef;
    int OriginalIndex = -1;

    std::map<unsigned int, double> wos;  // simevent -> wos
    std::map<unsigned int, double> wnt;  // simevent -> weighted number of truth matched tracks
    unsigned int wosmatch;               // index of the simevent providing the largest contribution to wos
    unsigned int wntmatch;               // index of the simevent providing the highest number of tracks
    double sumwos = 0;                   // total sum of wos of all truth matched tracks
    double sumwnt = 0;                   // total weighted number of truth matchted tracks
    double maxwos = 0;                   // largest wos sum from one sim event (wosmatch)
    double maxwnt = 0;                   // largest weighted number of tracks from one sim event (ntmatch)
    int maxwosnt = 0;                    // number of tracks from the simevt with highest wos
    unsigned int sim = NOT_MATCHED;      // best match (NO_MATCH if not matched)
    unsigned int matchQuality = 0;       // quality flag

    bool is_real() { return (matchQuality > 0) && (matchQuality < 99); }

    bool is_fake() { return (matchQuality <= 0) || (matchQuality >= 99); }

    bool is_signal() { return (sim == 0); }

    int split_from() {
      if (is_real())
        return -1;
      if ((maxwos > 0) && (maxwos > 0.3 * sumwos))
        return wosmatch;
      return -1;
    }
    bool other_fake() { return (is_fake() && (split_from() < 0)); }

    void addTrack(unsigned int iev, double twos, double wt) {
      sumwnt += wt;
      if (wnt.find(iev) == wnt.end()) {
        wnt[iev] = wt;
      } else {
        wnt[iev] += wt;
      }

      sumwos += twos;
      if (wos.find(iev) == wos.end()) {
        wos[iev] = twos;
      } else {
        wos[iev] += twos;
      }
    }
  };

public:
  explicit Primary4DVertexValidation(const edm::ParameterSet&);
  ~Primary4DVertexValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void bookHistograms(DQMStore::IBooker& i, edm::Run const&, edm::EventSetup const&) override;

private:
  void matchReco2Sim(std::vector<recoPrimaryVertex>&,
                     std::vector<simPrimaryVertex>&,
                     const edm::ValueMap<float>&,
                     const edm::ValueMap<float>&,
                     const edm::Handle<reco::BeamSpot>&);
  bool matchRecoTrack2SimSignal(const reco::TrackBaseRef&);
  std::pair<const edm::Ref<std::vector<TrackingParticle>>*, int> getMatchedTP(const reco::TrackBaseRef&,
                                                                              const TrackingVertexRef&);
  double timeFromTrueMass(double, double, double, double);
  bool select(const reco::Vertex&, int level = 0);
  void observablesFromJets(const std::vector<reco::Track>&,
                           const std::vector<double>&,
                           const std::vector<int>&,
                           const std::string&,
                           unsigned int&,
                           double&,
                           double&,
                           double&,
                           double&);
  void isParticle(const reco::TrackBaseRef&,
                  const edm::ValueMap<float>&,
                  const edm::ValueMap<float>&,
                  const edm::ValueMap<float>&,
                  unsigned int&,
                  bool&,
                  bool&,
                  bool&,
                  bool&);
  void getWosWnt(const reco::Vertex&,
                 const reco::TrackBaseRef&,
                 const edm::ValueMap<float>&,
                 const edm::Handle<reco::BeamSpot>&,
                 double&,
                 double&);
  std::vector<Primary4DVertexValidation::simPrimaryVertex> getSimPVs(const edm::Handle<TrackingVertexCollection>&);
  std::vector<Primary4DVertexValidation::recoPrimaryVertex> getRecoPVs(const edm::Handle<edm::View<reco::Vertex>>&);

  void printMatchedRecoTrackInfo(const reco::Vertex&,
                                 const reco::TrackBaseRef&,
                                 const TrackingParticleRef&,
                                 const unsigned int&);
  void printSimVtxRecoVtxInfo(const struct Primary4DVertexValidation::simPrimaryVertex&,
                              const struct Primary4DVertexValidation::recoPrimaryVertex&);
  const bool mvaTPSel(const TrackingParticle&);
  const bool mvaRecSel(const reco::TrackBase&, const reco::Vertex&, const double&, const double&);

  // ----------member data ---------------------------

  const std::string folder_;
  static constexpr unsigned int NOT_MATCHED = 66666;
  static constexpr double simUnit_ = 1e9;     // sim time in s while reco time in ns
  static constexpr double c_ = 2.99792458e1;  // c in cm/ns
  static constexpr double mvaL_ = 0.5;        // MVA cuts for MVA categories
  static constexpr double mvaH_ = 0.8;
  static constexpr double selNdof_ = 4.;
  static constexpr double maxRank_ = 8.;
  static constexpr double maxTry_ = 10.;
  static constexpr double zWosMatchMax_ = 1.;
  static constexpr double etacutGEN_ = 4.;   // |eta| < 4;
  static constexpr double etacutREC_ = 3.;   // |eta| < 3;
  static constexpr double pTcut_ = 0.7;      // PT > 0.7 GeV
  static constexpr double deltaZcut_ = 0.1;  // dz separation 1 mm
  static constexpr double trackMaxBtlEta_ = 1.5;
  static constexpr double trackMinEtlEta_ = 1.6;
  static constexpr double trackMaxEtlEta_ = 3.;
  static constexpr double tol_ = 1.e-4;          // tolerance on reconstructed track time, [ns]
  static constexpr double minThrSumWnt_ = 0.01;  // min threshold for filling histograms with logarithmic scale
  static constexpr double minThrSumWos_ = 0.1;
  static constexpr double minThrSumPt_ = 0.01;
  static constexpr double minThrSumPt2_ = 1.e-3;
  static constexpr double minThrMetPt_ = 1.e-3;
  static constexpr double minThrSumPz_ = 1.e-4;

  static constexpr float c_cm_ns = geant_units::operators::convertMmToCm(CLHEP::c_light);  // [mm/ns] -> [cm/ns]

  bool use_only_charged_tracks_;
  bool debug_;
  bool optionalPlots_;
  bool use3dNoTime_;

  const double minProbHeavy_;
  const double trackweightTh_;
  const double mvaTh_;
  const std::vector<double> lineDensityPar_;
  const reco::RecoToSimCollection* r2s_;
  const reco::SimToRecoCollection* s2r_;

  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;

  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> vecPileupSummaryInfoToken_;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::BeamSpot> RecBeamSpotToken_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> Rec4DVerToken_;

  // --- START cluster debugging --- //
  edm::EDGetTokenT<MtdRecoClusterToSimLayerClusterAssociationMap> r2sClusterAssociationMapToken_;
  edm::EDGetTokenT<MtdSimLayerClusterToRecoClusterAssociationMap> s2rClusterAssociationMapToken_;
  edm::EDGetTokenT<reco::TPToSimCollectionMtd> TP2simCluAssociationMapToken_;
  // --- END cluster debugging  --- //

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> momentumToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> timeToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPToken_;
  edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> pdtToken_;

  // ---------------------
  // --- NTUPLE SAVING ---
  // ---------------------

  void createBranch();
  void reset();

  edm::Service<TFileService> fs;
  TTree* tree_;

  /* --- variables --- */

  // EVENT-LEVEL
  int run;
  int event;

  int nRecoTracks;
  int nRecoVtxs;
  int nSimTracks;
  int nSimVtxs;
  std::map<int, int> trackToVertexMap;
  std::map<int, int> recoToSimTrackMap;
  std::map<int, int> recoToSimVertexMap;
  std::map<int, int> simTrackToSimVertexMap;
  
  // TODO: add MC <-> reco mapping for track <-> vtx
  // think abt it (should we make a separate collection for reco/sim tks/vtxs? And then also make a map for MC<->reco for them? Seems a bit convoluted)

  // TRACK-LEVEL, RECO
  std::vector<int> track_id;

  std::vector<float> track_pt;
  std::vector<float> track_eta;
  std::vector<float> track_phi;

  std::vector<float> track_PCAx;
  std::vector<float> track_PCAy;
  std::vector<float> track_PCAz;

  std::vector<float> track_MTDx;
  std::vector<float> track_MTDy;
  std::vector<float> track_MTDz;

  std::vector<float> track_tmtd;
  std::vector<float> track_sigmatmtd;
  std::vector<float> track_t0;
  std::vector<float> track_sigmat0;

  std::vector<float> track_tofPi;
  std::vector<float> track_tofK;
  std::vector<float> track_tofP;
  std::vector<float> track_sigmatofPi;
  std::vector<float> track_sigmatofK;
  std::vector<float> track_sigmatofP;

  std::vector<float> track_probPi;
  std::vector<float> track_probK;
  std::vector<float> track_probP;

  std::vector<float> track_mtdMVAflag;
  std::vector<float> track_vtx3Dwgt;
  std::vector<float> track_matchCategory;
  std::vector<bool>  track_MVAselected;

  // TRACK-LEVEL, GEN
  std::vector<int> mc_track_id;

  std::vector<float> mc_track_pt;
  std::vector<float> mc_track_eta;
  std::vector<float> mc_track_phi;

  std::vector<float> mc_track_PCAx;
  std::vector<float> mc_track_PCAy;
  std::vector<float> mc_track_PCAz;

  std::vector<float> mc_track_MTDx;
  std::vector<float> mc_track_MTDy;
  std::vector<float> mc_track_MTDz;

  std::vector<float> mc_track_t0;
  std::vector<float> mc_track_tmtd;

  std::vector<float> mc_track_pdgId;

  // VERTEX-LEVEL, RECO
  std::vector<bool> is_signal_vtx;
  std::vector<int> vtx_id;
  std::vector<int> n_tracks_per_vertex;
  std::vector<float> vtx_x;
  std::vector<float> vtx_y;
  std::vector<float> vtx_z;
  std::vector<float> vtx_sigmaz;
  std::vector<float> vtx_t;
  std::vector<float> vtx_sigmat;

  // VERTEX-LEVEL, GEN
  std::vector<int> mc_vtx_id;
  std::vector<bool> mc_is_signal_vtx;
  std::vector<int> mc_n_tracks_per_vertex;
  std::vector<float> mc_vtx_t;
  std::vector<float> mc_vtx_x;
  std::vector<float> mc_vtx_y;
  std::vector<float> mc_vtx_z;
};

// constructors and destructor
Primary4DVertexValidation::Primary4DVertexValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      use_only_charged_tracks_(iConfig.getParameter<bool>("useOnlyChargedTracks")),
      debug_(iConfig.getUntrackedParameter<bool>("debug")),
      optionalPlots_(iConfig.getUntrackedParameter<bool>("optionalPlots")),
      use3dNoTime_(iConfig.getParameter<bool>("use3dNoTime")),
      minProbHeavy_(iConfig.getParameter<double>("minProbHeavy")),
      trackweightTh_(iConfig.getParameter<double>("trackweightTh")),
      mvaTh_(iConfig.getParameter<double>("mvaTh")),
      lineDensityPar_(iConfig.getParameter<std::vector<double>>("lineDensityPar")),
      pdtToken_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>()) {
  vecPileupSummaryInfoToken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag(std::string("addPileupInfo")));
  trackingParticleCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  trackingVertexCollectionToken_ = consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  simToRecoAssociationToken_ =
      consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  recoToSimAssociationToken_ =
      consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));

  // --- START cluster debugging --- //
  r2sClusterAssociationMapToken_ = consumes<MtdRecoClusterToSimLayerClusterAssociationMap>(
      iConfig.getParameter<edm::InputTag>("r2sClusterAssociationMapTag"));
  s2rClusterAssociationMapToken_ = consumes<MtdSimLayerClusterToRecoClusterAssociationMap>(
      iConfig.getParameter<edm::InputTag>("r2sClusterAssociationMapTag"));
  TP2simCluAssociationMapToken_ = consumes<reco::TPToSimCollectionMtd>(
      iConfig.getParameter<edm::InputTag>("TP2simCluAssociationMapTag"));
  // --- END cluster debuggings --- //


  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("mtdTracks"));
  RecBeamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBS"));
  Rec4DVerToken_ = consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offline4DPV"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  momentumToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("momentumSrc"));
  timeToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("timeSrc"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  sigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaSrc"));
  tofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofPi"));
  tofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofK"));
  tofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofP"));
  sigmatofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofPi"));
  sigmatofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofK"));
  sigmatofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofP"));
  probPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPi"));
  probKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probK"));
  probPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probP"));

  tree_ = fs->make<TTree>( "tree", "tree" );
  createBranch();
}

Primary4DVertexValidation::~Primary4DVertexValidation() {}

//
// member functions
//
void Primary4DVertexValidation::bookHistograms(DQMStore::IBooker& ibook,
                                               edm::Run const& iRun,
                                               edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);
}

bool Primary4DVertexValidation::matchRecoTrack2SimSignal(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return false;

  //// reco track matched to some TP from signal vertex
  for (const auto& tp : found->val) {
    if (tp.first->eventId().bunchCrossing() == 0 && tp.first->eventId().event() == 0)
      return true;
  }

  // reco track not matched to any TP from signal vertex
  return false;
}

std::pair<const edm::Ref<std::vector<TrackingParticle>>*, int> Primary4DVertexValidation::getMatchedTP(
    const reco::TrackBaseRef& recoTrack, const TrackingVertexRef& vsim) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP (fake tracks)
  if (found == r2s_->end())
    return std::make_pair(nullptr, -1);

  // matched TP equal to any TP of a given sim vertex
  for (const auto& tp : found->val) {
    if (std::find_if(vsim->daughterTracks_begin(), vsim->daughterTracks_end(), [&](const TrackingParticleRef& vtp) {
          return tp.first == vtp;
        }) != vsim->daughterTracks_end())
      return std::make_pair(&tp.first, 0);
    // matched TP not associated to any daughter track of a given sim vertex but having the same eventID (track from secondary vtx)
    else if (tp.first->eventId().bunchCrossing() == vsim->eventId().bunchCrossing() &&
             tp.first->eventId().event() == vsim->eventId().event()) {
      return std::make_pair(&tp.first, 1);
    }
    // matched TP not associated to any sim vertex of a given simulated event (PU track)
    else {
      return std::make_pair(&tp.first, 2);
    }
  }

  // reco track not matched to any TP from vertex
  return std::make_pair(nullptr, -1);
}

double Primary4DVertexValidation::timeFromTrueMass(double mass, double pathlength, double momentum, double time) {
  if (time > 0 && pathlength > 0 && mass > 0) {
    double gammasq = 1. + momentum * momentum / (mass * mass);
    double v = c_ * std::sqrt(1. - 1. / gammasq);  // cm / ns
    double t_est = time - (pathlength / v);

    return t_est;
  } else {
    return -1;
  }
}

bool Primary4DVertexValidation::select(const reco::Vertex& v, int level) {
  /* level
   0  !isFake  && ndof>4  (default)
   1  !isFake  && ndof>4 && prob > 0.01
   2  !isFake  && ndof>4 && prob > 0.01 && ptmax2 > 0.4
   */
  if (v.isFake())
    return false;
  if ((level == 0) && (v.ndof() > selNdof_))
    return true;
  /*if ((level == 1) && (v.ndof() > selNdof_) && (vertex_pxy(v) > 0.01))
    return true;
  if ((level == 2) && (v.ndof() > selNdof_) && (vertex_pxy(v) > 0.01) && (vertex_ptmax2(v) > 0.4))
    return true;
  if ((level == 3) && (v.ndof() > selNdof_) && (vertex_ptmax2(v) < 0.4))
    return true;*/
  return false;
}

void Primary4DVertexValidation::observablesFromJets(const std::vector<reco::Track>& reco_Tracks,
                                                    const std::vector<double>& mass_Tracks,
                                                    const std::vector<int>& category_Tracks,
                                                    const std::string& skip_Tracks,
                                                    unsigned int& n_Jets,
                                                    double& sum_EtJets,
                                                    double& sum_Pt2Jets,
                                                    double& met_Pt,
                                                    double& sum_PzJets) {
  double sum_PtJets = 0;
  n_Jets = 0;
  sum_EtJets = 0;
  sum_Pt2Jets = 0;
  met_Pt = 0;
  sum_PzJets = 0;
  auto met = LorentzVector(0, 0, 0, 0);
  std::vector<fastjet::PseudoJet> fjInputs_;
  fjInputs_.clear();
  size_t countScale0 = 0;
  for (size_t i = 0; i < reco_Tracks.size(); i++) {
    const auto recotr = reco_Tracks[i];
    const auto mass = mass_Tracks[i];
    float scale = 1.;
    if (recotr.charge() == 0) {
      continue;
    }
    // skip PU tracks in jet definition if skip_PU is required
    if (skip_Tracks == "skip_PU" && category_Tracks[i] == 2) {
      continue;
    }
    // skip fake tracks in jet definition if skip_Fake is required
    if (skip_Tracks == "skip_Fake" && category_Tracks[i] == -1) {
      continue;
    }
    if (recotr.pt() != 0) {
      scale = (recotr.pt() - recotr.ptError()) / recotr.pt();
    }
    if (edm::isNotFinite(scale)) {
      edm::LogWarning("Primary4DVertexValidation") << "Scaling is NAN ignoring this recotrack" << std::endl;
      scale = 0;
    }
    if (scale < 0) {
      scale = 0;
      countScale0++;
    }
    if (scale != 0) {
      fjInputs_.push_back(fastjet::PseudoJet(recotr.px() * scale,
                                             recotr.py() * scale,
                                             recotr.pz() * scale,
                                             std::sqrt(recotr.p() * recotr.p() + mass * mass) * scale));
    }
  }
  fastjet::ClusterSequence sequence(fjInputs_, fastjet::JetDefinition(fastjet::antikt_algorithm, 0.4));
  auto jets = fastjet::sorted_by_pt(sequence.inclusive_jets(0));
  for (const auto& pj : jets) {
    auto p4 = LorentzVector(pj.px(), pj.py(), pj.pz(), pj.e());
    sum_EtJets += std::sqrt(p4.e() * p4.e() - p4.P() * p4.P() + p4.pt() * p4.pt());
    sum_PtJets += p4.pt();
    sum_Pt2Jets += (p4.pt() * p4.pt() * 0.8 * 0.8);
    met += p4;
    sum_PzJets += p4.pz();
    n_Jets++;
  }
  met_Pt = met.pt();
  double metAbove = met_Pt - 2 * std::sqrt(sum_PtJets);
  if (metAbove > 0) {
    sum_Pt2Jets += (metAbove * metAbove);
  }
  if (countScale0 == reco_Tracks.size()) {
    sum_Pt2Jets = countScale0 * 0.01;  //leave some epsilon value to sort vertices with unknown pt
  }
}

void Primary4DVertexValidation::isParticle(const reco::TrackBaseRef& recoTrack,
                                           const edm::ValueMap<float>& probPi,
                                           const edm::ValueMap<float>& probK,
                                           const edm::ValueMap<float>& probP,
                                           unsigned int& no_PIDtype,
                                           bool& no_PID,
                                           bool& is_Pi,
                                           bool& is_K,
                                           bool& is_P) {
  no_PIDtype = 0;
  no_PID = false;
  is_Pi = false;
  is_K = false;
  is_P = false;
  if (probPi[recoTrack] == -1) {
    no_PIDtype = 1;
  } else if (edm::isNotFinite(probPi[recoTrack])) {
    no_PIDtype = 2;
  } else if (probPi[recoTrack] == 1 && probK[recoTrack] == 0 && probP[recoTrack] == 0) {
    no_PIDtype = 3;
  }
  no_PID = no_PIDtype > 0;
  is_Pi = !no_PID && 1. - probPi[recoTrack] < minProbHeavy_;
  is_K = !no_PID && !is_Pi && probK[recoTrack] > probP[recoTrack];
  is_P = !no_PID && !is_Pi && !is_K;
}

void Primary4DVertexValidation::getWosWnt(const reco::Vertex& recoVtx,
                                          const reco::TrackBaseRef& recoTrk,
                                          const edm::ValueMap<float>& sigmat0,
                                          const edm::Handle<reco::BeamSpot>& BS,
                                          double& wos,
                                          double& wnt) {
  double dz2_beam = pow((*BS).BeamWidthX() * cos(recoTrk->phi()) / tan(recoTrk->theta()), 2) +
                    pow((*BS).BeamWidthY() * sin(recoTrk->phi()) / tan(recoTrk->theta()), 2);
  double dz2 =
      pow(recoTrk->dzError(), 2) + dz2_beam + pow(0.0020, 2);  // added 20 um, some tracks have crazy small resolutions
  wos = recoVtx.trackWeight(recoTrk) / dz2;
  wnt = recoVtx.trackWeight(recoTrk) * std::min(recoTrk->pt(), 1.0);

  if (sigmat0[recoTrk] > 0) {
    double sigmaZ = (*BS).sigmaZ();
    double sigmaT = sigmaZ / c_;  // c in cm/ns
    wos = wos / erf(sigmat0[recoTrk] / sigmaT);
  }
}

/* Extract information form TrackingParticles/TrackingVertex and fill
 * the helper class simPrimaryVertex with proper generation-level
 * information */
std::vector<Primary4DVertexValidation::simPrimaryVertex> Primary4DVertexValidation::getSimPVs(
    const edm::Handle<TrackingVertexCollection>& tVC) {
  std::vector<Primary4DVertexValidation::simPrimaryVertex> simpv;
  int current_event = -1;
  int s = -1;
  for (TrackingVertexCollection::const_iterator v = tVC->begin(); v != tVC->end(); ++v) {
    // We keep only the first vertex from all the events at BX=0.
    if (v->eventId().bunchCrossing() != 0)
      continue;
    if (v->eventId().event() != current_event) {
      current_event = v->eventId().event();
    } else {
      continue;
    }
    s++;
    if (std::abs(v->position().z()) > 1000)
      continue;  // skip junk vertices

    // could be a new vertex, check  all primaries found so far to avoid multiple entries
    simPrimaryVertex sv(v->position().x(), v->position().y(), v->position().z(), v->position().t());
    sv.eventId = v->eventId();
    sv.sim_vertex = TrackingVertexRef(tVC, std::distance(tVC->begin(), v));
    sv.OriginalIndex = s;

    for (TrackingParticleRefVector::iterator iTrack = v->daughterTracks_begin(); iTrack != v->daughterTracks_end();
         ++iTrack) {
      assert((**iTrack).eventId().bunchCrossing() == 0);
    }
    simPrimaryVertex* vp = nullptr;  // will become non-NULL if a vertex is found and then point to it
    for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      if ((sv.eventId == v0->eventId) && (std::abs(sv.x - v0->x) < 1e-5) && (std::abs(sv.y - v0->y) < 1e-5) &&
          (std::abs(sv.z - v0->z) < 1e-5)) {
        vp = &(*v0);
        break;
      }
    }
    if (!vp) {
      // this is a new vertex, add it to the list of sim-vertices
      simpv.push_back(sv);
      vp = &simpv.back();
    }

    // Loop over daughter track(s) as Tracking Particles
    for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
      auto momentum = (*(*iTP)).momentum();
      const reco::Track* matched_best_reco_track = nullptr;
      double match_quality = -1;
      if (use_only_charged_tracks_ && (**iTP).charge() == 0)
        continue;
      if (s2r_->find(*iTP) != s2r_->end()) {
        matched_best_reco_track = (*s2r_)[*iTP][0].first.get();
        match_quality = (*s2r_)[*iTP][0].second;
      }

      vp->ptot.setPx(vp->ptot.x() + momentum.x());
      vp->ptot.setPy(vp->ptot.y() + momentum.y());
      vp->ptot.setPz(vp->ptot.z() + momentum.z());
      vp->ptot.setE(vp->ptot.e() + (**iTP).energy());
      vp->pt += (**iTP).pt();
      vp->ptsq += ((**iTP).pt() * (**iTP).pt());
      vp->nGenTrk++;

      if (matched_best_reco_track) {
        vp->num_matched_reco_tracks++;
        vp->average_match_quality += match_quality;
      }
    }  // End of for loop on daughters sim-particles
    if (vp->num_matched_reco_tracks)
      vp->average_match_quality /= static_cast<float>(vp->num_matched_reco_tracks);
    if (debug_) {
      edm::LogPrint("Primary4DVertexValidation")
          << "average number of associated tracks: " << vp->num_matched_reco_tracks / static_cast<float>(vp->nGenTrk)
          << " with average quality: " << vp->average_match_quality;
    }
  }  // End of for loop on tracking vertices

  // In case of no simulated vertices, break here
  if (simpv.empty())
    return simpv;

  // Now compute the closest distance in z between all simulated vertex
  // first initialize
  auto prev_z = simpv.back().z;
  for (simPrimaryVertex& vsim : simpv) {
    vsim.closest_vertex_distance_z = std::abs(vsim.z - prev_z);
    prev_z = vsim.z;
  }
  // then calculate
  for (std::vector<simPrimaryVertex>::iterator vsim = simpv.begin(); vsim != simpv.end(); vsim++) {
    std::vector<simPrimaryVertex>::iterator vsim2 = vsim;
    vsim2++;
    for (; vsim2 != simpv.end(); vsim2++) {
      double distance = std::abs(vsim->z - vsim2->z);
      // need both to be complete
      vsim->closest_vertex_distance_z = std::min(vsim->closest_vertex_distance_z, distance);
      vsim2->closest_vertex_distance_z = std::min(vsim2->closest_vertex_distance_z, distance);
    }
  }
  return simpv;
}

/* Extract information form recoVertex and fill the helper class
 * recoPrimaryVertex with proper reco-level information */
std::vector<Primary4DVertexValidation::recoPrimaryVertex> Primary4DVertexValidation::getRecoPVs(
    const edm::Handle<edm::View<reco::Vertex>>& tVC) {
  std::vector<Primary4DVertexValidation::recoPrimaryVertex> recopv;
  int r = -1;
  for (auto v = tVC->begin(); v != tVC->end(); ++v) {
    r++;
    // Skip junk vertices
    if (std::abs(v->z()) > 1000)
      continue;
    if (v->isFake() || !v->isValid())
      continue;

    recoPrimaryVertex sv(v->position().x(), v->position().y(), v->position().z());
    sv.recVtx = &(*v);
    sv.recVtxRef = reco::VertexBaseRef(tVC, std::distance(tVC->begin(), v));

    sv.OriginalIndex = r;
    sv.ndof = v->ndof();
    // this is a new vertex, add it to the list of reco-vertices
    recopv.push_back(sv);
    Primary4DVertexValidation::recoPrimaryVertex* vp = &recopv.back();

    // Loop over daughter track(s)
    for (auto iTrack = v->tracks_begin(); iTrack != v->tracks_end(); ++iTrack) {
      auto momentum = (*(*iTrack)).innerMomentum();
      if (momentum.mag2() == 0)
        momentum = (*(*iTrack)).momentum();
      vp->pt += std::sqrt(momentum.perp2());
      vp->ptsq += (momentum.perp2());
      vp->nRecoTrk++;

      auto matched = r2s_->find(*iTrack);
      if (matched != r2s_->end()) {
        vp->num_matched_sim_tracks++;
      }

    }  // End of for loop on daughters reconstructed tracks
  }    // End of for loop on tracking vertices

  // In case of no reco vertices, break here
  if (recopv.empty())
    return recopv;

  // Now compute the closest distance in z between all reconstructed vertex
  // first initialize
  auto prev_z = recopv.back().z;
  for (recoPrimaryVertex& vreco : recopv) {
    vreco.closest_vertex_distance_z = std::abs(vreco.z - prev_z);
    prev_z = vreco.z;
  }
  for (std::vector<recoPrimaryVertex>::iterator vreco = recopv.begin(); vreco != recopv.end(); vreco++) {
    std::vector<recoPrimaryVertex>::iterator vreco2 = vreco;
    vreco2++;
    for (; vreco2 != recopv.end(); vreco2++) {
      double distance = std::abs(vreco->z - vreco2->z);
      // need both to be complete
      vreco->closest_vertex_distance_z = std::min(vreco->closest_vertex_distance_z, distance);
      vreco2->closest_vertex_distance_z = std::min(vreco2->closest_vertex_distance_z, distance);
    }
  }
  return recopv;
}

// ------------ method called to produce the data  ------------
void Primary4DVertexValidation::matchReco2Sim(std::vector<recoPrimaryVertex>& recopv,
                                              std::vector<simPrimaryVertex>& simpv,
                                              const edm::ValueMap<float>& sigmat0,
                                              const edm::ValueMap<float>& MVA,
                                              const edm::Handle<reco::BeamSpot>& BS) {
  for (auto vv : simpv) {
    vv.wnt.clear();
    vv.wos.clear();
  }
  for (auto rv : recopv) {
    rv.wnt.clear();
    rv.wos.clear();
  }

  for (unsigned int iv = 0; iv < recopv.size(); iv++) {
    const reco::Vertex* vertex = recopv.at(iv).recVtx;

    for (unsigned int iev = 0; iev < simpv.size(); iev++) {
      double wnt = 0;
      double wos = 0;
      double evwnt = 0;
      double evwos = 0;
      double evnt = 0;

      for (auto iTrack = vertex->tracks_begin(); iTrack != vertex->tracks_end(); ++iTrack) {
      
        if (vertex->trackWeight(*iTrack) < trackweightTh_)
          continue;
        if (MVA[(*iTrack)] < mvaTh_)
          continue;

        auto tp_info = getMatchedTP(*iTrack, simpv.at(iev).sim_vertex).first;
        int matchCategory = getMatchedTP(*iTrack, simpv.at(iev).sim_vertex).second;
        // matched TP equal to any TP of a given sim vertex
        if (tp_info != nullptr && matchCategory == 0) {
          getWosWnt(*vertex, *iTrack, sigmat0, BS, wos, wnt);
          simpv.at(iev).addTrack(iv, wos, wnt);
          recopv.at(iv).addTrack(iev, wos, wnt);
          evwos += wos;
          evwnt += wnt;
          evnt++;
        }
      }  // RecoTracks loop

      // require 2 tracks for a wos-match
      if ((evwos > 0) && (evwos > recopv.at(iv).maxwos) && (evnt > 1)) {
        recopv.at(iv).wosmatch = iev;
        recopv.at(iv).maxwos = evwos;
        recopv.at(iv).maxwosnt = evnt;

        simpv.at(iev).wos_dominated_recv.push_back(iv);
        simpv.at(iev).nwosmatch++;
      }

      // weighted track counting match, require at least one track
      if ((evnt > 0) && (evwnt > recopv.at(iv).maxwnt)) {
        recopv.at(iv).wntmatch = iev;
        recopv.at(iv).maxwnt = evwnt;
      }
    }  // TrackingVertex loop

  }  // RecoPrimaryVertex

  // after filling infos, goes for the sim-reco match
  for (auto& vrec : recopv) {
    vrec.sim = NOT_MATCHED;
    vrec.matchQuality = 0;
  }
  unsigned int iev = 0;
  for (auto& vv : simpv) {
    if (debug_) {
      edm::LogPrint("Primary4DVertexValidation") << "iev: " << iev;
      edm::LogPrint("Primary4DVertexValidation") << "wos_dominated_recv.size: " << vv.wos_dominated_recv.size();
    }
    for (unsigned int i = 0; i < vv.wos_dominated_recv.size(); i++) {
      auto recov = vv.wos_dominated_recv.at(i);
      if (debug_) {
        edm::LogPrint("Primary4DVertexValidation")
            << "index of reco vertex: " << recov << " that has a wos: " << vv.wos.at(recov) << " at position " << i;
      }
    }
    vv.rec = NOT_MATCHED;
    vv.matchQuality = 0;
    iev++;
  }
  // this tries a one-to-one match, taking simPV with highest wos if there are > 1 simPV candidates
  for (unsigned int rank = 1; rank < maxRank_; rank++) {
    for (unsigned int iev = 0; iev < simpv.size(); iev++) {  //loop on SimPV
      if (simpv.at(iev).rec != NOT_MATCHED)
        continue;
      if (simpv.at(iev).nwosmatch == 0)
        continue;
      if (simpv.at(iev).nwosmatch > rank)
        continue;
      unsigned int iv = NOT_MATCHED;
      for (unsigned int k = 0; k < simpv.at(iev).wos_dominated_recv.size(); k++) {
        unsigned int rec = simpv.at(iev).wos_dominated_recv.at(k);
        auto vrec = recopv.at(rec);
        if (vrec.sim != NOT_MATCHED)
          continue;  // already matched
        if (std::abs(simpv.at(iev).z - vrec.z) > zWosMatchMax_)
          continue;  // insanely far away
        if ((iv == NOT_MATCHED) || simpv.at(iev).wos.at(rec) > simpv.at(iev).wos.at(iv)) {
          iv = rec;
        }
      }
      if (iv !=
          NOT_MATCHED) {  // if the rec vertex has already been associated is possible that iv remains NOT_MATCHED at this point
        recopv.at(iv).sim = iev;
        simpv.at(iev).rec = iv;
        recopv.at(iv).matchQuality = rank;
        simpv.at(iev).matchQuality = rank;
      }
    }
  }
  //give vertices a chance that have a lot of overlap, but are still recognizably
  //caused by a specific simvertex (without being classified as dominating)
  //like a small peak sitting on the flank of a larger nearby peak
  unsigned int ntry = 0;
  while (ntry++ < maxTry_) {
    unsigned nmatch = 0;
    for (unsigned int iev = 0; iev < simpv.size(); iev++) {
      if ((simpv.at(iev).rec != NOT_MATCHED) || (simpv.at(iev).wos.empty()))
        continue;
      // find a rec vertex for the NOT_MATCHED sim vertex
      unsigned int rec = NOT_MATCHED;
      for (auto rv : simpv.at(iev).wos) {
        if ((rec == NOT_MATCHED) || (rv.second > simpv.at(iev).wos.at(rec))) {
          rec = rv.first;
        }
      }

      if (rec == NOT_MATCHED) {  // try with wnt match
        for (auto rv : simpv.at(iev).wnt) {
          if ((rec == NOT_MATCHED) || (rv.second > simpv.at(iev).wnt.at(rec))) {
            rec = rv.first;
          }
        }
      }

      if (rec == NOT_MATCHED)
        continue;
      if (recopv.at(rec).sim != NOT_MATCHED)
        continue;  // already gone

      // check if the recvertex can be  matched
      unsigned int rec2sim = NOT_MATCHED;
      for (auto sv : recopv.at(rec).wos) {
        if (simpv.at(sv.first).rec != NOT_MATCHED)
          continue;  // already used
        if ((rec2sim == NOT_MATCHED) || (sv.second > recopv.at(rec).wos.at(rec2sim))) {
          rec2sim = sv.first;
        }
      }
      if (iev == rec2sim) {
        // do the match and assign lowest quality (i.e. max rank)
        recopv.at(rec).sim = iev;
        recopv.at(rec).matchQuality = maxRank_;
        simpv.at(iev).rec = rec;
        simpv.at(iev).matchQuality = maxRank_;
        nmatch++;
      }
    }  // sim loop
    if (nmatch == 0) {
      break;
    }
  }  // ntry
}

void Primary4DVertexValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // for ntuplizer saved values
  reset();

  using edm::Handle;
  using edm::View;
  using std::cout;
  using std::endl;
  using std::vector;
  using namespace reco;

  std::vector<float> pileUpInfo_z;

  // get the pileup information
  edm::Handle<std::vector<PileupSummaryInfo>> puinfoH;
  if (iEvent.getByToken(vecPileupSummaryInfoToken_, puinfoH)) {
    for (auto const& pu_info : *puinfoH.product()) {
      if (pu_info.getBunchCrossing() == 0) {
        pileUpInfo_z = pu_info.getPU_zpositions();
        break;
      }
    }
  }

  edm::Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByToken(trackingParticleCollectionToken_, TPCollectionH);
  if (!TPCollectionH.isValid())
    edm::LogWarning("Primary4DVertexValidation") << "TPCollectionH is not valid";

  edm::Handle<TrackingVertexCollection> TVCollectionH;
  iEvent.getByToken(trackingVertexCollectionToken_, TVCollectionH);
  if (!TVCollectionH.isValid())
    edm::LogWarning("Primary4DVertexValidation") << "TVCollectionH is not valid";

  edm::Handle<reco::SimToRecoCollection> simToRecoH;
  iEvent.getByToken(simToRecoAssociationToken_, simToRecoH);
  if (simToRecoH.isValid())
    s2r_ = simToRecoH.product();
  else
    edm::LogWarning("Primary4DVertexValidation") << "simToRecoH is not valid";

  edm::Handle<reco::RecoToSimCollection> recoToSimH;
  iEvent.getByToken(recoToSimAssociationToken_, recoToSimH);
  if (recoToSimH.isValid())
    r2s_ = recoToSimH.product();
  else
    edm::LogWarning("Primary4DVertexValidation") << "recoToSimH is not valid";

  edm::Handle<reco::BeamSpot> BeamSpotH;
  iEvent.getByToken(RecBeamSpotToken_, BeamSpotH);
  if (!BeamSpotH.isValid())
    edm::LogWarning("Primary4DVertexValidation") << "BeamSpotH is not valid";

  std::vector<simPrimaryVertex> simpv;  // a list of simulated primary MC vertices
  simpv = getSimPVs(TVCollectionH);
  
  // this bool check if first vertex in that with highest pT
  bool signal_is_highest_pt =
      std::max_element(simpv.begin(), simpv.end(), [](const simPrimaryVertex& lhs, const simPrimaryVertex& rhs) {
        return lhs.ptsq < rhs.ptsq;
      }) == simpv.begin();

  std::vector<recoPrimaryVertex> recopv;  // a list of reconstructed primary MC vertices
  edm::Handle<edm::View<reco::Vertex>> recVtxs;
  iEvent.getByToken(Rec4DVerToken_, recVtxs);
  if (!recVtxs.isValid())
    edm::LogWarning("Primary4DVertexValidation") << "recVtxs is not valid";
  recopv = getRecoPVs(recVtxs);

  const auto& trackAssoc = iEvent.get(trackAssocToken_);
  const auto& pathLength = iEvent.get(pathLengthToken_);
  const auto& momentum = iEvent.get(momentumToken_);
  const auto& time = iEvent.get(timeToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& t0Safe = iEvent.get(t0SafePidToken_);
  const auto& sigmat0Safe = iEvent.get(sigmat0SafePidToken_);
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& sigmatMtd = iEvent.get(sigmatmtdToken_);
  const auto& tofPi = iEvent.get(tofPiToken_);
  const auto& tofK = iEvent.get(tofKToken_);
  const auto& tofP = iEvent.get(tofPToken_);
  const auto& sigmatofPi = iEvent.get(sigmatofPiToken_);
  const auto& sigmatofK = iEvent.get(sigmatofKToken_);
  const auto& sigmatofP = iEvent.get(sigmatofPToken_);
  const auto& probPi = iEvent.get(probPiToken_);
  const auto& probK = iEvent.get(probKToken_);
  const auto& probP = iEvent.get(probPToken_);
  const auto& fPDGTable = iSetup.getHandle(pdtToken_);

  // ---- BEGIN reco to sim cluster association map -----
  const auto& r2sClusterAssociationMap = iEvent.get(r2sClusterAssociationMapToken_);
  const auto& s2rClusterAssociationMap = iEvent.get(s2rClusterAssociationMapToken_);
  const auto& TP2simCluAssociationMap = iEvent.get(TP2simCluAssociationMapToken_);
  // ----  END reco to sim cluster association map  -----

  // I have simPV and recoPV collections
  matchReco2Sim(recopv, simpv, sigmat0Safe, mtdQualMVA, BeamSpotH);

  // SAVING TOTAL NUMBER OF TRACKS
  int vtx_id_counter = 0;
  int sim_track_id_counter = 0;
  int track_id_counter = 0;

  int cumulative_saved_tracks = 0;
  int cumulative_sim_saved_tracks = 0;

  // std::cout << "NEW EVENT" << std::endl;
  // std::cout << "RECO VTXS = " << recopv.size() << std::endl;
  // std::cout << "SIM PRIMARY VTXS = " << simpv.size() << std::endl;

  // saving flag
  bool has_saved_sim = 0;
  int saved_simVtxs = 0; // only save if saved_simVtxs < simpv.size()

  // std::cout << "NEW EVENT" << std::endl;

  // RECO VERTEX LOOP
  for (unsigned int iv = 0; iv < recopv.size(); iv++) {
    if (recopv.at(iv).ndof > selNdof_) {
      const reco::Vertex* vertex = recopv.at(iv).recVtx;

      // std::cout << "RECO VTX " << iv << std::endl;

      // ----- FILLING RECO VERTEX QUANTITIES ----- //

      // first save all vtx (both MC-matched and not)
      vtx_id.push_back(vtx_id_counter++);
      // n_tracks_per_vertex.push_back(vertex->tracksSize());

      vtx_x.push_back(vertex->position().x());
      vtx_y.push_back(vertex->position().y());
      vtx_z.push_back(vertex->position().z());
      vtx_sigmaz.push_back(vertex->zError());

      vtx_t.push_back(vertex->t());
      vtx_sigmat.push_back(vertex->tError());

      // reco-sim quantities booking
      bool selectedVtxMatching = 0, selectedLV = 0, selectedLVMatching = 0;

      // saving flags
      bool has_saved_tracks = 0;
      int saved_tracks = 0;
      cumulative_sim_saved_tracks = 0; //bc of recopv repetition, reset

      // SIM VERTEX LOOP
      for (unsigned int iev = 0; iev < simpv.size(); iev++) {

        auto vsim = simpv.at(iev).sim_vertex;

        selectedVtxMatching = recopv.at(iv).sim == iev && simpv.at(iev).rec == iv;
        selectedLV = simpv.at(iev).eventId.bunchCrossing() == 0 && simpv.at(iev).eventId.event() == 0 &&
                          recopv.at(iv).OriginalIndex == 0;
        selectedLVMatching = selectedVtxMatching && selectedLV;  // bool for reco vtx leading match
        
        if (selectedLVMatching && !recopv.at(iv).is_signal()) {
          edm::LogWarning("Primary4DVertexValidation")
              << "Reco vtx leading match inconsistent: BX/ID " << simpv.at(iev).eventId.bunchCrossing() << " "
              << simpv.at(iev).eventId.event();
        }
        if (selectedLVMatching && debug_) {
          printSimVtxRecoVtxInfo(simpv.at(iev), recopv.at(iv));
        }
        double vzsim = simpv.at(iev).z;
        double vtsim = simpv.at(iev).t * simUnit_;

        double wnt = 0, wos = 0;
        double PUsumWnt = 0, PUsumWos = 0, SecsumWos = 0, FakesumWos = 0, PUsumPt = 0, PUsumPt2 = 0;
        double sumWnt = 0, sumWos = 0, sumPt = 0, sumPt2 = 0;
        unsigned int nt = 0, PUnt = 0, Fakent = 0;

        std::vector<double> massVector;
        std::vector<reco::Track> recotracks;
        std::vector<int> categoryVector;

        // if(selectedVtxMatching){
        //   std::cout << "SIM VTX HAS BEEN MATCHED TO RECO" << std::endl;
        //   std::cout << "vtx index = " << vtx_id_counter - 1 << std::endl;
        //   std::cout << "sim vtx index = " << iev << std::endl;
        //   std::cout << "t_sim = " << vtsim << ", t_reco = " << vertex->t() << std::endl;
        //   std::cout << "z_sim = " << vzsim << ", z_reco = " << vertex->z() << std::endl;
        // }

        // --- SAVE GEN VTX QUANTITIES ---
        if(!has_saved_sim){
          saved_simVtxs++;
          mc_vtx_id.push_back(iev);
          mc_is_signal_vtx.push_back(simpv.at(iev).eventId.bunchCrossing() == 0 && simpv.at(iev).eventId.event() == 0);
          mc_n_tracks_per_vertex.push_back(vsim->nDaughterTracks());                      
          mc_vtx_x.push_back(vsim->position().x());
          mc_vtx_y.push_back(vsim->position().y());
          mc_vtx_z.push_back(vsim->position().z());
          mc_vtx_t.push_back(vtsim);

          // SAVE SIM TRACK QUANTITIES
          // std::cout << "Saving " << vsim->nDaughterTracks() << " sim tracks for vertex idx " << iev << std::endl;
          for (auto iTrack = vsim->daughterTracks_begin(); iTrack != vsim->daughterTracks_end(); ++iTrack) {
            simTrackToSimVertexMap.insert({sim_track_id_counter, iev});
            mc_track_id.push_back(sim_track_id_counter++);
            mc_track_pt.push_back((*iTrack)->pt());
            mc_track_eta.push_back((*iTrack)->eta());
            mc_track_phi.push_back((*iTrack)->phi());
            mc_track_t0.push_back((*iTrack)->parentVertex()->position().t() * simUnit_);
            mc_track_PCAx.push_back((*iTrack)->parentVertex()->position().x());
            mc_track_PCAy.push_back((*iTrack)->parentVertex()->position().y());
            mc_track_PCAz.push_back((*iTrack)->parentVertex()->position().z());
            
            // // RETRIEVE CLUSTER
            // std::cout << "Retrieving sim cluster associated to TP" << std::endl;
          
            double simClusterTime = -999.;
            double simClusterX = -999., simClusterY = -999., simClusterZ = -999;

            if(TP2simCluAssociationMap.find(*iTrack) != TP2simCluAssociationMap.end()){
              const auto& simClusterRefs = TP2simCluAssociationMap[*iTrack];
              // std::cout << "RETRIEVED" << std::endl;

              for(const auto& simCluster : simClusterRefs){
                if(simClusterTime < -50.){
                    simClusterTime = (*simCluster).simLCTime(); //only save first
                    simClusterX = (*simCluster).simLCPos().x();
                    simClusterY = (*simCluster).simLCPos().y();
                    simClusterZ = (*simCluster).simLCPos().z();
                }
              }
            }
            // else {
            //   std::cout << "NOT RETRIEVED" << std::endl;
            // }
            
            mc_track_MTDx.push_back(simClusterX);
            mc_track_MTDy.push_back(simClusterY);
            mc_track_MTDz.push_back(simClusterZ);
            mc_track_tmtd.push_back(simClusterTime);

            mc_track_pdgId.push_back((*iTrack)->pdgId());
          }
        }

        // --- SAVE RECO-TO-GEN MATCHED VTX QUANTITIES ---
        if(selectedVtxMatching){
          recoToSimVertexMap.insert({vtx_id_counter - 1, iev});
        }

        // --------------------------------------- //

        // TRACK (from vertex) LOOP
        int track_tmp_counter = 0;

        for (auto iTrack = vertex->tracks_begin(); iTrack != vertex->tracks_end(); ++iTrack) {

          if (trackAssoc[*iTrack] == -1) {
            LogTrace("mtdTracks") << "Extended track not associated";
            continue;
          }

          // if (vertex->trackWeight(*iTrack) < trackweightTh_)
          //   continue;

          bool noCrack = std::abs((*iTrack)->eta()) < trackMaxBtlEta_ || std::abs((*iTrack)->eta()) > trackMinEtlEta_;

          bool selectRecoTrk = mvaRecSel(**iTrack, *vertex, t0Safe[*iTrack], sigmat0Safe[*iTrack]);

          auto tp_info = getMatchedTP(*iTrack, vsim).first;
          int matchCategory = getMatchedTP(*iTrack, vsim).second;

          // PU, fake and secondary tracks
          if (selectedVtxMatching) {
            unsigned int no_PIDtype = 0;
            bool no_PID, is_Pi, is_K, is_P;
            int PartID = 211;  // pion
            isParticle(*iTrack, probPi, probK, probP, no_PIDtype, no_PID, is_Pi, is_K, is_P);
            if (!use3dNoTime_) {
              if (no_PID || is_Pi) {
                PartID = 211;
              } else if (is_K) {
                PartID = 321;
              } else if (is_P) {
                PartID = 2212;
              }
            }
            const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID)));
            double mass = PData->mass().value();
            massVector.push_back(mass);
            recotracks.push_back(**iTrack);
            getWosWnt(*vertex, *iTrack, sigmat0Safe, BeamSpotH, wos, wnt);
            // reco track matched to any TP
            if (tp_info != nullptr) {
              if (debug_ && selectedLV) {
                printMatchedRecoTrackInfo(*vertex, *iTrack, *tp_info, matchCategory);
              }
              // matched TP not associated to any daughter track of a given sim vertex but having the same eventID (track from secondary vtx)
              if (matchCategory == 1) {
                categoryVector.push_back(matchCategory);
                SecsumWos += wos;
              }
              // matched TP not associated to any sim vertex of a given simulated event (PU track)
              if (matchCategory == 2) {
                if (optionalPlots_) {
                  // mePUTrackWnt_->Fill(wnt);
                  if (selectedLV) {
                    // mePUTrackRecLVWnt_->Fill(wnt);
                  }
                }
                PUsumWnt += wnt;
                PUsumWos += wos;
                PUsumPt += (*iTrack)->pt();
                PUsumPt2 += ((*iTrack)->pt() * (*iTrack)->pt());
                PUnt++;
                categoryVector.push_back(2);
              }
            }
            // reco track not matched to any TP (fake tracks)
            else {
              categoryVector.push_back(matchCategory);
              FakesumWos += wos;
              Fakent++;
            }
            nt++;
            sumWnt += wnt;
            sumWos += wos;
            sumPt += (*iTrack)->pt();
            sumPt2 += ((*iTrack)->pt() * (*iTrack)->pt());
          }

          // CLUSTER INFORMATION (default -999 if track not matched to TP)
          std::vector<float> clusterRes, clusterResAlt, clusterEnergy;
          double simClusterTime = -999., recoClusterTime = -999.;

          // matched TP equal to any TP of a given sim vertex

          // ---- START TRACK TO TP/CLUSTER MATCHING ---- //
          if (tp_info != nullptr && matchCategory == 0) {
            categoryVector.push_back(matchCategory);
            double mass = (*tp_info)->mass();
            double tsim = (*tp_info)->parentVertex()->position().t() * simUnit_;
            double tEst = timeFromTrueMass(mass, pathLength[*iTrack], momentum[*iTrack], time[*iTrack]);

            double xsim = (*tp_info)->parentVertex()->position().x();
            double ysim = (*tp_info)->parentVertex()->position().y();
            double zsim = (*tp_info)->parentVertex()->position().z();
            double xPCA = (*iTrack)->vx();
            double yPCA = (*iTrack)->vy();
            double zPCA = (*iTrack)->vz();

            double dZ = zPCA - zsim;
            double d3D = std::sqrt((xPCA - xsim) * (xPCA - xsim) + (yPCA - ysim) * (yPCA - ysim) + dZ * dZ);
            // orient d3D according to the projection of RECO - SIM onto simulated momentum
            if ((xPCA - xsim) * ((*tp_info)->px()) + (yPCA - ysim) * ((*tp_info)->py()) + dZ * ((*tp_info)->pz()) <
                0.) {
              d3D = -d3D;
            }

            // select TPs associated to the signal event
            bool selectTP = mvaTPSel(**tp_info);

            if (sigmat0Safe[*iTrack] == -1)
              continue;

            // 2) SIM CLUSTERS

            if(TP2simCluAssociationMap.find(*tp_info) != TP2simCluAssociationMap.end()){

              // STRATEGY 1:
              // use sim clusters associated to TP directly (using track id)

              const auto& simClusterRefs = TP2simCluAssociationMap[*tp_info];
              // std::cout << "Found " << simClusterRefs.size() << " sim clusters associated to TP" << std::endl;
              
              for(const auto& simCluster : simClusterRefs){
                if(simClusterTime < -50.){
                  simClusterTime = (*simCluster).simLCTime(); //only save first
                }

                // 1) save (reco time of matched hit - true time of true hit)
                //    i.e. true MTD time - reco MTD time
                //    (overall cluster resolution = reconstruction + association)
                clusterRes.push_back(tMtd[*iTrack] - simClusterTime);
                
                auto itp = s2rClusterAssociationMap.equal_range(simCluster);
                if (itp.first != itp.second) {
                  std::vector<edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster>> recoClusterRefs = (*itp.first).second;
                  for (unsigned int i = 0; i < recoClusterRefs.size(); i++) {
                    if(recoClusterTime < -50.) recoClusterTime = (*recoClusterRefs[i]).time();

                    // 2) save (reco time of true hit - true time of true hit)
                    //    i.e. reco time of true hit - true MTD time
                    //    (just reconstruction resolution)
                    //    alto save cluster energy
                    clusterEnergy.push_back((*recoClusterRefs[i]).energy());
                    clusterResAlt.push_back(recoClusterTime - simClusterTime);
                  }
                }
              }
            }

            // --- END retrieve matched clusters   --- //
          }  // if tp_info != nullptr && MatchCategory == 0
          // ---- END TRACK TO TP/CLUSTER MATCHING ---- //

          // ------ FILLING TRACK QUANTITIES ------ //
          // if(selectRecoTrk || 1) { //also iterate if vertex was not matched to SIM (it's just 1 anyway, no double counting)
          if(!has_saved_tracks){
            saved_tracks++;

            trackToVertexMap.insert({track_id_counter, vtx_id_counter - 1}); //vtx_id_counter already ++ before
            track_id.push_back(track_id_counter++);

            track_pt.push_back((*iTrack)->pt());
            track_eta.push_back((*iTrack)->eta());
            track_phi.push_back((*iTrack)->phi());
						std::cout << "saving track with phi = " << (*iTrack)->phi() << std::endl;

            track_PCAx.push_back((*iTrack)->vx());
            track_PCAy.push_back((*iTrack)->vy());
            track_PCAz.push_back((*iTrack)->vz());
						std::cout << "track pca (x,y,z) = " << (*iTrack)->vx() << ", " << (*iTrack)->vy() << ", " << (*iTrack)->vz() << std::endl;

            // from tracks -> reco hit -> average

            track_tmtd.push_back(tMtd[*iTrack]);
            track_sigmatmtd.push_back(sigmatMtd[*iTrack]);
            track_t0.push_back(t0Safe[*iTrack]);
            track_sigmat0.push_back(sigmat0Safe[*iTrack]);

            track_tofPi.push_back(tofPi[*iTrack]);
            track_tofK.push_back(tofK[*iTrack]);
            track_tofP.push_back(tofP[*iTrack]);
            track_sigmatofPi.push_back(sigmatofPi[*iTrack]);
            track_sigmatofK.push_back(sigmatofK[*iTrack]);
            track_sigmatofP.push_back(sigmatofP[*iTrack]);

            track_probPi.push_back(probPi[*iTrack]);
            track_probK.push_back(probK[*iTrack]);
            track_probP.push_back(probP[*iTrack]);

            track_mtdMVAflag.push_back(mtdQualMVA[*iTrack]);
            track_vtx3Dwgt.push_back(vertex->trackWeight(*iTrack));
            track_matchCategory.push_back(matchCategory);
            track_MVAselected.push_back(selectRecoTrk);


            // RETRIEVE CLUSTER INFO
            std::vector<double> recoClustersX, recoClustersY, recoClustersZ;
            double recoClusterX = -999., recoClusterY = -999., recoClusterZ = -999;
  
            // // retrieve extended track from current trackref
            // const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(RecTrackToken_), trackAssoc[*iTrack]);
            // const reco::Track& track = *mtdTrackref;
            
            // for (const auto hit : track.recHits()) {
            //   if (hit->isValid() == false)
            //     continue;
            //   MTDDetId Hit = hit->geographicalId();
            //   if ((Hit.det() == 6) && (Hit.subdetId() == 1)) {
            //     // TODO: hit->globalPosition() crashes the program, info is probably not available.
            //     //       find a way to overcome
            //     recoClustersX.push_back(hit->globalPosition().x());
            //     recoClustersY.push_back(hit->globalPosition().y());
            //     recoClustersZ.push_back(hit->globalPosition().z());
            //   }
            // }
            // if(recoClustersX.size() == 1){
            //   recoClusterX = recoClustersX[0];
            //   recoClusterY = recoClustersY[0];
            //   recoClusterZ = recoClustersZ[0]; 
            // } else if (recoClustersX.size() > 1) {
            //   // take average
            //   std::cout << "taking average" << std::endl;
            //   recoClusterX = std::accumulate(recoClustersX.begin(), recoClustersX.end(), 0.0) / (float)recoClustersX.size();
            //   recoClusterY = std::accumulate(recoClustersY.begin(), recoClustersY.end(), 0.0) / (float)recoClustersY.size();
            //   recoClusterZ = std::accumulate(recoClustersZ.begin(), recoClustersZ.end(), 0.0) / (float)recoClustersZ.size();
            // }

            track_MTDx.push_back(recoClusterX);
            track_MTDy.push_back(recoClusterY);
            track_MTDz.push_back(recoClusterZ);
          }

          // -------------------------------------- //

          if(tp_info != nullptr && (matchCategory == 0 || matchCategory == 1)){

            // // find index of matched TP in sim vertex track collection
            // std::cout << "Matched TP found: " << tp_info << ", SEARCHING FOR MATCHED TP INDEX" << std::endl;
            int tp_index = -1;

            // auto res = std::find_if(vsim->daughterTracks_begin(), vsim->daughterTracks_end(), [&](const TrackingParticleRef& vtp) {
            //       return vtp == *tp_info;
            //     });
            // if(res != vsim->daughterTracks_end()){
            //   std::cout << "Found matched TP from find_if" << std::endl;
            //   tp_index = std::distance(vsim->daughterTracks_begin(), res);
            //   break;
            // }

            for(unsigned int i = 0; i < vsim->daughterTracks().size(); i++){
              if(vsim->daughterTracks().at(i)->pt() == (*tp_info)->pt() && vsim->daughterTracks().at(i)->eta() == (*tp_info)->eta()){
                // std::cout << "Found matched TP from pT and eta equality" << std::endl;
                tp_index = i;
                break;
              }
            }

            if(tp_index >= 0){
              std::cout << "\t index check [tp_info vs index]: " << (*tp_info)->pt() << " vs " << mc_track_pt[cumulative_sim_saved_tracks + tp_index] << std::endl;
              recoToSimTrackMap.insert({cumulative_saved_tracks + track_tmp_counter, cumulative_sim_saved_tracks + tp_index});
              // update matchCategory for matched tracks
              track_matchCategory[cumulative_saved_tracks + track_tmp_counter] = matchCategory;
            }

          }
          
          track_tmp_counter++;
        } // loop on tracks

        has_saved_tracks = 1; //stop saving tracks after first simPV iteration

        // increase total counter
        cumulative_sim_saved_tracks += int(vsim->daughterTracks().size());

      }  // loop on simpv

      // don't save anymore
      has_saved_sim = 1;

      // FINAL SAVING QUANTITIES
      is_signal_vtx.push_back(recopv.at(iv).is_signal());

      // save total number of saved tracks PER VTX (only update at the end of recopv iteration)
      n_tracks_per_vertex.push_back(saved_tracks);
      cumulative_saved_tracks += saved_tracks;

    }    // ndof
  }      // loop on recopv

  // ------- FILLING EVENT-LEVEL -------
  run = iEvent.id().run();
  event = iEvent.id().event();

  nRecoTracks = cumulative_saved_tracks;
  nSimTracks = cumulative_sim_saved_tracks;
  nRecoVtxs = recopv.size();
  nSimVtxs = simpv.size();
  // ------------------------------------

  tree_->Fill();

}  // end of analyze

void Primary4DVertexValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Vertices");
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));

  // --- START cluster debugging --- //
  desc.add<edm::InputTag>("r2sClusterAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<edm::InputTag>("TP2simCluAssociationMapTag", edm::InputTag("mtdSimLayerClusterToTPAssociation"));
  // --- END cluster debugging  --- //

  desc.add<edm::InputTag>("mtdTracks", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("offlineBS", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("offline4DPV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("momentumSrc", edm::InputTag("trackExtenderWithMTD:generalTrackp"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("timeSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmaSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("tofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
  desc.add<edm::InputTag>("tofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
  desc.add<edm::InputTag>("tofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
  desc.add<edm::InputTag>("sigmatofPi", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
  desc.add<edm::InputTag>("sigmatofK", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
  desc.add<edm::InputTag>("sigmatofP", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
  desc.add<edm::InputTag>("probPi", edm::InputTag("tofPID:probPi"));
  desc.add<edm::InputTag>("probK", edm::InputTag("tofPID:probK"));
  desc.add<edm::InputTag>("probP", edm::InputTag("tofPID:probP"));
  desc.add<bool>("useOnlyChargedTracks", true);
  desc.addUntracked<bool>("debug", false);
  desc.addUntracked<bool>("optionalPlots", false);
  desc.add<bool>("use3dNoTime", false);
  desc.add<double>("trackweightTh", 0.5);
  desc.add<double>("mvaTh", 0.01);
  desc.add<double>("minProbHeavy", 0.75);

  //lineDensity parameters have been obtained by fitting the distribution of the z position of the vertices,
  //using a 200k single mu ptGun sample (gaussian fit)
  std::vector<double> lDP;
  lDP.push_back(1.87);
  lDP.push_back(0.);
  lDP.push_back(42.5);
  desc.add<std::vector<double>>("lineDensityPar", lDP);
  descriptions.add("vertices4DValid", desc);
}

void Primary4DVertexValidation::printMatchedRecoTrackInfo(const reco::Vertex& vtx,
                                                          const reco::TrackBaseRef& trk,
                                                          const TrackingParticleRef& tp,
                                                          const unsigned int& categ) {
  std::string strTrk;
  switch (categ) {
    case 0:
      strTrk = "Reco_Track:";
      break;
    case 1:
      strTrk = "SecRecoTrk:";
      break;
    case 2:
      strTrk = "PU_RecoTrk:";
      break;
  }
  edm::LogPrint("Primary4DVertexValidation")
      << strTrk << " w =" << std::setw(6) << std::setprecision(2) << vtx.trackWeight(trk) << " pt =" << std::setw(6)
      << std::setprecision(2) << trk->pt() << " eta =" << std::setw(6) << std::setprecision(2) << trk->eta()
      << "  MatchedTP: Pt =" << std::setw(6) << std::setprecision(2) << tp->pt() << " eta =" << std::setw(6)
      << std::setprecision(2) << tp->eta() << "  Parent vtx: z =" << std::setw(8) << std::setprecision(4)
      << tp->parentVertex()->position().z() << " t =" << std::setw(8) << std::setprecision(4)
      << tp->parentVertex()->position().t() * simUnit_ << " BX =" << tp->parentVertex()->eventId().bunchCrossing()
      << " ev =" << tp->parentVertex()->eventId().event() << std::endl;
}

void Primary4DVertexValidation::printSimVtxRecoVtxInfo(
    const struct Primary4DVertexValidation::simPrimaryVertex& simpVtx,
    const struct Primary4DVertexValidation::recoPrimaryVertex& recopVtx) {
  edm::LogPrint("Primary4DVertexValidation")
      << "Sim vtx (x,y,z,t) = (" << std::setprecision(4) << simpVtx.x << "," << std::setprecision(4) << simpVtx.y << ","
      << std::setprecision(4) << simpVtx.z << "," << std::setprecision(4) << simpVtx.t * simUnit_ << ")";
  edm::LogPrint("Primary4DVertexValidation")
      << "Sim vtx: pt = " << std::setprecision(4) << simpVtx.pt << " ptsq = " << std::setprecision(6) << simpVtx.ptsq
      << " nGenTrk = " << simpVtx.nGenTrk << " nmatch recotrks = " << simpVtx.num_matched_reco_tracks;
  edm::LogPrint("Primary4DVertexValidation")
      << "Reco vtx (x,y,z) = (" << std::setprecision(4) << recopVtx.x << "," << std::setprecision(4) << recopVtx.y
      << "," << std::setprecision(4) << recopVtx.z << ")";
  edm::LogPrint("Primary4DVertexValidation")
      << "Reco vtx: pt = " << std::setprecision(4) << recopVtx.pt << " ptsq = " << std::setprecision(6) << recopVtx.ptsq
      << " nrecotrks = " << recopVtx.nRecoTrk << " nmatch simtrks = " << recopVtx.num_matched_sim_tracks;
  edm::LogPrint("Primary4DVertexValidation") << "wnt " << recopVtx.sumwnt << " wos = " << recopVtx.sumwos;
  for (auto iTP = simpVtx.sim_vertex->daughterTracks_begin(); iTP != simpVtx.sim_vertex->daughterTracks_end(); ++iTP) {
    if (use_only_charged_tracks_ && (**iTP).charge() == 0) {
      continue;
    }
    edm::LogPrint("Primary4DVertexValidation")
        << "Daughter track of sim vertex: pt =" << std::setw(6) << std::setprecision(2) << (*iTP)->pt()
        << "  eta =" << std::setw(6) << std::setprecision(2) << (*iTP)->eta();
  }
}

const bool Primary4DVertexValidation::mvaTPSel(const TrackingParticle& tp) {
  bool match = false;
  if (tp.status() != 1) {
    return match;
  }
  match = tp.charge() != 0 && tp.pt() > pTcut_ && std::abs(tp.eta()) < etacutGEN_;
  return match;
}

const bool Primary4DVertexValidation::mvaRecSel(const reco::TrackBase& trk,
                                                const reco::Vertex& vtx,
                                                const double& t0,
                                                const double& st0) {
  bool match = false;
  match = trk.pt() > pTcut_ && std::abs(trk.eta()) < etacutREC_ && std::abs(trk.vz() - vtx.z()) <= deltaZcut_;
  if (st0 > 0.) {
    match = match && std::abs(t0 - vtx.t()) < 3. * st0;
  }
  return match;
}


// branch title creation
void Primary4DVertexValidation::createBranch() { 

  // Event-level
  tree_->Branch("run", &run);
  tree_->Branch("event", &event);
  tree_->Branch("n_reco_tracks", &nRecoTracks);
  tree_->Branch("n_reco_vtxs", &nRecoVtxs);
  tree_->Branch("n_sim_tracks", &nSimTracks);
  tree_->Branch("n_sim_vtxs", &nSimVtxs);
  tree_->Branch("track_to_vertex_map", &trackToVertexMap);
  tree_->Branch("reco_to_sim_vertex_map", &recoToSimVertexMap);
  tree_->Branch("reco_to_sim_track_map", &recoToSimTrackMap);
  tree_->Branch("sim_track_to_sim_vertex_map", &simTrackToSimVertexMap);

  // Track quantities, reco
  tree_->Branch("track_id", &track_id);

  tree_->Branch("track_pt", &track_pt);
  tree_->Branch("track_eta", &track_eta);
  tree_->Branch("track_phi", &track_phi);
  
  tree_->Branch("track_PCAx", &track_PCAx);
  tree_->Branch("track_PCAy", &track_PCAy);
  tree_->Branch("track_PCAz", &track_PCAz);

  tree_->Branch("track_MTDx", &track_MTDx);
  tree_->Branch("track_MTDy", &track_MTDy);
  tree_->Branch("track_MTDz", &track_MTDz);

  tree_->Branch("track_tmtd", &track_tmtd);
  tree_->Branch("track_sigmatmtd", &track_sigmatmtd);
  tree_->Branch("track_t0", &track_t0);
  tree_->Branch("track_sigmat0", &track_sigmat0);

  tree_->Branch("track_tofPi", &track_tofPi);
  tree_->Branch("track_tofK", &track_tofK);
  tree_->Branch("track_tofP", &track_tofP);
  tree_->Branch("track_sigmatofPi", &track_sigmatofPi);
  tree_->Branch("track_sigmatofK", &track_sigmatofK);
  tree_->Branch("track_sigmatofP", &track_sigmatofP);
    
  tree_->Branch("track_probPi", &track_probPi);
  tree_->Branch("track_probK", &track_probK);
  tree_->Branch("track_probP", &track_probP);
  tree_->Branch("track_mtdMVAflag", &track_mtdMVAflag);
  tree_->Branch("track_vtx3Dwgt", &track_vtx3Dwgt);
  tree_->Branch("track_match_category", &track_matchCategory);
  tree_->Branch("track_MVA_selected", &track_MVAselected);

  // Track quantities, MC
  tree_->Branch("mc_track_id", &mc_track_id);

  tree_->Branch("mc_track_pt", &mc_track_pt);
  tree_->Branch("mc_track_eta", &mc_track_eta);
  tree_->Branch("mc_track_phi", &mc_track_phi);
  tree_->Branch("mc_track_PCAx", &mc_track_PCAx);
  tree_->Branch("mc_track_PCAy", &mc_track_PCAy);
  tree_->Branch("mc_track_PCAz", &mc_track_PCAz);

  tree_->Branch("mc_track_MTDx", &mc_track_MTDx);
  tree_->Branch("mc_track_MTDy", &mc_track_MTDy);
  tree_->Branch("mc_track_MTDz", &mc_track_MTDz);

  tree_->Branch("mc_track_t0", &mc_track_t0);
  tree_->Branch("mc_track_tmtd", &mc_track_tmtd);
  tree_->Branch("mc_track_pdgId", &mc_track_pdgId);
  
  // Vertex quantities, reco
  tree_->Branch("vtx_id", &vtx_id);
  tree_->Branch("n_tracks_per_vertex", &n_tracks_per_vertex);
  tree_->Branch("is_signal_vtx", &is_signal_vtx);
  tree_->Branch("vtx_x", &vtx_x);
  tree_->Branch("vtx_y", &vtx_y);
  tree_->Branch("vtx_z", &vtx_z);
  tree_->Branch("vtx_sigmaz", &vtx_sigmaz);
  tree_->Branch("vtx_t", &vtx_t);
  tree_->Branch("vtx_sigmat", &vtx_sigmat);

  // Vertex quantities, MC
  tree_->Branch("mc_vtx_id", &mc_vtx_id);
  tree_->Branch("mc_is_signal_vtx", &mc_is_signal_vtx);
  tree_->Branch("mc_n_tracks_per_vertex", &mc_n_tracks_per_vertex);
  tree_->Branch("mc_vtx_x", &mc_vtx_x);
  tree_->Branch("mc_vtx_y", &mc_vtx_y);
  tree_->Branch("mc_vtx_z", &mc_vtx_z);
  tree_->Branch("mc_vtx_t", &mc_vtx_t);
}

// branch reset
void Primary4DVertexValidation::reset() {

  // Event-level
  run = -1;
  event = -1;
  nRecoTracks = -1;
  nRecoVtxs = -1;
  nSimVtxs = -1;
  nSimTracks = -1;
  recoToSimTrackMap.clear();
  recoToSimVertexMap.clear();
  trackToVertexMap.clear();
  simTrackToSimVertexMap.clear();

  // Track quantities, reco
  track_id.clear();
  track_pt.clear();
  track_eta.clear();
  track_phi.clear();
  track_PCAx.clear();
  track_PCAy.clear();
  track_PCAz.clear();
  track_MTDx.clear();
  track_MTDy.clear();
  track_MTDz.clear();
  track_tmtd.clear();
  track_sigmatmtd.clear();
  track_t0.clear();
  track_sigmat0.clear();
  track_tofPi.clear();
  track_tofK.clear();
  track_tofP.clear();
  track_sigmatofPi.clear();
  track_sigmatofK.clear();
  track_sigmatofP.clear();
  track_probPi.clear();
  track_probK.clear();
  track_probP.clear();
  track_mtdMVAflag.clear();
  track_vtx3Dwgt.clear();
  track_matchCategory.clear();
  track_MVAselected.clear();

  // Track quantities, MC
  mc_track_id.clear();
  mc_track_pt.clear();
  mc_track_eta.clear();
  mc_track_phi.clear();
  mc_track_PCAx.clear();
  mc_track_PCAy.clear();
  mc_track_PCAz.clear();
  mc_track_MTDx.clear();
  mc_track_MTDy.clear();
  mc_track_MTDz.clear();
  mc_track_t0.clear();
  mc_track_tmtd.clear();
  mc_track_pdgId.clear();

  // Vertex quantities, reco
  is_signal_vtx.clear();
  vtx_id.clear();
  n_tracks_per_vertex.clear();
  vtx_x.clear();
  vtx_y.clear();
  vtx_z.clear();
  vtx_sigmaz.clear();
  vtx_t.clear();
  vtx_sigmat.clear();

  // Vertex quantities, MC
  mc_vtx_id.clear();
  mc_is_signal_vtx.clear();
  mc_n_tracks_per_vertex.clear();
  mc_vtx_x.clear();
  mc_vtx_y.clear();
  mc_vtx_z.clear();
  mc_vtx_t.clear();

}

DEFINE_FWK_MODULE(Primary4DVertexValidation);
