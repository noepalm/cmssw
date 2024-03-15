#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

using namespace std;
using namespace edm;

class TOFPIDProducer : public edm::stream::EDProducer<> {
public:
  TOFPIDProducer(const ParameterSet& pset);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  template <class H, class T>
  void fillValueMap(edm::Event& iEvent,
                    const edm::Handle<H>& handle,
                    const std::vector<T>& vec,
                    const std::string& name) const;

  void produce(edm::Event& ev, const edm::EventSetup& es) final;

private:
  static constexpr char t0Name[] = "t0";
  static constexpr char sigmat0Name[] = "sigmat0";
  static constexpr char t0safeName[] = "t0safe";
  static constexpr char sigmat0safeName[] = "sigmat0safe";
  static constexpr char probPiName[] = "probPi";
  static constexpr char probKName[] = "probK";
  static constexpr char probPName[] = "probP";
  static constexpr char dtSignPiName[] = "dtSignPi";
  static constexpr char dtSignKName[] = "dtSignK";
  static constexpr char dtSignPName[] = "dtSignP";
  static constexpr char betaName[] = "beta";
  static constexpr char sigmaBetaName[] = "sigmabeta";
  static constexpr char timeChisqPiName[] = "timeChisqPi";
  static constexpr char timeChisqKName[] = "timeChisqK";
  static constexpr char timeChisqPName[] = "timeChisqP";

  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofkToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> tofpToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofpiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofkToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmatofpToken_; 
  edm::EDGetTokenT<edm::ValueMap<float>> betapiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> betakToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> betapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmabetapiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmabetakToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmabetapToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxsToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trackMTDTimeQualityToken_;
  const double vtxMaxSigmaT_;
  const double maxDz_;
  const double maxDtSignificance_;
  const double minProbHeavy_;
  const double fixedT0Error_;
  const double probPion_;
  const double probKaon_;
  const double probProton_;
  const double minTrackTimeQuality_;
  const bool MVASel_;
  const bool vertexReassignment_;
};

TOFPIDProducer::TOFPIDProducer(const ParameterSet& iConfig)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksSrc"))),
      t0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"))),
      tmtdToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtdSrc"))),
      sigmat0Token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"))),
      sigmatmtdToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtdSrc"))),
      tofkToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofkSrc"))),
      tofpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tofpSrc"))),
      sigmatofpiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofpiSrc"))),
      sigmatofkToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofkSrc"))),
      sigmatofpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatofpSrc"))),
      betapiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("betapiSrc"))),
      betakToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("betakSrc"))),
      betapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("betapSrc"))),
      sigmabetapiToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmabetapiSrc"))),
      sigmabetakToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmabetakSrc"))),
      sigmabetapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmabetapSrc"))),
      vtxsToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxsSrc"))),
      trackMTDTimeQualityToken_(
          consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"))),
      vtxMaxSigmaT_(iConfig.getParameter<double>("vtxMaxSigmaT")),
      maxDz_(iConfig.getParameter<double>("maxDz")),
      maxDtSignificance_(iConfig.getParameter<double>("maxDtSignificance")),
      minProbHeavy_(iConfig.getParameter<double>("minProbHeavy")),
      fixedT0Error_(iConfig.getParameter<double>("fixedT0Error")),
      probPion_(iConfig.getParameter<double>("probPion")),
      probKaon_(iConfig.getParameter<double>("probKaon")),
      probProton_(iConfig.getParameter<double>("probProton")),
      minTrackTimeQuality_(iConfig.getParameter<double>("minTrackTimeQuality")),
      MVASel_(iConfig.getParameter<bool>("MVASel")),
      vertexReassignment_(iConfig.getParameter<bool>("vertexReassignment")) {
  produces<edm::ValueMap<float>>(t0Name);
  produces<edm::ValueMap<float>>(sigmat0Name);
  produces<edm::ValueMap<float>>(t0safeName);
  produces<edm::ValueMap<float>>(sigmat0safeName);
  produces<edm::ValueMap<float>>(probPiName);
  produces<edm::ValueMap<float>>(probKName);
  produces<edm::ValueMap<float>>(probPName);
  produces<edm::ValueMap<float>>(dtSignPiName);
  produces<edm::ValueMap<float>>(dtSignKName);
  produces<edm::ValueMap<float>>(dtSignPName);
  produces<edm::ValueMap<float>>(timeChisqPiName);
  produces<edm::ValueMap<float>>(timeChisqKName);
  produces<edm::ValueMap<float>>(timeChisqPName);
  produces<edm::ValueMap<float>>(betaName);
  produces<edm::ValueMap<float>>(sigmaBetaName);
}

// Configuration descriptions
void TOFPIDProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracksSrc", edm::InputTag("generalTracks"))->setComment("Input tracks collection");
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"))
      ->setComment("Input ValueMap for track time at beamline");
  desc.add<edm::InputTag>("tmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"))
      ->setComment("Input ValueMap for track time at MTD");
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"))
      ->setComment("Input ValueMap for track time uncertainty at beamline");
  desc.add<edm::InputTag>("sigmatmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"))
      ->setComment("Input ValueMap for track time uncertainty at MTD");
  desc.add<edm::InputTag>("tofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"))
      ->setComment("Input ValueMap for track tof as kaon");
  desc.add<edm::InputTag>("tofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"))
      ->setComment("Input ValueMap for track tof as proton");
  desc.add<edm::InputTag>("sigmatofpiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"))
      ->setComment("Input ValueMap for track sigma(tof) as pion");
  desc.add<edm::InputTag>("sigmatofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"))
      ->setComment("Input ValueMap for track sigma(tof) as kaon");
  desc.add<edm::InputTag>("sigmatofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"))
      ->setComment("Input ValueMap for track sigma(tof) as proton");
  desc.add<edm::InputTag>("betapiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackBetaPi"))
      ->setComment("Input ValueMap for track beta as pion");
  desc.add<edm::InputTag>("betakSrc", edm::InputTag("trackExtenderWithMTD:generalTrackBetaK"))
      ->setComment("Input ValueMap for track beta as kaon");
  desc.add<edm::InputTag>("betapSrc", edm::InputTag("trackExtenderWithMTD:generalTrackBetaP"))
      ->setComment("Input ValueMap for track beta as proton");
  desc.add<edm::InputTag>("sigmabetapiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaBetaPi"))
      ->setComment("Input ValueMap for track sigma(beta) as pion");
  desc.add<edm::InputTag>("sigmabetakSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaBetaK"))
      ->setComment("Input ValueMap for track sigma(beta) as kaon");
  desc.add<edm::InputTag>("sigmabetapSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaBetaP"))
      ->setComment("Input ValueMap for track sigma(beta) as proton");
  desc.add<edm::InputTag>("vtxsSrc", edm::InputTag("unsortedOfflinePrimaryVertices4DwithPID"))
      ->setComment("Input primary vertex collection");
  desc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"))
      ->setComment("Track MVA quality value");
  desc.add<double>("vtxMaxSigmaT", 0.025)
      ->setComment("Maximum primary vertex time uncertainty for use in particle id [ns]");
  desc.add<double>("maxDz", 0.1)
      ->setComment("Maximum distance in z for track-primary vertex association for particle id [cm]");
  desc.add<double>("maxDtSignificance", 5.0)
      ->setComment(
          "Maximum distance in time (normalized by uncertainty) for track-primary vertex association for particle id");
  desc.add<double>("minProbHeavy", 0.75)
      ->setComment("Minimum probability for a particle to be a kaon or proton before reassigning the timestamp");
  desc.add<double>("fixedT0Error", 0.)->setComment("Use a fixed T0 uncertainty [ns]");
  desc.add<double>("probPion", 1.)->setComment("A priori probability pions");
  desc.add<double>("probKaon", 1.)->setComment("A priori probability kaons");
  desc.add<double>("probProton", 1.)->setComment("A priori probability for protons");
  desc.add<double>("minTrackTimeQuality", 0.8)->setComment("Minimum MVA Quality selection on tracks");
  desc.add<bool>("MVASel", false)->setComment("Use MVA Quality selection");
  desc.add<bool>("vertexReassignment", true)->setComment("Track-vertex reassignment");

  descriptions.add("tofPIDProducer", desc);
}

template <class H, class T>
void TOFPIDProducer::fillValueMap(edm::Event& iEvent,
                                  const edm::Handle<H>& handle,
                                  const std::vector<T>& vec,
                                  const std::string& name) const {
  auto out = std::make_unique<edm::ValueMap<T>>();
  typename edm::ValueMap<T>::Filler filler(*out);
  filler.insert(handle, vec.begin(), vec.end());
  filler.fill();
  iEvent.put(std::move(out), name);
}

void TOFPIDProducer::produce(edm::Event& ev, const edm::EventSetup& es) {

  edm::Handle<reco::TrackCollection> tracksH;
  ev.getByToken(tracksToken_, tracksH);
  const auto& tracks = *tracksH;

  const auto& t0In = ev.get(t0Token_);

  const auto& tmtdIn = ev.get(tmtdToken_);

  const auto& sigmat0In = ev.get(sigmat0Token_);

  const auto& sigmatmtdIn = ev.get(sigmatmtdToken_);

  const auto& tofkIn = ev.get(tofkToken_);

  const auto& tofpIn = ev.get(tofpToken_);

  const auto& sigmatofpiIn = ev.get(sigmatofpiToken_);

  const auto& sigmatofkIn = ev.get(sigmatofkToken_);

  const auto& sigmatofpIn = ev.get(sigmatofpToken_);

  const auto& betapiIn = ev.get(betapiToken_);

  const auto& betakIn = ev.get(betakToken_);

  const auto& betapIn = ev.get(betapToken_);

  const auto& sigmabetapiIn = ev.get(sigmabetapiToken_);

  const auto& sigmabetakIn = ev.get(sigmabetakToken_);

  const auto& sigmabetapIn = ev.get(sigmabetapToken_);

  const auto& vtxs = ev.get(vtxsToken_);

  const auto& trackMVAQualIn = ev.get(trackMTDTimeQualityToken_);

  //output value maps (PID probabilities and recalculated time at beamline)
  std::vector<float> t0OutRaw;
  std::vector<float> sigmat0OutRaw;
  std::vector<float> t0safeOutRaw;
  std::vector<float> sigmat0safeOutRaw;
  std::vector<float> probPiOutRaw;
  std::vector<float> probKOutRaw;
  std::vector<float> probPOutRaw;
  std::vector<float> dtsignPiOutRaw;
  std::vector<float> dtsignKOutRaw;
  std::vector<float> dtsignPOutRaw;
  std::vector<float> chisqPiOutRaw;
  std::vector<float> chisqKOutRaw;
  std::vector<float> chisqPOutRaw;
  std::vector<float> betaOutRaw;
  std::vector<float> sigmaBetaOutRaw;

  //Do work here
  for (unsigned int itrack = 0; itrack < tracks.size(); ++itrack) {
    const reco::Track& track = tracks[itrack];
    const reco::TrackRef trackref(tracksH, itrack);
    float t0 = t0In[trackref];
    float t0safe = t0;
    float sigmat0safe = sigmat0In[trackref];
    float sigmatmtd = (sigmatmtdIn[trackref] > 0. && fixedT0Error_ > 0.) ? fixedT0Error_ : sigmatmtdIn[trackref];
    float sigmat0 = sigmatmtd;

    float sigmatofpi = sigmatofpiIn[trackref];
    float sigmatofk = sigmatofkIn[trackref];
    float sigmatofp = sigmatofpIn[trackref];

    float betapi = betapiIn[trackref];
    float betak = betakIn[trackref];
    float betap = betapIn[trackref];

    float sigmabetapi = sigmabetapiIn[trackref];
    float sigmabetak = sigmabetakIn[trackref];
    float sigmabetap = sigmabetapIn[trackref];

    float prob_pi = -1.;
    float prob_k = -1.;
    float prob_p = -1.;

    float dtsign_pi_out = -1.;
    float dtsign_k_out = -1.;
    float dtsign_p_out = -1.;

    float chisqmin_pi_out = -1.;
    float chisqmin_k_out = -1.;
    float chisqmin_p_out = -1.;

    float beta = betapi;
    float sigmabeta = betapi - betap;

    float trackMVAQual = trackMVAQualIn[trackref];

    if (sigmat0 > 0. && (!MVASel_ || (MVASel_ && trackMVAQual >= minTrackTimeQuality_))) {
      double rsigmazsq = 1. / track.dzError() / track.dzError();
      // double rsigmat = 1. / sigmatmtd;

      double rsigmat[3] = {1. / std::sqrt(sigmatmtd*sigmatmtd + sigmatofpi*sigmatofpi),
                           1. / std::sqrt(sigmatmtd*sigmatmtd + sigmatofk*sigmatofk),
                           1. / std::sqrt(sigmatmtd*sigmatmtd + sigmatofp*sigmatofp)};

      //find associated vertex
      int vtxidx = -1;
      int vtxidxmindz = -1;
      int vtxidxminchisq = -1;
      double mindz = maxDz_;
      double minchisq = std::numeric_limits<double>::max();
      //first try based on association weights, but keep track of closest in z and z-t as well
      for (unsigned int ivtx = 0; ivtx < vtxs.size(); ++ivtx) {
        const reco::Vertex& vtx = vtxs[ivtx];
        float w = vtx.trackWeight(trackref);
        if (w > 0.5) {
          vtxidx = ivtx;
          break;
        }
        double dz = std::abs(track.dz(vtx.position()));
        if (dz < mindz) {
          mindz = dz;
          vtxidxmindz = ivtx;
        }
        if (vtx.tError() > 0. && vtx.tError() < vtxMaxSigmaT_) {
          double dt = std::abs(t0 - vtx.t());
          double dtsig = dt * rsigmat[0]; //pion hp. uncertainty
          double chisq = dz * dz * rsigmazsq + dtsig * dtsig;
          if (dz < maxDz_ && dtsig < maxDtSignificance_ && chisq < minchisq) {
            minchisq = chisq;
            vtxidxminchisq = ivtx;
          }
        }
      }

      //if no vertex found based on association weights, fall back to closest in z or z-t
      if (vtxidx < 0) {
        //if closest vertex in z does not have valid time information, just use it,
        //otherwise use the closest vertex in z-t plane with timing info, with a fallback to the closest in z
        if (vtxidxmindz >= 0 && !(vtxs[vtxidxmindz].tError() > 0. && vtxs[vtxidxmindz].tError() < vtxMaxSigmaT_)) {
          vtxidx = vtxidxmindz;
        } else if (vtxidxminchisq >= 0) {
          vtxidx = vtxidxminchisq;
        } else if (vtxidxmindz >= 0) {
          vtxidx = vtxidxmindz;
        }
      }

      //testing mass hypotheses only possible if there is an associated vertex with time information
      if (vtxidx >= 0 && vtxs[vtxidx].tError() > 0. && vtxs[vtxidx].tError() < vtxMaxSigmaT_) {
        //compute chisq in z-t plane for nominal vertex and mass hypothesis (pion)
        const reco::Vertex& vtxnom = vtxs[vtxidx];
        double dznom = std::abs(track.dz(vtxnom.position()));
        double dtnom = std::abs(t0 - vtxnom.t());
        double dtsignom = dtnom * rsigmat[0];
        dtsign_pi_out = dtsignom;
        double chisqnom = dznom * dznom * rsigmazsq + dtsignom * dtsignom;

        //recompute t0 for alternate mass hypotheses
        double t0_best = t0;

        //reliable match, revert to raw mtd time uncertainty + tof uncertainty for pion hp
        if (dtsignom < maxDtSignificance_) {
          sigmat0safe = 1./rsigmat[0];
          sigmabeta = sigmabetapi;
        }

        double tmtd = tmtdIn[trackref];
        double t0_k = tmtd - tofkIn[trackref];
        double t0_p = tmtd - tofpIn[trackref];

        double chisqmin = chisqnom;

        double chisqmin_pi = chisqnom;
        double chisqmin_k = std::numeric_limits<double>::max();
        double chisqmin_p = std::numeric_limits<double>::max();
        //loop through vertices and check for better matches
        for (unsigned int ivtx = 0; ivtx < vtxs.size(); ++ivtx) {
          const reco::Vertex& vtx = vtxs[ivtx];
          if (!vertexReassignment_) {
            if (ivtx != (unsigned int)vtxidx)
              continue;
          }
          if (!(vtx.tError() > 0. && vtx.tError() < vtxMaxSigmaT_)) {
            continue;
          }

          double dz = std::abs(track.dz(vtx.position()));
          if (dz >= maxDz_) {
            continue;
          }

          double chisqdz = dz * dz * rsigmazsq;

          double dt_k = std::abs(t0_k - vtx.t());
          double dtsig_k = dt_k * rsigmat[1];
          dtsign_k_out = dtsig_k;
          double chisq_k = chisqdz + dtsig_k * dtsig_k;

          if (dtsig_k < maxDtSignificance_ && chisq_k < chisqmin_k) {
            chisqmin_k = chisq_k;
          }

          double dt_p = std::abs(t0_p - vtx.t());
          double dtsig_p = dt_p * rsigmat[2];
          dtsign_p_out = dtsig_p;
          double chisq_p = chisqdz + dtsig_p * dtsig_p;

          if (dtsig_p < maxDtSignificance_ && chisq_p < chisqmin_p) {
            chisqmin_p = chisq_p;
          }

          if (dtsig_k < maxDtSignificance_ && chisq_k < chisqmin) {
            chisqmin = chisq_k;
            t0_best = t0_k;
            t0safe = t0_k;
            sigmat0safe = 1./rsigmat[1];
            beta = betak;
            sigmabeta = sigmabetak;
          }

          if (dtsig_p < maxDtSignificance_ && chisq_p < chisqmin) {
            chisqmin = chisq_p;
            t0_best = t0_p;
            t0safe = t0_p;
            sigmat0safe = 1./rsigmat[2];
            beta = betap;
            sigmabeta = sigmabetap;
          }

        }

        //compute PID probabilities
        //*TODO* deal with heavier nucleons and/or BSM case here?
        double rawprob_pi = probPion_ * exp(-0.5 * chisqmin_pi);
        double rawprob_k = probKaon_ * exp(-0.5 * chisqmin_k);
        double rawprob_p = probProton_ * exp(-0.5 * chisqmin_p);

        double normprob = 1. / (rawprob_pi + rawprob_k + rawprob_p);

        prob_pi = rawprob_pi * normprob;
        prob_k = rawprob_k * normprob;
        prob_p = rawprob_p * normprob;

        double prob_heavy = 1. - prob_pi;

        if (prob_heavy > minProbHeavy_) {
          t0 = t0_best;
        }

        chisqmin_pi_out = chisqmin_pi;
        chisqmin_k_out = chisqmin_k;
        chisqmin_p_out = chisqmin_p;
      }
    }

    t0OutRaw.push_back(t0);
    sigmat0OutRaw.push_back(sigmat0);
    t0safeOutRaw.push_back(t0safe);
    sigmat0safeOutRaw.push_back(sigmat0safe);
    probPiOutRaw.push_back(prob_pi);
    probKOutRaw.push_back(prob_k);
    probPOutRaw.push_back(prob_p);
    dtsignPiOutRaw.push_back(dtsign_pi_out);
    dtsignKOutRaw.push_back(dtsign_k_out);
    dtsignPOutRaw.push_back(dtsign_p_out);
    chisqPiOutRaw.push_back(chisqmin_pi_out);
    chisqKOutRaw.push_back(chisqmin_k_out);
    chisqPOutRaw.push_back(chisqmin_p_out);
    betaOutRaw.push_back(beta);
    sigmaBetaOutRaw.push_back(sigmabeta);
  }

  fillValueMap(ev, tracksH, t0OutRaw, t0Name);
  fillValueMap(ev, tracksH, sigmat0OutRaw, sigmat0Name);
  fillValueMap(ev, tracksH, t0safeOutRaw, t0safeName);
  fillValueMap(ev, tracksH, sigmat0safeOutRaw, sigmat0safeName);
  fillValueMap(ev, tracksH, probPiOutRaw, probPiName);
  fillValueMap(ev, tracksH, probKOutRaw, probKName);
  fillValueMap(ev, tracksH, probPOutRaw, probPName);
  fillValueMap(ev, tracksH, dtsignPiOutRaw, dtSignPiName);
  fillValueMap(ev, tracksH, dtsignKOutRaw, dtSignKName);
  fillValueMap(ev, tracksH, dtsignPOutRaw, dtSignPName);
  fillValueMap(ev, tracksH, chisqPiOutRaw, timeChisqPiName);
  fillValueMap(ev, tracksH, chisqKOutRaw, timeChisqKName);
  fillValueMap(ev, tracksH, chisqPOutRaw, timeChisqPName);
  fillValueMap(ev, tracksH, betaOutRaw, betaName);
  fillValueMap(ev, tracksH, sigmaBetaOutRaw, sigmaBetaName);
}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TOFPIDProducer);
