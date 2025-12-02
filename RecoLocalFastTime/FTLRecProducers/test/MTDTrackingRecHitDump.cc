#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"
#include <iostream>

class MTDTrackingRecHitDump : public edm::one::EDAnalyzer<> {
public:
  explicit MTDTrackingRecHitDump(const edm::ParameterSet& cfg)
      : src_(cfg.getParameter<edm::InputTag>("src")), token_(consumes<MTDTrackingDetSetVector>(src_)) {}
  ~MTDTrackingRecHitDump() override = default;

  void analyze(const edm::Event& e, const edm::EventSetup&) override {
    edm::Handle<MTDTrackingDetSetVector> h;
    e.getByToken(token_, h);
    if (!h.isValid()) {
      std::cout << "MTDTrackingRecHitDump: no collection " << src_ << " in this event\n";
      return;
    }
    size_t total = 0;
    for (const auto& detset : *h) {
      const unsigned int detId = detset.detId();
      for (const auto& hit : detset) {
        ++total;
        auto lp = hit.localPosition();
        std::cout << "DetId " << detId << "  pos = (" << lp.x() << "," << lp.y() << ")"
                  << "  err_xx=" << hit.localPositionError().xx() << "  time=" << hit.time()
                  << "  timeErr=" << hit.timeError() << "  energy=" << hit.energy() << "\n";
      }
    }
    std::cout << "MTDTrackingRecHitDump: total hits = " << total << "\n";
  }

private:
  const edm::InputTag src_;
  const edm::EDGetTokenT<MTDTrackingDetSetVector> token_;
};

DEFINE_FWK_MODULE(MTDTrackingRecHitDump);