#ifndef RecoMuon_L3TrackFinder_MuonChamberMatch_h
#define RecoMuon_L3TrackFinder_MuonChamberMatch_h

// Matching of a tracker track to the muon chambers (segments and rec-hits)
// for the muon high-purity track selectors: the raw quantities of the
// candidates of every muon station.
//
// A muon of charge q, transverse momentum pT and direction (eta, phi) at its
// point of closest approach leaves in each station a segment (or rec-hits)
// displaced from the track direction by the bending in the magnetic field.
// For a (track, object) pair
//
//   bend    = q * dphi(position)  * pT   [rad GeV]
//   dirBend = q * dphi(direction) * pT   (segments)
//   dEta    = eta(object) - eta(track)
//
// bend and dirBend are negative for an object displaced on the side the
// charge bends to; the factor pT removes most of the momentum dependence.
// Where a muon's own objects lie in these quantities (per detector, station
// and |eta|) is neither tabulated nor fitted: the forest gets the raw
// quantities and learns it.
//
// Candidates of a track in a station: the objects inside the search window
// (|dEta| < kChamberMaxDEta, bend on the charge side and not beyond
// kChamberMaxBend; rec-hits within kChamberMaxHitBx bunch crossings) whose
// direction, for segments, is also bent on the charge side. They are ordered
// by |dEta|, the coordinate the field does not bend, and the first two are
// kept. The search window is a wide road, not a working point.
//
// This is the C++ counterpart of chamber_match.py of the training repository
// (muonHighPurityTrackSelection/production): same window, same ordering of
// the candidates, inputs rounded to float as the training n-tuple stores them
// and arithmetic in double. Any change here needs the same change there and
// a retraining.

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace muonhp {

  // Muon stations: ME0 segments, DT or CSC segments of stations 1-4, GE1/1 and GE2/1 rec-hits (the HLT has no
  // segments there), RPC rec-hits of stations 1-4 (barrel or endcap).
  enum ChamberStation : uint8_t {
    kME0 = 0,
    kSt1,
    kSt2,
    kSt3,
    kSt4,
    kGE11,
    kGE21,
    kRPC1,
    kRPC2,
    kRPC3,
    kRPC4,
    kNChamberStations
  };

  // Search window
  inline constexpr double kChamberMaxDEta = 0.3;
  inline constexpr double kChamberMaxBend = 8.0;  // |q * dphi * pT| [rad GeV]
  inline constexpr int kChamberMaxHitBx = 1;      // rec-hits: |bx| <= 1

  // Values of a station without a candidate: outside the range of any candidate (bend < 0, bendDiff < kChamberMaxBend,
  // |dEta| < kChamberMaxDEta), so their exact value is irrelevant to a forest.
  inline constexpr float kNoCandidateBend = 9.f;
  inline constexpr float kNoCandidateDEta = 0.5f;

  // A segment or rec-hit in global coordinates, as the training n-tuple stores it.
  struct ChamberPoint {
    float eta = 0.f;
    float phi = 0.f;
    float dirPhi = 0.f;  // segments only
    uint8_t station = kME0;
    bool isSegment = false;
  };

  // Transverse radius and pseudorapidity of a point from its float coordinates, every operation rounded to float.
  // The training n-tuple tables evaluate GlobalPoint::perp() and eta() through the ROOT interpreter, which does
  // not fuse x * x + y * y into a multiply-add as the compiler does here: the volatile products keep the two
  // roundings, so that the values are those the forest was trained on to the last bit.
  inline float chamberPerp(float x, float y) {
    volatile float xx = x * x;
    volatile float yy = y * y;
    return std::sqrt(xx + yy);
  }
  inline float chamberEta(float x, float y, float z) { return ::asinhf(z / chamberPerp(x, y)); }

  // phi(a) - phi(b), reduced to [-pi, pi] only when outside (the reduction of the training code)
  inline double chamberDeltaPhi(double a, double b) {
    const double d = a - b;
    return std::abs(d) <= M_PI ? d : d - std::rint(d * (1.0 / (2.0 * M_PI))) * (2.0 * M_PI);
  }

  // Track parameters as stored in the training n-tuple (float), in double.
  struct ChamberMatchTrack {
    double pt, eta, phi, charge;
  };

  struct ChamberCandidate {
    float bend = kNoCandidateBend;
    float dEta = kNoCandidateDEta;
    float bendDiff = kNoCandidateBend;  // dirBend - bend (segments): how much more the direction is bent than the position
  };

  // The two candidates nearest in eta of one station.
  struct StationCandidates {
    ChamberCandidate nearest;
    ChamberCandidate second;
  };

  using ChamberMatch = std::array<StationCandidates, kNChamberStations>;

  // `points`: the event's segments and hits in the order DT, CSC, ME0, GEM, RPC barrel, RPC endcap, each in the
  // order of its collection (the order of the training n-tuple tables): of candidates at the same |dEta| the first
  // one is the nearer.
  inline ChamberMatch matchChambers(const ChamberMatchTrack& t, const std::vector<ChamberPoint>& points) {
    ChamberMatch match;
    std::array<double, kNChamberStations> nearest, second;
    nearest.fill(std::numeric_limits<double>::max());
    second.fill(std::numeric_limits<double>::max());
    for (const auto& p : points) {
      const double dEta = static_cast<double>(p.eta) - t.eta;
      const double absDEta = std::abs(dEta);
      if (!(absDEta < kChamberMaxDEta))
        continue;
      const double bend = t.charge * chamberDeltaPhi(p.phi, t.phi) * t.pt;
      if (!(bend < 0. && bend > -kChamberMaxBend))
        continue;
      double dirBend = 0.;
      if (p.isSegment) {
        dirBend = t.charge * chamberDeltaPhi(p.dirPhi, t.phi) * t.pt;
        if (!(dirBend < 0.))
          continue;
      }
      const ChamberCandidate candidate{
          static_cast<float>(bend), static_cast<float>(dEta), static_cast<float>(dirBend - bend)};
      auto& station = match[p.station];
      if (absDEta < nearest[p.station]) {
        station.second = station.nearest;
        second[p.station] = nearest[p.station];
        station.nearest = candidate;
        nearest[p.station] = absDEta;
      } else if (absDEta < second[p.station]) {
        station.second = candidate;
        second[p.station] = absDEta;
      }
    }
    return match;
  }

}  // namespace muonhp

#endif
