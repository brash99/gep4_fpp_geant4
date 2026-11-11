#pragma once

#include "FPPGeometry.hh"

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <array>
#include <unordered_map>

class G4Event;

class EventAction final : public G4UserEventAction {
 public:
  explicit EventAction(bool singleAnalyzer = false,
                       FPPLayout layout = legacyFPPLayout(),
                       double analyzer1ThicknessCm = 55.88)
      : singleAnalyzer_(singleAnalyzer), layout_(layout),
        planes_(makeHallCFPPPlanes(layout)),
        analyzer1ThicknessCm_(analyzer1ThicknessCm) {}
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void setInput(long id, double weight, double momentumGeV,
                const G4ThreeVector& positionCm,
                const G4ThreeVector& direction);
  void registerTrack(int trackId, int parentId, int creatorProcessCategory);
  void recordHadronicInteraction(int analyzer, int trackId, int parentId,
                                 int processCategory, double zCm,
                                 double incomingMomentumGeV, double thetaDeg,
                                 double outgoingMomentumGeV,
                                 int chargedSecondaries, int particlePdg,
                                 double particleCharge);
  void recordPrimaryAnalyzerExit(const G4ThreeVector& positionCm,
                                 const G4ThreeVector& direction,
                                 double momentumGeV);
  void recordPrimaryKilledInAnalyzer(double momentumGeV);
  void recordChamberCrossing(int chamber, int trackId, int parentId,
                             bool primary, const G4ThreeVector& direction,
                             const G4ThreeVector& positionCm,
                             double momentumGeV);

 private:
  struct TrackObservation {
    int parentId = 0;
    int mask = 0;
    std::array<G4ThreeVector, 4> direction;
    std::array<G4ThreeVector, 4> positionCm;
    std::array<double, 4> momentumGeV{};
  };

  struct TrackMetadata {
    int parentId = 0;
    int rootTrackId = 0;
    int generation = 0;
    int lineageProcessCategory = 0;
    int elasticInteractions = 0;
    int inelasticInteractions = 0;
    int otherHadronicInteractions = 0;
  };

  struct Candidate {
    bool valid = false;
    bool conePass = false;
    double thetaDeg = 0.0;
    double phiRad = 0.0;
    double momentumGeV = 0.0;
    double zCloseCm = 0.0;
    double sCloseCm = 0.0;
    bool reconstructedValid = false;
    double reconstructedThetaDeg = 0.0;
    double reconstructedPhiRad = 0.0;
    double reconstructedReducedChi2 = 0.0;
    int reconstructionPermutations = 0;
  };

  Candidate candidateForTrack(const TrackObservation& track, int analyzer,
                              bool reconstruct) const;
  Candidate bestCandidate(int analyzer) const;
  static bool legacyConeTest(const G4ThreeVector& incomingPointCm,
                             const G4ThreeVector& incomingDirection,
                             double thetaDeg, double rearPlaneZCm);

  long eventId_ = -1;
  double weight_ = 1.0;
  double momentumGeV_ = 0.0;
  G4ThreeVector incidentPositionCm_{0., 0., 0.};
  G4ThreeVector incidentDirection_{0., 0., 1.};
  int analyzer_ = 0;
  int processCategory_ = 0;
  double interactionZCm_ = 0.0;
  double thetaDeg_ = 0.0;
  double outgoingMomentumGeV_ = 0.0;
  int chargedSecondaries_ = 0;
  int primaryMask_ = 0;
  int secondaryMask_ = 0;
  std::unordered_map<int, TrackObservation> tracks_;
  std::unordered_map<int, TrackMetadata> trackMetadata_;
  int interactionCount_ = 0;
  int primaryInteractionCount_ = 0;
  int secondaryInteractionCount_ = 0;
  int primaryElasticCount_ = 0;
  int primaryInelasticCount_ = 0;
  int secondaryElasticCount_ = 0;
  int secondaryInelasticCount_ = 0;
  int primaryExitState_ = 0;  // 1 downstream, 2 upstream, 3 side, 4 killed
  double primaryExitMomentumGeV_ = 0.0;
  double primaryExitThetaDeg_ = 0.0;
  bool singleAnalyzer_ = false;
  FPPLayout layout_;
  FPPPlaneArray planes_;
  double analyzer1ThicknessCm_ = 55.88;
};
