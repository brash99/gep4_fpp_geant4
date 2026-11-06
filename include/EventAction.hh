#pragma once

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <array>
#include <unordered_map>

class G4Event;

class EventAction final : public G4UserEventAction {
 public:
  explicit EventAction(bool singleAnalyzer = false)
      : singleAnalyzer_(singleAnalyzer) {}
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void setInput(long id, double weight, double momentumGeV,
                const G4ThreeVector& positionCm,
                const G4ThreeVector& direction);
  void recordHadronicInteraction(int analyzer, int processCategory, double zCm,
                                 double thetaDeg, double outgoingMomentumGeV,
                                 int chargedSecondaries);
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
  bool singleAnalyzer_ = false;
};
