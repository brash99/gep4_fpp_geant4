# FPP thickness-study selection contract

Confirmed by the user on 11 November 2026. This contract controls every
selected-track efficiency and FOM result in the analyzer-thickness study.

## Event denominator

Each thickness sample contains 100,000 unit-weight, monoenergetic 6.7 GeV/c
protons incident on the CH2 analyzer. Every efficiency and FOM is normalized
to this incident-proton denominator, not to a preselected event population.

## Candidate population and multiplicity

- A candidate is a distinct charged GEANT4 track crossing both FPP1 and FPP2.
- Raw multiplicity counts all such front-pair candidates before the angle
  filter, track ranking, or any other analysis cut.
- Zero, single, and multiple refer to zero, exactly one, and at least two
  front-pair candidates, respectively.

## Agreed selected-track algorithm

1. Remove every candidate with scattering angle theta less than or equal to
   1 degree. The lower bound is strictly theta greater than 1 degree.
2. Among the remaining candidates, choose the track with the smallest
   scattering angle. Track ID breaks an exact angular tie deterministically.
3. Do not use simulated outgoing momentum for ranking, filtering, tie-breaking,
   or any other selection decision.
4. Require the chosen track to cross FPP3 and FPP4, giving twelve tracking
   planes after the analyzer.
5. Require the all-four-chamber legacy cone test.
6. Require z-close to lie inside the physical CH2 analyzer interval for the
   thickness being analyzed.
7. Require s-close less than or equal to 1 cm.
8. Do not replace the chosen track if it fails a later requirement.
9. Require 0.07 <= pT <= 1.00 GeV/c, with pT calculated from the known
   incident-proton momentum as pT = p_incident sin(theta). The unmeasured
   outgoing-track momentum is never used.

## FOM bookkeeping

Apply the same event-by-event analyzing-power calculation used for the
smallest-angle-above-one-degree result reported on slide 15. Report separately:

- the useful efficiency and FOM contribution from raw single-candidate events;
- the useful efficiency and FOM contribution from raw multiple-candidate events,
  using the selected-track algorithm above; and
- the combined value as the sum of those two disjoint contributions.

The multi-track analyzing-power result remains provisional because the GEp-III
analyzing-power calibration was established for single-track events. That
physics caveat does not change the selection algorithm or normalization.
