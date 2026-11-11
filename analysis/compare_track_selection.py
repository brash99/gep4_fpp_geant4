#!/usr/bin/env python3
"""Compare baseline, exclusive single-track, and smallest-angle ideal FOM."""
import argparse
import csv
import json
import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path

import ROOT
from analyzing_power import analyzing_power, PT_MIN_GEV, PT_MAX_GEV
from calculate_fom import accepted_observable
from track_selection import Track, select


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root_file')
    parser.add_argument('output_stem')
    parser.add_argument('--min-angle-deg', type=float, default=1.0,
                        help='Strict lower bound applied before smallest-angle ranking only')
    parser.add_argument('--thickness', type=float, default=55.88)
    parser.add_argument('--beam-polarization', type=float, default=.85)
    parser.add_argument('--baseline-root', help='Verify every original event branch against old transport')
    args = parser.parse_args()
    source = ROOT.TFile.Open(args.root_file)
    if not source or source.IsZombie() or source.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Input ROOT file is missing, corrupt, or recovered')
    events, candidates = source.Get('fpp'), source.Get('candidates')
    if not events or not candidates:
        raise RuntimeError('Both fpp and candidates trees are required; rerun transport')
    old_file = ROOT.TFile.Open(args.baseline_root) if args.baseline_root else None
    old = old_file.Get('fpp') if old_file else None
    if args.baseline_root and (not old or old.GetEntries() != events.GetEntries()):
        raise AssertionError('Baseline entry mismatch')
    original_names = [b.GetName() for b in old.GetListOfBranches()] if old else []
    # The published sample predates the input-metadata persistence fix. Its
    # event_id/p fields are defaults; compare physics by entry, not those IDs.
    metadata_defaults = False
    if old and old.GetEntries():
        old.GetEntry(0)
        metadata_defaults = int(old.event_id) == -1 and float(old.incident_p_GeV) == 0.
    excluded = ['event_id', 'weight', 'incident_p_GeV'] if metadata_defaults else []
    original_names = [n for n in original_names if n not in excluded]
    tracks = defaultdict(list)
    seen = set()
    for c in candidates:
        key = (int(c.event_id), int(c.track_id))
        if key in seen:
            raise AssertionError(f'Duplicate candidate: {key}')
        seen.add(key)
        if (int(c.chamber_mask) & 3) != 3:
            raise AssertionError('Candidate does not cross the front pair')
        ancestry = (
            int(c.root_track_id), int(c.generation),
            int(c.lineage_process_category),
            int(c.track_elastic_interaction_count),
            int(c.track_inelastic_interaction_count),
            int(c.track_other_hadronic_interaction_count),
        ) if hasattr(c, "generation") else (0, 0, 0, 0, 0, 0)
        tracks[key[0]].append(Track(key[1], int(c.chamber_mask), float(c.theta_deg),
            float(c.phi_az_rad), float(c.outgoing_p_GeV), float(c.zclose_cm),
            float(c.sclose_cm), bool(c.front_cone_pass), bool(c.all4_cone_pass),
            *ancestry))
    policies = ('legacy', 'single', 'min_angle_uncut', 'min_angle')
    totals = {p: dict(accepted_events=0, useful_events=0, accepted_weight=0.,
                     useful_weight=0., ay2_weight=0., multi_useful_events=0,
                     multi_ay2_weight=0.) for p in policies}
    total_weight = 0.
    total_weight2 = 0.
    moments = {p: [0., 0., 0.] for p in policies}
    multiplicities = Counter()
    audit = []
    ids = set()
    overlap = Counter()
    single_above = dict(useful_events=0, ay2_weight=0.)
    for index, e in enumerate(events):
        eid, weight = int(e.event_id), float(e.weight)
        if eid in ids:
            raise AssertionError('Duplicate event ID')
        ids.add(eid)
        if old:
            old.GetEntry(index)
            if metadata_defaults and not (int(old.event_id) == -1 and
                    float(old.incident_p_GeV) == 0. and float(old.weight) == 1.):
                raise AssertionError('Legacy metadata are not uniformly the documented defaults')
            for name in original_names:
                if getattr(old, name) != getattr(e, name):
                    raise AssertionError(f'Original branch changed: {index} {name}')
        event_tracks = tracks.pop(eid, [])
        multiplicities[len(event_tracks)] += 1
        total_weight += weight
        total_weight2 += weight*weight
        chosen = {p: select(event_tracks, p, args.thickness, min_angle_deg=args.min_angle_deg) for p in policies}
        baseline = accepted_observable(e, args.thickness, 1., 12., 0.)
        legacy = chosen['legacy']
        if (baseline is None) != (legacy is None) or (legacy and not math.isclose(
                baseline[0], legacy.theta_deg, abs_tol=1e-10)):
            raise AssertionError(f'Candidate-table baseline mismatch: event {eid}')
        useful = {}
        row = dict(event_id=eid, weight=weight, front_track_multiplicity=len(event_tracks))
        for p, track in chosen.items():
            ay, pt = analyzing_power(float(e.incident_p_GeV), track.theta_deg) if track else (0., 0.)
            useful[p] = bool(track and PT_MIN_GEV <= pt <= PT_MAX_GEV)
            moments[p][0] += weight*weight*int(useful[p])
            moments[p][1] += weight*weight*ay*ay
            moments[p][2] += weight*weight*ay**4
            row.update({p+'_track_id': track.track_id if track else -1,
                        p+'_theta_deg': track.theta_deg if track else '',
                        p+'_useful': int(useful[p]), p+'_Ay2': ay*ay})
            if not track:
                continue
            t = totals[p]
            t['accepted_events'] += 1
            t['accepted_weight'] += weight
            t['ay2_weight'] += weight * ay*ay
            if useful[p]:
                t['useful_events'] += 1
                t['useful_weight'] += weight
                if len(event_tracks) > 1:
                    t['multi_useful_events'] += 1
                    t['multi_ay2_weight'] += weight*ay*ay
        for track in (chosen['min_angle'],):
            if track and not track.theta_deg > args.min_angle_deg:
                raise AssertionError('Smallest-angle track fails strict angular threshold')
        single_eligible = bool(useful['single'] and chosen['single'].theta_deg > args.min_angle_deg)
        if single_eligible:
            single_above['useful_events'] += 1
            single_above['ay2_weight'] += weight * row['single_Ay2']
        if single_eligible and not useful['min_angle']:
            raise AssertionError('Single-track useful sample must be contained in min-angle sample')
        overlap[f"legacy{int(useful['legacy'])}_min_angle{int(useful['min_angle'])}"] += 1
        audit.append(row)
    if tracks:
        raise AssertionError('Orphan candidate event IDs')
    if totals['min_angle']['useful_events'] != (single_above['useful_events'] +
            totals['min_angle']['multi_useful_events']):
        raise AssertionError('Recovered useful-event count does not close')
    if not math.isclose(totals['min_angle']['ay2_weight'],
            single_above['ay2_weight'] + totals['min_angle']['multi_ay2_weight'],
            rel_tol=1e-12, abs_tol=1e-12):
        raise AssertionError('Recovered multi-track FOM does not close')
    for p, t in totals.items():
        t.update(useful_efficiency=t['useful_weight']/total_weight,
                 intrinsic_fom=t['ay2_weight']/total_weight,
                 polarized_fom=args.beam_polarization**2*t['ay2_weight']/total_weight,
                 effective_Ay=math.sqrt(t['ay2_weight']/t['useful_weight']) if t['useful_weight'] else 0.)
        efficiency = t['useful_efficiency']
        fom = t['intrinsic_fom']
        m = moments[p]
        t['useful_efficiency_mc_se'] = math.sqrt(max(0.,
            m[0]*(1.-2.*efficiency)+efficiency**2*total_weight2))/total_weight
        t['intrinsic_fom_mc_se'] = math.sqrt(max(0.,
            m[2]-2.*fom*m[1]+fom*fom*total_weight2))/total_weight
    source_paths = ['src/EventAction.cc', 'src/RunAction.cc',
                    'analysis/track_selection.py', 'analysis/compare_track_selection.py',
                    'analysis/analyzing_power.py']
    repo = Path(__file__).resolve().parents[1]
    hashes = {p: hashlib.sha256((repo/p).read_bytes()).hexdigest() for p in source_paths}
    result = dict(min_angle_strict_lower_bound_deg=args.min_angle_deg,
                  single_above_threshold=single_above, source_sha256=hashes, input_file=str(Path(args.root_file).resolve()), events=len(ids),
                  incident_weight=total_weight, candidate_rows=len(seen),
                  front_multiplicity=dict(sorted(multiplicities.items())), policies=totals,
                  useful_overlap=dict(overlap), baseline_physics_branches_exact=bool(old),
                  baseline_excluded_metadata=excluded,
                  candidate_baseline_exact=True, beam_polarization=args.beam_polarization,
                  convention='Ideal truth-associated tracks; count all front-pair tracks before cuts; smallest angle strictly above configured threshold, then remaining cuts; multi-track Ay unvalidated')
    Path(args.output_stem+'.json').write_text(json.dumps(result, indent=2)+'\n')
    with open(args.output_stem+'.csv', 'w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=['policy']+list(totals['legacy']))
        writer.writeheader()
        for p,t in totals.items(): writer.writerow(dict(policy=p, **t))
    with open(args.output_stem+'_events.csv', 'w', newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(audit[0]))
        writer.writeheader(); writer.writerows(audit)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
