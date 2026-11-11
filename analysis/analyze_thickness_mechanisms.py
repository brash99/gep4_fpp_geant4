#!/usr/bin/env python3
"""Decompose useful FPP efficiency by track ancestry and hadronic process."""

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyzing_power import PT_MAX_GEV, PT_MIN_GEV, analyzing_power  # noqa: E402
from track_selection import Track, select  # noqa: E402


PROCESS_NAMES = {0: "none", 1: "elastic", 2: "inelastic", 3: "other"}
EXIT_NAMES = {0: "unknown", 1: "downstream", 2: "upstream", 3: "side", 4: "killed"}
CHANNELS = (
    "primary_elastic",
    "secondary_elastic",
    "primary_inelastic",
    "secondary_inelastic",
    "primary_other",
    "secondary_other",
)

MIN_ANGLE_DEG = 1.0
SCLOSE_MAX_CM = 1.0


def select_agreed_track(tracks, thickness):
    """Apply the confirmed slide-15 selection without outgoing momentum.

    Rank front-pair charged candidates strictly above one degree by increasing
    scattering angle.  After choosing, require that same track to cross all
    four chambers and pass the all-four-chamber cone, analyzer z-close, and
    s-close cuts.  A failing chosen track is not replaced.
    """
    return select(
        tracks,
        "min_angle",
        thickness=thickness,
        sclose_max=SCLOSE_MAX_CM,
        theta_switch=math.inf,
        min_angle_deg=MIN_ANGLE_DEG,
    )


def classify_channel(track, first_track_process):
    origin = "primary" if track.generation == 0 else "secondary"
    process = (first_track_process if track.generation == 0
               else track.lineage_process_category)
    return f"{origin}_{PROCESS_NAMES.get(process, 'other')}"


def load_tracks(tree):
    tracks = defaultdict(list)
    for row in tree:
        track = Track(
            int(row.track_id), int(row.chamber_mask), float(row.theta_deg),
            float(row.phi_az_rad), float(row.outgoing_p_GeV),
            float(row.zclose_cm), float(row.sclose_cm),
            bool(row.front_cone_pass), bool(row.all4_cone_pass),
            int(row.root_track_id), int(row.generation),
            int(row.lineage_process_category),
            int(row.track_elastic_interaction_count),
            int(row.track_inelastic_interaction_count),
            int(row.track_other_hadronic_interaction_count),
        )
        tracks[int(row.event_id)].append(track)
    return tracks


def first_processes(tree):
    first = {}
    particle_counts = Counter()
    for row in tree:
        key = (int(row.event_id), int(row.track_id))
        index = int(row.interaction_index)
        if key not in first or index < first[key][0]:
            first[key] = (index, int(row.process_category))
        generation = int(row.generation)
        pdg = int(row.particle_pdg)
        charge = float(row.particle_charge)
        if generation == 0 and pdg == 2212:
            particle_counts["primary_proton"] += 1
        elif generation > 0 and pdg == 2212:
            particle_counts["secondary_proton"] += 1
        elif generation > 0 and charge != 0.0:
            particle_counts["other_secondary_charged"] += 1
        else:
            particle_counts["neutral"] += 1
    return ({key: value[1] for key, value in first.items()}, particle_counts)


def analyze(label, thickness, root_path):
    source = ROOT.TFile.Open(root_path)
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open {root_path}")
    events = source.Get("fpp")
    candidates = source.Get("candidates")
    interactions = source.Get("interactions")
    if not events or not candidates or not interactions:
        raise RuntimeError(f"Missing required ntuple in {root_path}")

    tracks = load_tracks(candidates)
    first_process, particle_interactions = first_processes(interactions)
    totals = Counter()
    channels = Counter()
    multiplicity = Counter()
    useful_multiplicity = Counter()
    exit_states = Counter()
    interaction_bins = Counter()
    useful_interaction_bins = Counter()
    process_interactions = Counter()
    selected_histories = Counter()
    ay2_sums = Counter()
    ay4_weight2_sums = Counter()

    for event in events:
        event_id = int(event.event_id)
        weight = float(event.weight)
        event_tracks = tracks.get(event_id, [])
        count = len(event_tracks)
        mult = "zero" if count == 0 else ("single" if count == 1 else "multi")
        multiplicity[mult] += weight
        totals["incident"] += weight
        exit_states[EXIT_NAMES.get(int(event.primary_analyzer_exit_state), "unknown")] += weight
        n_interactions = int(event.hadronic_interaction_count)
        interaction_bin = "zero" if n_interactions == 0 else (
            "one" if n_interactions == 1 else "multiple")
        interaction_bins[interaction_bin] += weight
        process_interactions["primary_elastic"] += weight * int(event.primary_elastic_interaction_count)
        process_interactions["primary_inelastic"] += weight * int(event.primary_inelastic_interaction_count)
        process_interactions["secondary_elastic"] += weight * int(event.secondary_elastic_interaction_count)
        process_interactions["secondary_inelastic"] += weight * int(event.secondary_inelastic_interaction_count)

        chosen = select_agreed_track(event_tracks, thickness)
        if chosen is None:
            continue
        totals["accepted"] += weight
        ay, pt = analyzing_power(float(event.incident_p_GeV), chosen.theta_deg)
        if not (PT_MIN_GEV <= pt <= PT_MAX_GEV):
            continue
        totals["useful"] += weight
        useful_multiplicity[mult] += weight
        ay2_sums[mult] += weight * ay * ay
        ay4_weight2_sums[mult] += weight * weight * ay ** 4
        useful_interaction_bins[interaction_bin] += weight
        process = first_process.get((event_id, chosen.track_id), 0)
        channel = classify_channel(chosen, process)
        channels[channel if channel in CHANNELS else channel.split("_")[0] + "_other"] += weight
        if chosen.elastic_interactions and chosen.inelastic_interactions:
            selected_histories["mixed_elastic_inelastic"] += weight
        elif chosen.elastic_interactions:
            selected_histories["elastic_on_selected_track"] += weight
        elif chosen.inelastic_interactions:
            selected_histories["inelastic_on_selected_track"] += weight
        else:
            selected_histories["no_hadronic_on_selected_track"] += weight

    incident = totals["incident"]
    total_ay2 = ay2_sums["single"] + ay2_sums["multi"]
    fom_by_multiplicity = {}
    for key in ("single", "multi"):
        fom = ay2_sums[key] / incident
        fom_by_multiplicity[key] = {
            "useful_efficiency": useful_multiplicity[key] / incident,
            "intrinsic_fom": fom,
            "intrinsic_fom_mc_se": math.sqrt(max(
                0.0,
                ay4_weight2_sums[key] - 2.0 * fom * ay2_sums[key]
                + fom * fom * incident,
            )) / incident,
            "effective_Ay": math.sqrt(
                ay2_sums[key] / useful_multiplicity[key]
            ) if useful_multiplicity[key] else 0.0,
        }
    combined_fom = total_ay2 / incident
    fom_by_multiplicity["combined"] = {
        "useful_efficiency": totals["useful"] / incident,
        "intrinsic_fom": combined_fom,
        "intrinsic_fom_mc_se": math.sqrt(max(
            0.0,
            ay4_weight2_sums["single"] + ay4_weight2_sums["multi"]
            - 2.0 * combined_fom * total_ay2
            + combined_fom * combined_fom * incident,
        )) / incident,
        "effective_Ay": math.sqrt(total_ay2 / totals["useful"])
        if totals["useful"] else 0.0,
    }

    result = {
        "label": label,
        "thickness_cm": thickness,
        "root_file": root_path,
        "incident_weight": incident,
        "accepted_efficiency": totals["accepted"] / incident,
        "useful_efficiency": totals["useful"] / incident,
        "fom_by_multiplicity": fom_by_multiplicity,
        "multiplicity_fraction": {key: multiplicity[key] / incident
                                  for key in ("zero", "single", "multi")},
        "useful_multiplicity_efficiency": {
            key: useful_multiplicity[key] / incident
            for key in ("single", "multi")
        },
        "channel_efficiency": {key: channels[key] / incident for key in CHANNELS},
        "primary_exit_fraction": {key: exit_states[key] / incident
                                  for key in EXIT_NAMES.values()},
        "event_interaction_fraction": {key: interaction_bins[key] / incident
                                       for key in ("zero", "one", "multiple")},
        "useful_interaction_efficiency": {
            key: useful_interaction_bins[key] / incident
            for key in ("zero", "one", "multiple")
        },
        "mean_interactions_per_incident": {
            key: process_interactions[key] / incident for key in
            ("primary_elastic", "primary_inelastic", "secondary_elastic", "secondary_inelastic")
        },
        "mean_interactions_by_particle_per_incident": {
            key: particle_interactions[key] / incident for key in
            ("primary_proton", "secondary_proton",
             "other_secondary_charged", "neutral")
        },
        "useful_selected_track_history_efficiency": {
            key: selected_histories[key] / incident for key in
            ("no_hadronic_on_selected_track", "elastic_on_selected_track",
             "inelastic_on_selected_track", "mixed_elastic_inelastic")
        },
        "selection_contract": {
            "candidate_population": "charged tracks crossing FPP1 and FPP2",
            "strict_minimum_angle_deg": MIN_ANGLE_DEG,
            "ranking": "smallest scattering angle; track ID breaks exact ties",
            "outgoing_momentum_used": False,
            "chosen_track_must_cross_all_four_chambers": True,
            "replacement_after_failed_cut": False,
            "cone": "all-four-chamber legacy cone test",
            "zclose": "inside the physical CH2 analyzer interval",
            "sclose_max_cm": SCLOSE_MAX_CM,
            "pt_GeV": [PT_MIN_GEV, PT_MAX_GEV],
            "pt_momentum": "known incident-proton momentum",
            "multiplicity": "all front-pair candidates before angle or analysis cuts",
        },
    }
    source.Close()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", nargs=3, action="append", required=True,
                        metavar=("LABEL", "THICKNESS_CM", "ROOT_FILE"))
    parser.add_argument("--output-stem", required=True)
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)

    results = [analyze(label, float(thickness), path)
               for label, thickness, path in args.sample]
    output = Path(args.output_stem)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output) + ".json", "w") as stream:
        json.dump({"samples": results}, stream, indent=2, sort_keys=True)

    fields = [
        "label", "thickness_cm", "accepted_efficiency", "useful_efficiency",
        "zero_track_fraction", "single_track_fraction", "multi_track_fraction",
        *[f"{channel}_efficiency" for channel in CHANNELS],
        "primary_killed_fraction", "multiple_interaction_fraction",
    ]
    with open(str(output) + ".csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for result in results:
            row = {
                "label": result["label"],
                "thickness_cm": result["thickness_cm"],
                "accepted_efficiency": result["accepted_efficiency"],
                "useful_efficiency": result["useful_efficiency"],
                "zero_track_fraction": result["multiplicity_fraction"]["zero"],
                "single_track_fraction": result["multiplicity_fraction"]["single"],
                "multi_track_fraction": result["multiplicity_fraction"]["multi"],
                "primary_killed_fraction": result["primary_exit_fraction"]["killed"],
                "multiple_interaction_fraction": result["event_interaction_fraction"]["multiple"],
            }
            row.update({f"{key}_efficiency": value
                        for key, value in result["channel_efficiency"].items()})
            writer.writerow(row)
    print(json.dumps({"samples": results}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
