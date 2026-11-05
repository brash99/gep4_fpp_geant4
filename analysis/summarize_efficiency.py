#!/usr/bin/env python3
"""Report weighted geometry-level FPP efficiencies from the transport tree."""

import argparse
import math

import ROOT


def ratio_and_error(selected_weight, total_weight, total_weight2):
    efficiency = selected_weight / total_weight if total_weight else 0.0
    # Effective weighted-count approximation. This is a diagnostic statistical
    # uncertainty, not a transport-model systematic uncertainty.
    effective_total = (
        total_weight * total_weight / total_weight2
        if total_weight2 > 0.0
        else 0.0
    )
    error = (
        math.sqrt(max(0.0, efficiency * (1.0 - efficiency) / effective_total))
        if effective_total > 0.0
        else 0.0
    )
    return efficiency, error


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file")
    parser.add_argument("--theta-min", type=float, default=5.0)
    parser.add_argument("--theta-max", type=float, default=30.0)
    args = parser.parse_args()

    source = ROOT.TFile.Open(args.root_file)
    if not source or source.IsZombie():
        raise SystemExit(f"Cannot open {args.root_file}")
    tree = source.Get("fpp")
    if not tree:
        raise SystemExit("Missing fpp tree")

    total = 0.0
    total2 = 0.0
    selected = {1: 0.0, 2: 0.0}
    interactions = {1: 0.0, 2: 0.0}

    for event in tree:
        weight = float(event.weight)
        total += weight
        total2 += weight * weight
        analyzer = int(event.first_analyzer)
        if analyzer not in (1, 2):
            continue
        interactions[analyzer] += weight
        in_angle = args.theta_min <= event.theta_deg <= args.theta_max
        has_pair = (
            bool(event.same_track_front_pair)
            if analyzer == 1
            else bool(event.same_track_rear_pair)
        )
        if in_angle and has_pair:
            selected[analyzer] += weight

    print(f"input_weight {total:.9g}")
    print(f"theta_range_deg {args.theta_min:g} {args.theta_max:g}")
    for analyzer in (1, 2):
        efficiency, error = ratio_and_error(selected[analyzer], total, total2)
        print(f"analyzer{analyzer}_interaction_weight {interactions[analyzer]:.9g}")
        print(f"analyzer{analyzer}_selected_weight {selected[analyzer]:.9g}")
        print(f"analyzer{analyzer}_efficiency {efficiency:.9g} +/- {error:.3g}")


if __name__ == "__main__":
    main()
