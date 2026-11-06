#!/usr/bin/env python3
"""Scan geometry-level FPP efficiency versus closest-approach cuts."""

import argparse
import math

import ROOT


ANALYZERS = {
    1: ("front", 111.4215, 167.3015),
    2: ("rear", 209.7540, 265.6340),
}


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

    cuts = [1.0, 2.0, 3.0, 5.0, 10.0, math.inf]
    selected = {cut: {1: 0.0, 2: 0.0, "union": 0.0} for cut in cuts}
    total = 0.0

    for event in tree:
        weight = float(event.weight)
        total += weight
        for cut in cuts:
            accepted = {}
            for analyzer, (prefix, zmin, zmax) in ANALYZERS.items():
                valid = bool(getattr(event, f"{prefix}_track_valid"))
                theta = float(getattr(event, f"{prefix}_track_theta_deg"))
                zclose = float(getattr(event, f"{prefix}_zclose_cm"))
                sclose = float(getattr(event, f"{prefix}_sclose_cm"))
                accepted[analyzer] = (
                    valid
                    and args.theta_min <= theta <= args.theta_max
                    and zmin <= zclose <= zmax
                    and sclose <= cut
                )
                if accepted[analyzer]:
                    selected[cut][analyzer] += weight
            if accepted[1] or accepted[2]:
                selected[cut]["union"] += weight

    print(f"input_weight {total:.9g}")
    print(f"theta_range_deg {args.theta_min:g} {args.theta_max:g}")
    print("sclose_max_cm analyzer1_eff analyzer2_eff union_eff")
    for cut in cuts:
        label = "inf" if math.isinf(cut) else f"{cut:g}"
        print(
            f"{label:>13} "
            f"{selected[cut][1] / total:13.7f} "
            f"{selected[cut][2] / total:13.7f} "
            f"{selected[cut]['union'] / total:9.7f}"
        )


if __name__ == "__main__":
    main()
