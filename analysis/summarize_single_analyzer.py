#!/usr/bin/env python3
"""Summarize single-analyzer, four-chamber truth-level efficiency."""

import argparse
import ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file")
    parser.add_argument("--thickness", type=float, required=True)
    parser.add_argument("--theta-min", type=float, default=5.0)
    parser.add_argument("--theta-max", type=float, default=30.0)
    parser.add_argument("--sclose-max", type=float, default=1.0)
    parser.add_argument(
        "--no-cone-test",
        action="store_true",
        help="report the legacy efficiency without requiring full azimuthal acceptance",
    )
    args = parser.parse_args()
    zmin = 111.4215
    zmax = zmin + args.thickness

    source = ROOT.TFile.Open(args.root_file)
    tree = source.Get("fpp") if source else None
    if not tree:
        raise SystemExit("Cannot open fpp tree")
    total = selected = 0.0
    for event in tree:
        weight = float(event.weight)
        total += weight
        if (
            event.single_all4_valid
            and args.theta_min <= event.single_all4_theta_deg <= args.theta_max
            and zmin <= event.single_all4_zclose_cm <= zmax
            and event.single_all4_sclose_cm <= args.sclose_max
            and (args.no_cone_test or event.single_all4_cone_pass)
        ):
            selected += weight
    print(f"thickness_cm {args.thickness:g}")
    print(f"input_weight {total:.9g}")
    print(f"selected_weight {selected:.9g}")
    print(f"cone_test_required {0 if args.no_cone_test else 1}")
    print(f"efficiency {selected / total:.9g}")


if __name__ == "__main__":
    main()
