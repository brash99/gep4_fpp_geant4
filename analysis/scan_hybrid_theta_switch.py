#!/usr/bin/env python3
"""Scan a hybrid four-/two-chamber truth selection versus theta switch."""

import argparse

import ROOT


ZMIN_CM = 111.4215


def base(event, prefix, thickness, sclose_max):
    if prefix == "front":
        return (
            event.front_track_valid
            and event.front_cone_pass
            and ZMIN_CM <= event.front_zclose_cm <= ZMIN_CM + thickness
            and event.front_sclose_cm <= sclose_max
            and event.front_track_theta_deg >= 5.0
        )
    return (
        event.single_all4_valid
        and event.single_all4_cone_pass
        and ZMIN_CM <= event.single_all4_zclose_cm <= ZMIN_CM + thickness
        and event.single_all4_sclose_cm <= sclose_max
        and event.single_all4_theta_deg >= 5.0
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file")
    parser.add_argument("output_stem")
    parser.add_argument("--thickness", type=float, required=True)
    parser.add_argument("--sclose-max", type=float, default=1.0)
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)

    source = ROOT.TFile.Open(args.root_file)
    tree = source.Get("fpp") if source else None
    if not tree:
        raise SystemExit(f"Cannot open fpp tree in {args.root_file}")

    switches = [5.0 + 0.5 * index for index in range(28)]  # 5.0--18.5 deg
    selected = [0.0] * len(switches)
    total = 0.0
    four_total = 0.0
    for event in tree:
        weight = float(event.weight)
        total += weight
        two = base(event, "front", args.thickness, args.sclose_max)
        four = base(event, "all4", args.thickness, args.sclose_max)
        if four:
            four_total += weight
        theta = float(event.front_track_theta_deg) if two else -1.0
        for index, switch in enumerate(switches):
            # Below the switch require all four chambers; at and above it use
            # the first chamber pair and its shorter-lever-arm cone test.
            if (four and event.single_all4_theta_deg < switch) or (
                two and theta >= switch
            ):
                selected[index] += weight

    graph = ROOT.TGraph(len(switches))
    for index, switch in enumerate(switches):
        graph.SetPoint(index, switch, selected[index] / total)
    graph.SetTitle(
        "Hybrid FPP truth efficiency;#theta_{switch} (deg);accepted efficiency"
    )
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.1)
    graph.SetLineWidth(2)
    graph.SetMinimum(0.13)
    graph.SetMaximum(0.19)
    canvas = ROOT.TCanvas("hybrid_switch", "hybrid_switch", 850, 650)
    graph.Draw("APL")
    canvas.SaveAs(args.output_stem + ".pdf")
    canvas.SaveAs(args.output_stem + ".png")

    print(f"input_weight {total:.9g}")
    best_index = max(range(len(selected)), key=selected.__getitem__)
    print(
        f"maximum theta_switch_deg {switches[best_index]:.1f} "
        f"efficiency {selected[best_index] / total:.7f}"
    )
    print("theta_switch_deg efficiency relative_to_four")
    four_reference = four_total / total
    for switch, count in zip(switches, selected):
        if switch in (5.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 18.5):
            efficiency = count / total
            print(f"{switch:16.1f} {efficiency:10.7f} {efficiency/four_reference:16.5f}")


if __name__ == "__main__":
    main()
