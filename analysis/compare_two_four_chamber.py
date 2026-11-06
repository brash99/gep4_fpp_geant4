#!/usr/bin/env python3
"""Compare cone-tested truth efficiency using two or four rear chambers."""

import argparse
import math

import ROOT


ZMIN_CM = 111.4215


def accepted(event, prefix, thickness, theta_max, sclose_max):
    if prefix == "front":
        valid = event.front_track_valid and event.front_cone_pass
        theta = float(event.front_track_theta_deg)
        zclose = float(event.front_zclose_cm)
        sclose = float(event.front_sclose_cm)
    else:
        valid = event.single_all4_valid and event.single_all4_cone_pass
        theta = float(event.single_all4_theta_deg)
        zclose = float(event.single_all4_zclose_cm)
        sclose = float(event.single_all4_sclose_cm)
    return (
        valid
        and 5.0 <= theta <= theta_max
        and ZMIN_CM <= zclose <= ZMIN_CM + thickness
        and sclose <= sclose_max
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file")
    parser.add_argument("output_stem")
    parser.add_argument("--thickness", type=float, required=True)
    parser.add_argument("--sclose-max", type=float, default=1.0)
    parser.add_argument("--theta-plot-max", type=float, default=45.0)
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    source = ROOT.TFile.Open(args.root_file)
    tree = source.Get("fpp") if source else None
    if not tree:
        raise SystemExit(f"Cannot open fpp tree in {args.root_file}")

    bins = int(args.theta_plot_max)
    two_hist = ROOT.TH1D(
        "two_chamber_theta",
        ";#theta_{polar} (deg);accepted weighted events / input event / degree",
        bins,
        0.0,
        args.theta_plot_max,
    )
    four_hist = ROOT.TH1D(
        "four_chamber_theta",
        ";#theta_{polar} (deg);accepted weighted events / input event / degree",
        bins,
        0.0,
        args.theta_plot_max,
    )
    total = 0.0
    theta_maxima = list(range(6, int(args.theta_plot_max) + 1))
    two_cumulative = [0.0] * len(theta_maxima)
    four_cumulative = [0.0] * len(theta_maxima)
    accepted_theta = {"two": [], "four": []}

    for event in tree:
        weight = float(event.weight)
        total += weight
        two_base = accepted(
            event, "front", args.thickness, args.theta_plot_max, args.sclose_max
        )
        four_base = accepted(
            event, "all4", args.thickness, args.theta_plot_max, args.sclose_max
        )
        if two_base:
            theta = float(event.front_track_theta_deg)
            accepted_theta["two"].append(theta)
            two_hist.Fill(theta, weight)
            for index, theta_max in enumerate(theta_maxima):
                if theta <= theta_max:
                    two_cumulative[index] += weight
        if four_base:
            theta = float(event.single_all4_theta_deg)
            accepted_theta["four"].append(theta)
            four_hist.Fill(theta, weight)
            for index, theta_max in enumerate(theta_maxima):
                if theta <= theta_max:
                    four_cumulative[index] += weight

    for hist in (two_hist, four_hist):
        hist.Scale(1.0 / total, "width")
        hist.SetStats(False)
        hist.SetLineWidth(3)
    two_hist.SetLineColor(ROOT.kBlue + 1)
    four_hist.SetLineColor(ROOT.kRed + 1)

    two_graph = ROOT.TGraph(len(theta_maxima))
    four_graph = ROOT.TGraph(len(theta_maxima))
    for index, theta_max in enumerate(theta_maxima):
        two_graph.SetPoint(index, theta_max, two_cumulative[index] / total)
        four_graph.SetPoint(index, theta_max, four_cumulative[index] / total)
    for graph, color in ((two_graph, ROOT.kBlue + 1), (four_graph, ROOT.kRed + 1)):
        graph.SetLineColor(color)
        graph.SetLineWidth(3)

    canvas = ROOT.TCanvas("two_four", "two_four", 1500, 650)
    canvas.Divide(2, 1)
    canvas.cd(1)
    maximum = max(two_hist.GetMaximum(), four_hist.GetMaximum()) * 1.18
    two_hist.SetMaximum(maximum)
    two_hist.Draw("hist")
    four_hist.Draw("hist same")
    legend = ROOT.TLegend(0.52, 0.72, 0.88, 0.88)
    legend.AddEntry(two_hist, "FPP1+FPP2", "l")
    legend.AddEntry(four_hist, "FPP1--FPP4", "l")
    legend.Draw()

    canvas.cd(2)
    frame = ROOT.TH1D(
        "cumulative_frame",
        ";upper #theta_{polar} cut (deg);cumulative truth efficiency",
        100,
        5.0,
        args.theta_plot_max,
    )
    frame.SetStats(False)
    frame.SetMinimum(0.0)
    frame.SetMaximum(max(two_cumulative[-1], four_cumulative[-1]) / total * 1.15)
    frame.Draw()
    two_graph.Draw("L same")
    four_graph.Draw("L same")
    legend2 = ROOT.TLegend(0.52, 0.72, 0.88, 0.88)
    legend2.AddEntry(two_graph, "FPP1+FPP2", "l")
    legend2.AddEntry(four_graph, "FPP1--FPP4", "l")
    legend2.Draw()

    canvas.SaveAs(args.output_stem + ".pdf")
    canvas.SaveAs(args.output_stem + ".png")

    print(f"input_weight {total:.9g}")
    for label in ("two", "four"):
        values = accepted_theta[label]
        above20 = sum(value > 20.0 for value in values)
        maximum = max(values) if values else math.nan
        print(
            f"{label}_accepted_entries {len(values)} "
            f"above_20_deg {above20} max_theta_deg {maximum:.7f}"
        )
    for theta_max in (20, 25, 30, 35, 40, 45):
        if theta_max not in theta_maxima:
            continue
        index = theta_maxima.index(theta_max)
        two = two_cumulative[index] / total
        four = four_cumulative[index] / total
        ratio = two / four if four > 0.0 else math.inf
        print(f"theta_5_{theta_max}_deg two {two:.7f} four {four:.7f} ratio {ratio:.4f}")


if __name__ == "__main__":
    main()
