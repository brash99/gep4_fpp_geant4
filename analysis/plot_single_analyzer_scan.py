#!/usr/bin/env python3
"""Plot the single-analyzer thickness scan and endpoint diagnostics."""

import argparse
import math
import os

import ROOT


POINTS = [(40.0, "40"), (50.0, "50"), (55.88, "55p88"), (60.0, "60")]
ZMIN = 111.4215


def efficiency(tree, thickness):
    total = selected = 0.0
    for event in tree:
        total += event.weight
        if (
            event.single_all4_valid
            and 5.0 <= event.single_all4_theta_deg <= 30.0
            and ZMIN <= event.single_all4_zclose_cm <= ZMIN + thickness
            and event.single_all4_sclose_cm <= 1.0
            and event.single_all4_cone_pass
        ):
            selected += event.weight
    value = selected / total
    error = math.sqrt(value * (1.0 - value) / total)
    return value, error


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", help="Directory containing single_*cm run folders")
    parser.add_argument("output_stem")
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)

    values = []
    sources = {}
    for thickness, tag in POINTS:
        path = os.path.join(args.run_root, f"single_{tag}cm", "transport.root")
        source = ROOT.TFile.Open(path)
        if not source or source.IsZombie():
            raise SystemExit(f"Cannot open {path}")
        sources[tag] = source
        values.append((thickness, *efficiency(source.Get("fpp"), thickness)))

    graph = ROOT.TGraphErrors(len(values))
    for index, (thickness, value, error) in enumerate(values):
        graph.SetPoint(index, thickness, value)
        graph.SetPointError(index, 0.0, error)
    graph.SetTitle(
        "Single-analyzer four-chamber efficiency;CH_{2} thickness (cm);truth efficiency with cone test"
    )
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.3)
    graph.SetLineWidth(2)
    graph.SetMinimum(0.10)
    graph.SetMaximum(0.16)
    canvas = ROOT.TCanvas("thickness", "thickness", 850, 650)
    graph.Draw("APL")
    canvas.SaveAs(args.output_stem + "_efficiency.pdf")
    canvas.SaveAs(args.output_stem + "_efficiency.png")

    comparison = ROOT.TCanvas("comparison", "comparison", 1500, 650)
    comparison.Divide(2, 1)
    histograms = []
    for pad, (thickness, tag) in enumerate((POINTS[0], POINTS[-1]), start=1):
        comparison.cd(pad)
        hist = ROOT.TH2D(
            f"theta_z_{tag}",
            f"{thickness:g} cm CH_{{2}};z_{{close}} (cm);#theta_{{polar}} (deg)",
            120,
            90.0,
            190.0,
            90,
            0.0,
            45.0,
        )
        hist.SetStats(False)
        sources[tag].Get("fpp").Draw(
            f"single_all4_theta_deg:single_all4_zclose_cm>>{hist.GetName()}",
            "weight*single_all4_valid*single_all4_cone_pass",
            "colz",
        )
        for boundary in (ZMIN, ZMIN + thickness):
            line = ROOT.TLine(boundary, 0.0, boundary, 45.0)
            line.SetLineColor(ROOT.kRed + 1)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
            ROOT.SetOwnership(line, False)
        histograms.append(hist)
    comparison.SaveAs(args.output_stem + "_theta_zclose.pdf")
    comparison.SaveAs(args.output_stem + "_theta_zclose.png")
    _ = histograms


if __name__ == "__main__":
    main()
