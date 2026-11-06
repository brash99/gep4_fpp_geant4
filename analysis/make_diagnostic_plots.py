#!/usr/bin/env python3
"""Create first-pass interaction and angular diagnostics from FPP transport."""

import argparse

import ROOT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file")
    parser.add_argument("output_stem")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    source = ROOT.TFile.Open(args.root_file)
    if not source or source.IsZombie():
        raise SystemExit(f"Cannot open {args.root_file}")
    tree = source.Get("fpp")
    if not tree:
        raise SystemExit("Missing fpp tree")

    elastic = ROOT.TH1D("elastic", "", 90, 0.0, 45.0)
    inelastic = ROOT.TH1D("inelastic", "", 90, 0.0, 45.0)
    accepted = ROOT.TH1D("accepted", "", 90, 0.0, 45.0)
    all_interactions = ROOT.TH1D("all_interactions", "", 90, 0.0, 45.0)
    for hist in (elastic, inelastic, accepted, all_interactions):
        hist.Sumw2()

    for event in tree:
        if event.first_analyzer == 0:
            continue
        weight = float(event.weight)
        theta = float(event.theta_deg)
        all_interactions.Fill(theta, weight)
        if event.process_category == 1:
            elastic.Fill(theta, weight)
        elif event.process_category == 2:
            inelastic.Fill(theta, weight)
        if event.front_track_valid:
            accepted.Fill(float(event.front_track_theta_deg), weight)
        if event.rear_track_valid:
            accepted.Fill(float(event.rear_track_theta_deg), weight)

    canvas = ROOT.TCanvas("canvas", "canvas", 900, 700)
    all_interactions.SetTitle("6.7 GeV/c proton interactions in GEp-III FPP geometry;scattering angle (deg);weighted events / 0.5 deg")
    all_interactions.SetLineColor(ROOT.kBlack)
    all_interactions.SetLineWidth(2)
    elastic.SetLineColor(ROOT.kBlue + 1)
    elastic.SetLineWidth(2)
    inelastic.SetLineColor(ROOT.kRed + 1)
    inelastic.SetLineWidth(2)
    accepted.SetLineColor(ROOT.kGreen + 2)
    accepted.SetLineWidth(2)
    all_interactions.Draw("hist")
    elastic.Draw("hist same")
    inelastic.Draw("hist same")
    accepted.Draw("hist same")
    legend = ROOT.TLegend(0.57, 0.66, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.AddEntry(all_interactions, "all first hadronic interactions", "l")
    legend.AddEntry(elastic, "hadronic elastic", "l")
    legend.AddEntry(inelastic, "hadronic inelastic", "l")
    legend.AddEntry(accepted, "reconstructed chamber-pair tracks", "l")
    legend.Draw()
    canvas.SaveAs(args.output_stem + ".pdf")
    canvas.SaveAs(args.output_stem + ".png")


if __name__ == "__main__":
    main()
