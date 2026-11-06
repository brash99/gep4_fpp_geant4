#!/usr/bin/env python3
"""Plot reconstructed FPP polar angle versus two-track closest-approach z."""

import argparse

import ROOT


def analyzer_plot(tree, analyzer, canvas, pad):
    canvas.cd(pad)
    if analyzer == 1:
        theta = "front_track_theta_deg"
        zclose = "front_zclose_cm"
        valid = "front_track_valid"
        zmin, zmax = 111.4215, 167.3015
        title = "Analyzer 1: HMS track vs FPP1-FPP2 track"
    else:
        theta = "rear_track_theta_deg"
        zclose = "rear_zclose_cm"
        valid = "rear_track_valid"
        zmin, zmax = 209.754, 265.634
        title = "Analyzer 2: FPP1-FPP2 track vs FPP3-FPP4 track"

    hist = ROOT.TH2D(
        f"theta_z_analyzer{analyzer}",
        f"{title};z_{{close}} (cm);#theta_{{polar}} (deg)",
        160,
        80.0,
        320.0,
        90,
        0.0,
        45.0,
    )
    hist.SetStats(False)
    tree.Draw(f"{theta}:{zclose}>>{hist.GetName()}", f"weight*({valid})", "colz")

    for boundary in (zmin, zmax):
        line = ROOT.TLine(boundary, 0.0, boundary, 45.0)
        line.SetLineColor(ROOT.kRed + 1)
        line.SetLineWidth(2)
        line.SetLineStyle(2)
        line.Draw()
        # Preserve ownership until the canvas is saved.
        ROOT.SetOwnership(line, False)
    return hist


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

    canvas = ROOT.TCanvas("theta_zclose", "theta_zclose", 1500, 650)
    canvas.Divide(2, 1)
    histograms = [
        analyzer_plot(tree, 1, canvas, 1),
        analyzer_plot(tree, 2, canvas, 2),
    ]
    canvas.SaveAs(args.output_stem + ".pdf")
    canvas.SaveAs(args.output_stem + ".png")
    # Keep histograms alive through SaveAs.
    _ = histograms

    quality = ROOT.TCanvas("quality", "quality", 1500, 1000)
    quality.Divide(2, 2)
    quality_hists = []
    for analyzer, prefix, pad_offset in ((1, "front", 0), (2, "rear", 1)):
        quality.cd(1 + pad_offset)
        theta_sclose = ROOT.TH2D(
            f"theta_sclose_{analyzer}",
            f"Analyzer {analyzer};s_{{close}} (cm);#theta_{{polar}} (deg)",
            120,
            0.0,
            30.0,
            90,
            0.0,
            45.0,
        )
        theta_sclose.SetStats(False)
        tree.Draw(
            f"{prefix}_track_theta_deg:{prefix}_sclose_cm>>{theta_sclose.GetName()}",
            f"weight*({prefix}_track_valid)",
            "colz",
        )
        quality_hists.append(theta_sclose)

        quality.cd(3 + pad_offset)
        phi = ROOT.TH1D(
            f"phi_{analyzer}",
            f"Analyzer {analyzer};#phi_{{az}} (rad);weighted events / (#pi/36)",
            72,
            0.0,
            2.0 * ROOT.TMath.Pi(),
        )
        phi.SetStats(False)
        phi.SetLineWidth(2)
        tree.Draw(
            f"{prefix}_phi_az_rad>>{phi.GetName()}",
            f"weight*({prefix}_track_valid && {prefix}_track_theta_deg>=5 && {prefix}_track_theta_deg<=30)",
            "hist",
        )
        quality_hists.append(phi)

    quality.SaveAs(args.output_stem + "_sclose_phi.pdf")
    quality.SaveAs(args.output_stem + "_sclose_phi.png")
    _ = quality_hists


if __name__ == "__main__":
    main()
