#!/usr/bin/env python3
"""Compare FPP1+FPP2 wire-level reconstruction with truth polar angle."""

import argparse
import math

import ROOT


ZMIN_CM = 111.4215
ANGLE_BINS = ((5.0, 8.0), (8.0, 10.0), (10.0, 12.0),
              (12.0, 15.0), (15.0, 19.0))


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

    correlation = ROOT.TH2D(
        "reco_truth_theta",
        ";truth #theta_{polar} (deg);reconstructed #theta_{polar} (deg)",
        100, 0.0, 20.0, 100, 0.0, 20.0,
    )
    residual = ROOT.TH2D(
        "theta_residual",
        ";truth #theta_{polar} (deg);#theta_{reco}-#theta_{truth} (deg)",
        100, 0.0, 20.0, 160, -8.0, 8.0,
    )
    phi_residual = ROOT.TH2D(
        "phi_residual",
        ";truth #theta_{polar} (deg);wrapped #phi_{reco}-#phi_{truth} (deg)",
        100, 0.0, 20.0, 180, -180.0, 180.0,
    )
    statistics = {bounds: [0.0] * 7 for bounds in ANGLE_BINS}
    selected = reconstructed = 0.0

    for event in tree:
        if not (
            event.front_track_valid
            and event.front_cone_pass
            and ZMIN_CM <= event.front_zclose_cm <= ZMIN_CM + args.thickness
            and event.front_sclose_cm <= args.sclose_max
            and 5.0 <= event.front_track_theta_deg <= 19.0
        ):
            continue
        weight = float(event.weight)
        selected += weight
        if not event.front_reco_valid:
            continue
        reconstructed += weight
        truth = float(event.front_track_theta_deg)
        reco = float(event.front_reco_theta_deg)
        delta = reco - truth
        truth_phi = float(event.front_phi_az_rad)
        reco_phi = float(event.front_reco_phi_az_rad)
        delta_phi = math.degrees(
            math.atan2(math.sin(reco_phi - truth_phi),
                       math.cos(reco_phi - truth_phi))
        )
        correlation.Fill(truth, reco, weight)
        residual.Fill(truth, delta, weight)
        phi_residual.Fill(truth, delta_phi, weight)
        for bounds in ANGLE_BINS:
            low, high = bounds
            if low <= truth < high:
                values = statistics[bounds]
                values[0] += weight
                values[1] += weight * delta
                values[2] += weight * delta * delta
                if abs(delta) > 2.0:
                    values[3] += weight
                values[4] += weight * delta_phi
                values[5] += weight * delta_phi * delta_phi
                if abs(delta_phi) > 30.0:
                    values[6] += weight
                break

    for hist in (correlation, residual):
        hist.SetStats(False)
    canvas = ROOT.TCanvas("closure", "closure", 1500, 1200)
    canvas.Divide(2, 2)
    canvas.cd(1)
    correlation.Draw("colz")
    diagonal = ROOT.TLine(0.0, 0.0, 20.0, 20.0)
    diagonal.SetLineColor(ROOT.kRed + 1)
    diagonal.SetLineWidth(2)
    diagonal.Draw()
    canvas.cd(2)
    residual.Draw("colz")
    zero = ROOT.TLine(0.0, 0.0, 20.0, 0.0)
    zero.SetLineColor(ROOT.kRed + 1)
    zero.SetLineWidth(2)
    zero.Draw()
    canvas.cd(3)
    phi_residual.SetStats(False)
    phi_residual.Draw("colz")
    phi_zero = ROOT.TLine(0.0, 0.0, 20.0, 0.0)
    phi_zero.SetLineColor(ROOT.kRed + 1)
    phi_zero.SetLineWidth(2)
    phi_zero.Draw()
    canvas.cd(4)
    phi_distribution = ROOT.TH1D(
        "phi_residual_projection",
        ";wrapped #phi_{reco}-#phi_{truth} (deg);weighted events",
        180, -180.0, 180.0,
    )
    phi_distribution.SetStats(False)
    phi_distribution.SetLineWidth(2)
    for xbin in range(1, phi_residual.GetNbinsX() + 1):
        for ybin in range(1, phi_residual.GetNbinsY() + 1):
            phi_distribution.AddBinContent(
                ybin, phi_residual.GetBinContent(xbin, ybin)
            )
    phi_distribution.Draw("hist")
    canvas.SaveAs(args.output_stem + ".pdf")
    canvas.SaveAs(args.output_stem + ".png")

    print(f"truth_selected_weight {selected:.9g}")
    print(f"reconstructed_weight {reconstructed:.9g}")
    print(f"reconstruction_fraction {reconstructed / selected:.9g}")
    print(
        "theta_bin_deg weight theta_bias_deg theta_rms_deg "
        "theta_outlier_gt2deg phi_bias_deg phi_rms_deg phi_outlier_gt30deg"
    )
    for bounds, values in statistics.items():
        (weight, total_delta, total_delta2, outliers, total_phi,
         total_phi2, phi_outliers) = values
        bias = total_delta / weight if weight else math.nan
        rms = math.sqrt(total_delta2 / weight) if weight else math.nan
        fraction = outliers / weight if weight else math.nan
        phi_bias = total_phi / weight if weight else math.nan
        phi_rms = math.sqrt(total_phi2 / weight) if weight else math.nan
        phi_fraction = phi_outliers / weight if weight else math.nan
        print(
            f"{bounds[0]:g}-{bounds[1]:g} {weight:.9g} {bias:.7g} "
            f"{rms:.7g} {fraction:.7g} {phi_bias:.7g} {phi_rms:.7g} "
            f"{phi_fraction:.7g}"
        )


if __name__ == "__main__":
    main()
