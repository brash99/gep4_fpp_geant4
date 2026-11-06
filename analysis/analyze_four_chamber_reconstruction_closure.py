#!/usr/bin/env python3
"""Compare twelve-plane FPP reconstruction with truth below 10 degrees."""

import argparse
import math

import ROOT


ZMIN_CM = 111.4215
ANGLE_BINS = (
    (0.0, 0.5),
    (0.5, 1.0),
    (1.0, 2.0),
    (2.0, 3.0),
    (3.0, 4.0),
    (4.0, 5.0),
    (5.0, 7.0),
    (7.0, 10.0),
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

    correlation = ROOT.TH2D(
        "all4_reco_truth_theta",
        ";truth #theta_{polar} (deg);reconstructed #theta_{polar} (deg)",
        100, 0.0, 10.0, 100, 0.0, 10.0,
    )
    theta_residual = ROOT.TH2D(
        "all4_theta_residual",
        ";truth #theta_{polar} (deg);#theta_{reco}-#theta_{truth} (deg)",
        100, 0.0, 10.0, 200, -2.0, 2.0,
    )
    phi_residual = ROOT.TH2D(
        "all4_phi_residual",
        ";truth #theta_{polar} (deg);wrapped #phi_{reco}-#phi_{truth} (deg)",
        100, 0.0, 10.0, 180, -180.0, 180.0,
    )
    statistics = {bounds: [0.0] * 7 for bounds in ANGLE_BINS}
    selected = reconstructed = 0.0

    for event in tree:
        if not (
            event.single_all4_valid
            and event.single_all4_cone_pass
            and ZMIN_CM <= event.single_all4_zclose_cm <= ZMIN_CM + args.thickness
            and event.single_all4_sclose_cm <= args.sclose_max
            and 0.0 <= event.single_all4_theta_deg < 10.0
        ):
            continue
        weight = float(event.weight)
        selected += weight
        if not event.all4_reco_valid:
            continue
        reconstructed += weight
        truth = float(event.single_all4_theta_deg)
        reco = float(event.all4_reco_theta_deg)
        delta = reco - truth
        delta_phi = math.degrees(
            math.atan2(
                math.sin(event.all4_reco_phi_az_rad - event.single_all4_phi_az_rad),
                math.cos(event.all4_reco_phi_az_rad - event.single_all4_phi_az_rad),
            )
        )
        correlation.Fill(truth, reco, weight)
        theta_residual.Fill(truth, delta, weight)
        phi_residual.Fill(truth, delta_phi, weight)
        for bounds in ANGLE_BINS:
            if bounds[0] <= truth < bounds[1]:
                values = statistics[bounds]
                values[0] += weight
                values[1] += weight * delta
                values[2] += weight * delta * delta
                values[3] += weight if abs(delta) > 0.5 else 0.0
                values[4] += weight * delta_phi
                values[5] += weight * delta_phi * delta_phi
                values[6] += weight if abs(delta_phi) > 30.0 else 0.0
                break

    canvas = ROOT.TCanvas("all4_closure", "all4_closure", 1500, 1200)
    canvas.Divide(2, 2)
    for pad, histogram in enumerate((correlation, theta_residual, phi_residual), 1):
        canvas.cd(pad)
        histogram.SetStats(False)
        histogram.Draw("colz")
    canvas.cd(4)
    projection = phi_residual.ProjectionY("all4_phi_residual_projection")
    projection.SetTitle(";wrapped #phi_{reco}-#phi_{truth} (deg);weighted events")
    projection.SetStats(False)
    projection.SetLineWidth(2)
    projection.Draw("hist")
    canvas.SaveAs(args.output_stem + ".pdf")
    canvas.SaveAs(args.output_stem + ".png")

    print(f"truth_selected_weight {selected:.9g}")
    print(f"reconstructed_weight {reconstructed:.9g}")
    print(f"reconstruction_fraction {reconstructed / selected:.9g}")
    print(
        "theta_bin_deg weight theta_bias_deg theta_rms_deg "
        "theta_outlier_gt0p5deg phi_bias_deg phi_rms_deg phi_outlier_gt30deg"
    )
    for bounds, values in statistics.items():
        weight, sum_d, sum_d2, out_d, sum_p, sum_p2, out_p = values
        if not weight:
            print(f"{bounds[0]:g}-{bounds[1]:g} 0 nan nan nan nan nan nan")
            continue
        print(
            f"{bounds[0]:g}-{bounds[1]:g} {weight:.9g} {sum_d/weight:.7g} "
            f"{math.sqrt(sum_d2/weight):.7g} {out_d/weight:.7g} "
            f"{sum_p/weight:.7g} {math.sqrt(sum_p2/weight):.7g} {out_p/weight:.7g}"
        )


if __name__ == "__main__":
    main()
