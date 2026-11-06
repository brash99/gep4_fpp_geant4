#!/usr/bin/env python3
"""Calculate truth-level FPP efficiency and analyzing-power figure of merit."""

import argparse
import math
import sys
from pathlib import Path

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyzing_power import (  # noqa: E402
    AY_B,
    AY_C0,
    AY_C1,
    PT_MAX_GEV,
    PT_MIN_GEV,
    analyzing_power,
)

ZMIN_CM = 111.4215


def accepted_observable(
    event, thickness_cm, sclose_max_cm, theta_switch_deg, theta_min_deg
):
    """Return (theta, momentum) for the adopted hybrid topology, else None."""
    zmax_cm = ZMIN_CM + thickness_cm
    four = (
        event.single_all4_valid
        and event.single_all4_cone_pass
        and ZMIN_CM <= event.single_all4_zclose_cm <= zmax_cm
        and event.single_all4_sclose_cm <= sclose_max_cm
        and event.single_all4_theta_deg >= theta_min_deg
    )
    two = (
        event.front_track_valid
        and event.front_cone_pass
        and ZMIN_CM <= event.front_zclose_cm <= zmax_cm
        and event.front_sclose_cm <= sclose_max_cm
        and event.front_track_theta_deg >= theta_min_deg
    )
    if four and event.single_all4_theta_deg < theta_switch_deg:
        return float(event.single_all4_theta_deg), float(event.incident_p_GeV)
    if two and event.front_track_theta_deg >= theta_switch_deg:
        return float(event.front_track_theta_deg), float(event.incident_p_GeV)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file")
    parser.add_argument("output_stem")
    parser.add_argument("--thickness", type=float, required=True)
    parser.add_argument("--theta-switch", type=float, default=12.0)
    parser.add_argument("--theta-min", type=float, default=0.0)
    parser.add_argument("--sclose-max", type=float, default=1.0)
    parser.add_argument("--beam-polarization", type=float, default=0.85)
    parser.add_argument(
        "--momentum-gev",
        type=float,
        help="override incident_p_GeV (appropriate for a monoenergetic sample)",
    )
    parser.add_argument(
        "--incident-events",
        type=float,
        help="optional projected number of elastic protons incident on the analyzer",
    )
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)

    source = ROOT.TFile.Open(args.root_file)
    tree = source.Get("fpp") if source else None
    if not tree:
        raise SystemExit(f"Cannot open fpp tree in {args.root_file}")

    theta_hist = ROOT.TH1D(
        "accepted_theta", "Accepted events;#theta_{polar} (deg);weighted events", 50, 0, 25
    )
    fom_hist = ROOT.TH1D(
        "fom_theta", "Intrinsic FOM contribution;#theta_{polar} (deg);#Sigma w A_{y}^{2}", 50, 0, 25
    )
    ay_graph = ROOT.TGraph()
    ay_graph.SetName("ay_model")
    ay_graph.SetTitle("GEp-III analyzing-power model;#theta_{polar} (deg);A_{y}")

    total_weight = 0.0
    accepted_weight = 0.0
    useful_weight = 0.0
    ay2_weight = 0.0
    momentum_weight = 0.0
    for event in tree:
        weight = float(event.weight)
        total_weight += weight
        observable = accepted_observable(
            event,
            args.thickness,
            args.sclose_max,
            args.theta_switch,
            args.theta_min,
        )
        if observable is None:
            continue
        theta_deg, momentum_gev = observable
        if args.momentum_gev is not None:
            momentum_gev = args.momentum_gev
        accepted_weight += weight
        momentum_weight += weight * momentum_gev
        theta_hist.Fill(theta_deg, weight)
        ay, pt = analyzing_power(momentum_gev, theta_deg)
        if PT_MIN_GEV <= pt <= PT_MAX_GEV:
            useful_weight += weight
        ay2 = ay * ay
        ay2_weight += weight * ay2
        fom_hist.Fill(theta_deg, weight * ay2)

    mean_momentum = momentum_weight / accepted_weight if accepted_weight else 0.0
    for index in range(251):
        theta_deg = index * 0.1
        ay, _ = analyzing_power(mean_momentum, theta_deg)
        ay_graph.SetPoint(index, theta_deg, ay)

    efficiency = accepted_weight / total_weight
    useful_efficiency = useful_weight / total_weight
    intrinsic_fom = ay2_weight / total_weight
    effective_ay = math.sqrt(ay2_weight / accepted_weight) if accepted_weight else 0.0
    useful_effective_ay = math.sqrt(ay2_weight / useful_weight) if useful_weight else 0.0
    polarized_fom = args.beam_polarization ** 2 * intrinsic_fom

    cumulative_fom = fom_hist.Clone("cumulative_fom")
    cumulative_fom.Reset()
    running_fom = 0.0
    for bin_index in range(1, fom_hist.GetNbinsX() + 1):
        running_fom += fom_hist.GetBinContent(bin_index)
        cumulative_fom.SetBinContent(
            bin_index, running_fom / total_weight if total_weight else 0.0
        )
    cumulative_fom.SetTitle(
        "Cumulative intrinsic FOM;maximum #theta_{polar} (deg);"
        "#Sigma w A_{y}^{2}/#Sigma w_{incident}"
    )

    canvas = ROOT.TCanvas("fom", "fom", 1500, 1100)
    canvas.Divide(2, 2)
    canvas.cd(1)
    theta_hist.SetStats(False)
    theta_hist.SetLineWidth(2)
    theta_hist.Draw("hist")
    canvas.cd(2)
    fom_hist.SetStats(False)
    fom_hist.SetLineWidth(2)
    fom_hist.Draw("hist")
    canvas.cd(3)
    ay_graph.SetLineWidth(3)
    ay_graph.Draw("AL")
    canvas.cd(4)
    cumulative_fom.SetStats(False)
    cumulative_fom.SetLineWidth(3)
    cumulative_fom.Draw("hist")
    canvas.SaveAs(args.output_stem + ".pdf")
    canvas.SaveAs(args.output_stem + ".png")

    print("model gep3_high_energy")
    print(f"ay_formula ({AY_C0:.9g}+{AY_C1:.9g}/p)*pT*exp(-{AY_B:.9g}*pT^2)")
    print(f"pt_window_GeV {PT_MIN_GEV:g} {PT_MAX_GEV:g}")
    print(f"thickness_cm {args.thickness:g}")
    print(f"theta_switch_deg {args.theta_switch:g}")
    print(f"theta_min_deg {args.theta_min:g}")
    print(
        "momentum_source "
        + (f"override_{args.momentum_gev:g}_GeV" if args.momentum_gev is not None else "incident_p_GeV")
    )
    print(f"input_weight {total_weight:.9g}")
    print(f"accepted_weight {accepted_weight:.9g}")
    print(f"useful_pt_weight {useful_weight:.9g}")
    print(f"accepted_efficiency {efficiency:.9g}")
    print(f"useful_pt_efficiency {useful_efficiency:.9g}")
    print(f"effective_Ay_all_accepted {effective_ay:.9g}")
    print(f"effective_Ay_useful_pt {useful_effective_ay:.9g}")
    print(f"intrinsic_fom_per_incident {intrinsic_fom:.9g}")
    print(f"beam_polarization {args.beam_polarization:.9g}")
    print(f"polarized_fom_per_incident {polarized_fom:.9g}")
    if args.incident_events is not None:
        accepted_events = args.incident_events * efficiency
        useful_events = args.incident_events * useful_efficiency
        total_polarized_fom = args.incident_events * polarized_fom
        sigma_component = (
            math.sqrt(2.0 / total_polarized_fom)
            if total_polarized_fom > 0.0
            else math.inf
        )
        print(f"projected_incident_events {args.incident_events:.9g}")
        print(f"projected_accepted_events {accepted_events:.9g}")
        print(f"projected_useful_pt_events {useful_events:.9g}")
        print(f"projected_polarized_fom {total_polarized_fom:.9g}")
        print(f"projected_sigma_fpp_component {sigma_component:.9g}")


if __name__ == "__main__":
    main()
