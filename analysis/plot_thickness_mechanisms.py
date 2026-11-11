#!/usr/bin/env python3
"""Plot the five-thickness FPP mechanism decomposition from JSON summaries."""

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "accepted": "#4C78A8",
    "useful": "#E45756",
    "zero": "#BAB0AC",
    "single": "#59A14F",
    "multi": "#F28E2B",
    "primary": "#4C78A8",
    "secondary": "#E45756",
    "other": "#B279A2",
    "killed": "#E45756",
    "downstream": "#59A14F",
    "unknown": "#BAB0AC",
}


def pct(values):
    return 100.0 * np.asarray(values, dtype=float)


def style_axis(axis, thicknesses):
    axis.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.set_xlabel(r"CH$_2$ thickness (cm)")
    axis.set_xticks(thicknesses)
    axis.set_xticklabels([f"{value:g}" for value in thicknesses])


def save_figure(figure, stem):
    figure.savefig(str(stem) + ".pdf", bbox_inches="tight")
    figure.savefig(str(stem) + ".png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("combined_json")
    parser.add_argument("output_dir")
    args = parser.parse_args()

    with open(args.combined_json) as stream:
        samples = json.load(stream)["samples"]
    samples.sort(key=lambda sample: sample["thickness_cm"])

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    x = np.asarray([sample["thickness_cm"] for sample in samples])
    incident = np.asarray([sample["incident_weight"] for sample in samples])
    accepted = np.asarray([sample["accepted_efficiency"] for sample in samples])
    useful = np.asarray([sample["useful_efficiency"] for sample in samples])
    accepted_error = np.sqrt(accepted * (1.0 - accepted) / incident)
    useful_error = np.sqrt(useful * (1.0 - useful) / incident)

    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "legend.fontsize": 10,
        "figure.titlesize": 16,
    })

    figure, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), constrained_layout=True)
    figure.suptitle("Single-analyzer FPP thickness scan: event-level diagnostics")

    axis = axes[0, 0]
    axis.errorbar(x, pct(accepted), yerr=pct(accepted_error), marker="o",
                  linewidth=2.2, capsize=3, color=COLORS["accepted"],
                  label="Passes agreed track cuts")
    axis.errorbar(x, pct(useful), yerr=pct(useful_error), marker="o",
                  linewidth=2.2, capsize=3, color=COLORS["useful"],
                  label="Also in useful $p_T$ window")
    axis.axvline(70, color="#777777", linestyle="--", linewidth=1)
    axis.set_ylabel("Fraction of incident events (%)")
    axis.set_title("Accepted and useful efficiency")
    axis.legend(frameon=False)
    style_axis(axis, x)

    axis = axes[0, 1]
    for key, title in (("zero", "No front-pair track"),
                       ("single", "Exactly one track"),
                       ("multi", "Multiple tracks")):
        axis.plot(x, pct([sample["multiplicity_fraction"][key] for sample in samples]),
                  marker="o", linewidth=2.2, color=COLORS[key], label=title)
    axis.set_ylabel("Fraction of incident events (%)")
    axis.set_title("Charged-track multiplicity before cuts")
    axis.legend(frameon=False)
    style_axis(axis, x)

    axis = axes[1, 0]
    for key, title in (("downstream", "Primary exits downstream"),
                       ("killed", "Primary killed in analyzer"),
                       ("unknown", "Lost before/without classified exit")):
        axis.plot(x, pct([sample["primary_exit_fraction"][key] for sample in samples]),
                  marker="o", linewidth=2.2, color=COLORS[key], label=title)
    axis.set_ylabel("Fraction of incident events (%)")
    axis.set_title("Incident-primary fate")
    axis.legend(frameon=False)
    style_axis(axis, x)

    axis = axes[1, 1]
    for key, title, color in (("zero", "No analyzer interaction", "#4C78A8"),
                              ("one", "Exactly one", "#59A14F"),
                              ("multiple", "Multiple", "#F28E2B")):
        axis.plot(x, pct([sample["event_interaction_fraction"][key] for sample in samples]),
                  marker="o", linewidth=2.2, color=color, label=title)
    axis.set_ylabel("Fraction of incident events (%)")
    axis.set_title("All recorded hadronic interactions")
    axis.legend(frameon=False)
    style_axis(axis, x)
    save_figure(figure, output_dir / "thickness_event_diagnostics")

    figure, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), constrained_layout=True)
    figure.suptitle("Single-analyzer FPP thickness scan: useful-yield mechanisms")

    primary = np.asarray([sample["channel_efficiency"]["primary_elastic"] for sample in samples])
    secondary = np.asarray([sample["channel_efficiency"]["secondary_inelastic"] for sample in samples])
    other = useful - primary - secondary
    width = 4.0
    axis = axes[0, 0]
    axis.bar(x, pct(primary), width=width, color=COLORS["primary"], label="Primary elastic")
    axis.bar(x, pct(secondary), width=width, bottom=pct(primary),
             color=COLORS["secondary"], label="Secondary inelastic")
    axis.bar(x, pct(other), width=width, bottom=pct(primary + secondary),
             color=COLORS["other"], label="Other")
    axis.errorbar(x, pct(useful), yerr=pct(useful_error), fmt="ko", capsize=3,
                  label="Total useful")
    axis.set_ylabel("Efficiency per incident event (%)")
    axis.set_title("Useful-efficiency decomposition")
    axis.legend(frameon=False)
    style_axis(axis, x)

    axis = axes[0, 1]
    useful_single = np.asarray([
        sample["useful_multiplicity_efficiency"]["single"] for sample in samples])
    useful_multi = np.asarray([
        sample["useful_multiplicity_efficiency"]["multi"] for sample in samples])
    axis.bar(x, pct(useful_single), width=width, color=COLORS["single"],
             label="Useful single-track")
    axis.bar(x, pct(useful_multi), width=width, bottom=pct(useful_single),
             color=COLORS["multi"], label="Useful multi-track")
    axis.errorbar(x, pct(useful), yerr=pct(useful_error), fmt="ko", capsize=3)
    axis.set_ylabel("Efficiency per incident event (%)")
    axis.set_title("Useful yield by event multiplicity")
    axis.legend(frameon=False)
    style_axis(axis, x)

    axis = axes[1, 0]
    for key, title, color in (
        ("single", "Single-candidate events", COLORS["single"]),
        ("multi", "Multiple-candidate events", COLORS["multi"]),
        ("combined", "Combined", "#B279A2"),
    ):
        fom = np.asarray([
            sample["fom_by_multiplicity"][key]["intrinsic_fom"]
            for sample in samples
        ])
        fom_error = np.asarray([
            sample["fom_by_multiplicity"][key]["intrinsic_fom_mc_se"]
            for sample in samples
        ])
        axis.errorbar(x, 1.0e3 * fom, yerr=1.0e3 * fom_error, marker="o",
                      linewidth=2.2, capsize=3, color=color, label=title)
    axis.axvline(70, color="#777777", linestyle="--", linewidth=1)
    axis.set_ylabel(r"Intrinsic FOM contribution ($\times 10^{-3}$)")
    axis.legend(frameon=False)
    axis.set_title("FOM per incident proton by raw multiplicity")
    style_axis(axis, x)

    axis = axes[1, 1]
    particle_series = (
        ("primary_proton", "Primary proton", "#4C78A8"),
        ("secondary_proton", "Secondary protons", "#E45756"),
        ("other_secondary_charged", "Other charged secondaries", "#F28E2B"),
        ("neutral", "Neutral particles", "#BAB0AC"),
    )
    for key, title, color in particle_series:
        axis.plot(x, [sample["mean_interactions_by_particle_per_incident"][key]
                      for sample in samples], marker="o", linewidth=2.2,
                  color=color, label=title)
    axis.set_yscale("log")
    axis.set_ylabel("Mean interactions per incident event")
    axis.set_title("Interaction activity by particle class")
    axis.legend(frameon=False)
    style_axis(axis, x)
    save_figure(figure, output_dir / "thickness_mechanism_decomposition")


if __name__ == "__main__":
    main()
