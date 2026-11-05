#!/usr/bin/env python3
"""Generate deterministic monoenergetic proton input for transport closure tests."""

import argparse
import csv
import random


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("--events", type=int, required=True)
    parser.add_argument("--momentum", type=float, required=True, help="GeV/c")
    parser.add_argument("--sigma-x", type=float, default=0.0, help="cm")
    parser.add_argument("--sigma-y", type=float, default=0.0, help="cm")
    parser.add_argument("--sigma-theta", type=float, default=0.0, help="mrad per projection")
    parser.add_argument("--seed", type=int, default=20261105)
    args = parser.parse_args()

    if args.events <= 0 or args.momentum <= 0.0:
        parser.error("events and momentum must be positive")

    rng = random.Random(args.seed)
    with open(args.output, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["event_id", "weight", "x_cm", "y_cm", "z_cm", "px_GeV", "py_GeV", "pz_GeV"]
        )
        for event_id in range(args.events):
            x = rng.gauss(0.0, args.sigma_x)
            y = rng.gauss(0.0, args.sigma_y)
            tx = rng.gauss(0.0, args.sigma_theta * 1.0e-3)
            ty = rng.gauss(0.0, args.sigma_theta * 1.0e-3)
            norm = (1.0 + tx * tx + ty * ty) ** 0.5
            writer.writerow(
                [
                    event_id,
                    1.0,
                    f"{x:.9g}",
                    f"{y:.9g}",
                    0.0,
                    f"{args.momentum * tx / norm:.9g}",
                    f"{args.momentum * ty / norm:.9g}",
                    f"{args.momentum / norm:.9g}",
                ]
            )


if __name__ == "__main__":
    main()
