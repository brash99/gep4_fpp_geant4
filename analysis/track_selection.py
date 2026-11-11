"""Ideal charged-track policies; multiplicity is counted before analysis cuts."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Track:
    track_id: int
    chamber_mask: int
    theta_deg: float
    phi_rad: float
    momentum_gev: float  # Used only by the legacy reference policy.
    zclose_cm: float
    sclose_cm: float
    front_cone: bool
    all4_cone: bool
    root_track_id: int = 0
    generation: int = 0
    lineage_process_category: int = 0
    elastic_interactions: int = 0
    inelastic_interactions: int = 0
    other_hadronic_interactions: int = 0


def passes(track, thickness=55.88, sclose_max=1.0, theta_switch=12.0):
    four = track.theta_deg < theta_switch
    return (
        (not four or (track.chamber_mask & 15) == 15)
        and (track.all4_cone if four else track.front_cone)
        and 111.4215 <= track.zclose_cm <= 111.4215 + thickness
        and track.sclose_cm <= sclose_max
    )


def select(tracks, policy, thickness=55.88, sclose_max=1.0, theta_switch=12.0, min_angle_deg=1.0):
    """Choose before cuts; all inputs must cross FPP1 and FPP2.

    single requires exactly one such track, even if extra tracks fail cuts.
    min_angle ranks only tracks strictly above min_angle_deg, before other cuts.
    min_angle_uncut retains the earlier rule for regression comparisons.
    legacy reproduces the independent leading-momentum choices by topology.
    """
    if not tracks:
        return None
    if policy == "single":
        chosen = tracks[0] if len(tracks) == 1 else None
    elif policy in ("min_angle", "min_angle_uncut"):
        eligible = ([t for t in tracks if t.theta_deg > min_angle_deg]
                    if policy == "min_angle" else tracks)
        if not eligible:
            return None
        chosen = min(eligible, key=lambda t: (t.theta_deg, t.track_id))
    elif policy == "legacy":
        four = [t for t in tracks if (t.chamber_mask & 15) == 15]
        best_four = max(four, key=lambda t: t.momentum_gev) if four else None
        if best_four and best_four.theta_deg < theta_switch and passes(
            best_four, thickness, sclose_max, theta_switch
        ):
            return best_four
        chosen = max(tracks, key=lambda t: t.momentum_gev)
        if chosen.theta_deg < theta_switch:
            return None
    else:
        raise ValueError(f"Unknown policy: {policy}")
    return chosen if chosen and passes(chosen, thickness, sclose_max, theta_switch) else None
