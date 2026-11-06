"""GEp-III high-energy proton analyzing-power parameterization."""

import math


PT_MIN_GEV = 0.07
PT_MAX_GEV = 1.00
AY_B = (
    2.376 * 0.189 ** -2
    + 2.776 * 0.212 ** -2
    + 2.292 * 0.145 ** -2
) / (0.189 ** -2 + 0.212 ** -2 + 0.145 ** -2)
AY_SCALE = math.sqrt(2.0 * math.e * AY_B)
AY_C0 = 0.0342 * AY_SCALE
AY_C1 = 0.348 * AY_SCALE


def analyzing_power(momentum_gev, theta_deg):
    """Return GEp-III Ay and pT; zero outside its adopted useful pT range."""
    pt = momentum_gev * math.sin(math.radians(theta_deg))
    if momentum_gev <= 0.0 or not PT_MIN_GEV <= pt <= PT_MAX_GEV:
        return 0.0, pt
    ay = (AY_C0 + AY_C1 / momentum_gev) * pt * math.exp(-AY_B * pt * pt)
    return ay, pt
