#!/usr/bin/env python3
# about: transit times of known planets from TESS light curves, for timing studies (targets in targets/timing/)
"""timing.py - mid-transit times of a known planet from TESS light curves (Astro folder tool).

    bash run.sh timing --self-test              # model and error checks on synthetic data (seconds)
    bash run.sh timing --list-targets           # the planets in targets/timing/
    bash run.sh timing wasp22                   # download the light curves (first time) and time every transit
    bash run.sh timing wasp22 --fetch-only      # only download
    bash run.sh timing wasp22 --both            # also download the full-frame-image light curves where 2-min exist

What it does
  1. Reads TESS SPOC 2-min (*_lc.fits) or TESS-SPOC FFI (hlsp_tess-spoc_*_llc.fits) light curves:
     PDCSAP flux, quality flags removed, each sector normalised. Time = BTJD + 2457000 = BJD_TDB.
  2. Cuts a window of +-1.5 transit durations around each predicted transit and keeps only
     complete transits: ingress, mid-transit and egress covered, no data gap longer than 30 min
     between first and last contact, and at least 80% of the baseline on each side.
  3. Fits one transit shape to all of them together (Rp/R*, a/R*, b, quadratic limb darkening,
     plus a linear baseline per transit, solved exactly at each step).
  4. With that shape fixed, fits each transit's own mid-time (and its own baseline).
     Errors: from the chi-square curvature, scaled up by sqrt(reduced chi2) when > 1 and by the
     red-noise factor beta (Winn et al. 2008: rms of binned residuals vs white-noise expectation).
  5. Writes <out>/times.csv (BJD_TDB), <out>/fits.txt and a per-transit residual summary.

Data: SPOC 2-min light curves where they exist, otherwise TESS-SPOC full-frame-image light curves,
downloaded with the planet-search tool's own archive code (tools/transits.py: HTTPS only, MAST and
archive.stsci.edu, plus STScI's public S3 bucket after a MAST redirect). The FITS files are kept in
data/timing/<name>/ so the timing can be rerun offline. Only numpy, scipy and astropy are used.
The O-C analysis and its figure are done separately.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

VERSION = "0.5"           # 0.5: shape fit in offsets (converges), exact inner disc for z<p, MAST file names; 0.4: mid-times off the 15-s grid; 0.3: gap rule
BJD0 = 2457000.0
TOOL = "timing"
HERE = Path(__file__).resolve().parent      # Astro/tools
ROOT = HERE.parent                          # Astro
TARGET_DIR = ROOT / "targets" / TOOL
DATA_DIR = ROOT / "data" / TOOL
OUT_DIR = ROOT / "output" / TOOL
CREDITS = """This uses data collected by the TESS mission, obtained from the MAST data archive at the Space Telescope
Science Institute (STScI). Funding for the TESS mission is provided by the NASA Explorer Program. STScI is
operated by the Association of Universities for Research in Astronomy, Inc., under NASA contract NAS 5-26555.
Light curves: TESS SPOC (Jenkins et al. 2016) and TESS-SPOC full-frame-image light curves (Caldwell et al. 2020).
Transit times: timing.py (Astro folder), written for Stefano Ivan Stinga."""

# ----------------------------------------------------------------------------------------------
# Transit model: quadratic limb darkening, circular orbit, exact geometry by radial integration
# ----------------------------------------------------------------------------------------------
_GL_X, _GL_W = np.polynomial.legendre.leggauss(96)


def occulted_flux(z, p, u1, u2):
    """Fraction of the stellar flux hidden by a dark disc of radius p (stellar radii) whose centre
    is z stellar radii from the star's centre. Integrates I(r) times the length of the circle of
    radius r that lies inside the planet's disc, over r in [|z-p|, min(1, z+p)]; when z < p the fully
    covered inner disc r < p - z is added in closed form (no kink inside the quadrature)."""
    z = np.asarray(z, float)
    out = np.zeros_like(z)
    m = z < 1.0 + p
    if not m.any():
        return out
    zz = z[m]
    lo = np.abs(zz - p)          # for z < p, the annuli r < p - z are fully covered (added analytically below)
    hi = np.minimum(1.0, zz + p)
    # nodes mapped with r = lo + (hi-lo) * s, s = (1 - cos(pi*x'))/2 clusters nodes at both ends,
    # where the integrand has square-root behaviour
    x = 0.5 * (_GL_X + 1.0)
    s = 0.5 * (1.0 - np.cos(np.pi * x))
    ds = 0.5 * np.pi * np.sin(np.pi * x)
    r = lo[:, None] + (hi - lo)[:, None] * s[None, :]
    w = (hi - lo)[:, None] * ds[None, :] * 0.5 * _GL_W[None, :]
    rz = np.maximum(r * zz[:, None], 1e-300)
    cosang = (r * r + zz[:, None] ** 2 - p * p) / (2.0 * rz)
    theta = np.arccos(np.clip(cosang, -1.0, 1.0))
    mu = np.sqrt(np.clip(1.0 - r * r, 0.0, 1.0))
    inten = 1.0 - u1 * (1.0 - mu) - u2 * (1.0 - mu) ** 2
    blocked = np.sum(inten * 2.0 * r * theta * w, axis=1)
    # fully covered inner disc r < p - z (only when z < p): integral of I(r) 2 pi r dr, done in mu
    inner = zz < p
    if inner.any():
        def G(mu):
            return 2.0 * np.pi * ((1 - u1 - u2) * mu ** 2 / 2 + (u1 + 2 * u2) * mu ** 3 / 3 - u2 * mu ** 4 / 4)
        mu0 = np.sqrt(np.clip(1.0 - (p - zz[inner]) ** 2, 0.0, 1.0))
        blocked[inner] += G(1.0) - G(mu0)
    total = np.pi * (1.0 - u1 / 3.0 - u2 / 6.0)
    out[m] = blocked / total
    return out


def transit_model(t, tc, per, p, a, b, u1, u2, exptime=0.0, nsub=1):
    """Relative flux at times t (days) for a transit centred on tc; exposure-averaged over
    nsub sub-samples spanning exptime (days)."""
    t = np.asarray(t, float)
    if nsub > 1 and exptime > 0:
        off = (np.arange(nsub) + 0.5) / nsub - 0.5
        tt = t[:, None] + off[None, :] * exptime
    else:
        tt = t[:, None]
    ph = 2.0 * np.pi * (tt - tc) / per
    cosi = b / a
    x = a * np.sin(ph)
    y = a * cosi * np.cos(ph)
    z = np.sqrt(x * x + y * y)
    z = np.where(np.cos(ph) > 0, z, 10.0)          # planet behind the star: no transit
    f = 1.0 - occulted_flux(z.ravel(), p, u1, u2).reshape(z.shape)
    return f.mean(axis=1)


def duration_days(per, p, a, b):
    """Total duration T14 for a circular orbit."""
    arg = np.sqrt(max((1 + p) ** 2 - b * b, 0.0)) / (a * np.sqrt(1 - (b / a) ** 2))
    return per / np.pi * np.arcsin(min(arg, 1.0))


# ----------------------------------------------------------------------------------------------
# Light curves
# ----------------------------------------------------------------------------------------------
def read_lc(path):
    """TESS SPOC or TESS-SPOC light curve -> dict with t (BJD_TDB), f, e (normalised), sector,
    exptime (days)."""
    from astropy.io import fits
    with fits.open(path, memmap=False) as h:
        hd0, hd1 = h[0].header, h[1].header
        d = h[1].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d["PDCSAP_FLUX"], float)
        e = np.asarray(d["PDCSAP_FLUX_ERR"], float)
        q = np.asarray(d["QUALITY"], int)
        sector = int(hd0.get("SECTOR", -1))
        tid = int(hd0.get("TICID", 0) or hd0.get("OBJECT", "TIC 0").split()[-1])
        exptime = float(hd1.get("TIMEDEL", hd1.get("INT_TIME", 120.0 / 86400.0)))
        tsys = str(hd1.get("TIMESYS", "TDB")).strip()
        bjdref = float(hd1.get("BJDREFI", 2457000)) + float(hd1.get("BJDREFF", 0.0))
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (q == 0) & (e > 0)
    t, f, e = t[ok] + bjdref, f[ok], e[ok]
    med = np.median(f)
    return {"t": t, "f": f / med, "e": e / med, "sector": sector, "tic": tid, "exptime": exptime,
            "timesys": tsys, "path": str(path)}


def cut_windows(lcs, t0, per, t14, half=1.5, min_cover=0.8, max_gap_min=30.0):
    """Windows of +-half*T14 around each predicted transit. A transit is kept when the
    out-of-transit data on both sides cover at least min_cover of their span, ingress, mid-transit
    and egress all have points, and no gap between first and last contact exceeds max_gap_min."""
    wins = []
    for lc in lcs:
        t = lc["t"]
        n0 = int(np.ceil((t.min() - t0) / per - 0.5))
        n1 = int(np.floor((t.max() - t0) / per + 0.5))
        cad = np.median(np.diff(t))
        for n in range(n0, n1 + 1):
            tc = t0 + n * per
            m = np.abs(t - tc) < half * t14
            if m.sum() < 10:
                continue
            tw = t[m]
            pre = (tw < tc - 0.5 * t14)
            post = (tw > tc + 0.5 * t14)
            span_side = (half - 0.5) * t14
            cov_pre = pre.sum() * cad / span_side
            cov_post = post.sum() * cad / span_side
            ing = np.abs(tw - (tc - 0.5 * t14)) < 0.12 * t14
            egr = np.abs(tw - (tc + 0.5 * t14)) < 0.12 * t14
            mid = np.abs(tw - tc) < 0.3 * t14
            flags = []
            if cov_pre < min_cover:
                flags.append(f"pre-transit coverage {cov_pre:.2f}")
            if cov_post < min_cover:
                flags.append(f"post-transit coverage {cov_post:.2f}")
            if ing.sum() < 3:
                flags.append("ingress missing")
            if egr.sum() < 3:
                flags.append("egress missing")
            if mid.sum() < 5:
                flags.append("mid-transit missing")
            inside = np.sort(tw[np.abs(tw - tc) <= 0.5 * t14])
            edges = np.concatenate([[tc - 0.5 * t14], inside, [tc + 0.5 * t14]])
            gap = float(np.max(np.diff(edges))) * 1440.0
            if gap > max_gap_min:
                flags.append(f"{gap:.0f}-min gap in transit")
            wins.append({"n": n, "tc_pred": tc, "t": tw, "f": lc["f"][m], "e": lc["e"][m],
                         "sector": lc["sector"], "exptime": lc["exptime"], "flags": flags})
    return wins


def _nsub(exptime_days):
    return max(1, int(math.ceil(exptime_days * 86400.0 / 40.0)))   # sub-samples every <=40 s


def _baseline_solve(t, f, e, model, tc):
    """Best multiplicative linear baseline c0 + c1*(t-tc) for flux = model * baseline."""
    w = 1.0 / e
    A = np.column_stack([model * w, model * (t - tc) * w])
    c, *_ = np.linalg.lstsq(A, f * w, rcond=None)
    return model * (c[0] + c[1] * (t - tc)), c


# ----------------------------------------------------------------------------------------------
# Fits
# ----------------------------------------------------------------------------------------------
def q_to_u(q1, q2):
    s = math.sqrt(max(q1, 0.0))
    return 2.0 * s * q2, s * (1.0 - 2.0 * q2)


def fit_shape(wins, per, t0, guess):
    """Joint fit of the shape (and a linear ephemeris) to all usable windows."""
    good = [w for w in wins if not w["flags"]]
    if len(good) < 2:
        raise SystemExit("fewer than 2 complete transits: cannot fit the shape")
    nmid = np.median([w["n"] for w in good])

    tc_mid = t0 + nmid * per

    def resid(th):
        dtc_min, dper_s, p, a, b, q1, q2 = th
        tc0 = tc_mid + dtc_min / 1440.0          # offsets, not absolute BJD: the optimiser's relative
        dper = dper_s / 86400.0                  # steps and tolerances are meaningless at JD 2.46e6
        u1, u2 = q_to_u(q1, q2)
        out = []
        for w in good:
            tc = tc0 + (w["n"] - nmid) * (per + dper)
            mod = transit_model(w["t"], tc, per, p, a, b, u1, u2, w["exptime"], _nsub(w["exptime"]))
            full, _ = _baseline_solve(w["t"], w["f"], w["e"], mod, tc)
            out.append((w["f"] - full) / w["e"])
        return np.concatenate(out)

    x0 = [0.0, 0.0, guess["p"], guess["a"], guess["b"], guess.get("q1", 0.36), guess.get("q2", 0.3)]
    lo = [-72.0, -86.4, 0.01, 1.5, 0.0, 0.0, 0.0]
    hi = [72.0, 86.4, 0.5, 50.0, 0.95, 1.0, 1.0]
    # central differences with a relative step of 1e-5: the default tiny forward steps let the model's
    # ~1e-9 quadrature noise swamp the gradient along the weak (b, a, limb-darkening) directions, and the
    # fit stopped early. Restart from the result until chi2 stops improving.
    x, best, r = np.array(x0, float), None, None
    for _ in range(4):
        r = least_squares(resid, x, bounds=(lo, hi), x_scale=[0.5, 0.5, 1e-3, 0.1, 0.05, 0.05, 0.05],
                          jac="3-point", diff_step=1e-5, xtol=1e-12, ftol=1e-12, gtol=1e-10, max_nfev=2000)
        c = float(np.sum(r.fun ** 2))
        improved = best is None or best - c > 1e-3
        x, best = r.x, c if best is None else min(best, c)
        if not improved:
            break
    dtc_min, dper_s, p, a, b, q1, q2 = x
    tc0, dper = tc_mid + dtc_min / 1440.0, dper_s / 86400.0
    u1, u2 = q_to_u(q1, q2)
    chi2 = float(np.sum(r.fun ** 2))
    return {"tc_ref": tc0, "n_ref": nmid, "per": per + dper, "p": p, "a": a, "b": b, "q1": q1, "q2": q2,
            "u1": u1, "u2": u2, "chi2": chi2, "npts": int(r.fun.size), "ntransits": len(good),
            "t14": duration_days(per + dper, p, a, b), "fit_status": int(r.status), "fit_nfev": int(r.nfev)}


def beta_factor(t, res, e, bins_min=(10, 15, 20, 30, 40)):
    """Red-noise factor: rms of binned residuals over the white-noise expectation (>= 1)."""
    sig1 = np.std(res)
    betas = []
    dt = np.median(np.diff(t))
    for bm in bins_min:
        nb = max(1, int(round(bm / (dt * 1440.0))))
        if nb < 2 or len(res) < 4 * nb:
            continue
        k = len(res) // nb
        rb = res[:k * nb].reshape(k, nb).mean(axis=1)
        expect = sig1 / math.sqrt(nb) * math.sqrt(k / max(k - 1, 1))
        betas.append(np.std(rb) / expect)
    return max(1.0, float(np.max(betas))) if betas else 1.0      # conservative: worst bin size


def time_one(w, shape, search_min=20.0):
    """Mid-time of one transit with the shape fixed."""
    per, p, a, b, u1, u2 = (shape[k] for k in ("per", "p", "a", "b", "u1", "u2"))
    nsub = _nsub(w["exptime"])

    def chi2(tc):
        mod = transit_model(w["t"], tc, per, p, a, b, u1, u2, w["exptime"], nsub)
        full, _ = _baseline_solve(w["t"], w["f"], w["e"], mod, tc)
        return float(np.sum(((w["f"] - full) / w["e"]) ** 2))

    tc0 = w["tc_pred"]
    # Work in offsets from the predicted time (minutes). Minimising over the absolute BJD would make
    # the optimiser's relative tolerance ~20 s at JD 2.46e6 and leave the result on the grid (v0.3 bug).
    off = lambda x_min: chi2(tc0 + x_min / 1440.0)
    xs = np.linspace(-search_min, search_min, 161)
    c = np.array([off(x) for x in xs])
    k = int(np.argmin(c))
    if k in (0, len(xs) - 1):
        return None
    r = minimize_scalar(off, bounds=(xs[k - 1], xs[k + 1]), method="bounded", options={"xatol": 1e-5})
    tc = tc0 + float(r.x) / 1440.0
    cmin = chi2(tc)
    h = 20.0 / 86400.0
    curv = (chi2(tc + h) - 2 * cmin + chi2(tc - h)) / h ** 2
    sig = math.sqrt(2.0 / curv) if curv > 0 else float("nan")
    dof = max(len(w["t"]) - 3, 1)
    red = cmin / dof
    mod = transit_model(w["t"], tc, per, p, a, b, u1, u2, w["exptime"], nsub)
    full, base = _baseline_solve(w["t"], w["f"], w["e"], mod, tc)
    res = w["f"] - full
    beta = beta_factor(w["t"], res, w["e"])
    sig_final = sig * math.sqrt(max(red, 1.0)) * beta
    return {"n": w["n"], "sector": w["sector"], "tc": tc, "err": sig_final, "err_white": sig,
            "redchi2": red, "beta": beta, "npts": len(w["t"]), "rms_ppm": float(np.std(res) * 1e6),
            "depth_scale": float(base[0])}


# ----------------------------------------------------------------------------------------------
# Self-test on synthetic data
# ----------------------------------------------------------------------------------------------
def self_test(verbose=True):
    rng = np.random.default_rng(1)
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        if verbose:
            print(f"  [{'ok' if cond else 'FAIL'}] {name} {detail}")

    # 1. uniform source against the analytic overlap area
    p = 0.1
    for z in (0.0, 0.5, 0.95, 1.0, 1.05, 1.09):
        if z <= 1 - p:
            exact = p * p
        elif z >= 1 + p:
            exact = 0.0
        else:
            k0 = math.acos((p * p + z * z - 1) / (2 * p * z))
            k1 = math.acos((1 - p * p + z * z) / (2 * z))
            exact = (p * p * k0 + k1 - 0.5 * math.sqrt(4 * z * z - (1 + z * z - p * p) ** 2)) / math.pi
        got = occulted_flux(np.array([z]), p, 0.0, 0.0)[0]
        check(f"uniform disc z={z}", abs(got - exact) < 2e-6, f"{got:.7f} vs {exact:.7f}")
    # 2. limb-darkened centre depth: small-planet limit p^2 * I(0)/<I>
    p = 0.01
    u1, u2 = 0.4, 0.25
    got = occulted_flux(np.array([0.0]), p, u1, u2)[0]
    exact = p * p * 1.0 / (1 - u1 / 3 - u2 / 6)
    check("limb-darkened centre (small planet)", abs(got / exact - 1) < 1e-3, f"{got:.4e} vs {exact:.4e}")
    # 3. symmetric in time
    tt = np.linspace(-0.1, 0.1, 201)
    m = transit_model(tt, 0.0, 3.5, 0.1, 8.0, 0.3, 0.4, 0.25)
    check("symmetric light curve", np.max(np.abs(m - m[::-1])) < 1e-9)
    # 2b. planet covering the stellar centre (z < p): compare with adaptive quadrature
    from scipy.integrate import quad
    p_, u1_, u2_ = 0.1, 0.3, 0.25
    for zt in (0.02, 0.05, 0.09):
        def integrand(rr):
            mu_ = math.sqrt(max(1 - rr * rr, 0.0))
            th = math.pi if rr < p_ - zt else math.acos(max(-1.0, min(1.0, (rr * rr + zt * zt - p_ * p_) / (2 * rr * zt))))
            return (1 - u1_ * (1 - mu_) - u2_ * (1 - mu_) ** 2) * 2 * rr * th
        ref = quad(integrand, 0, p_ + zt, points=[p_ - zt], limit=200, epsabs=1e-14)[0] / (math.pi * (1 - u1_ / 3 - u2_ / 6))
        got = occulted_flux(np.array([zt]), p_, u1_, u2_)[0]
        check(f"limb-darkened, centre covered z={zt}", abs(got - ref) < 1e-9, f"{got:.10f} vs {ref:.10f}")
    # 3b. precision: a noiseless transit placed between grid points must come back exactly
    per, p, a, b, u1, u2 = 3.5327, 0.1, 8.5, 0.15, 0.3, 0.25
    exptime = 120.0 / 86400.0
    tt = np.arange(2461200.0, 2461200.6, exptime)
    true_tc = 2461200.3 + 7.3 / 86400.0
    ff = transit_model(tt, true_tc, per, p, a, b, u1, u2, exptime, 3)
    wq = {"n": 0, "tc_pred": 2461200.3, "t": tt, "f": ff, "e": np.full(tt.size, 1e-4), "sector": 0,
          "exptime": exptime, "flags": []}
    rq = time_one(wq, {"per": per, "p": p, "a": a, "b": b, "u1": u1, "u2": u2})
    dq = (rq["tc"] - true_tc) * 86400.0
    check("timing precision (off-grid injection)", abs(dq) < 0.3, f"recovered offset error {dq:+.3f} s")
    # 4. recovery of injected times with TESS-like noise
    per, p, a, b, u1, u2 = 3.5327, 0.1035, 8.2, 0.20, 0.38, 0.24
    t14 = duration_days(per, p, a, b)
    exptime = 120.0 / 86400.0
    true_shift = rng.normal(0, 60, 12) / 86400.0
    lcs = []
    t0 = 2461205.0
    for sec in range(2):
        t = np.arange(t0 + sec * 27.0, t0 + sec * 27.0 + 26.0, exptime)
        f = np.ones_like(t)
        for n in range(int((t[-1] - t0) / per) + 2):
            tc = t0 + n * per + true_shift[n % 12]
            mm = np.abs(t - tc) < 0.3
            f[mm] = transit_model(t[mm], tc, per, p, a, b, u1, u2, exptime, 3)
        sig = 1800e-6
        red = np.convolve(rng.normal(0, 400e-6, t.size), np.ones(15) / 15, mode="same") * np.sqrt(15)
        f = f * (1 + 2e-4 * np.sin(2 * np.pi * (t - t0) / 5.0)) + rng.normal(0, sig, t.size) + red
        lcs.append({"t": t, "f": f, "e": np.full(t.size, sig), "sector": 900 + sec, "exptime": exptime})
    wins = cut_windows(lcs, t0, per, t14)
    shape = fit_shape(wins, per, t0, {"p": 0.1, "a": 8.0, "b": 0.3})
    shape2 = fit_shape(wins, per, t0, {"p": 0.11, "a": 7.5, "b": 0.1})
    check("shape fit converges to the same solution from two starts",
          abs(shape["tc_ref"] - shape2["tc_ref"]) * 86400 < 1.0 and abs(shape["chi2"] - shape2["chi2"]) < 0.5,
          f"dtc {abs(shape['tc_ref'] - shape2['tc_ref']) * 86400:.3f} s, dchi2 {abs(shape['chi2'] - shape2['chi2']):.4f} "
          f"(b {shape['b']:.3f} vs {shape2['b']:.3f}: low b is weakly constrained), status {shape['fit_status']}")
    check("shape: Rp/R*", abs(shape["p"] - p) < 0.003, f"{shape['p']:.4f} (true {p})")
    check("shape: period", abs(shape["per"] - per) < 2e-4, f"{shape['per']:.6f}")
    pulls, errs = [], []
    for w in wins:
        if w["flags"]:
            continue
        r = time_one(w, shape)
        true = t0 + w["n"] * per + true_shift[w["n"] % 12]
        pulls.append((r["tc"] - true) / r["err"])
        errs.append(r["err"] * 86400)
    pulls = np.array(pulls)
    check("timings found", len(pulls) >= 10, f"{len(pulls)} transits")
    check("no bias", abs(np.mean(pulls)) < 3 / math.sqrt(len(pulls)), f"mean pull {np.mean(pulls):+.2f}")
    check("errors calibrated", 0.5 < np.std(pulls) < 1.5, f"pull std {np.std(pulls):.2f}, median error {np.median(errs):.0f} s")
    print("self-test", "passed" if ok else "FAILED")
    return ok


# ----------------------------------------------------------------------------------------------
# Download (through tools/transits.py)
# ----------------------------------------------------------------------------------------------
def transits_module():
    """The planet-search tool, imported as a library (its main() only runs as a script)."""
    import importlib.util
    path = HERE / "transits.py"
    if not path.exists():
        raise SystemExit("tools/transits.py is missing: the download step uses its archive code")
    spec = importlib.util.spec_from_file_location("astro_transits", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fetch(cfg, both=False, log=print):
    """Download the light curves of cfg['tic'] for cfg['sectors'] into data/timing/<name>/.
    SPOC 2-min first; TESS-SPOC (full-frame images) when there is no 2-min file (or also, with both)."""
    T = transits_module()
    d = DATA_DIR / cfg["name"]
    d.mkdir(parents=True, exist_ok=True)
    cache_path = DATA_DIR / "_cache" / "spoc_constants.json"
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        cache = json.loads(cache_path.read_text())
    except (OSError, ValueError):
        cache = {}
    tic = int(cfg["tic"])
    got = {}
    for sector in cfg["sectors"]:
        have = []
        for prod in ("spoc-2min", "tess-spoc"):
            if prod == "tess-spoc" and have and not both:
                break
            dest = d / f"{tic}_s{int(sector):04d}_{prod}_lc.fits"
            if dest.exists():
                have.append(prod)
                continue
            try:
                if prod == "spoc-2min":
                    c = T.spoc_constants(sector, cache)
                    cache_path.write_text(json.dumps(cache, indent=1))
                    if c is None:
                        log(f"  sector {sector}: no public 2-min light curves for this sector")
                        continue
                    url = T.spoc_url(tic, sector, c)
                else:
                    url = T.ffi_url(tic, sector)
                T.http_file(url, dest, max_bytes=250_000_000)
                T.read_lc_fits(dest, expect_tic=tic, expect_sector=sector)      # right star, right sector
                have.append(prod)
                log(f"  sector {sector}: {prod} light curve downloaded ({dest.stat().st_size / 1e6:.1f} MB)")
            except T.Missing:
                log(f"  sector {sector}: no {prod} light curve for TIC {tic}")
            except Exception as exc:
                if dest.exists():
                    dest.unlink()
                log(f"  sector {sector}: {prod} failed: {exc}")
        got[sector] = have
    (d / "downloads.json").write_text(json.dumps({"tic": tic, "products": got}, indent=1))
    return got


def list_targets():
    for p in sorted(TARGET_DIR.glob("*.json")):
        try:
            c = json.loads(p.read_text())
            print(f"  {p.stem:12s} {c.get('planet', '')}  TIC {c.get('tic')}  sectors {c.get('sectors')}")
        except (OSError, ValueError) as exc:
            print(f"  {p.stem:12s} unreadable: {exc}")


def pick_files(d: Path, prefer: str):
    """One light curve per sector: the preferred product if present."""
    files = {}
    for p in sorted(set(d.glob("*_lc.fits")) | set(d.glob("*_llc.fits"))):
        m = p.stem.split("_")
        if len(m) >= 3 and m[1].startswith("s") and m[1][1:].isdigit():      # this tool's own names
            sector, prod = int(m[1][1:]), m[2]
        else:                                                                  # MAST names: read the header
            from astropy.io import fits
            try:
                sector = int(fits.getheader(p, 0).get("SECTOR"))
            except Exception:
                continue
            prod = "tess-spoc" if "tess-spoc" in p.name else "spoc-2min"
        if sector not in files or prod == prefer:
            files[sector] = p
    return [files[k] for k in sorted(files)]


# ----------------------------------------------------------------------------------------------
def main(argv=None):
    import datetime as dt
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", help="name of a file in targets/timing/ (without .json)")
    ap.add_argument("--lc-dir", help="time the light curves in this folder instead of data/timing/<name>/")
    ap.add_argument("--fetch-only", action="store_true", help="only download")
    ap.add_argument("--both", action="store_true", help="also download TESS-SPOC where 2-min files exist")
    ap.add_argument("--product", choices=("spoc-2min", "tess-spoc"), default="spoc-2min",
                    help="which light curve to time when a sector has both (default spoc-2min)")
    ap.add_argument("--list-targets", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        return 0 if self_test() else 1
    if args.list_targets:
        list_targets()
        return 0
    if not args.target:
        ap.error("give a target name (bash run.sh timing --list-targets) or --self-test")
    tpath = TARGET_DIR / f"{args.target}.json"
    if not tpath.exists():
        raise SystemExit(f"no target file {tpath}")
    cfg = json.loads(tpath.read_text())
    cfg.setdefault("name", args.target)
    if args.lc_dir:
        lc_dir = Path(args.lc_dir)
    else:
        lc_dir = DATA_DIR / cfg["name"]
        print(f"light curves for {cfg.get('planet', cfg['name'])} (TIC {cfg['tic']}), sectors {cfg['sectors']}:")
        fetch(cfg, both=args.both)
        if args.fetch_only:
            return 0
    files = pick_files(lc_dir, args.product)
    if not files:
        raise SystemExit(f"no light curves in {lc_dir}")
    lcs = [read_lc(p) for p in files]
    lcs = [lc for lc in lcs if lc["sector"] in cfg["sectors"]]
    for lc in lcs:
        if lc["timesys"] != "TDB":
            raise SystemExit(f"{lc['path']}: TIMESYS {lc['timesys']} (expected TDB)")
    eph = cfg["ephemeris"]
    t0, per = eph["t0"], eph["per"]
    t14 = duration_days(per, cfg["guess"]["p"], cfg["guess"]["a"], cfg["guess"]["b"])
    wins = cut_windows(lcs, t0, per, t14)
    shape = fit_shape(wins, per, t0, cfg["guess"])
    out = OUT_DIR / cfg["name"] / dt.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    out.mkdir(parents=True, exist_ok=True)
    rows = [(w, time_one(w, shape) if not w["flags"] else None) for w in wins]
    with open(out / "times.csv", "w") as fh:
        fh.write("epoch_tool,sector,t_mid_bjd_tdb,err_d,err_white_d,redchi2,beta,npts,rms_ppm,flags\n")
        for w, r in rows:
            if r is None:
                why = "; ".join(w["flags"]) or "fit failed (minimum at the edge of the search)"
                fh.write(f"{w['n']},{w['sector']},,,,,,{len(w['t'])},,\"{why}\"\n")
            else:
                fh.write(f"{r['n']},{r['sector']},{r['tc']:.6f},{r['err']:.6f},{r['err_white']:.6f},"
                         f"{r['redchi2']:.3f},{r['beta']:.3f},{r['npts']},{r['rms_ppm']:.0f},\n")
    timed = [r for _, r in rows if r is not None]
    with open(out / "fits.txt", "w") as fh:
        fh.write(f"timing.py {VERSION}, {dt.datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"target {cfg['name']} ({cfg.get('planet', '')}), TIC {cfg['tic']}\n")
        fh.write(f"reference ephemeris: T0 = {t0}, P = {per} ({eph.get('source', '')})\nfiles:\n")
        for lc in lcs:
            fh.write(f"  {Path(lc['path']).name}: sector {lc['sector']}, cadence {lc['exptime'] * 86400:.0f} s, "
                     f"{lc['t'].size} good points, TIMESYS {lc['timesys']}\n")
        fh.write("shape fit (all complete transits together):\n")
        for k, v in shape.items():
            fh.write(f"  {k} = {v}\n")
        fh.write(f"transits: {len(wins)} windows, {len(timed)} timed, {len(wins) - len(timed)} skipped (see times.csv)\n")
        if timed:
            fh.write(f"median error {np.median([r['err'] for r in timed]) * 86400:.0f} s, "
                     f"median red-noise factor {np.median([r['beta'] for r in timed]):.2f}\n")
    (out / "CREDITS.txt").write_text(CREDITS + "\n")
    print(f"{len(timed)} of {len(wins)} transits timed -> {out / 'times.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
