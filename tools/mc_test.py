"""Monte Carlo: error calibration of timing.py at 120 s and 200 s cadence, white + red noise,
and the FITS reader on synthetic SPOC / TESS-SPOC style files."""
import sys, math, os, tempfile, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import timing as T
from astropy.io import fits
rng = np.random.default_rng(7)
per, p, a, b, u1, u2 = 3.5327294, 0.1035, 8.2, 0.20, 0.38, 0.24
t14 = T.duration_days(per, p, a, b)
for cad_s, sig in ((120, 1800e-6), (200, 1400e-6)):
    exptime = cad_s / 86400
    pulls = []
    for trial in range(6):
        t0 = 2461205.0 + rng.uniform(0, 1)
        shifts = rng.normal(0, 90, 20) / 86400
        lcs = []
        for sec in range(3):
            t = np.arange(t0 + sec * 27.3, t0 + sec * 27.3 + 25.5, exptime)
            f = np.ones_like(t)
            for n in range(int((t[-1] - t0) / per) + 2):
                tc = t0 + n * per + shifts[n % 20]
                m = np.abs(t - tc) < 0.3
                f[m] = T.transit_model(t[m], tc, per, p, a, b, u1, u2, exptime, 7)
            red = np.convolve(rng.normal(0, 500e-6, t.size), np.ones(9) / 9, mode="same") * 3
            f = f * (1 + 3e-4 * np.sin(2 * np.pi * (t - t0) / 4.3 + trial)) + rng.normal(0, sig, t.size) + red
            lcs.append({"t": t, "f": f, "e": np.full(t.size, sig), "sector": 105 + sec, "exptime": exptime})
        wins = T.cut_windows(lcs, t0, per, t14)
        shape = T.fit_shape(wins, per, t0, {"p": 0.1, "a": 8.0, "b": 0.3})
        for w in wins:
            if w["flags"]:
                continue
            r = T.time_one(w, shape)
            if r is None:
                continue
            pulls.append((r["tc"] - (t0 + w["n"] * per + shifts[w["n"] % 20])) / r["err"])
    pulls = np.array(pulls)
    print(f"cadence {cad_s} s: {len(pulls)} transits, mean pull {pulls.mean():+.2f}, pull std {pulls.std():.2f}, "
          f"|pull|>3: {np.sum(np.abs(pulls)>3)}")

# FITS reader: SPOC-like file
t = np.arange(3500.0, 3510.0, 120 / 86400)
cols = [fits.Column(name="TIME", format="D", array=t), fits.Column(name="PDCSAP_FLUX", format="E", array=np.full(t.size, 1000.0)),
        fits.Column(name="PDCSAP_FLUX_ERR", format="E", array=np.full(t.size, 2.0)), fits.Column(name="QUALITY", format="J", array=np.zeros(t.size, int))]
h0 = fits.PrimaryHDU(); h0.header["SECTOR"] = 105; h0.header["TICID"] = 257567854
h1 = fits.BinTableHDU.from_columns(cols); h1.header["TIMEDEL"] = 120 / 86400; h1.header["TIMESYS"] = "TDB"; h1.header["BJDREFI"] = 2457000; h1.header["BJDREFF"] = 0.0
fits.HDUList([h0, h1]).writeto(os.path.join(tempfile.gettempdir(), "timing_test_lc.fits"), overwrite=True)
lc = T.read_lc(os.path.join(tempfile.gettempdir(), "timing_test_lc.fits"))
print("reader:", lc["sector"], lc["tic"], round(lc["t"][0], 3), lc["timesys"], round(lc["exptime"] * 86400), "s", round(float(np.median(lc["f"])), 3))
