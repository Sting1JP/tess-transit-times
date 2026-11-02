# tess-transit-times

This tool measures the mid-transit times of a known transiting planet in TESS light curves. The tool is one Python file, `tools/timing.py`.

The `timing.py` in this version is the code used for the TESS transit times of WASP-22 b in Stinga (2026, Research Notes of the AAS, 10, 297) and of HD 1397 b, WASP-161 b and WASP-193 b in Stinga (2026, doi:10.5281/zenodo.23855016). It is the same file as in the WASP-22 b data package (doi:10.5281/zenodo.23395079), and it reports itself as `timing.py 0.5`.

## What it does

For one planet and a list of TESS sectors, `timing.py`:

1. reads SPOC 2-min or TESS-SPOC full-frame-image light curves: the PDCSAP flux, without the points that have a quality flag, each sector divided by its median. Times are BTJD + 2457000 = BJD_TDB, and a file whose `TIMESYS` is not TDB is refused;
2. cuts a window of ±1.5 transit durations around each transit predicted by the reference ephemeris. A transit is kept only if ingress, mid-transit and egress have data, no gap between first and last contact is longer than 30 min, and the data cover at least 80% of the baseline on each side;
3. fits one transit shape to all the kept transits together: Rp/R★, a/R★, the impact parameter and quadratic limb darkening (the q1, q2 parameters of Kipping 2013) on a circular orbit, with small corrections to the reference time and period. The model is averaged over each exposure, and each transit has its own linear baseline, solved exactly at every step;
4. with that shape fixed, fits each transit's mid-time and baseline. The error is the χ² curvature error, multiplied by √(reduced χ²) when that is above 1 and by the red-noise factor β of Winn et al. (2008), the largest over 10–40-min bins.

It writes `times.csv`, `fits.txt` (the files used and the shape fit) and `CREDITS.txt` to a new folder `output/timing/<name>/<date>_<time>/`.

## Install

You need Python 3.9 or newer with numpy, scipy and astropy:

```
python3 -m pip install -r requirements.txt
```

Tested on Linux with Python 3.9.25 (numpy 1.26.4, scipy 1.13.1, astropy 6.0.1) and Python 3.13.18 (numpy 2.5.3, scipy 1.18.1, astropy 8.0.1), and on macOS with Python 3.9.6 (numpy 1.26.4, scipy 1.13.1, astropy 6.0.1).

## Use

Run the commands from the package folder. To check the install, run the self-test. It uses synthetic light curves and needs no internet:

```
python3 tools/timing.py --self-test
```

Then download the light curves of a target and time them, for example WASP-22 b in TESS Sectors 4 and 105–107:

```
python3 tools/fetch_lc.py wasp22
python3 tools/timing.py wasp22 --lc-dir data/timing/wasp22
```

`fetch_lc.py` downloads each sector's SPOC 2-min light curve from MAST, or the TESS-SPOC light curve when a sector has no 2-min one; `--both` downloads both. It also works without a target file: `--tic N --sectors S1 S2 ... --out DIR`. It stops if the folder already holds light curves of another star, and it exits with status 1 if a sector got no light curve or a download failed.

The folder given with `--lc-dir` must hold only the target's SPOC 2-min and TESS-SPOC light curves, directly in it (not in subfolders). `timing.py` reads every `*_lc.fits` and `*_llc.fits` file there, without checking the star, so another star's file, or a QLP or other high-level light curve, can take the place of the right one. Files downloaded some other way work under the names MAST gives them.

Always give `--lc-dir`. Without it, `timing.py` tries to download the files with a module from the author's own setup, which is not included, and stops; its `--fetch-only` and `--both` options act only in that step. Its help text shows the author's wrapper command, `bash run.sh timing ...`, which does the same as `python3 tools/timing.py ...`. The help text also names TESS-SPOC files `*_llc.fits` (MAST names them `*_lc.fits`, which the tool reads) and mentions a residual summary that the tool does not write. `CREDITS.txt` lists both light-curve products whichever was used.

When a folder holds both kinds of light curve for a sector, `--product tess-spoc` times the TESS-SPOC one instead of the 2-min one.

`python3 tools/mc_test.py` checks the error bars on synthetic light curves with 120-s and 200-s cadence. It takes a few minutes and prints, for each cadence, the standard deviation of the timing errors divided by their stated uncertainties, which should be close to 1. It leaves a small test file in the system's temporary folder.

## A new planet

Copy `targets/timing/wasp22.json` to `targets/timing/<name>.json` and edit it:

| Key | What it sets |
|---|---|
| `name` | the folder names under `data/timing/` and `output/timing/` |
| `planet` | the planet's name in the outputs |
| `tic` | the TESS Input Catalog number |
| `sectors` | the TESS sectors to use |
| `ephemeris` | `t0` (BJD_TDB) and `per` (days), which place the windows, and `source` |
| `guess` | starting values of `p` (Rp/R★), `a` (a/R★) and `b`, and optionally `q1` and `q2` |

The window width comes from the transit duration implied by `guess`. Each mid-time is searched within ±20 min of the reference ephemeris, so the ephemeris has to predict every transit in the sectors to a few minutes. The shape fit needs at least two complete transits and keeps Rp/R★ between 0.01 and 0.5, a/R★ between 1.5 and 50, the impact parameter between 0 and 0.95, and the corrections within ±72 min (reference time) and ±86.4 s (period); a starting value outside these limits stops the run, so transits with an impact parameter above 0.95 cannot be fitted. The `lc_dir` key in the example file is not read.

## Output

`times.csv` has one row per transit window with at least 10 points:

| Column | Meaning |
|---|---|
| `epoch_tool` | epoch number, counted from `t0` of the target file |
| `sector` | TESS sector |
| `t_mid_bjd_tdb` | mid-transit time (BJD_TDB) |
| `err_d` | adopted 1σ error in days: curvature error × √(reduced χ², if above 1) × β |
| `err_white_d` | the curvature error alone, in days |
| `redchi2` | reduced χ² of the transit fit |
| `beta` | red-noise factor |
| `npts` | number of points in the window |
| `rms_ppm` | residual rms in ppm |
| `flags` | why the transit was not timed (empty if it was) |

## Limits

- One planet per target, on a circular orbit. Transits of other planets in a window are not modeled.
- The tool measures times only. Ephemeris fits and O−C analysis are left to the user.
- Light curves with 30-min cadence (the full-frame images of TESS Sectors 1–26) are in practice not usable. The 30-min gap limit equals the cadence, so one missing point, or the slight stretch of the time stamps in BJD, fails a transit. The contact rule (at least three points within 0.12 transit durations of each contact) also needs transits longer than about 4–6 hours, and the red-noise factor, with its 10–40-min bins, does not act on 30-min data.
- With 10-min cadence (Sectors 27–55) the contact rule passes only for transits longer than about 1.5–2 hours.

## Known answer

`known_answers/wasp22/` holds `times.csv` and `fits.txt` from the published WASP-22 b record (doi:10.5281/zenodo.23395079). They were made with Python 3.13, numpy 2.5.3, scipy 1.18.1 and astropy 8.0.1, from the SPOC light curves of data release 42 (Sector 4) and 143–145 (Sectors 105–107). To compare a rerun:

```
python3 tools/compare_times.py known_answers/wasp22/times.csv output/timing/wasp22/DATE_TIME/times.csv
```

Put the run folder that `timing.py` prints at the end in place of `DATE_TIME`.

With those versions (Python 3.13.18 on Linux) the new `times.csv` is identical. With Python 3.9, numpy 1.26.4, scipy 1.13.1 and astropy 6.0.1 the shape fit settles slightly differently, and the times differ by less than 0.3 s: 0.07 s rms on Linux and 0.12 s rms on macOS.

## Files

- `tools/timing.py`: the tool.
- `tools/fetch_lc.py`: downloads the light curves from MAST.
- `tools/compare_times.py`: compares two `times.csv` files.
- `tools/mc_test.py`: checks the error bars on synthetic light curves.
- `targets/timing/wasp22.json`: the WASP-22 b target.
- `known_answers/wasp22/`: the published WASP-22 b times.
- `requirements.txt`, `LICENSE`, `CITATION.cff`, `SHA256SUMS`.

## Credits

This uses data collected by the TESS mission, obtained from the MAST data archive at the Space Telescope Science Institute (STScI). Funding for the TESS mission is provided by the NASA Explorer Program. STScI is operated by the Association of Universities for Research in Astronomy, Inc., under NASA contract NAS 5-26555. SPOC light curves: Jenkins et al. (2016), doi:10.17909/t9-nmc8-f686. TESS-SPOC light curves: Caldwell et al. (2020), doi:10.17909/t9-wpz1-8s54.

References: Caldwell, D. A., et al. 2020, Research Notes of the AAS, 4, 201. Jenkins, J. M., et al. 2016, Proc. SPIE, 9913, 99133E. Kipping, D. M. 2013, MNRAS, 435, 2152. Winn, J. N., et al. 2008, ApJ, 683, 1076.

## License

The code is under the MIT License (`LICENSE`). The files in `known_answers/` come from doi:10.5281/zenodo.23395079 and are under CC BY 4.0.

## How to cite

Stinga, S. I. 2026, tess-transit-times (version 0.5.0), Zenodo, doi:10.5281/zenodo.23875330

All versions: doi:10.5281/zenodo.23875329

Please also cite Stinga, S. I. 2026, Research Notes of the AAS, 10, 297, doi:10.3847/2515-5172/aeb5e5, which used this code.

Stefano Ivan Stinga, Independent researcher (ORCID 0009-0005-3860-3944), stefano@stinga.jp
