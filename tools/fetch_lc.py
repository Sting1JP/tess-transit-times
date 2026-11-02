#!/usr/bin/env python3
"""Download the TESS light curves of a timing.py target from MAST.

    python3 tools/fetch_lc.py wasp22                 # the sectors listed in targets/timing/wasp22.json
    python3 tools/fetch_lc.py wasp22 --both          # also the TESS-SPOC full-frame-image light curves
    python3 tools/fetch_lc.py --tic 257567854 --sectors 4 105 --out data/timing/wasp22

For each sector it downloads the SPOC 2-min light curve. When a sector has none, it downloads the
TESS-SPOC full-frame-image light curve instead (with --both, always). The files are saved as
<TIC>_s<sector>_spoc-2min_lc.fits and <TIC>_s<sector>_tess-spoc_lc.fits, the names timing.py reads,
in data/timing/<name>/ unless --out is given. Files already there are checked and kept. The script
stops if the folder holds light curves of another star, and it exits with status 1 if a sector got no
light curve or a download failed. Then:

    python3 tools/timing.py wasp22 --lc-dir data/timing/wasp22

The SPOC file name is looked up in STScI's public S3 copy of the TESS archive; every file is
downloaded through the MAST download service and checked for the right star and sector.
"""
import argparse
import json
import re
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

SSL_CONTEXT = ssl.create_default_context()
try:                                   # some Python installs (e.g. python.org builds on macOS) have no
    import certifi                     # certificate store of their own: add certifi's when it is installed
    SSL_CONTEXT.load_verify_locations(cafile=certifi.where())
except ImportError:
    pass

VERSION = "0.5.0"
ROOT = Path(__file__).resolve().parent.parent
USER_AGENT = f"tess-transit-times/{VERSION} (+https://github.com/Sting1JP/tess-transit-times)"
S3 = "https://stpubdata.s3.amazonaws.com"
MAST_DOWNLOAD = "https://mast.stsci.edu/api/v0.1/Download/file?uri="


class Missing(Exception):
    """The archive has no such light curve."""


def get(url, timeout=300):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=SSL_CONTEXT) as r:
            return r.read()
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise Missing(url) from None
        raise


def tic_path(tic):
    s = f"{int(tic):016d}"
    return "/".join(s[i:i + 4] for i in range(0, 16, 4))


def spoc_uri(tic, sector):
    prefix = f"tess/public/tid/s{sector:04d}/{tic_path(tic)}/"
    listing = get(f"{S3}/?list-type=2&prefix={urllib.parse.quote(prefix)}").decode()
    keys = [k for k in re.findall(r"<Key>([^<]+)</Key>", listing) if k.endswith("-s_lc.fits")]
    if not keys:
        raise Missing(prefix)
    return "mast:TESS/product/" + sorted(keys)[-1].rsplit("/", 1)[1]


def tess_spoc_uri(tic, sector):
    return (f"mast:HLSP/tess-spoc/s{sector:04d}/target/{tic_path(tic)}/"
            f"hlsp_tess-spoc_tess_phot_{int(tic):016d}-s{sector:04d}_tess_v1_lc.fits")


def check_file(path, tic, sector):
    """Header keys of a light curve; raises ValueError if it is incomplete or for another star or sector."""
    from astropy.io import fits
    with fits.open(path, memmap=False) as h:
        h.readall()
        last = h[-1].fileinfo()
        if h[-1].name != "APERTURE" or Path(path).stat().st_size < last["datLoc"] + last["datSpan"]:
            raise ValueError("incomplete file")
        h0, h1 = h[0].header, h[1].header
        got_tic, got_sector = int(h0.get("TICID", -1)), int(h0.get("SECTOR", -1))
        if got_tic != int(tic) or got_sector != int(sector):
            raise ValueError(f"file is for TIC {got_tic}, sector {got_sector}")
        for col in ("TIME", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR", "QUALITY"):
            if col not in h[1].columns.names:
                raise ValueError(f"no {col} column")
        if len(h[1].data["TIME"]) == 0:
            raise ValueError("empty light curve")
        return {"data_rel": h0.get("DATA_REL"), "procver": h0.get("PROCVER"), "timesys": h1.get("TIMESYS")}


def download(uri, dest, tic, sector):
    part = dest.with_name(dest.name + ".part")
    part.write_bytes(get(MAST_DOWNLOAD + urllib.parse.quote(uri, safe=":/")))
    try:
        info = check_file(part, tic, sector)
    except Exception:
        part.unlink()
        raise
    part.replace(dest)
    return info


def other_stars(folder, tic):
    """Light-curve files in folder whose names do not show this TIC (timing.py would read them too)."""
    names = sorted({p.name for pattern in ("*_lc.fits", "*_llc.fits") for p in folder.glob(pattern)})
    bad = []
    for name in names:
        m = re.match(r"(\d+)_s\d{4}_", name) or re.search(r"(?<!\d)(\d{16})(?!\d)", name)
        if not m or int(m.group(1)) != int(tic):
            bad.append(name)
    return bad


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", help="name of a file in targets/timing/ (without .json)")
    ap.add_argument("--tic", type=int, help="TIC number (instead of a target file)")
    ap.add_argument("--sectors", type=int, nargs="+", help="sectors (instead of a target file)")
    ap.add_argument("--out", help="folder for the light curves (default data/timing/<name>/)")
    ap.add_argument("--both", action="store_true", help="also download TESS-SPOC where SPOC 2-min files exist")
    args = ap.parse_args(argv)
    if args.target:
        if args.tic or args.sectors:
            ap.error("give either a target name or --tic and --sectors, not both")
        tpath = ROOT / "targets" / "timing" / f"{args.target}.json"
        if not tpath.exists():
            ap.error(f"no target file {tpath}")
        cfg = json.loads(tpath.read_text())
        name, tic, sectors = cfg.get("name", args.target), int(cfg["tic"]), [int(s) for s in cfg["sectors"]]
    elif args.tic and args.sectors and args.out:
        name, tic, sectors = None, args.tic, args.sectors
    else:
        ap.error("give a target name, or --tic, --sectors and --out")
    out = Path(args.out) if args.out else ROOT / "data" / "timing" / name
    out.mkdir(parents=True, exist_ok=True)
    log = out / "fetch_lc.json"
    try:
        record = json.loads(log.read_text())
    except (OSError, ValueError):
        record = {"tic": tic, "files": {}}
    record.setdefault("files", {})
    bad = other_stars(out, tic)
    if record.get("tic") != tic and not record["files"] and not bad:
        record = {"tic": tic, "files": {}}             # an earlier run for another TIC found nothing
    if record.get("tic") != tic or bad:
        print(f"{out} already holds light curves of another star "
              f"({', '.join(bad) if bad else 'TIC ' + str(record.get('tic')) + ' in fetch_lc.json'}). Use another folder.")
        return 1
    errors, missing = 0, []
    for sector in sectors:
        have, spoc_failed = [], False
        for prod, find in (("spoc-2min", spoc_uri), ("tess-spoc", tess_spoc_uri)):
            if prod == "tess-spoc" and not args.both:
                if have:
                    break
                if spoc_failed:
                    print(f"sector {sector}: TESS-SPOC not tried because the SPOC 2-min step failed; run again")
                    break
            dest = out / f"{tic}_s{sector:04d}_{prod}_lc.fits"
            if dest.exists():
                try:
                    info = check_file(dest, tic, sector)
                except Exception as exc:
                    errors += 1
                    spoc_failed = spoc_failed or prod == "spoc-2min"
                    print(f"sector {sector}: {dest.name} is there but cannot be used ({exc}); delete it and run again")
                    continue
                record["files"].setdefault(dest.name, {"uri": None, **info})
                log.write_text(json.dumps(record, indent=1) + "\n")
                have.append(prod)
                print(f"sector {sector}: {prod} already there ({dest.name})")
                continue
            try:
                uri = find(tic, sector)
                info = download(uri, dest, tic, sector)
            except Missing:
                print(f"sector {sector}: no {prod} light curve for TIC {tic}")
                continue
            except Exception as exc:
                errors += 1
                spoc_failed = spoc_failed or prod == "spoc-2min"
                print(f"sector {sector}: {prod} failed: {exc}")
                continue
            record["files"][dest.name] = {"uri": uri, **info}
            log.write_text(json.dumps(record, indent=1) + "\n")
            have.append(prod)
            print(f"sector {sector}: {prod} downloaded ({dest.stat().st_size / 1e6:.1f} MB, "
                  f"data release {info['data_rel']}, {info['procver']})")
        if not have:
            missing.append(sector)
    log.write_text(json.dumps(record, indent=1) + "\n")
    if missing:
        print("no light curve for sector " + " ".join(map(str, missing)))
    print(f"light curves in {out}")
    return 1 if errors or missing else 0


if __name__ == "__main__":
    sys.exit(main())
