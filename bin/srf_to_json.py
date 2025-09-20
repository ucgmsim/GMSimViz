#!/usr/bin/env python3
import json, math, argparse
from pathlib import Path
from source_modelling import srf as sm_srf

def deg_per_km_lat(): return 1.0/111.1
def deg_per_km_lon(lat): return deg_per_km_lat()/max(0.2, math.cos(math.radians(lat)))

def estimate_bounds_from_header(planes):
    polys = []
    for p in planes:
        lon0, lat0 = p["centre"]
        L, W = max(0.1, float(p["length"])), max(0.1, float(p["width"]))
        strike = float(p.get("strike", 0.0))
        dtop = float(p.get("dtop", 0.0))
        lat_k = deg_per_km_lat(); lon_k = deg_per_km_lon(lat0)
        hx, hy = 0.5*L*lon_k, 0.5*W*lat_k
        th = math.radians((90.0 - strike) % 360.0)
        corners = [(-hx,-hy),(hx,-hy),(hx,hy),(-hx,hy)]
        poly=[]
        for x,y in corners:
            xr = x*math.cos(th)-y*math.sin(th)
            yr = x*math.sin(th)+y*math.cos(th)
            poly.append((lon0+xr, lat0+yr, dtop))
        polys.append(poly)
    return polys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("srf", help="Path to new-format SRF")
    ap.add_argument("-o", "--out", required=True, help="Output JSON path")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    srf_file = sm_srf.read_srf(args.srf)
    header = srf_file.header
    planes = []
    for _, row in header.iterrows():
        def coalesce(*names, default=None):
            for n in names:
                if n in row and not (row[n] is None or (isinstance(row[n], float) and math.isnan(row[n]))):
                    return row[n]
            return default

        elon = float(coalesce("elon", "lon", "elng"))
        elat = float(coalesce("elat", "lat"))
        nstk = int(coalesce("nstk", "nstrike", "n_strike", default=1))
        ndip = int(coalesce("ndip", "n_dip", "n_dip_cells", default=1))

        # Try all the usual suspects for fault length/width (km)
        L = coalesce("flen", "flen_km", "fault_length", "length", "length_km", "len_km")
        W = coalesce("fwid", "fwid_km", "fault_width", "width", "width_km", "wid_km")

        # If not present, derive from spacings (km) × counts
        # Typical names for along-strike and along-dip cell size:
        dl = coalesce("dl", "dl_km", "dlen", "dlen_km", "dx_km", "strike_dx_km")
        dd = coalesce("dd", "dd_km", "ddip", "ddip_km", "dy_km", "dip_dy_km")

        if L is None and dl is not None:
            L = float(dl) * float(nstk)
        if W is None and dd is not None:
            W = float(dd) * float(ndip)

        # Final guard: if still missing, hard fail so you don’t silently get 1.0
        if L is None or W is None:
            raise ValueError(
                f"Could not determine plane length/width (km). "
                f"Have nstk={nstk}, ndip={ndip}, dl={dl}, dd={dd}. "
                f"Header columns: {list(header.columns)}"
            )

        planes.append({
            "centre": [elon, elat],
            "nstrike": nstk,
            "ndip": ndip,
            "length": float(L),
            "width": float(W),
            "strike": float(coalesce("strike", "stk", "strike_deg", default=0.0)),
            "dip": float(coalesce("dip", "dip_deg", default=45.0)),
            "shyp": float(coalesce("shyp", default=0.5)),
            "dhyp": float(coalesce("dhyp", default=0.5)),
            "dtop": float(coalesce("dtop", "depth_top", "top_depth", default=0.0)),
        })

    print("Header LxW (km):", planes[0]["length"], "x", planes[0]["width"])

    # hypocentre: earliest tinit if present, else plane[0] top
    pts = srf_file.points
    if "tinit" in pts.columns:
        idx = pts["tinit"].idxmin()
        lon = float(pts.loc[idx, "lon"] if "lon" in pts.columns else pts.loc[idx, "elon"])
        lat = float(pts.loc[idx, "lat"] if "lat" in pts.columns else pts.loc[idx, "elat"])
        dep = float(pts.loc[idx, "dep"] if "dep" in pts.columns else pts.loc[idx, "depth"])
        hypo = [lon, lat, dep]
    else:
        p0 = planes[0]
        hypo = [p0["centre"][0], p0["centre"][1], float(p0.get("dtop",0.0))]

    bounds = estimate_bounds_from_header(planes)

    out = {"planes": planes, "hypo": hypo, "bounds": bounds, "meta": {}}
    if args.title:
        out["meta"]["title"] = args.title

    Path(args.out).write_text(json.dumps(out, indent=2))
    print(f"Wrote {args.out} with {len(planes)} plane(s).")

if __name__ == "__main__":
    main()
