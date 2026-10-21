import json
import matplotlib
import matplotlib.colors
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.patheffects as pe
import geopandas as gpd, pandas as pd

import os
SRC = "data/china_provinces_datav.json"
TAB = "tables/table_3.3_region.csv"
OUT = "figures/fig3_region_map_5class.png"

gdf = gpd.GeoDataFrame.from_features(json.load(open(SRC, encoding="utf-8"))["features"], crs="EPSG:4326")
nine = gdf[gdf["adcode"] == "100000_JD"]
prov = gdf[gdf["adcode"] != "100000_JD"].copy()
prov["name"] = prov["name"].astype(str)

EAST = ["北京市","天津市","河北省","上海市","江苏省","浙江省","福建省","山东省","广东省","海南省","台湾省","香港特别行政区","澳门特别行政区"]
CENTRAL = ["山西省","安徽省","江西省","河南省","湖北省","湖南省"]
WEST = ["内蒙古自治区","广西壮族自治区","重庆市","四川省","贵州省","云南省","西藏自治区","陕西省","甘肃省","青海省","宁夏回族自治区","新疆维吾尔自治区"]
NORTHEAST = ["辽宁省","吉林省","黑龙江省"]
def reg(n):
    return "East" if n in EAST else "Central" if n in CENTRAL else "West" if n in WEST else "Northeast" if n in NORTHEAST else None
prov["region"] = prov["name"].map(reg)
prov = prov.dropna(subset=["region"])
prov["geometry"] = prov["geometry"].simplify(0.05, preserve_topology=True)
ALBERS = "+proj=aea +lat_1=25 +lat_2=47 +lon_0=105 +datum=WGS84 +units=m +no_defs"
prov_alb = prov.to_crs(ALBERS)
regions = prov.dissolve(by="region").reset_index().to_crs(ALBERS)
nine_p = nine.to_crs(ALBERS)
hainan_p = prov[prov["name"] == "海南省"].to_crs(ALBERS)

tab = pd.read_csv(TAB)
K = {"东北":"Northeast","东部":"East","中部":"Central","西部":"West"}
shares = {r["group"]: {K[k]: float(r[k]) for k in K} for _, r in tab.iterrows()}
gp = pd.read_csv("tables/table_3.3_group_profile.csv")
ns = {r["group"]: int(r["n"]) for _, r in gp.iterrows()}
SHORT = {"G1":"low pollution, gradual decline","G2":"moderate pollution, steady decline",
         "G3":"high pollution, partial rebound","G4":"highest pollution, partial rebound",
         "G5":"lowest pollution, steep decline"}

norm = Normalize(0, 80)
import numpy as np
_bp = plt.get_cmap("Blues")
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("Blues_trunc", _bp(np.linspace(0, 0.68, 256)))
fig = plt.figure(figsize=(7.2, 6.2)); fig.patch.set_facecolor("white")
gs = fig.add_gridspec(2, 3, hspace=0.30, wspace=0.05)
axes = [fig.add_subplot(gs[i, j]) for i, j in [(0,0),(0,1),(0,2),(1,0),(1,1)]]
cax = fig.add_subplot(gs[1, 2])

minx, miny, maxx, maxy = regions.total_bounds
POS = {"West": (0.22, 0.55), "Northeast": (0.72, 0.94), "Central": (0.50, 0.40), "East": (0.80, 0.42)}
for ax, g in zip(axes, ["G1","G2","G3","G4","G5"]):
    for _, row in regions.iterrows():
        v = shares[g][row["region"]]
        c = "#F7F7F7" if v == 0 else cmap(norm(v))
        geoms = row.geometry.geoms if row.geometry.geom_type == "MultiPolygon" else [row.geometry]
        for gg in geoms:
            if gg.is_empty: continue
            xs, ys = gg.exterior.xy
            ax.fill(xs, ys, facecolor=c, edgecolor="white", linewidth=0.25)
    for _, prow in prov_alb.iterrows():
        pg = prow.geometry
        pgeoms = pg.geoms if pg.geom_type == "MultiPolygon" else [pg]
        for pgg in pgeoms:
            if pgg.is_empty: continue
            pxs, pys = pgg.exterior.xy
            ax.plot(pxs, pys, color="#AAAAAA", linewidth=0.25, zorder=1)
    for rn, (fx, fy) in POS.items():
        _v = shares[g][rn]
        ax.text(fx, fy, f"{rn}\n{_v:.1f}%", transform=ax.transAxes, ha="center", va="center",
                fontsize=6.8, color="#777777" if _v == 0 else "#333333",
                linespacing=1.15, path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])
    ax.set_title(f"{g}  (n = {ns[g]:,}, {100*ns[g]/sum(ns.values()):.1f}%)", fontsize=8.5, loc="left", pad=2)
    ax.text(0.0, 1.005, SHORT[g], transform=ax.transAxes, fontsize=7, color="#666666", va="bottom")
    ax.set_xlim(minx, maxx); ax.set_ylim(miny, maxy); ax.set_aspect("equal"); ax.axis("off")
    ins = ax.inset_axes([0.73, 0.05, 0.23, 0.25])
    hainan_p.plot(ax=ins, color="#EFEFEF", edgecolor="#999999", linewidth=0.3)
    nine_p.plot(ax=ins, color="#555555", linewidth=0.9)
    ins.text(0.5, -0.08, "South China Sea", transform=ins.transAxes, ha="center", va="top", fontsize=5.5, color="#666666")
    bx = nine_p.total_bounds
    ins.set_xlim(bx[0]-1e5, bx[2]+1e5); ins.set_ylim(bx[1]-1e5, bx[3]+1e5)
    ins.set_xticks([]); ins.set_yticks([])
    for s in ins.spines.values(): s.set_linewidth(0.4); s.set_color("#999999")

cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_label("Share of class members in the region (%)", fontsize=7.2)
cb.ax.tick_params(labelsize=7)
pos = cax.get_position()
cax.set_position([pos.x0 + 0.02, pos.y0 + 0.30, pos.width - 0.06, 0.018])
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print("saved", OUT)
