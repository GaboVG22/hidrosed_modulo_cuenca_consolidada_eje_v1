
"""
HidroSed · Módulo Cuenca Consolidada + Eje Hidráulico Activo v1
----------------------------------------------------------------
Lee un KMZ consolidado con:
- cuenca de descarga/soporte;
- cuenca hidrológica;
- eje del cauce;
- curvas de nivel existentes;
- puntos de control.

Luego:
- calcula intercuenca = cuenca descarga - cuenca hidrológica;
- asocia el área incremental al tramo activo del eje;
- recorta curvas al corredor del eje;
- genera curvas auxiliares incrementadas/interpoladas a lo largo del eje;
- genera tabla de aporte incremental distribuido por km;
- exporta KMZ, Excel, CSV y JSON.

Nota técnica:
Las curvas auxiliares NO son topografía levantada. Son una densificación/interpolación basada en curvas existentes.
"""

from __future__ import annotations

import io
import json
import math
import re
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from shapely.geometry import (
    Point, LineString, Polygon, MultiPolygon, MultiLineString, GeometryCollection,
    mapping
)
from shapely.ops import unary_union, transform as shp_transform
try:
    from shapely.ops import substring as shp_substring
except Exception:
    shp_substring = None

try:
    from shapely.strtree import STRtree
except Exception:
    STRtree = None

from pyproj import CRS as PyCRS, Transformer

APP_TITLE = "HidroSed · Cuenca Consolidada + Eje Activo v1"

st.set_page_config(page_title=APP_TITLE, page_icon="🌊", layout="wide")


# =============================================================================
# Dataclasses
# =============================================================================

@dataclass
class PolygonItem:
    name: str
    geom_ll: object
    geom_proj: object
    area_km2: float


@dataclass
class PointItem:
    name: str
    lon: float
    lat: float
    geom_proj: Point


@dataclass
class AxisItem:
    name: str
    geom_ll: LineString
    geom_proj: LineString
    length_km: float


@dataclass
class ContourItem:
    name: str
    elev: float
    geom_ll: LineString
    geom_proj: LineString
    source: str = "original"


# =============================================================================
# KML/KMZ parsing
# =============================================================================

def extract_kml_from_upload(uploaded_file) -> str:
    raw = uploaded_file.getvalue()
    name = uploaded_file.name.lower()
    if name.endswith(".kmz"):
        with zipfile.ZipFile(io.BytesIO(raw), "r") as zf:
            kml_names = [n for n in zf.namelist() if n.lower().endswith(".kml")]
            if not kml_names:
                raise ValueError("El KMZ no contiene un archivo KML interno.")
            preferred = next((n for n in kml_names if n.lower().endswith("doc.kml")), kml_names[0])
            return zf.read(preferred).decode("utf-8", errors="ignore")
    if name.endswith(".kml"):
        return raw.decode("utf-8", errors="ignore")
    raise ValueError("El archivo debe ser KMZ o KML.")


def _coord_tokens_to_xy(text: str) -> List[Tuple[float, float]]:
    coords = []
    if not text:
        return coords
    raw = text.strip().replace("\n", " ").replace("\t", " ")
    for token in raw.split():
        parts = [p for p in token.split(",") if p != ""]
        if len(parts) >= 2:
            try:
                coords.append((float(parts[0]), float(parts[1])))
            except Exception:
                pass
    return coords


def parse_kml_items(kml_text: str) -> Dict[str, Any]:
    root = ET.fromstring(kml_text.encode("utf-8"))
    ns = {"k": "http://www.opengis.net/kml/2.2"}

    polygons_ll: List[Tuple[str, object]] = []
    points_ll: List[Tuple[str, float, float]] = []
    axes_ll: List[Tuple[str, LineString]] = []
    contours_ll: List[Tuple[str, float, LineString]] = []

    for pm in root.findall(".//k:Placemark", ns):
        name = pm.findtext("k:name", default="", namespaces=ns).strip()

        # Points
        for pnode in pm.findall(".//k:Point", ns):
            cnode = pnode.find(".//k:coordinates", ns)
            if cnode is not None and cnode.text:
                coords = _coord_tokens_to_xy(cnode.text)
                if coords:
                    lon, lat = coords[0]
                    points_ll.append((name or "Punto", lon, lat))

        # Polygons
        poly_geoms = []
        for pnode in pm.findall(".//k:Polygon", ns):
            outer = pnode.find(".//k:outerBoundaryIs/k:LinearRing/k:coordinates", ns)
            if outer is None or not outer.text:
                continue
            coords = _coord_tokens_to_xy(outer.text)
            if len(coords) >= 3:
                try:
                    poly = Polygon(coords)
                    if not poly.is_valid:
                        poly = poly.buffer(0)
                    if not poly.is_empty:
                        poly_geoms.append(poly)
                except Exception:
                    pass
        if poly_geoms:
            geom = unary_union(poly_geoms)
            polygons_ll.append((name or "Cuenca", geom))

        # Lines
        for lnode in pm.findall(".//k:LineString", ns):
            cnode = lnode.find(".//k:coordinates", ns)
            if cnode is None or not cnode.text:
                continue
            coords = _coord_tokens_to_xy(cnode.text)
            if len(coords) < 2:
                continue
            try:
                line = LineString(coords)
            except Exception:
                continue

            m = re.search(r"curva\s+(-?\d+(?:[.,]\d+)?)\s*m", name, re.I)
            if m:
                elev = float(m.group(1).replace(",", "."))
                contours_ll.append((name, elev, line))
            else:
                # Eje, cauce o cualquier línea no curva.
                lname = name.lower()
                if "eje" in lname or "cauce" in lname or "axis" in lname or "thalweg" in lname:
                    axes_ll.append((name or "Eje cauce", line))
                else:
                    # Respaldo: si no es curva pero no hay otro eje, puede seleccionarse luego.
                    axes_ll.append((name or "Línea", line))

    return {
        "polygons_ll": polygons_ll,
        "points_ll": points_ll,
        "axes_ll": axes_ll,
        "contours_ll": contours_ll,
    }


# =============================================================================
# CRS helpers
# =============================================================================

def choose_local_utm_epsg(lons: List[float], lats: List[float]) -> int:
    lon = float(np.nanmean(lons))
    lat = float(np.nanmean(lats))
    zone = int(math.floor((lon + 180.0) / 6.0) + 1)
    if lat >= 0:
        return 32600 + zone
    return 32700 + zone


def project_items(parsed: Dict[str, Any]) -> Dict[str, Any]:
    all_lons, all_lats = [], []
    for _, geom in parsed["polygons_ll"]:
        c = geom.centroid
        all_lons.append(c.x); all_lats.append(c.y)
    for _, lon, lat in parsed["points_ll"]:
        all_lons.append(lon); all_lats.append(lat)
    for _, line in parsed["axes_ll"][:3]:
        c = line.centroid
        all_lons.append(c.x); all_lats.append(c.y)

    if not all_lons:
        raise ValueError("No se encontraron coordenadas válidas en el KMZ.")
    epsg = choose_local_utm_epsg(all_lons, all_lats)
    crs_ll = PyCRS.from_epsg(4326)
    crs_proj = PyCRS.from_epsg(epsg)
    tr = Transformer.from_crs(crs_ll, crs_proj, always_xy=True)
    inv = Transformer.from_crs(crs_proj, crs_ll, always_xy=True)

    polygons = []
    for name, geom_ll in parsed["polygons_ll"]:
        geom_proj = shp_transform(tr.transform, geom_ll)
        polygons.append(PolygonItem(name, geom_ll, geom_proj, float(geom_proj.area / 1e6)))

    points = []
    for name, lon, lat in parsed["points_ll"]:
        x, y = tr.transform(lon, lat)
        points.append(PointItem(name, lon, lat, Point(x, y)))

    axes = []
    for name, geom_ll in parsed["axes_ll"]:
        geom_proj = shp_transform(tr.transform, geom_ll)
        if geom_proj.length > 0:
            axes.append(AxisItem(name, geom_ll, geom_proj, float(geom_proj.length / 1000.0)))

    contours = []
    for name, elev, geom_ll in parsed["contours_ll"]:
        geom_proj = shp_transform(tr.transform, geom_ll)
        if geom_proj.length > 0:
            contours.append(ContourItem(name, elev, geom_ll, geom_proj, "original"))

    return {
        "epsg": epsg,
        "crs_proj": crs_proj,
        "to_proj": tr,
        "to_ll": inv,
        "polygons": polygons,
        "points": points,
        "axes": axes,
        "contours": contours,
    }


# =============================================================================
# Geometry logic
# =============================================================================

def auto_roles(polygons: List[PolygonItem]) -> Tuple[Optional[int], Optional[int]]:
    if len(polygons) < 2:
        return (0 if polygons else None, None)
    order = sorted(range(len(polygons)), key=lambda i: polygons[i].area_km2, reverse=True)
    return order[0], order[1]


def nearest_point_index(points: List[PointItem], target_geom) -> Optional[int]:
    if not points:
        return None
    dists = [p.geom_proj.distance(target_geom) for p in points]
    return int(np.argmin(dists))


def select_active_axis(axis: LineString, p1: Point, p2: Point) -> Tuple[LineString, float, float]:
    d1 = axis.project(p1)
    d2 = axis.project(p2)
    a, b = (d1, d2) if d1 <= d2 else (d2, d1)
    if shp_substring is not None:
        try:
            sub = shp_substring(axis, a, b)
            if isinstance(sub, LineString) and sub.length > 0:
                return sub, float(a), float(b)
        except Exception:
            pass

    # Fallback: sample points along interval.
    n = max(2, int((b - a) / 10.0) + 1)
    coords = [axis.interpolate(float(d)).coords[0] for d in np.linspace(a, b, n)]
    return LineString(coords), float(a), float(b)


def clean_lines_from_intersection(geom) -> List[LineString]:
    lines = []
    if geom is None or geom.is_empty:
        return lines
    if isinstance(geom, LineString):
        if geom.length > 0 and len(geom.coords) >= 2:
            lines.append(geom)
    elif isinstance(geom, MultiLineString):
        for g in geom.geoms:
            lines.extend(clean_lines_from_intersection(g))
    elif isinstance(geom, GeometryCollection):
        for g in geom.geoms:
            lines.extend(clean_lines_from_intersection(g))
    return lines


def geom_to_kml_coords(geom, to_ll: Transformer) -> str:
    geom_ll = shp_transform(to_ll.transform, geom)
    return " ".join([f"{x:.8f},{y:.8f},0" for x, y in geom_ll.coords])


def polygon_outer_kml(poly, to_ll: Transformer) -> str:
    poly_ll = shp_transform(to_ll.transform, poly)
    if isinstance(poly_ll, MultiPolygon):
        poly_ll = max(poly_ll.geoms, key=lambda g: g.area)
    return " ".join([f"{x:.8f},{y:.8f},0" for x, y in poly_ll.exterior.coords])


def sample_axis_stations(axis: LineString, step_m: float) -> np.ndarray:
    if axis.length <= 0:
        return np.array([0.0])
    n = max(2, int(math.floor(axis.length / max(step_m, 1.0))) + 1)
    stations = np.linspace(0.0, axis.length, n)
    if stations[-1] < axis.length:
        stations = np.r_[stations, axis.length]
    return stations


def normal_transect(axis: LineString, dist: float, half_width_m: float) -> Optional[LineString]:
    if axis.length <= 0:
        return None
    eps = max(1.0, min(10.0, axis.length / 50.0))
    d0 = max(0.0, dist - eps)
    d1 = min(axis.length, dist + eps)
    p = axis.interpolate(dist)
    p0 = axis.interpolate(d0)
    p1 = axis.interpolate(d1)
    dx = p1.x - p0.x
    dy = p1.y - p0.y
    norm = math.hypot(dx, dy)
    if norm <= 0:
        return None
    # Normal izquierda/derecha
    nx = -dy / norm
    ny = dx / norm
    a = (p.x - nx * half_width_m, p.y - ny * half_width_m)
    b = (p.x + nx * half_width_m, p.y + ny * half_width_m)
    return LineString([a, b])


def point_like_from_intersection(geom) -> List[Point]:
    pts = []
    if geom is None or geom.is_empty:
        return pts
    if isinstance(geom, Point):
        pts.append(geom)
    elif geom.geom_type == "MultiPoint":
        pts.extend(list(geom.geoms))
    elif isinstance(geom, LineString):
        # Si una curva coincide parcialmente con transecta, tomar punto medio.
        pts.append(geom.interpolate(0.5, normalized=True))
    elif isinstance(geom, GeometryCollection):
        for g in geom.geoms:
            pts.extend(point_like_from_intersection(g))
    return pts


def generate_auxiliary_contours(
    active_axis: LineString,
    contours: List[ContourItem],
    corridor: object,
    interval_m: float = 1.0,
    station_step_m: float = 25.0,
    half_width_m: float = 50.0,
    min_points: int = 4,
    simplify_m: float = 0.5,
    max_output_lines: int = 5000,
) -> Tuple[List[ContourItem], List[ContourItem], pd.DataFrame]:
    """Genera curvas auxiliares longitudinales a partir de intersecciones de curvas existentes con transectas."""
    clipped: List[ContourItem] = []
    for c in contours:
        if not c.geom_proj.intersects(corridor):
            continue
        inter = c.geom_proj.intersection(corridor)
        for line in clean_lines_from_intersection(inter):
            clipped.append(ContourItem(c.name, c.elev, c.geom_ll, line, "original_recortada"))

    if not clipped:
        return [], [], pd.DataFrame()

    min_e = math.floor(min(c.elev for c in clipped) / interval_m) * interval_m
    max_e = math.ceil(max(c.elev for c in clipped) / interval_m) * interval_m
    targets = np.arange(min_e, max_e + interval_m * 0.1, interval_m)
    targets = [float(t) for t in targets]

    stations = sample_axis_stations(active_axis, station_step_m)

    # Build optional tree
    geoms = [c.geom_proj for c in clipped]
    tree = STRtree(geoms) if STRtree is not None and geoms else None
    geom_to_items: Dict[int, List[ContourItem]] = {}
    for c in clipped:
        geom_to_items.setdefault(id(c.geom_proj), []).append(c)

    # Dict: (target elevation, side) -> list of (station, point)
    acc: Dict[Tuple[float, str], List[Tuple[float, Point]]] = {}
    station_rows = []

    for s in stations:
        tran = normal_transect(active_axis, float(s), half_width_m)
        if tran is None:
            continue
        candidates = []
        if tree is not None:
            try:
                for g in tree.query(tran):
                    candidates.extend(geom_to_items.get(id(g), []))
            except Exception:
                candidates = clipped
        else:
            candidates = clipped

        crossings: List[Tuple[float, float, Point]] = []
        for c in candidates:
            if not c.geom_proj.intersects(tran):
                continue
            inter = c.geom_proj.intersection(tran)
            for pt in point_like_from_intersection(inter):
                off = tran.project(pt) - half_width_m
                if -half_width_m - 1e-6 <= off <= half_width_m + 1e-6:
                    crossings.append((float(off), float(c.elev), pt))

        # Remove rough duplicates
        crossings.sort(key=lambda x: (x[0], x[1]))
        uniq = []
        for off, elev, pt in crossings:
            if not uniq or abs(off - uniq[-1][0]) > 0.75 or abs(elev - uniq[-1][1]) > 0.01:
                uniq.append((off, elev, pt))
        crossings = sorted(uniq, key=lambda x: x[0])

        if len(crossings) < 2:
            station_rows.append({"station_m": float(s), "n_crossings": len(crossings), "n_interpolated": 0})
            continue

        n_interp = 0
        # Interpolate between consecutive known crossings
        for (off1, e1, _), (off2, e2, _) in zip(crossings[:-1], crossings[1:]):
            if abs(e2 - e1) < 1e-9:
                continue
            lo, hi = sorted([e1, e2])
            for t in targets:
                if t <= lo or t >= hi:
                    continue
                f = (t - e1) / (e2 - e1)
                off = off1 + f * (off2 - off1)
                pt = tran.interpolate(off + half_width_m)
                side = "izquierda" if off < 0 else "derecha"
                acc.setdefault((round(t, 6), side), []).append((float(s), pt))
                n_interp += 1
        station_rows.append({"station_m": float(s), "n_crossings": len(crossings), "n_interpolated": n_interp})

    aux: List[ContourItem] = []
    for (elev, side), pts in sorted(acc.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if len(aux) >= max_output_lines:
            break
        pts = sorted(pts, key=lambda x: x[0])
        # Split by station gaps
        chunks = []
        cur = []
        last_s = None
        max_gap = station_step_m * 2.5
        for s, pt in pts:
            if last_s is not None and (s - last_s) > max_gap:
                if len(cur) >= min_points:
                    chunks.append(cur)
                cur = []
            cur.append((s, pt))
            last_s = s
        if len(cur) >= min_points:
            chunks.append(cur)

        for chunk in chunks:
            if len(aux) >= max_output_lines:
                break
            line = LineString([pt.coords[0] for _, pt in chunk])
            if simplify_m > 0:
                line = line.simplify(simplify_m, preserve_topology=False)
            if line.length <= 0 or len(line.coords) < min_points:
                continue
            aux.append(ContourItem(f"Curva auxiliar {elev:g} m · {side}", float(elev), line, line, "auxiliar_interpolada"))

    qc_df = pd.DataFrame(station_rows)
    return clipped, aux, qc_df


def build_distributed_q_table(
    axis_len_km: float,
    area_inc_km2: float,
    q_hidro: float,
    q_inc: float,
    step_km: float,
) -> pd.DataFrame:
    if axis_len_km <= 0:
        return pd.DataFrame()
    n = max(2, int(math.floor(axis_len_km / max(step_km, 0.001))) + 1)
    kms = np.linspace(0, axis_len_km, n)
    if kms[-1] < axis_len_km:
        kms = np.r_[kms, axis_len_km]
    rows = []
    for km in kms:
        frac = min(1.0, max(0.0, km / axis_len_km))
        area_acum = area_inc_km2 * frac
        q_acum = q_inc * frac
        rows.append({
            "Km desde PC-HIDRO hacia PC-DESCARGA": km,
            "% longitud activa": frac * 100.0,
            "Área incremental acumulada (km²)": area_acum,
            "Q incremental acumulado (m³/s)": q_acum,
            "Q total en tramo (m³/s)": q_hidro + q_acum,
        })
    return pd.DataFrame(rows)


# =============================================================================
# Export
# =============================================================================

def kml_escape(v: Any) -> str:
    return escape(str(v), quote=True)


def style_block() -> str:
    return """
    <Style id="cuenca_desc"><LineStyle><color>ffff0000</color><width>3</width></LineStyle><PolyStyle><color>330000ff</color></PolyStyle></Style>
    <Style id="cuenca_hidro"><LineStyle><color>ff00ffff</color><width>3</width></LineStyle><PolyStyle><color>3300ffff</color></PolyStyle></Style>
    <Style id="intercuenca"><LineStyle><color>ff00aa00</color><width>3</width></LineStyle><PolyStyle><color>3300aa00</color></PolyStyle></Style>
    <Style id="corredor"><LineStyle><color>ff9900ff</color><width>2</width></LineStyle><PolyStyle><color>209900ff</color></PolyStyle></Style>
    <Style id="eje"><LineStyle><color>ff0000ff</color><width>4</width></LineStyle></Style>
    <Style id="curva_original"><LineStyle><color>ff996633</color><width>1</width></LineStyle></Style>
    <Style id="curva_aux"><LineStyle><color>ff00ff00</color><width>1.4</width></LineStyle></Style>
    <Style id="pc_h"><IconStyle><color>ff00ffff</color><scale>1.1</scale><Icon><href>http://maps.google.com/mapfiles/kml/paddle/ylw-circle.png</href></Icon></IconStyle></Style>
    <Style id="pc_d"><IconStyle><color>ffff0000</color><scale>1.1</scale><Icon><href>http://maps.google.com/mapfiles/kml/paddle/blu-circle.png</href></Icon></IconStyle></Style>
    """


def placemark_polygon(name: str, geom, style: str, to_ll: Transformer, desc: str = "") -> str:
    if geom is None or geom.is_empty:
        return ""
    polys = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
    out = []
    for i, p in enumerate(polys, start=1):
        if p.is_empty or not hasattr(p, "exterior"):
            continue
        coords = polygon_outer_kml(p, to_ll)
        out.append(f"""
        <Placemark>
          <name>{kml_escape(name)} {i}</name>
          <description>{kml_escape(desc)}</description>
          <styleUrl>#{style}</styleUrl>
          <Polygon><outerBoundaryIs><LinearRing><coordinates>{coords}</coordinates></LinearRing></outerBoundaryIs></Polygon>
        </Placemark>""")
    return "\n".join(out)


def placemark_line(name: str, line: LineString, style: str, to_ll: Transformer, desc: str = "") -> str:
    if line is None or line.is_empty or line.length <= 0:
        return ""
    coords = geom_to_kml_coords(line, to_ll)
    return f"""
        <Placemark>
          <name>{kml_escape(name)}</name>
          <description>{kml_escape(desc)}</description>
          <styleUrl>#{style}</styleUrl>
          <LineString><tessellate>1</tessellate><coordinates>{coords}</coordinates></LineString>
        </Placemark>"""


def placemark_point(name: str, pt: PointItem, style: str) -> str:
    return f"""
        <Placemark>
          <name>{kml_escape(name)}</name>
          <styleUrl>#{style}</styleUrl>
          <Point><coordinates>{pt.lon:.8f},{pt.lat:.8f},0</coordinates></Point>
        </Placemark>"""


def make_kmz(
    desc_poly: PolygonItem,
    hidro_poly: PolygonItem,
    inter_geom,
    axis_active: LineString,
    corridor_geom,
    pc_h: Optional[PointItem],
    pc_d: Optional[PointItem],
    clipped_contours: List[ContourItem],
    aux_contours: List[ContourItem],
    to_ll: Transformer,
    include_original_clipped: bool,
) -> bytes:
    placemarks = []
    placemarks.append(placemark_polygon("Cuenca descarga / soporte", desc_poly.geom_proj, "cuenca_desc", to_ll, f"Área {desc_poly.area_km2:.3f} km²"))
    placemarks.append(placemark_polygon("Cuenca hidrológica", hidro_poly.geom_proj, "cuenca_hidro", to_ll, f"Área {hidro_poly.area_km2:.3f} km²"))
    placemarks.append(placemark_polygon("Intercuenca incremental", inter_geom, "intercuenca", to_ll, "Área incremental = descarga - hidrológica"))
    placemarks.append(placemark_polygon("Corredor hidráulico del eje activo", corridor_geom, "corredor", to_ll, "Área usada para recortar/incrementar curvas"))
    placemarks.append(placemark_line("Eje hidráulico activo PC-HIDRO a PC-DESCARGA", axis_active, "eje", to_ll, "Tramo activo usado para aporte incremental y curvas"))
    if pc_h:
        placemarks.append(placemark_point("PC-HIDRO", pc_h, "pc_h"))
    if pc_d:
        placemarks.append(placemark_point("PC-DESCARGA", pc_d, "pc_d"))

    if include_original_clipped:
        for c in clipped_contours:
            placemarks.append(placemark_line(f"Original recortada {c.elev:g} m", c.geom_proj, "curva_original", to_ll, "Curva original del KMZ consolidado recortada al corredor"))

    for c in aux_contours:
        placemarks.append(placemark_line(c.name, c.geom_proj, "curva_aux", to_ll, "Curva auxiliar interpolada; no representa topografía levantada"))

    kml = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>HidroSed Cuenca Consolidada + Eje Activo</name>
    {style_block()}
    {''.join(placemarks)}
  </Document>
</kml>"""
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("doc.kml", kml.encode("utf-8"))
    return bio.getvalue()


def make_excel_bytes(summary_df, q_df, qc_df, contour_summary_df) -> bytes:
    bio = io.BytesIO()
    with pd.ExcelWriter(bio, engine="openpyxl") as writer:
        summary_df.to_excel(writer, index=False, sheet_name="Resumen")
        q_df.to_excel(writer, index=False, sheet_name="Q_distribuido")
        contour_summary_df.to_excel(writer, index=False, sheet_name="Curvas")
        if qc_df is not None and not qc_df.empty:
            qc_df.to_excel(writer, index=False, sheet_name="QC_transectas")
    return bio.getvalue()


def make_result_zip(kmz, excel, summary_json, q_df, contour_summary_df) -> bytes:
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("01_cuenca_consolidada_eje_activo.kmz", kmz)
        zf.writestr("02_resultados.xlsx", excel)
        zf.writestr("03_q_distribuido.csv", q_df.to_csv(index=False).encode("utf-8-sig"))
        zf.writestr("04_resumen_curvas.csv", contour_summary_df.to_csv(index=False).encode("utf-8-sig"))
        zf.writestr("05_resumen_integracion_hidrosed.json", json.dumps(summary_json, ensure_ascii=False, indent=2).encode("utf-8"))
    return bio.getvalue()


# =============================================================================
# UI
# =============================================================================

st.title("🌊 HidroSed · Cuenca Consolidada + Eje Activo v1")
st.caption("Lee una cuenca consolidada, asocia intercuenca al eje del cauce e incrementa curvas solo en el corredor hidráulico.")

with st.expander("Sensibilización: ¿qué significan buffer y bbox?", expanded=True):
    st.markdown(
        """
        **Buffer** es una franja de trabajo alrededor de una línea o punto.  
        En la práctica, si el eje del cauce es el centro, un buffer de **50 m** significa:  
        **tomar 50 m hacia la izquierda y 50 m hacia la derecha del cauce**.  
        Esa franja se usa para recortar curvas, revisar riberas y preparar secciones.

        **BBox** viene de *bounding box* o “caja envolvente”.  
        En la práctica es el **rectángulo mínimo** que encierra una zona.  
        Si una app descarga DEM usando un bbox pequeño, es como cortar un mapa demasiado chico:  
        la cuenca puede quedar incompleta.  

        En este módulo se evita depender del bbox para redelimitar cuencas: se parte desde una **cuenca consolidada ya validada**.
        """
    )

left, right = st.columns([0.9, 1.1])

with left:
    st.subheader("1. Cargar KMZ consolidado")
    file = st.file_uploader("Cuenca consolidada KMZ/KML", type=["kmz", "kml"])
    st.caption("Debe contener al menos dos polígonos de cuenca, un eje de cauce y curvas de nivel.")

    st.subheader("2. Parámetros del eje y corredor")
    corridor_m = st.number_input("Buffer/corredor a cada lado del eje (m)", min_value=5.0, max_value=1000.0, value=80.0, step=5.0)
    station_step_m = st.number_input("Separación de transectas de análisis (m)", min_value=5.0, max_value=200.0, value=25.0, step=5.0)
    contour_interval_m = st.number_input("Equidistancia de curvas auxiliares (m)", min_value=0.5, max_value=20.0, value=1.0, step=0.5)
    simplify_m = st.number_input("Simplificación de curvas auxiliares (m)", min_value=0.0, max_value=20.0, value=1.0, step=0.5)
    min_pts = st.number_input("Mínimo de puntos para aceptar una curva auxiliar", min_value=2, max_value=20, value=4, step=1)
    include_original = st.checkbox("Incluir curvas originales recortadas en KMZ", value=True)

with right:
    st.subheader("3. Aporte / caudal incremental")
    q_hidro = st.number_input("Q en PC-HIDRO (m³/s)", min_value=0.0, value=0.0, step=0.1)
    q_mode = st.radio("Método para Q incremental", ["Sin Q incremental", "Manual", "Proporcional por área"], index=0)
    q_inc_manual = 0.0
    if q_mode == "Manual":
        q_inc_manual = st.number_input("Q incremental manual (m³/s)", min_value=0.0, value=0.0, step=0.1)
    q_specific = 0.0
    if q_mode == "Proporcional por área":
        q_specific = st.number_input("Caudal específico incremental (m³/s/km²)", min_value=0.0, value=0.0, step=0.01)
    q_step_km = st.number_input("Paso de tabla Q distribuido (km)", min_value=0.05, max_value=5.0, value=0.5, step=0.05)

run = st.button("Procesar cuenca consolidada", type="primary", use_container_width=True)

if run:
    if file is None:
        st.error("Debe cargar un KMZ/KML consolidado.")
        st.stop()

    try:
        with st.status("Procesando KMZ consolidado...", expanded=True) as status:
            status.write("Leyendo geometrías...")
            kml_text = extract_kml_from_upload(file)
            parsed = parse_kml_items(kml_text)
            data = project_items(parsed)

            polygons: List[PolygonItem] = data["polygons"]
            points: List[PointItem] = data["points"]
            axes: List[AxisItem] = data["axes"]
            contours: List[ContourItem] = data["contours"]

            if len(polygons) < 2:
                raise ValueError("El KMZ debe contener al menos dos polígonos: cuenca descarga y cuenca hidrológica.")
            if not axes:
                raise ValueError("El KMZ debe contener un eje del cauce como LineString.")
            if not contours:
                raise ValueError("El KMZ no contiene curvas de nivel reconocibles con nombres tipo 'Curva 340 m'.")

            idx_desc, idx_hidro = auto_roles(polygons)
            desc_poly = polygons[idx_desc]
            hidro_poly = polygons[idx_hidro]

            # Eje: escoger el más largo por defecto
            axis = max(axes, key=lambda a: a.length_km)

            # PC: escoger puntos más cercanos a los cierres de cuenca. Priorizamos puntos cerca del borde de cada polígono.
            pc_d_idx = nearest_point_index(points, desc_poly.geom_proj.boundary)
            pc_h_idx = nearest_point_index(points, hidro_poly.geom_proj.boundary)
            pc_d = points[pc_d_idx] if pc_d_idx is not None else None
            pc_h = points[pc_h_idx] if pc_h_idx is not None else None

            if pc_h is None or pc_d is None:
                raise ValueError("No se pudieron identificar puntos de control en el KMZ.")

            status.write("Recortando eje activo entre PC-HIDRO y PC-DESCARGA...")
            active_axis, d0, d1 = select_active_axis(axis.geom_proj, pc_h.geom_proj, pc_d.geom_proj)
            if active_axis.length <= 0:
                raise ValueError("No se pudo construir el eje activo entre los puntos.")

            status.write("Calculando intercuenca...")
            inter_geom = desc_poly.geom_proj.difference(hidro_poly.geom_proj)
            inter_area_km2 = float(inter_geom.area / 1e6) if not inter_geom.is_empty else 0.0
            area_simple = desc_poly.area_km2 - hidro_poly.area_km2
            pct_inside = float(desc_poly.geom_proj.intersection(hidro_poly.geom_proj).area / hidro_poly.geom_proj.area * 100.0)

            corridor = active_axis.buffer(float(corridor_m), cap_style=2, join_style=2)
            corridor_area_km2 = float(corridor.area / 1e6)

            status.write("Recortando e incrementando curvas dentro del corredor del eje...")
            clipped, aux, qc_df = generate_auxiliary_contours(
                active_axis,
                contours,
                corridor,
                interval_m=float(contour_interval_m),
                station_step_m=float(station_step_m),
                half_width_m=float(corridor_m),
                min_points=int(min_pts),
                simplify_m=float(simplify_m),
            )

            if q_mode == "Sin Q incremental":
                q_inc = 0.0
            elif q_mode == "Manual":
                q_inc = float(q_inc_manual)
            else:
                q_inc = float(q_specific) * max(inter_area_km2, 0.0)

            q_df = build_distributed_q_table(
                axis_len_km=active_axis.length / 1000.0,
                area_inc_km2=inter_area_km2,
                q_hidro=float(q_hidro),
                q_inc=float(q_inc),
                step_km=float(q_step_km),
            )

            contour_summary_df = pd.DataFrame([
                {"Tipo": "Curvas originales totales en KMZ", "Cantidad": len(contours), "Observación": "Curvas leídas desde archivo consolidado"},
                {"Tipo": "Curvas originales recortadas al corredor", "Cantidad": len(clipped), "Observación": "Curvas reales del KMZ dentro del corredor"},
                {"Tipo": "Curvas auxiliares interpoladas", "Cantidad": len(aux), "Observación": "Curvas generadas por interpolación entre curvas existentes"},
                {"Tipo": "Equidistancia original estimada", "Cantidad": None, "Observación": "Depende del KMZ; usualmente 20 m en el archivo de respaldo"},
                {"Tipo": "Equidistancia auxiliar solicitada", "Cantidad": float(contour_interval_m), "Observación": "m"},
            ])

            area_per_km = inter_area_km2 / (active_axis.length / 1000.0) if active_axis.length > 0 else float("nan")
            summary_df = pd.DataFrame([
                {"Parámetro": "CRS proyectado usado", "Valor": f"EPSG:{data['epsg']}"},
                {"Parámetro": "Cuenca descarga / soporte", "Valor": desc_poly.name},
                {"Parámetro": "Área cuenca descarga (km²)", "Valor": desc_poly.area_km2},
                {"Parámetro": "Cuenca hidrológica", "Valor": hidro_poly.name},
                {"Parámetro": "Área cuenca hidrológica (km²)", "Valor": hidro_poly.area_km2},
                {"Parámetro": "Área incremental geométrica (km²)", "Valor": inter_area_km2},
                {"Parámetro": "Diferencia simple de áreas (km²)", "Valor": area_simple},
                {"Parámetro": "% cuenca hidrológica dentro de descarga", "Valor": pct_inside},
                {"Parámetro": "Eje usado", "Valor": axis.name},
                {"Parámetro": "Longitud eje activo (km)", "Valor": active_axis.length / 1000.0},
                {"Parámetro": "Área incremental por km de eje (km²/km)", "Valor": area_per_km},
                {"Parámetro": "Buffer/corredor a cada lado del eje (m)", "Valor": float(corridor_m)},
                {"Parámetro": "Área corredor hidráulico (km²)", "Valor": corridor_area_km2},
                {"Parámetro": "Q PC-HIDRO (m³/s)", "Valor": float(q_hidro)},
                {"Parámetro": "Q incremental (m³/s)", "Valor": q_inc},
                {"Parámetro": "Q total en PC-DESCARGA (m³/s)", "Valor": float(q_hidro) + q_inc},
            ])

            summary_json = {
                "generated_at": datetime.utcnow().isoformat() + "Z",
                "app": APP_TITLE,
                "epsg": data["epsg"],
                "area_descarga_km2": desc_poly.area_km2,
                "area_hidrologica_km2": hidro_poly.area_km2,
                "area_incremental_km2": inter_area_km2,
                "longitud_eje_activo_km": active_axis.length / 1000.0,
                "area_incremental_por_km": area_per_km,
                "q_hidro_m3s": float(q_hidro),
                "q_incremental_m3s": q_inc,
                "q_total_descarga_m3s": float(q_hidro) + q_inc,
                "buffer_corredor_m": float(corridor_m),
                "curvas_originales_recortadas": len(clipped),
                "curvas_auxiliares": len(aux),
                "nota": "Curvas auxiliares interpoladas; no son topografía levantada.",
            }

            kmz_bytes = make_kmz(
                desc_poly, hidro_poly, inter_geom, active_axis, corridor,
                pc_h, pc_d, clipped, aux, data["to_ll"], bool(include_original)
            )
            excel_bytes = make_excel_bytes(summary_df, q_df, qc_df, contour_summary_df)
            zip_bytes = make_result_zip(kmz_bytes, excel_bytes, summary_json, q_df, contour_summary_df)

            st.session_state["result"] = {
                "summary_df": summary_df,
                "q_df": q_df,
                "qc_df": qc_df,
                "contour_summary_df": contour_summary_df,
                "kmz_bytes": kmz_bytes,
                "excel_bytes": excel_bytes,
                "zip_bytes": zip_bytes,
                "desc_poly": desc_poly,
                "hidro_poly": hidro_poly,
                "inter_geom": inter_geom,
                "active_axis": active_axis,
                "corridor": corridor,
                "clipped_n": len(clipped),
                "aux_n": len(aux),
                "epsg": data["epsg"],
                "q_inc": q_inc,
            }
            status.update(label="Proceso terminado", state="complete")

    except Exception as exc:
        st.exception(exc)
        st.stop()

if "result" in st.session_state:
    r = st.session_state["result"]
    st.divider()
    st.header("Resultados")

    sd = r["summary_df"].set_index("Parámetro")["Valor"].to_dict()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Área descarga", f"{float(sd['Área cuenca descarga (km²)']):.3f} km²")
    c2.metric("Área hidrológica", f"{float(sd['Área cuenca hidrológica (km²)']):.3f} km²")
    c3.metric("Intercuenca", f"{float(sd['Área incremental geométrica (km²)']):.3f} km²")
    c4.metric("Eje activo", f"{float(sd['Longitud eje activo (km)']):.3f} km")

    c5, c6, c7 = st.columns(3)
    c5.metric("Curvas originales recortadas", f"{r['clipped_n']}")
    c6.metric("Curvas auxiliares", f"{r['aux_n']}")
    c7.metric("Q incremental", f"{r['q_inc']:.3f} m³/s")

    st.subheader("Resumen técnico")
    st.dataframe(r["summary_df"], use_container_width=True)

    st.subheader("Caudal distribuido a lo largo del eje")
    st.dataframe(r["q_df"], use_container_width=True, height=300)

    st.subheader("Resumen de curvas")
    st.dataframe(r["contour_summary_df"], use_container_width=True)

    with st.expander("Control de calidad por transecta", expanded=False):
        st.dataframe(r["qc_df"], use_container_width=True, height=300)

    st.subheader("Vista técnica simplificada")
    try:
        fig, ax = plt.subplots(figsize=(10, 8))
        for geom, label, lw in [
            (r["desc_poly"].geom_proj, "Cuenca descarga", 2.5),
            (r["hidro_poly"].geom_proj, "Cuenca hidrológica", 2.5),
            (r["inter_geom"], "Intercuenca", 1.5),
            (r["corridor"], "Corredor eje", 1.0),
        ]:
            geoms = list(geom.geoms) if isinstance(geom, MultiPolygon) else [geom]
            first = True
            for g in geoms:
                if hasattr(g, "exterior"):
                    x, y = g.exterior.xy
                    ax.plot(x, y, linewidth=lw, label=label if first else None)
                    first = False
        x, y = r["active_axis"].xy
        ax.plot(x, y, linewidth=3, label="Eje activo")
        ax.set_aspect("equal", adjustable="box")
        ax.set_title("Cuenca consolidada + eje activo + corredor")
        ax.legend()
        st.pyplot(fig)
    except Exception as exc:
        st.warning(f"No se pudo graficar: {exc}")

    st.subheader("Descargas")
    d1, d2, d3 = st.columns(3)
    d1.download_button(
        "ZIP completo",
        r["zip_bytes"],
        file_name="hidrosed_cuenca_consolidada_eje_activo_v1_resultados.zip",
        mime="application/zip",
        use_container_width=True,
    )
    d2.download_button(
        "KMZ preparado",
        r["kmz_bytes"],
        file_name="cuenca_consolidada_eje_activo_curvas_incrementadas.kmz",
        mime="application/vnd.google-earth.kmz",
        use_container_width=True,
    )
    d3.download_button(
        "Excel resultados",
        r["excel_bytes"],
        file_name="cuenca_consolidada_eje_activo_resultados.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        use_container_width=True,
    )

st.divider()
st.markdown(
    """
    **Alcance v1:** usa una cuenca consolidada validada, no redelimita cuencas desde cero.  
    **Advertencia:** las curvas auxiliares interpoladas son apoyo para preparación hidráulica; no reemplazan levantamiento topográfico, DEM de alta resolución ni secciones de terreno.
    """
)
