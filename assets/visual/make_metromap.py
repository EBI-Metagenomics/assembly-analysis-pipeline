#!/usr/bin/env python3
"""Render a pipeline metro map (nf-core workflow-schematic style) from a TOML spec.

The spec gives the lines, the lanes, the stages and the stations. The script
calculates all the positions: the columns from the label widths, the lanes from
the space that the labels need, and the curves between the lanes.

The style is the style of the MGnify pipeline maps (miassembler, metaViraVerse):

- one coloured line for each input or sequence type, in parallel tracks;
- grey boxes for the stages;
- the tool name (bold) above each station, a short description below;
- a black dot for a mandatory step, an empty station for an optional step;
- file icons, input triangles, output bars, a title and a legend.

This script comes from the mi-metro-map skill of EBI-Metagenomics/mi-agentic
(``skills/mi-metro-map/scripts/metromap.py``). The spec format is in
``skills/mi-metro-map/references/spec.md`` of that repository. Change the spec, not
this copy: make script changes in mi-agentic and copy the script again.

The script needs Python 3.11 or later and no other packages.

Usage::

    make_metromap.py SPEC.toml [--out MAP.svg] [--font MavenPro.ttf] [--png MAP.png]
"""

import argparse
import base64
import math
import shutil
import subprocess
import sys
import tomllib
from dataclasses import dataclass, field
from itertools import pairwise
from pathlib import Path

FONT = "'Maven Pro', 'Helvetica Neue', Arial, sans-serif"

# Lines and stations
GAP = 14  # distance between parallel tracks
LW = 9  # line width
R = 13  # station radius
DOT = 5.5  # radius of the black dot of a mandatory station

# Labels: font size, line height and colour of each kind of text
SIZE = {"name": 28, "desc": 26, "note": 22}
LH = {"name": 32, "desc": 30, "note": 28}
FILL = {"name": "#000", "desc": "#000", "note": "#555"}
NAME_DY = 30  # last tool-name baseline, above the highest track in use at the station
DESC_DY = 44  # first description baseline, below the lowest track in use at the station
FORK_DY = 45  # branch offset of a fork with two branches
FORK_PITCH = 80  # minimum branch distance of a fork with three or more branches

# Horizontal layout
COL_GAP = 24  # minimum space between the labels of two columns
MIN_PITCH = 150  # minimum distance between two columns (room for a curve)
CURVE = 90  # width of the curve from one lane to another
STAGE_PAD, STAGE_GAP = 20, 15
TERMINUS_GAP = 20
MARGIN = 16

# Vertical layout
LANE_GAP = 120  # minimum distance between the last track of a lane and the next lane
CLEAR = 16  # minimum space between the labels and tracks of different lanes
TITLE_LH = 32

# File icons: (width, height) of the map inputs and outputs, and of the station inputs
ICON = (62 * 1.25, 80 * 1.25)
SMALL_ICON = (62, 80)


class SpecError(Exception):
    """The spec has errors. The message lists all of them."""


@dataclass
class Station:
    names: list  # tool name lines (stacked above the station)
    branches: list  # for a fork: [(name lines, desc lines)], else []
    desc: list
    note: list
    lane: str
    stage: str
    col: int
    lines: list
    optional: bool
    input: tuple | None = None  # (label, caption lines, note lines)
    x: float = 0.0
    top: float = 0.0  # highest track in use at the station (set by the layout)
    bottom: float = 0.0  # lowest track in use at the station (set by the layout)
    icon_y: float = 0.0  # top of the station input icon (set by the layout)

    @property
    def fork(self):
        return bool(self.branches)


@dataclass
class Model:
    title: dict
    style: dict
    legend: dict
    lines: dict  # id -> {"label", "color", "start", "end"}, in legend order
    lanes: dict  # id -> list of tracks (each a list of line ids), lanes from top to bottom
    stages: list  # [(id, title lines)], from left to right
    stations: list
    inputs: list  # [(label, caption lines, lane or None)]
    outputs: list
    weight: int = field(init=False)

    def __post_init__(self):
        self.weight = self.style["tool_weight"]


def text_lines(value):
    """Split a label into lines. A label is a string with ``\\n`` or a list."""
    if not value:
        return []
    if isinstance(value, list):
        return [str(v) for v in value]
    return str(value).split("\n")


def text_width(s, size, weight=400):
    """Estimate the width of a text in Maven Pro (no font metrics needed)."""
    return len(s) * size * (0.53 if weight < 700 else 0.6 if weight < 800 else 0.63)


def lane_lines(tracks):
    return [ln for t in tracks for ln in t]


def parse(spec):
    """Check the spec and make the model.

    :param spec: The spec, as a dict from ``tomllib``.
    :return: A :class:`Model`.
    :raises SpecError: If the spec has errors.
    """
    errors = []
    lines = {}
    for ln in spec.get("lines", []):
        if not ln.get("id") or not ln.get("color"):
            errors.append(f"line {ln}: needs 'id' and 'color'")
            continue
        if ln["id"] in lines:
            errors.append(f"line '{ln['id']}': duplicate id")
        start, end = ln.get("start", "input"), ln.get("end", "output")
        if start not in ("input", "station") or end not in ("output", "station"):
            errors.append(f"line '{ln['id']}': 'start' must be input or station, 'end' must be output or station")
        lines[ln["id"]] = {"label": ln.get("label", ln["id"]), "color": ln["color"], "start": start, "end": end}
    if not lines:
        errors.append("the spec has no [[lines]]")

    lanes = {}
    for lane in spec.get("lanes", []):
        lid = lane.get("id")
        if not lid or lid in lanes:
            errors.append(f"lane {lane}: needs a unique 'id'")
            continue
        tracks = [t if isinstance(t, list) else [t] for t in lane.get("lines", [])]
        flat = lane_lines(tracks)
        unknown = [x for x in flat if x not in lines]
        if unknown or not flat:
            errors.append(f"lane '{lid}': unknown or no lines {unknown}")
        elif len(flat) != len(set(flat)):
            errors.append(f"lane '{lid}': a line is on two tracks")
        lanes[lid] = tracks
    if not lanes:
        errors.append("the spec has no [[lanes]]")

    style = {"numbered_stages": False, "tool_weight": 700}
    style.update(spec.get("style", {}))
    stages = []
    for st in spec.get("stages", []):
        if not st.get("id") or st["id"] in [s[0] for s in stages]:
            errors.append(f"stage {st}: needs a unique 'id'")
            continue
        title = text_lines(st.get("title", st["id"]))
        if style["numbered_stages"]:
            title[0] = f"{len(stages) + 1}. {title[0]}"
        stages.append((st["id"], title))
    stage_ids = [s[0] for s in stages]

    stations = []
    last_col = {}
    for i, raw in enumerate(spec.get("stations", []), 1):
        branches = [(text_lines(b.get("tool")), text_lines(b.get("desc"))) for b in raw.get("branches", [])]
        names = text_lines(raw.get("tool"))
        label = " / ".join(" ".join(n) for n in [names] + [b[0] for b in branches] if n) or "no tool"
        where = f"station {i} ({label})"
        if bool(names) == bool(branches):
            errors.append(f"{where}: needs 'tool' or 'branches', not both")
        if branches and (len(branches) < 2 or not all(b[0] for b in branches)):
            errors.append(f"{where}: 'branches' needs two or more tables, each with 'tool'")
        lane, stage = raw.get("lane"), raw.get("stage")
        if lane not in lanes:
            errors.append(f"{where}: unknown lane '{lane}'")
            continue
        if stage not in stage_ids:
            errors.append(f"{where}: unknown stage '{stage}'")
            continue
        # Default column: the next column of this lane in this stage.
        col = raw.get("col", last_col.get((stage, lane), 0) + 1)
        if not isinstance(col, int) or col < 1:
            errors.append(f"{where}: 'col' must be a positive integer")
            continue
        last_col[(stage, lane)] = col
        st_lines = raw.get("lines", lane_lines(lanes[lane]))
        outside = [x for x in st_lines if x not in lane_lines(lanes[lane])]
        if outside or not st_lines:
            errors.append(f"{where}: lines {outside or st_lines} are not in lane '{lane}'")
            continue
        icon = raw.get("input")
        if icon is not None and not icon.get("label"):
            errors.append(f"{where}: 'input' needs 'label'")
        icon = icon and (icon.get("label", ""), text_lines(icon.get("caption")), text_lines(icon.get("note")))
        stations.append(Station(names=names, branches=branches, desc=text_lines(raw.get("desc")),
                                note=text_lines(raw.get("note")), lane=lane, stage=stage, col=col,
                                lines=list(st_lines), optional=bool(raw.get("optional")), input=icon))

    seen = {}
    for s in stations:
        key = (s.stage, s.col, s.lane)
        if key in seen:
            errors.append(f"two stations share lane '{s.lane}', stage '{s.stage}' and col {s.col}")
        seen[key] = True
    for lid in lines:
        cols = [(stage_ids.index(s.stage), s.col) for s in stations if lid in s.lines]
        if not cols:
            errors.append(f"line '{lid}' has no stations")
        elif len(cols) != len(set(cols)):
            errors.append(f"line '{lid}' has two stations in the same column")
    for sid, _ in stages:
        if not any(s.stage == sid for s in stations):
            errors.append(f"stage '{sid}' has no stations")

    def files(key):
        out = []
        for f in spec.get(key, []):
            if not f.get("label"):
                errors.append(f"{key} {f}: needs 'label'")
            if f.get("lane") is not None and f["lane"] not in lanes:
                errors.append(f"{key} {f}: unknown lane '{f['lane']}'")
            out.append((f.get("label", ""), text_lines(f.get("caption")), f.get("lane")))
        return out

    inputs, outputs = files("inputs"), files("outputs")
    if errors:
        raise SpecError("\n".join(f"- {e}" for e in errors))
    legend = {"mandatory": "Mandatory", "optional": "Optional", "input": "Input", "output": "Output", "note": ""}
    legend.update(spec.get("legend", {}))
    return Model(title=spec.get("title", {}), style=style, legend=legend, lines=lines, lanes=lanes,
                 stages=stages, stations=stations, inputs=inputs, outputs=outputs)


# ---------------------------------------------------------------------------
# Geometry relative to a lane (y = 0 is the first track of the lane)
# ---------------------------------------------------------------------------


def track(m, lane, line):
    return next(i for i, t in enumerate(m.lanes[lane]) if line in t) * GAP


def station_span(m, s):
    """y of the first and the last track of the station lines."""
    ys = [track(m, s.lane, ln) for ln in s.lines]
    return min(ys), max(ys)


def block_height(items):
    """Height from the first to the last baseline of a block of (kind, text) items."""
    return sum(LH[k] for k, _ in items[1:])


def fork_offsets(m, s):
    """y offset of each branch of a station. A normal station has one branch."""
    n = len(s.branches)
    if n < 2:
        return [0.0]
    y0, y1 = station_span(m, s)
    if n == 2:  # the station description and note go between the two branches
        middle = [("desc", d) for d in s.desc] + [("note", t) for t in s.note]
        need = block_height(middle) + (LH[middle[0][0]] if middle else 0) + 2 * R + 24
        pitch = max(2 * FORK_DY, need)
    else:  # the name and description of each branch go above it
        above = max(len(nm) * LH["name"] + len(d) * LH["desc"] for nm, d in s.branches[1:])
        pitch = max(FORK_PITCH, above + 2 * R + 16)
    pitch += y1 - y0
    return [(i - (n - 1) / 2) * pitch for i in range(n)]


def stack(x, base, items, up):
    """Baselines for a block of (kind, text) items.

    :param up: True: the last item is on ``base`` and the block goes up. False: the first
        item is on ``base`` and the block goes down.
    """
    out, y = [], base
    seq = reversed(items) if up else items
    for i, (kind, t) in enumerate(seq):
        if i:
            y += -LH[kind] if up else LH[kind]
        out.append((x, y, t, kind))
    return out


def station_labels(m, s):
    """Labels of a station: (x, baseline y, text, kind). kind is name, desc or note."""
    below = [("desc", d) for d in s.desc] + [("note", t) for t in s.note]
    if not s.fork:
        return (stack(s.x, s.top - NAME_DY, [("name", n) for n in s.names], up=True)
                + stack(s.x, s.bottom + DESC_DY, below, up=False))
    y0, y1 = station_span(m, s)
    offs = fork_offsets(m, s)

    def block(b):
        names, desc = s.branches[b]
        return [("name", n) for n in names] + [("desc", d) for d in desc]

    if len(s.branches) == 2:  # labels outside the branches, station description between them
        out = stack(s.x, offs[0] + y0 - NAME_DY, block(0), up=True)
        out += stack(s.x, offs[1] + y1 + DESC_DY + 2, block(1), up=False)
        if below:
            first = (y0 + y1) / 2 - block_height(below) / 2 + 0.3 * SIZE[below[0][0]]
            out += stack(s.x, first, below, up=False)
        return out
    out = []
    for b, off in enumerate(offs):  # labels above each branch, station description below the last one
        out += stack(s.x, off + y0 - NAME_DY, block(b), up=True)
    return out + stack(s.x, offs[-1] + y1 + DESC_DY, below, up=False)


def label_rect(m, x, y, t, kind):
    w = text_width(t, SIZE[kind], m.weight if kind == "name" else 400)
    return (x - w / 2, x + w / 2, y - 0.8 * SIZE[kind], y + 0.25 * SIZE[kind])


def station_width(m, s):
    return max([2 * R] + [r[1] - r[0] for r in (label_rect(m, *lb) for lb in station_labels(m, s))])


def icon_block(m, s):
    """x range and height of the input icon of a station, with its caption on the right."""
    w, h = SMALL_ICON
    texts = [text_width(c, SIZE["desc"]) for c in s.input[1]] + [text_width(c, SIZE["note"]) for c in s.input[2]]
    return s.x - w / 2 - 8, s.x + w / 2 + 18 + max(texts, default=0), h


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


class Layout:
    """All the positions of a map."""

    def __init__(self, m):
        self.m = m
        # Pass 1 finds the gaps between columns that have curves. Pass 2 makes these gaps wide
        # enough for a curve between the labels. The order of the columns does not change.
        self.place_columns({})
        self.make_routes()
        need = {i: (self.col_w[i] + self.col_w[i + 1]) / 2 + CURVE + 2 * CLEAR for i in self.curve_gaps}
        self.place_columns(need)
        self.make_routes()
        self.find_tracks_in_use()
        self.place_station_icons()
        self.place_lanes()

    # -- x ------------------------------------------------------------------
    def place_columns(self, need):
        """Set the x of each column and stage box.

        :param need: Minimum distance between column i and column i + 1, for the gaps with curves.
        """
        m = self.m
        widths = {}
        for s in m.stations:
            key = (s.stage, s.col)
            widths[key] = max(widths.get(key, 0), station_width(m, s))

        self.in_w = self.file_column_width(m.inputs)
        self.x_start = MARGIN + (self.in_w + TERMINUS_GAP if m.inputs else 0)
        cursor = self.x_start + TERMINUS_GAP
        col_x, boxes = {}, []
        prev_x = prev_w = None
        gap = 1  # gap 0 is between the input terminus and the first column
        for sid, title in m.stages:
            x0 = cursor
            for i, c in enumerate(sorted(c for (st, c) in widths if st == sid)):
                w = widths[(sid, c)]
                if prev_x is None:
                    x = x0 + STAGE_PAD + w / 2
                elif i == 0:
                    x = max(x0 + STAGE_PAD + w / 2, prev_x + MIN_PITCH, prev_x + need.get(gap, 0))
                    gap += 1
                else:
                    x = prev_x + max((prev_w + w) / 2 + COL_GAP, MIN_PITCH, need.get(gap, 0))
                    gap += 1
                col_x[(sid, c)] = x
                prev_x, prev_w = x, w
            title_w = max(text_width(t, SIZE["name"], 700) for t in title) + 48
            x1 = max(prev_x + prev_w / 2 + STAGE_PAD, x0 + title_w)
            boxes.append([x0, x1, title])
            cursor = x1 + STAGE_GAP
        for s in m.stations:
            s.x = col_x[(s.stage, s.col)]
        self.boxes = boxes
        self.x_end = boxes[-1][1] + 12
        keys = sorted(col_x, key=lambda k: col_x[k])
        self.columns = [self.x_start] + [col_x[k] for k in keys] + [self.x_end]
        self.col_w = [0] + [widths[k] for k in keys] + [0]
        self.out_w = self.file_column_width(m.outputs)
        self.width = self.x_end + 11 + (TERMINUS_GAP + self.out_w if m.outputs else 0) + MARGIN

    @staticmethod
    def file_column_width(files):
        if not files:
            return 0
        return max([ICON[0] + 16] + [text_width(c, SIZE["desc"]) for _, caps, _ in files for c in caps])

    def curve(self, gap):
        """x range of a curve in the gap between column gap and column gap + 1.

        The curve is in the middle of the free space between the labels of the two columns.
        """
        self.curve_gaps.add(gap)
        a = self.columns[gap] + self.col_w[gap] / 2
        b = self.columns[gap + 1] - self.col_w[gap + 1] / 2
        mid, w = (a + b) / 2, min(CURVE, self.columns[gap + 1] - self.columns[gap] - 2 * R - 20)
        return mid - w / 2, mid + w / 2

    # -- routes -------------------------------------------------------------
    def make_routes(self):
        """For each line, one route per fork branch: a list of (x, lane, y in lane)."""
        m = self.m
        self.routes = {}
        self.curve_gaps = set()
        for ln, info in m.lines.items():
            stations = sorted((s for s in m.stations if ln in s.lines), key=lambda s: s.x)
            n_branches = max(len(fork_offsets(m, s)) for s in stations)
            routes = []
            for b in range(n_branches):
                pts = []
                prev = None
                for s in stations:
                    offs = fork_offsets(m, s)
                    y = track(m, s.lane, ln) + offs[min(b, len(offs) - 1)]
                    if prev is None:
                        if info["start"] == "input":
                            pts.append((self.x_start, s.lane, y))
                    elif (prev[1], prev[2]) != (s.lane, y):
                        # After a fork, the curve is just before the next station; else just after the last one.
                        gap = self.columns.index(s.x) - 1 if prev[3].fork else self.columns.index(prev[0])
                        t0, t1 = self.curve(gap)
                        pts += [(t0, prev[1], prev[2]), (t1, s.lane, y)]
                    pts.append((s.x, s.lane, y))
                    prev = (s.x, s.lane, y, s)
                if info["end"] == "output":
                    pts.append((self.x_end, prev[1], prev[2]))
                routes.append(pts)
            self.routes[ln] = routes

    def segments(self, lane):
        """Straight parts of the routes in a lane: (x0, x1, y)."""
        for routes in self.routes.values():
            for pts in routes:
                for (xa, la, ya), (xb, lb, yb) in pairwise(pts):
                    if la == lb == lane and ya == yb:
                        yield xa, xb, ya

    def find_tracks_in_use(self):
        """Set the highest and lowest track in use at each station: the labels go outside them."""
        m = self.m
        for s in self.m.stations:
            tracks = {i * GAP for i in range(len(m.lanes[s.lane]))}
            used = {y for x0, x1, y in self.segments(s.lane) if x0 <= s.x <= x1 and y in tracks}
            used.update(station_span(m, s))
            s.top, s.bottom = min(used), max(used)

    def place_station_icons(self):
        """Put each station input icon above the labels of the stations under it."""
        m = self.m
        for s in m.stations:
            if not s.input:
                continue
            x0, x1, h = icon_block(m, s)
            tops = [r[2] for o in m.stations if o.lane == s.lane
                    for r in (label_rect(m, *lb) for lb in station_labels(m, o)) if r[0] < x1 and r[1] > x0]
            tops.append(s.top - R)
            s.icon_y = min(tops) - 14 - h

    # -- y ------------------------------------------------------------------
    def lane_rects(self, lane):
        """Rectangles (x0, x1, y0, y1) of labels, stations, icons and tracks, relative to the lane."""
        m = self.m
        rects = []
        for s in m.stations:
            if s.lane != lane:
                continue
            rects += [label_rect(m, *lb) for lb in station_labels(m, s)]
            y0, y1 = station_span(m, s)
            for off in fork_offsets(m, s):
                rects.append((s.x - R, s.x + R, y0 + off - R, y1 + off + R))
            if s.input:
                x0, x1, h = icon_block(m, s)
                rects.append((x0, x1, s.icon_y, s.icon_y + h))
        rects += [(x0, x1, y - LW / 2, y + LW / 2) for x0, x1, y in self.segments(lane)]
        return rects

    def place_lanes(self):
        m = self.m
        rects = {lane: self.lane_rects(lane) for lane in m.lanes}
        self.lane_y = {}
        prev = None
        for lane in m.lanes:
            y = 0.0 if prev is None else self.lane_y[prev] + (len(m.lanes[prev]) - 1) * GAP + LANE_GAP
            for done, yd in self.lane_y.items():
                for a in rects[done]:
                    for b in rects[lane]:
                        if a[0] < b[1] + CLEAR and b[0] < a[1] + CLEAR:
                            y = max(y, yd + a[3] - b[2] + CLEAR)
            self.lane_y[lane] = y
            prev = lane

        # Stage boxes: the same top and bottom for all; room for the titles at the top.
        absolute = [(r[0], r[1], r[2] + self.lane_y[ln], r[3] + self.lane_y[ln]) for ln in m.lanes for r in rects[ln]]
        top = math.inf
        for x0, x1, title in self.boxes:
            inside = [r[2] for r in absolute if r[0] < x1 and r[1] > x0]
            top = min(top, min(inside) - (46 + (len(title) - 1) * TITLE_LH + 20))
        bottom = max(r[3] for r in absolute) + 30
        shift = 20 - top
        for ln in self.lane_y:
            self.lane_y[ln] += shift
        self.box_top, self.box_bottom = 20, bottom + shift

    def y(self, lane, rel):
        return self.lane_y[lane] + rel


# ---------------------------------------------------------------------------
# SVG
# ---------------------------------------------------------------------------


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def fmt(v):
    return f"{v:.1f}".rstrip("0").rstrip(".")


def text(x, y, s, size=SIZE["desc"], weight=400, anchor="middle", fill="#000"):
    return (f'<text x="{fmt(x)}" y="{fmt(y)}" font-family="{FONT}" font-size="{fmt(size)}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}">{esc(s)}</text>')


def path_through(points):
    """Polyline through points; a change of y is a smooth S-curve."""
    d = f"M{fmt(points[0][0])},{fmt(points[0][1])}"
    for (x0, y0), (x1, y1) in pairwise(points):
        if y0 == y1:
            d += f" L{fmt(x1)},{fmt(y1)}"
        else:
            xm = (x0 + x1) / 2
            d += f" C{fmt(xm)},{fmt(y0)} {fmt(xm)},{fmt(y1)} {fmt(x1)},{fmt(y1)}"
    return d


def pill(x, y0, y1, r, **attrs):
    extra = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in attrs.items())
    return (f'<rect x="{fmt(x - r)}" y="{fmt(y0 - r)}" width="{fmt(2 * r)}" height="{fmt(y1 - y0 + 2 * r)}" '
            f'rx="{fmt(r)}" {extra}/>')


def draw_station(lay, s):
    m = lay.m
    y0, y1 = station_span(m, s)
    out = []
    for off in fork_offsets(m, s):
        a, b = lay.y(s.lane, y0 + off), lay.y(s.lane, y1 + off)
        out.append(pill(s.x, a, b, R, fill="#fff", stroke="#000", stroke_width=4))
        if not s.optional:
            out.append(pill(s.x, a, b, DOT, fill="#000"))
    for x, y, t, kind in station_labels(m, s):
        out.append(text(x, lay.y(s.lane, y), t, SIZE[kind], m.weight if kind == "name" else 400, fill=FILL[kind]))
    if s.input:
        label, caption, note = s.input
        out += file_icon(s.x - SMALL_ICON[0] / 2, lay.y(s.lane, s.icon_y), label, caption, note, SMALL_ICON, True)
    return out


def file_icon(x, y, label, caption, note=(), size=ICON, caption_right=False):
    """nf-core style document icon; (x, y) is the top-left corner.

    The caption goes below the icon, or on its right (for a station input). The note
    lines (grey) follow the caption.
    """
    w, h = size
    sc = w / 62
    fold = 18 * sc
    page = f"M{fmt(x)},{fmt(y)} H{fmt(x + w - fold)} L{fmt(x + w)},{fmt(y + fold)} V{fmt(y + h)} H{fmt(x)} Z"
    corner = f"M{fmt(x + w - fold)},{fmt(y)} V{fmt(y + fold)} H{fmt(x + w)}"
    out = [f'<path d="{page}" fill="#fff" stroke="#000" stroke-width="5" stroke-linejoin="round"/>',
           f'<path d="{corner}" fill="none" stroke="#000" stroke-width="4" stroke-linejoin="round"/>']
    bw = max(w + 16, text_width(label, 18 * sc, 700) + 16)
    out.append(f'<rect x="{fmt(x + w / 2 - bw / 2)}" y="{fmt(y + h * 0.42)}" width="{fmt(bw)}" '
               f'height="{fmt(26 * sc)}" rx="4" fill="#000"/>')
    out.append(text(x + w / 2, y + h * 0.42 + 19 * sc, label, 18 * sc, 700, fill="#fff"))
    items = [("desc", c) for c in caption] + [("note", n) for n in note]
    if caption_right:
        first = y + h / 2 - block_height(items) / 2 + 8
        for _, yy, t, kind in stack(0, first, items, up=False):
            out.append(text(x + w + 18, yy, t, SIZE[kind], anchor="start", fill=FILL[kind]))
    else:
        for _, yy, t, kind in stack(0, y + h + 32, items, up=False):
            out.append(text(x + w / 2, yy, t, SIZE[kind], fill=FILL[kind]))
    return out


def file_column(lay, files, centre_x):
    """Icons in a column. An icon with a lane is level with the lane; the others are stacked.

    :return: A tuple ``(svg parts, bottom y)``.
    """
    gap = 40
    heights = [ICON[1] + 40 + LH["desc"] * len(c) for _, c, _ in files]
    y = lay.box_top
    if not any(lane for _, _, lane in files):  # centre the stack on the boxes if it fits
        y = max(y, (lay.box_top + lay.box_bottom - sum(heights) - gap * (len(files) - 1)) / 2)
    out = []
    for (label, caption, lane), h in zip(files, heights):
        if lane:
            centre = lay.y(lane, (len(lay.m.lanes[lane]) - 1) * GAP / 2)
            y = max(y, centre - h / 2 + 20)
        out += file_icon(centre_x - ICON[0] / 2, y, label, caption)
        y += h + gap
    return out, y - gap


def termini(lay, end):
    """Input triangles or output bars: one for each lane where lines start or end."""
    groups = {}
    for ln, routes in lay.routes.items():
        if lay.m.lines[ln]["end" if end else "start"] == "station":
            continue
        _, lane, y = routes[0][-1 if end else 0]
        groups.setdefault(lane, []).append(lay.y(lane, y))
    out = []
    for ys in groups.values():
        y0, y1 = min(ys) - 12, max(ys) + 12
        if end:
            out.append(f'<rect x="{fmt(lay.x_end)}" y="{fmt(y0)}" width="11" height="{fmt(y1 - y0)}" fill="#000"/>')
        else:
            xs = lay.x_start
            out.append(f'<path d="M{fmt(xs)},{fmt(y0)} L{fmt(xs + 26)},{fmt((y0 + y1) / 2)} '
                       f'L{fmt(xs)},{fmt(y1)} Z" fill="#000"/>')
    return out


def legend(lay, top):
    """Title (bottom left) and legend box (bottom right). Returns (svg parts, bottom y)."""
    m = lay.m
    out = []
    left = lay.boxes[0][0]
    t = m.title
    title_w = max(text_width(t.get("org", ""), 46, 700), text_width(t.get("name", ""), 86, 800),
                  text_width(t.get("subtitle", ""), 28))
    if t.get("org"):
        out.append(text(left, top + 70, t["org"], 46, 700, "start", "#555"))
    if t.get("name"):
        out.append(text(left, top + 140, t["name"], 86, 800, "start"))
    if t.get("subtitle"):
        out.append(text(left, top + 190, t["subtitle"], 28, 400, "start", "#444"))

    symbols = [("mandatory", m.legend["mandatory"])]
    if any(s.optional for s in m.stations):
        symbols.append(("optional", m.legend["optional"]))
    symbols += [("in", m.legend["input"]), ("out", m.legend["output"])]
    lines_w = 120 + max(text_width(v["label"], 28) for v in m.lines.values())
    sym_w = 52 + max(text_width(lab, 28) for _, lab in symbols)
    note = m.legend["note"]
    bw = max(30 + lines_w + 40 + sym_w + 30, text_width(note, 24) + 60)
    rows = max(len(m.lines), len(symbols))
    bh = 100 + (rows - 1) * 42 + 40 + (40 if note else 0)
    bx, by = lay.boxes[-1][1] - bw, top
    if bx < left + title_w + 60:  # no room on the right of the title: put the legend below it
        bx, by = left, top + 230
    out.append(f'<rect x="{fmt(bx)}" y="{fmt(by)}" width="{fmt(bw)}" height="{fmt(bh)}" rx="24" fill="#E9E9E9"/>')
    out.append(text(bx + 30, by + 50, "Legend", 34, 700, "start"))
    for i, v in enumerate(m.lines.values()):
        y = by + 100 + i * 42
        out.append(f'<line x1="{fmt(bx + 30)}" y1="{fmt(y)}" x2="{fmt(bx + 100)}" y2="{fmt(y)}" '
                   f'stroke="{v["color"]}" stroke-width="{LW + 3}"/>')
        out.append(text(bx + 120, y + 8, v["label"], 28, 400, "start"))
    sx = bx + 30 + lines_w + 40
    for i, (kind, lab) in enumerate(symbols):
        y, cx = by + 100 + i * 42, sx + 20
        if kind in ("mandatory", "optional"):
            out.append(f'<circle cx="{fmt(cx)}" cy="{fmt(y)}" r="{R}" fill="#fff" stroke="#000" stroke-width="4"/>')
            if kind == "mandatory":
                out.append(f'<circle cx="{fmt(cx)}" cy="{fmt(y)}" r="{DOT}" fill="#000"/>')
        elif kind == "in":
            out.append(f'<path d="M{fmt(cx - 10)},{fmt(y - 14)} L{fmt(cx + 12)},{fmt(y)} '
                       f'L{fmt(cx - 10)},{fmt(y + 14)} Z" fill="#000"/>')
        else:
            out.append(f'<rect x="{fmt(cx - 5)}" y="{fmt(y - 15)}" width="10" height="30" fill="#000"/>')
        out.append(text(sx + 52, y + 9, lab, 28, 400, "start"))
    if note:
        out.append(text(bx + 30, by + bh - 25, note, 24, 400, "start", "#555"))
    return out, max(by + bh, top + 200)


def render(m, font_b64=None):
    """Make the SVG.

    :param m: The :class:`Model`.
    :param font_b64: A TTF font as base64, to embed in the SVG, or ``None``.
    :return: A tuple ``(svg text, width, height)``.
    """
    lay = Layout(m)
    body = []
    for x0, x1, title in lay.boxes:
        body.append(f'<rect x="{fmt(x0)}" y="{fmt(lay.box_top)}" width="{fmt(x1 - x0)}" '
                    f'height="{fmt(lay.box_bottom - lay.box_top)}" rx="24" fill="#E9E9E9"/>')
        body += [text(x0 + 24, lay.box_top + 46 + i * TITLE_LH, t, SIZE["name"], 700, "start", "#8A8A8A")
                 for i, t in enumerate(title)]
    for ln, routes in lay.routes.items():
        for pts in routes:
            d = path_through([(x, lay.y(lane, y)) for x, lane, y in pts])
            body.append(f'<path d="{d}" fill="none" stroke="{m.lines[ln]["color"]}" stroke-width="{LW}" '
                        f'stroke-linejoin="round"/>')
    body += termini(lay, end=False) + termini(lay, end=True)
    for s in m.stations:
        body += draw_station(lay, s)
    # The legend is below the stage boxes. The file icons are on the left and right of the boxes,
    # so they can go lower than the boxes.
    parts, bottom = legend(lay, lay.box_bottom + 40)
    body += parts
    for files, x in [(m.inputs, MARGIN + lay.in_w / 2), (m.outputs, lay.x_end + 11 + TERMINUS_GAP + lay.out_w / 2)]:
        if files:
            parts, icons_bottom = file_column(lay, files, x)
            body += parts
            bottom = max(bottom, icons_bottom)

    w, h = math.ceil(lay.width), math.ceil(bottom + 40)
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">']
    if font_b64:
        head.append("<defs><style>@font-face{font-family:'Maven Pro';font-weight:400 900;"
                    f"src:url(data:font/ttf;base64,{font_b64}) format('truetype');}}</style></defs>")
    head.append(f'<rect width="{w}" height="{h}" fill="#fff"/>')
    return "\n".join(head + body + ["</svg>"]), w, h


CHROME = ["/Applications/Google Chrome.app/Contents/MacOS/Google Chrome", "google-chrome", "chromium",
          "chromium-browser", "chrome"]


def render_png(svg, png, w, h):
    """Take a screenshot of the SVG with headless Chrome (it shows the embedded font)."""
    chrome = next((c if Path(c).exists() else shutil.which(c) for c in CHROME
                   if Path(c).exists() or shutil.which(c)), None)
    if not chrome:
        sys.exit("no Chrome or Chromium found: open the SVG in a browser and export the PNG by hand")
    subprocess.run([chrome, "--headless", f"--screenshot={png.resolve()}", f"--window-size={w},{h}",
                    "--hide-scrollbars", svg.resolve().as_uri()], check=True, capture_output=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("spec", type=Path, help="TOML spec of the map")
    ap.add_argument("--out", type=Path, help="SVG file (default: the spec path with .svg)")
    ap.add_argument("--font", type=Path, help="TTF to embed in the SVG (Maven Pro)")
    ap.add_argument("--png", type=Path, help="also write a PNG (needs Chrome or Chromium)")
    a = ap.parse_args(argv)
    try:
        with open(a.spec, "rb") as f:
            model = parse(tomllib.load(f))
    except (SpecError, tomllib.TOMLDecodeError) as e:
        sys.exit(f"{a.spec}: errors in the spec:\n{e}")
    font = base64.b64encode(a.font.read_bytes()).decode() if a.font else None
    svg, w, h = render(model, font)
    out = a.out or a.spec.with_suffix(".svg")
    out.write_text(svg)
    print(f"wrote {out} ({w} x {h})")
    if a.png:
        render_png(out, a.png, w, h)
        print(f"wrote {a.png}")


if __name__ == "__main__":
    main()
