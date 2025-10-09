#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ---------------- Helpers ----------------

def detect_gff_gtf(file_path):
    with open(file_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                print("Detected GFF3 file format.")
                return "GFF3"
            elif '"' in line:
                print("Detected GTF file format.")
                return "GTF"
    print(f"Cannot detect GFF3 or GTF format in {file_path}. Exiting.", file=sys.stderr)
    sys.exit(1)

def parse_gff_line(line):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 9:
        return None
    return {
        "seqid": parts[0],
        "source": parts[1],
        "type": parts[2],
        "start": int(parts[3]),
        "end": int(parts[4]),
        "score": parts[5],
        "strand": parts[6],
        "phase": parts[7],
        "attributes": parts[8],
    }

def parse_attributes(attr_str, file_format="GFF3"):
    attrs = {}
    parts = [p.strip() for p in attr_str.split(";") if p.strip()]
    if file_format == "GFF3":
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                attrs[k.strip()] = v.strip().strip('"')
    elif file_format == "GTF":
        for part in parts:
            if " " in part:
                k, v = part.split(" ", 1)
                attrs[k.strip()] = v.strip().strip('"')
    return attrs

def pack_intervals(intervals):
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    track_ends = []
    packed_intervals = []
    for start, end, payload in sorted_intervals:
        placed = False
        for idx, last_end in enumerate(track_ends):
            if start > last_end:
                track_ends[idx] = end
                packed_intervals.append((start, end, payload, idx))
                placed = True
                break
        if not placed:
            track_ends.append(end)
            packed_intervals.append((start, end, payload, len(track_ends)-1))
    return packed_intervals, len(track_ends)

def attr_hover_text(attr_dict, file_format="GFF3"):
    hover_lines = []
    feature_id = attr_dict.get("locus_tag") or attr_dict.get("ID") or attr_dict.get("id")
    if feature_id:
        hover_lines.append(f"Feature: {feature_id}")
    gene_name = attr_dict.get("gene") or attr_dict.get("Name") or attr_dict.get("Gene") or attr_dict.get("name")
    if gene_name:
        hover_lines.append(f"Gene: {gene_name}")
    product = attr_dict.get("product") or attr_dict.get("description") or attr_dict.get("note")
    if product:
        hover_lines.append(f"Description: {product}")
    if not hover_lines:
        if file_format.upper() == "GFF3":
            kv_pairs = [f"{k}={v}" for k,v in attr_dict.items()]
            if kv_pairs: hover_lines.append("; ".join(kv_pairs))
        elif file_format.upper() == "GTF":
            kv_pairs = [f'{k}="{v}"' for k,v in attr_dict.items()]
            if kv_pairs: hover_lines.append("; ".join(kv_pairs))
    return "<br>".join(hover_lines)

def strand_color(s):
    if s=="+": return "blue"
    elif s=="-": return "red"
    else: return "grey"

def find_overlaps(intervals):
    overlaps = 0
    intervals_sorted = sorted(intervals, key=lambda x: (x[0], x[1]))
    for i in range(len(intervals_sorted)):
        s1, e1, p1 = intervals_sorted[i]
        chrom1, start1, end1, strand1, attrs1 = p1
        for j in range(i+1, len(intervals_sorted)):
            s2, e2, p2 = intervals_sorted[j]
            chrom2, start2, end2, strand2, attrs2 = p2
            if s2 > e1: break
            if chrom1 == chrom2 and s2 <= e1:
                id1 = attrs1.get("locus_tag") or attrs1.get("ID") or f"{chrom1}:{start1}-{end1}"
                id2 = attrs2.get("locus_tag") or attrs2.get("ID") or f"{chrom2}:{start2}-{end2}"
                print(f"Overlap on {chrom1}: {id1} ({start1}-{end1}) overlaps {id2} ({start2}-{end2})")
                overlaps += 1
    return overlaps

# ---------------- Main ----------------
def main():
    p = argparse.ArgumentParser(description="Interactive genome plot: value panel + feature panel")
    p.add_argument("-v","--values_tsv", required=True)
    p.add_argument("-f","--features_gff", required=True)
    p.add_argument("-t","--feature_type", required=True)
    p.add_argument("--output", default="genome_plot.html")
    p.add_argument("--feature-height", type=float, default=0.2)
    p.add_argument("--max-feature-tracks", type=int, default=5)
    p.add_argument("--title", default="Genome Browser")
    p.add_argument("--value-ytitle", default="Value")
    p.add_argument("--feature-ytitle", default="Feature")
    args = p.parse_args()

    # Detect file format
    file_format = detect_gff_gtf(args.features_gff)
    if file_format == "GTF":
        print(
            "Warning: Detected GTF. Consider adding '##sequence-region' lines to determine chromosome lengths.",
            file=sys.stderr
        )

    # Parse ##sequence-region
    chrom_lengths = {}
    chrom_order = []
    with open(args.features_gff) as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    chrom = parts[1]
                    end = int(parts[3])
                    chrom_lengths[chrom] = end
                    if chrom not in chrom_order:
                        chrom_order.append(chrom)
    if not chrom_order:
        print("No chromosomes found in ##sequence-region. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Read features
    chrom_feature_map = defaultdict(list)
    with open(args.features_gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parsed = parse_gff_line(line)
            if parsed and parsed["type"]==args.feature_type and parsed["seqid"] in chrom_order:
                chrom_feature_map[parsed["seqid"]].append(parsed)

    # Load values TSV
    df_vals = pd.read_csv(args.values_tsv, sep="\t", dtype={"chr": str})
    mandatory_cols = {"chr","start","end"}
    missing = mandatory_cols - set(df_vals.columns)
    if missing:
        print(f"Missing mandatory column(s): {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)
    df_vals = df_vals[df_vals["chr"].isin(chrom_order)].copy()
    if df_vals.empty:
        print("No values for chromosomes in ##sequence-region. Exiting.", file=sys.stderr)
        sys.exit(1)
    df_vals["mid"] = (df_vals["start"] + df_vals["end"])/2.0

    # Chromosome offsets
    offsets = {}
    cum = 0
    for c in chrom_order:
        offsets[c] = cum
        cum += chrom_lengths[c]
    total_length = cum
    df_vals["x"] = df_vals["mid"] + df_vals["chr"].map(offsets)

    # Sample columns
    sample_cols = [c for c in df_vals.columns if c not in ["chr","start","end","mid","x"]]
    for col in sample_cols:
        df_vals[col] = pd.to_numeric(df_vals[col], errors="coerce")
        non_numeric = df_vals[col][df_vals[col].apply(lambda x: not pd.api.types.is_number(x))]
        if not non_numeric.empty:
            print(f"Warning: Column {col} contains non-numeric values. They will be skipped.", file=sys.stderr)

    # Prepare features
    all_feature_intervals = []
    for chrom in chrom_order:
        for f in chrom_feature_map.get(chrom, []):
            attrs = parse_attributes(f["attributes"], file_format)
            x0 = f["start"] + offsets[chrom]
            x1 = f["end"] + offsets[chrom]
            all_feature_intervals.append((x0, x1, (chrom, f["start"], f["end"], f["strand"], attrs)))

    packed, n_tracks = pack_intervals(all_feature_intervals)
    n_tracks = min(n_tracks, args.max_feature_tracks)
    track_y = lambda track: (-0.95 + track*(1.9/n_tracks), -0.95 + (track+1)*(1.9/n_tracks))

    # Line segments
    lines_by_color = defaultdict(lambda: {"x":[],"y":[],"hover":[]})
    for s,e,payload,track in packed:
        if track >= args.max_feature_tracks:
            track = args.max_feature_tracks-1
        chrom,start_,end_,strand,attrs = payload
        y0,y1 = track_y(track)
        y_mid = (y0+y1)/2
        lines_by_color[strand_color(strand)]["x"] += [s,e,np.nan]
        lines_by_color[strand_color(strand)]["y"] += [y_mid,y_mid,np.nan]
        hover_text = "<br>".join([f"<b>{chrom}:{start_}-{end_}</b>", attr_hover_text(attrs, file_format)])
        lines_by_color[strand_color(strand)]["hover"] += [hover_text]*3

    print("Checking overlaps...")
    nb_overlaps = find_overlaps(all_feature_intervals)
    print("Total overlaps:", nb_overlaps)
    print("Tracks to draw:", n_tracks)

    # ---------------- Build figure ----------------
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True,
        row_heights=[1-args.feature_height, args.feature_height],
        vertical_spacing=0.02,
        specs=[[{"type":"xy"}],[{"type":"xy"}]]
    )

    # Value panel
    for col in sample_cols:
        valid = df_vals[~df_vals[col].isna()]
        fig.add_trace(go.Scatter(
            x=valid["x"], y=valid[col],
            mode="lines+markers",
            name=col,
            hovertemplate="<b>%{text}</b><br>x: %{x}<br>y: %{y}<extra></extra>",
            text=valid.apply(lambda r: f"{r['chr']}:{int(r['start'])}-{int(r['end'])}", axis=1)
        ), row=1, col=1)

    # Feature panel
    for color,d in lines_by_color.items():
        fig.add_trace(go.Scatter(
            x=d["x"], y=d["y"],
            mode="lines",
            line=dict(color=color, width=6),
            hoverinfo="text",
            text=d["hover"],
            showlegend=False
        ), row=2, col=1)

    # Chromosome separators & labels
    chrom_midpoints = []
    chrom_ranges = {}
    for i,c in enumerate(chrom_order):
        if i>0:
            fig.add_vline(x=offsets[c], line=dict(dash="dash", width=1), row=1,col=1)
            fig.add_vline(x=offsets[c], line=dict(dash="dash", width=1), row=2,col=1)
        chrom_midpoints.append(offsets[c]+chrom_lengths[c]/2.0)
        chrom_ranges[c] = (offsets[c], offsets[c]+chrom_lengths[c])

    fig.update_yaxes(title_text=args.value_ytitle, row=1,col=1)
    fig.update_yaxes(title_text=args.feature_ytitle, row=2,col=1, showticklabels=False, range=[-1,1])
    fig.update_xaxes(tickmode="array", tickvals=chrom_midpoints, ticktext=chrom_order, row=2,col=1)
    fig.update_xaxes(title_text="Genomic position (concatenated chromosomes)", row=2,col=1)
    fig.update_xaxes(rangeslider=dict(visible=True), row=2,col=1)

    # ---------------- Dropdown menu ----------------
    buttons = [dict(label="All", method="relayout",
                    args=[{"xaxis.range":[0,total_length], "xaxis2.range":[0,total_length]}])]
    for chrom, (x0,x1) in chrom_ranges.items():
        buf = 0.01*(x1-x0)
        buttons.append(dict(label=chrom, method="relayout",
                            args=[{"xaxis.range":[x0-buf,x1+buf], "xaxis2.range":[x0-buf,x1+buf]}]))
    fig.update_layout(
        updatemenus=[dict(buttons=buttons, direction="down", showactive=True, x=1, y=1.1,
                          xanchor="left", yanchor="top")]
    )

    # ---------------- Layout ----------------
    fig.update_layout(title=args.title,
                      height=800,
                      hovermode="closest",
                      legend=dict(orientation="v", x=1, y=1.02, xanchor="left", yanchor="top"))

    fig.write_html(args.output, include_plotlyjs=True)
    print(f"Plot saved to {args.output}")

if __name__=="__main__":
    main()
