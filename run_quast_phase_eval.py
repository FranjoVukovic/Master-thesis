#!/usr/bin/env python3
"""Run QUAST on all phase FASTAs vs paternal/maternal references, then summarize."""

import subprocess
import os
import sys
from pathlib import Path

BASE = Path("/mnt/share1_Jabba/fvukovic/chromosome_18")
QUAST = "/home/fvukovic/miniconda3/envs/FreeHiC/bin/quast"
REF_PAT = "/mnt/share1_Jabba/ftomas/reads/chr18/reference/chr18_PATERNAL.fasta"
REF_MAT = "/mnt/share1_Jabba/ftomas/reads/chr18/reference/chr18_MATERNAL.fasta"
OUT_BASE = BASE / "quast_phase_eval"
THREADS = "8"

TOOLS = [
    {
        "name": "GFAse_Shasta",
        "dir": BASE / "gfase_results",
        "phase0": "phase_0.fasta",
        "phase1": "phase_1.fasta",
    },
    {
        "name": "GFAse_Homology",
        "dir": BASE / "gfase_results_homology",
        "phase0": "phase_0.fasta",
        "phase1": "phase_1.fasta",
    },
    {
        "name": "MyTool_Shasta",
        "dir": BASE / "results_shasta",
        "phase0": "chr18_csv.phase_0.fasta",
        "phase1": "chr18_csv.phase_1.fasta",
    },
    {
        "name": "MyTool_MinHash",
        "dir": BASE / "results_minhash",
        "phase0": "chr18.phase_0.fasta",
        "phase1": "chr18.phase_1.fasta",
    },
    {
        "name": "MyTool_GFAsePairs",
        "dir": BASE / "results_gfase_pairs",
        "phase0": "chr18.phase_0.fasta",
        "phase1": "chr18.phase_1.fasta",
    },
]

REFS = [
    ("PATERNAL", REF_PAT),
    ("MATERNAL", REF_MAT),
]

TSV_KEYS = {
    "# contigs":             "contigs",
    "Total length":          "total_length",
    "Largest contig":        "largest_contig",
    "N50":                   "N50",
    "NG50":                  "NG50",
    "# misassemblies":       "misassemblies",
    "Genome fraction (%)":   "genome_fraction",
    "Duplication ratio":     "duplication_ratio",
    "# mismatches per 100 kbp": "mismatches_per_100k",
    "# indels per 100 kbp":  "indels_per_100k",
    "Total aligned length":  "aligned_length",
    "Reference length":      "ref_length",
}


def run_quast(fasta: Path, ref: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / "report.tsv"
    if report.exists():
        print(f"  [skip] {out_dir.name} (already done)")
        return
    cmd = [QUAST, "-r", str(ref), "-o", str(out_dir), "--threads", THREADS, str(fasta)]
    print(f"  [run] {' '.join(cmd)}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print("  [ERROR] quast failed:\n" + result.stderr.decode()[-1000:], file=sys.stderr)


def parse_report(tsv: Path) -> dict:
    metrics = {}
    if not tsv.exists():
        return metrics
    with open(tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            key = parts[0]
            if key in TSV_KEYS:
                try:
                    val = float(parts[1])
                except ValueError:
                    val = parts[1]
                metrics[TSV_KEYS[key]] = val
    return metrics


def fmt_len(n) -> str:
    try:
        n = float(n)
        return f"{n/1e6:.2f} Mb"
    except Exception:
        return str(n)


def fmt_val(v) -> str:
    try:
        f = float(v)
        if f == int(f):
            return str(int(f))
        return f"{f:.3f}"
    except Exception:
        return str(v)


def score_phase(m: dict) -> float:
    """Higher = better match to this reference."""
    gf = float(m.get("genome_fraction", 0))
    mm = float(m.get("mismatches_per_100k", 9999))
    ind = float(m.get("indels_per_100k", 9999))
    return gf - 0.1 * mm - 0.05 * ind


def determine_assignment(p0_pat, p0_mat, p1_pat, p1_mat):
    """Return (phase0_parent, phase1_parent)."""
    s_p0_pat = score_phase(p0_pat)
    s_p0_mat = score_phase(p0_mat)
    s_p1_pat = score_phase(p1_pat)
    s_p1_mat = score_phase(p1_mat)

    opt_A = s_p0_pat + s_p1_mat  # phase0=PAT, phase1=MAT
    opt_B = s_p0_mat + s_p1_pat  # phase0=MAT, phase1=PAT

    if opt_A >= opt_B:
        return "PATERNAL", "MATERNAL"
    else:
        return "MATERNAL", "PATERNAL"


def print_metrics_block(label: str, m: dict, indent: str = "    "):
    lines = []
    lines.append(f"{indent}{label}")
    lines.append(f"{indent}  Ukupna duljina sklopa : {fmt_len(m.get('total_length', 'N/A'))}  (referenca: {fmt_len(m.get('ref_length', 'N/A'))})")
    lines.append(f"{indent}  Broj kontiga         : {fmt_val(m.get('contigs', 'N/A'))}")
    lines.append(f"{indent}  Largest contig       : {fmt_len(m.get('largest_contig', 'N/A'))}")
    lines.append(f"{indent}  N50                  : {fmt_len(m.get('N50', 'N/A'))}")
    lines.append(f"{indent}  NG50                 : {fmt_len(m.get('NG50', 'N/A'))}")
    lines.append(f"{indent}  Genome fraction      : {fmt_val(m.get('genome_fraction', 'N/A'))}%")
    lines.append(f"{indent}  Total aligned length : {fmt_len(m.get('aligned_length', 'N/A'))}")
    lines.append(f"{indent}  Misassemblies        : {fmt_val(m.get('misassemblies', 'N/A'))}")
    lines.append(f"{indent}  Mismatches/100kbp    : {fmt_val(m.get('mismatches_per_100k', 'N/A'))}")
    lines.append(f"{indent}  Indels/100kbp        : {fmt_val(m.get('indels_per_100k', 'N/A'))}")
    lines.append(f"{indent}  Duplication ratio    : {fmt_val(m.get('duplication_ratio', 'N/A'))}")
    return lines


def main():
    # ── Step 1: run QUAST ──────────────────────────────────────────────────────
    print("=== Running QUAST (16 jobs) ===")
    for tool in TOOLS:
        for phase_key, phase_label in [("phase0", "phase_0"), ("phase1", "phase_1")]:
            fasta = tool["dir"] / tool[phase_key]
            if not fasta.exists():
                print(f"  [WARN] missing: {fasta}")
                continue
            for ref_label, ref_path in REFS:
                out_dir = OUT_BASE / tool["name"] / phase_label / ref_label
                run_quast(fasta, Path(ref_path), out_dir)

    # ── Step 2: collect results ────────────────────────────────────────────────
    print("\n=== Parsing results ===")
    results = {}
    for tool in TOOLS:
        tname = tool["name"]
        results[tname] = {}
        for phase_key, phase_label in [("phase0", "phase_0"), ("phase1", "phase_1")]:
            results[tname][phase_label] = {}
            for ref_label, _ in REFS:
                tsv = OUT_BASE / tname / phase_label / ref_label / "report.tsv"
                results[tname][phase_label][ref_label] = parse_report(tsv)

    # ── Step 3: build report ───────────────────────────────────────────────────
    lines = []
    lines.append("=" * 72)
    lines.append("QUAST PHASE EVALUATION — chr18")
    lines.append("=" * 72)
    lines.append("")

    for tool in TOOLS:
        tname = tool["name"]
        r = results[tname]

        p0_pat = r["phase_0"].get("PATERNAL", {})
        p0_mat = r["phase_0"].get("MATERNAL", {})
        p1_pat = r["phase_1"].get("PATERNAL", {})
        p1_mat = r["phase_1"].get("MATERNAL", {})

        if all([p0_pat, p0_mat, p1_pat, p1_mat]):
            assign0, assign1 = determine_assignment(p0_pat, p0_mat, p1_pat, p1_mat)
        else:
            assign0, assign1 = "?", "?"

        lines.append(f"{'─' * 72}")
        lines.append(f"TOOL: {tname}")
        lines.append(f"  Phase 0 → {assign0}")
        lines.append(f"  Phase 1 → {assign1}")
        lines.append("")

        for phase_label, phase_key in [("phase_0", "phase0"), ("phase_1", "phase1")]:
            lines.append(f"  [ {phase_label} ]")
            for ref_label, _ in REFS:
                m = r[phase_label].get(ref_label, {})
                if not m:
                    lines.append(f"    vs {ref_label}: NO DATA")
                    continue
                gf    = fmt_val(m.get("genome_fraction", "N/A"))
                al    = fmt_len(m.get("aligned_length", "N/A"))
                mm    = fmt_val(m.get("mismatches_per_100k", "N/A"))
                ind   = fmt_val(m.get("indels_per_100k", "N/A"))
                misas = fmt_val(m.get("misassemblies", "N/A"))
                lines.append(f"    vs {ref_label:<10}  Genome fraction: {gf:>7}%  Aligned: {al:>10}  Mismatches: {mm:>7}/100k  Indels: {ind:>7}/100k  Misassemblies: {misas}")
            lines.append("")

        lines.append("  Detailed metrics:")
        lines.append("")

        for phase_label in ["phase_0", "phase_1"]:
            for ref_label, _ in REFS:
                m = r[phase_label].get(ref_label, {})
                if not m:
                    continue
                label = f"{phase_label} vs {ref_label}"
                lines += print_metrics_block(label, m)
                lines.append("")

        lines.append("")

    report_text = "\n".join(lines)
    report_path = OUT_BASE / "phase_eval_report.txt"
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report_text)
    print(f"\nReport written to: {report_path}")
    print(report_text)


if __name__ == "__main__":
    main()
