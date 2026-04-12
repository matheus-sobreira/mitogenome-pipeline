#!/usr/bin/env python3
"""
generate_genbank.py — Gera GenBank Flat File (.gbk) e mapa circular (SVG/PDF)
a partir da saída do MITOS2 (GFF3 + FASTA).

Produz:
  - arquivo .gbk no formato GenBank Flat File (LOCUS, FEATURES, ORIGIN)
  - mapa circular do mitogenoma em SVG e PDF (via Biopython GenomeDiagram)

O .gbk pode ser:
  - Carregado no OGDRAW para mapa de publicação
  - Submetido ao GenBank via BankIt
  - Usado como entrada para análises downstream

Uso:
    python3 generate_genbank.py \\
        --gff result.gff \\
        --fasta assembly.fasta \\
        --organism "Anodorhynchus leari" \\
        --outdir genbank_submission/
"""

import argparse
import os
import re
import sys
from datetime import date

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# ── Nomes padronizados para GenBank ──────────────────────────────────────────

GENE_PRODUCT = {
    "nad1":  ("ND1",  "NADH dehydrogenase subunit 1"),
    "nad2":  ("ND2",  "NADH dehydrogenase subunit 2"),
    "nad3":  ("ND3",  "NADH dehydrogenase subunit 3"),
    "nad4":  ("ND4",  "NADH dehydrogenase subunit 4"),
    "nad4l": ("ND4L", "NADH dehydrogenase subunit 4L"),
    "nad5":  ("ND5",  "NADH dehydrogenase subunit 5"),
    "nad6":  ("ND6",  "NADH dehydrogenase subunit 6"),
    "cox1":  ("COX1", "cytochrome c oxidase subunit I"),
    "cox2":  ("COX2", "cytochrome c oxidase subunit II"),
    "cox3":  ("COX3", "cytochrome c oxidase subunit III"),
    "atp6":  ("ATP6", "ATP synthase F0 subunit 6"),
    "atp8":  ("ATP8", "ATP synthase F0 subunit 8"),
    "cob":   ("CYTB", "cytochrome b"),
}

TRNA_PRODUCT = {
    "trnF":  ("trnF",  "tRNA-Phe"),
    "trnV":  ("trnV",  "tRNA-Val"),
    "trnL2": ("trnL2", "tRNA-Leu"),
    "trnL1": ("trnL1", "tRNA-Leu"),
    "trnI":  ("trnI",  "tRNA-Ile"),
    "trnQ":  ("trnQ",  "tRNA-Gln"),
    "trnM":  ("trnM",  "tRNA-Met"),
    "trnW":  ("trnW",  "tRNA-Trp"),
    "trnA":  ("trnA",  "tRNA-Ala"),
    "trnN":  ("trnN",  "tRNA-Asn"),
    "trnC":  ("trnC",  "tRNA-Cys"),
    "trnY":  ("trnY",  "tRNA-Tyr"),
    "trnS2": ("trnS2", "tRNA-Ser"),
    "trnS1": ("trnS1", "tRNA-Ser"),
    "trnD":  ("trnD",  "tRNA-Asp"),
    "trnK":  ("trnK",  "tRNA-Lys"),
    "trnG":  ("trnG",  "tRNA-Gly"),
    "trnR":  ("trnR",  "tRNA-Arg"),
    "trnH":  ("trnH",  "tRNA-His"),
    "trnT":  ("trnT",  "tRNA-Thr"),
    "trnP":  ("trnP",  "tRNA-Pro"),
    "trnE":  ("trnE",  "tRNA-Glu"),
}

RRNA_PRODUCT = {
    "rrnS": ("rrnS", "s-rRNA", "12S ribosomal RNA"),
    "rrnL": ("rrnL", "l-rRNA", "16S ribosomal RNA"),
}

# Cores para o mapa circular (por categoria funcional)
COLOR_CDS_NADH  = "#4CAF50"   # Verde — Complex I (NADH dehydrogenase)
COLOR_CDS_COX   = "#2196F3"   # Azul — Complex IV (cytochrome c oxidase)
COLOR_CDS_ATP   = "#FF9800"   # Laranja — ATP synthase
COLOR_CDS_CYTB  = "#F44336"   # Vermelho — Complex III (cytochrome b)
COLOR_TRNA      = "#9C27B0"   # Roxo — tRNAs
COLOR_RRNA      = "#795548"   # Marrom — rRNAs
COLOR_DLOOP     = "#607D8B"   # Cinza — D-loop / control region

def organism_to_seqid(organism):
    """Convert organism name to a short sequence ID. e.g. 'Anodorhynchus leari' → 'A_leari'."""
    parts = organism.strip().split()
    if len(parts) >= 2:
        return f"{parts[0][0]}_{parts[1]}"
    return organism.replace(" ", "_")


def get_cds_color(gene_symbol):
    """Return color based on gene function."""
    if gene_symbol.startswith("ND"):
        return COLOR_CDS_NADH
    elif gene_symbol.startswith("COX"):
        return COLOR_CDS_COX
    elif gene_symbol.startswith("ATP"):
        return COLOR_CDS_ATP
    elif gene_symbol == "CYTB":
        return COLOR_CDS_CYTB
    return "#9E9E9E"


def parse_gff(gff_path):
    """Parse MITOS2 GFF3, returning a list of feature dicts."""
    features = []
    with open(gff_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs_str = cols
            attrs = {}
            for pair in attrs_str.split(";"):
                if "=" in pair:
                    k, v = pair.split("=", 1)
                    attrs[k] = v
            features.append({
                "seqid":  seqid,
                "source": source,
                "type":   ftype,
                "start":  int(start),
                "end":    int(end),
                "strand": strand,
                "phase":  phase,
                "attrs":  attrs,
            })
    return features


def read_fasta(fasta_path):
    """Read FASTA and return (seqid, sequence)."""
    seqid = None
    seq_parts = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                seqid = line[1:].split()[0]
            else:
                seq_parts.append(line)
    return seqid, "".join(seq_parts)


def build_genbank_record(gff_features, seqid, sequence, organism, transl_table, topology):
    """Build a BioPython SeqRecord with GenBank features from MITOS2 GFF3."""

    record = SeqRecord(
        Seq(sequence),
        id=seqid,
        name=seqid,
        description=f"{organism} mitochondrion, complete genome",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = topology
    record.annotations["data_file_division"] = "VRT"
    record.annotations["date"] = date.today().strftime("%d-%b-%Y").upper()
    record.annotations["organism"] = organism
    record.annotations["source"] = f"mitochondrion {organism}"
    record.annotations["taxonomy"] = []
    record.annotations["references"] = []

    # Source feature
    source_loc = FeatureLocation(0, len(sequence), strand=1)
    source_feat = SeqFeature(source_loc, type="source")
    source_feat.qualifiers["organism"] = [organism]
    source_feat.qualifiers["organelle"] = ["mitochondrion"]
    source_feat.qualifiers["mol_type"] = ["genomic DNA"]
    record.features.append(source_feat)

    # Collect gene-level features
    gene_features = [f for f in gff_features
                     if f["type"] in ("gene", "ncRNA_gene", "origin_of_replication")]

    # Detect nad3 frameshift
    nad3_parts = {}
    for f in gene_features:
        name = f["attrs"].get("Name", "")
        if name in ("nad3_0", "nad3_1"):
            nad3_parts[name] = f

    # Process each feature
    processed_nad3 = False
    for f in sorted(gene_features, key=lambda x: x["start"]):
        name = f["attrs"].get("Name", "")
        base_name = re.sub(r"_\d+$", "", name)
        strand_val = -1 if f["strand"] == "-" else 1

        # GFF is 1-based inclusive; BioPython is 0-based half-open
        start_0 = f["start"] - 1
        end_0 = f["end"]

        # ── nad3 frameshift ──
        if name in ("nad3_0", "nad3_1"):
            if processed_nad3:
                continue
            if "nad3_0" in nad3_parts and "nad3_1" in nad3_parts:
                processed_nad3 = True
                p0 = nad3_parts["nad3_0"]
                p1 = nad3_parts["nad3_1"]
                first = p1 if p1["start"] <= p0["start"] else p0
                second = p0 if p1["start"] <= p0["start"] else p1

                # Gene spanning full range
                gene_loc = FeatureLocation(first["start"] - 1, second["end"], strand=strand_val)
                gene_feat = SeqFeature(gene_loc, type="gene")
                gene_feat.qualifiers["gene"] = ["ND3"]
                record.features.append(gene_feat)

                # CDS with join (CompoundLocation)
                loc1 = FeatureLocation(first["start"] - 1, first["end"], strand=strand_val)
                loc2 = FeatureLocation(second["start"] - 1, second["end"], strand=strand_val)
                if strand_val == -1:
                    cds_loc = CompoundLocation([loc2, loc1])
                else:
                    cds_loc = CompoundLocation([loc1, loc2])
                cds_feat = SeqFeature(cds_loc, type="CDS")
                cds_feat.qualifiers["gene"] = ["ND3"]
                cds_feat.qualifiers["product"] = ["NADH dehydrogenase subunit 3"]
                cds_feat.qualifiers["transl_table"] = [transl_table]
                cds_feat.qualifiers["exception"] = ["ribosomal slippage"]
                cds_feat.qualifiers["note"] = [
                    "programmed frameshift; frameshift mechanism unknown "
                    "(Mindell et al., 1998, Mol. Biol. Evol., 15:1568-1571)"
                ]
                record.features.append(cds_feat)
                continue

        # ── CDS (protein-coding) ──
        if base_name in GENE_PRODUCT:
            gene_sym, product = GENE_PRODUCT[base_name]
            loc = FeatureLocation(start_0, end_0, strand=strand_val)

            gene_feat = SeqFeature(loc, type="gene")
            gene_feat.qualifiers["gene"] = [gene_sym]
            record.features.append(gene_feat)

            cds_feat = SeqFeature(loc, type="CDS")
            cds_feat.qualifiers["gene"] = [gene_sym]
            cds_feat.qualifiers["product"] = [product]
            cds_feat.qualifiers["transl_table"] = [transl_table]
            record.features.append(cds_feat)

        # ── tRNA ──
        elif base_name in TRNA_PRODUCT:
            gene_sym, product = TRNA_PRODUCT[base_name]
            loc = FeatureLocation(start_0, end_0, strand=strand_val)

            gene_feat = SeqFeature(loc, type="gene")
            gene_feat.qualifiers["gene"] = [gene_sym]
            record.features.append(gene_feat)

            trna_feat = SeqFeature(loc, type="tRNA")
            trna_feat.qualifiers["gene"] = [gene_sym]
            trna_feat.qualifiers["product"] = [product]
            record.features.append(trna_feat)

        # ── rRNA ──
        elif base_name in RRNA_PRODUCT:
            gene_sym, display, product = RRNA_PRODUCT[base_name]
            loc = FeatureLocation(start_0, end_0, strand=strand_val)

            gene_feat = SeqFeature(loc, type="gene")
            gene_feat.qualifiers["gene"] = [gene_sym]
            record.features.append(gene_feat)

            rrna_feat = SeqFeature(loc, type="rRNA")
            rrna_feat.qualifiers["gene"] = [gene_sym]
            rrna_feat.qualifiers["product"] = [product]
            record.features.append(rrna_feat)

        # ── D-loop / origin of replication ──
        elif name.startswith("OH"):
            loc = FeatureLocation(start_0, end_0, strand=strand_val)
            dloop_feat = SeqFeature(loc, type="D-loop")
            dloop_feat.qualifiers["note"] = ["control region"]
            record.features.append(dloop_feat)

    return record


def draw_circular_map(record, outdir, organism):
    """Draw a circular genome map using Biopython GenomeDiagram."""
    from Bio.Graphics import GenomeDiagram
    from reportlab.lib import colors
    from reportlab.lib.units import cm

    def hex_to_color(hex_str):
        hex_str = hex_str.lstrip("#")
        r, g, b = int(hex_str[0:2], 16), int(hex_str[2:4], 16), int(hex_str[4:6], 16)
        return colors.Color(r/255, g/255, b/255)

    diagram = GenomeDiagram.Diagram(organism)

    # Track 1: Genes (CDS + tRNA + rRNA + D-loop)
    gene_track = diagram.new_track(1, name="Genes", greytrack=False,
                                    scale=True, scale_ticks=True,
                                    scale_smalltick_interval=1000,
                                    scale_largetick_interval=5000,
                                    scale_smalltick_labels=False,
                                    scale_largetick_labels=True,
                                    scale_fontsize=6)
    gene_set = gene_track.new_set()

    for feature in record.features:
        if feature.type == "source":
            continue

        # Only draw functional features (skip gene type, use CDS/tRNA/rRNA/D-loop)
        if feature.type == "gene":
            continue

        name = feature.qualifiers.get("gene", feature.qualifiers.get("product", [""]))[0]
        label_text = name

        if feature.type == "CDS":
            gene_sym = feature.qualifiers.get("gene", [""])[0]
            col = hex_to_color(get_cds_color(gene_sym))
            label_text = gene_sym
        elif feature.type == "tRNA":
            col = hex_to_color(COLOR_TRNA)
            product = feature.qualifiers.get("product", [""])[0]
            label_text = product  # e.g. "tRNA-Phe"
        elif feature.type == "rRNA":
            col = hex_to_color(COLOR_RRNA)
            product = feature.qualifiers.get("product", [""])[0]
            if "12S" in product:
                label_text = "12S"
            elif "16S" in product:
                label_text = "16S"
        elif feature.type == "D-loop":
            col = hex_to_color(COLOR_DLOOP)
            label_text = "D-loop"
        else:
            continue

        gene_set.add_feature(
            feature,
            color=col,
            border=colors.black,
            label=True,
            label_size=5,
            label_angle=0,
            name=label_text,
            label_position="middle",
        )

    # Draw circular
    diagram.draw(
        format="circular",
        circular=True,
        pagesize=(25*cm, 25*cm),
        start=0,
        end=len(record),
        circle_core=0.3,
        tracklines=False,
    )

    svg_path = os.path.join(outdir, "circular_map.svg")
    pdf_path = os.path.join(outdir, "circular_map.pdf")
    diagram.write(svg_path, "SVG")
    diagram.write(pdf_path, "PDF")

    return svg_path, pdf_path


def main():
    parser = argparse.ArgumentParser(
        description="Gera GenBank Flat File (.gbk) e mapa circular a partir do MITOS2"
    )
    parser.add_argument("--gff", required=True, help="MITOS2 result.gff")
    parser.add_argument("--fasta", required=True, help="Assembly FASTA")
    parser.add_argument("--organism", required=True,
                        help='Nome do organismo (ex: "Anodorhynchus leari")')
    parser.add_argument("--seqid", default=None,
                        help="Sequence ID (default: from FASTA header)")
    parser.add_argument("--transl-table", type=int, default=2,
                        help="Genetic code (default: 2, vertebrate mitochondrial)")
    parser.add_argument("--topology", default="circular",
                        choices=["circular", "linear"])
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--no-map", action="store_true",
                        help="Skip circular map generation")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Read sequence
    _, sequence = read_fasta(args.fasta)
    seqid = args.seqid if args.seqid else organism_to_seqid(args.organism)

    # Parse GFF
    gff_features = parse_gff(args.gff)

    # Build GenBank record
    record = build_genbank_record(
        gff_features, seqid, sequence,
        args.organism, args.transl_table, args.topology
    )

    # Write GenBank flat file
    gbk_path = os.path.join(args.outdir, f"{seqid}.gbk")
    with open(gbk_path, "w") as out:
        SeqIO.write(record, out, "genbank")

    # Count features
    cds_count = sum(1 for f in record.features if f.type == "CDS")
    trna_count = sum(1 for f in record.features if f.type == "tRNA")
    rrna_count = sum(1 for f in record.features if f.type == "rRNA")
    dloop_count = sum(1 for f in record.features if f.type == "D-loop")

    print(f"\n  GenBank Flat File: {gbk_path}")
    print(f"  Organismo:    {args.organism}")
    print(f"  Comprimento:  {len(sequence):,} bp ({args.topology})")
    print(f"  CDS:          {cds_count}")
    print(f"  tRNAs:        {trna_count}")
    print(f"  rRNAs:        {rrna_count}")
    print(f"  D-loop:       {dloop_count}")

    # Draw circular map
    if not args.no_map:
        try:
            svg_path, pdf_path = draw_circular_map(record, args.outdir, args.organism)
            print(f"  Mapa circular: {svg_path}")
            print(f"  Mapa circular: {pdf_path}")
        except Exception as e:
            print(f"  AVISO: Não foi possível gerar mapa circular: {e}", file=sys.stderr)

    print()


if __name__ == "__main__":
    main()
