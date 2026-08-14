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

# ── Anotação corrigida ───────────────────────────────────────────────────────
#
# Os nomes padronizados e — sobretudo — a lógica de correção da saída bruta do
# MITOS2 (frameshift do ND3, fragmentos espúrios, extremidades truncadas,
# região controle) vivem em `gff2genbank.py`. Este módulo NÃO deve reimplementar
# nada disso: manter duas cópias da mesma lógica já fez o `.gbk` divergir do
# `.tbl` de submissão.

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gff2genbank import (  # noqa: E402
    GENE_PRODUCT,
    TRNA_PRODUCT,
    RRNA_PRODUCT,
    base_gene_name,
    build_annotation,
    organism_to_seqid,
    validate_cds,
)

# Cores para o mapa circular (por categoria funcional)
COLOR_CDS_NADH  = "#4CAF50"   # Verde — Complex I (NADH dehydrogenase)
COLOR_CDS_COX   = "#2196F3"   # Azul — Complex IV (cytochrome c oxidase)
COLOR_CDS_ATP   = "#FF9800"   # Laranja — ATP synthase
COLOR_CDS_CYTB  = "#F44336"   # Vermelho — Complex III (cytochrome b)
COLOR_TRNA      = "#9C27B0"   # Roxo — tRNAs
COLOR_RRNA      = "#795548"   # Marrom — rRNAs
COLOR_DLOOP     = "#607D8B"   # Cinza — D-loop / control region

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


def build_genbank_record(ann, seqid, organism, transl_table, topology):
    """Monta o SeqRecord a partir da anotação já corrigida (`build_annotation`).

    Recebe o dicionário devolvido por `gff2genbank.build_annotation`, de modo
    que o `.gbk` e o mapa circular reflitam exatamente as mesmas coordenadas do
    `.tbl` de submissão.
    """
    sequence = ann["seq"]

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

    # ── D-loop / região controle ──
    # Pode vir em duas partes quando cruza a origem das coordenadas.
    dloop_parts = ann["dloop_parts"]
    if dloop_parts:
        locs = [FeatureLocation(a - 1, b, strand=1) for (a, b) in dloop_parts]
        dloop_loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)
        dloop_feat = SeqFeature(dloop_loc, type="D-loop")
        dloop_feat.qualifiers["note"] = ["control region"]
        record.features.append(dloop_feat)

    # ── CDS ──
    for f in ann["cds_feats"]:
        gene_sym, product = GENE_PRODUCT[base_gene_name(f["attrs"].get("Name", ""))]
        strand_val = -1 if f["strand"] == "-" else 1
        parts = f["parts"]

        # GFF é 1-based inclusivo; BioPython é 0-based semiaberto
        gene_loc = FeatureLocation(parts[0][0] - 1, parts[-1][1], strand=strand_val)
        gene_feat = SeqFeature(gene_loc, type="gene")
        gene_feat.qualifiers["gene"] = [gene_sym]
        record.features.append(gene_feat)

        locs = [FeatureLocation(a - 1, b, strand=strand_val) for (a, b) in parts]
        if strand_val == -1:
            locs = list(reversed(locs))
        cds_loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)

        cds_feat = SeqFeature(cds_loc, type="CDS")
        cds_feat.qualifiers["gene"] = [gene_sym]
        cds_feat.qualifiers["product"] = [product]
        cds_feat.qualifiers["transl_table"] = [transl_table]

        notes = []
        if f.get("partial_stop"):
            a, b = f["partial_stop"]
            pos = f"{a}" if a == b else f"{a}..{b}"
            if strand_val == -1:
                pos = f"complement({pos})"
            cds_feat.qualifiers["transl_except"] = [f"(pos:{pos},aa:TERM)"]
            notes.append("TAA stop codon is completed by the addition of "
                         "3' A residues to the mRNA")
        if f.get("frameshift"):
            cds_feat.qualifiers["exception"] = ["ribosomal slippage"]
            notes.append("programmed frameshift; frameshift mechanism unknown "
                         "(Mindell et al., 1998, Mol. Biol. Evol., 15:1568-1571)")
        if notes:
            cds_feat.qualifiers["note"] = notes

        record.features.append(cds_feat)

    # ── tRNA e rRNA ──
    for f in ann["rna_feats"]:
        name = base_gene_name(f["attrs"].get("Name", ""))
        strand_val = -1 if f["strand"] == "-" else 1
        loc = FeatureLocation(f["start"] - 1, f["end"], strand=strand_val)

        if name in TRNA_PRODUCT:
            gene_sym, product, anticodon = TRNA_PRODUCT[name]
            rna_type = "tRNA"
        elif name in RRNA_PRODUCT:
            gene_sym, product = RRNA_PRODUCT[name]
            anticodon = None
            rna_type = "rRNA"
        else:
            continue

        gene_feat = SeqFeature(loc, type="gene")
        gene_feat.qualifiers["gene"] = [gene_sym]
        record.features.append(gene_feat)

        rna_feat = SeqFeature(loc, type=rna_type)
        rna_feat.qualifiers["gene"] = [gene_sym]
        rna_feat.qualifiers["product"] = [product]
        if anticodon:
            rna_feat.qualifiers["note"] = [f"anticodon:{anticodon}"]
        record.features.append(rna_feat)

    record.features.sort(key=lambda x: (int(x.location.start), x.type != "gene"))
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

    seqid = args.seqid if args.seqid else organism_to_seqid(args.organism)

    # Anotação corrigida — mesma fonte que alimenta o .tbl de submissão
    ann = build_annotation(args.gff, args.fasta, args.transl_table)
    sequence = ann["seq"]

    record = build_genbank_record(
        ann, seqid, args.organism, args.transl_table, args.topology
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

    if ann["log"]:
        print(f"\n  Correções aplicadas sobre a anotação bruta do MITOS2:")
        for item in ann["log"]:
            print(f"    · {item}")

    problems = validate_cds(ann["cds_feats"], sequence)
    if problems:
        print(f"\n  ATENÇÃO — problemas remanescentes:")
        for p in problems:
            print(f"    ! {p}")

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
