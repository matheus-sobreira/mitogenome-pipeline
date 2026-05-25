#!/usr/bin/env python3
"""
gff2genbank.py — Converte saída do MITOS2 (GFF3 + FASTA) para formato GenBank
(feature table .tbl + FASTA .fsa) pronto para submissão via table2asn.

Trata automaticamente o frameshift do ND3 em mitogenomas aviários,
unificando nad3_0 e nad3_1 em um único CDS com join().

Uso:
    python3 gff2genbank.py \\
        --gff result.gff \\
        --fasta assembly.fasta \\
        --organism "Anodorhynchus leari" \\
        --outdir genbank_submission/
"""

import argparse
import os
import re
import sys


def organism_to_seqid(organism):
    """Convert organism name to a short sequence ID. e.g. 'Anodorhynchus leari' → 'A_leari'."""
    parts = organism.strip().split()
    if len(parts) >= 2:
        return f"{parts[0][0]}_{parts[1]}"
    return organism.replace(" ", "_")
from collections import OrderedDict

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
    "trnL2": ("trnL2", "tRNA-Leu"),  # UUR
    "trnL1": ("trnL1", "tRNA-Leu"),  # CUN
    "trnI":  ("trnI",  "tRNA-Ile"),
    "trnQ":  ("trnQ",  "tRNA-Gln"),
    "trnM":  ("trnM",  "tRNA-Met"),
    "trnW":  ("trnW",  "tRNA-Trp"),
    "trnA":  ("trnA",  "tRNA-Ala"),
    "trnN":  ("trnN",  "tRNA-Asn"),
    "trnC":  ("trnC",  "tRNA-Cys"),
    "trnY":  ("trnY",  "tRNA-Tyr"),
    "trnS2": ("trnS2", "tRNA-Ser"),  # UCN
    "trnS1": ("trnS1", "tRNA-Ser"),  # AGY
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
    "rrnS": ("rrnS", "12S ribosomal RNA"),
    "rrnL": ("rrnL", "16S ribosomal RNA"),
}


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


def merge_nad3(features, transl_table=2):
    """Detect nad3_0 and nad3_1, merge into single nad3 with frameshift join.

    Em vertebrados (transl_table=2), o MITOS2 frequentemente reporta as duas
    ORFs do ND3 com sobreposição de poucos nucleotídeos. A convenção do
    GenBank para o programmed ribosomal frameshift descrito por Mindell et
    al. (1998) --- observado em aves, tartarugas, crocodilianos e demais
    vertebrados não-mamíferos --- é representar a junção como skip de 1 nt
    (gap de 1 base entre os dois exons). Quando essa sobreposição é
    detectada e o código genético é o vertebrado mitocondrial, as
    coordenadas são ajustadas para o skip canônico antes da emissão do .tbl.
    """
    nad3_parts = {}
    other = []
    for f in features:
        name = f["attrs"].get("Name", "")
        if name in ("nad3_0", "nad3_1") and f["type"] == "gene":
            nad3_parts[name] = f
        elif name.startswith("nad3_") and f["type"] == "exon":
            nad3_parts[f"exon_{name}"] = f
        else:
            other.append(f)

    if "nad3_0" not in nad3_parts or "nad3_1" not in nad3_parts:
        return features, None

    p0 = nad3_parts["nad3_0"]
    p1 = nad3_parts["nad3_1"]

    # nad3_1 is typically the first (upstream) part
    first = p1 if p1["start"] <= p0["start"] else p0
    second = p0 if p1["start"] <= p0["start"] else p1

    raw_coords = (first["start"], first["end"], second["start"], second["end"])

    # Ajuste para skip-1 canônico (Mindell 1998) em vertebrados quando há overlap.
    overlap = first["end"] - second["start"] + 1
    if transl_table == 2 and overlap > 0:
        adjusted_first_end = second["start"] - 2
        join_coords = (first["start"], adjusted_first_end,
                       second["start"], second["end"])
        print(f"  ND3 overlap de {overlap} bp detectado; ajustado para skip-1 "
              f"canônico (Mindell 1998): "
              f"join({join_coords[0]}..{join_coords[1]},"
              f"{join_coords[2]}..{join_coords[3]})")
    else:
        join_coords = raw_coords

    merged = {
        "seqid":  first["seqid"],
        "source": "mitos",
        "type":   "gene",
        "start":  first["start"],
        "end":    second["end"],
        "strand": first["strand"],
        "phase":  "0",
        "attrs":  {"Name": "nad3", "ID": "gene_nad3"},
        "frameshift": True,
        "join_coords": join_coords,
    }

    # Insert merged nad3 and sort all by start position
    result = other + [merged]
    result.sort(key=lambda x: x["start"])

    return result, merged


def write_tbl(features, seqid, genome_length, transl_table, outpath):
    """Write GenBank feature table (.tbl) in NCBI 5-column format."""
    with open(outpath, "w") as out:
        out.write(f">Feature {seqid}\n")

        for f in features:
            name = f["attrs"].get("Name", "")
            ftype = f["type"]

            if ftype not in ("gene", "ncRNA_gene", "origin_of_replication"):
                continue

            base_name = re.sub(r"_\d+$", "", name)

            # ── Protein-coding gene ──
            if base_name in GENE_PRODUCT:
                gene_sym, product = GENE_PRODUCT[base_name]
                is_complement = f["strand"] == "-"

                if f.get("frameshift"):
                    c = f["join_coords"]  # (start1, end1, start2, end2)

                    # gene spans full range
                    if is_complement:
                        out.write(f"{f['end']}\t{f['start']}\tgene\n")
                    else:
                        out.write(f"{f['start']}\t{f['end']}\tgene\n")
                    out.write(f"\t\t\tgene\t{gene_sym}\n")

                    # CDS with join — first line has feature key, continuation lines don't
                    if is_complement:
                        out.write(f"{c[3]}\t{c[2]}\tCDS\n")
                        out.write(f"{c[1]}\t{c[0]}\n")
                    else:
                        out.write(f"{c[0]}\t{c[1]}\tCDS\n")
                        out.write(f"{c[2]}\t{c[3]}\n")
                    out.write(f"\t\t\tgene\t{gene_sym}\n")
                    out.write(f"\t\t\tproduct\t{product}\n")
                    out.write(f"\t\t\ttransl_table\t{transl_table}\n")
                    out.write(f"\t\t\texception\tribosomal slippage\n")
                    out.write(f'\t\t\tnote\tprogrammed frameshift; frameshift mechanism unknown (Mindell et al., 1998, Mol. Biol. Evol., 15:1568-1571)\n')

                else:
                    # Normal CDS
                    if is_complement:
                        out.write(f"{f['end']}\t{f['start']}\tgene\n")
                        out.write(f"\t\t\tgene\t{gene_sym}\n")
                        out.write(f"{f['end']}\t{f['start']}\tCDS\n")
                    else:
                        out.write(f"{f['start']}\t{f['end']}\tgene\n")
                        out.write(f"\t\t\tgene\t{gene_sym}\n")
                        out.write(f"{f['start']}\t{f['end']}\tCDS\n")
                    out.write(f"\t\t\tproduct\t{product}\n")
                    out.write(f"\t\t\ttransl_table\t{transl_table}\n")

            # ── tRNA ──
            elif base_name in TRNA_PRODUCT:
                gene_sym, product = TRNA_PRODUCT[base_name]
                is_complement = f["strand"] == "-"
                if is_complement:
                    out.write(f"{f['end']}\t{f['start']}\tgene\n")
                    out.write(f"\t\t\tgene\t{gene_sym}\n")
                    out.write(f"{f['end']}\t{f['start']}\ttRNA\n")
                else:
                    out.write(f"{f['start']}\t{f['end']}\tgene\n")
                    out.write(f"\t\t\tgene\t{gene_sym}\n")
                    out.write(f"{f['start']}\t{f['end']}\ttRNA\n")
                out.write(f"\t\t\tproduct\t{product}\n")

            # ── rRNA ──
            elif base_name in RRNA_PRODUCT:
                gene_sym, product = RRNA_PRODUCT[base_name]
                out.write(f"{f['start']}\t{f['end']}\tgene\n")
                out.write(f"\t\t\tgene\t{gene_sym}\n")
                out.write(f"{f['start']}\t{f['end']}\trRNA\n")
                out.write(f"\t\t\tproduct\t{product}\n")

            # ── D-loop / origin of replication ──
            elif name.startswith("OH"):
                out.write(f"{f['start']}\t{f['end']}\tD-loop\n")


def write_fsa(fasta_path, seqid, organism, topology, outpath):
    """Write GenBank-formatted FASTA (.fsa) with proper header."""
    seq_lines = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_lines.append(line.strip())
    sequence = "".join(seq_lines)

    with open(outpath, "w") as out:
        header = (
            f">{seqid} "
            f"[organism={organism}] "
            f"[location=mitochondrion] "
            f"[topology={topology}] "
            f"[completeness=complete] "
            f"[molecule_type=genomic DNA] "
            f"[gcode=2]"
        )
        out.write(header + "\n")
        # Write 70-char lines
        for i in range(0, len(sequence), 70):
            out.write(sequence[i:i+70] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Converte MITOS2 GFF3 → GenBank feature table (.tbl) + FASTA (.fsa)"
    )
    parser.add_argument("--gff", required=True, help="MITOS2 result.gff")
    parser.add_argument("--fasta", required=True, help="Assembly FASTA (montagem original)")
    parser.add_argument("--organism", required=True, help='Nome do organismo (ex: "Anodorhynchus leari")')
    parser.add_argument("--seqid", default=None, help="Sequence ID (default: from FASTA header)")
    parser.add_argument("--transl-table", type=int, default=2, help="Genetic code (default: 2, vertebrate mitochondrial)")
    parser.add_argument("--topology", default="circular", choices=["circular", "linear"])
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Get seqid from organism name if not provided
    seqid = args.seqid
    if not seqid:
        seqid = organism_to_seqid(args.organism)

    # Parse GFF
    features = parse_gff(args.gff)

    # Get genome length from region feature
    genome_length = 0
    for f in features:
        if f["type"] == "region":
            genome_length = f["end"]
            break
    if genome_length == 0:
        # Fallback: read FASTA
        with open(args.fasta) as fh:
            seq = ""
            for line in fh:
                if not line.startswith(">"):
                    seq += line.strip()
            genome_length = len(seq)

    # Keep only gene-level and origin_of_replication features
    gene_features = [f for f in features if f["type"] in ("gene", "ncRNA_gene", "origin_of_replication")]

    # Merge nad3 frameshift
    gene_features, nad3_merged = merge_nad3(gene_features, transl_table=args.transl_table)
    if nad3_merged:
        print(f"  ND3 frameshift detectado: join({nad3_merged['join_coords'][0]}..{nad3_merged['join_coords'][1]},"
              f"{nad3_merged['join_coords'][2]}..{nad3_merged['join_coords'][3]})")

    # Write outputs
    tbl_path = os.path.join(args.outdir, f"{seqid}.tbl")
    fsa_path = os.path.join(args.outdir, f"{seqid}.fsa")

    write_tbl(gene_features, seqid, genome_length, args.transl_table, tbl_path)
    write_fsa(args.fasta, seqid, args.organism, args.topology, fsa_path)

    # Summary
    cds_count = sum(1 for f in gene_features
                    if f["type"] in ("gene",) and re.sub(r"_\d+$", "", f["attrs"].get("Name", "")) in GENE_PRODUCT)
    trna_count = sum(1 for f in gene_features
                     if f["type"] in ("gene", "ncRNA_gene") and re.sub(r"_\d+$", "", f["attrs"].get("Name", "")) in TRNA_PRODUCT)
    rrna_count = sum(1 for f in gene_features
                     if f["type"] in ("gene", "ncRNA_gene") and re.sub(r"_\d+$", "", f["attrs"].get("Name", "")) in RRNA_PRODUCT)

    print(f"\n  Organismo:    {args.organism}")
    print(f"  Comprimento: {genome_length} bp ({args.topology})")
    print(f"  Genes CDS:   {cds_count}")
    print(f"  tRNAs:        {trna_count}")
    print(f"  rRNAs:        {rrna_count}")
    print(f"  Código gen.:  {args.transl_table} (mitocôndria de vertebrados)")
    print(f"\n  Arquivos gerados:")
    print(f"    {tbl_path}")
    print(f"    {fsa_path}")
    print(f"\n  Para gerar .sqn (submissão):")
    print(f"    table2asn -i {fsa_path} -f {tbl_path} -o {seqid}.sqn")


if __name__ == "__main__":
    main()
