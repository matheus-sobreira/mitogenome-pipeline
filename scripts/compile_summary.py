#!/usr/bin/env python3
"""
compile_summary.py — Compila os resultados do pipeline em uma pasta
organizada com todos os entregáveis para apresentação/submissão.

Entregáveis gerados:
  01_genome_assembly.fasta        — Genoma mitocondrial circularizado
  02_coding_genes_nt.fasta        — Genes codificantes (nucleotídeos)
  03_coding_genes_aa.fasta        — Genes codificantes (aminoácidos)
  04_ribosomal_genes.fasta        — Genes ribossomais (12S, 16S)
  05_transport_genes.fasta        — Genes transportadores (tRNAs)
  06_gene_positions.tsv           — Tabela de posição física dos genes
  07_start_stop_codons.tsv        — Tabela de códons start/stop (CDS)
  08_trna_anticodons.tsv          — Tabela de anticódons dos tRNAs
  09_genome_map.png               — Mapa linear do genoma (MITOS2)
  10_annotation.gff               — Anotação GFF3 completa
  structure_svgs/                  — Estruturas secundárias (tRNA + rRNA)
  genbank_submission/              — Arquivos para submissão ao GenBank

Uso:
    python3 compile_summary.py \\
        --assembly  assembly.fasta \\
        --mitos-dir mitos2_output/ \\
        --genbank-dir genbank_submission/ \\
        --organism "Anodorhynchus leari" \\
        --outdir summary/
"""

import argparse
import os
import re
import shutil
import sys
from collections import OrderedDict

# ── Mapeamento de nomes MITOS2 → nomes padronizados ─────────────────────────

CDS_GENES = {
    "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6",
    "cox1", "cox2", "cox3", "atp6", "atp8", "cob",
}

# nad3_0 e nad3_1 = frameshift do mesmo gene nad3
CDS_ALIASES = {"nad3_0": "nad3", "nad3_1": "nad3"}

RRNA_GENES = {"rrnS", "rrnL"}

TRNA_GENES = {
    "trnF", "trnV", "trnL2", "trnL1", "trnI", "trnQ", "trnM",
    "trnW", "trnA", "trnN", "trnC", "trnY", "trnS2", "trnS1",
    "trnD", "trnK", "trnG", "trnR", "trnH", "trnT", "trnP", "trnE",
}

GENE_DISPLAY_NAMES = {
    "nad1": "ND1", "nad2": "ND2", "nad3": "ND3", "nad4": "ND4",
    "nad4l": "ND4L", "nad5": "ND5", "nad6": "ND6",
    "cox1": "COX1", "cox2": "COX2", "cox3": "COX3",
    "atp6": "ATP6", "atp8": "ATP8", "cob": "CYTB",
    "rrnS": "12S rRNA", "rrnL": "16S rRNA",
}

TRNA_AA_NAMES = {
    "trnF": "Phe", "trnV": "Val", "trnL2": "Leu(UUR)", "trnL1": "Leu(CUN)",
    "trnI": "Ile", "trnQ": "Gln", "trnM": "Met", "trnW": "Trp",
    "trnA": "Ala", "trnN": "Asn", "trnC": "Cys", "trnY": "Tyr",
    "trnS2": "Ser(UCN)", "trnS1": "Ser(AGY)", "trnD": "Asp",
    "trnK": "Lys", "trnG": "Gly", "trnR": "Arg", "trnH": "His",
    "trnT": "Thr", "trnP": "Pro", "trnE": "Glu",
}

# Tabela genética mitocondrial de vertebrados (transl_table=2)
MT_CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Met", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "*",   "TAG": "*",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Trp", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "*",   "AGG": "*",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

MT_START_CODONS = {"ATG", "ATA", "ATT", "ATC", "GTG"}
MT_STOP_CODONS = {"TAA", "TAG", "AGA", "AGG", "TA", "T"}


def parse_fasta_multi(path):
    """Parse multi-FASTA, return list of (header, sequence) tuples."""
    records = []
    header = None
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def parse_mitos_tab(path):
    """Parse result.mitos tabular file."""
    entries = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 10:
                continue
            entries.append({
                "seqid":     cols[0],
                "feat_type": cols[1],
                "name":      cols[2],
                "source":    cols[3],
                "start":     int(cols[4]) + 1,  # BED is 0-based, convert to 1-based
                "end":       int(cols[5]),
                "strand":    "+" if int(cols[6]) == 1 else "-",
                "score":     cols[7],
                "bit_score": cols[8],
                "anticodon": cols[9] if cols[9] != "-" else None,
                "length":    cols[10] if len(cols) > 10 else None,
            })
    return entries


def parse_bed(path):
    """Parse result.bed."""
    entries = []
    with open(path) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            name_raw = cols[3]
            name = name_raw.split("(")[0]  # Remove anticodon e.g. trnF(gaa) → trnF
            entries.append({
                "name_raw": name_raw,
                "name":     name,
                "start":    int(cols[1]) + 1,  # BED is 0-based
                "end":      int(cols[2]),
                "strand":   cols[5],
                "length":   int(cols[2]) - int(cols[1]),
            })
    return entries


def get_gene_name(header):
    """Extract gene name from MITOS2 FASTA header like '>Contig1; 4181-5154; +; nad1'."""
    parts = header.split(";")
    if len(parts) >= 4:
        return parts[-1].strip()
    return header.split()[-1]


def extract_start_stop(seq, strand):
    """Extract start and stop codons from a nucleotide CDS sequence."""
    seq = seq.upper().replace("U", "T")
    if len(seq) < 6:
        return ("?", "?")

    start_codon = seq[:3]
    # Stop codon: check for truncated (T, TA) or full (TAA, TAG, AGA, AGG)
    stop_codon = seq[-3:]
    if stop_codon in MT_STOP_CODONS:
        pass
    elif seq[-2:] in ("TA",):
        stop_codon = seq[-2:] + "(+A)"  # Polyadenylation completes TAA
    elif seq[-1:] == "T":
        stop_codon = "T(+AA)"  # Polyadenylation completes TAA
    return (start_codon, stop_codon)


def write_fasta(records, outpath, wrap=70):
    """Write list of (header, sequence) to FASTA."""
    with open(outpath, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                out.write(seq[i:i+wrap] + "\n")


def compile_summary(args):

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    mitos_dir = args.mitos_dir
    assembly  = args.assembly

    # ── 01. Genoma circularizado ─────────────────────────────────────────
    out_assembly = os.path.join(outdir, "01_genome_assembly.fasta")
    shutil.copy2(assembly, out_assembly)
    print(f"  [01] Genoma circularizado → {out_assembly}")

    # ── Parse saídas MITOS2 ──────────────────────────────────────────────
    fas_path = os.path.join(mitos_dir, "result.fas")
    faa_path = os.path.join(mitos_dir, "result.faa")
    bed_path = os.path.join(mitos_dir, "result.bed")
    mitos_path = os.path.join(mitos_dir, "result.mitos")
    gff_path = os.path.join(mitos_dir, "result.gff")
    png_path = os.path.join(mitos_dir, "result.png")

    fas_records = parse_fasta_multi(fas_path)
    faa_records = parse_fasta_multi(faa_path)
    bed_entries = parse_bed(bed_path)
    mitos_entries = parse_mitos_tab(mitos_path)

    # ── 02. Genes codificantes (nucleotídeos) ────────────────────────────
    cds_nt = []
    seen_nad3 = False
    for header, seq in fas_records:
        gene = get_gene_name(header)
        gene_base = CDS_ALIASES.get(gene, gene)
        if gene_base in CDS_GENES:
            # Para nad3 com frameshift, concatenar nad3_1 + nad3_0
            if gene_base == "nad3":
                if not seen_nad3:
                    seen_nad3 = True
                    # Buscar ambas as partes e concatenar
                    nad3_parts = [(h, s) for h, s in fas_records
                                  if get_gene_name(h) in ("nad3_0", "nad3_1", "nad3")]
                    if len(nad3_parts) >= 2:
                        # Ordenar por posição (nad3_1 geralmente é upstream)
                        nad3_parts.sort(key=lambda x: int(x[0].split(";")[1].strip().split("-")[0]))
                        combined_seq = nad3_parts[0][1] + nad3_parts[1][1]
                        display = GENE_DISPLAY_NAMES.get("nad3", "ND3")
                        pos0 = nad3_parts[0][0].split(";")[1].strip().split("-")
                        pos1 = nad3_parts[1][0].split(";")[1].strip().split("-")
                        new_header = f"{display} {pos0[0]}..{pos1[1]}(+) [frameshift join]"
                        cds_nt.append((new_header, combined_seq))
                    elif len(nad3_parts) == 1:
                        display = GENE_DISPLAY_NAMES.get("nad3", "ND3")
                        cds_nt.append((f"{display} {nad3_parts[0][0]}", nad3_parts[0][1]))
                continue
            display = GENE_DISPLAY_NAMES.get(gene_base, gene_base.upper())
            # Reconstruir header limpo: >GENE start..end(strand)
            parts = header.split(";")
            pos = parts[1].strip() if len(parts) > 1 else ""
            strand = parts[2].strip() if len(parts) > 2 else "+"
            cds_nt.append((f"{display} {pos}({strand})", seq))

    out_cds_nt = os.path.join(outdir, "02_coding_genes_nt.fasta")
    write_fasta(cds_nt, out_cds_nt)
    print(f"  [02] Genes codificantes (nt) → {out_cds_nt} ({len(cds_nt)} genes)")

    # ── 03. Genes codificantes (aminoácidos) ─────────────────────────────
    cds_aa = []
    seen_nad3_aa = False
    for header, seq in faa_records:
        gene = get_gene_name(header)
        gene_base = CDS_ALIASES.get(gene, gene)
        if gene_base in CDS_GENES:
            if gene_base == "nad3":
                if not seen_nad3_aa:
                    seen_nad3_aa = True
                    nad3_aa = [(h, s) for h, s in faa_records
                               if get_gene_name(h) in ("nad3_0", "nad3_1", "nad3")]
                    if len(nad3_aa) >= 2:
                        nad3_aa.sort(key=lambda x: int(x[0].split(";")[1].strip().split("-")[0]))
                        # Protein sequences: keep both reading frames separated by X
                        combined = nad3_aa[0][1].rstrip("*") + "X" + nad3_aa[1][1]
                        display = GENE_DISPLAY_NAMES.get("nad3", "ND3")
                        cds_aa.append((f"{display} [frameshift - two reading frames joined with X]", combined))
                    elif len(nad3_aa) == 1:
                        display = GENE_DISPLAY_NAMES.get("nad3", "ND3")
                        cds_aa.append((f"{display}", nad3_aa[0][1]))
                continue
            display = GENE_DISPLAY_NAMES.get(gene_base, gene_base.upper())
            parts = header.split(";")
            pos = parts[1].strip() if len(parts) > 1 else ""
            strand = parts[2].strip() if len(parts) > 2 else "+"
            cds_aa.append((f"{display} {pos}({strand})", seq))

    out_cds_aa = os.path.join(outdir, "03_coding_genes_aa.fasta")
    write_fasta(cds_aa, out_cds_aa)
    print(f"  [03] Genes codificantes (aa) → {out_cds_aa} ({len(cds_aa)} genes)")

    # ── 04. Genes ribossomais ────────────────────────────────────────────
    rrna = []
    for header, seq in fas_records:
        gene = get_gene_name(header)
        if gene in RRNA_GENES:
            display = GENE_DISPLAY_NAMES.get(gene, gene)
            parts = header.split(";")
            pos = parts[1].strip() if len(parts) > 1 else ""
            strand = parts[2].strip() if len(parts) > 2 else "+"
            rrna.append((f"{display} {pos}({strand})", seq))

    out_rrna = os.path.join(outdir, "04_ribosomal_genes.fasta")
    write_fasta(rrna, out_rrna)
    print(f"  [04] Genes ribossomais → {out_rrna} ({len(rrna)} genes)")

    # ── 05. Genes transportadores (tRNAs) ────────────────────────────────
    trna = []
    for header, seq in fas_records:
        gene_full = get_gene_name(header)
        gene_base = gene_full.split("(")[0]  # trnF(gaa) → trnF
        if gene_base in TRNA_GENES:
            aa = TRNA_AA_NAMES.get(gene_base, gene_base)
            parts = header.split(";")
            pos = parts[1].strip() if len(parts) > 1 else ""
            strand = parts[2].strip() if len(parts) > 2 else "+"
            trna.append((f"tRNA-{aa} {pos}({strand})", seq))

    out_trna = os.path.join(outdir, "05_transport_genes.fasta")
    write_fasta(trna, out_trna)
    print(f"  [05] Genes transportadores → {out_trna} ({len(trna)} tRNAs)")

    # ── 06. Tabela de posição física dos genes ───────────────────────────
    out_positions = os.path.join(outdir, "06_gene_positions.tsv")
    with open(out_positions, "w") as out:
        out.write("Region\tStart\tStop\tStrand\tLength\n")
        for entry in bed_entries:
            name = entry["name"]
            # Skip nad3_0/nad3_1 duplicates — show unified nad3
            if name in ("nad3_0", "nad3_1"):
                continue
            if name.startswith("OH"):
                display = f"D-loop ({name})"
            elif name in GENE_DISPLAY_NAMES:
                display = GENE_DISPLAY_NAMES[name]
            elif name.split("(")[0] in TRNA_AA_NAMES:
                base = name.split("(")[0]
                anticodon = name.split("(")[1].rstrip(")") if "(" in name else ""
                display = f"tRNA-{TRNA_AA_NAMES[base]}"
            else:
                display = name
            strand_display = "H" if entry["strand"] == "+" else "L"
            out.write(f"{display}\t{entry['start']}\t{entry['end']}\t{strand_display}\t{entry['length']}\n")

        # Add unified nad3 entry
        nad3_entries = [e for e in bed_entries if e["name"] in ("nad3_0", "nad3_1")]
        if len(nad3_entries) >= 2:
            nad3_entries.sort(key=lambda x: x["start"])
            start = nad3_entries[0]["start"]
            end = nad3_entries[-1]["end"]
            out.write(f"ND3 [frameshift]\t{start}\t{end}\tH\t{end - start + 1}\n")

    print(f"  [06] Posição dos genes → {out_positions}")

    # ── 07. Tabela de start/stop codons (CDS) ────────────────────────────
    out_codons = os.path.join(outdir, "07_start_stop_codons.tsv")
    with open(out_codons, "w") as out:
        out.write("Gene\tStart_Codon\tStop_Codon\tLength_bp\tLength_aa\n")
        for header, seq in fas_records:
            gene = get_gene_name(header)
            gene_base = CDS_ALIASES.get(gene, gene)
            if gene_base not in CDS_GENES:
                continue
            # Skip nad3 parts individually, handle unified
            if gene in ("nad3_0", "nad3_1"):
                continue
            display = GENE_DISPLAY_NAMES.get(gene_base, gene_base.upper())
            start_c, stop_c = extract_start_stop(seq, "+")
            aa_len = len(seq) // 3
            out.write(f"{display}\t{start_c}\t{stop_c}\t{len(seq)}\t{aa_len}\n")

        # nad3 unified
        nad3_nt_parts = [(h, s) for h, s in fas_records
                         if get_gene_name(h) in ("nad3_0", "nad3_1")]
        if len(nad3_nt_parts) >= 2:
            nad3_nt_parts.sort(key=lambda x: int(x[0].split(";")[1].strip().split("-")[0]))
            combined = nad3_nt_parts[0][1] + nad3_nt_parts[1][1]
            start_c = combined[:3].upper()
            stop_c = combined[-3:].upper()
            out.write(f"ND3 [frameshift]\t{start_c}\t{stop_c}\t{len(combined)}\t{len(combined)//3}\n")

    print(f"  [07] Start/Stop codons → {out_codons}")

    # ── 08. Tabela de anticódons dos tRNAs ───────────────────────────────
    out_anticodons = os.path.join(outdir, "08_trna_anticodons.tsv")
    with open(out_anticodons, "w") as out:
        out.write("tRNA\tLength\tAnticodon\tCodon\tAmino_Acid\n")
        for entry in mitos_entries:
            if entry["feat_type"] != "tRNA":
                continue
            name = entry["name"]
            if name not in TRNA_AA_NAMES:
                continue
            aa = TRNA_AA_NAMES[name]
            anticodon = entry["anticodon"]
            if anticodon:
                # Anticodon → codon (reverse complement)
                complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
                codon = "".join(complement.get(b, b) for b in anticodon.upper()[::-1])
            else:
                codon = "-"
                anticodon = "-"
            length = entry["end"] - entry["start"] + 1
            out.write(f"tRNA-{aa}\t{length}\t{anticodon.upper()}\t{codon}\t{aa}\n")

    print(f"  [08] Anticódons tRNA → {out_anticodons}")

    # ── 09. Mapa do genoma ───────────────────────────────────────────────
    if os.path.exists(png_path):
        out_map = os.path.join(outdir, "09_genome_map.png")
        shutil.copy2(png_path, out_map)
        print(f"  [09] Mapa do genoma → {out_map}")

    # ── 10. Anotação GFF3 ────────────────────────────────────────────────
    if os.path.exists(gff_path):
        out_gff = os.path.join(outdir, "10_annotation.gff")
        shutil.copy2(gff_path, out_gff)
        print(f"  [10] Anotação GFF3 → {out_gff}")

    # ── 11. Estruturas secundárias (SVGs) ────────────────────────────────
    svg_src = os.path.join(mitos_dir, "structure_svgs")
    if os.path.isdir(svg_src):
        svg_dst = os.path.join(outdir, "structure_svgs")
        if os.path.exists(svg_dst):
            shutil.rmtree(svg_dst)
        shutil.copytree(svg_src, svg_dst)
        svg_count = sum(1 for _, _, files in os.walk(svg_dst) for f in files if f.endswith(".svg"))
        print(f"  [11] Estruturas secundárias → {svg_dst} ({svg_count} SVGs)")

    # ── 12. Arquivos para GenBank ────────────────────────────────────────
    if args.genbank_dir and os.path.isdir(args.genbank_dir):
        gb_dst = os.path.join(outdir, "genbank_submission")
        if os.path.exists(gb_dst):
            shutil.rmtree(gb_dst)
        shutil.copytree(args.genbank_dir, gb_dst)
        print(f"  [12] GenBank submission → {gb_dst}")

    # ── 13. Ordem gênica ─────────────────────────────────────────────────
    geneorder_path = os.path.join(mitos_dir, "result.geneorder")
    if os.path.exists(geneorder_path):
        out_go = os.path.join(outdir, "13_gene_order.txt")
        shutil.copy2(geneorder_path, out_go)
        print(f"  [13] Ordem gênica → {out_go}")

    # ── 14. Plots de qualidade MITOS2 ────────────────────────────────────
    plots_src = os.path.join(mitos_dir, "plots")
    if os.path.isdir(plots_src):
        plots_dst = os.path.join(outdir, "quality_plots")
        if os.path.exists(plots_dst):
            shutil.rmtree(plots_dst)
        shutil.copytree(plots_src, plots_dst)
        print(f"  [14] Plots de qualidade → {plots_dst}")

    # ── Resumo final ─────────────────────────────────────────────────────
    # Read assembly to get size
    genome_size = 0
    with open(assembly) as fh:
        for line in fh:
            if not line.startswith(">"):
                genome_size += len(line.strip())

    print(f"\n{'='*60}")
    print(f"  RESUMO — {args.organism}")
    print(f"{'='*60}")
    print(f"  Genoma:       {genome_size:,} bp (circular)")
    print(f"  Genes CDS:    {len(cds_nt)}")
    print(f"  tRNAs:        {len(trna)}")
    print(f"  rRNAs:        {len(rrna)}")
    print(f"  Pasta saída:  {os.path.abspath(outdir)}")
    print(f"{'='*60}")


def main():
    parser = argparse.ArgumentParser(
        description="Compila resultados do pipeline em pasta organizada de entregáveis"
    )
    parser.add_argument("--assembly", required=True,
                        help="FASTA do genoma circularizado")
    parser.add_argument("--mitos-dir", required=True,
                        help="Diretório de saída do MITOS2 (contém result.*)")
    parser.add_argument("--genbank-dir", default=None,
                        help="Diretório com arquivos GenBank (.tbl, .fsa)")
    parser.add_argument("--organism", default="Unknown organism",
                        help='Nome do organismo (ex: "Anodorhynchus leari")')
    parser.add_argument("--outdir", default="summary",
                        help="Diretório de saída para entregáveis")
    args = parser.parse_args()

    if not os.path.isfile(args.assembly):
        print(f"ERRO: Arquivo de montagem não encontrado: {args.assembly}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(args.mitos_dir):
        print(f"ERRO: Diretório MITOS2 não encontrado: {args.mitos_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"\n  Compilando resultados para: {args.organism}")
    print(f"  {'─'*55}")
    compile_summary(args)


if __name__ == "__main__":
    main()
