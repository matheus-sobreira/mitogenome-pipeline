#!/usr/bin/env python3
"""
gff2genbank.py — Converte saída do MITOS2 (GFF3 + FASTA) para formato GenBank
(feature table .tbl + FASTA .fsa) pronto para submissão via table2asn.

Correções aplicadas sobre a anotação bruta do MITOS2 (todas registradas no
relatório de validação impresso ao final):

  1. Frameshift do ND3 (Mindell et al., 1998) — o MITOS2 representa o gene de
     DUAS formas conforme a montagem: (a) dois genes irmãos `nad3_0`/`nad3_1`,
     ou (b) um único gene com dois exons `nad3-a`/`nad3-b`. Ambas são
     reconhecidas e unificadas no join() com skip de 1 nt.
  2. Fragmentos espúrios — o MITOS2 pode reportar cópias curtas e de score
     baixíssimo do mesmo gene (ex.: um `nad5_1` de 60 nt com score 355 ao lado
     do `nad5_0` real com score 2e9). Cópias não-adjacentes são descartadas
     mantendo-se a de maior score.
  3. Extensão 5' — o início do CDS é estendido até o códon de iniciação em
     fase mais a montante que não invada a feature anterior.
  4. Extensão 3' — se os 3 nt seguintes em fase formam um códon de parada, ele
     é incorporado ao CDS.
  5. Códon de parada incompleto — CDS cujo comprimento não é múltiplo de 3
     termina em T ou TA completado por poliadenilação do mRNA; emite-se
     transl_except (pos:..,aa:TERM) + note, conforme exigido pelo NCBI.
  6. Região controle — anotada como D-loop cobrindo TODA a maior região
     não-codificante do genoma circular (e não apenas os fragmentos OH que o
     MITOS2 reporta), com join() quando cruza a origem das coordenadas.
  7. Anticódons — emitidos como note nos tRNAs, desambiguando os isoaceptores
     duplicados (tRNA-Leu UUR/CUN e tRNA-Ser UCN/AGY).

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
from collections import OrderedDict, defaultdict

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

# (símbolo, produto, anticódon) — anticódons canônicos do mtDNA de vertebrados
TRNA_PRODUCT = {
    "trnF":  ("trnF",  "tRNA-Phe", "GAA"),
    "trnV":  ("trnV",  "tRNA-Val", "TAC"),
    "trnL2": ("trnL2", "tRNA-Leu", "TAA"),  # UUR
    "trnL1": ("trnL1", "tRNA-Leu", "TAG"),  # CUN
    "trnI":  ("trnI",  "tRNA-Ile", "GAT"),
    "trnQ":  ("trnQ",  "tRNA-Gln", "TTG"),
    "trnM":  ("trnM",  "tRNA-Met", "CAT"),
    "trnW":  ("trnW",  "tRNA-Trp", "TCA"),
    "trnA":  ("trnA",  "tRNA-Ala", "TGC"),
    "trnN":  ("trnN",  "tRNA-Asn", "GTT"),
    "trnC":  ("trnC",  "tRNA-Cys", "GCA"),
    "trnY":  ("trnY",  "tRNA-Tyr", "GTA"),
    "trnS2": ("trnS2", "tRNA-Ser", "TGA"),  # UCN
    "trnS1": ("trnS1", "tRNA-Ser", "GCT"),  # AGY
    "trnD":  ("trnD",  "tRNA-Asp", "GTC"),
    "trnK":  ("trnK",  "tRNA-Lys", "TTT"),
    "trnG":  ("trnG",  "tRNA-Gly", "TCC"),
    "trnR":  ("trnR",  "tRNA-Arg", "TCG"),
    "trnH":  ("trnH",  "tRNA-His", "GTG"),
    "trnT":  ("trnT",  "tRNA-Thr", "TGT"),
    "trnP":  ("trnP",  "tRNA-Pro", "TGG"),
    "trnE":  ("trnE",  "tRNA-Glu", "TTC"),
}

RRNA_PRODUCT = {
    "rrnS": ("rrnS", "12S ribosomal RNA"),
    "rrnL": ("rrnL", "16S ribosomal RNA"),
}

# ── Código genético 2 (mitocôndria de vertebrados) ───────────────────────────

_BASES = "TCAG"
_AA2 = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
CODON_TABLE2 = {
    a + b + c: _AA2[i * 16 + j * 4 + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}

START_CODONS = {"ATG", "ATA", "ATT", "ATC", "GTG"}
STOP_CODONS = {"TAA", "TAG", "AGA", "AGG"}

_COMP = str.maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
                      "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")


def revcomp(s):
    """Reverse complement, preservando códigos IUPAC ambíguos."""
    return s.translate(_COMP)[::-1]


def translate2(nt):
    """Traduz com o código genético 2; códons incompletos/ambíguos viram 'X'."""
    return "".join(CODON_TABLE2.get(nt[i:i + 3].upper(), "X")
                   for i in range(0, len(nt) - len(nt) % 3, 3))


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
            try:
                score_val = float(score)
            except ValueError:
                score_val = 0.0
            features.append({
                "seqid":  seqid,
                "source": source,
                "type":   ftype,
                "start":  int(start),
                "end":    int(end),
                "score":  score_val,
                "strand": strand,
                "phase":  phase,
                "attrs":  attrs,
            })
    return features


def base_gene_name(name):
    """`nad3_1` → `nad3`; `nad3-a` → `nad3`; `OH_0` → `OH`."""
    return re.sub(r"[-_][0-9a-z]+$", "", name)


def resolve_cds_groups(gene_feats, exons_by_parent, transl_table, log):
    """Resolve cópias múltiplas de um mesmo gene codificante.

    Retorna uma lista de features de CDS, cada uma com `parts` = lista de
    intervalos (start, end) em coordenadas ascendentes.

    Duas situações são distinguidas:

    * **Frameshift programado** — as cópias são adjacentes ou se sobrepõem
      (gap < ADJACENT_MAX). É o caso do ND3 dos vertebrados não-mamíferos.
      As partes são unificadas num join() com skip de 1 nt (Mindell 1998).
    * **Fragmento espúrio** — as cópias estão distantes entre si. O MITOS2
      atribui score muito baixo ao fragmento falso; mantém-se apenas a cópia
      de maior score.
    """
    ADJACENT_MAX = 100  # nt de folga entre partes de um mesmo gene real

    def rank(f):
        """Critério de desempate entre cópias do mesmo gene.

        O campo `score` das linhas `gene` do GFF do MITOS2 é sempre '.'; o
        score de verdade fica nas linhas `exon` filhas. Sem herdar esse valor,
        todas as cópias empatam em zero e a escolha vira arbitrária — foi o que
        fez um fragmento espúrio de 60 nt vencer o ND5 real de 1809 nt.
        """
        exons = (exons_by_parent.get(f["attrs"].get("ID", "")) or
                 exons_by_parent.get(f["attrs"].get("Name", "")) or [])
        score = max([e["score"] for e in exons], default=f["score"])
        return (score, f["end"] - f["start"] + 1)

    groups = defaultdict(list)
    for f in gene_feats:
        groups[base_gene_name(f["attrs"].get("Name", ""))].append(f)

    resolved = []
    for gname, members in groups.items():
        if gname not in GENE_PRODUCT:
            continue

        # Um gene com múltiplos exons no GFF já traz o frameshift explícito.
        if len(members) == 1:
            f = members[0]
            exons = exons_by_parent.get(f["attrs"].get("ID", ""), [])
            if len(exons) < 2:
                exons = exons_by_parent.get(f["attrs"].get("Name", ""), [])
            if len(exons) >= 2:
                parts = sorted([(e["start"], e["end"]) for e in exons])
                f = dict(f)
                f["parts"] = apply_skip1(parts, gname, transl_table, log)
                f["frameshift"] = True
                resolved.append(f)
                continue
            f = dict(f)
            f["parts"] = [(f["start"], f["end"])]
            resolved.append(f)
            continue

        # Múltiplas cópias: adjacentes → frameshift; distantes → espúrio.
        members.sort(key=lambda x: x["start"])
        gap = members[1]["start"] - members[0]["end"]
        if gap <= ADJACENT_MAX:
            parts = [(m["start"], m["end"]) for m in members]
            merged = dict(members[0])
            merged["end"] = members[-1]["end"]
            merged["attrs"] = dict(merged["attrs"], Name=gname)
            merged["parts"] = apply_skip1(parts, gname, transl_table, log)
            merged["frameshift"] = True
            resolved.append(merged)
        else:
            best = max(members, key=rank)
            for m in members:
                if m is not best:
                    log.append(
                        f"fragmento espúrio descartado: {m['attrs'].get('Name')} "
                        f"{m['start']}..{m['end']} ({m['end'] - m['start'] + 1} nt, "
                        f"score {rank(m)[0]:.4g}) — mantido "
                        f"{best['attrs'].get('Name')} {best['start']}..{best['end']} "
                        f"({best['end'] - best['start'] + 1} nt, score {rank(best)[0]:.4g})")
            best = dict(best)
            best["attrs"] = dict(best["attrs"], Name=gname)
            best["parts"] = [(best["start"], best["end"])]
            resolved.append(best)

    resolved.sort(key=lambda x: x["start"])
    return resolved


def apply_skip1(parts, gname, transl_table, log):
    """Converte partes sobrepostas no join() canônico com skip de 1 nt."""
    if len(parts) < 2 or transl_table != 2:
        return parts
    (s1, e1), (s2, e2) = parts[0], parts[1]
    overlap = e1 - s2 + 1
    if overlap > 0:
        e1_adj = s2 - 2
        log.append(f"{gname.upper()}: sobreposição de {overlap} nt entre as partes; "
                   f"ajustado para o skip-1 canônico (Mindell 1998) → "
                   f"join({s1}..{e1_adj},{s2}..{e2})")
        return [(s1, e1_adj), (s2, e2)]
    return parts


def refine_cds(parts, strand, seq, blocked, gname, log):
    """Ajusta as bordas do CDS e detecta códon de parada incompleto.

    Devolve `(parts, partial_stop)`, onde `partial_stop` é a lista de posições
    (1-based, no genoma) do códon de parada incompleto, ou None.

    * Extensão 5': recua o início até o códon de iniciação em fase mais a
      montante que não invada `blocked` (posições ocupadas por outras features).
    * Extensão 3': incorpora o códon de parada seguinte, se houver.
    """
    L = len(seq)
    parts = [list(p) for p in parts]

    def nt_of(ps):
        out = ""
        for (a, b) in ps:
            out += seq[a - 1:b] if strand == "+" else revcomp(seq[a - 1:b])
        return out if strand == "+" else out

    def coding(ps):
        """Sequência codificante na orientação de leitura."""
        if strand == "+":
            return "".join(seq[a - 1:b] for (a, b) in ps)
        return "".join(revcomp(seq[a - 1:b]) for (a, b) in reversed(ps))

    # ── extensão 5' ──
    # Varre os códons em fase a montante e guarda o início válido mais distante.
    # A varredura para ao encontrar a feature vizinha (não se pode invadi-la) ou
    # um códon de parada em fase (não se pode ler através dele).
    origin = parts[0][0] if strand == "+" else parts[-1][1]
    furthest = None
    for step in range(1, 31):
        if strand == "+":
            new = origin - 3 * step
            if new < 1 or not blocked.isdisjoint(range(new, origin)):
                break
            codon = seq[new - 1:new + 2].upper()
        else:
            new = origin + 3 * step
            if new > L or not blocked.isdisjoint(range(origin + 1, new + 1)):
                break
            codon = revcomp(seq[new - 3:new]).upper()
        if codon in STOP_CODONS:
            break
        if codon in START_CODONS:
            furthest = (new, codon)

    if furthest:
        new, codon = furthest
        if strand == "+":
            parts[0][0] = new
        else:
            parts[-1][1] = new
        log.append(f"{gname.upper()}: início estendido para {new} "
                   f"(códon de iniciação {codon} em fase, sem invadir a feature anterior)")

    # ── extensão 3' ──
    #
    # O MITOS2 frequentemente trunca a extremidade 3' do CDS (no ND1 e no ND4
    # das araras, por 7 e 15 nt respectivamente, conferido contra OR209186.1).
    # Portanto a borda anotada não é confiável e o fim correto é procurado
    # varrendo em fase até o limite imposto pela feature seguinte:
    #
    #   · códon de parada completo — pode invadir a feature vizinha em até 2 nt
    #     (sobreposição CDS/tRNA é comum no mtDNA; o ND1 sobrepõe o trnI em 2 nt);
    #   · códon de parada incompleto (T ou TA completado por poliadenilação) —
    #     só é aceito encostando na vizinha, sem invadi-la.
    #
    # Havendo as duas possibilidades, a parada completa tem precedência.
    OVERLAP_TOL = 2

    original_end = parts[-1][1] if strand == "+" else parts[0][0]
    if strand == "+":
        after = [x for x in blocked if x > original_end]
        boundary = min(after) if after else L + 1
        hard_limit = min(boundary + OVERLAP_TOL, L)
    else:
        before = [x for x in blocked if x < original_end]
        boundary = max(before) if before else 0
        hard_limit = max(boundary - OVERLAP_TOL, 1)

    def with_end(e):
        p = [list(x) for x in parts]
        if strand == "+":
            p[-1][1] = e
        else:
            p[0][0] = e
        return [tuple(x) for x in p]

    complete_end = None
    partial_end = None
    steps = range(original_end, hard_limit + 1) if strand == "+" \
        else range(original_end, hard_limit - 1, -1)
    for e in steps:
        cds = coding(with_end(e))
        rem = len(cds) % 3
        if rem == 0:
            if cds[-3:].upper() in STOP_CODONS and complete_end is None:
                complete_end = e
                break
        elif cds[-rem:].upper() in ("T", "TA"):
            touching = (e < boundary) if strand == "+" else (e > boundary)
            if touching:
                partial_end = e

    partial_stop = None
    if complete_end is not None:
        if complete_end != original_end:
            parts = [list(p) for p in with_end(complete_end)]
            log.append(f"{gname.upper()}: extremidade 3' estendida de {original_end} "
                       f"para {complete_end} até o códon de parada em fase "
                       f"(o MITOS2 havia truncado o CDS)")
    elif partial_end is not None:
        parts = [list(p) for p in with_end(partial_end)]
        cds = coding([tuple(p) for p in parts])
        rem = len(cds) % 3
        if strand == "+":
            partial_stop = (partial_end - rem + 1, partial_end)
        else:
            partial_stop = (partial_end, partial_end + rem - 1)
        extended = "" if partial_end == original_end else \
            f" (extremidade 3' estendida de {original_end} para {partial_end})"
        log.append(f"{gname.upper()}: comprimento {len(cds)} nt (resto {rem}) — "
                   f"códon de parada incompleto completado por poliadenilação; "
                   f"emitido transl_except em {partial_stop[0]}..{partial_stop[1]}{extended}")

    return [tuple(p) for p in parts], partial_stop


def find_control_region(occupied, genome_length):
    """Maior região não-codificante do genoma circular → região controle.

    Devolve lista de intervalos (1 se contígua, 2 se cruza a origem).
    """
    if not occupied:
        return []
    ivs = sorted(occupied)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    gaps = []
    for i in range(len(merged) - 1):
        gaps.append(([(merged[i][1] + 1, merged[i + 1][0] - 1)],
                     merged[i + 1][0] - merged[i][1] - 1))
    # gap que cruza a origem
    tail, head = merged[-1][1], merged[0][0]
    wrap = []
    if tail < genome_length:
        wrap.append((tail + 1, genome_length))
    if head > 1:
        wrap.append((1, head - 1))
    if wrap:
        gaps.append((wrap, sum(b - a + 1 for a, b in wrap)))

    if not gaps:
        return []
    return max(gaps, key=lambda g: g[1])[0]


def write_tbl(cds_feats, rna_feats, dloop_parts, seqid, transl_table, seq, outpath):
    """Write GenBank feature table (.tbl) in NCBI 5-column format."""
    rows = []
    for f in cds_feats:
        rows.append((f["parts"][0][0], "cds", f))
    for f in rna_feats:
        rows.append((f["start"], "rna", f))
    rows.sort(key=lambda r: r[0])

    with open(outpath, "w") as out:
        out.write(f">Feature {seqid}\n")

        if dloop_parts:
            (a, b) = dloop_parts[0]
            out.write(f"{a}\t{b}\tD-loop\n")
            for (a, b) in dloop_parts[1:]:
                out.write(f"{a}\t{b}\n")

        for _, kind, f in rows:
            name = base_gene_name(f["attrs"].get("Name", ""))
            comp = f["strand"] == "-"

            if kind == "cds":
                gene_sym, product = GENE_PRODUCT[name]
                parts = f["parts"]
                gs, ge = parts[0][0], parts[-1][1]

                # feature `gene` cobre a extensão total
                out.write(f"{ge}\t{gs}\tgene\n" if comp else f"{gs}\t{ge}\tgene\n")
                out.write(f"\t\t\tgene\t{gene_sym}\n")

                # feature `CDS` — join() quando há frameshift
                ordered = list(reversed(parts)) if comp else parts
                first = True
                for (a, b) in ordered:
                    lo, hi = (b, a) if comp else (a, b)
                    out.write(f"{lo}\t{hi}\tCDS\n" if first else f"{lo}\t{hi}\n")
                    first = False
                out.write(f"\t\t\tgene\t{gene_sym}\n")
                out.write(f"\t\t\tproduct\t{product}\n")
                out.write(f"\t\t\ttransl_table\t{transl_table}\n")

                if f.get("partial_stop"):
                    a, b = f["partial_stop"]
                    pos = f"{a}" if a == b else f"{a}..{b}"
                    if comp:
                        pos = f"complement({pos})"
                    out.write(f"\t\t\ttransl_except\t(pos:{pos},aa:TERM)\n")
                    out.write("\t\t\tnote\tTAA stop codon is completed by the "
                              "addition of 3' A residues to the mRNA\n")

                if f.get("frameshift"):
                    out.write("\t\t\texception\tribosomal slippage\n")
                    out.write("\t\t\tnote\tprogrammed frameshift; frameshift mechanism "
                              "unknown (Mindell et al., 1998, Mol. Biol. Evol., "
                              "15:1568-1571)\n")

            elif name in TRNA_PRODUCT:
                gene_sym, product, anticodon = TRNA_PRODUCT[name]
                s, e = f["start"], f["end"]
                out.write(f"{e}\t{s}\tgene\n" if comp else f"{s}\t{e}\tgene\n")
                out.write(f"\t\t\tgene\t{gene_sym}\n")
                out.write(f"{e}\t{s}\ttRNA\n" if comp else f"{s}\t{e}\ttRNA\n")
                out.write(f"\t\t\tproduct\t{product}\n")
                out.write(f"\t\t\tnote\tanticodon:{anticodon}\n")

            elif name in RRNA_PRODUCT:
                gene_sym, product = RRNA_PRODUCT[name]
                s, e = f["start"], f["end"]
                out.write(f"{e}\t{s}\tgene\n" if comp else f"{s}\t{e}\tgene\n")
                out.write(f"\t\t\tgene\t{gene_sym}\n")
                out.write(f"{e}\t{s}\trRNA\n" if comp else f"{s}\t{e}\trRNA\n")
                out.write(f"\t\t\tproduct\t{product}\n")


def write_fsa(fasta_path, seqid, organism, topology, transl_table, outpath):
    """Write GenBank-formatted FASTA (.fsa) with proper header.

    Os modificadores de origem precisam usar os nomes que o table2asn
    reconhece: 'moltype' (não 'molecule_type') e 'mgcode' (código genético
    MITOCONDRIAL, não 'gcode', que é o nuclear). Com os nomes errados o
    table2asn trata o colchete como modificador desconhecido e deixa TODO o
    header sem parsear — sem organismo, sem location=mitochondrion, código
    genético 0/1 — o que gera ~50 erros espúrios de 'internal stop' (TGA lido
    como stop em vez de Trp). Verificado na submissão GenBank das araras
    (17/08/2026). mgcode usa o transl_table do organismo (2=vertebrado,
    5=invertebrado), não um valor fixo.
    """
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
            f"[moltype=genomic DNA] "
            f"[mgcode={transl_table}]"
        )
        out.write(header + "\n")
        # Write 70-char lines
        for i in range(0, len(sequence), 70):
            out.write(sequence[i:i+70] + "\n")


def validate_cds(cds_feats, seq):
    """Traduz cada CDS e devolve uma lista de problemas remanescentes."""
    problems = []
    for f in cds_feats:
        gene_sym, _ = GENE_PRODUCT[base_gene_name(f["attrs"].get("Name", ""))]
        parts = f["parts"]
        if f["strand"] == "+":
            nt = "".join(seq[a - 1:b] for (a, b) in parts)
        else:
            nt = "".join(revcomp(seq[a - 1:b]) for (a, b) in reversed(parts))
        prot = translate2(nt)
        body = prot[:-1] if prot.endswith("*") else prot
        if "*" in body:
            problems.append(f"{gene_sym}: {body.count('*')} códon(s) de parada interno(s)")
        # GTG/ATT/ATC iniciam tradução no mtDNA de vertebrados sem codificar Met
        # (o COX1 de aves tipicamente inicia em GTG), logo a checagem é sobre o
        # códon e não sobre o aminoácido traduzido.
        if nt[:3].upper() not in START_CODONS:
            problems.append(f"{gene_sym}: códon de iniciação inválido ({nt[:3]})")
        if not prot.endswith("*") and not f.get("partial_stop"):
            problems.append(f"{gene_sym}: sem códon de parada terminal")
        if "X" in prot:
            problems.append(f"{gene_sym}: base ambígua dentro do CDS "
                            f"(posição(ões) {[i + 1 for i, c in enumerate(nt.upper()) if c not in 'ACGT']})")
    return problems


def build_annotation(gff_path, fasta_path, transl_table=2):
    """Anotação corrigida do mitogenoma a partir da saída bruta do MITOS2.

    É o **ponto único de verdade** da anotação: tanto o `.tbl`/`.fsa` de
    submissão (este módulo) quanto o `.gbk` e o mapa circular
    (`generate_genbank.py`) devem consumir esta função, sob pena de os
    entregáveis divergirem entre si.

    Devolve um dicionário com:
      `seq`, `genome_length`, `cds_feats`, `rna_feats`, `dloop_parts`, `log`.

    Cada item de `cds_feats` traz `parts` (lista de intervalos ascendentes; mais
    de um quando há frameshift), `partial_stop` e `frameshift`.
    """
    seq = "".join(l.strip() for l in open(fasta_path) if not l.startswith(">"))

    features = parse_gff(gff_path)

    genome_length = 0
    for f in features:
        if f["type"] == "region":
            genome_length = f["end"]
            break
    if genome_length == 0:
        genome_length = len(seq)

    log = []

    # Exons indexados pelo transcrito pai (forma alternativa do frameshift)
    exons_by_parent = defaultdict(list)
    for f in features:
        if f["type"] == "exon":
            parent = f["attrs"].get("Parent", "")
            exons_by_parent[parent.replace("transcript_", "gene_")].append(f)
            exons_by_parent[parent.replace("transcript_", "")].append(f)

    gene_feats = [f for f in features if f["type"] in ("gene", "ncRNA_gene")]

    coding = [f for f in gene_feats
              if base_gene_name(f["attrs"].get("Name", "")) in GENE_PRODUCT]
    rna_feats = [f for f in gene_feats
                 if base_gene_name(f["attrs"].get("Name", "")) in TRNA_PRODUCT
                 or base_gene_name(f["attrs"].get("Name", "")) in RRNA_PRODUCT]

    cds_feats = resolve_cds_groups(coding, exons_by_parent, transl_table, log)

    # Posições ocupadas por features vizinhas (impedem a extensão do CDS)
    blocked = set()
    for f in rna_feats:
        blocked.update(range(f["start"], f["end"] + 1))
    for f in cds_feats:
        for (a, b) in f["parts"]:
            blocked.update(range(a, b + 1))

    for f in cds_feats:
        own = set()
        for (a, b) in f["parts"]:
            own.update(range(a, b + 1))
        gname = base_gene_name(f["attrs"].get("Name", ""))
        f["parts"], f["partial_stop"] = refine_cds(
            f["parts"], f["strand"], seq, blocked - own, gname, log)

    # Região controle = maior região não-codificante
    occupied = [(f["start"], f["end"]) for f in rna_feats]
    for f in cds_feats:
        occupied.extend(f["parts"])
    dloop_parts = find_control_region(occupied, genome_length)
    if dloop_parts:
        total = sum(b - a + 1 for a, b in dloop_parts)
        span = "+".join(f"{a}..{b}" for a, b in dloop_parts)
        log.append(f"Região controle (D-loop) anotada como {span} ({total} bp)")
    for f in features:
        if f["type"] == "origin_of_replication":
            log.append(f"fragmento OH do MITOS2 em {f['start']}..{f['end']} "
                       f"absorvido pela região controle")

    return {
        "seq": seq,
        "genome_length": genome_length,
        "cds_feats": cds_feats,
        "rna_feats": rna_feats,
        "dloop_parts": dloop_parts,
        "log": log,
    }


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

    ann = build_annotation(args.gff, args.fasta, args.transl_table)
    seq = ann["seq"]
    genome_length = ann["genome_length"]
    cds_feats = ann["cds_feats"]
    rna_feats = ann["rna_feats"]
    dloop_parts = ann["dloop_parts"]
    log = ann["log"]

    # Write outputs
    tbl_path = os.path.join(args.outdir, f"{seqid}.tbl")
    fsa_path = os.path.join(args.outdir, f"{seqid}.fsa")

    write_tbl(cds_feats, rna_feats, dloop_parts, seqid, args.transl_table, seq, tbl_path)
    write_fsa(args.fasta, seqid, args.organism, args.topology, args.transl_table, fsa_path)

    # Summary
    cds_count = len(cds_feats)
    trna_count = sum(1 for f in rna_feats
                     if base_gene_name(f["attrs"].get("Name", "")) in TRNA_PRODUCT)
    rrna_count = sum(1 for f in rna_feats
                     if base_gene_name(f["attrs"].get("Name", "")) in RRNA_PRODUCT)

    print(f"\n  Organismo:    {args.organism}")
    print(f"  Comprimento: {genome_length} bp ({args.topology})")
    print(f"  Genes CDS:   {cds_count}")
    print(f"  tRNAs:        {trna_count}")
    print(f"  rRNAs:        {rrna_count}")
    print(f"  Código gen.:  {args.transl_table} (mitocôndria de vertebrados)")

    if log:
        print(f"\n  Correções aplicadas sobre a anotação bruta do MITOS2:")
        for item in log:
            print(f"    · {item}")

    problems = validate_cds(cds_feats, seq)
    if problems:
        print(f"\n  ATENÇÃO — problemas remanescentes (revisar antes de submeter):")
        for p in problems:
            print(f"    ! {p}")
    else:
        print(f"\n  Validação: todos os {cds_count} CDS traduzem sem stop interno, "
              f"iniciam em Met e terminam em códon de parada.")

    print(f"\n  Arquivos gerados:")
    print(f"    {tbl_path}")
    print(f"    {fsa_path}")
    print(f"\n  Para gerar .sqn (submissão):")
    print(f"    table2asn -t template.sbt -i {fsa_path} -f {tbl_path} "
          f"-o {seqid}.sqn -V vb -Z")


if __name__ == "__main__":
    main()
