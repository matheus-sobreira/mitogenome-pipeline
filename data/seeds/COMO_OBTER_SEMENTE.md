# Como Obter a Semente para o NOVOPlasty

O NOVOPlasty precisa de uma sequência **semente (seed)** — tipicamente o gene **cox1** — para iniciar a montagem.

## Espécie: *Caenorhabditis elegans*

- **Referência mtDNA completo:** NC_001328.1 (13.794 bp)
- **Gene cox1:** posição ~4.758–6.287 do mtDNA

---

## Opção 1 — Download via NCBI (recomendado)

Acesse a página do gene cox1 de *C. elegans* no NCBI:

```
https://www.ncbi.nlm.nih.gov/gene/3565218
```

Ou baixe diretamente a region cox1 do mtDNA completo:

1. Acesse: https://www.ncbi.nlm.nih.gov/nuccore/NC_001328.1
2. Clique em **"Send to"** → **"File"** → FORMAT: **FASTA**
3. Marque **"Selected Features"** e selecione o gene `cox1`
4. Salve como `c_elegans_cox1.fasta` nesta pasta

---

## Opção 2 — Script Python (download automático via Entrez)

Salve e execute o script abaixo (`baixar_semente.py`) na pasta `data/seeds/`:

```python
from Bio import Entrez, SeqIO

Entrez.email = "seu_email@exemplo.com"   # Obrigatório pelo NCBI

# Baixa o mtDNA completo de C. elegans (NC_001328.1)
handle = Entrez.efetch(
    db="nucleotide", id="NC_001328.1",
    rettype="gb", retmode="text"
)
record = SeqIO.read(handle, "genbank")
handle.close()

# Extrai apenas a feature do gene cox1
for feature in record.features:
    if feature.type == "CDS" and "cox1" in feature.qualifiers.get("gene", []):
        seq = feature.extract(record.seq)
        with open("c_elegans_cox1.fasta", "w") as f:
            f.write(f">cox1_Caenorhabditis_elegans_NC_001328.1\n{seq}\n")
        print(f"Semente salva: {len(seq)} bp")
        break
```

Instale a dependência (se necessário):
```bash
pip install biopython
```

---

## Formato esperado do arquivo `c_elegans_cox1.fasta`

```
>cox1_Caenorhabditis_elegans_NC_001328.1
ATGTTTTTGATTTTTTCCAAGGCCATCCCGTTGTAATCGGAGGTTTCGGAAATTGACTTC
GTCCGTTAATACTTGGAGCCCCAGATATAGCTTTCCCTCGTATAAATAATATAAGCTTCT
...
```

---

## Nota para outras espécies futuramente

| Espécie | Acesso mtDNA | Cox1 aprox. (bp) |
|---|---|---|
| *C. elegans* | NC_001328.1 | ~1.530 |
| *Anodorhynchus leari* | A obter/publicar | ~1.560 |
| Espécie comparação (aves) | A definir | ~1.560 |
