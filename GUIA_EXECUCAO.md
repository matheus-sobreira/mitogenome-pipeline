# Guia de Execução — mitogenome-pipeline

## Pré-requisitos

- **Docker** instalado e rodando (`docker info`)
- **Nextflow** ≥ 23.x (`nextflow -version`)
- **Java** 17+ (`java -version`)
- ~50 GB de disco livre (download SRA + trabalho temporário)

---

## 1. Construir as imagens Docker (apenas uma vez)

```bash
cd mitogenome-pipeline

# Linux/macOS
for dir in sra-tools fastqc trim-galore novoplasty mitos2; do
    docker build -t "mitogenome-pipeline/${dir}:1.0" "docker/${dir}"
done

# Windows (PowerShell)
.\build_images.ps1
```

Verificar se foram criadas:

```bash
docker images | grep mitogenome-pipeline
```

Saída esperada:

```
mitogenome-pipeline/sra-tools      1.0
mitogenome-pipeline/fastqc         1.0
mitogenome-pipeline/trim-galore    1.0
mitogenome-pipeline/novoplasty     1.0
mitogenome-pipeline/mitos2         1.0
```

---

## 2. Baixar o banco de dados MITOS2 (apenas uma vez)

```bash
mkdir -p data/databases
wget https://zenodo.org/record/4284483/files/refseq89m.tar.bz2
tar -xjf refseq89m.tar.bz2 -C data/databases/
rm refseq89m.tar.bz2
```

Verificar:

```bash
ls data/databases/refseq89m/
# Esperado: featureProt/  ncRNA/  prot_db.fas  ...
```

---

## 3. Executar o pipeline

### Perfil A. leari (arara-azul-de-lear)

```bash
nextflow run main.nf -profile a_leari
```

### Perfil de teste (D. bifasciatus)

```bash
nextflow run main.nf -profile test
```

### Execução resiliente (recomendado para sessões longas)

O script `run_pipeline.sh` executa o pipeline em background com `nohup`.
O terminal pode ser fechado sem interromper a execução:

```bash
# Executar com qualquer perfil
./run_pipeline.sh -profile a_leari

# Com parâmetros extras
./run_pipeline.sh -profile a_leari --sra_max_reads 20000000

# Retomar run anterior
./run_pipeline.sh -profile a_leari -resume -w work/a_leari/2026-04-12_14h30
```

Após iniciar, o script informa:
- **Log**: `logs/a_leari_2026-04-12_14h30.log`
- **PID**: para cancelar com `kill <PID>`

Acompanhar em tempo real:

```bash
tail -f logs/a_leari_2026-04-12_14h30.log
```

Ver progresso resumido:

```bash
grep -E 'process|Completed|ERROR' logs/a_leari_2026-04-12_14h30.log
```

### Modo direto (terminal aberto)

Se preferir acompanhar diretamente no terminal (não resiliente a desconexão):

```bash
nextflow run main.nf -profile a_leari
```

---

## 4. Fluxo do pipeline

```
SRA_PILOT (se sra_max_reads não definido)
    │ analisa 500K reads → recomenda max_reads
    ▼
SRA_DOWNLOAD → FASTQC → TRIM_GALORE → NOVOPLASTY → MITOS2 → COMPILE_SUMMARY
```

- **SRA_PILOT**: baixa amostra de 500K reads, analisa Q30%, adaptadores e fração mitocondrial, calcula o número ideal de reads. Roda automaticamente quando `sra_max_reads` não está definido no config.
- **SRA_DOWNLOAD**: baixa e trunca as reads do SRA
- **FASTQC**: controle de qualidade das reads brutas
- **TRIM_GALORE**: remoção de adaptadores e reads de baixa qualidade
- **NOVOPLASTY**: montagem do mitogenoma por seed-and-extend
- **MITOS2**: anotação funcional (13 CDS, 22 tRNA, 2 rRNA)
- **COMPILE_SUMMARY**: compila todos os entregáveis em `<outdir>/summary/`

---

## 5. Forçar número de reads (pular pilot)

Se você já sabe quantas reads usar:

```bash
nextflow run main.nf -profile a_leari --sra_max_reads 20000000
```

---

## 6. Retomar execução interrompida

O diretório `work/` é organizado por espécie e timestamp:

```
work/
├── a_leari/
│   ├── 2026-04-12_14h30/    ← run que falhou
│   └── 2026-04-12_20h15/    ← run que funcionou
└── d_bifasciatus/
    └── 2026-04-10_09h00/
```

Para retomar um run específico, use `-resume` com `-w` apontando para o diretório:

```bash
# Retomar o último run (Nextflow detecta automaticamente)
nextflow run main.nf -profile a_leari -resume -w work/a_leari/2026-04-12_20h15

# Sem -w, uma nova pasta com timestamp é criada (não encontra o cache anterior)
```

> **Dica**: use `ls work/a_leari/` para ver os runs anteriores e escolher qual retomar.

---

## 7. Saídas

Após a execução, os resultados ficam em `<outdir>/` (ex: `results/a_leari/`):

```
results/a_leari/
├── assembly/                          # Montagem NOVOPlasty
│   ├── Circularized_assembly_*.fasta
│   ├── config.txt
│   └── log_*.txt
├── annotation/
│   └── mitos2/                        # Saídas brutas MITOS2
├── qc/
│   ├── fastqc_raw/                    # Relatórios FastQC
│   └── trim_galore/                   # Relatórios Trim Galore
└── summary/                           # ★ Entregáveis compilados
    ├── 01_genome_assembly.fasta       # Genoma circularizado
    ├── 02_coding_genes_nt.fasta       # 13 CDS (nucleotídeos)
    ├── 03_coding_genes_aa.fasta       # 13 CDS (aminoácidos)
    ├── 04_ribosomal_genes.fasta       # 12S + 16S rRNA
    ├── 05_transport_genes.fasta       # 22 tRNAs
    ├── 06_gene_positions.tsv          # Posição dos genes
    ├── 07_start_stop_codons.tsv       # Start/Stop codons
    ├── 08_trna_anticodons.tsv         # Anticódons tRNA
    ├── 09_circular_map.svg            # Mapa circular
    ├── 09_circular_map.pdf            # Mapa circular (PDF)
    ├── 09b_genome_map_linear.png      # Mapa linear MITOS2
    ├── 10_annotation.gff              # Anotação GFF3
    ├── 11_structure_svgs/             # Estruturas secundárias tRNA/rRNA
    ├── 12_quality_plots/              # Plots de qualidade MITOS2
    ├── 13_gene_order.txt              # Ordem gênica
    └── genbank_submission/            # Arquivos para submissão GenBank
        ├── <especie>.gbk              # GenBank Flat File
        ├── <especie>.tbl              # Feature table
        └── <especie>.fsa              # FASTA formatado
```

---

## 8. Limpeza

```bash
# Remover cache de trabalho de uma espécie inteira
rm -rf work/a_leari/

# Remover apenas um run específico
rm -rf work/a_leari/2026-04-12_14h30/

# Remover resultados de uma espécie
sudo rm -rf results/a_leari/

# Remover logs do Nextflow
rm -rf .nextflow/ .nextflow.log*

# Limpeza completa (tudo exceto código e banco de dados)
rm -rf work/ results/ .nextflow/ .nextflow.log*
```

> **Nota**: `sudo` pode ser necessário pois Docker cria arquivos como root.

---

## 9. Executar para uma nova espécie

1. **Encontrar o SRA**: buscar no [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) por WGS paired-end da espécie
2. **Obter semente COX1**: baixar do GenBank o gene cox1 de uma espécie próxima
3. **Criar config** em `conf/nova_especie.config`:

```groovy
params {
    sra_accession = 'SRRXXXXXXX'
    seed          = "${projectDir}/data/seeds/especie_proxima_cox1.fasta"
    outdir        = 'results/nova_especie'
    genome_range  = '15000-18500'     // Ajustar ao tamanho esperado do mtDNA
    novoplasty_kmers = '39,33'
    mitos2_db    = "${projectDir}/data/databases/refseq89m"
    genetic_code = 2                  // 2=vertebrado, 5=invertebrado
    organism     = 'Genus species'
}
```

4. **Adicionar perfil** em `nextflow.config`:

```groovy
profiles {
    nova_especie {
        includeConfig 'conf/nova_especie.config'
    }
}
```

5. **Executar**:

```bash
nextflow run main.nf -profile nova_especie
```
