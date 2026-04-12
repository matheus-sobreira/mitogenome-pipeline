/*
 * Módulo: MITOS2
 * Anotação funcional do genoma mitocondrial montado pelo NOVOPlasty.
 *
 * Detecta e anota:
 *   - 13 genes proteicos (nad1-6, nad4L, cox1-3, atp6, atp8, cytb)
 *   - 2 rRNAs (rrnL, rrnS)
 *   - 22 tRNAs
 *   - Região de controle (D-loop / CR)
 *
 * Entradas : tuple (sample_id, assembly.fasta)  ← montagem circularizada (NOVOPlasty)
 *            path(db_dir)                        ← banco de dados RefSeq (download externo)
 *
 * Saídas   : diretório completo da anotação
 *            arquivo GFF com genes anotados
 *            proteínas preditas (FAA)
 *
 * Parâmetros relevantes (nextflow.config):
 *   params.genetic_code — código genético mitocondrial (2 = vertebrado, 5 = invertebrado)
 *   params.mitos2_db    — caminho completo para o banco mitos2-refseq89m.db
 *
 * Referência: https://gitlab.com/Bernt/MITOS
 * Banco de dados: https://zenodo.org/record/4284483
 */

process MITOS2 {

    tag "${sample_id}"

    publishDir "${params.outdir}/annotation/mitos2", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(assembly)
    path  db_dir

    output:
    tuple val(sample_id), path("${sample_id}/"),     emit: annotation
    path "${sample_id}/*.gff",                        emit: gff,      optional: true
    path "${sample_id}/*.faa",                        emit: proteins, optional: true
    path "${sample_id}/*.bed",                        emit: bed,      optional: true
    path "${sample_id}/structure_svgs/**/*.svg",      emit: svgs,     optional: true

    script:
    def db_name = db_dir.name
    """
    mkdir -p ${sample_id}
    runmitos.py \\
        -i ${assembly} \\
        -c ${params.genetic_code} \\
        -o ${sample_id} \\
        -r ${db_name} \\
        -R \${PWD}

    # ── Gerar SVGs de estrutura secundária (tRNA + rRNA) via RNAplot ──
    WORKDIR=\$(pwd)
    MITFI="\${WORKDIR}/${sample_id}/mitfi-global"
    SVGDIR="\${WORKDIR}/${sample_id}/structure_svgs"
    mkdir -p "\${SVGDIR}/tRNA" "\${SVGDIR}/rRNA"
    TMPDIR=\$(mktemp -d) && cd "\${TMPDIR}"

    for NC_TYPE in tRNA rRNA; do
        if [[ "\${NC_TYPE}" == "tRNA" ]]; then
            NC_FILE="\${MITFI}/sequence.fas-0_tRNAout.nc"
            DEST="\${SVGDIR}/tRNA"
        else
            NC_FILE="\${MITFI}/sequence.fas-0_rRNAout.nc"
            DEST="\${SVGDIR}/rRNA"
        fi
        [[ ! -f "\${NC_FILE}" ]] && continue

        name="" seq="" struct=""
        while IFS= read -r line; do
            if [[ "\${line}" == ">"* ]]; then
                if [[ -n "\${name}" && -n "\${seq}" && -n "\${struct}" ]]; then
                    echo -e ">\${name}\\n\${seq}\\n\${struct}" | RNAplot -f svg 2>/dev/null
                    svgout=\$(ls -t *_ss.svg 2>/dev/null | head -1)
                    [[ -n "\${svgout}" ]] && mv "\${svgout}" "\${DEST}/\${name}.svg"
                fi
                is_primary=\$(echo "\${line}" | awk -F'|' '{print \$(NF-2)"|"\$(NF-1)}')
                if [[ "\${is_primary}" != "1|1" ]]; then name=""; continue; fi
                name=\$(echo "\${line}" | awk -F'|' '{print \$(NF-3)}')
                seq="" struct=""
            elif [[ -n "\${name}" ]]; then
                if [[ -z "\${seq}" ]]; then seq="\${line}"; else struct="\${line}"; fi
            fi
        done < "\${NC_FILE}"
        if [[ -n "\${name}" && -n "\${seq}" && -n "\${struct}" ]]; then
            echo -e ">\${name}\\n\${seq}\\n\${struct}" | RNAplot -f svg 2>/dev/null
            svgout=\$(ls -t *_ss.svg 2>/dev/null | head -1)
            [[ -n "\${svgout}" ]] && mv "\${svgout}" "\${DEST}/\${name}.svg"
        fi
    done

    cd - > /dev/null
    rm -rf "\${TMPDIR}"
    """
}
