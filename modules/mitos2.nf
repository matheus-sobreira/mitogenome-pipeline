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

    publishDir "${params.outdir}/annotation/mitos2", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path  db_dir

    output:
    tuple val(sample_id), path("${sample_id}/"),     emit: annotation
    path "${sample_id}/*.gff",                        emit: gff,      optional: true
    path "${sample_id}/*.faa",                        emit: proteins, optional: true
    path "${sample_id}/*.bed",                        emit: bed,      optional: true

    script:
    def db_name   = db_dir.name
    def db_parent = db_dir.parent
    """
    runmitos.py \\
        -i ${assembly} \\
        -c ${params.genetic_code} \\
        -o ${sample_id} \\
        -r ${db_name} \\
        -R ${db_parent}
    """
}
