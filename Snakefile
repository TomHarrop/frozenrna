#!/usr/bin/env python3


#############
# FUNCTIONS #
#############

def get_reads(wildcards):
    input_keys = ['r1', 'r2']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}


def get_bc(wildcards):
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return my_pep['bc']


def pick_trinity_input(wildcards):
    if wildcards.run == 'merged':
        return {
            'r1': 'output/000_tmp/{sample}.r1_joined_with_merged.fastq',
            'r2': 'output/020_merged/{sample}_R2.fastq'}
    elif wildcards.run == 'raw':
        return {
            'r1': 'output/010_reads/{sample}_R1.fastq',
            'r2': 'output/010_reads/{sample}_R2.fastq'}
    else:
        raise ValueError(f'wtf run {wildcards.run}')


###########
# GLOBALS #
###########

pepfile: 'config/config.yaml'
all_samples = pep.sample_table['sample_name']

# this is the bbduk adaptors file with the illumina stranded mrna added
adaptors = 'data/adaptors.fasta'

# containers
bbduk = 'shub://TomHarrop/seq-utils:bbmap_38.76'
bioconductor = 'shub://TomHarrop/r-containers:bioconductor_3.10'
biopython = 'shub://TomHarrop/singularity-containers:biopython_1.73'
busco = 'docker://ezlabgva/busco:v4.0.4_cv1'
pandas_container = 'shub://TomHarrop/py-containers:pandas_0.25.3'
r = 'shub://TomHarrop/r-containers:r_3.6.2'
trinity = 'shub://TomHarrop/assemblers:trinity_2.11.0'
trinotate = 'shub://TomHarrop/trinotate_pipeline:v0.0.12'

########
# MAIN #
########

# read the csv of raw file locations
# raw_read_df = pandas.read_csv(raw_read_csv,
#                               index_col=['indiv', 'lane', 'read'])

# all_indivs = sorted(set(x[0] for x in raw_read_df.index.values))
# all_lanes = sorted(set(x[1] for x in raw_read_df.index.values))
# indiv_with_r2 = sorted(set(x[0] for x in
#                        raw_read_df.loc[(
#                         slice(None),
#                         slice(None),
#                         'R2'), :].index.values))
# indiv_r1_only = [x for x in all_indivs if x not in indiv_with_r2]
# ordered_indivs = indiv_with_r2 + indiv_r1_only

#########
# RULES #
#########

wildcard_constraints:
    sample = '|'.join(all_samples),


# rule target:
#     input:
#         expand('output/040_trinity-abundance/{run}/ExN50.pdf',
#                run=['merged', 'raw']),
#         expand('output/060_de/{run}/pca.pdf',
#                run=['merged', 'raw']),
#         expand('output/060_de/{run}/de_wald_all_pfam.csv',
#                run=['merged', 'raw']),
#         expand('output/070_trinotate/{run}/trinotate/trinotate_annotation_report.txt',
#                run=['merged', 'raw']),
#         expand('output/099_busco/{run}.{filter}/fixed_full_table.csv',
#                run=['merged', 'raw'],
#                filter=['length', 'expr']),
#         'output/099_busco/busco_plot.pdf',
#         expand('output/045_transcript-length/{run}/{result}',
#                run=['merged', 'raw'],
#                result=['blastx.outfmt6.grouped.w_pct_hit_length.txt',
#                        'blastx.outfmt6.w_pct_hit_length.txt'])

# # omg wtf busco
# rule plot_busco_results:
#     input:
#         busco_files = expand(
#             'output/099_busco/{run}.{filter}/fixed_full_table.csv',
#             run=['merged', 'raw'],
#             filter=['length', 'expr'])
#     output:
#         plot = 'output/099_busco/busco_plot.pdf'
#     log:
#         'output/logs/plot_busco_results.log'
#     singularity:
#         r
#     script:
#         'src/plot_busco_results.R'


# # busco is borked, have to fix the table
# rule fix_busco_table:
#     input:
#         'output/099_busco/{run}.{filter}/busco/run_hymenoptera_odb10/full_table.tsv'
#     output:
#         'output/099_busco/{run}.{filter}/fixed_full_table.csv'
#     log:
#         'output/logs/fixed_full_table.{run}.{filter}.log'
#     singularity:
#         r
#     script:
#         'src/fix_busco_table.R'

# rule busco:
#     input:
#         fasta = 'output/050_transcripts-by/{run}/{filter}.fasta',
#         lineage = 'data/hymenoptera_odb10'
#     output:
#         ('output/099_busco/{run}.{filter}/'
#          'busco/run_hymenoptera_odb10/'
#          'full_table.tsv'),
#     log:
#         Path(('output/logs/'
#               'busco.{run}.{filter}.log')).resolve()
#     params:
#         wd = 'output/099_busco/{run}.{filter}',
#         name = 'busco',
#         fasta = lambda wildcards, input: Path(input.fasta).resolve(),
#         lineage = lambda wildcards, input:
#             Path(input.lineage).resolve()
#     threads:
#         multiprocessing.cpu_count()
#     singularity:
#         busco
#     shell:
#         'cd {params.wd} || exit 1 ; '
#         'busco '
#         '--force '
#         '--in {params.fasta} '
#         '--out {params.name} '
#         '--lineage_dataset {params.lineage} '
#         '--cpu {threads} '
#         '--augustus_species honeybee1 '
#         '--mode transcriptome '
#         '&> {log}'


# # filter transcripts
# rule filter_isoforms:
#     input:
#         transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta',
#         names = 'output/050_transcripts-by/{run}/transcripts_by_{filter}.txt'
#     output:
#         'output/050_transcripts-by/{run}/{filter}.fasta'
#     log:
#         'output/logs/filter_isoforms.{run}.{filter}.log'
#     singularity:
#         bbduk
#     shell:
#         'filterbyname.sh '
#         'in={input.transcripts} '
#         'include=t '
#         'names={input.names} '
#         'out={output} '
#         '2> {log}'

# rule transcripts_by:
#     input:
#         qf = expand('output/040_trinity-abundance/{{run}}/{indiv}/quant.sf',
#                     indiv=ordered_indivs),
#         expr = ('output/040_trinity-abundance/{run}/'
#                 'salmon.isoform.TMM.EXPR.matrix'),
#         gtm = 'output/030_trinity/trinity_{run}/Trinity.fasta.gene_trans_map'
#     output:
#         ln = ('output/050_transcripts-by/{run}/'
#               'transcripts_by_length.txt'),
#         expr = ('output/050_transcripts-by/{run}/'
#                 'transcripts_by_expr.txt')
#     log:
#         'output/logs/transcripts_by.{run}.log'
#     singularity:
#         r
#     script:
#         'src/transcripts_by.R'

rule bc_target:
    input:
        expand('output/070_trinotate/{sample}/{run}/trinotate/Trinotate.sqlite',
               sample=all_samples,
               run=['raw', 'merged'])


# # annotation
rule trinotate:
    input:
        transcripts = 'output/030_trinity/trinity.{sample}.{run}/Trinity.fasta',
        blast_db = 'data/db/uniprot_sprot.pep',
        hmmer_db = 'data/db/Pfam-A.hmm',
        sqlite_db = 'data/db/Trinotate.sqlite',
    output:
        'output/070_trinotate/{sample}/{run}/trinotate/trinotate_annotation_report.txt',
        'output/070_trinotate/{sample}/{run}/blastx/blastx.outfmt6',
        'output/070_trinotate/{sample}/{run}/trinotate/Trinotate.sqlite'
    params:
        outdir = 'output/070_trinotate/{sample}/{run}'
    log:
        'output/logs/trinotate.{sample}.{run}.log'
    threads:
        10
    singularity:
        trinotate
    shell:
        'trinotate_pipeline '
        '--trinity_fasta {input.transcripts} '
        '--blast_db {input.blast_db} '
        '--hmmer_db {input.hmmer_db} '
        '--sqlite_db {input.sqlite_db} '
        '--outdir {params.outdir} '
        '--threads {threads} '
        '&> {log}'


# # differential expression
# rule annotate_de_results:
#     input:
#         de = 'output/060_de/{run}/de_wald_all.csv',
#         db = 'output/070_trinotate/{run}/trinotate/Trinotate.sqlite',
#         report = ('output/070_trinotate/{run}/trinotate/'
#                   'trinotate_annotation_report.txt')
#     output:
#         de_annot = 'output/060_de/{run}/de_wald_all_annot.csv',
#         de_pfam = 'output/060_de/{run}/de_wald_all_pfam.csv'
#     log:
#         'output/logs/annotate_de_results.{run}.log'
#     singularity:
#         bioconductor
#     script:
#         'src/annotate_de_results.R'

# rule de_wald_all:
#     input:
#         dds = 'output/060_de/{run}/dds.Rds'
#     output:
#         de_all = 'output/060_de/{run}/de_wald_all.csv'
#     params:
#         lfc_cutoff = 2,
#         alpha = 0.01
#     log:
#         'output/logs/de_wald_all.{run}.log'
#     threads:
#         multiprocessing.cpu_count()
#     singularity:
#         bioconductor
#     script:
#         'src/de_wald_all.R'

# rule plot_pca:
#     input:
#         vst_blind = 'output/060_de/{run}/vst_blind.Rds'
#     output:
#         plot = 'output/060_de/{run}/pca.pdf'
#     log:
#         'output/logs/plot_pca.{run}.log'
#     singularity:
#         bioconductor
#     script:
#         'src/plot_pca.R'


# rule generate_deseq_object:
#     input:
#         tx2gene = ('output/030_trinity/trinity_{run}/'
#                    'Trinity.fasta.gene_trans_map'),
#         quant = expand('output/040_trinity-abundance/{{run}}/{indiv}/quant.sf',
#                        indiv=ordered_indivs)
#     output:
#         dds = 'output/060_de/{run}/dds.Rds',
#         vst_blind = 'output/060_de/{run}/vst_blind.Rds'
#     log:
#         'output/logs/generate_deseq_object.{run}.log'
#     singularity:
#         bioconductor
#     script:
#         'src/generate_deseq_object.R'

# # analyse Trinity output
# rule group_blast_hits:
#     input:
#         blastx = 'output/070_trinotate/{run}/blastx/blastx.outfmt6',
#         db = 'data/db/uniprot_sprot.pep',
#         transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta'
#     output:
#         'output/045_transcript-length/{run}/blastx.outfmt6.grouped',
#         ('output/045_transcript-length/{run}/'
#          'blastx.outfmt6.grouped.w_pct_hit_length.txt')
#     params:
#         wd = 'output/045_transcript-length/{run}',
#         transcripts = lambda wildcards, input:
#             Path(input.transcripts).resolve(),
#         blastx = lambda wildcards, input:
#             Path(input.blastx).resolve(),
#         db = lambda wildcards, input:
#             Path(input.db).resolve()
#     log:
#         Path('output/logs/group_blast_hits.{run}.log').resolve()
#     singularity:
#         trinity
#     shell:
#         'cd {params.wd} || exit 1 ; '
#         'ln -s {params.blastx} blastx-group.outfmt6 ; '
#         'blast_outfmt6_group_segments.pl '
#         'blastx-group.outfmt6 '
#         '{params.transcripts} '
#         '{params.db} '
#         '> blastx.outfmt6.grouped '
#         '2> {log} ;'
#         'blast_outfmt6_group_segments.tophit_coverage.pl '
#         'blastx.outfmt6.grouped '
#         '> blastx.outfmt6.grouped.w_pct_hit_length.txt '
#         '2>> {log}'

# rule add_blast_coverage:
#     input:
#         blastx = 'output/070_trinotate/{run}/blastx/blastx.outfmt6',
#         db = 'data/db/uniprot_sprot.pep',
#         transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta'
#     output:
#         ('output/045_transcript-length/{run}/'
#          'blastx.outfmt6.w_pct_hit_length.txt'),
#         ('output/045_transcript-length/{run}/'
#          'blastx-coverage.outfmt6.w_pct_hit_length')
#     params:
#         wd = 'output/045_transcript-length/{run}',
#         transcripts = lambda wildcards, input:
#             Path(input.transcripts).resolve(),
#         blastx = lambda wildcards, input:
#             Path(input.blastx).resolve(),
#         db = lambda wildcards, input:
#             Path(input.db).resolve()
#     log:
#         Path('output/logs/add_blast_coverage.{run}.log').resolve()
#     singularity:
#         trinity
#     shell:
#         'cd {params.wd} || exit 1 ; '
#         'ln -s {params.blastx} blastx-coverage.outfmt6 ; '
#         'analyze_blastPlus_topHit_coverage.pl '
#         'blastx-coverage.outfmt6 '
#         '{params.transcripts} '
#         '{params.db} '
#         '> blastx.outfmt6.w_pct_hit_length.txt '
#         '2> {log}'

# rule plot_exn50:
#     input:
#         exn50 = 'output/040_trinity-abundance/{run}/ExN50.stats'
#     output:
#         plot = 'output/040_trinity-abundance/{run}/ExN50.pdf'
#     log:
#         'output/logs/plot_exn50.{run}.log'
#     singularity:
#         r
#     script:
#         'src/plot_exn50.R'

# rule contig_exn50:
#     input:
#         transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta',
#         expr = ('output/040_trinity-abundance/{run}/'
#                 'salmon.isoform.TMM.EXPR.matrix')
#     output:
#         inputs = ('output/040_trinity-abundance/{run}/'
#                   'salmon.isoform.TMM.EXPR.matrix.E-inputs'),
#         stats = ('output/040_trinity-abundance/{run}/'
#                  'ExN50.stats')
#     params:
#         outdir = 'output/040_trinity-abundance/{run}',
#         transcripts = lambda wildcards, input:
#             Path(input.transcripts).resolve(),
#         expr = lambda wildcards, input:
#             Path(input.expr).resolve()
#     log:
#         Path('output/logs/abundance_to_matrix.{run}.log').resolve()
#     singularity:
#         trinity
#     shell:
#         'cd {params.outdir} || exit 1 ; '
#         'contig_ExN50_statistic.pl '
#         '{params.expr} '
#         '{params.transcripts} '
#         '> ExN50.stats '
#         '2> {log}'


# rule abundance_to_matrix:
#     input:
#         qf = expand('output/040_trinity-abundance/{{run}}/{indiv}/quant.sf',
#                     indiv=ordered_indivs),
#         gtm = 'output/030_trinity/trinity_{run}/Trinity.fasta.gene_trans_map'
#     output:
#         'output/040_trinity-abundance/{run}/salmon.isoform.counts.matrix',
#         'output/040_trinity-abundance/{run}/salmon.isoform.TMM.EXPR.matrix'
#     params:
#         outdir = 'output/040_trinity-abundance/{run}',
#         qf = lambda wildcards, input:
#             sorted(set(Path(x).resolve() for x in input.qf)),
#         gtm = lambda wildcards, input:
#             Path(input.gtm).resolve()
#     log:
#         Path('output/logs/abundance_to_matrix.{run}.log').resolve()
#     singularity:
#         trinity
#     shell:
#         'cd {params.outdir} || exit 1 ; '
#         'abundance_estimates_to_matrix.pl '
#         '--est_method salmon '
#         '--gene_trans_map {params.gtm} '
#         '--name_sample_by_basedir '
#         '{params.qf} '
#         '&> {log}'

# rule trinity_abundance:
#     input:
#         'output/030_trinity/trinity_{run}/Trinity.fasta.salmon.idx',
#         transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta',
#         r1 = expand('output/010_reads/{indiv}_R1.fastq.gz',
#                     indiv=ordered_indivs),
#         r2 = expand('output/010_reads/{indiv}_R2.fastq.gz',
#                     indiv=indiv_with_r2),
#         samples_txt = 'output/030_trinity/samples.txt'
#     output:
#         expand('output/040_trinity-abundance/{{run}}/{indiv}/quant.sf',
#                indiv=ordered_indivs)
#     params:
#         outdir = 'output/040_trinity-abundance/{run}',
#         transcripts = lambda wildcards, input:
#             Path(input.transcripts).resolve(),
#         samples_txt = lambda wildcards, input:
#             Path(input.samples_txt).resolve()
#     log:
#         Path('output/logs/trinity_abundance.{run}.log').resolve()
#     threads:
#         multiprocessing.cpu_count()
#     singularity:
#         trinity
#     shell:
#         'cd {params.outdir} || exit 1 ; '
#         'align_and_estimate_abundance.pl '
#         '--transcripts {params.transcripts} '
#         '--seqType fq '
#         '--samples_file {params.samples_txt} '
#         '--est_method salmon '
#         '--SS_lib_type RF '
#         '--thread_count {threads} '
#         '--trinity_mode '
#         '&> {log}'


# rule trinity_abundance_prep:
#     input:
#         transcripts = 'output/030_trinity/trinity_{run}/Trinity.fasta'
#     output:
#         directory('output/030_trinity/trinity_{run}/Trinity.fasta.salmon.idx')
#     log:
#         'output/logs/trinity_abundance_prep.{run}.log'
#     singularity:
#         trinity
#     shell:
#         'align_and_estimate_abundance.pl '
#         '--transcripts {input.transcripts} '
#         '--est_method salmon '
#         '--trinity_mode '
#         '--prep_reference '
#         '&> {log}'

# rule generate_samples_text:
#     input:
#         raw_read_csv = raw_read_csv
#     output:
#         samples_txt = 'output/030_trinity/samples.txt'
#     singularity:
#         pandas_container
#     script:
#         'src/generate_samples_text.py'


# trinity
rule trinity:
    input:
        unpack(pick_trinity_input)
    output:
        'output/030_trinity/trinity.{sample}.{run}/Trinity.fasta',
        'output/030_trinity/trinity.{sample}.{run}/Trinity.fasta.gene_trans_map'
    params:
        outdir = 'output/030_trinity/trinity.{sample}.{run}'
    log:
        'output/logs/trinity.{sample}.{run}.log'
    threads:
        workflow.cores
    singularity:
        trinity
    shell:
        'Trinity '
        # '--FORCE '
        '--seqType fq '
        '--max_memory 800G '
        '--left {input.r1} '
        '--right {input.r2} '
        '--SS_lib_type RF '
        '--CPU {threads} '
        '--output {params.outdir} '
        '&> {log}'


# merge the input reads, try with and without
rule join_merged_with_r1:
    input:
        r1 = 'output/020_merged/{sample}_R1.fastq',
        merged = 'output/020_merged/{sample}_merged.fastq'
    output:
        temp('output/000_tmp/{sample}.r1_joined_with_merged.fastq')
    singularity:
        bbduk
    shell:
        'cat {input.r1} {input.merged} > {output}'

rule merge:
    input:
        r1 = 'output/010_reads/{sample}_R1.fastq',
        r2 = 'output/010_reads/{sample}_R2.fastq'
    output:
        merged = 'output/020_merged/{sample}_merged.fastq',
        r1 = 'output/020_merged/{sample}_R1.fastq',
        r2 = 'output/020_merged/{sample}_R2.fastq',
        ihist = 'output/020_merged/{sample}.ihist.txt'
    params:
        adaptors = adaptors
    log:
        'output/logs/merge.{sample}.log'
    threads:
        workflow.cores
    singularity:
        bbduk
    shell:
        'bbmerge.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.merged} '
        'outu1={output.r1} '
        'outu2={output.r2} '
        'ihist={output.ihist} '
        'verystrict=t '
        'adapters={params.adaptors} '
        '2>{log}'

# # trinity doesn't do interleaved
# rule dont_split:
#     input:
#         'output/000_tmp/{indiv}.trim.fastq'
#     output:
#         r1 = 'output/010_reads/{indiv}_R1.fastq.gz'
#     wildcard_constraints:
#         indiv = '|'.join(indiv_r1_only),
#     log:
#         'output/logs/split.{indiv}.log'
#     singularity:
#         bbduk
#     shell:
#         'reformat.sh '
#         'in={input} '
#         'int=f '
#         'out={output.r1} '
#         'zl=9 '
#         '2> {log}'

rule split:
    input:
        'output/000_tmp/{sample}.trim.fastq'
    output:
        r1 = 'output/010_reads/{sample}_R1.fastq',
        r2 = 'output/010_reads/{sample}_R2.fastq'
    log:
        'output/logs/split.{sample}.log'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'reformat.sh '
        'in={input} '
        'int=t '
        'out={output.r1} '
        'out2={output.r2} '
        '2> {log}'

rule trim:
    input:
        'output/000_tmp/{sample}.decon.fastq'
    output:
        temp('output/000_tmp/{sample}.trim.fastq')
    params:
        trim = adaptors
    log:
        log = 'output/logs/trim.{sample}.log',
        stats = 'output/logs/trim.{sample}.stats'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={log.stats} '
        '>> {output} '
        '2> {log.log} '

rule decon:
    input:
        'output/000_tmp/{sample}.repair.fastq'
    output:
        pipe('output/000_tmp/{sample}.decon.fastq')
    params:
        filter = '/phix174_ill.ref.fa.gz'
    log:
        log = 'output/logs/decon.{sample}.log',
        stats = 'output/logs/decon.{sample}.stats'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={log.stats} '
        '>> {output} '
        '2> {log.log} '

rule repair:
    input:
        'output/000_tmp/{sample}.barcode.fastq'
    output:
        # this should be a pipe but it won't let me
        temp('output/000_tmp/{sample}.repair.fastq')
    log:
        'output/logs/repair.{sample}.log'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'repair.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'repair=t '
        '>> {output} '
        '2> {log}'

rule check_barcodes:
    input:
        unpack(get_reads)
    output:
        pipe('output/000_tmp/{sample}.barcode.fastq')
    params:
        bc = lambda wildcards: get_bc(wildcards)
    log:
        'output/logs/check_barcodes.{sample}.log'
    group:
        'process'
    singularity:
        bbduk
    shell:
        'reformat.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'int=f '
        'out=stdout.fastq '
        'barcodefilter=t '
        'barcodes={params.bc} '
        '>> {output} '
        '2> {log}'
