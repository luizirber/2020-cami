#  * compute sigs for all samples (in /data/cami/rhizosphere_reads)
#  * run gather with LCA for each sample
#  * generate profile for each gather run
#  * concatenate profiles

rule all_cami2_rhimg:
  input:
    expand('outputs/cami_ii_rhimg/results/{exp}-k{ksize}.profile',
           ksize = (21, 31, 51),
           exp = ('short_read', 'long_read_nano', 'long_read_pacbio'))

rule compute_cami_sigs:
  output: 'outputs/cami_ii_rhimg/sigs/{exp}/{sample}.sig'
  input: '/data/cami/rhizosphere_reads/reads%2F{sample}.fq.gz'
  params:
    sample = "{sample}"
  shell: '''
    sourmash compute -k 21,31,51 \
                     --scaled 2000 \
                     --track-abundance \
                     --name {params.sample} \
                     -o {output} \
                     {input}
  '''

SCALED_FOR_EXP = {
  'short_read': 10000,
  'long_read_nano': 10000,
  'long_read_pacbio': 10000,
}

rule gather_cami_sample:
  output: 'outputs/cami_ii_rhimg/gather/{exp}-k{ksize}/{sample}.csv'
  input:
    sig = 'outputs/cami_ii_rhimg/sigs/{exp}/{sample}.sig',
    db = lambda w: 'outputs/lca/refseq-k{{ksize}}-s{scaled}.lca.json.gz'.format(scaled=SCALED_FOR_EXP[w.exp])
  threads: 2
  params:
    scaled = lambda w: SCALED_FOR_EXP[w.exp],
    ksize = "{ksize}"
  shell: '''
    sourmash gather \
        --scaled {params.scaled} \
        -k {params.ksize} \
        --output {output} \
        {input.sig} \
        {input.db}
  '''

rule gather_to_profile:
  output: 'outputs/cami_ii_rhimg/profile/{exp}-k{ksize}/{sample}.profile'
  input: 
    profile = 'outputs/cami_ii_rhimg/gather/{exp}-k{ksize}/{sample}.csv'
  params:
    sample = lambda w: w.sample.split("_")[-2]
  shell: '''
    python scripts/gather_to_opal.py profile \
      --output {output} \
      --taxdump inputs/ncbi_taxonomy \
      --taxid4index outputs/lca/taxid4index.csv \
      {params.sample} \
      {input.profile}
  '''

rule merge_profiles:
  output: 'outputs/cami_ii_rhimg/results/{exp}-k{ksize}.profile'
  input: expand('outputs/cami_ii_rhimg/profile/{{exp}}-k{{ksize}}/rhimgCAMI2_{{exp}}_sample_{i}_reads.profile', i=range(21))
  shell: '''
    cat {input} > {output}
  '''
