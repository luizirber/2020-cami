from glob import glob
import os
from pathlib import PosixPath
import urllib

include: 'rules/common.smk'
include: 'rules/rhizo.smk'

rule all:
  input:
    "outputs/cami_i_low/opal_output/results.html",
    "outputs/cami_i_low/opal_output_all/results.html",
    "outputs/cami_ii_mg/opal_output/results.html",
    "outputs/cami_ii_mg/opal_output_all/results.html",
    #expand("outputs/lca/refseq-k{k}-s{scaled}.lca.json.gz", k=(21,31,51), scaled=(2000, 10000)),


### Include external rules for index construction
# Uncomment this for initial run and rebuilding indices,
# keeping it separated to avoid taking too long to recompute snakemake DAG
#include: "rules/build_indices.smk"

### Links for data

LCA_URLS = {
  'genbank': {
    '21': 'https://osf.io/d7rv8/download',
    '31': 'https://osf.io/4f8n3/download',
    '51': 'https://osf.io/nemkw/download'
  }
}

SBT_URLS = {
  'genbank': {
    '21': 'https://osf.io/nqs7k/download',
    '31': 'https://osf.io/h7k8a/download',
    '51': 'https://osf.io/jznwe/download'
  }
}

PROFILES = {
  # Data from OPAL repo: convenient, but not official
  "cami_i_low": {
    "data/metalign_profiles/metalign_default_cami1_low1.tsv": "Metalign",
    "data/opal/grave_wright_13": "Quikr",
    "data/opal/focused_archimedes_13": "MetaPhyler",
    "data/opal/furious_elion_13": "MP2.0",
    "data/opal/cranky_wozniak_13": "TIPP",
    "data/opal/agitated_blackwell_7": "CLARK",
    "data/opal/jolly_pasteur_3": "FOCUS",
    "data/opal/evil_darwin_13": "mOTU",
    #"data/metalign_profiles/metalign_sensitive_cami1_low1.tsv": "Metalign sensitive",
    #"data/metalign_profiles/metalign_precise_cami1_low1.tsv": "Metalign precise",
  },
  #'cami_i_low': {},
  'cami_ii_mg': {},
}

# I ran metalign myself, NOT OFFICIAL
PROFILES['cami_ii_mg']['data/metalign_profiles/metalign_default_cami_ii_mg.profile'] = "metalign"

# cami_ii_mg available profiles
for tool in ("motus2.5.1", "metaphlan2.9.21", "metaphlan2.2.0", "metapalette1.0.0",
             "motus1.1", "metaphyler1.25", "tipp2.0.0", "bracken2.5",
             "camiarkquikr1.0.0", "focus0.31"):
  path = f"data/cami_ii_mg_profiles/cami2_mouse_gut_{tool}.profile"
  PROFILES['cami_ii_mg'][path] = tool

# cami_i_low available profiles from official results

CAMI_I_LOW_PROFILES = {
  "CLARK_v1.1.3": "serene_almeida_reformatted",
#  "Common_Kmers_": "RL_S001__insert_270.fq-QC-default.profile",
#  "Common_Kmers_Sensitive_Unnormalized": "all.fq-QC-sensitive-unnormalized.profile",
  "commonkmers_sjanssen": "result_3.profile",
  "DUDes_": "RL_diginorm_0.1_k60-t1m0a0.000005-strain.out",
#  "DUDes_old": "RL_diginorm_subset10M_k10_profile.out",
#  "FOCUS_cfk7b": "cfk7b.out",
#  "FOCUS_cfk7bd": "cfk7bd.out",
#  "FOCUS_cfk7d": "cfk7d.out",
#  "FOCUS_cfk8b": "cfk8b.out",
#  "FOCUS_cfk8bd": "cfk8bd.out",
#  "FOCUS_cfk8d": "cfk8d.out",
  "FOCUS_sjanssen": "result_3.profile",
  "MetaPhlAn2.0_db_v20": "result_3_pairedend.txt.profile",
  "MetaPhyler_V1.25": "result_3.profile",
  "mOTU_1.1.1": "result_3.profile",
  "Quickr_sjanssen": "result_3.profile",
#  "Taxy-Pro_": "cami_low 2.profile",
  "Taxy-Pro_sjanssen": "result_3.profile",
  "TIPP_1.1": "result_3.profile",
}

#for (tool, path) in CAMI_I_LOW_PROFILES.items():
#  path = f"data/program_results/profiling/1st_CAMI_Challenge_Dataset_1_CAMI_low/Low/{tool}/{path}"
#  PROFILES['cami_i_low'][path] = tool

### Download CAMI databases

rule download_camiClient:
  output: "bin/camiClient.jar"
  shell: "wget -qO {output} https://data.cami-challenge.org/camiClient.jar"

rule download_taxonomy:
  output: "inputs/taxdump_cami2_toy.tar.gz"
  shell: "wget -qO {output} https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_DATABASES/taxdump_cami2_toy.tar.gz"

rule download_taxonomy_cami_i:
  output: "inputs/taxonomy_cami_i.tar.gz"
  shell: "wget -qO {output} ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100344/databases.dir/taxonomy.tar.gz"

rule extract_cami2_toy_taxonomy:
  output: expand("inputs/taxdump/{file}.dmp", file=('names', 'nodes'))
  input: "inputs/taxdump_cami2_toy.tar.gz"
  params:
    outdir = lambda w, input, output: os.path.dirname(output[0])
  shell: """
    mkdir -p {params.outdir}
    cd {params.outdir} && tar xf ../../{input}
  """

rule download_refseq_genomic:
  output: "inputs/refseq/RefSeq_genomic_20190108.tar"
  shell: """
    wget -qO {output[0]} https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_2_DATABASES/RefSeq_genomic_20190108.tar
    mkdir -p inputs/refseq/sequences
    cd inputs/refseq/sequences && tar xf ../RefSeq_genomic_20190108.tar
  """

rule download_cami2_taxonomy:
  output: "inputs/ncbi_taxonomy.tar"
  shell: """
    wget -qO {output[0]} https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_2_DATABASES/ncbi_taxonomy.tar
  """

rule extract_cami2_taxonomy:
  output: expand("inputs/ncbi_taxonomy/{file}.dmp", file=('names', 'nodes'))
  input: "inputs/ncbi_taxonomy.tar"
  params:
    infile = lambda w, input: os.path.basename(input[0])
  shell: """
    cd inputs && tar xf {params.infile}
    cd ncbi_taxonomy && tar xf taxdump.tar.gz
  """

rule download_and_extract_cami2_acc2taxid:
  output:
    compressed = "inputs/ncbi_taxonomy/accession2taxid/ncbi_taxonomy_accession2taxid.tar",
    gb = "inputs/ncbi_taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
    wgs = "inputs/ncbi_taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
  input:
    "inputs/ncbi_taxonomy/names.dmp",
    wgs = "db/nucl_wgs.accession2taxid.gz",
  shell: """
    wget -qO {output.compressed} https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_2_DATABASES/ncbi_taxonomy_accession2taxid.tar
    tar xf {output.compressed} -C inputs/ncbi_taxonomy/accession2taxid --strip-components 1
    cp {input.wgs} {output.wgs}
  """

### Metalign profiles

# available precisions: default, sensitive, precise
rule download_metalign_profiles:
  output: "data/metalign_profiles/metalign_{precision}_cami1_low1.tsv"
  shell: "wget -qO {output[0]} https://github.com/nlapier2/Metalign/raw/378262587455fd5899a1de25995e28b09bc8d03b/paper_data/raw_results/cami1/metalign_results/metalign_{wildcards.precision}_cami1_low1.tsv"

### Current acc-to-taxid mappings (genbank and wgs)

rule download_gb_acc2taxid:
  output: "db/nucl_gb.accession2taxid.gz"
  shell: "wget -qO {output} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"

rule download_wgs_acc2taxid:
  output: "db/nucl_wgs.accession2taxid.gz"
  shell: "wget -qO {output} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"

### Download CAMI gold standards, other tools profiles and biobox definitions

#### CAMI 2 data

rule download_gs_cami_ii_mg:
  output: "data/gs_cami_ii_mg.profile"
  shell: "wget -qO {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/gs_cami_ii_mg.profile"

rule download_biobox_cami_ii_mg:
  output: "data/biobox_cami_ii_mg.yaml"
  shell: "wget -qO {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/biobox_cami_ii_mg.yaml"

rule download_cami_ii_mg_summary:
  output: "data/summary_cami_ii_mg.tsv"
  shell: "wget -qO {output} https://github.com/CAMI-challenge/data/raw/256460db625384a561aaf67ca168efd9ae070a52/CAMI2/toy/mouse_gut/taxonomic_profiling.tsv"

rule download_cami_ii_mg_profiles:
  output:
    expand("data/cami_ii_mg_profiles/cami2_mouse_gut_{tool}.profile",
           tool=("bracken2.5", "camiarkquikr1.0.0", "focus0.31", "metapalette1.0.0",
                 "metaphlan2.2.0", "metaphlan2.9.21", "metaphyler1.25", "motus1.1",
                 "motus2.5.1", "tipp2.0.0")
    )
  input: "data/summary_cami_ii_mg.tsv"
  run:
    import pandas as pd
    dirpath = os.path.dirname(output[0])

    t = pd.read_table(input[0])
    for link in t['DirectLink']:
      outfile = os.path.basename(urllib.parse.urlparse(link).path)
      shell(f"wget -qO {dirpath}/{outfile} {link}")

rule download_cami_ii_mg:
  output: expand("inputs/cami_ii_mg/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_{n}/reads/anonymous_reads.fq.gz", n=range(0,64))
  input: "bin/camiClient.jar"
  shell: """
    mkdir -p inputs/cami_ii_mg && \
    java -jar {input} \
      -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMISIM_MOUSEGUT \
      inputs/cami_ii_mg \
      -p 'mousegut_scaffolds.*anonymous_reads.fq.gz' \
      -t 64
  """

rule download_cami_ii_marine:
  output:
    expand("inputs/CAMI_2_MARINE/marmgCAMI2_{type}_read_sample_{n}_reads.fq.gz",
        n=range(0,10),
        type=("short", "long"))
  input:
    client = "bin/camiClient.jar",
    linkfile = "data/marine.linkfile",
  shell: """
    mkdir -p inputs/CAMI_2_MARINE && \
    java -jar {input} \
      -d {input.linkfile} \
      inputs/CAMI_2_MARINE \
      -t 64
  """

rule download_cami_ii_strain:
  output:
    expand("inputs/CAMI_2_STRAIN/reads/strmgCAMI2_{type}_read_sample_{n}_reads.fq.gz",
        n=range(0,100),
        type=("short", "long"))
  input:
    client = "bin/camiClient.jar",
    linkfile = "data/strain.linkfile",
  shell: """
    mkdir -p inputs/CAMI_2_STRAIN && \
    java -jar {input} \
      -d {input.linkfile} \
      inputs/CAMI_2_STRAIN \
      -t 64
  """

### CAMI I data

rule download_gs_cami_i_hc:
  output: "data/gs_cami_i_hc.profile"
  shell: "wget -qO {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/gs_cami_i_hc.profile"

rule download_biobox_cami_i_hc:
  output: "data/biobox_cami_i_hc.yaml"
  shell: "wget -qO {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/biobox_cami_i_hc.yaml"

rule download_cami_i_hc:
  output: expand("inputs/cami_i_hc/RH_S00{n}__insert_270.fq.gz", n=range(1,6))
  input: "bin/camiClient.jar"
  params:
    output_dir = lambda w, output: os.path.dirname(output[0])
  shell: """
    mkdir -p {params.output_dir} && \
    java -jar bin/camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_HIGH {params.output_dir}
  """

rule download_cami_i_medium:
  output:
    expand("inputs/cami_i_medium/RM1_S00{n}__insert_5000.fq.gz", n=(1, 2)),
    expand("inputs/cami_i_medium/RM2_S00{n}__insert_270.fq.gz", n=(1, 2)),
  input: "bin/camiClient.jar"
  params:
    output_dir = lambda w, output: os.path.dirname(output[0])
  shell: """
    mkdir -p {params.output_dir} && \
    java -jar bin/camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_MEDIUM {params.output_dir}
  """

rule download_gs_cami_i_low:
  output: "data/gs_cami_i_low.profile"
  shell: "wget -qO {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/goldstandard_low_1.bin"

rule download_biobox_cami_i_low:
  output: "data/biobox_cami_i_low.yaml"
  # TODO: fix this
  shell: "wget -qO {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/biobox_cami_i_hc.yaml"

rule download_program_results_cami_i:
  output: "data/program_results.tar.gz"
  shell: "wget -qO {output} ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100344/program_results.tar.gz"

#rule extract_program_results_cami_i:
#  output: PROFILES['cami_i_low'].keys()
#  input: "data/program_results.tar.gz"
#  shell: """
#    mkdir -p data
#    cd data && tar xf ../{input} program_results/profiling
#  """

# TODO: fix sample names in cami_i_low profiles

rule download_cami_i_low:
  output: "inputs/cami_i_low/RL_S001__insert_270.fq.gz"
  input: "bin/camiClient.jar"
  params:
    output_dir = lambda w, output: os.path.dirname(output[0])
  shell: """
    mkdir -p {params.output_dir} && \
    java -jar bin/camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_LOW {params.output_dir} -p fq.gz
  """

### Download prepared sourmash databases

rule download_lca_database:
  output: "db/{db}-k{ksize}.lca.json.gz"
  params:
    url = lambda w: LCA_URLS[w.db][w.ksize]
  shell: "wget -qO {output[0]} {params.url}"

rule download_sbt_database:
  output:
    db="db/{db}-d2-k{ksize}.sbt.json",
    compressed="db/{db}-k{ksize}.tar.gz"
  params:
    url = lambda w: SBT_URLS[w.db][w.ksize],
    basename_compressed = lambda w, output: os.path.basename(output.compressed)
  shell: """
    wget -qO {output.compressed} {params.url} && \
    cd db && tar xf {params.basename_compressed}
  """

### snakemake rules to mirror opal workflow
### This is needed because opal workflow runs everything serially,
### and that takes too long...

rule profile_for_challenge_sample:
  output:
    "outputs/{challenge}/profiles/{sample}",
  input:
    #data = input_for_sample,
    biobox = "data/biobox_{challenge}.yaml",
    #db = "db/genbank-k51.lca.json.gz",
    #db = "outputs/sbt/refseq-k51.sbt.json",
    db = "outputs/lca/refseq-k51-s10000.lca.json.gz",
    acc2taxid_gb = "inputs/ncbi_taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
    acc2taxid_wgs = "inputs/ncbi_taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
    taxonomy = "inputs/ncbi_taxonomy/names.dmp"
  params:
    datadir = datadir_for_sample,
    outputdir = lambda w: f"outputs/{w.challenge}",
    dbdir = lambda w, input: os.path.dirname(input.db),
    taxdir = lambda w, input: os.path.dirname(input.taxonomy),
  run:
    ## TODO: make a new yaml just for the sample

    shell(f"""opal_stats.py \
      --input_dir $(pwd)/{params.datadir} \
      --output_dir $(pwd)/{params.outputdir} \
      --yaml $(pwd)/{input.biobox} \
      --volume $(pwd)/{params.taxdir}:/biobox/share/taxonomy/:ro \
      --volume $(pwd)/{params.dbdir}:/exchange/db:ro \
      quay.io/sourmash.bio/sourmash:latest
    """)

### Rules for opal workflow (calculating sourmash profiles)

rule run_opal_workflow:
  output:
    "outputs/{sample}/opal_output/results.html",
    "outputs/{sample}/quay.io-sourmash.bio-sourmash-latest/all_results.profile",
  input:
    data = input_for_sample,
    biobox = "data/biobox_{sample}.yaml",
    gs = "data/gs_{sample}.profile",
    #db = "db/genbank-k51.lca.json.gz",
    #db = "outputs/sbt/refseq-k51.sbt.json",
    db = "outputs/lca/refseq-k51-s10000.lca.json.gz",
    acc2taxid_gb = "inputs/ncbi_taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
    acc2taxid_wgs = "inputs/ncbi_taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
    taxid4index = "outputs/lca/taxid4index.csv",
    taxonomy = "inputs/ncbi_taxonomy/names.dmp",
  params:
    datadir = datadir_for_sample,
    outputdir = lambda w: f"outputs/{w.sample}",
    dbdir = lambda w, input: os.path.dirname(input.db),
    taxdir = lambda w, input: os.path.dirname(input.taxonomy),
  shell: """
    opal_workflow.py \
    quay.io/sourmash.bio/sourmash:latest \
    --labels "sourmash" \
    --input_dir $(pwd)/{params.datadir} \
    --output_dir $(pwd)/{params.outputdir} \
    --yaml $(pwd)/{input.biobox} \
    --volume $(pwd)/{params.taxdir}:/biobox/share/taxonomy/:ro \
    --volume $(pwd)/{params.dbdir}:/exchange/db:ro \
    --gold_standard_file $(pwd)/{input.gs} \
    --plot_abundances \
    --desc "{wildcards.sample}"
  """

### Generate an OPAL report comparing with other tools

rule run_opal_report:
  output:
    "outputs/{sample}/opal_output_all/results.html",
  input:
    "outputs/{sample}/quay.io-sourmash.bio-sourmash-latest/all_results.profile",
    gs = "data/gs_{sample}.profile",
    profiles = lambda w: list(PROFILES[w.sample].keys())
  params:
    outputdir = lambda w: f"outputs/{w.sample}/opal_output_all/",
    desc = "{sample}",
    labels = lambda w, input: "sourmash," + ",".join(PROFILES[w.sample][path] for path in input.profiles)
  shell: """
    opal.py \
    --gold_standard_file $(pwd)/{input.gs} \
    --output_dir $(pwd)/{params.outputdir} \
    --plot_abundances \
    --desc='{params.desc}' \
    -l '{params.labels}' \
    {input[0]} \
    {input.profiles}
  """
