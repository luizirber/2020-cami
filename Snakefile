import os

LCA_URLS = {                                                                                                                                                                                                                                             
  'genbank': {                                                                                                                                                                                                                                           
    '21': 'https://osf.io/d7rv8/download',                                                                                                                                                                                                               
    '31': 'https://osf.io/4f8n3/download',                                                                                                                                                                                                               
    '51': 'https://osf.io/nemkw/download'                                                                                                                                                                                                                
  }                                                                                                                                                                                                                                                      
}  

rule all:
  input:
    "outputs/cami_i_low/opal_output/results.html",
    "outputs/cami_i_hc/opal_output/results.html"

rule download_taxonomy:
  output: "inputs/taxdump_cami2_toy.tar.gz"
  shell: "wget -O {output} https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_DATABASES/taxdump_cami2_toy.tar.gz"

rule extract_taxonomy:
  output: expand("inputs/taxdump/{file}.dmp", file=('names', 'nodes'))
  input: "inputs/taxdump_cami2_toy.tar.gz"
  params:
    outdir = lambda w, input, output: os.path.dirname(output[0])
  shell: """
    mkdir -p {params.outdir}
    cd {params.outdir} && tar xf ../../{input}
  """

rule download_gb_acc2taxid:
  output: "db/nucl_gb.accession2taxid.gz"
  shell: "wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"

rule download_wgs_acc2taxid:
  output: "db/nucl_wgs.accession2taxid.gz"
  shell: "wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz"

rule download_camiClient:
  output: "bin/camiClient.jar"
  shell: "wget -O {output} https://data.cami-challenge.org/camiClient.jar"

rule download_gs_cami_ii_mg:
  output: "data/gs_cami_ii_mg.profile"
  shell: "wget -O {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/gs_cami_ii_mg.profile"

rule download_biobox_cami_ii_mg:
  output: "data/biobox_cami_ii_mg.yaml"
  shell: "wget -O {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/biobox_cami_ii_mg.yaml"

rule download_gs_cami_i_hc:
  output: "data/gs_cami_i_hc.profile"
  shell: "wget -O {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/gs_cami_i_hc.profile"

rule download_biobox_cami_i_hc:
  output: "data/biobox_cami_i_hc.yaml"
  shell: "wget -O {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/biobox_cami_i_hc.yaml"

rule download_cami_i_hc:
  output: expand("inputs/cami_i_hc/RH_S00{n}__insert_270.fq.gz", n=range(1,6))
  input: "bin/camiClient.jar"
  params:
    output_dir = lambda w, output: os.path.dirname(output[0])
  shell: """
    mkdir -p {params.output_dir} && \
    java -jar bin/camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_HIGH {params.output_dir}
  """

rule download_gs_cami_i_low:
  output: "data/gs_cami_i_low.profile"
  shell: "wget -O {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/goldstandard_low_1.bin"

rule download_biobox_cami_i_low:
  output: "data/biobox_cami_i_low.yaml"
  # TODO: fix this
  shell: "wget -O {output} https://github.com/CAMI-challenge/OPAL/raw/master/data/biobox_cami_i_hc.yaml"

rule download_cami_i_low:
  output: "inputs/cami_i_low/RL_S001__insert_270.fq.gz"
  input: "bin/camiClient.jar"
  params:
    output_dir = lambda w, output: os.path.dirname(output[0])
  shell: """
    mkdir -p {params.output_dir} && \
    java -jar bin/camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_LOW {params.output_dir} -p fq.gz
  """

rule download_lca_database:                                                                                                                                                                                                                              
  output: "db/{db}-k{ksize}.lca.json.gz"                                                                                                                                                                                                         
  params:                                                                                                                                                                                                                                                
    url = lambda w: LCA_URLS[w.db][w.ksize]                                                                                                                                                                                                              
  shell: "curl -L -o {output[0]} {params.url}"   

rule run_opal_workflow:
  output:
    "outputs/{sample}/opal_output/results.html",
    "outputs/{sample}/sourmash_bioboxes-latest/all_results.profile"
  input:
    data = "inputs/{sample}/",
    biobox = "data/biobox_{sample}.yaml",
    gs = "data/gs_{sample}.profile",
    db = "db/genbank-k51.lca.json.gz",
    acc2taxid_gb = "db/nucl_gb.accession2taxid.gz",
    acc2taxid_wgs = "db/nucl_wgs.accession2taxid.gz",
    taxonomy = "inputs/taxdump/names.dmp"
  params:
    datadir = lambda w, input: os.path.dirname(input.data),
    outputdir = lambda w: f"outputs/{w.sample}",
    dbdir = lambda w, input: os.path.dirname(input.db),
    taxdir = lambda w, input: os.path.dirname(input.taxonomy),
  shell: """
    opal_workflow.py \
    sourmash_bioboxes:latest \
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

rule run_opal_report:
  output:
    "outputs/{sample}/opal_output_all/results.html",
  input:
    "outputs/{sample}/sourmash_bioboxes-latest/all_results.profile",
    gs = "data/gs_{sample}.profile",
  params:
    outputdir = lambda w: f"outputs/{w.sample}_all",
  shell: """
    opal.py \
    --gold_standard_file $(pwd)/{input.gs} \
    --output_dir $(pwd)/{params.outputdir} \
    --plot_abundances \
    --desc='1st CAMI low' \
    -l "sourmash, TIPP, Quikr, MP2.0, MetaPhyler, mOTU, CLARK, FOCUS"
    {input[0]} \
    data/opal/cranky_wozniak_13 \
    data/opal/grave_wright_13 \
    data/opal/furious_elion_13 \
    data/opal/focused_archimedes_13 \
    data/opal/evil_darwin_13 \
    data/opal/agitated_blackwell_7 \
    data/opal/jolly_pasteur_3
  """
