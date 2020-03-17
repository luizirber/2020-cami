### Prepare sourmash SBT using cami 2 refseq database

rule sourmash_compute:
  output: "inputs/refseq/sigs/{sequence}.sig"
  input: "inputs/refseq/sequences/{sequence}"
  shell: """
    sourmash compute -k 21,31,51 \
                     --scaled 2000 \
                     --track-abundance \
                     --name-from-first \
                     -o {output} \
                     {input}
  """

def refseq_sigs(w):
  return expand("inputs/refseq/sigs/{sequence}.sig",
                sequence=(os.path.basename(g) for g in glob("inputs/refseq/sequences/*.fna.gz")))

rule sbt_index:
  output: "outputs/sbt/refseq-k{ksize}.sbt.json"
  input:
    tar_file = "inputs/refseq/RefSeq_genomic_20190108.tar",
    sigs = refseq_sigs
  params:
    sigs_dir = lambda w, input: os.path.dirname(input.sigs[0]),
    ksize = lambda w: w.ksize,
    bfsize = 500000
  shell: """
    sourmash index -k {params.ksize} \
                   -d 2 \
                   -x {params.bfsize} \
                   --traverse-directory \
                   {output} {params.sigs_dir}
  """

### Prepare sourmash LCA using the SBT

rule download_lca_scripts:
  output:
    get_accession = "scripts/get-accessions-from-sbt.py",
    make_acc_taxid_mapping = "scripts/make-acc-taxid-mapping.py",
    make_lineage_csv = "scripts/make-lineage-csv.py",
    taxdump_utils = "scripts/ncbi_taxdump_utils.py",
  shell: """
    wget -qO {output.get_accession} https://raw.githubusercontent.com/dib-lab/2018-ncbi-lineages/63e8dc784af092293362f2e8e671ae03d1a84d1d/get-accessions-from-sbt.py
    chmod +x {output.get_accession}
    wget -qO {output.make_acc_taxid_mapping} https://raw.githubusercontent.com/dib-lab/2018-ncbi-lineages/63e8dc784af092293362f2e8e671ae03d1a84d1d/make-acc-taxid-mapping.py
    chmod +x {output.make_acc_taxid_mapping}
    wget -qO {output.make_lineage_csv} https://raw.githubusercontent.com/dib-lab/2018-ncbi-lineages/63e8dc784af092293362f2e8e671ae03d1a84d1d/make-lineage-csv.py
    chmod +x {output.make_lineage_csv}
    wget -qO {output.taxdump_utils} https://raw.githubusercontent.com/dib-lab/2018-ncbi-lineages/63e8dc784af092293362f2e8e671ae03d1a84d1d/ncbi_taxdump_utils.py
    chmod +x {output.taxdump_utils}
  """

rule lca_lineage_csv:
  output: "outputs/lca/refseq_lineage.csv"
  input:
    sbt = "outputs/sbt/refseq-k51.sbt.json",
    get_accession = "scripts/get-accessions-from-sbt.py",
    make_acc_taxid_mapping = "scripts/make-acc-taxid-mapping.py",
    make_lineage_csv = "scripts/make-lineage-csv.py",
    acc2taxid_gb = "inputs/ncbi_taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
    acc2taxid_wgs = "inputs/ncbi_taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
    taxonomy = expand("inputs/ncbi_taxonomy/{file}.dmp", file=('nodes', 'names')),
  shell: """
    {input.get_accession} {input.sbt} -o refseq.acc.txt
    {input.make_acc_taxid_mapping} refseq.acc.txt {input.acc2taxid_gb} {input.acc2taxid_wgs}
    {input.make_lineage_csv} {input.taxonomy} refseq.acc.txt.taxid -o {output}
    rm refseq.acc.txt refseq.acc.txt.taxid
  """

rule lca_index:
  output: "outputs/lca/refseq-k{ksize}-s{scaled}.lca.json.gz"
  input:
    sbt = "outputs/sbt/refseq-k{ksize}.sbt.json",
    lineage_csv = "outputs/lca/refseq_lineage.csv",
  params:
    sigs_dir = "inputs/refseq/sigs",
    ksize = "{ksize}",
    scaled = "{scaled}",
  shell: """
    sourmash lca index \
      {input.lineage_csv} \
      {output} \
      -k {params.ksize} \
      --scaled={params.scaled} \
      -C 3 \
      -f \
      --traverse-directory {params.sigs_dir} \
      --split-identifiers
  """

