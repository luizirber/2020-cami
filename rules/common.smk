### Functions for opal rules

def input_for_sample(w):
  if w.sample == "cami_i_low":
    return ["inputs/cami_i_low/RL_S001__insert_270.fq.gz"]
  elif w.sample == "cami_ii_mg":
    return expand("inputs/cami_ii_mg/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_{n}/reads/anonymous_reads.fq.gz", n=range(0,64))
  else:
    # TODO: fix
    return []

def datadir_for_sample(w, input):
  if w.sample == "cami_ii_mg":
    return os.path.join(*PosixPath(input.data[0]).parts[:3]),
  else:
    return os.path.join(*PosixPath(input.data[0]).parts[:2]),

def refseq_sigs(w):
  return expand("inputs/refseq/sigs/{sequence}.sig",
                sequence=(os.path.basename(g) for g in glob("inputs/refseq/sequences/*.fna.gz")))
