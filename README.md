# Formatting sourmash for CAMI evaluations

![.github/workflows/bench.yml](https://github.com/luizirber/2020-cami/workflows/.github/workflows/bench.yml/badge.svg)

This repo has the snakemake pipeline for running CAMI 2 evaluations (with [OPAL]),
including indices preparation based on the CAMI 2 databases.

The [sourmash biobox] definition is currently in [sourmash-bio/bioboxes](https://github.com/sourmash-bio/bioboxes)

[OPAL]: https://github.com/CAMI-challenge/OPAL/
[sourmash biobox]: https://quay.io/repository/sourmash.bio/sourmash

## Examples

The snakemake pipeline is continuously built by [GitHub Actions] on every change,
and OPAL reports are published to github pages. Here are the links to access them:
  - [CAMI I low complexity challenge dataset]: [Metalign results] from their repo,
    results for other tools come from the [OPAL repo].
  - [CAMI II mouse gut toy dataset]: I ran Metalign based on the [setup instructions]
    and [running with example data] wiki pages, and other tools results are from the
    official [data repository] for the CAMI 2 challenges.

[GitHub Actions]: https://github.com/luizirber/2020-cami/actions
[CAMI I low complexity challenge dataset]: https://luizirber.github.io/2020-cami/cami_i_low/opal_output_all/results.html
[CAMI II mouse gut toy dataset]: https://luizirber.github.io/2020-cami/cami_ii_mg/opal_output_all/results.html
[OPAL repo]: https://github.com/CAMI-challenge/OPAL/tree/3b24e5dc44ad3e6018d2381e6bf37ab89ba1a00a/data
[Metalign results]: https://github.com/nlapier2/Metalign/blob/378262587455fd5899a1de25995e28b09bc8d03b/paper_data/raw_results/cami1/metalign_results/metalign_default_cami1_low1.tsv
[setup instructions]: https://github.com/nlapier2/Metalign/wiki/Home/c1f422e246df20c11cab531ed412531bbe2b7028#overview
[running with example data]: https://github.com/nlapier2/Metalign/wiki/Running-some-example-data/ae546a62dc7adf384d2a9cb848d0dd8772cd9d1c#instructions
[data repository]: https://github.com/CAMI-challenge/data/blob/256460db625384a561aaf67ca168efd9ae070a52/CAMI2/toy/mouse_gut/taxonomic_profiling.tsv

## Other similar CAMI evaluation repos

These are other tools which participated in CAMI before and made their
submission process public.
 - [CCMetagen](https://github.com/vrmarcelino/CriticalAss2)
