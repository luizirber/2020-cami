# Formatting sourmash for CAMI evaluations

This repo has the snakemake pipeline for running CAMI 2 evaluations (with [OPAL]),
including indices preparation based on the CAMI 2 databases.

The [sourmash biobox] definition is currently in [Pull Request] to [docker_profiling_tools].

[OPAL]: https://github.com/CAMI-challenge/OPAL/
[sourmash biobox]: https://quay.io/repository/luizirber/sourmash_biobox
[Pull Request]: https://github.com/CAMI-challenge/docker_profiling_tools/pull/10
[docker_profiling_tools]: https://github.com/CAMI-challenge/docker_profiling_tools/

## Examples

The snakemake pipeline is continuously built by [GitHub Actions] on every change,
and OPAL reports are published to github pages. Here are the links to access them:
  - [CAMI I low complexity challenge dataset]
  - [CAMI II mouse gut toy dataset]

[CAMI I low complexity challenge dataset]: https://luizirber.github.io/2020-cami/cami_i_low/opal_output_all/results.html
[CAMI II mouse gut toy dataset]: https://luizirber.github.io/2020-cami/cami_ii_mg/opal_output_all/results.html
