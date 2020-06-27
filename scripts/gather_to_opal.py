#! /usr/bin/env python3
"""
Take a gather CSV and one or more NCBI 'accession2taxid' files
and create
  1) csv containing accessions, taxid, and
  2) csv with lineage

run:
```
python gather-to-opal.py example_output example_output.csv \
    --acc2taxid nucl_gb.accession2taxid.gz \
    --acc2taxid nucl_wgs.accession2taxid.gz \
    --output example_output.profile
```
"""

from __future__ import print_function

import argparse
import gzip
import os

import pandas as pd
import taxonomy

__version__ = "0.1.0"


def taxid4index(index, acc2taxids, output):
    index = load_index(index)

    acc_set = set()
    for sig in index.signatures():
        acc = sig.name().split()[0].split(".")[0]
        acc_set.add(acc)

        if acc.startswith("NZ_") or acc.startswith("NC_"):
            acc_set.add(acc[3:])

    matches = matches_for_acc2taxids(acc_set, acc2taxids)

    if output:
        sr = pd.Series(matches)
        sr.columns = ("acc", "taxid")
        sr.to_csv(output)

    return matches


def load_index(filename):
    import sourmash

    index = None
    try:
        index = sourmash.load_sbt_index(filename)
    except (ValueError, EnvironmentError):
        pass

    if index is None:
        try:
            index = sourmash.lca.lca_utils.LCA_Database()
            index.load(filename)
        except (ValueError, EnvironmentError, TypeError):
            pass

    if index is None:
        # TODO: raise error
        pass

    return index


def matches_for_acc2taxids(acc_set, acc2taxids):
    m = 0
    matches = {}
    for filename in acc2taxids:
        if not acc_set:
            break

        xopen = open
        if filename.endswith(".gz"):
            xopen = gzip.open

        with xopen(filename, "rt") as fp:  # ruun through acc2taxid files
            next(fp)  # skip headers
            for n, line in enumerate(fp):
                if not acc_set:
                    break

                if n and n % 1000000 == 0:
                    print("\r\033[K", end="")
                    print(
                        "... read {} lines of {}; found {} of {}".format(
                            n, filename, m, m + len(acc_set)
                        ),
                        end="\r",
                    )

                try:
                    acc, acc_version, taxid, _ = [l.strip() for l in line.split()]
                except ValueError:
                    print("ignoring line", (line,))
                    continue

                if acc in acc_set:
                    m += 1
                    matches[acc] = str(taxid)
                    acc_set.remove(acc)

                    if not acc_set:
                        break

    if acc_set:
        print("failed to find {} acc: {}".format(len(acc_set), acc_set))
    else:
        print("found all {} accessions!".format(m))

    return matches


def get_taxid(gather_csv, acc2taxid_files, taxid4index):
    gather_info = pd.read_csv(gather_csv)
    # grab the acc from gather column `name`
    gather_info["accession"] = gather_info["name"].str.replace(r"\..*", "")
    gather_info["percentage"] = gather_info["f_unique_weighted"] * 100

    # init opal_info df
    opal_info = gather_info.loc[:, ["accession", "percentage"]].copy()
    opal_info.set_index("accession", inplace=True)
    acc_set = set(opal_info.index)

    if taxid4index:
        taxids = pd.read_csv(taxid4index, names=("acc", "taxid"))
        taxids.set_index("acc", inplace=True)
        matches = {}
        not_found = set()
        for acc in acc_set:
            taxid = None
            try:
                taxid = taxids.loc[acc]
            except KeyError:
                not_found.add(acc)
            else:
                matches[acc] = taxid
        if not_found:
            print("failed to find {} acc: {}".format(len(not_found), not_found))
        else:
            print("found all {} accessions!".format(len(matches)))
    else:
        matches = matches_for_acc2taxids(acc_set, acc2taxid_files)

    for acc, taxid in matches.items():
        opal_info.loc[acc, "taxid"] = str(taxid["taxid"])

    return opal_info


def get_row_taxpath(row, taxo, ranks):
    current_taxid = str(row["taxid"])
    try:
        current_rank = taxo.rank(current_taxid)
    except Exception as e:
        # pyo3 doesn't export Exceptions properly, so doing this for now...
        if "TaxonomyError" in repr(e):
            print(e)
            # a TaxonomyError is raised when the taxid is not in the taxonomy.
            # Skipping row for now, and not setting 'taxpath' means it needs to
            # be filtered out later.
            return
        else:
            raise

    if current_rank == "no rank":
        # need to figure out based on parent
        parent = taxo.parent(current_taxid)
        parent_rank = taxo.rank(parent)

        if parent_rank == "species":
            # we have a strain
            row["rank"] = "strain"
        elif parent_rank == "genus":
            # it might be a species-like rank,
            # but we should leave empty for OPAL
            row["rank"] = "no rank"
        elif parent_rank == "subspecies":
            # we have a strain
            row["rank"] = "strain"
        elif parent_rank == "no rank":
            # seen this with taxid 1620419, which has lineage
            # no rank|no rank|subspecies|species|genus|family|order|class|phylum|superkingdom
            parent = taxo.parent(parent)
            parent_rank = taxo.rank(parent)
            if parent_rank in ("species", "subspecies"):
                row["rank"] = "strain"
            elif parent_rank == "no rank":
                # taxid 1194420 has lineage
                # no rank|no rank|no rank|species|species subgroup|species group|genus|family|order|class|phylum|superkingdom
                parent = taxo.parent(parent)
                parent_rank = taxo.rank(parent)
                if parent_rank == "species":
                    row["rank"] = "strain"
                else:
                    raise Exception("TODO other ranks testing for parent")
            else:
                raise Exception("TODO other ranks testing for parent")
    elif current_rank == "species":
        row["rank"] = "species"
    elif current_rank == "subspecies":
        # TODO: should probably drop it from profile if subspecies not in ranks,
        # but need to do it later? (can't remove row during apply)
        row["rank"] = "no rank"
    else:
        raise Exception("TODO other ranks testing for current")

    # uses taxonomy pkg
    lineage = taxo.lineage(current_taxid)
    valid_ranks = {}

    if row["rank"] == "strain" and "strain" in ranks:
        valid_ranks["strain"] = current_taxid

    for l in lineage:
        current_rank = taxo.rank(l)
        if current_rank in ranks:
            valid_ranks[current_rank] = l

    final_lineage = []
    for rank in ranks:
        if rank in valid_ranks:
            final_lineage.append(valid_ranks[rank])
        else:
            final_lineage.append("")

    row["taxpath"] = "|".join(final_lineage)
    return row


def summarize_all_levels(df, ranks):
    new_rows = []
    for (percentage, tax_id, rank, taxpath) in df.itertuples(index=False, name=None):
        if not taxpath:
            # no taxpath means tax_id couldn't be found in the taxonomy,
            # and so it's empty. Dropping it for now.
            continue

        lineage_values = taxpath.split("|")
        for i, (rank, tax_id) in enumerate(zip(ranks, lineage_values), 1):
            if not tax_id:
                continue

            taxpath = "|".join(lineage_values[:i])
            new_rows.append([percentage, tax_id, rank, taxpath])

    new_df = pd.DataFrame(new_rows, columns=df.columns)
    return new_df.groupby(["taxid", "rank", "taxpath"], as_index=False).sum()


def gen_report(sample_id, ranks, taxons, *, taxonomy_id=None, program=None):
    output_lines = f"""# Taxonomic Profiling Output
@SampleID:{sample_id}
@Version:0.9.1
@Ranks:{ranks}
""".splitlines()

    if taxonomy_id is not None:
        output_lines.append(f"@TaxonomyID:{taxonomy_id}")
#    if program is not None:
#        output_lines.append(f"@__program__: {program}")
    output_lines.append(f"@@TAXID\tRANK\tTAXPATH\tPERCENTAGE")

    for tax in taxons.itertuples(index=False, name=None):
        tax_line = "\t".join(str(t) for t in tax)
        output_lines.append(tax_line)

    return "\n".join(output_lines)


def gather_to_opal(
    sample_id,
    gather_csv,
    taxdump,
    tax_ranks,
    *,
    acc2taxid=None,
    taxid4index=None,
    opal_csv=None,
    taxonomy_id=None,
    program=None,
):
    opal_info = get_taxid(gather_csv, acc2taxid, taxid4index)

    # Drop tax_ids not found. There is a warning already on the `get_taxid`
    # function, but might want to be more eloquent...
    opal_info.dropna(subset=["taxid"], inplace=True)

    # load ncbi taxonomy info
    taxo = taxonomy.Taxonomy.from_ncbi(
        os.path.join(taxdump, "nodes.dmp"), os.path.join(taxdump, "names.dmp")
    )

    # get lineage using taxid
    tax_df = opal_info.astype(object).apply(
        lambda row: get_row_taxpath(row, taxo, tax_ranks), axis=1
    )

    # summarize taxonomic ranks
    rank_df = summarize_all_levels(tax_df, tax_ranks)
    if not opal_csv:
        opal_csv = gather_csv.rsplit(".csv")[0] + "_opal.csv"

    out = gen_report(
        sample_id,
        "|".join(tax_ranks),
        rank_df,
        taxonomy_id=taxonomy_id,
        program=program,
    )
    with open(opal_csv, "w") as f:
        f.write(out)
        f.write('\n')


def gather_to_cami_main(args):
    gather_to_opal(
        args.sample_id,
        args.gather_csv,
        args.taxdump,
        args.ranks.split("|"),
        taxid4index=args.taxid4index,
        acc2taxid=args.acc2taxid,
        opal_csv=args.output,
        taxonomy_id=args.taxonomy_id,
        program=args.program,
    )


def tx4idx_main(args):
    taxid4index(args.index, args.acc2taxid, args.output)


def main():
    parser = argparse.ArgumentParser()
    subp = parser.add_subparsers()

    profile = subp.add_parser("profile")
    profile.add_argument("sample_id")
    profile.add_argument("gather_csv")
    profile.add_argument("--acc2taxid", action="append")
    profile.add_argument("--taxid4index", type=str)
    profile.add_argument("--taxdump", default="taxdump")
    profile.add_argument(
        "--ranks", default="superkingdom|phylum|class|order|family|genus|species|strain"
    )
    profile.add_argument("-o", "--output")
    profile.add_argument("--taxonomy_id")
    profile.add_argument("--program", default="sourmash gather")
    profile.set_defaults(func=gather_to_cami_main)

    tx4idx = subp.add_parser("taxid4index")
    tx4idx.add_argument("index")
    tx4idx.add_argument("--acc2taxid", action="append")
    tx4idx.add_argument("--output")
    tx4idx.set_defaults(func=tx4idx_main)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
