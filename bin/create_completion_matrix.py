#!/usr/bin/env python

import argparse

import pythoncyc
from pythoncyc import PTools as PTools
from pythoncyc.PTools import PToolsError as PToolsError
from pythoncyc.PTools import PythonCycError as PythonCycError


def get_genes_of_reaction(rxn, pgdb):
    genes = set()
    for gene in pgdb.genes_of_reaction(rxn):
        genes.add(gene.split("|")[1])
    for enz in pgdb.enzymes_of_reaction(rxn):
        for gene in pgdb.genes_of_protein(enz):
            genes.add(gene.split("|")[1])
    return genes


def get_pathways_none_spontaneous_reactions(pgdb):
    pathways = dict()
    for path in pgdb.all_pathways(selector="all", base=True):
        pathways[path] = dict()
        pathways[path]["Name"] = pgdb[path].common_name
        pathways[path]["Reactions"] = set()
        for rxn in pgdb[path].reaction_list:
            if pgdb[rxn].spontaneous_p is not None:
                if not pgdb[rxn].spontaneous_p:
                    pathways[path]["Reactions"].add(rxn)
            else:
                pathways[path]["Reactions"].add(rxn)
    return pathways


def get_pathways_none_spontaneous_orphan_reactions(pgdb):
    pathways = dict()
    for path in pgdb.all_pathways(selector="all", base=True):
        pathways[path] = dict()
        pathways[path]["Name"] = pgdb[path].common_name
        pathways[path]["Reactions"] = set()
        for rxn in pgdb[path].reaction_list:
            is_spontaneous = False
            is_orphan = None
            is_orphan_in_metacyc = False

            if pgdb[rxn].spontaneous_p is not None and pgdb[rxn].spontaneous_p:
                is_spontaneous = True

            if not is_spontaneous:
                if len(get_genes_of_reaction(rxn, pgdb)) == 0:
                    is_orphan_in_metacyc = True

                if pgdb[rxn].orphan_p is not None and pgdb[rxn].orphan_p[0] == "|NO|":
                    is_orphan = False
                else:
                    is_orphan = True

            if (
                not (is_orphan_in_metacyc and (is_orphan is None or is_orphan))
                and not is_spontaneous
            ):
                pathways[path]["Reactions"].add(rxn)

    return pathways


def get_reactions_in_pangenome(pgdb):
    reactions = set()
    for rxn in pgdb.all_rxns(type_of_reactions=":all"):
        if len(get_genes_of_reaction(rxn, pgdb)) != 0:
            reactions.add(rxn)
    return reactions


def get_reactions_with_families(pgdb):
    reactions_with_families = dict()
    for rxn in pgdb.all_rxns(type_of_reactions=":all"):
        families = get_genes_of_reaction(rxn, pgdb)
        if len(families) != 0:
            reactions_with_families[rxn] = set()
            reactions_with_families[rxn].update(families)
    return reactions_with_families


def get_pathways(pgdb):
    pathways = set()
    for path in pgdb.all_pathways(selector="all", base=True):
        pathways.add(path)
    return pathways


def write_pgdb_reaction_presence_absence_by_strain(
    pgdb, strains_with_families, partitions_with_families, meta, pgdb_name
):
    file_name = pgdb_name + "_reaction_presence_absence.tsv"

    _pgdb_reactions = get_reactions_in_pangenome(pgdb)
    reactions_with_families = get_reactions_with_families(pgdb)

    with open(file_name, "w") as pgdb_write:
        header = "PGDB\tReaction\tFamilies\tPersistent %\tShell %\tCloud %"
        for strain in strains_with_families:
            header += "\t" + strain
        pgdb_write.write(header + "\n")

        for rxn in reactions_with_families:
            families = reactions_with_families[rxn]
            families_to_str = ",".join(families)

            persistent_pct = str(
                len(families.intersection(partitions_with_families["persistent"]))
                * 100
                / len(families)
            )
            shell_pct = str(
                len(families.intersection(partitions_with_families["shell"]))
                * 100
                / len(families)
            )
            cloud_pct = str(
                len(families.intersection(partitions_with_families["cloud"]))
                * 100
                / len(families)
            )

            line = (
                pgdb_name
                + "\t"
                + rxn.split("|")[1]
                + "\t"
                + families_to_str
                + "\t"
                + persistent_pct
                + "\t"
                + shell_pct
                + "\t"
                + cloud_pct
            )
            for strain in strains_with_families:
                if len(families.intersection(strains_with_families[strain])) >= 1:
                    line += "\t1"
                else:
                    line += "\t0"
            pgdb_write.write(line + "\n")


def write_pgdb_pathway_completion_by_strain(
    pgdb,
    strains_with_families,
    partitions_with_families,
    modules_with_families,
    meta,
    pgdb_name,
    use_orphan,
):
    if use_orphan:
        meta_pathways = get_pathways_none_spontaneous_reactions(meta)
        file_name = pgdb_name + "_pathway_completion_by_strain.tsv"
    else:
        meta_pathways = get_pathways_none_spontaneous_orphan_reactions(meta)
        file_name = pgdb_name + "_pathway_completion_wo_orphan_by_strain.tsv"

    pgdb_reactions = get_reactions_in_pangenome(pgdb)
    pgdb_pathways = get_pathways(pgdb)
    reactions_with_families = get_reactions_with_families(pgdb)

    with open(file_name, "w") as pgdb_write:
        header = "PGDB\tPathway\tPathway name\tNb reactions\tFamilies\tPersistent %\tShell %\tCloud %\tModules (reaction cov.)\tGlobal completion\tMax completion"
        for strain in strains_with_families:
            header += "\t" + strain
        pgdb_write.write(header + "\n")
        for path in pgdb_pathways:
            global_completion = len(
                meta_pathways[path]["Reactions"].intersection(pgdb_reactions)
            ) / len(meta_pathways[path]["Reactions"])

            families = set()
            for rxn in meta_pathways[path]["Reactions"]:
                if rxn in reactions_with_families:
                    families.update(reactions_with_families[rxn])
            families_to_str = ",".join(families)

            persistent_pct = str(
                len(families.intersection(partitions_with_families["persistent"]))
                * 100
                / len(families)
            )
            shell_pct = str(
                len(families.intersection(partitions_with_families["shell"]))
                * 100
                / len(families)
            )
            cloud_pct = str(
                len(families.intersection(partitions_with_families["cloud"]))
                * 100
                / len(families)
            )

            modules_to_str = ""
            for module in modules_with_families:
                nb_common_fam = len(
                    families.intersection(modules_with_families[module])
                )
                if nb_common_fam >= 2:
                    nb_rxn_in_module = 0
                    nb_rxn_with_families = 0
                    for rxn in meta_pathways[path]["Reactions"]:
                        if rxn in reactions_with_families:
                            nb_rxn_with_families += 1
                            if (
                                len(
                                    reactions_with_families[rxn].intersection(
                                        modules_with_families[module]
                                    )
                                )
                                >= 1
                            ):
                                nb_rxn_in_module += 1
                    modules_to_str += (
                        module
                        + " ("
                        + str(nb_rxn_in_module / nb_rxn_with_families)
                        + ") ,"
                    )
            modules_to_str = modules_to_str.rstrip(" ,")

            line = (
                pgdb_name
                + "\t"
                + path.split("|")[1]
                + "\t"
                + meta_pathways[path]["Name"]
                + "\t"
                + str(len(meta_pathways[path]["Reactions"]))
                + "\t"
                + families_to_str
                + "\t"
                + persistent_pct
                + "\t"
                + shell_pct
                + "\t"
                + cloud_pct
                + "\t"
                + modules_to_str
                + "\t"
                + str(global_completion)
            )

            strain_completion = dict()
            for strain in strains_with_families:
                nbreactions_in_strain = 0
                for rxn in meta_pathways[path]["Reactions"]:
                    if (
                        rxn in reactions_with_families
                        and len(
                            reactions_with_families[rxn].intersection(
                                strains_with_families[strain]
                            )
                        )
                        >= 1
                    ):
                        nbreactions_in_strain += 1
                strain_completion[strain] = nbreactions_in_strain / len(
                    meta_pathways[path]["Reactions"]
                )
            max_completion = max(strain_completion.values())
            line += "\t" + str(max_completion)
            for strain in strain_completion:
                line += "\t" + str(strain_completion[strain])
            pgdb_write.write(line + "\n")


def read_pangenome_rtab(rtab):
    rtabfile = open(rtab, "r")
    header = rtabfile.readline().rstrip()
    strains = header.split()[1:]
    strains_with_families = dict()

    for strain in strains:
        strains_with_families[strain] = set()

    for line in rtabfile:
        i = 0
        for field in line.rstrip().split("\t"):
            if i == 0:
                famid = field.upper()
            elif field == "1":
                strains_with_families[strains[i - 1]].add(famid)
            i += 1

    return strains_with_families


def read_functional_module_file(module_file_name):

    modules_with_families = dict()

    modulefile = open(module_file_name, "r")

    _header = modulefile.readline().rstrip()

    for line in modulefile:
        (moduleid, famid) = line.rstrip().split("\t")
        if moduleid not in modules_with_families:
            modules_with_families[moduleid] = set()
        modules_with_families[moduleid].add(famid.upper())

    return modules_with_families


def read_partition_files(persistent_file_name, shell_file_name, cloud_file_name):

    partitions_with_families = dict()

    with open(persistent_file_name) as fh:
        partitions_with_families["persistent"] = {line.rstrip().upper() for line in fh}

    with open(shell_file_name) as fh:
        partitions_with_families["shell"] = {line.rstrip().upper() for line in fh}

    with open(cloud_file_name) as fh:
        partitions_with_families["cloud"] = {line.rstrip().upper() for line in fh}

    return partitions_with_families


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-pgdb", help="Name of the PGDB (example : eco)", required=True, type=str
    )
    parser.add_argument(
        "-ppanggolin_dir",
        help="Directory of PPanGGOLiN result files",
        required=True,
        type=str,
    )
    args = parser.parse_args()

    pgdb_name = args.pgdb
    ppanggolin_dir = args.ppanggolin_dir
    rtab = ppanggolin_dir + "/gene_presence_absence.Rtab"
    module_file_name = ppanggolin_dir + "/modules/functional_modules.tsv"
    persistent_file_name = ppanggolin_dir + "/partitions/persistent.txt"
    shell_file_name = ppanggolin_dir + "/partitions/shell.txt"
    cloud_file_name = ppanggolin_dir + "/partitions/cloud.txt"

    strains_with_families = read_pangenome_rtab(rtab)

    partitions_with_families = read_partition_files(
        persistent_file_name, shell_file_name, cloud_file_name
    )

    modules_with_families = read_functional_module_file(module_file_name)

    pgdb = pythoncyc.select_organism(pgdb_name)
    meta = pythoncyc.select_organism("meta")

    write_pgdb_reaction_presence_absence_by_strain(
        pgdb, strains_with_families, partitions_with_families, meta, pgdb_name
    )

    write_pgdb_pathway_completion_by_strain(
        pgdb,
        strains_with_families,
        partitions_with_families,
        modules_with_families,
        meta,
        pgdb_name,
        True,
    )
    write_pgdb_pathway_completion_by_strain(
        pgdb,
        strains_with_families,
        partitions_with_families,
        modules_with_families,
        meta,
        pgdb_name,
        False,
    )


if __name__ == "__main__":
    main()
