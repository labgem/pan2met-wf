#!/usr/bin/env python

# -*- coding: utf-8 -*-

import time
import re
import io
from collections import defaultdict
import logging

from urllib.request import urlopen

import pythoncyc
from pythoncyc import PTools as PTools
from pythoncyc.PTools import PToolsError as PToolsError
from pythoncyc.PTools import PythonCycError as PythonCycError

logger = logging.getLogger("panmeta:crossref")

def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate

@static_vars(previous=0)
def query_kegg(resource):
    # Extracted from biopython.Bio.KEGG.REST
    # Rate limit to at most 3 query per second.
    delay = 0.333333333  # one third of a second
    current = time.time()
    wait = query_kegg.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        query_kegg.previous = current + wait
    else:
        query_kegg.previous = current

    URL = "https://rest.kegg.jp/%s"
    resp = urlopen(URL % resource)
    handle = io.TextIOWrapper(resp, encoding="UTF-8")
    handle.url = resp.url
    return handle


def get_pgdb_reactions(pgdb):
    return pgdb.all_rxns(type_of_reactions=":all")


def get_kegg_reactions_to_pgdb(pgdb):
    reactions = get_pgdb_reactions(pgdb)
    kegg_reactions_to_pgdb = dict()
    for rxn in reactions:
        if pgdb[rxn].dblinks is not None and "|LIGAND-RXN|" in pgdb[rxn].dblinks:
            for kegg_rxn in pgdb[rxn].dblinks["|LIGAND-RXN|"]:
                if isinstance(kegg_rxn, str) and re.match("^R[0-9]+$", kegg_rxn):
                    if kegg_rxn not in kegg_reactions_to_pgdb:
                        kegg_reactions_to_pgdb[kegg_rxn] = set()
                    kegg_reactions_to_pgdb[kegg_rxn].add(rxn.split("|")[1])
    return kegg_reactions_to_pgdb


def write_monomers_to_reactions(pgdb, outputfile, reactions_tokeep=list()):
    reactions = get_pgdb_reactions(pgdb)
    if len(reactions_tokeep):
        reactions = list(set(reactions) & set(reactions_tokeep))
    monomerswithrxn = dict()
    for rxn in reactions:
        enzymes = pgdb.enzymes_of_reaction(rxn)
        for enz in enzymes:
            monomers = pgdb.monomers_of_protein(enz, unmodify=True)
            for mon in monomers:
                if mon is not None:
                    if isinstance(mon, list):
                        for mon_of_mon in mon:
                          if isinstance(mon_of_mon, str):
                            if mon_of_mon not in monomerswithrxn:
                                monomerswithrxn[mon_of_mon] = set()
                            monomerswithrxn[mon_of_mon].add(rxn.split("|")[1])
                    elif isinstance(mon, str):
                        if mon not in monomerswithrxn:
                            monomerswithrxn[mon] = set()
                        monomerswithrxn[mon].add(rxn.split("|")[1])
                    else:
                        logger.warning(f"{mon} neither list nor str")
    with open(outputfile, "w") as monomers_to_reactions_file:
        for mon in monomerswithrxn:
            monomers_to_reactions_file.write(
                mon.split("|")[1] + "\t" + ",".join(monomerswithrxn[mon]) + "\n"
            )


def write_kegg_kos_to_metacyc_reactions(outputfile, kegg_reactions_to_metacyc):
    reactions = query_kegg("list/reaction").read()
    reactions = [line.split("\t")[0] for line in reactions.split("\n")]
    ko_to_rxn = defaultdict(set)
    rxn_to_ko = [
        line.split('\t')
        for line in query_kegg("link/orthology/rn").read().split('\n')[:-1]
    ]
    for rxn, ko in rxn_to_ko:
        rxn = rxn.split(":")[1]
        ko = ko.split(":")[1]
        ko_to_rxn[ko].add(rxn)

    with open(outputfile, "w") as kegg_kos_to_metacyc_reactions_file:
        for ko in ko_to_rxn:
            metacyc_rxn = set()
            for rxn in ko_to_rxn[ko]:
                if rxn in kegg_reactions_to_metacyc:
                    metacyc_rxn.update(kegg_reactions_to_metacyc[rxn])
            kegg_kos_to_metacyc_reactions_file.write(
                ko
                + "\t"
                + ",".join(ko_to_rxn[ko])
                + "\t"
                + ",".join(metacyc_rxn)
                + "\n"
            )


def write_kegg_kos_to_ec_numbers(outputfile):
    response = query_kegg("link/ec/ko").read()
    ko_to_ec = defaultdict(list)
    for row in response.split("\n"):
        parts = row.split("\t")
        if len(parts) != 2:
            continue
        ko, ec = parts
        ko = ko.replace("ko:", "")
        ec = ec.replace("ec:", "")
        ko_to_ec[ko].append(ec)
    with open(outputfile, "w") as kegg_kos_to_ec_numbers_file:
        for ko, ec_numbers in ko_to_ec.items():
            kegg_kos_to_ec_numbers_file.write(
                "\t".join([ko, ",".join(ec_numbers)]) + "\n"
            )


def main():
    meta = pythoncyc.select_organism("meta")
    metacyc_reactions = get_pgdb_reactions(meta)
    write_monomers_to_reactions(meta, "metacyc_monomers_to_reactions.tsv")
    kegg_reactions_to_metacyc = get_kegg_reactions_to_pgdb(meta)
    eco = pythoncyc.select_organism("eco")
    write_monomers_to_reactions(
        eco, "ecocyc_monomers_to_metacyc_reactions.tsv", metacyc_reactions
    )
    write_kegg_kos_to_metacyc_reactions(
        "kegg_kos_to_metacyc_reactions.tsv", kegg_reactions_to_metacyc
    )
    write_kegg_kos_to_ec_numbers("kegg_kos_to_ec_numbers.tsv")


if __name__ == "__main__":
    main()
