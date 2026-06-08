# Prepare reference files

## MetaCyc to KEGG reaction cross-references

`pan2met-wf` requires multiple tabulated files with cross-references between reference databases on metabolism.
Prepare the required KEGG / MetaCyc / EC-numpbers crossref files:

You will need a local installation of PathwayTools with MetaCyc, and a working internet connection to be able to query the KEGG REST API.

In a shell, launch the PathwayTools' python API:

```bash
pathway-tools -lisp -python-local-only-non-strict
```

In another shell, run `extract_kegg_metacyc_data.py`, on a environment with [pythoncyc](https://pypi.org/project/PythonCyc) installed.

```bash
python3 -m venv .venv/pythoncyc
source .venv/pythoncyc/bin/activate
pip install PythonCyc
```

```bash
mkdir KEGG_MetaCyc_crossref
pushd KEGG_MetaCyc_crossref
python3 ../bin/extract_kegg_metacyc_data.py
popd
```

Then, do not forget to adapt the configuration file, you should be good to go updating the `params.KEGG_MetaCyc_crossref` full path, the other crossreference file paths will depend on it.

```nextflow
params {
    KEGG_MetaCyc_crossref = '/path/to/KEGG_MetaCyc_crossref/'
    ecocyc_monomers_to_metacyc_reactions = params.KEGG_MetaCyc_crossref + 'ecocyc_monomers_to_metacyc_reactions.tsv'
    metacyc_monomers_to_reactions = params.KEGG_MetaCyc_crossref + 'metacyc_monomers_to_reactions.tsv'
    kegg_kos_to_metacyc_reactions = params.KEGG_MetaCyc_crossref + 'kegg_kos_to_metacyc_reactions.tsv'
    kegg_kos_to_ec_numbers = params.KEGG_MetaCyc_crossref + 'kegg_kos_to_ec_numbers.tsv'
    
}
```

## KOfam HMM database

The KOFamScan run on the gene families reference proteins require a local installation of the KOfam HMM profiles database.
You can find this database on the KEGG FTP server from <https://www.genome.jp/tools/kofamkoala/>: <https://www.genome.jp/ftp/db/kofam/> you will need both `ko_list.gz` and `profiles.tar.gz` files. 
You will need to unzip these files to a folder, e.g. `KOfam/`.
You will also need to adapt the configuration file to update the parameter `params.kofam_db` to the full path of the `KOfam` folder:

```nextflow
params {
    kofam_db = '/path/to/KOfam/'
}
```

The KOfamScan software enables to associate protein sequence to protein families, which can further be associated to EC-number or KEGG Orthology groups, which can be associated themselves to KEGG reactions and via cross-references to MetaCyc reactions.

## NCBIfam HMM database

NCBIfam (former TIGRfam) is another protein families HMM profile database.

NCBIfam are available on NCBI website at <https://www.ncbi.nlm.nih.gov/refseq/annotation_prok/tigrfams/>. 
One can download `hmm_PGAP.HMM.tgz` file from NCBI FTP server at <https://ftp.ncbi.nlm.nih.gov/hmm/current/>.

You will need to uncompress the `tgz` archive:
```bash
tar xvzf hmm_PGAP.HMM.tgz
```

You will also need to prepare the HMM for `hmmerscan`:

```bash
cat hmm_PGAP.HMM/*.hmm pgap.hmm
hmmpress pgap.hmm
```

