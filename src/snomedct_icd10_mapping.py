"""
Project of SNOMED-CT and ICD10 Mapping by keywords you interest in
PyMedTermino2 is integrated with Owlready2, and store medical terminologies as OWL ontlogies.
This allows relating medical terms from terminologies with user created concepts.
"""

from owlready2.pymedtermino2.umls import *
from pathlib import Path
from typing import Union, List
import pandas as pd
import collections

PathLike = Union[Path | str]


def init_db(dir_umls_data: PathLike, dir_sqldb: PathLike, terminologies: List[str] = None):
    """
    Build a database from UMLS data.
    please download UMLS data here: https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html
    only need to run this function once.

    :param dir_umls_data: Absolute path (do not use ~/) to the UMLS data (zip format). e.g.,/PATH/YOUR/DIRECTORY/umls-2025AA-metathesaurus-full.zip
    :param dir_sqldb: Absolute path for ontology database with database name e.g.,/PATH/YOUR/DIRECTORY
    :param terminologies: by default, terminologies includes ICD10, SNOMEDCT_US, and CUI
    """
    try:
        # Check if directory exists before using it
        if not dir_sqldb.exists():
            raise FileNotFoundError(f"Directory not found: {dir_sqldb}")

        default_world.set_backend(filename=dir_sqldb / 'med_terminology.sqlite3',
                                  exclusive=False)  # data are stored in path_sqldb which is the quadstore file.
        print(f'SQL database has been created at {dir_sqldb}')

    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        print("Please check your directory path or database file.")

    print('Setting backend...')
    if terminologies is None:
        # CUI correspond to concepts: a given concept gathers equivalent terms or codes from various terminologies.
        terminologies = ["ICD10", "SNOMEDCT_US", "CUI"]

    try:
        # check if directory exists before using it
        if not dir_umls_data.exists():
            raise FileNotFoundError(f"Directory not found:{dir_umls_data}")
        umls_zip = dir_umls_data / 'umls-2025AA-metathesaurus-full.zip'
        if not umls_zip.exists():
            print(
                f"Please see contents on this directory {[x.name for x in dir_umls_data.iterdir() if x.suffix == '.zip']}")

            raise FileNotFoundError(f"UMLS data file not found: {umls_zip}")

        import_umls(dir_umls_data / 'umls-2025AA-metathesaurus-full.zip', terminologies)
        default_world.save()  # for using persistence
        print(f'Medical ontology from UMLS data is successfully saved at {dir_sqldb} as a SQLite Database')
    except (ValueError, FileNotFoundError) as e:
        print(f"Error:{e}")
        print("Please check your directory path or UMLS data file.")


def open_db(dir_sqldb: PathLike):
    """
    Loading the medical ontology from path_sqldb (e.g., db.sqlite3) built from init_db
    :param dir_sqldb: Directory of db.sqlite3
    """
    owlready2 = World()  # build a clean world

    try:
        # Check if directory exists before using it
        if not dir_sqldb.exists():
            raise FileNotFoundError(f"Directory not found: {dir_sqldb}")

        owlready2.set_backend(filename=dir_sqldb / 'med_terminology.sqlite3')  # load .sqlite3 database
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        print("Please check your directory path or database file.")
    med_ontology = owlready2.get_ontology("http://PYM/").load()
    return owlready2, med_ontology


def get_hierarchy(keywords: str, ancestor: bool = False, descendant: bool = True,
                  output_path: PathLike = None) -> pd.DataFrame:
    """
    Using Keywords to search SNOMED-CT concept with all hierarchy
    Both methods ("ancestor_concepts()" and "descendant_concepts()") remove duplicates automatically.
    :param keywords: keywords you interest in (e.g., breast cancer)
    :param ancestor: concept of all ancestor includes itself, otherwise, setting include_self = False
    :param descendant: concept of all descendant includes itself, otherwise, setting include_self = False
    :param output_path: path of output file
    :return: all hierarchy dataframe
    """
    ret = collections.defaultdict(list)
    concepts = list(snomedct.search(keywords))

    if ancestor:
        seen_anc_code = set()  # remove duplicate
        for concept in concepts:
            for parents in concept.ancestor_concepts(include_self=True):
                anc_code = parents.name
                if anc_code in seen_anc_code:
                    continue
                seen_anc_code.add(anc_code)
                ret['SNOMED-CT Code'].append(anc_code)
                ret['Term'].append(parents.label.first())

    if descendant:
        seen_dec_code = set()  # remove duplicate
        for concept in concepts:
            for child in concept.descendant_concepts(include_self=True):
                child_code = child.name
                if child_code in seen_dec_code:
                    continue
                seen_dec_code.add(child_code)
                ret['SNOMED-CT Code'].append(
                    child_code)  # as concept includes itself, so it doesn't need to append concept name or label
                ret['Term'].append(
                    child.label.first())

    all_hierarchy = pd.DataFrame(ret)
    all_hierarchy.sort_values(by='Term', ascending=True, inplace=True)  # sort alphabetically
    # print(all_hierarchy.to_markdown(index=False)) # for README
    if output_path is not None:
        all_hierarchy.to_csv(f'{keywords}_{output_path}', index=False)
    return all_hierarchy


def mapping_icd10(keywords: str, refining_term: str = None, map_child: bool = True,
                  output_path: PathLike = None) -> pd.DataFrame:
    """
    Mapping SNOMED-CT to ICD10 by keywords and/or refining term
    :param keywords: keyword
    :param refining_term: searching for specific term
    :param map_child: looking for terms in all descendant
    :param output_path: output path
    :return: df
    """
    ret = collections.defaultdict(list)
    concepts = list(snomedct.search(keywords))

    if map_child:
        seen_code = set()  # remove duplicate
        for concept in concepts:
            for child in concept.descendant_concepts(include_self=True):
                for icd in child >> icd10:
                    key = (child.name, icd.name)

                    if key in seen_code:
                        continue
                    seen_code.add(icd)
                    ret['SNOMED-CT Code'].append(child.name)
                    ret['SNOMED-CT Term'].append(child.label.first())
                    ret['ICD10 Code'].append(icd.name)
                    ret['ICD10 Term'].append(icd.label.first())

    snomed_icd10_report = pd.DataFrame(ret)
    snomed_icd10_report.sort_values(by='SNOMED-CT Term', ascending=True, inplace=True)  # sort alphabetically

    if refining_term is not None:
        mask = snomed_icd10_report['SNOMED-CT Term'].apply(
            lambda x: refining_term.lower() in str(x).lower()
        )
        refining_report = snomed_icd10_report[mask]
        refining_report.to_csv(f'{refining_term}_{output_path}', index=False)
        # print(refining_report.to_markdown(index=False)) # for README

    if output_path is not None:
        snomed_icd10_report.to_csv(f'{keywords}_{output_path}', index=False)

    return snomed_icd10_report


if __name__ == '__main__':
    dir_umls = Path('/PATH/YOUR/DIRECTORY/')
    dir_sqldb = Path('/PATH/YOUR/DIRECTORY/')
    init_db(dir_umls, dir_sqldb)
    owl, med_ontology = open_db(dir_sqldb)
    snomedct = med_ontology["SNOMEDCT_US"]
    icd10 = med_ontology["ICD10"]
    keywords = 'breast cancer'
    all_hierarchy = get_hierarchy(keywords=keywords, ancestor=False, descendant=True, output_path='snomed.csv')
    mapping_icd10(keywords=keywords, refining_term='hormone receptor positive', output_path='snomed_icd_mapping.csv')
