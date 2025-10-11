"""
Project of SNOMED-CT and ICD10 Mapping

PyMedTermino - Medical Terminologies, including SNOMED CT, ICD10, MedDRA, CUI for Python
PyMedTermino2 is integrated with Owlready2, and store medical terminologies as OWL ontlogies. This allows relating medical terms from terminologies with user created concepts.
"""

from owlready2.pymedtermino2.umls import *
from pathlib import Path
from typing import Union, List
import pandas as pd
import collections

PathLike = Union[Path | str]


def init_db(path_umls_data: PathLike, path_sqldb: PathLike, terminologies: List[str] = None):
    """
    Build a database from UMLS data.
    please download UMLS data here: https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html
    only need to run this function once.

    :param path_umls_data: Absolute path (do not use ~/) to the UMLS data (zip format) dir path of UMLS data. e.g.,/PATH/YOUR/DIRECTORY/umls-2025AA-metathesaurus-full.zip
    :param path_sqldb: Absolute path for ontology database with database name e.g.,/PATH/YOUR/DIRECTORY/db.sqlite3
    :param terminologies: by default, terminologies includes ICD10, SNOMEDCT_US, and CUI
    """

    default_world.set_backend(filename=path_sqldb,
                              exclusive=False)  # data are stored in path_sqldb which is the quadstore file.
    print('Setting backend...')
    if terminologies is None:
        # CUI correspond to concepts: a given concept gathers equivalent terms or codes from various terminologies.
        terminologies = ["ICD10", "SNOMEDCT_US", "CUI"]
    print(f'SQL database has been created at {path_sqldb}')
    import_umls(path_umls_data, terminologies)
    default_world.save()  # for using persistence
    print(f'Medical ontology from UMLS data is successfully saved at {path_sqldb} as a SQLite Database')


def open_db(path_sqldb: PathLike):
    """
    Loading the medical ontology from path_sqldb (e.g., db.sqlite3) built from init_db
    :param path_sqldb: path of db.sqlite3
    """
    owlready2 = World()  # build a clean world
    owlready2.set_backend(filename=path_sqldb)  # load .sqlite3 database
    med_ontology = owlready2.get_ontology("http://PYM/").load()
    return owlready2, med_ontology


def get_hierarchy(keywords: str, ancestor: bool = False, descendant: bool = True,
                  output_path: PathLike = None) -> pd.DataFrame:
    # todo: create hierarchy plot for visualisation
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
    if output_path is not None:
        all_hierarchy.to_csv(f'{keywords}_{output_path}', index=False)
    return all_hierarchy


def mapping_icd10(keywords: str, map_child: bool = True, output_path: PathLike = None):
    ret = collections.defaultdict(list)
    concepts = list(snomedct.search(keywords))

    if map_child:
        seen_icd_code = set()  # remove duplicate
        for concept in concepts:
            for child in concept.descendant_concepts(include_self=True):
                for code_mapping in child >> icd10:
                    if code_mapping in seen_icd_code:
                        continue
                    seen_icd_code.add(code_mapping)

                    ret['ICD10 Code'].append(code_mapping.name)
                    ret['ICD10 Term'].append(code_mapping.label.first())
    snomed_icd10_report = pd.DataFrame(ret)

    if output_path is not None:
        snomed_icd10_report.to_csv(f'{keywords}_{output_path}', index=False)

    return snomed_icd10_report


if __name__ == '__main__':
    umls_path = '/PATH/YOUR/DIRECTORY/umls-2025AA-metathesaurus-full.zip'
    path_sqldb = '/PATH/YOUR/DIRECTORY/med_terminology.sqlite3'
    init_db(umls_path, path_sqldb)
    owl, med_ontology = open_db(path_sqldb)
    snomedct = med_ontology["SNOMEDCT_US"]
    icd10 = med_ontology["ICD10"]
    keywords = 'breast cancer'
    all_hierarchy = get_hierarchy(keywords=keywords, ancestor=False, descendant=True,
                                  output_path='snomed.csv')
    mapping_icd10(keywords=keywords, output_path='snomed_icd_mapping.csv')
