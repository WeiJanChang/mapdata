"""
Project of SNOMED-CT and ICD10 Mapping

PyMedTermino - Medical Terminologies, including SNOMED CT, ICD10, MedDRA, CUI for Python
PyMedTermino2 is integrated with Owlready2, and store medical terminologies as OWL ontlogies. This allows relating medical terms from terminologies with user created concepts.
"""

from owlready2.pymedtermino2.umls import *
from pathlib import Path
from typing import Union, List

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


if __name__ == '__main__':
    umls_path = '/PATH/YOUR/DIRECTORY/umls-2025AA-metathesaurus-full.zip'
    path_sqldb = '/PATH/YOUR/DIRECTORY/med_terminology.sqlite3'
    init_db(umls_path, path_sqldb)
    owl, med_ontology = open_db(path_sqldb)
    snomedct = med_ontology["SNOMEDCT_US"]
    icd10 = med_ontology["ICD10"]
