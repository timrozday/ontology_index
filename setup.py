from setuptools import setup

setup(
   name='ontology_index',
   version='0.6',
   description='Allows parsing and saving of EFO and MeSH ontologies to BSDDB, and SPARQL querying using `rdflib`. \
                For finding related entities within the ontologies (either by measuring distance between multiple entities, or by retrieving all related entities). \
                This involves the construction, storage and querying of indexes.',
   author='Tim Rozday',
   author_email='timrozday@ebi.ac.uk',
   packages=['ontology_index'],  #same as name
)
