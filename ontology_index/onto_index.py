import rdflib
import itertools as it
from collections import defaultdict
import pickle
import json
from tqdm.auto import tqdm
import zipfile

class EfoIndex():
    equivalent_rels = {
        "http://www.w3.org/2002/07/owl#equivalentClass",
        "http://purl.obolibrary.org/obo/mondo#exactMatch",
    }
    close_rels = {
        "http://purl.obolibrary.org/obo/mondo#narrowMatch",
        "http://purl.obolibrary.org/obo/mondo#broadMatch",
        "http://purl.obolibrary.org/obo/mondo#closeMatch",
    }
    xref_rels = {
        "http://www.geneontology.org/formats/oboInOwl#hasDbXref",
    }
    child_rels = set()
    parent_rels = {
        'http://www.w3.org/2000/01/rdf-schema#subClassOf',
    }
    name_labels = {
        'http://www.w3.org/2000/01/rdf-schema#label',
        'http://www.w3.org/2004/02/skos/core#prefLabel',
        'http://purl.obolibrary.org/obo/ArrayExpress_label',
        'http://www.geneontology.org/formats/oboInOwl#hasExactSynonym',
        'http://www.geneontology.org/formats/oboInOwl#hasRelatedSynonym',
        'http://www.geneontology.org/formats/oboInOwl#hasBroadSynonym',
        'http://www.geneontology.org/formats/oboInOwl#hasNarrowSynonym',
        'http://www.geneontology.org/formats/oboInOwl#shorthand',
        'http://purl.obolibrary.org/obo/IAO_0100001',
        'http://purl.obolibrary.org/obo/SRA_label'
    }
    pref_label = 'http://www.w3.org/2000/01/rdf-schema#label'
    
    name_ranks = {
        'http://www.w3.org/2000/01/rdf-schema#label': 1,
        'http://www.w3.org/2004/02/skos/core#prefLabel': 2,
        'http://www.geneontology.org/formats/oboInOwl#hasExactSynonym': 3,
        'http://www.geneontology.org/formats/oboInOwl#shorthand': 4,
        'http://www.geneontology.org/formats/oboInOwl#hasRelatedSynonym': 5,
        'http://www.geneontology.org/formats/oboInOwl#hasBroadSynonym': 6,
        'http://www.geneontology.org/formats/oboInOwl#hasNarrowSynonym': 7,
        'http://purl.obolibrary.org/obo/SRA_label': 8,
        'http://purl.obolibrary.org/obo/ArrayExpress_label': 9,
        'http://purl.obolibrary.org/obo/IAO_0100001': 10,
    }
    
    disease_root_iris = {
        "http://www.ebi.ac.uk/efo/EFO_0000408",   # EFO disease
        "http://purl.obolibrary.org/obo/EFO_0000408",  # EFO disease
        "http://www.ebi.ac.uk/efo/EFO_0000651",  # EFO phenotype
        "http://purl.obolibrary.org/obo/EFO_0000651",  # EFO phenotype
        "http://purl.obolibrary.org/obo/MONDO_0000001",  # MONDO disease
        "http://purl.obolibrary.org/obo/OMIT_0005457",  # OMIT disease
        "http://purl.obolibrary.org/obo/OMIT_0002893",  # OMIT Mental disorders
        "http://purl.obolibrary.org/obo/DOID_4",  # DOID disease
        "http://purl.obolibrary.org/obo/NCIT_C7057",  # NCIT disease, disorder or finding
        "http://www.orpha.net/ORDO/Orphanet_377788",  # Orphanet disease
        "http://www.orpha.net/ORDO/Orphanet_C001",  # Orphanet clinical entity
        "http://purl.obolibrary.org/obo/OMIT_0001003",  # OMIT diseases category
        "http://purl.obolibrary.org/obo/MONDO_0042489",  # EFO/MONDO disease susceptibility
        "http://purl.obolibrary.org/obo/OMIT_0005461",  # OMIT disease susceptibility
#         "http://purl.obolibrary.org/obo/GO_0008150",  # EFO Biological process
        "http://purl.obolibrary.org/obo/OBI_1110122",  # EFO Pathologic process
        "http://purl.obolibrary.org/obo/NCIT_C16956",  # NCIT Pathologic process
        "http://purl.obolibrary.org/obo/HP_0000118",  # HP Phenotypic abnormality
        "http://purl.obolibrary.org/obo/SYMP_0000462",  # SYMP Symptom
#         "http://purl.obolibrary.org/obo/OMIT_0011629",  # OMIT phenotype
#         "http://purl.obolibrary.org/obo/NCIT_C16977",  # NCIT phenotype
#         "http://purl.obolibrary.org/obo/UPHENO_0001001",  # DOID phenotype
    }
    
    relevant_root_nodes = {
        "http://purl.obolibrary.org/obo/EFO_0000408",
        "http://purl.obolibrary.org/obo/EFO_0000651",
        "http://purl.obolibrary.org/obo/MONDO_0042489",
    }
    
    def __init__(self, data_dir='.'):
        self.data_dir = data_dir
        
        self.rel_dict = {
            **{k:'equivalent' for k in self.equivalent_rels}, 
            **{k:'close' for k in self.close_rels},
            **{k:'xref' for k in self.xref_rels},
            **{k:'child' for k in self.child_rels},
            **{k:'parent' for k in self.parent_rels},
        }
        
        self.rev_rel_dict = {
            **{k:'equivalent' for k in self.equivalent_rels}, 
            **{k:'close' for k in self.close_rels},
            **{k:'xref' for k in self.xref_rels},
            **{k:'parent' for k in self.child_rels},
            **{k:'child' for k in self.parent_rels},
        }
        
        self.efo_graph = rdflib.ConjunctiveGraph(store="Sleepycat")
        
        
        try:
            r = self.efo_graph.open(f"{self.data_dir}/efo.db", create=False)
            assert r == rdflib.store.VALID_STORE, "Invalid EFO store"
        except:
            pass
        
        try:
            self.load_indexes()
        except:
            pass
        
        self.cache = {}

    def is_disease(self, iri):     
        if (iri in self.iri2name) or (iri in self.rels_index) or (iri in self.rev_rels_index):        
            return iri in self.disease_iris
        
    def get_children(self, iri, equivalents=True):
        rels = set()

        if iri in self.rels_index:
            rels.update(self.rels_index[iri])
        if iri in self.rev_rels_index:
            rels.update(self.rev_rels_index[iri])

        allowed_p = {'child'}
        if equivalents:
            allowed_p.add('equivalent')

        return {i for p,i in rels if p in allowed_p}

    def get_descendents(self, iris, covered_iris=None, jumps=1, equivalents=True):
        if isinstance(iris, str):
            iris = {iris}
        iris = set(iris)

        if covered_iris is None:
            covered_iris = iris

        xrefs = set()

        for iri in iris:
            xrefs.update(self.get_children(iri, equivalents=equivalents))  # ontology xrefs

        new_xrefs = xrefs - covered_iris

        if new_xrefs and not (jumps==0 or jumps==1):
            xrefs.update(
                self.get_descendents(
                    new_xrefs, 
                    covered_iris=covered_iris|xrefs, 
                    jumps=jumps-1
                )
            )

        return xrefs
        
    def get_distant_efo_relatives(self, iri, distance=2, distant_rels={'close', 'child', 'parent'}, equivalent_rels={'equivalent'}):

        def rec_f(iri, distance=2, related_iris={}):

            def get_efo_relatives(iri):
                rels = set()

                if self.rels_index:
                    if iri in self.rels_index:
                        rels.update(self.rels_index[iri])
                    if iri in self.rev_rels_index:
                        rels.update(self.rev_rels_index[iri])
                    
                    return rels
                
                else:
                    if not iri in self.cache:
                        p_str = ','.join(f"<{i}>" for i in self.equivalent_rels|self.close_rels|self.xref_rels|self.child_rels|self.parent_rels)
                        
                        query = self.efo_graph.query(f"SELECT ?p ?o WHERE {{ ?q ?p ?o . FILTER ( ?p IN({p_str}) )}}", initBindings={'q': rdflib.URIRef(iri)})
                        rels.update({(self.rel_dict[str(p)],o) for p,o in query if isinstance(o, rdflib.term.URIRef)})

                        query = self.efo_graph.query(f"SELECT ?p ?s WHERE {{ ?s ?p ?q . FILTER ( ?p IN({p_str}) )}}", initBindings={'q': rdflib.URIRef(iri)})
                        rels.update({(self.rev_rel_dict[str(p)],s) for p,s in query if isinstance(s, rdflib.term.URIRef)})
                        
                        self.cache[iri] = rels
                        
                    return self.cache[iri]

            for predicate, related_iri in get_efo_relatives(iri):

                new_d = None
                if predicate in distant_rels:
                    new_d = distance-1
                if predicate in equivalent_rels:
                    new_d = distance

                if (not new_d is None) and (new_d >= 0):
                    if not ((related_iri in related_iris) and (new_d < related_iris[related_iri])):
                        related_iris[related_iri] = new_d
                        related_iris = rec_f(related_iri, distance=new_d, related_iris=related_iris)

            return related_iris

        r = rec_f(iri, distance=distance, related_iris={})
        r = {str(k):distance-d for k,d in r.items() if not str(k) == str(iri)}  # adjust distances

        return r

    def get_efo_links(self, iris, distance=2):
        mappings = {}
        for iri in iris:
            rels = self.get_distant_efo_relatives(rdflib.URIRef(iri), distance=distance)
            mappings[iri] = set(rels.items())

        links = defaultdict(set)
        for s_iri in iris:
            if s_iri in mappings:
                for t_iri,d in mappings[s_iri]:
                    if t_iri in iris:
                        links[tuple(sorted([s_iri,t_iri]))].add(d)

        return {(k1,k2,min(vs)) for (k1,k2),vs in links.items()}
    
    def get_name(self, iri):
#         try:
#             query = self.efo_graph.query(f"SELECT ?o WHERE {{ ?q <http://www.w3.org/2000/01/rdf-schema#label> ?o }}", initBindings={'q': rdflib.URIRef(iri)})
#             return str(list(query)[0][0])
#         except:
#             return None
        
        if iri in self.iri2pref_name and self.iri2pref_name[iri]:
            return self.iri2pref_name[iri]
        
        if iri in self.iri2name and self.iri2name[iri]:
            return sorted([(p,n) for p,n in self.iri2name[iri]], key=lambda x:(self.name_ranks[x[0]], len(x[1])) )[0][1]
        
    def get_names(self, iri):
        return {(n,p,self.name_ranks[p]) for p,n in self.iri2name[iri]}
    
    def get_xrefs(self, iri):
        xrefs = set()
        if iri in self.xref_index:
            xrefs.update(self.xref_index[iri])
        if iri in self.rev_xref_index:
            xrefs.update(self.rev_xref_index[iri])
        return xrefs
    
    def gen_rel_indexes(self):
        p_str = ','.join(f"<{i}>" for i in self.equivalent_rels|self.close_rels|self.child_rels|self.parent_rels)  # ['owl:equivalentClass', ':exactMatch', ':closeMatch', ':narrowMatch', ':broadMatch', 'rdfs:subClassOf', 'oboInOwl:inSubset']
        
        self.rels_index = defaultdict(set)
        self.rev_rels_index = defaultdict(set)
        for p in self.equivalent_rels|self.close_rels|self.child_rels|self.parent_rels:
            p_iri = rdflib.URIRef(p)
            for s,o in tqdm(self.efo_graph.query(f"SELECT ?s ?o WHERE {{ ?s ?p ?o }}", initBindings={'p': p_iri}), leave=True, position=0, desc=str(p)):
                if isinstance(s, rdflib.term.URIRef) and isinstance(o, rdflib.term.URIRef):
                    self.rels_index[str(s)].add((self.rel_dict[str(p)],str(o)))
                    self.rev_rels_index[str(o)].add((self.rev_rel_dict[str(p)],str(s)))
    
    def gen_xref_indexes(self):
        def efo_norm_xref(iri, \
                          prefix_source_map = {'MESH': 'mesh',
                                               'MSH': 'mesh',
                                               'MeSH': 'mesh',
                                               'SCTID': 'snomed',
                                               'SCTID_2010_1_31': 'snomed',
                                               'SNOMEDCT': 'snomed',
                                               'SNOMEDCT_2010_1_31': 'snomed',
                                               'SNOMEDCT_US': 'snomed',
                                               'SNOMEDCT_US_2018_03_01': 'snomed',
                                               'UMLS': 'umls',
                                               'UMLS CUI': 'umls',
                                               'UMLS_CUI': 'umls'
                                              },\
                          source_ns_map = {'snomed': 'snomed:',
                                           'mesh': 'http://id.nlm.nih.gov/mesh/2021/',
                                           'umls': 'UMLS:'
                                          }):
            prefix = iri.split(':')[0]
            code = ':'.join(iri.split(':')[1:])
            if prefix in prefix_source_map:
                source = prefix_source_map[prefix]
                ns = source_ns_map[source]
                return f"{ns}{code}"
        
        self.xref_index = defaultdict(set)
        self.rev_xref_index = defaultdict(set)
        
        p = "http://www.geneontology.org/formats/oboInOwl#hasDbXref"
        p_iri = rdflib.URIRef(p)
        for s,o in tqdm(self.efo_graph.query(f"SELECT ?s ?o WHERE {{ ?s ?p ?o }}", initBindings={'p': p_iri}), leave=True, position=0):
            if isinstance(s, rdflib.term.URIRef):
                o = efo_norm_xref(o)
                self.xref_index[str(s)].add((self.rel_dict[str(p)],str(o)))
                self.rev_xref_index[str(o)].add((self.rev_rel_dict[str(p)],str(s)))
    
        self.xref_index = dict(self.xref_index)
        self.rev_xref_index = dict(self.rev_xref_index)
    
    def gen_disease_indexes(self):
        self.disease_iris = self.get_descendents(self.disease_root_iris, covered_iris=None, jumps=-1, equivalents=True)
    
    def gen_name_indexes(self):
        self.iri2name = defaultdict(set)
        self.iri2pref_name = {}
        for p in self.name_labels:
            p_iri = rdflib.URIRef(p)
            for s,o in tqdm(self.efo_graph.query(f"SELECT ?s ?o WHERE {{ ?s ?p ?o }}", initBindings={'p': p_iri}), leave=True, position=0, desc=str(p)):
                if isinstance(s, rdflib.term.URIRef):
                    self.iri2name[str(s)].add((str(p), str(o)))
                    if p == self.pref_label:
                        self.iri2pref_name[str(s)] = str(o)
                        
        self.iri2name = dict(self.iri2name)
    
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        
        with open(f"{data_dir}/efo_disease_iris.json", 'wt') as f:
            json.dump(list(self.disease_iris), f)
        with open(f"{data_dir}/efo_rels_index.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.rels_index.items()}, f)
        with open(f"{data_dir}/efo_rev_rels_index.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.rev_rels_index.items()}, f)
        with open(f"{data_dir}/efo_xref_index.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.xref_index.items()}, f)
        with open(f"{data_dir}/efo_rev_xref_index.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.rev_xref_index.items()}, f)
        with open(f"{data_dir}/efo_iri2name.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.iri2name.items()}, f)
        with open(f"{data_dir}/efo_iri2pref_name.json", 'wt') as f:
            json.dump(self.iri2pref_name, f)
        
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        
        with open(f"{data_dir}/efo_disease_iris.json", 'rt') as f:
            self.disease_iris = set(json.load(f))
        with open(f"{data_dir}/efo_rels_index.json", 'rt') as f:
            self.rels_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/efo_rev_rels_index.json", 'rt') as f:
            self.rev_rels_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/efo_xref_index.json", 'rt') as f:
            self.xref_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/efo_rev_xref_index.json", 'rt') as f:
            self.rev_xref_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/efo_iri2name.json", 'rt') as f:
            self.iri2name = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/efo_iri2pref_name.json", 'rt') as f:
            self.iri2pref_name = json.load(f)
    

class MeshIndex():
    term_rels = {
        'http://id.nlm.nih.gov/mesh/vocab#term',
        'http://id.nlm.nih.gov/mesh/vocab#preferredTerm',
#         'http://id.nlm.nih.gov/mesh/vocab#useInstead', 
#         'http://id.nlm.nih.gov/mesh/vocab#mappedTo', 
#         'http://id.nlm.nih.gov/mesh/vocab#preferredMappedTo', 
#         'http://id.nlm.nih.gov/mesh/vocab#preferredConcept', 
#         'http://id.nlm.nih.gov/mesh/vocab#concept'
    }  
    concept_rels = {
        'http://id.nlm.nih.gov/mesh/vocab#preferredConcept',
        'http://id.nlm.nih.gov/mesh/vocab#concept',
    }
    child_rels = set()  # 'mesh_vocab:narrowerConcept'
    parent_rels = {
        'http://id.nlm.nih.gov/mesh/vocab#parentTreeNumber',
#         'mesh_vocab:broaderConcept', 
#         'mesh_vocab:broaderDescriptor', 
#         'mesh_vocab:hasDescriptor'
    }  
    name_labels = {
        'http://www.w3.org/2000/01/rdf-schema#label',
        'http://id.nlm.nih.gov/mesh/vocab#altLabel',
        'http://id.nlm.nih.gov/mesh/vocab#prefLabel'
    }
    pref_label = 'http://id.nlm.nih.gov/mesh/vocab#prefLabel'
    
    
    name_ranks = {
        'http://id.nlm.nih.gov/mesh/vocab#prefLabel': 1,
        'http://www.w3.org/2000/01/rdf-schema#label': 2,
        'http://id.nlm.nih.gov/mesh/vocab#altLabel': 3,
    }
    
    relevant_root_treenumbers = {
        "C01", "C02", "C03", "C04", "C05", 
        "C06", "C07", "C08", "C09", "C10", 
        "C11", "C12", "C13", "C14", "C15", 
        "C16", "C17", "C18", "C19", "C20", 
        "C21", "C22", "C23", "C24", "C25", 
        "C26", "F03", 
    }
    
    def __init__(self, data_dir='.'):
        self.data_dir = data_dir
        
        self.mesh_graph = rdflib.ConjunctiveGraph(store="Sleepycat")
        
        try:
            rt = self.mesh_graph.open(f"{data_dir}/mesh.db", create=False)
            assert rt == rdflib.store.VALID_STORE, "Invalid MeSH store"

            self.mesh_graph.bind('mesh2021', rdflib.URIRef("http://id.nlm.nih.gov/mesh/2021/"))
            self.mesh_graph.bind('vocab', rdflib.URIRef("http://id.nlm.nih.gov/mesh/vocab#"))
        except:
            pass
        
        try:
            self.load_indexes()
        except:
            pass
        
        self.cache = {}
    
    def is_disease(self, iri):
        for tn in self.get_treenumber(iri):
            if tn.split('.')[0] in self.relevant_root_treenumbers:
                return True
        return False
    
    def get_mesh_treenumbers(self, mesh_descriptor_id):
        if not mesh_descriptor_id in self.cache:
            query = self.mesh_graph.query(f"SELECT ?o WHERE {{ mesh2021:{mesh_descriptor_id} vocab:treeNumber ?o }}")
            self.cache[mesh_descriptor_id] = {o.split('/')[-1] for o, in query}
        return self.cache[mesh_descriptor_id]

#     def get_mesh_links(self, iris, distance=2):

#         # get treenumbers
#         iri_treenumbers = {}
#         for iri in iris:
#             if self.iri2treenumber:
#                 r = self.iri2treenumber[iri]
#             else:
#                 r = self.get_mesh_treenumbers(iri)
#             iri_treenumbers[iri] = r

#         treenumber_iris = defaultdict(set)
#         for iri, treenumbers in iri_treenumbers.items():
#             for t in treenumbers:
#                 treenumber_iris[t].add(iri)

#         # parse treenumbers and create distances
#         treenumbers = list(treenumber_iris.keys())
#         treenumber_distances = defaultdict(set)
#         for i,t in enumerate(treenumbers):
#             treenumber_distances[t].add((i,0))
#             t = t.split('.')
#             for d in range(1,distance+1):
#                 trimmed_t = '.'.join(t[:-d])
#                 treenumber_distances[trimmed_t].add((i,d))

#         # find links between treenumbers
#         links = defaultdict(set)
#         for t,ds in treenumber_distances.items():
#             for (i1,d1),(i2,d2) in it.combinations(ds,2):
#                 if i1==i2:
#                     continue
#                 links[tuple(sorted([i1,i2]))].add(d1+d2)

#         # convert to iri links
#         iri_links = set()
#         for (i1,i2),ds in links.items():
#             d = min(ds)
#             if d <= distance:
#                 for iri1,iri2 in it.product(treenumber_iris[treenumbers[i1]], treenumber_iris[treenumbers[i2]]):
#                     if iri1 == iri2:
#                         continue
#                     iri_links.add((iri1, iri2, d))
#         return iri_links

    def get_descendents(self, iri, distance=None):
        xrefs = set()
        for tn in self.get_treenumber(iri):
            if (distance is None) or (distance==-1):
                xrefs.update({f"http://id.nlm.nih.gov/mesh/2021/{i}" for d,i in self.treenumber_index[tn]})
            else:
                xrefs.update({f"http://id.nlm.nih.gov/mesh/2021/{i}" for d,i in self.treenumber_index[tn] if d<=distance})
        
        return xrefs
    
    def get_distant_mesh_relatives(self, iri, distance=2, search_up=True, search_down=True):

        def rec_f(tn, distance=2, related_iris=set()):
            if tn in self.treenumber_index:
                for d, related_iri in self.treenumber_index[tn]:
                    if (d <= distance) and (not related_iri in related_iris):
                        if search_down or d==0:
                            related_iris.add(related_iri)
                            for related_tn in self.iri2treenumber[related_iri]:
                                related_iris = rec_f(related_tn, distance=distance-d, related_iris=related_iris)

            return related_iris

        related_iris = set()
        for original_tn in self.get_treenumber(iri):
            tn_split = original_tn.split('.')

            for i,idx in enumerate(reversed(range(len(tn_split)))):
                if i > distance:
                    break
                sub_tn = '.'.join(tn_split[:idx+1])
                related_iris.update(rec_f(sub_tn, distance=distance-i, related_iris=related_iris))
                
                # if config dictates, terminate the upwards search before it begins
                if not search_up:
                    break

        return related_iris - {iri}

    def get_iri(self, iri):
        if iri in self.concept2iri:
            return self.concept2iri[iri]
        if iri in self.term2iri:
            return self.term2iri[iri]
        return iri
    
    def get_treenumber(self, iri):
        iri = self.get_iri(iri)
        return self.iri2treenumber[iri.split('/')[-1]]
    
    def get_name(self, iri):
#         try:
#             query = self.mesh_graph.query(f"SELECT ?o WHERE {{ ?q <http://www.w3.org/2000/01/rdf-schema#label> ?o }}", initBindings={'q': rdflib.URIRef(iri)})
#             return str(list(query)[0][0])
#         except:
#             return None
        
        if iri in self.iri2pref_name and self.iri2pref_name[iri]:
            return self.iri2pref_name[iri]
        
        if iri in self.iri2name and self.iri2name[iri]:
            return sorted([(p,n) for p,n in self.iri2name[iri]], key=lambda x:(self.name_ranks[x[0]], len(x[1])) )[0][1]

    def get_names(self, iri, all_terms=True):
        iri = self.get_iri(iri)
        
        names = {(n,p,self.name_ranks[p]) for p,n in self.iri2name[iri]}
        
        if all_terms:
            if iri in self.iri2term:
                for p,o in self.iri2term[iri]:
                    names.update((n,p,self.name_ranks[p]) for p,n in self.iri2name[o])
        
        return names
        
    def gen_treenumber_indexes(self):
        self.treenumber_index = defaultdict(set)
        self.iri2treenumber = defaultdict(set)
        for iri,tn in tqdm(self.mesh_graph.query(f"SELECT ?s ?o WHERE {{ ?s vocab:treeNumber ?o }}"), leave=True, position=0):
            iri = str(iri).split('/')[-1]
            tn = tn.split('/')[-1]
            self.iri2treenumber[iri].add(tn)
            tn_split = tn.split('.')
            for i,idx in enumerate(reversed(range(len(tn_split)))):
                self.treenumber_index['.'.join(tn_split[:idx+1])].add((i,iri))
                
        self.treenumber_index = dict(self.treenumber_index)
        self.iri2treenumber = dict(self.iri2treenumber)
    
    def gen_name_indexes(self):
        self.iri2name = defaultdict(set)
        self.iri2pref_name = {}
        for p in self.name_labels:
            p_iri = rdflib.URIRef(p)
            for s,o in tqdm(self.mesh_graph.query(f"SELECT ?s ?o WHERE {{ ?s ?p ?o }}", initBindings={'p': p_iri}), leave=True, position=0, desc=str(p)):
                if isinstance(s, rdflib.term.URIRef):
                    self.iri2name[str(s)].add((str(p), str(o)))
                    if p == self.pref_label:
                        self.iri2pref_name[str(s)] = str(o)
        
        self.iri2name = dict(self.iri2name)
        
    def gen_term_indexes(self):
        def convert_concept(iri):
            if iri in self.concept2iri:
                return self.concept2iri[iri]
            else:
                return iri
        
        self.iri2term = defaultdict(set)
        self.term2iri = {}
        for p in self.term_rels:
            p_iri = rdflib.URIRef(p)
            for s,o in tqdm(self.mesh_graph.query(f"SELECT ?s ?o WHERE {{ ?s ?p ?o }}", initBindings={'p': p_iri}), leave=True, position=0, desc=str(p)):
                if isinstance(s, rdflib.term.URIRef) and isinstance(o, rdflib.term.URIRef):
                    self.iri2term[str(s)].add((str(p), str(o)))
                    self.term2iri[str(o)] = convert_concept(str(s))
        
        self.iri2term = dict(self.iri2term)
        self.term2iri = dict(self.term2iri)
    
    def gen_concept_indexes(self):
        self.iri2concept = defaultdict(set)
        self.concept2iri = {}
        for p in self.concept_rels:
            p_iri = rdflib.URIRef(p)
            for s,o in tqdm(self.mesh_graph.query(f"SELECT ?s ?o WHERE {{ ?s ?p ?o }}", initBindings={'p': p_iri}), leave=True, position=0, desc=str(p)):
                if isinstance(s, rdflib.term.URIRef) and isinstance(o, rdflib.term.URIRef):
                    self.iri2concept[str(s)].add((str(p), str(o)))
                    self.concept2iri[str(o)] = str(s)
        
        self.iri2concept = dict(self.iri2concept)
        self.concept2iri = dict(self.concept2iri)
        
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/treenumber_index.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.treenumber_index.items()}, f)
        with open(f"{data_dir}/iri2treenumber.json", 'wt') as f:
            json.dump({k:list(vs) for k,vs in self.iri2treenumber.items()}, f)
        with open(f"{data_dir}/mesh_iri2name.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.iri2name.items()}, f)
        with open(f"{data_dir}/mesh_iri2pref_name.json", 'wt') as f:
            json.dump(self.iri2pref_name, f)
        
        with open(f"{data_dir}/mesh_iri2term.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.iri2term.items()}, f)
        with open(f"{data_dir}/mesh_term2iri.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.term2iri.items()}, f)
        with open(f"{data_dir}/mesh_iri2concept.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.iri2concept.items()}, f)
        with open(f"{data_dir}/mesh_concept2iri.json", 'wt') as f:
            json.dump(self.concept2iri, f)
        
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/treenumber_index.json", 'rt') as f:
            self.treenumber_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/iri2treenumber.json", 'rt') as f:
            self.iri2treenumber = {k:set(vs) for k,vs in json.load(f).items()}
        with open(f"{data_dir}/mesh_iri2name.json", 'rt') as f:
            self.iri2name = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/mesh_iri2pref_name.json", 'rt') as f:
            self.iri2pref_name = json.load(f)
            
        with open(f"{data_dir}/mesh_iri2term.json", 'rt') as f:
            self.iri2term = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/mesh_term2iri.json", 'rt') as f:
            self.term2iri = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/mesh_iri2concept.json", 'rt') as f:
            self.iri2concept = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/mesh_concept2iri.json", 'rt') as f:
            self.concept2iri = json.load(f)
        
class UmlsIndex():
    name = "umls"
    pref_label = 'umls:cui_pref_string'
    name_labels = {
        'umls:cui_string', 
        'umls:cui_pref_string'
    }
    relevant_root_nodes = set()
    child_rels = set()
    parent_rels = set()
    equivalent_rels = {
        'umls:same_cui'
    }
    xref_rels = set()
    
    pred_types = {
        'umls:same_cui': {'equivalent'},
        'umls:cui_string': {'name'},
        'umls:cui_pref_string': {'name', 'pref_label'}
    }
    
    name_ranks = {
        'PF': 1,
        'VC': 2,
        'VW': 3,
        'VCW': 4,
        'VO': 5
    }
    
    name_types = {
        'PF': 'umls:pref_term',
        'VCW': 'umls:case_word_order_variant',
        'VC': 'umls:case_variant',
        'VO': 'umls:variant',
        'VW': 'umls:word_order_variant',
    }
    
    good_semantic_types = {
        'T019',  # Congenital Abnormality
        'T020',  # Acquired Abnormality
        'T190',  # Anatomical Abnormality
        'T033',  # Finding 
        'T037',  # Injury or Poisoning
        'T041',  # Mental Process
        'T046',  # Pathologic Function
        'T047',  # Disease or Syndrome
        'T048',  # Mental or Behavioral Dysfunction
        'T049',  # Cell or Molecular Dysfunction
        'T050',  # Experimental Model of Disease
        'T184',  # Sign or Symptom
        'T191',  # Neoplastic Process
        'T201',  # Clinical Attribute
    }
    
    def __init__(self, filepath=None, data_dir='.'):
        self.data_dir = data_dir
        self.filepath = filepath
        
        try:
            self.load_indexes()
        except:
            pass
    
    def is_disease(self, iri):
        semantic_types = self.iri2semantic_types[iri]
        return bool(semantic_types & self.good_semantic_types)
    
    def gen_terms_and_rel_indexes(self, filepath=None):
        if filepath is None:
            filepath = self.filepath
            
        def gen_iri(source, code, source_name_map={'MSH': 'http://id.nlm.nih.gov/mesh/2021/', 'SNOMEDCT_US': 'snomed:'}):
            if source in source_name_map:
                prefix = source_name_map[source]
                return f"{prefix}{code}"
        
        def clean_line(l):
            line = l.decode('utf-8')
            while True:
                if line[-1] in {'\n', '\r'}: 
                    line = line[:-1]
                else:
                    break
            return line

#         cui_terms = defaultdict(list)
        sources = set()
        equivalent_entities = defaultdict(set)
        cui_terms = defaultdict(list)
        iri2semantic_types = defaultdict(set)
        
#         with tarfile.open(self.filepath, 'r:') as tf:
#             for member in tf.getmembers():
#                 if re.match('^.*?META/MRCONSO\.RRF$', member.name):
#                     for line in (clean_line(l) for l in tf.extractfile(member)):

        with zipfile.ZipFile(filepath) as f:
            with f.open(name='umls-2020AB-data/MRCONSO.RRF') as df:
                cols = ['CUI','LAT','TS','LUI','STT','SUI','ISPREF','AUI','SAUI','SCUI','SDUI','SAB','TTY','CODE','STR','SRL','SUPPRESS','CVF']
                for line in tqdm((clean_line(l) for l in df), leave=True, position=0, desc='Extracting file'):
                    row_dict = {cols[i]:v for i,v in enumerate(line.split('|')) if i < len(cols)}
                    cui = f"UMLS:{row_dict['CUI']}"
                    if all([
                        str(row_dict['LAT']) == 'ENG',  # only english
                        str(row_dict['SUPPRESS']) == 'N',
                    ]):
#                             sources.add(row_dict['SAB'])  # keep track of what sources are present
                        cui_terms[cui].append({'string': row_dict['STR'], \
                                               'source': row_dict['SAB'], \
                                               'string_type': row_dict['STT'], \
                                               'is_pref': row_dict['ISPREF'], \
                                               'term_status': row_dict['TS'], \
                                               'term_type_in_source': row_dict['TTY']})

                        iri = gen_iri(row_dict['SAB'], row_dict['CODE'])
                        if iri:
                            equivalent_entities[cui].add(iri)
            
            with f.open(name='umls-2020AB-data/MRSTY.RRF') as df:
                cols = ['CUI','STY','?1','?2','?3','?4',]
                for line in tqdm((clean_line(l) for l in df), leave=True, position=0, desc='Extracting file'):
                    row_dict = {cols[i]:v for i,v in enumerate(line.split('|')[:-1])}
                    cui = f"UMLS:{row_dict['CUI']}"
                    iri2semantic_types[cui].add(row_dict['STY'])
        
        self.iri2semantic_types = {k:set(vs) for k,vs in iri2semantic_types.items()}

        self.entity_rels = defaultdict(set)
        self.iri2name = defaultdict(set)
        self.iri2pref_name = {}
        for cui,iris in tqdm(equivalent_entities.items(), leave=True, position=0, desc='Processing rels'):
            for iri1, iri2 in it.permutations(iris|{cui},2):
                self.entity_rels[iri1].add(('umls:same_cui', iri2))
        
        for cui,vs in tqdm(cui_terms.items(), leave=True, position=0, desc='Processing names'):
            for v in vs:
                if v['is_pref']=='Y':
                    self.iri2name[cui].add(('umls:cui_pref_string', v['string_type'], v['string']))
                    self.iri2pref_name[cui] = v['string']
                else:
                    self.iri2name[cui].add(('umls:cui_string', v['string_type'], v['string']))
        
        self.entity_rels = dict(self.entity_rels)
        self.iri2name = dict(self.iri2name)
    
    def get_name(self, iri):
        if iri in self.iri2pref_name and self.iri2pref_name[iri]:
            return self.iri2pref_name[iri]
        
        if iri in self.iri2name and self.iri2name[iri]:
            return sorted([(p,n) for _,p,n in self.iri2name[iri]], key=lambda x:(self.name_ranks[x[0]], len(x[1])) )[0][1]
    
    def get_names(self, iri):
        return {(n,p,self.name_ranks[p]) for _,p,n in self.iri2name[iri]}
    
    def get_xrefs(self, iri):
        if iri in self.entity_rels:
            return self.entity_rels[iri]
    
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        
        with open(f"{data_dir}/umls_iri2semantic_types.json", 'wt') as f:
            json.dump({k:list(vs) for k,vs in self.iri2semantic_types.items()}, f)
        with open(f"{data_dir}/umls_entity_rels.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.entity_rels.items()}, f)
        with open(f"{data_dir}/umls_iri2name.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.iri2name.items()}, f)
        with open(f"{data_dir}/umls_iri2pref_name.json", 'wt') as f:
            json.dump(self.iri2pref_name, f)
        
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/umls_iri2semantic_types.json", 'rt') as f:
            self.iri2semantic_types = {k:set(vs) for k,vs in json.load(f).items()}
        with open(f"{data_dir}/umls_entity_rels.json", 'rt') as f:
            self.entity_rels = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/umls_iri2name.json", 'rt') as f:
            self.iri2name = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/umls_iri2pref_name.json", 'rt') as f:
            self.iri2pref_name = json.load(f)
            
