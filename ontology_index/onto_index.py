import rdflib
import itertools as it
from collections import defaultdict
import pickle
import json
from tqdm.auto import tqdm

class EfoIndex():
    equivalent_rels = {
        "http://www.w3.org/2002/07/owl#equivalentClass",
        "http://purl.obolibrary.org/obo/http://www.ebi.ac.uk/efo/efo.owl#exactMatch",
        "http://purl.obolibrary.org/obo/http://www.ebi.ac.uk/efo/efo.owl#closeMatch",
    }
    close_rels = {
        "http://purl.obolibrary.org/obo/http://www.ebi.ac.uk/efo/efo.owl#narrowMatch",
        "http://purl.obolibrary.org/obo/http://www.ebi.ac.uk/efo/efo.owl#broadMatch",
    }
    xref_rels = {
        "http://www.geneontology.org/formats/oboInOwl#hasDbXref",
    }
    child_rels = set()
    parent_rels = {
        'http://www.w3.org/2000/01/rdf-schema#subClassOf',
        'http://purl.obolibrary.org/obo/http://www.ebi.ac.uk/efo/efo.owl#excluded_subClassOf',
        'http://www.geneontology.org/formats/oboInOwl#inSubset'
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

    def get_distant_efo_relatives(self, iri, distance=2):

        def rec_f(iri, distance=2, related_iris={}):

            def get_efo_relatives(iri):
                p_str = ','.join(f"<{i}>" for i in equivalent_rels|close_rels|xref_rels|child_rels|parent_rels)  # ['owl:equivalentClass', ':exactMatch', ':closeMatch', ':narrowMatch', ':broadMatch', 'rdfs:subClassOf', 'oboInOwl:inSubset']
                rels = set()

                if self.rels_index:
                    rels.update(self.rels_index[iri])
                    rels.update(self.rev_rels_index[iri])
                    
                    return rels
                
                else:
                    if not iri in self.cache:
                        query = self.efo_graph.query(f"SELECT ?p ?o WHERE {{ ?q ?p ?o . FILTER ( ?p IN({p_str}) )}}", initBindings={'q': rdflib.URIRef(iri)})
                        rels.update({(self.rel_dict[str(p)],o) for p,o in query if isinstance(o, rdflib.term.URIRef)})

                        query = self.efo_graph.query(f"SELECT ?p ?s WHERE {{ ?s ?p ?q . FILTER ( ?p IN({p_str}) )}}", initBindings={'q': rdflib.URIRef(iri)})
                        rels.update({(self.rel_dict[str(p)],s) for p,s in query if isinstance(s, rdflib.term.URIRef)})
                        
                        self.cache[iri] = rels
                        
                    return self.cache[iri]

            for predicate, related_iri in get_efo_relatives(iri):

                new_d = None
                if predicate in {'close', 'child', 'parent'}:
                    new_d = distance-1
                if predicate in {'equivalent'}:
                    new_d = distance

                if (not new_d is None) and (new_d >= 0):
                    if not ((related_iri in related_iris) and (new_d < related_iris[related_iri])):
                        related_iris[related_iri] = new_d
                        related_iris = rec_f(related_iri, distance=new_d, related_iris=related_iris)

            return related_iris

        r = rec_f(iri, distance=distance, related_iris={})
        r = {str(k):distance-d for k,d in r.items() if not str(k) == str(iri)}

        return r  # adjust distances

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

    def get_efo_name(self, iri):
        query = self.efo_graph.query(f"SELECT ?o WHERE {{ ?q <http://www.w3.org/2000/01/rdf-schema#label> ?o }}", initBindings={'q': rdflib.URIRef(iri)})
        try:
            return str(list(query)[0][0])
        except IndexError as e:
            return None
        
    def gen_rel_indexes(self):
        p_str = ','.join(f"<{i}>" for i in self.equivalent_rels|self.close_rels|self.xref_rels|self.child_rels|self.parent_rels)  # ['owl:equivalentClass', ':exactMatch', ':closeMatch', ':narrowMatch', ':broadMatch', 'rdfs:subClassOf', 'oboInOwl:inSubset']
        
        self.rels_index = defaultdict(set)
        self.rev_rels_index = defaultdict(set)
        for s,p,o in tqdm(self.efo_graph.query(f"SELECT ?s ?p ?o WHERE {{ ?s ?p ?o . FILTER ( ?p IN({p_str}) )}}"), leave=True, position=0):
            if isinstance(s, rdflib.term.URIRef) and isinstance(o, rdflib.term.URIRef):
                self.rels_index[str(s)].add((self.rel_dict[str(p)],str(o)))
                self.rev_rels_index[str(o)].add((self.rel_dict[str(p)],str(s)))
    
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/efo_rels_index.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.rels_index.items()}, f)
        with open(f"{data_dir}/efo_rev_rels_index.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.rev_rels_index.items()}, f)
        
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/efo_rels_index.json", 'rt') as f:
            self.rels_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/efo_rev_rels_index.json", 'rt') as f:
            self.rev_rels_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
    

class MeshIndex():
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
    
    def get_mesh_treenumbers(self, mesh_descriptor_id):
        if not mesh_descriptor_id in self.cache:
            query = self.mesh_graph.query(f"SELECT ?o WHERE {{ mesh2021:{mesh_descriptor_id} vocab:treeNumber ?o }}")
            self.cache[mesh_descriptor_id] = {o.split('/')[-1] for o, in query}
        return self.cache[mesh_descriptor_id]

    def get_mesh_links(self, iris, distance=2):

        # get treenumbers
        iri_treenumbers = {}
        for iri in iris:
            if self.iri2treenumber:
                r = self.iri2treenumber[iri]
            else:
                r = self.get_mesh_treenumbers(iri)
            iri_treenumbers[iri] = r

        treenumber_iris = defaultdict(set)
        for iri, treenumbers in iri_treenumbers.items():
            for t in treenumbers:
                treenumber_iris[t].add(iri)

        # parse treenumbers and create distances
        treenumbers = list(treenumber_iris.keys())
        treenumber_distances = defaultdict(set)
        for i,t in enumerate(treenumbers):
            treenumber_distances[t].add((i,0))
            t = t.split('.')
            for d in range(1,distance+1):
                trimmed_t = '.'.join(t[:-d])
                treenumber_distances[trimmed_t].add((i,d))

        # find links between treenumbers
        links = defaultdict(set)
        for t,ds in treenumber_distances.items():
            for (i1,d1),(i2,d2) in it.combinations(ds,2):
                if i1==i2:
                    continue
                links[tuple(sorted([i1,i2]))].add(d1+d2)

        # convert to iri links
        iri_links = set()
        for (i1,i2),ds in links.items():
            d = min(ds)
            if d <= distance:
                for iri1,iri2 in it.product(treenumber_iris[treenumbers[i1]], treenumber_iris[treenumbers[i2]]):
                    if iri1 == iri2:
                        continue
                    iri_links.add((iri1, iri2, d))
        return iri_links

    def get_distant_mesh_relatives(self, iri, distance=2):

        def rec_f(tn, distance=2, related_iris=set()):
            if tn in self.treenumber_index:
                for d, related_iri in self.treenumber_index[tn]:
                    if (d <= distance) and (not related_iri in related_iris):
                        related_iris.add(related_iri)
                        for related_tn in self.iri2treenumber[related_iri]:
                            related_iris = rec_f(related_tn, distance=distance-d, related_iris=related_iris)

            return related_iris

        related_iris = set()
        for original_tn in self.iri2treenumber[iri]:
            tn_split = original_tn.split('.')

            for i,idx in enumerate(reversed(range(len(tn_split)))):
                if i > distance:
                    break
                sub_tn = '.'.join(tn_split[:idx+1])
                related_iris.update(rec_f(sub_tn, distance=distance-i, related_iris=related_iris))

        return related_iris

    def get_mesh_name(self, iri):
        query = self.mesh_graph.query(f"SELECT ?o WHERE {{ mesh2021:{iri} <http://www.w3.org/2000/01/rdf-schema#label> ?o }}")
        try:
            return str(list(query)[0][0])
        except IndexError as e:
            return None

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
        
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/treenumber_index.json", 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.treenumber_index.items()}, f)
        with open(f"{data_dir}/iri2treenumber.json", 'wt') as f:
            json.dump({k:list(vs) for k,vs in self.iri2treenumber.items()}, f)
        
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/treenumber_index.json", 'rt') as f:
            self.treenumber_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
        with open(f"{data_dir}/iri2treenumber.json", 'rt') as f:
            self.iri2treenumber = {k:set(vs) for k,vs in json.load(f).items()}
        
        
