from .onto_index import EfoIndex, MeshIndex, UmlsIndex
from .name_index import NameIndex, QualifierIndex

from collections import defaultdict

class XrefIndex():
    
    def __init__(self, data_dir='.', efo_index=None, mesh_index=None, umls_index=None, name_index=None, qualifier_index=None):
        self.data_dir = data_dir
        
        if efo_index:
            self.efo_index = efo_index
        else:
            self.efo_index = EfoIndex(data_dir=self.data_dir)
            
        if mesh_index:
            self.mesh_index = mesh_index
        else:
            self.mesh_index = MeshIndex(data_dir=self.data_dir)
            
        if umls_index:
            self.umls_index = umls_index
        else:
            self.umls_index = UmlsIndex(filepath=None, data_dir=self.data_dir)
        
        if name_index:
            self.name_index = name_index
        else:
            self.name_index = NameIndex(data_dir=self.data_dir, efo_index=efo_index, mesh_index=mesh_index, umls_index=umls_index)
        
        if qualifier_index:
            self.qualifier_index = qualifier_index
        else:
            self.qualifier_index = QualifierIndex(data_dir=self.data_dir)
        
        try:
            self.load_indexes()
        except:
            pass
    
    def name_xref(self, iri, min_length=4, extract_qualifiers=True):
        def get_names(iri, min_length=4):
            r = self.name_index.get_names(iri)
            if r:
                iri_names = {filtered_name for name, filtered_name, tokens in r}
            else:
                iri_names = set()

            if min_length:
                iri_names = {n for n in iri_names if len(n)>=min_length}

            return iri_names
        
        iri_names = get_names(iri, min_length=min_length)
            
        candidates = defaultdict(set)
        for n in iri_names:
            r = self.name_index.query(n)
            if r:
                candidates[None].update(r)
            if extract_qualifiers:
                new_n, quals = self.qualifier_index.extract_qualifiers(n)
                r = self.name_index.query(new_n)
                if r:
                    candidates[tuple(quals)].update(r)
        
        filtered_iri_names = {self.name_index.filter_name(n) for n in iri_names}
        if extract_qualifiers:
            filtered_iri_names = {self.qualifier_index.extract_qualifiers(n)[0] for n in filtered_iri_names}
        
        for quals, qual_candidates in candidates.items():
            for c in qual_candidates:
                c_names = get_names(c, min_length=min_length)
                if c_names:
                    filtered_c_names = {self.name_index.filter_name(n) for n in c_names}
                    if extract_qualifiers:
                        filtered_c_names = {self.qualifier_index.extract_qualifiers(n)[0] for n in filtered_c_names}
                        
                    overlap = filtered_c_names & filtered_iri_names
                    scores = [len(overlap)/len(filtered_c_names), len(overlap)/len(filtered_iri_names)]
                    yield (
                        c,
                        max(scores),
                        min(scores),
                        iri_names,
                        c_names,
                        overlap,
                        quals
                    )
    
    def ontology_xref(self, iri, equivalents=True):
        xrefs = set()
        
        r = self.efo_index.get_xrefs(iri)
        if r:
            xrefs.update({o for p,o in r})
        
        r = self.umls_index.get_xrefs(iri)
        if r:
            xrefs.update({o for p,o in r})
        
        if equivalents:
            try:
                r = self.efo_index.get_distant_efo_relatives(iri, distance=0, distant_rels={'close', 'child', 'parent'}, equivalent_rels={'equivalent'})
                if r:
                    xrefs.update(r.keys())
            except:
                pass
            
            try:
                r = self.mesh_index.get_distant_mesh_relatives(iri.split('/')[-1], distance=0, search_up=True)
                if r:
                    xrefs.update({f"http://id.nlm.nih.gov/mesh/2021/{i}" for i in r})
            except:
                pass
            
        return xrefs
        
    
    def get_xrefs(self, iris, covered_iris=None, jumps=1, ontology_based=True, name_based=True, extract_qualifiers=True, name_xref_score_threshold=0.05, equivalents=True, min_name_length=4):
        if isinstance(iris, str):
            iris = {iris}
        iris = set(iris)
        
        if covered_iris is None:
            covered_iris = iris
        
        xrefs = set()
        
        for iri in iris:
            if ontology_based:
                xrefs.update(self.ontology_xref(iri, equivalents=equivalents))  # ontology xrefs
            if name_based:
                for m, max_score, min_score, _, _, _, quals in self.name_xref(iri, min_length=min_name_length, extract_qualifiers=extract_qualifiers):
                    if max_score >= name_xref_score_threshold:
                        xrefs.add(m)  # name-based xrefs
            
        new_xrefs = xrefs - covered_iris
        
        if new_xrefs and not (jumps==0 or jumps==1):
            xrefs.update(
                self.get_xrefs(
                    new_xrefs, 
                    covered_iris=covered_iris|xrefs, 
                    jumps=jumps-1, 
                    ontology_based=ontology_based, 
                    name_based=name_based, 
                    name_xref_score_threshold=name_xref_score_threshold,
                    min_name_length=min_name_length
                )
            )
            
        return xrefs
        
        
