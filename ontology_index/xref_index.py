from .onto_index import EfoIndex, MeshIndex, UmlsIndex
from .name_index import NameIndex

class XrefIndex():
    
    def __init__(self, data_dir='.', efo_index=None, mesh_index=None, umls_index=None, name_index=None):
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
            self.name_index = NameIndex(filepath=None, data_dir=self.data_dir, efo_index=efo_index, mesh_index=mesh_index, umls_index=umls_index)
        
        try:
            self.load_indexes()
        except:
            pass
    
    def name_xref(self, iri):
        r = self.name_index.get_names(iri)
        if r:
            iri_names = {n for n,_,_ in r}
        else:
            iri_names = set()
            
        candidates = set()
        for n in iri_names:
            candidates.update(self.name_index.query(n))
        
        filtered_iri_names = {self.name_index.filter_name(n) for n in iri_names}
        
        for c in candidates:
            r = self.name_index.get_names(c)
            if r:
                c_names = {n for n,_,_ in r}
            else:
                c_names = set()
                
            filtered_c_names = {self.name_index.filter_name(n) for n in c_names}
            overlap = filtered_c_names & filtered_iri_names
            scores = [len(overlap)/len(filtered_c_names), len(overlap)/len(filtered_iri_names)]
            yield (
                c,
                max(scores),
                min(scores),
                iri_names,
                c_names,
                overlap
            )
    
    def ontology_xref(self, iri):
        xrefs = set()
        
        r = self.efo_index.get_xrefs(iri)
        if r:
            xrefs.update(r)
        
        r = self.umls_index.get_xrefs(iri)
        if r:
            xrefs.update(r)
        
        return xrefs
        
    
    def get_xrefs(self, iris, covered_iris=None, jumps=1, ontology_based=True, name_based=True, name_xref_score_threshold=0.2):
        if isinstance(iris, str):
            iris = {iris}
        iris = set(iris)
        
        if covered_iris is None:
            covered_iris = iris
        
        xrefs = set()
        
        for iri in iris:
            if ontology_based:
                xrefs.update(self.ontology_xref(iri))  # ontology xrefs
            if name_based:
                for m, max_score, min_score, _, _, _ in self.name_xref(iri):
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
                    name_xref_score_threshold=name_xref_score_threshold
                )
            )
            
        return xrefs
        
        
