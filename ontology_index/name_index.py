import re
import json
from collections import defaultdict
from .onto_index import EfoIndex, MeshIndex, UmlsIndex

class NameIndex():
    
    umls_name_types = {
        'PF': 'umls:pref_term',
        'VCW': 'umls:case_word_order_variant',
        'VC': 'umls:case_variant',
        'VO': 'umls:variant',
        'VW': 'umls:word_order_variant',
    }

    def __init__(self, data_dir='.'):
        self.data_dir = data_dir
        
        self.efo_index = EfoIndex(data_dir=self.data_dir)
        self.mesh_index = MeshIndex(data_dir=self.data_dir)
        self.umls_index = UmlsIndex(filepath=None, data_dir=self.data_dir)
        
        try:
            self.load_indexes()
        except:
            pass
        
    
    def filter_name(s):
        def normalise_whitespace(s):
            s = re.sub('\s+', ' ', s)
            s = re.sub('\s+$', '', s)
            s = re.sub('^\s+', '', s)
            return s

        s = s.lower()
        if re.match('.*(\s|^)\(.*\)(\s|$).*', s):
            s = re.sub('(\s|^)\(.*\)(\s|$)', ' ', s)

            # tidy up problems that occur due to removing brackets
            s = normalise_whitespace(s)

            if bool(s):
                if s[-1] == ',':
                    s = normalise_whitespace(s[:-1])
                if s[0] == ',':
                    s = normalise_whitespace(s[1:])

        else:
            s = normalise_whitespace(s)

        if re.match('.*(\s|^)\[.*\](\s|$).*', s):
            s = re.sub('(\s|^)\[.*\](\s|$)', ' ', s)

            # tidy up problems that occur due to removing brackets
            s = normalise_whitespace(s)

            if bool(s):
                if s[-1] == ',':
                    s = normalise_whitespace(s[:-1])
                if s[0] == ',':
                    s = normalise_whitespace(s[1:])

        else:
            s = normalise_whitespace(s)

        s = re.sub('(\s|^)&(\s|$)', 'and', s)
        s = re.sub('[/\-]', ' ', s)
        s = re.sub('[^a-z., ]', '', s)
        s = re.sub('(\s|^)\.\s', ' ', s)

        return s
    
    def gen_query_index(self):
        self.name_index = defaultdict(set)
        
        for iri, d in tqdm(self.efo_index.iri2name.items(), leave=True, position=0):
            for name_type, name in d:
                filtered_name = self.filter_name(name)
                self.name_index[filtered_name].add(iri)  # (name, name_type, iri)
                
        for iri, d in tqdm(self.mesh_index.iri2name.items(), leave=True, position=0):
            for name_type, name in d:
                filtered_name = self.filter_name(name)
                self.name_index[filtered_name].add(iri)  # (name, name_type, iri)

        for iri, d in tqdm(self.umls_index.iri2name.items(), leave=True, position=0):
            for name_type, umls_name_type, name in d:
                filtered_name = self.filter_name(name)
                name_type = self.umls_name_types[umls_name_type]
                self.name_index[filtered_name].add(iri)  # (name, name_type, iri)
                
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        
        with open(f'{data_dir}/data/name_index.json', 'wt') as f:
            json.dump({k:[list(v) for v in vs] for k,vs in self.name_index.items()}, f)
            
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f'{data_dir}/name_index.json', 'rt') as f:
            self.name_index = {k:{tuple(v) for v in vs} for k,vs in json.load(f).items()}
            
    def query(self, q, filter_query=True):
        if filter_query:
            q  = self.filter_name(q)
        if q in self.name_index:
            return self.name_index[q]
