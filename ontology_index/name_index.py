import re
import json
from collections import defaultdict
from .onto_index import EfoIndex, MeshIndex, UmlsIndex
from tqdm.auto import tqdm

class NameIndex():

    def __init__(self, data_dir='.', efo_index=None, mesh_index=None, umls_index=None):
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
        
        try:
            self.load_indexes()
        except:
            pass
    
    exclude_suffixes = {
        ', nos',
        ',nos',
        ', not otherwise specified',
        ', unspecified',
        ', not elsewhere classified'
    }
    
    def filter_name(self,s):
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
        
        for suffix in self.exclude_suffixes:
            if s[-len(suffix):] == suffix:
                s = normalise_whitespace(s[:-len(suffix)])

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
                name_type = self.umls_index.name_types[umls_name_type]
                self.name_index[filtered_name].add(iri)  # (name, name_type, iri)
                
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        
        with open(f'{data_dir}/name_index.json', 'wt') as f:
            json.dump({k:list(vs) for k,vs in self.name_index.items()}, f)
            
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f'{data_dir}/name_index.json', 'rt') as f:
            self.name_index = {k:set(vs) for k,vs in json.load(f).items()}
            
    def query(self, q, filter_query=True):
        if filter_query:
            q  = self.filter_name(q)
        if q in self.name_index:
            return self.name_index[q]
        
    def get_name(self, iri):
        try:
            name = self.efo_index.get_name(iri)
            if name:
                return name
        except:
            pass

        try:
            name = self.mesh_index.get_name(iri)
            if name:
                return name
        except:
            pass

        try:
            name = self.umls_index.get_name(iri)
            if name:
                return name
        except:
            pass
