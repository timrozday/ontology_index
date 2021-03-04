import re
import json
import pickle
from collections import defaultdict
from .onto_index import EfoIndex, MeshIndex, UmlsIndex
from tqdm.auto import tqdm

import requests
import urllib
import re

class TextFilter():
    def __init__(self):
        pass
    
    def normalise_whitespace(self, s):
        s = re.sub('\s+', ' ', s)
        s = re.sub('\s+$', '', s)
        s = re.sub('^\s+', '', s)
        return s

    exclude_suffixes = {
        'nos',
        'nec',
        'not otherwise specified',
        'unspecified',
        'not elsewhere classified'
    }
    
    def filter_name(self, s):

        s = s.lower()
        if re.match('.*(\s|^)\([^)]*\)(\s|$).*', s):
            s = re.sub('(\s|^)\([^)]*\)(\s|$)', ' ', s)

            # tidy up problems that occur due to removing brackets
            s = self.normalise_whitespace(s)

            if bool(s):
                if s[-1] == ',':
                    s = self.normalise_whitespace(s[:-1])
            if bool(s):
                if s[0] == ',':
                    s = self.normalise_whitespace(s[1:])

        else:
            s = self.normalise_whitespace(s)

#         if re.match('.*(\s|^)\[[^\]]*\](\s|$).*', s):
#             s = re.sub('(\s|^)\[[^\]]*\](\s|$)', ' ', s)

#             # tidy up problems that occur due to removing brackets
#             s = self.normalise_whitespace(s)

#             if bool(s):
#                 if s[-1] == ',':
#                     s = self.normalise_whitespace(s[:-1])
#             if bool(s):
#                 if s[0] == ',':
#                     s = self.normalise_whitespace(s[1:])

#         else:
#             s = self.normalise_whitespace(s)

        s = re.sub('(\s|^)&(\s|$)', 'and', s)  # normalise ands
        s = re.sub('[/\-]', ' ', s)
        s = re.sub('[\']', '', s)
        s = re.sub('[^a-z0-9., ]', ' ', s)
        s = re.sub('(\s|^)\.(\s|$)', ' ', s)
        s = self.normalise_whitespace(s)
        
        for suffix in self.exclude_suffixes:
            if s[-len(suffix):] == suffix:
                s = self.normalise_whitespace(s[:-len(suffix)])
        
        if bool(s):
            if s[-1] == ',':
                s = self.normalise_whitespace(s[:-1])
        if bool(s):
            if s[0] == ',':
                s = self.normalise_whitespace(s[1:])
        
        return s

    def remove_punctuation(self, s):
        s = re.sub('[^a-z0-9 ]', '', s)  # remove punctuation
        s = self.normalise_whitespace(s)

        return s
    
    def trim(self, s):
        change = True
        while change:
            change = False
            if bool(s):
                if s[-1] in {',','.'}:
                    s = s[:-1]
                    change = True
            if bool(s):
                if s[0] == {',','.'}:
                    s = s[1:]
                    change = True
                    
        return s
    
    def tokenize(self, s):
        return [self.trim(t) for t in re.split('(?<=\S)[\s](?=\S)', self.normalise_whitespace(s))]
    

class NameIndex(TextFilter):

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
    
    def gen_query_index(self):
        self.name_index = defaultdict(set)
        self.iri_name_index = defaultdict(set)
        
        for iri, d in tqdm(self.efo_index.iri2name.items(), leave=True, position=0):
            for name_type, name in d:
                filtered_name = self.filter_name(name)
                if filtered_name:
                    if name_type in {
                        'http://www.w3.org/2000/01/rdf-schema#label',
                        'http://www.w3.org/2004/02/skos/core#prefLabel',
                        'http://www.geneontology.org/formats/oboInOwl#hasExactSynonym',
                        'http://www.geneontology.org/formats/oboInOwl#shorthand'
                    }:
                        self.name_index[filtered_name].add(iri)  # (name, name_type, iri)
                        tokens = self.tokenize(self.remove_punctuation(filtered_name))  # remove all punctuation
                        self.iri_name_index[iri].add((name, filtered_name, tuple(tokens)))
                
        for iri, d in tqdm(self.mesh_index.iri2name.items(), leave=True, position=0):
            for name_type, name in d:
                filtered_name = self.filter_name(name)
                if filtered_name:
                    if name_type in {
                        'http://id.nlm.nih.gov/mesh/vocab#prefLabel',
                        'http://www.w3.org/2000/01/rdf-schema#label',
                        'http://id.nlm.nih.gov/mesh/vocab#altLabel'
                    }:
                        self.name_index[filtered_name].add(iri)  # (name, name_type, iri)
                        tokens = self.tokenize(self.remove_punctuation(filtered_name))  # remove all punctuation
                        self.iri_name_index[iri].add((name, filtered_name, tuple(tokens)))

        for iri, d in tqdm(self.umls_index.iri2name.items(), leave=True, position=0):
            for name_type, umls_name_type, name in d:
                filtered_name = self.filter_name(name)
                if filtered_name:
                    name_type = self.umls_index.name_types[umls_name_type]
                    if name_type in {
                        'umls:pref_term',
                        'umls:case_word_order_variant',
                        'umls:case_variant',
                        'umls:variant',
                        'umls:word_order_variant'
                    }:
                        self.name_index[filtered_name].add(iri)  # (name, name_type, iri)
                        tokens = self.tokenize(self.remove_punctuation(filtered_name))  # remove all punctuation
                        self.iri_name_index[iri].add((name, filtered_name, tuple(tokens)))
                        
        self.gen_kmer_index()
    
    def gen_kmer_index(self, size_limit=None):
        def gen_kmers(l, k=3):
            if len(l) < k:
                yield tuple(sorted(l))
            for i in range(len(l)-k+1):
                yield tuple(sorted(l[i:i+k]))

        token_index = defaultdict(lambda :defaultdict(set))
        for iri,data in tqdm(self.iri_name_index.items(), position=0, leave=True, desc="Generating kmer index"):
            for name, filtered_name, tokens in data:
                for kmer in gen_kmers(tokens):
                    token_index[len(kmer)][kmer].add(iri)
        
        if size_limit:
            self.token_index = {k1:{k2:v2 for k2,v2 in v1.items() if len(v2)<=size_limit} for k1,v1 in token_index.items()}
        else:
            self.token_index = {k1:{k2:v2 for k2,v2 in v1.items()} for k1,v1 in token_index.items()}
    
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        
        with open(f'{data_dir}/name_index.json', 'wt') as f:
            json.dump({k:list(vs) for k,vs in self.name_index.items()}, f)
        with open(f'{data_dir}/iri_name_index.json', 'wt') as f:
            json.dump({k:list(vs) for k,vs in self.iri_name_index.items()}, f)
            
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f'{data_dir}/name_index.json', 'rt') as f:
            self.name_index = {k:set(vs) for k,vs in json.load(f).items()}
        with open(f'{data_dir}/iri_name_index.json', 'rt') as f:
            self.iri_name_index = {k:{(n,f,tuple(t)) for n,f,t in vs} for k,vs in json.load(f).items()}
            
        self.gen_kmer_index()
            
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
        
    def get_names(self, iri):
#         try:
#             names = self.efo_index.get_names(iri)
#             if names:
#                 return names, 'efo'
#         except:
#             pass

#         try:
#             names = self.mesh_index.get_names(iri)
#             if names:
#                 return names, 'mesh'
#         except:
#             pass

#         try:
#             names = self.umls_index.get_names(iri)
#             if names:
#                 return names, 'umls'
#         except:
#             pass
        
        if iri in self.iri_name_index:
            return self.iri_name_index[iri]

        
        
class QualifierIndex(TextFilter):
    """For extraction of allowed qualifiers from indications

From `NCIT` and `HPO`.

`NCIT`:
* `http://purl.obolibrary.org/obo/NCIT_C41009` Qualifier
  * `http://purl.obolibrary.org/obo/NCIT_C13442` Anatomical qualifier
  * `http://purl.obolibrary.org/obo/NCIT_C21514` Temporal qualifier
  * `http://purl.obolibrary.org/obo/NCIT_C73706` Spatial qualifier
  * `http://purl.obolibrary.org/obo/NCIT_C27992` Disease qualifier
    * `http://purl.obolibrary.org/obo/NCIT_C28102` Disease clinical qualifier
  * `http://purl.obolibrary.org/obo/NCIT_C27993` General qualifier (contains 'primary')
  
`HPO`:
* `http://purl.obolibrary.org/obo/HP_0031797` Clinical course
  * `http://purl.obolibrary.org/obo/HP_0011008` Temporal pattern
  * `http://purl.obolibrary.org/obo/HP_0003679` Pace of progression
  * `http://purl.obolibrary.org/obo/HP_0003674` Onset
* `http://purl.obolibrary.org/obo/HP_0012823` Clinical modifier
  * `http://purl.obolibrary.org/obo/HP_0025280` Pain characteristic
  * `http://purl.obolibrary.org/obo/HP_0012824` Severity
* `http://purl.obolibrary.org/obo/HP_0040279` Frequency

Miscellaneous
* `http://purl.obolibrary.org/obo/NCIT_C159950` Relapse
* `http://purl.obolibrary.org/obo/NCIT_C43623` Persistent
* `http://purl.obolibrary.org/obo/NCIT_C25251` Primary
* `http://purl.obolibrary.org/obo/NCIT_C14174` Metastatic
* `http://purl.obolibrary.org/obo/NCIT_C25378` Partial
* `http://purl.obolibrary.org/obo/NCIT_C28242` Idiopathic
* `http://purl.obolibrary.org/obo/NCIT_C34340` Accidental
  """
    
    def __init__(self, data_dir='.'):
        self.data_dir = data_dir
        
        try:
            self.load_indexes()
        except:
            pass

    def get_iri_terms(self, q):
        iri, source = q
        iri_str = urllib.parse.quote(iri, safe='')
        iri_str = urllib.parse.quote(iri_str, safe='')
        r = requests.get(f"https://www.ebi.ac.uk/ols/api/ontologies/{source}/terms/{iri_str}", params={'size':10, 'page':1})
        return r.json()

    def get_iri_decendents(self, q):
        iri, source = q
        iri_str = urllib.parse.quote(iri, safe='')
        iri_str = urllib.parse.quote(iri_str, safe='')
        data = []
        p = 0
        while True:
            r = requests.get(f"https://www.ebi.ac.uk/ols/api/ontologies/{source}/terms/{iri_str}/hierarchicalDescendants", params={'size': 50, 'page': p}).json()
            if '_embedded' in r:
                data.extend(r['_embedded']['terms'])
            if r['page']['number'] >= r['page']['totalPages'] -1:
                break
            else:
                p+=1
        return data
    
    def gen_indexes(self, ncit=True, hpo=True, miscellaneous=True):

        hp_qualifiers = {}
        if hpo:
            for iri in tqdm({
                'http://purl.obolibrary.org/obo/HP_0031797',
                'http://purl.obolibrary.org/obo/HP_0011008', 
                'http://purl.obolibrary.org/obo/HP_0003679', 
                'http://purl.obolibrary.org/obo/HP_0003674', 
                'http://purl.obolibrary.org/obo/HP_0012823',
                'http://purl.obolibrary.org/obo/HP_0025280', 
                'http://purl.obolibrary.org/obo/HP_0012824', 
                'http://purl.obolibrary.org/obo/HP_0040279',
            }, leave=True, position=0, desc="HPO qualifiers"):
                for r in self.get_iri_decendents((iri, 'hp')):
                    hp_qualifiers[(r['iri'], r['ontology_name'])] = {r['label']}
                    if 'synonyms' in r:
                        if r['synonyms']:
                            hp_qualifiers[(r['iri'], r['ontology_name'])].update(r['synonyms'])

        ncit_qualifiers = {}
        if ncit:
            for iri in tqdm({
            #     'http://purl.obolibrary.org/obo/NCIT_C41009', # Qualifier
                'http://purl.obolibrary.org/obo/NCIT_C13442', # Anatomical qualifier
                'http://purl.obolibrary.org/obo/NCIT_C21514', # Temporal qualifier
            #     'http://purl.obolibrary.org/obo/NCIT_C73706', # Spatial qualifier
                'http://purl.obolibrary.org/obo/NCIT_C27992', # Disease qualifier
                'http://purl.obolibrary.org/obo/NCIT_C28102', # Disease clinical qualifier
                'http://purl.obolibrary.org/obo/NCIT_C27993', # General qualifier (contains 'primary')
            }, leave=True, position=0, desc="NCIT qualifiers"):
                for r in self.get_iri_decendents((iri, 'ncit')):
                    ncit_qualifiers[(r['iri'], r['ontology_name'])] = {r['label']}
                    if 'synonyms' in r:
                        if r['synonyms']:
                            ncit_qualifiers[(r['iri'], r['ontology_name'])].update(r['synonyms'])

        if miscellaneous:
            custom_qualifiers = {
                ('http://purl.obolibrary.org/obo/NCIT_C159950', 'ncit'): {'relapse', 'relapsed'},
                ('http://purl.obolibrary.org/obo/NCIT_C43623', 'ncit'): {'persistent'},
                ('http://purl.obolibrary.org/obo/NCIT_C25251', 'ncit'): {'primary'},
                ('http://purl.obolibrary.org/obo/NCIT_C14174', 'ncit'): {'metastatic'},
                ('http://purl.obolibrary.org/obo/NCIT_C25378', 'ncit'): {'partial'},
                ('http://purl.obolibrary.org/obo/NCIT_C28242', 'ncit'): {'idiopathic'},
                ('http://purl.obolibrary.org/obo/NCIT_C34340', 'ncit'): {'accidental'},
            }
        else:
            custom_qualifiers = {}

        self.ols_qualifiers = {**hp_qualifiers, **custom_qualifiers, **ncit_qualifiers}

        token_qualifier_index = defaultdict(set)
        for (iri, source), terms in self.ols_qualifiers.items():
            for s in terms:
                for t in self.tokenize(self.filter_name(s)):
                    token_qualifier_index[t].add((s, (iri, source)))

        self.token_qualifier_index = dict(token_qualifier_index)

        
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f'{data_dir}/ols_token_qualifier_index.pkl', 'wb') as f:
            pickle.dump(self.token_qualifier_index, f)
        with open(f'{data_dir}/ols_qualifiers.pkl', 'wb') as f:
            pickle.dump(self.ols_qualifiers, f)
            
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f'{data_dir}/ols_token_qualifier_index.pkl', 'rb') as f:
            self.token_qualifier_index = pickle.load(f)
        with open(f'{data_dir}/ols_qualifiers.pkl', 'rb') as f:
            self.ols_qualifiers = pickle.load(f)

            
    def extract_qualifiers(self, q):
        q_tokens = set(self.tokenize(self.filter_name(q)))
        candidate_matches = set()
        for t in q_tokens:
            if t in self.token_qualifier_index:
                candidate_matches.update(self.token_qualifier_index[t])

    #     return candidate_matches
        candidate_matches_2 = set()
        for m,(iri,source) in candidate_matches:
            m_tokens = set(self.tokenize(self.filter_name(m)))
            if not m_tokens - q_tokens:
                candidate_matches_2.add((m,(iri,source)))

        filter_q = re.sub('-', ' ', q.lower())
        qualifiers = set()
        for m,(iri,source) in sorted(candidate_matches, key=lambda x:(len(x[0]),x[0],x), reverse=True):
            matches = list(re.finditer(self.filter_name(m), filter_q))
            if matches:
                qualifiers.add((m,(iri,source)))
                match_spans = [match.span() for match in matches]
                for (start, end) in sorted(match_spans, reverse=True):
                    q = self.normalise_whitespace(q[:start] + q[end:])
                    filter_q = re.sub('-', ' ', q.lower())


        return q, tuple(qualifiers)
