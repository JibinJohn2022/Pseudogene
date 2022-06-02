
import pandas as pd
#from itertools import Predicate
from rdflib import Namespace
from rdflib.namespace import CSVW, DC, DCAT, DCTERMS, DOAP, FOAF, ODRL2, ORG, OWL, \
                           PROF, PROV, RDF, RDFS, SDO, SH, SKOS, SOSA, SSN, TIME, \
                           VOID, XMLNS, XSD
from rdflib import Graph, URIRef, Literal, BNode



Schema=pd.read_csv("/mnt/c/Users/ikror/OneDrive/Desktop/Projects/PseudoGene/PseudoGeneGraph/PseudoGeneSchema_Long.csv",encoding='cp1252')
#Schema.head()
Data=pd.read_csv("/mnt/c/Users/ikror/OneDrive/Desktop/Projects/PseudoGene/PseudoGeneGraph/psiDR.v0_ForGraph.csv",sep=",", dtype=str)

#Data["TaxID"]="NCBITaxon_9606"


ListcOlumns=list(Schema["Subjectcolumns"].unique())+list(Schema["ObjectColumn"].unique())
ListcOlumns=[x for x in ListcOlumns if x != 'NotApplicable']

#Data=Data.head(n=10)

Data=Data[ListcOlumns]
g = Graph()
Subjeclist=set(Data[",".join(Schema['Subjectcolumns'].unique())].to_list())

for subject in Subjeclist:
    for i, j in Schema.iterrows():
        #Subject
        SPrefix=j['Subject']
        SubjectValue=subject
        Subject = URIRef(SPrefix+SubjectValue)
        
        #Predicate
        Predicate=j['Predicate']
        
        
        #object_ObjectPrefix
        ObPrefix=j['object_ObjectPrefix']
        #object_ObjectPrefix Value
        if  j['object_ObjectPrefix']=="NotApplicable":
            object_ObjectPrefix=ObPrefix
            g.add((Subject, URIRef(Predicate), URIRef(object_ObjectPrefix)))
        
        elif j['object_ObjectPrefix']=="Literal()":
            object_ObjectPrefixSufixs=set(Data[Data[j['Subjectcolumns']]==subject][j['ObjectColumn']].to_list())
            for object_ObjectPrefixSufix in object_ObjectPrefixSufixs:
                object_ObjectPrefix=Literal(object_ObjectPrefixSufix)
                g.add((Subject, URIRef(Predicate), object_ObjectPrefix))
        elif "http" in j['object_ObjectPrefix'] :
            if j['ObjectColumn']=="NotApplicable":
                object_ObjectPrefix=ObPrefix
                g.add((Subject, URIRef(Predicate), URIRef(object_ObjectPrefix)))
            else:
                object_ObjectPrefixSufixs=[x for x in set(Data[Data[j['Subjectcolumns']]==subject][j['ObjectColumn']].to_list()) if str(x) != 'nan']
                #print(object_ObjectPrefixSufixs)
                for object_ObjectPrefixSufix in object_ObjectPrefixSufixs:
                    object_ObjectPrefix=URIRef(j['object_ObjectPrefix']+str(object_ObjectPrefixSufix))
                    g.add((Subject, URIRef(Predicate), object_ObjectPrefix))



for s, p, o in g:
    print(s, p, o)


s = g.serialize(format="ttl",destination="/mnt/c/Users/ikror/OneDrive/Desktop/Projects/PseudoGene/PseudoGeneGraph/PseudoGenepsiDR-v0.ttl")
s = g.serialize(format="nt",destination="/mnt/c/Users/ikror/OneDrive/Desktop/Projects/PseudoGene/PseudoGeneGraph/PseudoGenepsiDR-v0.nt")

print(s)






Same.
ip is swiadmin@18.119.94.81
pass: chanakyapuri@108
      chanakyapuri@108















