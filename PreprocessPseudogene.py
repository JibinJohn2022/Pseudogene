

import numpy as np
import pandas as pd
from SPARQLWrapper import SPARQLWrapper, JSON

Pgene="/mnt/c/Users/ikror/OneDrive/Desktop/Projects/PseudoGene/RawData/psiDR.v0.txt"  # DFrom http://pseudogene.org/psidr/psiDR.v0.txt ; Downloaded date: 09pwd-05-2022


#Prepare pseudogene data
pgenedf=pd.read_csv(Pgene,sep="\t",skiprows=22)
pgeneDcolumns=pgenedf[[0,8]].dropna()
pgeneDcolumns.columns=["PseudogeneID","EnsembleGeneID"]

pgeneDcolumns['Pseudogenes'] = pgeneDcolumns.groupby(['EnsembleGeneID'])['PseudogeneID'].transform(lambda x : ','.join(x))
PgeneGrouped=pgeneDcolumns[["EnsembleGeneID","Pseudogenes"]].drop_duplicates().reset_index().drop("index",axis=1)



## Extract gene name from http://ENSEMBL.org graph database
query = """
      PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> 
      PREFIX owl: <http://www.w3.org/2002/07/owl#>  
      PREFIX ensembl: <http://ensembl.org/id/>

      SELECT ?EID ?SYMBOL ?HGNCID FROM <http://ENSEMBL.org/> where 
    {
	    ?EID  rdfs:label ?SYMBOL .
	    ?EID  owl:equivalentClass ?HGNCID.
    } 
"""
sparql = SPARQLWrapper("http://18.223.39.126:8890/sparql")
sparql.setQuery(query)
sparql.setReturnFormat(JSON)
result = sparql.query().convert()

## Convert to df
Ensembldf=pd.json_normalize(result['results']['bindings'])[['EID.value','SYMBOL.value']]
Ensembldf.columns=["EID",	"SYMBOL"]
Ensembldf["EID"]=Ensembldf["EID"].str.split("/").str[-1]



## Merge Psudogene and Gene name with df
Pgene_Ensembl=pd.merge(PgeneGrouped,Ensembldf,left_on="EnsembleGeneID",right_on="EID",how="left").drop("EID",axis=1)
PseudogeneWithSymbol=Pgene_Ensembl[~Pgene_Ensembl["SYMBOL"].isnull()]
PseudogeneWithSymbol=PseudogeneWithSymbol.set_index("EnsembleGeneID")

#Since Gene symbol has duplicate values we proform group by to make gene symbol uniq
PseudogeneWithSymbol1=PseudogeneWithSymbol.reset_index().drop("EnsembleGeneID",axis=1)
PseudogeneWithSymbol1['PseudogeneID']=PseudogeneWithSymbol1.groupby(['SYMBOL'])['Pseudogenes'].transform(lambda x : ','.join(x))
PseudogeneWithSymbol1=PseudogeneWithSymbol1[["SYMBOL","PseudogeneID"]].drop_duplicates().reset_index().drop("index",axis=1)
PseudogeneWithSymbol.to_csv("PseudogenesWithGeneSymbols_EnsemblGraph.csv",index=None)


#Convert Pseudogene with Gene symbol to json file 
dictionary=PseudogeneWithSymbol1.set_index('SYMBOL').to_dict()["PseudogeneID"]
jsonString =json.dumps(dictionary) 
jsonFile = open("PseudoGenewithGeneSymbol.json", "w")
jsonFile.write(jsonString)
jsonFile.close()



#Genes without Ensembl ID at Ensembl graph database
PseudogeneWithoutSymbol=Pgene_Ensembl[Pgene_Ensembl["SYMBOL"].isnull()]
PseudogeneWithoutSymbol=PseudogeneWithoutSymbol.set_index("EnsembleGeneID").drop(["EID","SYMBOL"],axis=1)

PseudogeneWithoutSymbol.to_csv("PseudogenesWithoutGeneSymbols_in_EnsemblGraph.csv",index=None)


sudo echo "export PATH=$PATH:/home/jibin/.local/bin" >> /home/jibin/.bash_profile







