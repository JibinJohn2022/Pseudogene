

import numpy as np
import pandas as pd
import json
from SPARQLWrapper import SPARQLWrapper, JSON
import sparql_dataframe  #https://github.com/lawlesst/sparql-dataframe
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


#df=pd.read_csv("13585_S28.hard-filtered",sep="\t")


def check_pseudogene(self,df):
    Replace={'PseudoGeneID':"Pseudogene ID",'Transcript':"Parent Gene Transcript id",'chr':"Pseudogene location chromosome",
    'start':"Pseudogene location start","end":"Pseudogene location end"}
    Columns=['Pseudogene ID','Parent Gene Transcript id','Pseudogene location chromosome', 'Pseudogene location start','Pseudogene location end']
    
    FinalDf=pd.DataFrame(columns = ["Gene Symbol","PSEUDOGENE"])
    GeneNames=list(set(df['SYMBOL'].unique()))
    
    for i in range(0, len(GeneNames), 100):
        inp=tuple(GeneNames[i:i+100])
        query = """
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX obo: <http://purl.obolibrary.org/obo/>
        PREFIX gvaont: <http://semanticwebindia.org/GVA/ont/>
        
        SELECT * FROM <http://PSEUDOGENE.org/> FROM <http://ENSEMBL.org/> where 
        {
            ?ParentGene rdfs:label ?ParentGene_Symbol . FILTER (?ParentGene_Symbol IN %(inp)s) .
            ?PseudoGeneID gvaont:hasGene ?ParentGene ;
                          gvaont:hasGeneTranscript ?Transcript ; 
                          obo:RO_0002231 ?start ;
                          obo:RO_0002232 ?end ;
                          obo:RO_0001025 ?chr .
        }
        """
        query = query % {'inp':inp}
        endpoint = "http://18.223.39.126:8890/sparql"
        PseudogeneDf= sparql_dataframe.get(endpoint, query)
        PseudogeneDf.rename(columns=Replace,inplace=True)
        PseudogeneDf['Pseudogene ID']=PseudogeneDf['Pseudogene ID'].str.replace("http://semanticwebindia.in/GVA/res/pseudogene/id/",'')
        for gene in PseudogeneDf['ParentGene_Symbol'].unique():
            result={}
            tempdf=PseudogeneDf[PseudogeneDf['ParentGene_Symbol']==gene]
            Jsondf=tempdf[Columns].reset_index().drop('index',axis=1)
            Json=json.loads(Jsondf.to_json(orient="index"))
            result["PSEUDOGENE"]=str(Json)
            result["Gene Symbol"]=gene
            result=pd.DataFrame([result])
            FinalDf=pd.concat([FinalDf, result], ignore_index=True)
    VEP_Transcript_Psedo_df=pd.merge(df,FinalDf,left_on='SYMBOL',right_on="Gene Symbol",how="left").drop("Gene Symbol",axis=1)

    return VEP_Transcript_Psedo_df





