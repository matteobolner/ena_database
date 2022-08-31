import os
import requests
import pandas as pd
from tqdm import tqdm
import re
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
import numpy as np

tqdm.pandas()


class ena_db:
    """
    Access ENA api
    """
    def __init__(self, taxon, create_db=True, include_subordinate_taxa=False, breeds='', sexes=['male','female']):
        self.taxon = taxon
        self.sexes= self.get_words(sexes)
        self.breeds = self.get_words(breeds)
        if create_db==True:
            self.ena_df = self.get_ena_records(include_subordinate_taxa)

    def get_ena_records(self, include_subordinate_taxa=False):
        if include_subordinate_taxa==True:
            query_tax = 'tax_tree'
        else:
            query_tax='tax_eq'

        headers = {
          'Content-Type': 'application/x-www-form-urlencoded',
        }

        data = {
        'result': 'read_run',
        'query': f'{query_tax}({self.taxon})',
        'fields':'accession,altitude,assembly_quality,assembly_software,base_count,binning_software,bio_material,broker_name,cell_line,cell_type,center_name,checklist,collected_by,collection_date,collection_date_submitted,completeness_score,contamination_score,country,cram_index_aspera,cram_index_ftp,cram_index_galaxy,cultivar,culture_collection,depth,description,dev_stage,ecotype,elevation,environment_biome,environment_feature,environment_material,environmental_package,environmental_sample,experiment_accession,experiment_alias,experiment_title,experimental_factor,fastq_aspera,fastq_bytes,fastq_ftp,fastq_galaxy,fastq_md5,first_created,first_public,germline,host,host_body_site,host_genotype,host_gravidity,host_growth_conditions,host_phenotype,host_sex,host_status,host_tax_id,identified_by,instrument_model,instrument_platform,investigation_type,isolate,isolation_source,last_updated,lat,library_construction_protocol,library_layout,library_name,library_selection,library_source,library_strategy,location,lon,mating_type,nominal_length,nominal_sdev,parent_study,ph,project_name,protocol_label,read_count,run_accession,run_alias,salinity,sample_accession,sample_alias,sample_capture_status,sample_collection,sample_description,sample_material,sample_title,sampling_campaign,sampling_platform,sampling_site,scientific_name,secondary_sample_accession,secondary_study_accession,sequencing_method,serotype,serovar,sex,specimen_voucher,sra_aspera,sra_bytes,sra_ftp,sra_galaxy,sra_md5,strain,study_accession,study_alias,study_title,sub_species,sub_strain,submission_accession,submission_tool,submitted_aspera,submitted_bytes,submitted_format,submitted_ftp,submitted_galaxy,submitted_host_sex,submitted_md5,submitted_sex,target_gene,tax_id,taxonomic_classification,taxonomic_identity_marker,temperature,tissue_lib,tissue_type,variety',
        'format':'json'
        }

        r = requests.post('https://www.ebi.ac.uk/ena/portal/api/search', headers=headers, data=data)

        r_json =r.json()
        df = pd.DataFrame.from_dict(r_json)
        df.dropna(how='all', axis=1)
        emptycols = []
        for col in df.columns:
            if len(df[col].unique())==1 and df[col].unique()[0]=='':
                emptycols.append(col)

        df=df.drop(emptycols, axis=1)
        return(df)

    def init_pbar(self,length):
        pbar = tqdm(total=length)
        return(pbar)

    def split_df_in_rows(self, df):
        rows_list = [df.loc[i] for i in df.index]
        return(rows_list)

    def get_words(self, input_parameters):
        if isinstance(input_parameters, list):
            words=set(input_parameters)
            return(words)
        elif os.path.exists(input_parameters):
            with open(input_parameters) as f:
                words = f.read().splitlines()
            return(words)
        else:
            raise ValueError(f"Wrong input: {input_parameters}")

    def search_row(self, row, words):
        words = set(w.lower() for w in words)
        words_found=set()
        row_set = set(x.lower() for x in set(row))
        for i in words:
            if str(row_set).find(i) != -1:
                words_found.add(i)
        if len(words_found)>0:
            return(",".join(list(words_found)))
        else:
            return('unknown')

    def get_sex_words(self, df):

        sexwords = set(self.sexes)
        for col in df.columns:
            if 'sex' in col:
                sexwords=sexwords.union(set(df[col].unique()))
        sexwords.remove('')
        return(list(sexwords))

    def get_biosample_from_ena_id(self, id):

        biosamples_server = 'https://www.ebi.ac.uk/biosamples/samples/'
        r = requests.get(biosamples_server+id, headers={ "Content-Type" : "application/json"})
        if r.ok:
            decoded = r.json()
            decoded=pd.json_normalize(decoded).transpose().to_dict()[0]
            #decoded = pd.DataFrame.from_dict([decoded])
            return(decoded)
        else:
            return()

    def build_biosamples_df(self, df, thread_number=30):

        ids = df['sample_accession'].unique().tolist()

        def biosample_with_progress(x):
            biosample=self.get_biosample_from_ena_id(x)
            pbar.update(1)
            return(biosample)

        pbar = tqdm(total=len(ids))
        pool=ThreadPool(thread_number)

        results = pool.map(biosample_with_progress, ids)
        pool.close()
        pool.join()

        biosamples_df = pd.DataFrame(results)
        biosamples_df['domain']=biosamples_df['domain'].str.replace("self.BiosampleImport","")
        biosamples_df.columns=biosamples_df.columns.str.replace("characteristics.","")

        return(biosamples_df)

    def clean_biosamples_df(self, df, thread_number=30):
        pbar = tqdm(total=len(df))
        def clean_row(row):
            acceptables_list = [str, bool, float, int, np.int64, np.bool, np.bool_]
            pbar.update(1)
            newrow=[]
            for col in df.columns:
                temptype = type(row[col])
                if temptype==list:
                    if len(row[col])==1:
                        try:
                            newrow.append(row[col][0]['text'])
                        except:
                            newrow.append(np.nan)
                    else:
                        templist=[]
                        for i in row[col]:
                            try:
                                templist.append(i['text'])
                            except:
                                templist.append(np.nan)
                        newrow.append(templist)
                elif temptype in acceptables_list:
                    newrow.append(row[col])
                else:
                    print(temptype)
                    break
                    #print("ocio")
                    #print(row[col])
            return(newrow)

        biosamples_rows = self.split_df_in_rows(df)

        pool=ThreadPool(thread_number)
        clean_rows = pool.map(clean_row, biosamples_rows)
        pool.close()
        pool.join()
        clean_df = pd.DataFrame(clean_rows,columns=df.columns)

        return(clean_df)

    def all_df_xrefs(self, df, thread_number=30):
        pbar = tqdm(total=len(df['study_accession'].unique()))
        def get_xref(project_id):
            pbar.update(1)
            headers={ "Content-Type" : "application/json"}
            r = requests.get('https://www.ebi.ac.uk/ena/xref/rest/json/search?accession=' + project_id, headers=headers)
            if r.ok:
                decoded = r.json()
            urls = []
            for i in decoded:
                if i['Source']=='EuropePMC' or i['Source']=='PubMed':
                    urls.append(i['Source URL'])
                else:
                    continue
            if len(urls)==0:
                return('NONE')
            else:
                return({project_id:";".join(urls)})
        pool=ThreadPool(thread_number)
        project_ids = df['study_accession'].unique()
        xref_list = pool.map(get_xref, project_ids)
        pool.close()
        pool.join()
        xref_dict = {}
        for i in xref_list:
            if i!='NONE':
                xref_dict.update(i)
        xref_df = pd.DataFrame.from_dict([xref_dict]).transpose()
        xref_df.columns=['xref']
        return(xref_df)
