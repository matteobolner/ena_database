import requests
import pandas as pd
from tqdm import tqdm
import io
import re
import numpy as np
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
import unicodedata
tqdm.pandas()

class EnaBiosamplesXrefDB:
    """
    Setup ena, biosamples and xref DB for taxon
    """
    def __init__(
        self,
        taxon,
        species,
        include_subordinate_taxa=True,
        ena_reads_db_path=None,
        biosamples_db_path=None,
        xrefs_db_path=None,
        n_samples='all', #for testing
    ):
        self.taxon = taxon
        self.species= species
        if ena_reads_db_path:
            print("Reading ENA db from file...")
            self.ena_reads_db=pd.read_table(ena_reads_db_path)
        else:
            self.ena_reads_db = self.build_ena_reads_db(taxon, include_subordinate_taxa, n_samples)

        self.samples=self.ena_reads_db['sample_accession'].dropna().unique().tolist()
        self.studies=self.ena_reads_db['study_accession'].dropna().unique().tolist()
        if biosamples_db_path:
            print("Reading BioSamples db from file...")
            self.biosamples_db=pd.read_table(biosamples_db_path)
        else:
            self.biosamples_db = self.build_biosamples_db(self.samples)
        if xrefs_db_path:
            print("Reading Xrefs db from file...")
            self.xrefs_db=pd.read_table(xrefs_db_path)
            print("Finished database setup")
        else:
            self.xrefs_db = self.build_study_xrefs_db(self.studies)

    def build_ena_reads_db(self, taxon, include_subordinate_taxa, n_samples):
        """
        Query ENA for all reads belonging to the input taxa
        """

        if include_subordinate_taxa == True:
            query_tax = "tax_tree"
            print(f"Searching ENA for taxon {self.taxon} and subordinate taxa...")
        else:
            query_tax = "tax_eq"
            print(f"Searching ENA for taxon {self.taxon}...")

        headers = {
            "Content-Type": "application/x-www-form-urlencoded",
        }

        result='read_run'
        query=f"{query_tax}({self.taxon})"
        format='tsv'
        fields="run_accession%2Cexperiment_title%2Ctax_id%2Cage%2Caligned%2Caltitude%2Cassembly_quality%2Cassembly_software%2Cbam_aspera%2Cbam_bytes%2Cbam_ftp%2Cbam_galaxy%2Cbam_md5%2Cbase_count%2Cbinning_software%2Cbio_material%2Cbisulfite_protocol%2Cbroad_scale_environmental_context%2Cbroker_name%2Ccage_protocol%2Ccell_line%2Ccell_type%2Ccenter_name%2Cchecklist%2Cchip_ab_provider%2Cchip_protocol%2Cchip_target%2Ccollected_by%2Ccollection_date%2Ccollection_date_end%2Ccollection_date_start%2Ccompleteness_score%2Ccontamination_score%2Ccontrol_experiment%2Ccountry%2Ccultivar%2Cculture_collection%2Cdatahub%2Cdepth%2Cdescription%2Cdev_stage%2Cdisease%2Cdnase_protocol%2Cecotype%2Celevation%2Cenvironment_biome%2Cenvironment_feature%2Cenvironment_material%2Cenvironmental_medium%2Cenvironmental_sample%2Cexperiment_accession%2Cexperiment_alias%2Cexperiment_target%2Cexperimental_factor%2Cexperimental_protocol%2Cextraction_protocol%2Cfaang_library_selection%2Cfastq_aspera%2Cfastq_bytes%2Cfastq_ftp%2Cfastq_galaxy%2Cfastq_md5%2Cfile_location%2Cfirst_created%2Cfirst_public%2Cgermline%2Chi_c_protocol%2Chost%2Chost_body_site%2Chost_genotype%2Chost_gravidity%2Chost_growth_conditions%2Chost_phenotype%2Chost_scientific_name%2Chost_sex%2Chost_status%2Chost_tax_id%2Cidentified_by%2Cinstrument_model%2Cinstrument_platform%2Cinvestigation_type%2Cisolate%2Cisolation_source%2Clast_updated%2Clat%2Clibrary_construction_protocol%2Clibrary_gen_protocol%2Clibrary_layout%2Clibrary_max_fragment_size%2Clibrary_min_fragment_size%2Clibrary_name%2Clibrary_pcr_isolation_protocol%2Clibrary_prep_date%2Clibrary_prep_date_format%2Clibrary_prep_latitude%2Clibrary_prep_location%2Clibrary_prep_longitude%2Clibrary_selection%2Clibrary_source%2Clibrary_strategy%2Clocal_environmental_context%2Clocation%2Clocation_end%2Clocation_start%2Clon%2Cmarine_region%2Cmating_type%2Cncbi_reporting_standard%2Cnominal_length%2Cnominal_sdev%2Cpcr_isolation_protocol%2Cph%2Cproject_name%2Cprotocol_label%2Cread_count%2Cread_strand%2Crestriction_enzyme%2Crestriction_enzyme_target_sequence%2Crestriction_site%2Crna_integrity_num%2Crna_prep_3_protocol%2Crna_prep_5_protocol%2Crna_purity_230_ratio%2Crna_purity_280_ratio%2Crt_prep_protocol%2Crun_alias%2Crun_date%2Csalinity%2Csample_accession%2Csample_alias%2Csample_capture_status%2Csample_collection%2Csample_description%2Csample_material%2Csample_prep_interval%2Csample_prep_interval_units%2Csample_storage%2Csample_storage_processing%2Csample_title%2Csampling_campaign%2Csampling_platform%2Csampling_site%2Cscientific_name%2Csecondary_project%2Csecondary_sample_accession%2Csecondary_study_accession%2Csequencing_date%2Csequencing_date_format%2Csequencing_location%2Csequencing_longitude%2Csequencing_method%2Csequencing_primer_catalog%2Csequencing_primer_lot%2Csequencing_primer_provider%2Cserotype%2Cserovar%2Csex%2Cspecimen_voucher%2Csra_aspera%2Csra_bytes%2Csra_ftp%2Csra_galaxy%2Csra_md5%2Cstatus%2Cstrain%2Cstudy_accession%2Cstudy_alias%2Cstudy_title%2Csub_species%2Csub_strain%2Csubmission_accession%2Csubmission_tool%2Csubmitted_aspera%2Csubmitted_bytes%2Csubmitted_format%2Csubmitted_ftp%2Csubmitted_galaxy%2Csubmitted_host_sex%2Csubmitted_md5%2Csubmitted_read_type%2Ctag%2Ctarget_gene%2Ctaxonomic_classification%2Ctaxonomic_identity_marker%2Ctemperature%2Ctissue_lib%2Ctissue_type%2Ctransposase_protocol%2Cvariety"
        data = f'result={result}&query={query}&fields={fields}&format={format}'

        r = requests.post(
            "https://www.ebi.ac.uk/ena/portal/api/search", headers=headers, data=data
        )
        print("Data obtained, now formatting...")
        df = pd.read_table(io.StringIO(r.text))
        df=df.dropna(how="all", axis=1)

        emptycols = []
        for col in df.columns:
            if len(df[col].unique()) == 1 and df[col].unique()[0] == "":
                emptycols.append(col)
        df = df.drop(emptycols, axis=1)
        print(f"All ENA reads for taxon {self.taxon} obtained (n={len(df)})")
        if n_samples=='all':
            return(df)
        else:
            return df.sample(n_samples)

    def get_biosample_entry_from_ena_id(self, id):
        biosamples_server = "https://www.ebi.ac.uk/biosamples/samples/"
        r = requests.get(
            biosamples_server + id, headers={"Content-Type": "application/json"}
        )
        if r.ok:
            decoded = r.json()
            decoded = pd.json_normalize(decoded).transpose().to_dict()[0]
            decoded = pd.Series(decoded)
            decoded.index = [i.replace("characteristics.", "") for i in decoded.index]
            decoded=decoded.dropna()
            decoded=decoded.rename(index={'accession':'sample_accession'})
            return decoded
        else:
            return {"sample_accession":id}

    def clean_biosamples_entry(self, row):
        acceptables_list = [str, bool, float, int, np.int64, np.bool_, np.float64]
        newrow = pd.Series()
        for index, cell in row.items():
            if index in newrow.index:
                index=index+"_2"
            cell_type = type(cell)
            if isinstance(cell, list):
                if len(cell) == 1:
                    try:
                        newrow[index]=cell[0]["text"]
                    except:
                        pass
            elif cell_type in acceptables_list:
                newrow[index]=cell
            else:
                #print(cell_type)
                raise TypeError(f"Cell type: {cell_type}")

        return newrow


    def build_biosamples_db(self, input, thread_number=4):
        if isinstance(input, list):
            ids = set(input)
        elif isinstance(input, pd.DataFrame):
            ids=set(input['sample_accession'])
        else:
            raise Exception("Wrong format for biosamples IDs: needs list or pandas Dataframe with column sample_accession")

        def biosample_with_progress(x):
            biosample = self.get_biosample_entry_from_ena_id(x)
            biosample_clean = self.clean_biosamples_entry(biosample)
            pbar.update(1)
            return biosample_clean

        print("Building biosamples db from ENA samples...")
        pbar = tqdm(total=len(ids))
        pool = ThreadPool(thread_number)
        results = pool.map(biosample_with_progress, ids)
        pool.close()
        pool.join()
        biosamples_db = pd.concat(results, axis=1).transpose()
        print(f"Biosamples db completed (n={len(biosamples_db)})")

        return biosamples_db

    def split_df_in_rows(self, df):
        rows_list = [row for _, row in df.iterrows()]
        return rows_list

    def build_study_xrefs_db(self, input, thread_number=4):
        if isinstance(input, list):
            study_ids = set(input)
        elif isinstance(input, pd.DataFrame):
            study_ids=set(input['study_accession'])
        else:
            raise Exception("Wrong format for Xrefs IDs: needs list or pandas Dataframe with column study_accession")

        pbar = tqdm(total=len(study_ids))
        study_ids=set(study_ids)
        def get_xref(study_id):
            pbar.update(1)
            headers = {"Content-Type": "application/json"}
            r = requests.get(
                "https://www.ebi.ac.uk/ena/xref/rest/json/search?accession="
                + study_id,
                headers=headers,
            )
            if r.ok:
                decoded = r.json()
                return(pd.DataFrame(decoded))
            else:
                return(np.nan)
        pool = ThreadPool(thread_number)
        xref_list = pool.map(get_xref, study_ids)
        pool.close()
        pool.join()
        xrefs_df=pd.concat(xref_list)
        return xrefs_df


    def remove_diacritics(self, input_string):
        input_string=str(input_string)
        return ''.join(
            char for char in unicodedata.normalize('NFD', input_string)
            if unicodedata.category(char) != 'Mn'  # Filter out combining characters
        )

    def strip_string_down_to_letters(self, input_string):
        stripped_string=self.remove_diacritics(input_string)
        #stripped_string=re.sub(r'\W+', '', stripped_string)
        stripped_string = re.sub(r'[^a-zA-Z0-9\s]+', '', stripped_string)
        stripped_string=stripped_string.lower()
        return(stripped_string)

    def get_fao_breeds_dataset(self):
        fao_dataset=pd.read_csv(f"https://raw.githubusercontent.com/matteobolner/livestock_breeds_database/main/species/{self.species}/all_breeds.csv")
        fao_dataset['breed_no_diacritics']=fao_dataset['Breed/Most common name'].dropna().apply(self.remove_diacritics)
        fao_dataset['breed_letters_only']=fao_dataset['Breed/Most common name'].apply(self.strip_string_down_to_letters)
        fao_dataset=fao_dataset[fao_dataset['Breed/Most common name'].apply(len)>=3]
        return(fao_dataset)

    def get_fao_breeds_stripped_list(self):
        fao_breeds_list=self.get_fao_breeds_dataset()['breed_letters_only'].unique().tolist()
        return(fao_breeds_list)

    def count_breed_occurrences_in_string(self, sample, breed_list, input_string):
        counts = {substring: len(re.findall(substring, input_string)) for substring in breed_list}
        counts = {k:v for k,v in counts.items() if v!=0}
        return(counts)

    def annotate_biosamples_with_breeds(self, breed_list='FAO', db='default', thread_number=4):
        if breed_list=='FAO':
            breed_list=self.get_fao_breeds_stripped_list()
        if isinstance(db, pd.DataFrame):
            biosamples=db
        elif (isinstance(db, string))&(db=='default'):
            biosamples=self.biosamples_db
        else:
            raise Exception("Unclear format for BioSamples DB")
        #biosamples=test
        datetimes_and_urls=[]
        counter=0
        for col in biosamples.columns:
            counter+=1
            if col.startswith("_links"):
                datetimes_and_urls.append(col)
            if biosamples[col].dtype == 'object':
                try:
                    pd.to_datetime(biosamples[col], format='mixed')
                    datetimes_and_urls.append(col)
                except:
                    pass
        biosamples=biosamples.drop(columns=datetimes_and_urls)
        biosamples=biosamples.drop(columns=['INSDC center name'])
        biosamples=biosamples.dropna(how='all', axis=1)
        biosamples=biosamples.dropna(how='all', axis=0)
        biosamples_stripped=biosamples.set_index('sample_accession').apply(lambda row: row.dropna().astype(str).str.cat(sep=' '), axis=1).apply(self.strip_string_down_to_letters)

        def annotate_sample_breed(sample_stripped_info):
            pbar.update(1)
            id=sample_stripped_info[0]
            info=sample_stripped_info[1]
            counts={id:self.count_breed_occurrences_in_string(sample=id, breed_list=breed_list, input_string=info)}
            return(counts)
        pbar = tqdm(total=len(biosamples))
        pool = ThreadPool(thread_number)
        annotated_samples = pool.map(annotate_sample_breed, list(biosamples_stripped.items()))
        pool.close()
        pool.join()
        annotated_samples_joined={}
        for i in annotated_samples:
            annotated_samples_joined.update(i)
        biosamples['breed_found_counts']=biosamples['sample_accession'].apply(lambda x:annotated_samples_joined[x])
        sample_breeds=biosamples[['sample_accession','organism','breed','breed_found_counts']]

        def remove_duplicates(input_dict):
            # Sort the dictionary by values in descending order
            sorted_dict_by_value = {k: v for k, v in sorted(input_dict.items(), key=lambda item: item[1], reverse=True)}
            sorted_dict_by_key = {}
            for k in sorted(sorted_dict_by_value, key=len, reverse=True):
                sorted_dict_by_key[k] = sorted_dict_by_value[k]

            # Create a list of keys from the sorted dictionary
            dict_keys = list(sorted_dict_by_key.keys())

            # Loop through the keys to identify and remove duplicates or overlapping keys
            for index, key in enumerate(dict_keys):
                for index2, key2 in enumerate(dict_keys[index+1:]):
                    if " ".join(sorted(key2.split())) in " ".join(sorted(key.split())):
                        try:
                            del sorted_dict_by_key[key2]
                        except KeyError:
                            pass
                    #elif sorted(key.split()) == sorted(key2.split()):
                    #    del input_dict[key2]
                    else:
                        if set(key2.split()).issubset(set(key.split())):
                            try:
                                del sorted_dict_by_key[key2]
                            except KeyError:
                                pass
            return sorted_dict_by_key

        def remove_dict_if_only_1_key(input_dict):
            if len(input_dict)==1:
                return(list(input_dict)[0])
            else:
                return(input_dict)

        sample_breeds['breed_found_counts']=sample_breeds['breed_found_counts'].apply(remove_duplicates)
        sample_breeds['breed_found_counts']=sample_breeds['breed_found_counts'].apply(remove_dict_if_only_1_key)
        sample_breeds['inferred_breed']=sample_breeds['breed'].where(~sample_breeds['breed'].isna(), np.nan)
        return(sample_breeds)
