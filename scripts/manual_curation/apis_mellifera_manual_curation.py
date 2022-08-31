import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pycountry_convert as pc


def move_column_inplace(df, col, pos):
    col = df.pop(col)
    df.insert(pos, col.name, col)


df = pd.read_csv("ena_databases/apis/genomic_data/manually_curated_db.csv")
df=df.set_index("Unnamed: 0")
df.index.name=None

#other_df = pd.read_csv("ena_databases/apis/merged_databases.csv")
#other_df=other_df[['run_accession','base_count','library_strategy']]
#df=df.merge(other_df)
#df.to_csv("ena_databases/apis/genomic_data/manually_curated_db_2.csv")
df=df[df['library_strategy'].isin(['WGS','WGA','OTHER'])]

df['notes']=df['notes'].fillna("")

terms_to_remove = ['crispr','methylation','metagenome','rna','rna-seq', 'RNASeq', 'Rna-seq']
indexes_to_remove = []
for term in terms_to_remove:
    indexes_to_remove.append(df[df['notes'].str.contains(term)].index.tolist())


indexes_to_remove=sum(indexes_to_remove, [])

df = df.drop(indexes_to_remove)

df=df.drop(columns=['country_lookup','state_by_coordinates','country_by_coordinates','sex_ena_db','sex_bs_db','sex','ploidy','sub_species'])

non_mellifera = ['dorsata','cerana','florea','andreniformis','laboriosa']

species_list = []
subspecies_list = []
for index,row in df.iterrows():
    if row['subspecies_manually_curated']=='cerana japonica':
        species_list.append('cerana')
        subspecies_list.append('japonica')

    elif row['subspecies_manually_curated'] in non_mellifera:
        species_list.append(row['subspecies_manually_curated'])
        subspecies_list.append(np.nan)

    else:
        species_list.append("mellifera")
        subspecies_list.append(row['subspecies_manually_curated'])

df['apis_species']=species_list
df['apis_subspecies']=subspecies_list

df=df.drop(columns=['species','scientific_name','subspecies_manually_curated'])


move_column_inplace(df, 'accession',0)
move_column_inplace(df, 'study_accession',1)
move_column_inplace(df, 'run_accession',2)

move_column_inplace(df, 'pmid/doi',2)
move_column_inplace(df, 'apis_species',3)
move_column_inplace(df, 'apis_subspecies',4)
move_column_inplace(df, 'geographic_location_manually_curated',5)
move_column_inplace(df, 'sex_manually_curated',6)
move_column_inplace(df, 'class_manually_curated',7)
move_column_inplace(df, 'sampling_year_manually_curated',8)

move_column_inplace(df, 'notes',9)
move_column_inplace(df, 'read_count',10)
move_column_inplace(df, 'base_count',11)
move_column_inplace(df, 'depth',12)
move_column_inplace(df, 'sample_average_depth',13)
move_column_inplace(df, 'library_strategy',14)
move_column_inplace(df, 'library_source',15)
move_column_inplace(df, 'library_layout',16)




df=df.drop(columns=['collection_date'])

pooled_verdict=[]
for index,row in df.iterrows():
    if "pool" in row['notes']:
        pooled_verdict.append('pooled')
    else:
        if 'pooled' in row['pooled_or_not']:
            pooled_verdict.append('pooled')
        else:
            pooled_verdict.append('unknown')

df['pooled_samples']=pooled_verdict

df=df.drop(columns=['pooled_or_not'])

move_column_inplace(df, 'pooled_samples',9)

move_column_inplace(df, 'coverage',16)


#df.to_csv("ena_databases/apis/genomic_data/manually_curated_db_fixed.csv")


mellifera=df[df['apis_species']=='mellifera'].copy()

len(mellifera.groupby('accession'))

mellifera['apis_subspecies']=mellifera['apis_subspecies'].fillna('unknown')

len(mellifera[mellifera['apis_subspecies']!='unknown'].groupby('accession'))
len(mellifera.groupby('accession'))

mellifera['pmid/doi']=mellifera['pmid/doi'].fillna("unknown")
len(mellifera[mellifera['pmid/doi']=='unknown'].groupby('accession'))
len(mellifera[mellifera['pmid/doi']!='unknown'].groupby('accession'))

mellifera['geographic_location_manually_curated']=mellifera['geographic_location_manually_curated'].fillna('Unknown')

len(mellifera[mellifera['geographic_location_manually_curated']!='Unknown'].groupby('accession'))
mellifera[mellifera['geographic_location_manually_curated']!='Unknown']
mellifera['geographic_location_manually_curated'].unique().tolist()
len(mellifera.groupby('accession'))


mellifera['sex_manually_curated']=mellifera['sex_manually_curated'].fillna("unknown")
len(mellifera[mellifera['sex_manually_curated']!='unknown'].groupby('accession'))
mellifera['class_manually_curated']=mellifera['class_manually_curated'].fillna("unknown")



plt.clf()
sns.set_theme()
sns.countplot()

sns.set(rc={'figure.figsize':(15,10)})


samplesdf = pd.DataFrame(index=df['accession'].unique(), columns=['species','subspecies','sex','geographic_location','pooled_samples'])


for name,group in df.groupby(by='accession'):
    templist = []
    templist.append(",".join(group['apis_species'].fillna("unknown").unique().tolist()))
    templist.append(",".join(group['apis_subspecies'].fillna("unknown").unique().tolist()))
    templist.append(",".join(group['sex_manually_curated'].fillna("unknown").unique().tolist()))
    templist.append(",".join(group['geographic_location_manually_curated'].fillna("unknown").unique().tolist()))
    templist.append(",".join(group['pooled_samples'].fillna("unknown").unique().tolist()))

    samplesdf.loc[name]=templist

samplesdf_mellifera = samplesdf[samplesdf['species']=='mellifera'].copy()

plt.clf()
ax = sns.countplot(x='species', data=samplesdf)
for item in ax.get_xticklabels():
    item.set_rotation(90)
ax.set_xlabel("Apis species")
ax.bar_label(ax.containers[0])
plt.savefig("ena_databases/apis/stats/species.png")

plt.clf()
ax = sns.countplot(data=samplesdf_mellifera, x='subspecies',order = samplesdf_mellifera['subspecies'].value_counts().index, palette='muted')
for item in ax.get_xticklabels():
    item.set_rotation(90)
ax.set_xlabel("Apis mellifera subspecies")
ax.bar_label(ax.containers[0])
plt.savefig("ena_databases/apis/stats/subspecies.png")

plt.clf()
ax = sns.countplot(data=samplesdf_mellifera, x='sex',order = samplesdf_mellifera['sex'].value_counts().index, palette='muted')
ax.set_xlabel("Sex")
ax.bar_label(ax.containers[0])
plt.savefig("ena_databases/apis/stats/sex.png")


country_code = pc.country_name_to_country_alpha2("Switzerland", cn_name_format="default")

countries = samplesdf_mellifera['geographic_location']
countries=countries.str.replace("Corsica","France")
countries=countries.str.replace("Scotland","GBR")

countries=countries.str.replace("armenia","Armenia")
countries=countries.str.replace("kenya","Kenya")
countries=countries.str.replace("France imported from Georgia","Georgia")
countries=countries.str.replace("Rodrigues","Mauritius")
countries=countries.str.replace("Reunion Island","Mauritius")
countries=countries.str.replace("One of two Baltic Sea Islands","Sweden")


country_codes = []
for country in countries:
    a = country.split(",")
    if len(a)>1:
        i=a[1].strip().replace("?","")
        country_codes.append(pc.country_name_to_country_alpha2(i, cn_name_format="default"))
    else:
        a=a[0]
        #print(a)
        if "Bern" in a:
            i='Switzerland'
        elif "Spain" in a:
            i='Spain'
        else:
            i=a.strip().replace("?","")
        if a=='unknown' or a=='Unknown':
            country_codes.append(np.nan)
        else:
            country_codes.append(pc.country_name_to_country_alpha2(i, cn_name_format="default"))

continents = []
for i in country_codes:
    if i==i:
        continent_name = pc.country_alpha2_to_continent_code(i)
        continents.append(pc.convert_continent_code_to_continent_name(continent_name))
    else:
        continents.append("Unknown")

samplesdf_mellifera['continent']=continents
plt.clf()
ax = sns.countplot(data=samplesdf_mellifera, x='continent',order = samplesdf_mellifera['continent'].value_counts().index, palette='muted')
ax.set_xlabel("Continent")
ax.bar_label(ax.containers[0])
plt.savefig("ena_databases/apis/stats/continent.png")


plt.clf()
ax = sns.countplot(data=samplesdf_mellifera, x='pooled_samples',order = samplesdf_mellifera['pooled_samples'].value_counts().index, palette='muted')
ax.set_xlabel("Pooled samples")
ax.bar_label(ax.containers[0])
plt.savefig("ena_databases/apis/stats/pooled.png")


mellifera=mellifera.merge(samplesdf_mellifera['continent'], left_on='accession', right_index=True)

move_column_inplace(mellifera, "continent", 6)

len(mellifera[mellifera['sample_average_depth']<5]['accession'].unique())
mellifera.to_csv("ena_databases/apis/genomic_data/mellifera_db.csv")
