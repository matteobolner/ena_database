import pandas as pd
import pycountry
import pycountry_convert as pc

def get_country_code_from_country_name(country_name):
    if country_name!=country_name:
        return(None)
    country = pycountry.countries.get(name=country_name)
    if country is None:
        return None
    return country.alpha_2

def get_country_name_from_country_code(country_code):
    if country_code!=country_code:
        return(None)
    country = pycountry.countries.get(alpha_2=country_code)
    if country is None:
        return None
    return country.name

def get_continent(country_name_or_code):
    if country_name_or_code != country_name_or_code:
        return None
    try:
        country_code = pycountry.countries.get(name=country_name_or_code).alpha_2
    except:
        country_code = country_name_or_code
    try:
        continent_code = pc.country_alpha2_to_continent_code(country_code)
    except:
        return None
    return continent_code


df=pd.read_csv("/home/pelmo/work/workspace/Metadata/Sus_scrofa_breeds.csv")

df['Country code'] = df['Country code'].fillna(df['Country'].astype(str).map(lambda x: get_country_code_from_country_name(x)))

df['Country'] = df['Country'].fillna(df['Country code'].astype(str).map(lambda x: get_country_name_from_country_code(x)))
df['Continent'] = df['Continent'].fillna(df['Continent code'].astype(str).map(lambda x: get_continent(x)))

df
