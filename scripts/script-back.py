import pandas as pd
import geopandas as gpd
import os

os.chdir(
    os.path.join(
        os.path.expanduser('~'),
         'Documents', 'Manuscripts', 'GermplasmCurationTools', 'scripts')
         )

ca_popgen = pd.read_csv(
    os.path.join('..', 'data', 'CaliPopGen_dataset_1_population_genetic_diversity_TSV.tsv'),
    sep = '\t'
)

# select only core columns, no need for lifeform characteristics
ca_popgen=ca_popgen.loc[:, 'CitationID':'HabitatType'] 

# only interested in seed plants
ca_popgen= ca_popgen[ ca_popgen['Phylum']=='Magnoliophyta']

# and we focus on terrestrial habitat
ca_popgen= ca_popgen[ ca_popgen['HabitatType']=='terrestrial']


#################
# restrict to populations with coordinate uncertatinty less than 20 KM ??? (verify distance measure, assuming km)
ca_popgen[ca_popgen['CoordError']>20].CitationFull.unique()

coordErrors=ca_popgen[ca_popgen['CoordError']>20].CitationID.unique() # grab study IDs
# and remove all based on them. 
ca_popgen = ca_popgen[~ca_popgen['CitationID'].isin(coordErrors)]

del coordErrors
################


#################
## determine how many populations were sampled in the study
ca_popgen['PopulationsCt'] = ca_popgen.groupby(
    ['CitationID', 'SpeciesID']
)[ca_popgen.columns[0]].transform('count')


## some of these just seem like they were outgroups for a focal taxon. 
ca_popgen[ca_popgen['PopulationsCt']<4].CitationFull.unique()

outgroups = ca_popgen[ca_popgen['PopulationsCt']<=4].CitationID.unique()

ca_popgen = ca_popgen[
      # find studies with 4 or more populations for the species. 
      # OR find the species in studies with 4 or more populations
    (~ca_popgen['CitationID'].isin(outgroups)) | 
    (
        (ca_popgen['CitationID'].isin(outgroups)) &
        (ca_popgen['PopulationsCt'] >= 4)
    )
] ## alternatively, and straight to the point, ca_popgen[ca_popgen['PopulationsCt'] >= 4]

del outgroups

#######################
## determine if species have repeat treatments across studies

possible_dupes = (
    ca_popgen[['CitationID', 'SpeciesID']]
    .drop_duplicates()
    .groupby(['SpeciesID'])
    .size()
    .reset_index(name='n')
) 

# compare the studies 
possible_dupes = pd.merge(
    possible_dupes[possible_dupes['n']>1], 
    ca_popgen[['CitationFull', 'CitationID','SpeciesID']].drop_duplicates(),
    on='SpeciesID',
    how='inner')

## most of these are just on rye and non-native species anyways, not rare taxa. 
# we will just drop them all. 

ca_popgen=ca_popgen[~ca_popgen['CitationID'].isin(possible_dupes['CitationID'].tolist())]

del possible_dupes

#######################
## remove noxious / invasive species from data set



#########################
## select variables to remove 
cols2remove = ['EffectivePopSize', 'PercentPolyLoci', 'AllelesPerLocus', 'InbreedingCoefType', 'InbreedingCoefValue', 'HaploDiv', 'NucDiversity']
ca_popgen.drop(cols2remove, axis=1, inplace=True)

del cols2remove

## create a spatial object 
gdf = gpd.GeoDataFrame(
    ca_popgen, geometry=gpd.points_from_xy(ca_popgen.LongitudeDD, ca_popgen.LatitudeDD), crs="EPSG:4326")

gdf.columns
