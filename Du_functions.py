import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multitest import fdrcorrection


def Hawk_smash(List):
    '''Flattern lists within a list '''
    return [item for sublist in List for item in sublist]

def ANOVA_test(meta_data,df_otu_counts,Group_factor,Group_list,P_cut):
    '''Conducting one way ANOVA test to find significantly enriched OTUs in each group
       dict_enriched_OTU_group: enriched OTUs in each group
       df_ANOVA_results_group: F, p, and p_adj for each enriched OTU within each group
    '''
    
    meta_data_filter = meta_data.loc[df_otu_counts.columns]
    dict_group = dict(zip(Group_list,[list(meta_data_filter[meta_data_filter[Group_factor] == x].index) for x in Group_list]))

    ANOVA_results = {}
    for group in Group_list:
        Others = list(set(df_otu_counts.columns) - set(dict_group[group]))
        ANOVA_results[group] = [f_oneway(list(df_otu_counts[dict_group[group]].loc[x]), list(df_otu_counts[Others].loc[x]))[0:2] for x in df_otu_counts.index] 
            
    df_ANOVA_results_group = {}
    for group in Group_list:
        df_ANOVA_results_group[group] = pd.DataFrame(ANOVA_results[group],index = df_otu_counts.index, columns=['F_ratio','p_value'])
        df_ANOVA_results_group[group]['P_adj'] = list(fdrcorrection(df_ANOVA_results_group[group]['p_value'])[1])
        df_ANOVA_results_group[group] = df_ANOVA_results_group[group].sort_values(by = 'F_ratio',ascending = False)        
    

    # Target the enriched ones for each group of interest
    Enrich_ratio = pd.DataFrame(index = df_otu_counts.index)
    for group in Group_list:
        Other_group = list(set(df_otu_counts.columns) - set(dict_group[group]))
        Enrich_ratio[group] = df_otu_counts[dict_group[group]].mean(axis = 1) - df_otu_counts[Other_group].mean(axis = 1)
    
    dict_enriched_OTU_group = {}
    for group in Group_list:
        df_temp = df_ANOVA_results_group[group]
        pos_list = set(df_temp.index).intersection(Enrich_ratio[Enrich_ratio[group]>0][group].index)
        dict_enriched_OTU_group[group] = df_temp.loc[pos_list][df_temp.P_adj<P_cut].index
    
    return dict_enriched_OTU_group,df_ANOVA_results_group

def Count_OTU_by_reference(df_candidate,df_Annotation_reference,group):

    '''df_candidate: index-OTU, col - D(1)- D(n) taxa levels
       df_Annotation_reference: Reference Taxa structure with custom defined group names
       dict_OTU_reference: OTUs that belong to the defined group in reference
       dict_OTU_reference_counts: OTU counts for each group'''
    dict_OTU_reference = {}
    Index_pool = list(df_Annotation_reference.index)
    for col in list(reversed(df_Annotation_reference.columns)):
        Pick_list = df_Annotation_reference.loc[Index_pool][col].dropna().index
        Index_pool = list(set(Index_pool)-set(Pick_list))
        for Pick_index in Pick_list:
            index_temp = list(df_candidate[df_candidate[col] == df_Annotation_reference.loc[Pick_index][col]].index)
            dict_OTU_reference[Pick_index] = list(set(index_temp)-set(Hawk_smash(dict_OTU_reference.values())))
    df_OTU_reference_counts = pd.DataFrame({group:[len(dict_OTU_reference[name]) for name in dict_OTU_reference.keys()]},index = dict_OTU_reference.keys())        
    return dict_OTU_reference,df_OTU_reference_counts