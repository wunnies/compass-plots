# -*- coding: utf-8 -*-
#Created on Mon Aug  2 20:01:42 2021

#@author: Wunna
def cohens_d(x, y):
    pooled_std = np.sqrt(((len(x)-1) * np.var(x, ddof=1) 
                          + (len(y)-1) * np.var(y, ddof=1)) / 
                             (len(x) + len(y) - 2))
    return (np.mean(x) - np.mean(y)) / pooled_std
    


def wilcoxon_test(consistencies_matrix, group_A_cells, group_B_cells):
    """
        Performs an unpaired wilcoxon test (or mann-whitney U test) for each reaction between group_A and group_B
    """
    #per reaction/meta-reaction, conduct wilcoxon test between group_A and group_B
    group_A = consistencies_matrix.loc[:,group_A_cells]
    group_B = consistencies_matrix.loc[:,group_B_cells]
    results = pd.DataFrame(index = consistencies_matrix.index, columns = ['wilcox_stat', 'wilcox_pval', 'cohens_d'], dtype='float64')
    warning_count=0
    for rxn in consistencies_matrix.index:
        try:
            A, B = group_A.loc[rxn].to_numpy().ravel(), group_B.loc[rxn].to_numpy().ravel()
            stat, pval = mannwhitneyu(A, B, alternative='two-sided')
            c_d = cohens_d(A, B)
            results.loc[rxn, ['wilcox_stat', 'wilcox_pval', 'cohens_d']] = stat, pval, c_d
        except:
            print("Exception occurred for reaction: " + str(rxn) + ". Not included in analysis")
            warning_count=warning_count+1
    results['adjusted_pval'] = np.array(multipletests(results['wilcox_pval'], method='fdr_bh')[1], dtype='float64')
    print("Warnings: "+ str(warning_count) + " reactions out of "+ str(len(consistencies_matrix.index)) + " not included")
    return results

def get_reaction_consistencies(compass_reaction_scores, min_range=1e-3):
    df = -np.log(compass_reaction_scores + 1)
    df = df[df.max(axis=1) - df.min(axis=1) >= min_range]
    df = df - df.min().min()
    return df

def plot_overview(data, celltype_a_name, celltype_b_name, interactive=True, save=None, xlims=None):
    plt.figure(figsize=(12,12))
    axs = plt.gca()
    #Sorts the reactions for plotting
    d = data[data['adjusted_pval'] < 0.1].groupby('subsystem')['cohens_d'].median().abs()
    if xlims:
        plt.xlim(xlims)   # set the xlim to left, right

    axs.scatter(d[d.argsort], d[d.argsort].index, alpha=0)
    color = data['cohens_d'].map(lambda x: 'r' if x >= 0 else 'b')
    alpha = data['adjusted_pval'].map(lambda x: 1.0 if x < 0.1 else 0.25)
    axs.scatter(data['cohens_d'], data['subsystem'], c=color, alpha=alpha ) # dots
    
    # static labels for a .png
    # for row in data.iterrows():
    #     rxn=row[1]
    #     if abs(rxn["cohens_d"]) > 0.5 and rxn["adjusted_pval"] < 0.1:
    #             offset = (-10, 20)
    #             #annot=str(rxn["associated_genes"]) + "\n" + rxn["reaction_name"]
    #             annot=str(rxn["associated_genes"])
    #             axs.annotate(annot, (rxn['cohens_d'],rxn['subsystem']), xytext = offset, 
    #                          textcoords='offset pixels', arrowprops={'arrowstyle':"-"})
    
    
    # big axis arrows at bottom of graph
    axs.annotate('', xy=(0.5, -0.08), xycoords='axes fraction', xytext=(0, -0.08), 
                arrowprops=dict(arrowstyle="<-", color='#348C73', linewidth=4))
    axs.annotate(celltype_a_name, xy=(0.60, -0.12), xycoords='axes fraction', fontsize=16)
    axs.annotate('', xy=(0.5, -0.08), xycoords='axes fraction', xytext=(1, -0.08), 
                arrowprops=dict(arrowstyle="<-", color='#E92E87', linewidth=4))
    axs.annotate(celltype_b_name, xy=(0.10, -0.12), xycoords='axes fraction', fontsize=16)
    axs.set_xlabel("Cohen's d")
    
    # add interactive labels for ipynb
    if (interactive==True):
        mplcursors.cursor(multiple = True).connect(
            "add", lambda sel: sel.annotation.set_text(     
                  str("Reaction: " + str(data["reaction_name"][sel.target.index]) + "\n" + "Genes: " +  str(data["associated_genes"][sel.target.index]))
        ))
    if save:
        plt.savefig( save + ".png", dpi=300, bbox_inches='tight')

### metareaction analysis

def plot_differential_scores(data, title, celltype_a_name, celltype_b_name, label_limits = (0.5, 1.0), interactive=False):
    plt.figure(figsize=(10,10))
    axs = plt.gca()
    color = data['cohens_d'].map(lambda x: 'r' if x >= 0 else 'b')
    axs.scatter(data['cohens_d'], -np.log10(data['adjusted_pval']), c=color)
    axs.set_xlabel("Cohen's d", fontsize=16)
    axs.set_ylabel("-log10 (Wilcoxon-adjusted p)", fontsize=16)
    #Everything after this should be tweaked depending on your application
    axs.set_xlim(-2.2, 2.2)
    axs.axvline(0, dashes=(3,3), c='black')
    axs.axhline(1, dashes=(3,3), c='black')
    axs.set_title(title, fontdict={'fontsize':20})
    axs.annotate('', xy=(0.5, -0.08), xycoords='axes fraction', xytext=(0, -0.08), 
            arrowprops=dict(arrowstyle="<-", color='#348C73', linewidth=4))
    axs.annotate(celltype_a_name, xy=(0.75, -0.12), xycoords='axes fraction', fontsize=16)
    axs.annotate('', xy=(0.5, -0.08), xycoords='axes fraction', xytext=(1, -0.08), 
            arrowprops=dict(arrowstyle="<-", color='#E92E87', linewidth=4))
    axs.annotate(celltype_b_name, xy=(0.25, -0.12), xycoords='axes fraction', fontsize=16)
    if interactive == True: # cursor labels for ipynb
        mplcursors.cursor(multiple = True).connect(
            "add", lambda sel: sel.annotation.set_text(
                "Reaction: " + str(data["reaction_name"][sel.target.index]) + "\n" + "Genes: " +  str(data["associated_genes"][sel.target.index]))
            )
    
    else: # static labels (for saving a .png)
        for row in data.iterrows():
            rxn=row[1]
            log_pval=-np.log10(rxn["adjusted_pval"])
            if abs(rxn["cohens_d"]) > label_limits[0] and log_pval > label_limits[1]:
                 offset = (-10, 20)
                 annot=str(rxn["associated_genes"]) + "\n" + rxn["reaction_name"]
                 #annot=str(rxn["associated_genes"])
                 axs.annotate(annot, (rxn['cohens_d'],log_pval), xytext = offset, 
                              textcoords='offset pixels', arrowprops={'arrowstyle':"-"})

        
def get_metareactions(reactions, height=0.02):
    """
        Returns an array of metareaction labels for each reaction
        Index k in the returned array has the metareaction label for reaction k.
    """
    #pairwise_reaction_correlations = reactions.T.corr(method='spearman') #Pandas method here is orders of magnitude slower
    pairwise_reaction_correlations = np.corrcoef(reactions.rank(axis=1))
    #Unfortunately due to floating point issues, these matrices are not always perfectly symmetric and the diagonal may be slightly off from 1
    pairwise_reaction_correlations[np.arange(reactions.shape[0]), np.arange(reactions.shape[0])] = 1.0
    pairwise_reaction_correlations = (pairwise_reaction_correlations + pairwise_reaction_correlations.T)/2
    assert(np.all(pairwise_reaction_correlations == pairwise_reaction_correlations.T))

    Z = hcluster.complete(squareform(1 - pairwise_reaction_correlations))
    return hcluster.fcluster(Z, height, criterion='distance')

import pandas as pd
import numpy as np
from scipy.stats import wilcoxon, mannwhitneyu, ranksums
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hcluster
from scipy.spatial.distance import squareform
import mplcursors

reactions = pd.read_csv("reactions_OC.tsv", sep="\t", index_col = 0)
cell_metadata = pd.read_csv("seurat_metadata_osteomorph.csv", index_col = 0)
cell_metadata.loc[cell_metadata["type"]=="Other", "type"] = "Other_single_positive"


celltype_a_name="Osteomorph"
celltype_b_name="Osteoclast"
celltype_a = cell_metadata[cell_metadata['type'] == celltype_a_name]["Sample"].values
celltype_b = cell_metadata[cell_metadata['type'] == celltype_b_name]["Sample"].values

reaction_metadata = pd.read_csv("reaction_metadata.txt", index_col = 0)

reaction_consistencies = get_reaction_consistencies(reactions)
reaction_consistencies.T.to_csv("consistencies_for_nmf.csv") # export to R for nmf

wilcox_results = wilcoxon_test(reaction_consistencies, celltype_a, celltype_b)
wilcox_results['metadata_r_id'] = ""
for r in wilcox_results.index:
    if r in reaction_metadata.index:
        wilcox_results.loc[r, 'metadata_r_id'] = r
    elif r[:-4] in reaction_metadata.index:
        wilcox_results.loc[r, 'metadata_r_id'] = r[:-4]
    else:
        print("Should not occur")
        
W = wilcox_results.merge(reaction_metadata, how='left', 
                         left_on='metadata_r_id', right_index=True, validate='m:1')
W = W[W['confidence'].isin([0,4])]
W = W[~W['EC_number'].isna()]

#### Filtering out citric acid cycle - Don't understand reason they gave "its outside mitochondria"
#W.loc[(W['formula'].map(lambda x: '[m]' not in x)) & (W['subsystem'] == "Citric acid cycle"), 'subsystem'] = 'Other'

data=W
#data = W[~W['subsystem'].isin(["Miscellaneous", "Unassigned"])]

##### Filtering certain metabolic pathways out - why?
#data = data[~data['subsystem'].map(lambda x: "Transport" in x or "Exchange" in x or x == "Other")]

items, counts = np.unique(data['subsystem'], return_counts=True)
items = [items[i] for i in range(len(items)) if counts[i] > 5] #filter(n() > 5) %>%
data = data[data['subsystem'].isin(items)]

plot_overview(data, celltype_a_name, celltype_b_name, interactive=True, save=None)

reactions[reactions <= 1e-4] = 0
reactions = reactions[np.all(reactions != 0, axis=1)]
reactions = reactions[reactions.max(axis=1) - reactions.min(axis=1) != 0]

meta_rxns_map = get_metareactions(reactions)
meta_rxns = reactions.join(pd.DataFrame(meta_rxns_map, columns=["meta_rxn_id"], index = reactions.index)).groupby("meta_rxn_id").mean()
meta_rxn_consistencies = get_reaction_consistencies(meta_rxns)
wilcox_meta_rxn_results = wilcoxon_test(meta_rxn_consistencies, celltype_a, celltype_b)

wilcox_meta_rxn_expanded = pd.DataFrame(index=reactions.index, columns=wilcox_meta_rxn_results.columns)
for i in range(len(wilcox_meta_rxn_expanded.index)):
    if (meta_rxns_map[i] in wilcox_meta_rxn_results.index):
        wilcox_meta_rxn_expanded.loc[wilcox_meta_rxn_expanded.index[i]] = wilcox_meta_rxn_results.loc[meta_rxns_map[i]]
wilcox_meta_rxn_expanded = wilcox_meta_rxn_expanded.dropna().astype('float64')

wilcox_meta_rxn_expanded['metadata_r_id'] = ""
for r in wilcox_meta_rxn_expanded.index:
    if r in reaction_metadata.index:
        wilcox_meta_rxn_expanded.loc[r, 'metadata_r_id'] = r
    elif r[:-4] in reaction_metadata.index:
        wilcox_meta_rxn_expanded.loc[r, 'metadata_r_id'] = r[:-4]
    else:
        print("Should not occur")
        
        
W = wilcox_meta_rxn_expanded.merge(reaction_metadata, how='left', 
                        left_on='metadata_r_id', right_index=True, validate='m:1')
W = W[W['confidence'].isin([0,4])]
W = W[~W['EC_number'].isna()]
#W.loc[(W['formula'].map(lambda x: '[m]' not in x)) & (W['subsystem'] == "Citric acid cycle"), 'subsystem'] = 'Other'

subsystem_name="Citric acid cycle"
subsystem_data = W[W['subsystem'] == subsystem_name]
plot_differential_scores(subsystem_data, title=subsystem_name, celltype_a_name=celltype_a_name, celltype_b_name=celltype_b_name, interactive=True,)

#system_data = W[W['subsystem'] == "Nucleotide interconversion"]
#plot_differential_scores(subsystem_data, title='Nucleotide interconversion', label_limits=(0.5, 1.0), interactive_labels=True)

# Lower the cohen's D, the more associated with osteoclasts. Left=osteoclasts, right=osteomorphs
