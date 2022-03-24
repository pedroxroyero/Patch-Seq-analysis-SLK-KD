#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pylab as plt
import pandas as pn
from glmnet import glmnet
import scipy
import seaborn as sns; sns.set()
import sparseRRR

#%% Preprocessing of the data for sparse RRR analysis

# Load count and ephys matrices, gene information and sample table
counts = pn.read_csv('count_exons_introns_shSLK_vs_hrGFP_QC_Norm.csv')
genes = pn.read_csv('genes_list_shSLK_vs_hrGFP_QC.csv')
samples = pn.read_csv('sample_table_shSLK_vs_hrGFP_QC.csv')
ephys = pn.read_csv('Ephys_PSC+PasPro_QC_OnlyshSLK_hrGFP.csv')

# Set row names for gene and sample information tables
z = genes.copy()
z = z.set_index(z.iloc[:,0])

s = samples.copy()
s = s.set_index(s.iloc[:,0])

# Set row names and columns for count matrix
# Sort rows and columns according to gene and sample information tables
# Convert to numpy array
x = counts.transpose().copy()
x.columns = x.iloc[0,:]
x = x.iloc[1:,:]
x = x.set_index(x.index.astype(float))
x = x.reindex(s.index)
x = x.reindex(columns=z.index)
x = x.astype(np.float)
x = np.array(x)

# Set row names for ephys matrix
# Sort rows according to sample information table
y = ephys.copy()
y = y.set_index(y.iloc[:,0])
y = y.reindex(s.index)
y = y.iloc[:,1:]

# Drop non-altered ephys properties from ephys matrix                       
ephys_features = ['mIPSC frequency', 'Input resistance', 'Capacitance']
y = y.drop(['I.E.balance.amplitude', 'I.E.balance.frequency','EPSC.amplitude', 'EPSC.frequency', 'IPSC.amplitude'], axis=1) 

# Convert missing values to numpy nan in ephys matrix
# Convert to numpy array
y[y=='NA'] = np.nan
y = y.astype(np.float)
y = np.array(y)

# Remove samples with missing values in the ephys variables from the count and ephys matrices
remove_nan = ~np.isnan(y).any(axis=1)
y = y[remove_nan]
x = x[remove_nan]

# Select the most abundant genes in the count matrix
genesForRRR = np.isin(np.arange(0,x.shape[1]), (np.argsort(x.copy().sum(axis=0), axis=0)[::-1][:3005]))
x = x[:,genesForRRR]

# Standarize gene expression and ephys variables to zero mean and unit variance
x = x - np.mean(x, axis=0)
x = x / np.std(x, axis=0)

y = y - np.nanmean(y, axis=0)
y = y / np.nanstd(y, axis=0)

print('Shape of X:', x.shape, '\nShape of Y:', y.shape)

#%% Cross-validation for hyperparameter tuning

# Cross-validated R2 with different alpha values and number of selected genes
lambdas = np.concatenate((np.arange(.01,1.01,.1), np.arange(2,10)))
alphas = np.array([0.1, 0.25, .5, .75, 1])

cvresults = elastic_rrr_cv(x, y, rank=2, reps=1, folds=2, alphas=alphas, lambdas=lambdas)

r2, r2_relaxed, nonzeros, corrs, corrs_relaxed = cvresults
with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=RuntimeWarning)
    nonzeros = np.nanmean(nonzeros, axis=(0,1))
    r2_relaxed = np.nanmean(r2_relaxed, axis=(0,1))
    r2 = np.nanmean(r2, axis=(0,1))
    corrs_relaxed = np.nanmean(corrs_relaxed, axis=(0,1))

#
plt.plot(nonzeros, r2_relaxed, '.-', linewidth=1.5, markersize=4)
plt.ylim([-0.1,0.1])
plt.legend(['$\\alpha='+str(a)+'$' for a in [0.1, 0.25, .5, .75, 1]], frameon=True)
plt.show()

# Cross-validated R2 with different rank values and number of selected genes
lambdas = np.arange(.2,5,.1)
alphas = np.array([0.25])
ranks = np.arange(1, y.shape[1]+1)

cvresults_rank = {}
for r in ranks:
    cvresults_rank[r] = elastic_rrr_cv(x, y, rank=r, reps=1, folds=2, alphas=alphas, lambdas=lambdas)

#     
maxRank = len(cvresults_rank)
colA = np.array([100, 114, 176]) / 256
colB = np.array([100, 40, 50]) / 256
perf25 = np.zeros(maxRank+1)
for rank in range(1,maxRank+1):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        n = np.nanmean(cvresults_rank[rank][2], axis=(0,1))
        c = np.nanmean(cvresults_rank[rank][0], axis=(0,1))
        cr = np.nanmean(cvresults_rank[rank][1], axis=(0,1))
        f = scipy.interpolate.interp1d(n.squeeze(), cr.squeeze(), kind='linear')
        if n.squeeze().max() > 150 > n.squeeze().min():
            perf25[rank] = f(150)
        else:
            aux = []
            for valor in n.squeeze():
                aux.append(abs(150-valor))
                closer_valur = n.squeeze()[aux.index(min(aux))]
            perf25[rank] = f(closer_valur)        
    col = (colA * (rank-1)/(maxRank-1) + colB * (maxRank-rank)/(maxRank-1))
    plt.plot(n[:,0], cr[:,0], '.-', linewidth=1.5, markersize=4, label=rank)
plt.legend(frameon=True, loc='lower right')
plt.ylim([-0.06,0.1])

#
plt.plot(np.arange(1,maxRank+1), perf25[1:], 'k.-', linewidth=1.5, markersize=4, color='black')
plt.ylim([0,0.15])
sns.despine()
plt.tight_layout()

#%% Sparse RRR model plots

# Run the sparse RRR model. Hyperparameter values obtained from cross-validation analyses
w,v = relaxed_elastic_rrr(x, y, rank=2, lambdau=0.47, alpha=0.25)
print('\nGenes selected: {}'.format(np.sum(w[:,0]!=0)))
print(', '.join(z[genesForRRR].iloc[w[:,0]!=0,3].astype(str)))

color_samples = {'hrGFP':'black', 'shSLK':'red'}

bibiplot(x, y, w, v, 
                   titles = ['RNA expression', 'Electrophysiology'],
                   cellTypes = s.iloc[remove_nan,1], 
                   cellTypeColors = color_samples, 
                   YdimsNames = ephys_features,
                   XdimsNames = z.iloc[genesForRRR,3],                 
                   )

# For visualization purposes, the transcriptomic space was plotted labeling different relevant gene ontology terms separately
# Enriched ontology terms were obtained using ShinyGo and setting the false discovery rate to 0.05 (Ge et al., 2020). 

genes_NervousSystemDevelopment = ['Nav1', 'Adcy1', 'Kif2a', 'Prkn', 'Epha4', 'Srgap2', 'Sema6d', 'Elavl4', 'Sema3a', 'Foxp1', 'Cntn3', 'Cadm1', 'Herc1', 'Atp1a3', 'Sox5', 'Grm5', 'Nyap2', 'Slit3', 'Ntrk3', 'Agbl4', 'Ctnna2', 'Unc5d', 'Macrod2', 'Syndig1', 'Scn2a', 'Syne1']
genes_TranssynapticSignaling = ['Cacnb4', 'Adcy1', 'Prkn', 'Epha4', 'Elavl4', 'Cadm1', 'Cdh8', 'Plcl2', 'Chrm3', 'Grm5', 'Prkcb', 'Cacng3','Syne1', 'Sema6d', 'Sema3a', 'Slit3']
genes_cytoskeletalDynamics = ['Snx6', 'Klhl3', 'Kifa2', 'Prkn', 'Ptpn4', 'Vps16', 'Gnas', 'Mical2', 'Myrip', 'Ctnna3', 'Agbl4', 'Ctnna2', 'Syne1']
genes_IonChannels = ['Cacnb4', 'Trpm4', 'Atp1a3', 'Grm5', 'Kcnq3', 'Scn2a', 'Gabbr2', 'Chrm3', 'Cacng3', 'Gpr85', 'Atp11b']

Gene_terms = [genes_NervousSystemDevelopment, genes_TranssynapticSignaling, genes_cytoskeletalDynamics, genes_IonChannels]

for term in Gene_terms:
    Genes_plot = np.where(z[genesForRRR]['mgi_symbol'].isin(term))[0]
    bibiplot(x, y, w, v, 
                       titles = ['RNA expression', 'Electrophysiology'],
                       cellTypes = s.iloc[remove_nan,1], 
                       cellTypeColors = color_samples, 
                       YdimsNames = ephys_features,             
                       XdimsNames = z.iloc[genesForRRR,3],
                       YdimsToShow=None,                       
                       XdimsToShow = Genes_plot
                       )
    
#%% Cross-validation of single hyperparameters
lambdas = np.array([0.47])
alphas = np.array([0.25])

r2, r2_relaxed, nonzero, corrs, corrs_relaxed = elastic_rrr_cv(
    x, y, folds=2, reps=1, alphas=alphas, lambdas=lambdas, rank=2)

print(np.nanmean(r2_relaxed).squeeze())

