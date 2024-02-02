from mudata import MuData
import numpy as np
from tqdm import tqdm
import pandas as pd
import anndata
import scanpy
import csv
import sys
import loompy
import logging
import os
from matplotlib import pyplot as plt
import plotly.graph_objects as go

# function reads the loomfile downloaded from Tapestri portal
def read_tapestri_loom(filename):
    """
    Read data from MissionBio's formatted loom file.

    Parameters
    ----------
    filename : str
        Path to the loom file (.loom)

    Returns
    -------
    anndata.AnnData
        An anndata object with the following layers:
        adata.X: GATK calls
        adata.layers['e']: reads with evidence of mutation
        adata.layers['no_e']: reads without evidence of mutation
    """
    loom_file = loompy.connect(filename)

    variant_names, amplicon_names, chromosome, location = (
        loom_file.ra['id'], loom_file.ra['amplicon'], loom_file.ra['CHROM'], loom_file.ra['POS'])

    barcodes = [barcode.split('-')[0] for barcode in loom_file.ca['barcode']]
    adata = anndata.AnnData(np.transpose(loom_file[:, :]), dtype=loom_file[:, :].dtype)
    adata.layers['e'] = np.transpose(loom_file.layers['AD'][:, :])
    adata.layers['no_e'] = np.transpose(loom_file.layers['RO'][:, :])
    adata.var_names = variant_names
    adata.obs_names = barcodes
    adata.varm['amplicon'] = amplicon_names
    adata.varm['chrom'] = chromosome
    adata.varm['loc'] = location

    loom_file.close()

    return adata

def read_HyPR_file(HyPR_file):

    HyPR_array = np.loadtxt(HyPR_file, dtype=str, skiprows=1)
    barcodes = HyPR_array[:, 0]
    with open(HyPR_file, 'r') as f:
        first_line = f.readline().strip()
        values = first_line.split('\t')
        var_names = np.array(values, dtype=str)[:]
    X = np.asarray(HyPR_array[:, 1:], dtype=int)
    if np.size(var_names) != X.shape[1]: # sloppy I shouldn't need this
        var_names = var_names[1:]
    df = pd.DataFrame(X, columns=var_names)
    df.columns = df.columns.str.split('_').str[0]
    collapsed_df = df.groupby(axis=1, level=0).sum()
    collapsed_X = collapsed_df.to_numpy()
    collapsed_var_names = collapsed_df.columns.to_numpy()

    adata_hypr = anndata.AnnData(collapsed_X)
    adata_hypr.obs_names = barcodes

    try:
        adata_hypr.var_names = collapsed_var_names
    except:
        adata_hypr.var_names = collapsed_var_names[1:]

    return adata_hypr

# Function for finding the intersecting barcodes between the modalities
def find_intersecting(mdata):
    """
    Find intersecting barcodes and add 'intersecting' to obs with boolean values.

    Parameters
    ----------
    mdata : object
        The mudata object.

    Raises
    ------
    AssertionError
        If the length of mod_names is not equal to 2.
    """

    mod_names = list(mdata.mod.keys())
    assert len(mod_names) == 2, 'Function not implemented for mod_names with length different from 2'

    obs_1, obs_2 = mdata[mod_names[0]].obs_names, mdata[mod_names[1]].obs_names

    cmn_barcodes, idx_1, idx_2 = np.intersect1d(obs_1, obs_2, return_indices=True)

    # Initializing 'intersecting' columns to False
    mdata[mod_names[0]].obs['intersecting'] = False
    mdata[mod_names[1]].obs['intersecting'] = False

    # Updating 'intersecting' columns where True
    mdata[mod_names[0]].obs['intersecting'].iloc[idx_1] = True
    mdata[mod_names[1]].obs['intersecting'].iloc[idx_2] = True

    logger.info(f"Found {idx_1.size} intersecting barcodes")


def annotate_genotype(mdata, variants, window=4, obs_key='mutant_type',
                      genotype_key='tapestri', phenotype_key='HyPR',
                      ignore_bystanders=False, write_all_as=False,
                      ignore_mixed=False, only_type=None):

    adata = mdata.mod[genotype_key]

    # loop through the variants in the variant list
    for variant in variants:

        if obs_key in mdata.mod[phenotype_key].obs:
            pass
        else:
            mdata.mod[phenotype_key].obs[obs_key] = 'unannotated'
            mdata.mod[genotype_key].obs[obs_key] = 'unannotated'

        try:
            idx = adata.var_names.get_loc(variant)
        except KeyError:
            logger.info(str(variant) + ' was not found in the genotype modality anndata')
            continue

        # identify which cells have the variant
        bkg_idx, het_idx, hom_idx = [np.flatnonzero(adata.X[:, idx] == i) for i in (0, 1, 2)]

        if only_type == 'het':
            hom_idx = None
        if only_type == 'hom':
            het_idx = None

        if ignore_bystanders == True:

            if het_idx is not None:
                cmn_barcodes, idx_genotype, idx_phenotype = np.intersect1d(mdata.mod[genotype_key].obs_names[het_idx],
                                                                        mdata.mod[phenotype_key].obs_names, return_indices=True)

                if write_all_as == False:
                    mdata.mod[phenotype_key].obs[obs_key][idx_phenotype] = str(variant) + '_het'
                    mdata.mod[genotype_key].obs[obs_key][het_idx]= str(variant) + '_het'

                else:
                    mdata.mod[phenotype_key].obs[obs_key][idx_phenotype] = write_all_as

            if hom_idx is not None:
                cmn_barcodes, idx_genotype, idx_phenotype = np.intersect1d(mdata.mod[genotype_key].obs_names[hom_idx],
                                                                       mdata.mod[phenotype_key].obs_names, return_indices=True)
                if write_all_as == False:
                    mdata.mod[phenotype_key].obs[obs_key][idx_phenotype]= str(variant) + '_hom'
                    mdata.mod[genotype_key].obs[obs_key][hom_idx]= str(variant) + '_hom'
                else:
                    mdata.mod[phenotype_key].obs[obs_key][idx_phenotype] = write_all_as

        else:

          # identify indices of nearby SNPs that are within the bystander window
            chromosome, location = adata.varm['chrom'][idx], adata.varm['loc'][idx]
            nearby_ = np.where(adata.varm['chrom'] == chromosome)
            idx_ = np.where(np.logical_and(np.asarray(adata.varm['loc'], dtype=int) < int(location + window), np.asarray(adata.varm['loc'], dtype=int) > int(location - window)))
            nearby = np.intersect1d(nearby_, idx_)

            mut_vals = adata.X[:, nearby]
            mut_vals[mut_vals == 3] = 0
            mut_vals[mut_vals > 1] = 1
            sum_muts = np.sum(mut_vals, axis=1)

            if het_idx is not None:
                pure_het_idx = np.intersect1d(np.argwhere(sum_muts == 1), het_idx)
                not_pure_het_idx = np.intersect1d(np.argwhere(sum_muts != 1), het_idx)
                mdata.mod[genotype_key].obs[obs_key][pure_het_idx] = str(variant) + '_het_pure'
                mdata.mod[genotype_key].obs[obs_key][not_pure_het_idx] = str(variant) + '_het_bystander'

                cmn_barcodes, idx_genotype, idx_phenotype = np.intersect1d(mdata.mod[genotype_key].obs_names[not_pure_het_idx],
                                                                          mdata.mod[phenotype_key].obs_names, return_indices=True)
                mdata.mod[phenotype_key].obs[obs_key][idx_phenotype]= str(variant) + '_het_bystander'
                logger.info(f"Annotated {idx_genotype.size} heterozygous mutants with bystanders")


                cmn_barcodes, idx_genotype, idx_phenotype = np.intersect1d(mdata.mod[genotype_key].obs_names[pure_het_idx],
                                                                      mdata.mod[phenotype_key].obs_names, return_indices=True)
                mdata.mod[phenotype_key].obs[obs_key][idx_phenotype] = str(variant) + '_het_pure'
                logger.info(f"Annotated {idx_genotype.size} pure heterozygous mutants")


            if hom_idx is not None:

                pure_hom_idx = np.intersect1d(np.argwhere(sum_muts == 1), hom_idx)
                not_pure_hom_idx = np.intersect1d(np.argwhere(sum_muts != 1), hom_idx)

                mdata.mod[genotype_key].obs[obs_key][pure_hom_idx] = str(variant) + '_hom_pure'
                mdata.mod[genotype_key].obs[obs_key][not_pure_hom_idx] = str(variant) + '_hom_bystander'

                cmn_barcodes, idx_genotype, idx_phenotype = np.intersect1d(mdata.mod[genotype_key].obs_names[not_pure_hom_idx],
                                                                        mdata.mod[phenotype_key].obs_names, return_indices=True)
                mdata.mod[phenotype_key].obs[obs_key][idx_phenotype] = str(variant) + '_hom_bystander'
                logger.info(f"Annotated {idx_genotype.size} homozygous mutants with bystanders")

                cmn_barcodes, idx_genotype, idx_phenotype = np.intersect1d(mdata.mod[genotype_key].obs_names[pure_hom_idx],
                                                                        mdata.mod[phenotype_key].obs_names, return_indices=True)
                mdata.mod[phenotype_key].obs[obs_key][idx_phenotype] = str(variant) + '_hom_pure'
                logger.info(f"Annotated {idx_genotype.size} pure homozygous mutants")


    if ignore_mixed == False:
      # first, lets find the cells that have mixed edits of the intended variants
      indices = []
        for variant in variants:
            try:
                indices.append(adata.var_names.get_loc(variant))
            except KeyError:
                logger.info(str(variant) + ' was not found in the genotype modality anndata')
                continue
        idxs = np.asarray(indices).flatten()

        edit_array = np.zeros(adata.X[:, idxs].shape)
        edit_array[adata.X[:, idxs] == 1] = 1
        edit_array[adata.X[:, idxs] == 2] = 1
        mixed_idxs = np.argwhere(np.sum(edit_array, axis=1) >= 2).ravel()

        mdata.mod[genotype_key].obs[obs_key][mixed_idxs] = 'mixed_mutant'
        cmn_barcodes, idx_genotype, idx_phenotype = np.intersect1d(mdata.mod[genotype_key].obs_names[mixed_idxs],
                                                                      mdata.mod[phenotype_key].obs_names, return_indices=True)
        mdata.mod[phenotype_key].obs[obs_key][idx_phenotype]= 'mixed_mutant'