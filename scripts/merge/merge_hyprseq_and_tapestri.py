import sys
import os
import re
import argparse
import numpy as np
import pandas as pd
import loompy
import logging
import anndata
import seaborn as sns
from mudata import MuData
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")


# Parse the arguments.
def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge HyPRseq and Tapestri output.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-r", "--hyprseq", type=str, required=True, help="HyPRseq output file."
    )  
    parser.add_argument(
        "-d", "--tapestri", type=str, required=True, help="Tapestri output file."
    ) 
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output folder."
    )
    parser.add_argument(
        "-n", "--name", type=str, required=True, help="Variant name table."
    )
    parser.add_argument(
        "-p", "--plot", type=str, default="False", help="Plot the output."
    )
    parser.add_argument(
        "-g", "--germline", type=str, default=None, help="Germline mutations file."
    )
    args = parser.parse_args()
    return args


def read_HyPR_file(HyPR_file):
    """
    Read in the HyPRseq output file. This is a matrix of counts, where every row is a cell barcode and every column is a probe.
    """
    HyPR_array = np.loadtxt(HyPR_file, dtype=str, skiprows=1)
    barcodes = HyPR_array[:, 0]
    with open(HyPR_file, "r") as f:
        first_line = f.readline().strip()
        values = first_line.split("\t")
        var_names = np.array(values, dtype=str)[:]
    X = np.asarray(HyPR_array[:, 1:], dtype=int)
    if np.size(var_names) != X.shape[1]:  # sloppy I shouldn't need this
        var_names = var_names[1:]
    df = pd.DataFrame(X, columns=var_names)
    df.columns = df.columns.str.split("_").str[0]
    collapsed_df = df.groupby(axis=1, level=0).sum()
    collapsed_X = collapsed_df.to_numpy()
    collapsed_var_names = collapsed_df.columns.to_numpy()

    adata_hypr = anndata.AnnData(collapsed_X)
    adata_hypr.obs_names = barcodes

    try:
        adata_hypr.var_names = collapsed_var_names
    except Exception as e:
        print("Error message:")
        print(e)
        adata_hypr.var_names = collapsed_var_names[1:]

    return adata_hypr


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

    variant_names, amplicon_names, chromosome, location, ref, alt = (
        loom_file.ra["id"],
        loom_file.ra["amplicon"],
        loom_file.ra["CHROM"],
        loom_file.ra["POS"],
        loom_file.ra["REF"],
        loom_file.ra["ALT"],
    )

    barcodes = [barcode.split("-")[0] for barcode in loom_file.ca["barcode"]]
    adata = anndata.AnnData(np.transpose(loom_file[:, :]), dtype=loom_file[:, :].dtype)
    adata.layers["e"] = np.transpose(loom_file.layers["AD"][:, :])
    adata.layers["no_e"] = np.transpose(loom_file.layers["RO"][:, :])
    # Also add the GQ and DP layers.
    adata.layers["GQ"] = np.transpose(loom_file.layers["GQ"][:, :])
    adata.layers["DP"] = np.transpose(loom_file.layers["DP"][:, :])

    adata.var_names = variant_names
    adata.obs_names = barcodes
    adata.varm["amplicon"] = amplicon_names
    adata.varm["chrom"] = chromosome
    adata.varm["loc"] = location
    adata.varm["ref"] = ref
    adata.varm["alt"] = alt

    loom_file.close()

    return adata


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
    assert (
        len(mod_names) == 2
    ), "Function not implemented for mod_names with length different from 2"

    obs_1, obs_2 = mdata[mod_names[0]].obs_names, mdata[mod_names[1]].obs_names

    _, idx_1, idx_2 = np.intersect1d(obs_1, obs_2, return_indices=True)

    # Initializing 'intersecting' columns to False
    mdata[mod_names[0]].obs["intersecting"] = False
    mdata[mod_names[1]].obs["intersecting"] = False

    # Updating 'intersecting' columns where True
    mdata[mod_names[0]].obs["intersecting"].iloc[idx_1] = True
    mdata[mod_names[1]].obs["intersecting"].iloc[idx_2] = True
    # Print how many barcodes in each modality to begin with.
    logger.info(f"Found {obs_1.size} barcodes in modality {mod_names[0]}")
    logger.info(f"Found {obs_2.size} barcodes in modality {mod_names[1]}")
    logger.info(f"Found {idx_1.size} intersecting barcodes")


def process_mutant_types(
    mdata,
    variants,
    false_amplicons,
    germline_amplicons,
    rename_rules_df,
    window_size=10,
    plot=False,
    save_path=None):

    def get_new_name(old_name, rename_rules_df=rename_rules_df):
        if old_name in rename_rules_df["old_name"].values:
            return rename_rules_df.loc[rename_rules_df["old_name"] == old_name, "new_name"].iloc[0]


    def get_nearby_variants(variant_names, target_loc, germline_amplicons, window_size):
        """
        Get nearby variants within a window size of the target_loc.
        Args:
            variant_names (list): List of variant names.
            target_loc (str): Target loci (e.g. chr1:123456:A/G).
            germline_amplicons (list): List of germline amplicons.
            window_size (int): Window size.
        Returns:
            list: List of nearby variants.
        """
        chrom_t, loc_t, _ = target_loc.split(":")
        nearby_variants = []
        for variant in variant_names:
            chrom, loc, edit_type = variant.split(":")
            is_valid_edit_type = not edit_type.endswith("*") and re.match(
                r"^[A-Za-z]/[A-Za-z]$", edit_type
            )
            is_not_germline = variant not in germline_amplicons
            if (
                chrom == chrom_t
                and (int(loc_t) - window_size) <= int(loc) <= (int(loc_t) + window_size)
                and is_valid_edit_type
                and is_not_germline
            ):
                nearby_variants.append(variant)

        return nearby_variants

    def process_bystander_cells(target_loc, cells, nearby_variants):
        """
        Identify which cells have bystander editing.
        Args:
            target_loc (str): Target loci (e.g. chr1:123456:A/G).
            cells (list): List of cells (that supposedly have the mutation).
            nearby_variants (list): List of nearby variants.
        Returns:
            list: List of pure cells.
            pd.DataFrame: New matrix.
            pd.DataFrame: Matrix with false amplicons.
        """
        nearby_adata = adata[cells, nearby_variants].copy()
        new_matrix = pd.DataFrame(
            nearby_adata.X, index=nearby_adata.obs_names, columns=nearby_adata.var_names
        )

        # Change all 3 to 0 in the new matrix, but excluding the target loci column.
        for col in new_matrix.columns:
            if col != target_loc:
                new_matrix[col] = new_matrix[col].apply(
                    lambda val: 0 if val == 3 else val
                )

        # Remove the target loci column to generate a bystander position only matrix.
        matrix_without_target_loci = new_matrix.loc[
            :, new_matrix.columns != target_loc
        ].copy()
        bystander_sum_for_every_cell = matrix_without_target_loci.sum(axis=1)
        # Remove cells that have a row sum of 0, or there's no bystander editing.
        matrix_without_target_loci = matrix_without_target_loci.loc[
            bystander_sum_for_every_cell != 0
        ]

        # If there's no bystander editing, return the original cells.
        if len(matrix_without_target_loci) == 0:
            return set(cells), new_matrix

        # Everything that's left is a bystander cell.
        pure_cells = set(cells) - set(matrix_without_target_loci.index)
        return pure_cells, new_matrix

    def plot_mutant_types(target_loc, new_matrix, mt_type, save_path, index=""):
        mutation_counts = new_matrix.apply(np.count_nonzero, axis=0)
        plt.figure(figsize=(10, 6))
        bars = plt.bar(mutation_counts.index, mutation_counts.values)
        for i, idx in enumerate(mutation_counts.index):
            plt.text(
                bars[i].get_x() + bars[i].get_width() / 2,
                bars[i].get_height(),
                str(mutation_counts.values[i]),
                ha="center",
                va="bottom",
            )
        if target_loc in mutation_counts.index:
            bars[mutation_counts.index.tolist().index(target_loc)].set_color("r")
        red_patch = mpatches.Patch(color="red", label="Target Loci")
        plt.legend(handles=[red_patch])

        plot_title = (
            f"Mutation Counts: {mt_type} - {get_new_name(target_loc)} ({index}) with total cells {len(new_matrix)}"
        )
        plt.title(plot_title)
        plt.xticks(rotation=90)
        plt.xlabel("")
        plt.ylabel("Number of Mutated Cells")

        save_index = f"{get_new_name(target_loc)}_{mt_type}_mutant_counts_{index}.pdf"
        plt.savefig(os.path.join(save_path, save_index), bbox_inches="tight")
        plt.close()

        # Plot a cluster map of matrix but cluster rows only.
        sns.clustermap(
            new_matrix,
            cmap="YlGnBu",
            row_cluster=True,
            col_cluster=False,
            figsize=(10, 10),
        )
        plt.title(
            f"Mutation Matrix: {mt_type} - {get_new_name(target_loc)} aka {target_loc} with shape {new_matrix.shape}"
        )
        plt.xlabel("Variant")
        plt.ylabel("Cell Barcode")
        save_index = f"{get_new_name(target_loc)}_{mt_type}_mutant_matrix_{index}.pdf"
        plt.savefig(os.path.join(save_path, save_index), bbox_inches="tight")
        plt.close()

    def plot_output_df(df, save_path):
        """
        Plot the df as a heatmap.
        """
        df_to_plot = df.copy()
        # Filter for rows where the sum is not 0.
        df_to_plot = df_to_plot.loc[df_to_plot.sum(axis=1) != 0]

        n_colors = 5
        cmap = plt.get_cmap("Paired", n_colors)  
        sns.clustermap(df_to_plot, cmap=cmap, row_cluster=True, col_cluster=False, figsize=(10, 10))
        plt.title("Mutation Matrix")
        plt.savefig(os.path.join(save_path, "cell_barcode_all_mutant_calls.pdf"))
        plt.close()

    
    def process_single(
        target_loc,
        adata,
        variant_names,
        window_size,
        plot,
        save_path,
        false_amplicons,
        germline_amplicons):

        """
        Process a single target loci. Output is a df with cell barcodes as index and target loci as column.
        The values represent the mutant type:
            0: Background
            1: Pure heterozygous
            2: Heterozygous with bystander editing
            3: Pure homozygous
            4: Homozygous with bystander editing
        """
        idx = adata.var_names.get_loc(target_loc)
        bkg_cells, het_cells, hom_cells = [
            adata.obs_names[np.flatnonzero(adata.X[:, idx] == i)] for i in (0, 1, 2)
        ]
        nearby_variants = get_nearby_variants(
            variant_names, target_loc, germline_amplicons, window_size
        )

        # print(f"Processing: {target_loc} heterozygous")
        het_pure, het_matrix = process_bystander_cells(
            target_loc, het_cells, nearby_variants
        )
        # print(f"Processing: {target_loc} homozygous")
        hom_pure, hom_matrix = process_bystander_cells(
            target_loc, hom_cells, nearby_variants
        )

        if plot is True:
            plot_mutant_types(target_loc, het_matrix, "Heterozygous", save_path, "1-0")
            plot_mutant_types(target_loc, hom_matrix, "Homozygous", save_path, "1-0")

        df = pd.DataFrame(0, index=adata.obs_names, columns=[target_loc])
        df.loc[list(het_pure), target_loc] = 1
        df.loc[list(set(het_cells) - het_pure), target_loc] = 2
        df.loc[list(hom_pure), target_loc] = 3
        df.loc[list(set(hom_cells) - hom_pure), target_loc] = 4
        df.loc[list(bkg_cells), target_loc] = 0
        return df

    def call_AAVS_cells(variant_names):
        aavs4 = []
        aavs1 = []
        for variant in variant_names:
            chrom, loc, edit_type = variant.split(":")
            if (
                chrom == "chr19"
                and 55115745 <= int(loc) <= 55115764
                and (edit_type == "A/G" or edit_type == "T/C")
            ):
                aavs4.append(variant)
            if (
                chrom == "chr19"
                and 55115752 <= int(loc) <= 55115771
                and (edit_type == "C/T" or edit_type == "G/A")
            ):
                aavs1.append(variant)
        return set(aavs1), set(aavs4)

    def process_aavs(
        aavs,
        adata,
        variant_names,
        window_size,
        col_index,
        plot,
        save_path,
        false_amplicons,
        germline_amplicons=germline_amplicons):

        variant_dfs = []
        for target_loc in aavs:
            df_aavs = process_single(
                target_loc,
                adata,
                variant_names,
                window_size,
                plot=plot,
                save_path=save_path,
                false_amplicons=false_amplicons,
                germline_amplicons=germline_amplicons,
            )
            variant_dfs.append(df_aavs)
        if len(variant_dfs) == 0:
            return pd.DataFrame(0, index=adata.obs_names, columns=[col_index])
        else:
            combined_df_aavs = pd.concat(variant_dfs, axis=1)
            df_aavs = combined_df_aavs.sum(axis=1)
            df_aavs_final = pd.DataFrame(
                df_aavs.apply(lambda x: 1 if x != 0 else 0), columns=[col_index]
            )
            return df_aavs_final

    def annotate_mutant_types(row):
        non_zero_count = (row != 0).sum()
        total_sum = row.sum()

        if total_sum == 0:
            return "Unedited", "Unedited"
        elif row["AAVS"] == 1 and total_sum == 1:
            return "AAVS", "AAVS"

        long_format = []
        short_format = []
        for col, val in row.items():
            if col != "AAVS" and val != 0:
                long_format.append(f"{col}-{val}")
                mutant_type = "Het" if val in [1, 2] else "Hom"
                short_format.append(f"{col}-{mutant_type}")

                if non_zero_count == 1:
                    return f"{col}-{val}", f"{col}-{mutant_type}"
        
        long_format_string = ", ".join(long_format)
        short_format_string = ", ".join(short_format)

        return "Mixed;" + long_format_string, "Mixed;" + short_format_string

    def replace_name(row, rename_col, rename_rules_df):
        parts = row[rename_col].split("-")
        old_name = parts[0]
        suffix = "-".join(parts[1:]) if len(parts) > 1 else ""

        new_name_row = rename_rules_df.loc[rename_rules_df["old_name"] == old_name]
        if not new_name_row.empty:
            new_name = new_name_row["new_name"].values[0]
            return f"{new_name}-{suffix}" if suffix else new_name
        else:
            return row[rename_col]

    adata = mdata.mod["tapestri"].copy()
    variant_names = np.asarray(mdata.mod["tapestri"].var_names.values)

    # AAVS
    aavs1, aavs4 = call_AAVS_cells(variant_names)
    df_aavs_final_1 = process_aavs(
        aavs1,
        adata,
        variant_names,
        window_size,
        "AAVS1",
        plot=False,
        save_path=None,
        false_amplicons=false_amplicons,
    )
    df_aavs_final_4 = process_aavs(
        aavs4,
        adata,
        variant_names,
        window_size,
        "AAVS4",
        plot=False,
        save_path=None,
        false_amplicons=false_amplicons,
    )
    df_aavs_final = process_aavs(
        aavs4,
        adata,
        variant_names,
        window_size,
        "AAVS",
        plot=False,
        save_path=None,
        false_amplicons=false_amplicons,
    )
    df_aavs_sub = pd.concat([df_aavs_final_1, df_aavs_final_4], axis=1)

    # Variant
    variant_dfs = []
    for target_loc in variants:
        try:
            # print(adata.var_names)
            # Save the var names to a txt file. 
            with open(os.path.join(save_path, "all_variant_names.txt"), "w") as f:
                for item in adata.var_names:
                    f.write("%s\n" % item)
            df = process_single(
                target_loc,
                adata,
                variant_names,
                window_size,
                plot=plot,
                save_path=save_path,
                false_amplicons=false_amplicons,
                germline_amplicons=germline_amplicons,
            )
            variant_dfs.append(df)
        except KeyError:
            logger.info(
                str(target_loc) + " was not found in the genotype modality anndata"
            )
        continue

    if len(variant_dfs) == 0:
        return pd.DataFrame(0, index=adata.obs_names, columns=variant_names), df_aavs_sub
    combined_df = pd.concat(variant_dfs, axis=1)
    df = pd.concat([df_aavs_final, combined_df], axis=1)

    # Write this df to output.
    df_output_filename = os.path.join(save_path, "cell_barcode_raw_mutant_calls.csv")
    df.to_csv(df_output_filename, sep=",", index=True)
    # Plot this as heatmap.
    plot_output_df(df, save_path)

    # Annotation
    df[["Mutant_type", "Mutant_type_short"]] = df.apply(annotate_mutant_types, axis=1, result_type="expand")
    df = df[["Mutant_type", "Mutant_type_short"]]
    df["Name_Mutant_type"] = df.apply(replace_name, rename_col="Mutant_type", rename_rules_df=rename_rules_df, axis=1)
    df["Name_Mutant_type_short"] = df.apply(replace_name,rename_col="Mutant_type_short",rename_rules_df=rename_rules_df,axis=1,)

    return df, df_aavs_sub


def read_germline_file(germline_file):
    germline_set = set()
    with open(germline_file, 'r') as file:
        for line in file:
            item = line.strip()
            germline_set.add(item)
    return germline_set


def plot_non_intersecting_hyprseq_barcodes(intersecting_barcodes_df, non_intersecting_barcodes_df, save_path):
    # Plot histogram of counts per barcode for intersecting and non-intersecting barcodes.
    plt.figure(figsize=(10, 6))
    plt.hist(intersecting_barcodes_df.sum(axis=1), bins=50, alpha=0.5, label="Intersecting")
    plt.hist(non_intersecting_barcodes_df.sum(axis=1), bins=50, alpha=0.5, label="Non-intersecting")
    plt.title("Histogram of UMI per cell barcode for HyPR-seq")
    plt.xlabel("Counts")
    plt.ylabel("Number of barcodes")
    plt.legend()
    plt.savefig(os.path.join(save_path, "counts_per_barcode.pdf"), bbox_inches="tight")
    return


if __name__ == "__main__":
    args = parse_args()
    hyprseq_path = args.hyprseq
    tapestri_path = args.tapestri
    germline_path = args.germline
    output_path = args.output
    name_table = args.name
    plot = args.plot

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("####################################################################################################")
    logger.info("Running script with the following arguments:")
    logger.info(f"HyPRseq file: {hyprseq_path}")
    logger.info(f"Tapestri file: {tapestri_path}")
    logger.info(f"Output path: {output_path}")
    logger.info(f"Variant name table: {name_table}")
    if germline_path is None:
        logger.info("Germline mutations file: Default (germline_mutations.txt)")
    else:
        logger.info(f"Germline mutations file: {germline_path}")
    logger.info("####################################################################################################")


    if germline_path is not None:
        germline_amplicons_list = read_germline_file(germline_path)
    else:
        # This is the list that is filtered from the P500/P1000/P2000 dataset.
        germline_amplicons_list = read_germline_file('./germline_mutations.txt')

    # Read in the variant name table.
    rename_rules_df = pd.read_csv(name_table, names=['old_name', 'new_name'])
    variant_list = rename_rules_df["old_name"].tolist()

    # Create output directory if it doesn't exist.
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    fig_path = os.path.join(output_path, "figures")
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    intermediate_fig_path = os.path.join(output_path, "intermediate_figures")
    if not os.path.exists(intermediate_fig_path):
        os.makedirs(intermediate_fig_path)

    # Read in the HyPRseq output file.
    adata_hypr = read_HyPR_file(hyprseq_path)
    # Read in the Tapestri output file.
    adata_tapestri = read_tapestri_loom(tapestri_path)

    mdata = MuData({"HyPR": adata_hypr, "tapestri": adata_tapestri})
    # Find intersecting barcodes and add 'intersecting' to obs with boolean values.
    find_intersecting(mdata)

    # Find the non-intersecting barcodes in the HyPRseq modality.
    non_intersecting_barcodes = mdata.mod["HyPR"][mdata.mod["HyPR"].obs["intersecting"] == False].copy()
    
    # Get the raw counts for the non-intersecting barcodes and convert to pandas dataframe.
    non_intersecting_barcodes_df = pd.DataFrame(
        non_intersecting_barcodes.X,
        index=non_intersecting_barcodes.obs_names,
        columns=non_intersecting_barcodes.var_names,
    )
    
    # Save the non-intersecting barcodes to output.
    non_intersecting_barcodes_filename = os.path.join(output_path, "HyPR_cells_non_intersect_with_tapestri.tsv")
    logger.info(f"Writing non-intersecting barcodes to {non_intersecting_barcodes_filename}")
    non_intersecting_barcodes_df.to_csv(
        non_intersecting_barcodes_filename, sep="\t", index=True
    )

    # Get the raw counts for the intersecting_barcodes and convert to pandas dataframe.
    intersecting_barcodes = mdata.mod["HyPR"][mdata.mod["HyPR"].obs["intersecting"] == True].copy()
    intersecting_barcodes_df = pd.DataFrame(
        intersecting_barcodes.X,
        index=intersecting_barcodes.obs_names,
        columns=intersecting_barcodes.var_names,
    )

    # Plot some stats about the non-intersecting barcodes.
    if plot is True:
        plot_non_intersecting_hyprseq_barcodes(intersecting_barcodes_df, non_intersecting_barcodes_df, intermediate_fig_path)

    # Filter to data that are intersecting.
    mdata.mod["tapestri"] = mdata.mod["tapestri"][mdata.mod["tapestri"].obs["intersecting"] == True].copy()
    mdata.mod["HyPR"] = mdata.mod["HyPR"][mdata.mod["HyPR"].obs["intersecting"] == True].copy()

    # Sort the mdata.mod['HyPR'] alphabetically by cell barcodes.
    mdata.mod["HyPR"] = mdata.mod["HyPR"][mdata.mod["HyPR"].obs_names.sort_values()].copy()
    # Sort the mdata.mod['tapestri'] by the same order as the mdata.mod['HyPR'].
    mdata.mod["tapestri"] = mdata.mod["tapestri"][mdata.mod["HyPR"].obs_names].copy()

    # Process the mutant types.
    # Removed the false amplicons list for now.
    false_amplicons = []
    df, _ = process_mutant_types(
        mdata,
        variant_list,
        false_amplicons,
        germline_amplicons_list,
        rename_rules_df,
        plot=plot,
        save_path=fig_path
    )

    # Write df to output.
    logger.info(f"Writing df to {output_path}")
    df_output_filename = os.path.join(output_path, "cell_barcode_DNA_mutations_calls.tsv")
    df.to_csv(df_output_filename, sep="\t", index=True)

    # Add the df to the mdata.
    df = df.reindex(mdata.mod["HyPR"].obs_names)
    mdata.mod["HyPR"].obs = df.copy()
    mdata.mod["tapestri"].obs = df.copy()

    # Write mdata to output as a h5mu file.
    mdata_filename = os.path.join(output_path, "merged_annotation.h5mu")
    logger.info(f"Writing mdata to {mdata_filename}")
    mdata.write(mdata_filename)
