import scanpy as sc
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import ktplotspy as kpy
import pandas as pd
from IPython.display import display, HTML
from pathlib import Path

cluster_key = 'seurat_clusters'
celltype_key = 'celltype_orig'
cellphonedb_celltype_key = 'celltype_cellphonedb'
samples = ['UMM059', 'UMM063', 'UMM064', 'UMM065', 'UMM066']
hybrid_clusters_list = [[8], [3, 12], [11], [11], [6]]

# Shared code for all samples
def create_celltype_for_cellphonedb(adata: sc.AnnData, 
                                    hybrid_clusters: list,
                                    celltype_key: str = celltype_key,
                                    cluster_key: str = cluster_key,
                                    cellphonedb_celltype_key: str = cellphonedb_celltype_key,
                                    save_file: str = None):
    def map_func(row):
        if (row[cluster_key] in hybrid_clusters):
            return 'Hybrid'
        if 'tumor' in row[celltype_key].lower():
            return 'Tumor'
        return row[celltype_key]
    adata.obs[cellphonedb_celltype_key] = adata.obs.apply(map_func, axis = 1)
    if save_file is not None:
        adata.obs.to_csv(save_file, columns=[cellphonedb_celltype_key], index_label='cell')
    return adata


def filter_cell_types(adata: sc.AnnData,
                      cellphonedb_celltype_key: str):
    # To control the display, list cell types have at least 20 cells
    celltype_counts = adata.obs.groupby(cellphonedb_celltype_key)['orig.ident'].count()
    print(celltype_counts)
    selected_celltypes = celltype_counts > 20 # Filter out fFollicular B-Cells too
    selected_celltypes = celltype_counts.index[selected_celltypes].to_list()
    # Filter the adata so that we care about what we need
    adata = adata[adata.obs[cellphonedb_celltype_key].isin(selected_celltypes)]
    # This format is required by the following function call to plot_cpdb 
    selected_celltypes = '|'.join(selected_celltypes)
    return selected_celltypes


def plot_cpdb(adata: sc.AnnData,
              cell_type1: str,
              selected_celltypes: str,
               means: pd.DataFrame,
               pvalues: pd.DataFrame,
               cellphonedb_celltype_key: str,
               title: str,
               figsize = (10, 5)):
    ggplot = kpy.plot_cpdb(
        adata = adata,
        cell_type1 = cell_type1,
        cell_type2 = selected_celltypes, 
        means = means,
        pvals = pvalues,
        celltype_key = cellphonedb_celltype_key,
        # genes = ["TGFB2", "CSF1R"],
        figsize = figsize,
        title = title,
        max_size = 6,
        highlight_size = 0.75,
        standard_scale = True,
        special_character_regex_pattern='/',
        # gene_family="costimulatory"
    )
    print(ggplot)
    # Print out the table
    cpdb_table = kpy.plot_cpdb(
        adata = adata,
        cell_type1 = cell_type1,
        cell_type2 = selected_celltypes, 
        means = means,
        pvals = pvalues,
        celltype_key = cellphonedb_celltype_key,
        special_character_regex_pattern='/',
        standard_scale = True,
        return_table=True
    )
    # print_cpdb_table(cpdb_table)
    return cpdb_table


def print_cpdb_table(table: any):
    if isinstance(table, dict):
        for key, value in table.items():
            print(key)
            display(HTML(value.to_html(max_rows=25)))
    elif isinstance(table, pd.DataFrame):
        display(HTML(table.to_html(max_rows=25)))


def run_cellphonedb(sample: str, 
                    hybrid_clusters: list,
                    cellphonedb_dir: str,
                    cpdb_file_path: str,
                    base_dir: str,
                    selected_celltypes: list = None,
                    iterations: int = 1000):
    """Here the default 1,000 samplings are used. This is different from Chris' analysis. He used a much
    larger sampling number, 1,000,000. Not sure if it is really necessary for this type of analysis.

    Args:
        sample (str): _description_
        hybrid_clusters (list): _description_
        cellphonedb_dir (str): _description_
        cpdb_file_path (str): _description_
        base_dir (str): _description_

    Returns:
        _type_: _description_
    """
    cellphonedb_meta_file = cellphonedb_dir + 'meta_{}.csv'.format(sample)
    counts_h5ad_file = base_dir + 'rna_{}.h5ad'.format(sample)
    cellpohone_out_dir = cellphonedb_dir + sample
    if not Path(cellpohone_out_dir).exists:
        Path(cellpohone_out_dir).mkdir()
    adata = sc.read_h5ad(counts_h5ad_file)
    # Need to create cell types for cellphone db. Here we merge all tumor cells together as a single tumor cell type
    # Also based on clusters information to create hybrid clusters
    create_celltype_for_cellphonedb(adata, 
                                    hybrid_clusters, 
                                    save_file=cellphonedb_meta_file,
                                    celltype_key=celltype_key,
                                    cluster_key=cluster_key,
                                    cellphonedb_celltype_key=cellphonedb_celltype_key)
    # Run the random sampling based method. For the time being, using 1,000 permuations
    deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cpdb_file_path,                 # mandatory: CellPhoneDB database zip file.
        meta_file_path = cellphonedb_meta_file,                 # mandatory: tsv file defining barcodes to cell label.
        counts_file_path = counts_h5ad_file,             # mandatory: normalized count matrix.
        counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
        microenvs_file_path = None,       # optional (default: None): defines cells per microenvironment.
        iterations = iterations,                               # denotes the number of shufflings performed in the analysis.
        threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
        threads = 4,                                     # number of threads to use in the analysis.
        debug_seed = 42,                                 # debug randome seed. To disable >=0.
        result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
        pvalue = 0.05,                                   # P-value threshold to employ for significance.
        subsampling = False,                             # To enable subsampling the data (geometri sketching).
        subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
        subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
        subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
        separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
        debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
        output_path = cellpohone_out_dir,                          # Path to save results.
        output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
        )
    # Output some results
    print("Heatmap without filtering out any cell types")
    kpy.plot_cpdb_heatmap(
        adata = adata,
        pvals = pvalues,
        celltype_key = cellphonedb_celltype_key,
        figsize = (5,5),
        title = "{}: Number of significant interactions".format(sample),
        symmetrical = False,
        return_tables=False,
    )
    print("Heatmap after filtering cell types having less than 20 cells")
    # First filter
    if selected_celltypes is None: # Do filtering based on cell types. This may give us different cell types across different samples
        selected_celltypes = filter_cell_types(adata=adata, 
                                               cellphonedb_celltype_key=cellphonedb_celltype_key)
    print('Selected cell types: {}'.format(selected_celltypes))
    # Filter the pvalues for the plot
    print("Filtering pvalues based on selected cell types...")
    print("Before filtering: {}".format(pvalues.shape))
    def filter_columns_based_on_cell_types(column: str):
        if '|' in column:
            # Make sure both of the cell types are selected
            cell_type_names = column.split('|')
            return (cell_type_names[0] in selected_celltypes) and (cell_type_names[1] in selected_celltypes)
        else:
            return True
    filter_pvalues = pvalues.filter(items=filter(filter_columns_based_on_cell_types, pvalues.columns), axis=1)
    print('After filtering: {}'.format(filter_pvalues.shape))

    kpy.plot_cpdb_heatmap(
        adata = adata,
        pvals = filter_pvalues,
        celltype_key = cellphonedb_celltype_key,
        figsize = (5,5),
        title = "{}: Number of significant interactions".format(sample),
        symmetrical = False,
        return_tables=False,
    )

    print("Dotplot for hybrid")
    hybrid_cpdb_table = plot_cpdb(
               adata,
               'Hybrid',
               selected_celltypes,
               means,
               pvalues,
               cellphonedb_celltype_key,
               title='{}: Interactions for Hybrid Cells'.format(sample))
    print("Dot plot for tumor")
    tumor_cpdb_table = plot_cpdb(
               adata,
               'Tumor',
               selected_celltypes,
               means,
               pvalues,
               cellphonedb_celltype_key,
               title='{}: Interactions for Tumor Cells'.format(sample),
               figsize=(10, 17))
    print("Dot plot for macrophage")    
    macrophage_cpdb_table = plot_cpdb(
               adata,
               'Macrophages/Monocytes',
               selected_celltypes,
               means,
               pvalues,
               cellphonedb_celltype_key,
               title='{}: Interactions for Macrophage Cells'.format(sample),
               figsize=(10, 27))
    combined_cpdb_table = plot_cpdb(
               adata,
               'Hybrid|Tumor|Macrophages/Monocytes',
               selected_celltypes,
               means,
               pvalues,
               cellphonedb_celltype_key,
               title='{}: Interactions for Hybrid, Tumor and Macrophage Cells'.format(sample),
               figsize=(10, 30))
    # Save these tables
    filenames = ['hybrid_cpdb_table_{}.csv'.format(sample), 
                 'tumor_cpdb_table_{}.csv'.format(sample),
                 'macrophage_cpdb_table_{}.csv'.format(sample),
                 'combined_cpdb_table_{}.csv'.format(sample)]
    tables = [hybrid_cpdb_table, tumor_cpdb_table, macrophage_cpdb_table, combined_cpdb_table]
    for i in range(len(filenames)):
        filename = cellphonedb_dir + filenames[i]
        print('Saving ' + filename)
        # signficant interactions only
        table_sig = tables[i]
        table_sig = table_sig[table_sig['significant'] == 'yes']
        table_sig.to_csv(filename)
    return hybrid_cpdb_table, tumor_cpdb_table, macrophage_cpdb_table, combined_cpdb_table


def plot_overlapped_interactions(selected_celltypes: str,
                                 dir: str,
                                 base_dir: str='/Volumes/ssd/results/missy_sc_rna/uvm/cellphonedb/',
                                 overlap_interaction_file: str = 'combined_cpdb_table_overlap.csv'):
    # Try to dotplot for shared interactions only
    shared_for_all_df = pd.read_csv(base_dir + overlap_interaction_file, index_col=False)
    shared_for_all_df = shared_for_all_df[shared_for_all_df['overlap'] == 5]
    selected_interactions = shared_for_all_df['interaction_group'].to_list()
    selected_interactions = [i.replace('-', '_') for i in selected_interactions]
    file_patterns = [
        'statistical_analysis_{}_08_03_2023_11:49:50.txt',
        'statistical_analysis_{}_08_03_2023_11:51:18.txt',
        'statistical_analysis_{}_08_03_2023_11:53:38.txt',
        'statistical_analysis_{}_08_03_2023_11:56:16.txt',
        'statistical_analysis_{}_08_03_2023_11:58:22.txt'
    ]
    for i in range(len(samples)):
        sample = samples[i]
        hybrid_cluster = hybrid_clusters_list[i]
        file_pattern = file_patterns[i]
        means = pd.read_csv('{}{}'.format(base_dir, sample) + '/' + file_pattern.format('means'), sep='\t', index_col=False)
        means = means[means['interacting_pair'].isin(selected_interactions)]
        pvalues = pd.read_csv('{}{}'.format(base_dir, sample) + '/' + file_pattern.format('pvalues'), sep='\t', index_col=False)
        pvalues = pvalues[pvalues['interacting_pair'].isin(selected_interactions)]
        counts_h5ad_file = dir + 'rna_{}.h5ad'.format(sample)
        adata = sc.read_h5ad(counts_h5ad_file)
        adata = create_celltype_for_cellphonedb(adata, hybrid_cluster)

        ggplot = kpy.plot_cpdb(
            adata = adata,
            cell_type1 = 'Hybrid|Macrophages/Monocytes|Tumor',
            cell_type2 = selected_celltypes, 
            means = means,
            pvals = pvalues,
            celltype_key = cellphonedb_celltype_key,
            figsize = (10, 12),
            title = 'Selected Interactions for {}'.format(sample),
            max_size = 6,
            highlight_size = 0.75,
            standard_scale = True,
            special_character_regex_pattern='/',
        )
        print(ggplot)
        filename = '{}{}_Selected_Interactions.pdf'.format(base_dir, sample)
        # Save the file
        ggplot.save(filename)