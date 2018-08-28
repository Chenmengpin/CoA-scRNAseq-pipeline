if __name__ == '__main__':  # required to run outside jupyter notebook

    import os  # used to interface with the operating system on a basic level
    import glob  # finds pathnames for UNIX
    import pickle # allows serialization of objects
    import pandas as pd  # required for data array manipulation

    from dask.diagnostics import ProgressBar  # creates progress bar for cisTarget run

    from arboreto.utils import load_tf_names   # needed to import TFs in
    from arboreto.algo import grnboost2  # needed to run GENIE3

    from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase # imports cisTarget database metadata
    from pyscenic.utils import modules_from_adjacencies, load_motifs # creates modules from GENIE3 adjacencies
    from pyscenic.prune import prune2df, df2regulons
    from pyscenic.aucell import aucell
    
    # load paths for repeatedly invoked files
    motifs_filename = os.path.join("data/motifs.csv")
    regulons_filename = os.path.join("data/regulons.p")
    TF_list_filename = os.path.join("resources/mm_tfs.txt")
    
    # this cell creates the TF list if not made yet
    tf_raw = pd.read_csv("resources/TF_import.txt", delimiter = "\t",     # imports in the raw gene annotation information from TFCat 
                         error_bad_lines = False, encoding = "ISO-8859-1")  
    tfs = tf_raw[["Gene ID", "Evidence Strength"]].drop_duplicates().dropna()   # cleans raw TF annotations
    tfs["ID"] = list(map(int, tfs["Gene ID"]))      # lists the genes by ID
    conv_tfs = pd.read_csv("resources/TF_conversion.txt", delimiter = "\t")   # imports in MGI gene names from DAVID
    def extract_symbol(name):       # extracts the abbreviated ID from full name
        s_idx = name.rfind('(')
        e_idx = name.rfind(')')
        return name[s_idx+1:e_idx]
    conv_tfs["Gene Name"].apply(extract_symbol).to_csv(TF_list_filename, index = False)     # turns list into CSV for future import
    
    # this cell loads TF list, expression matrix, and databases for downstream SCENIC analysis
    tf_names = load_tf_names(TF_list_filename)      # if TF list has been made, this imports it
    ex_matrix = pd.read_csv("data/GENIE3_import.csv", header = 0, index_col = 0).T    # loads expression matrix, make sure you transpose back
    databases_glob = os.path.join("databases/mm10__*.feather") # loads cisTarget databases into memory
    db_fnames = glob.glob(databases_glob)
    def name(fname):
        return os.path.basename(fname).split(".")[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    
    # GENIE3 process: returns co-expression modules
    adjacencies = grnboost2(ex_matrix, tf_names = tf_names, verbose = True)     # runs improved GRNBoost instance of GENIE3
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))    # identifies modules from GENIE3
    
    # cisTarget process: IDs cis-regulatory footprints from motifs around the TSS
    with ProgressBar(): # calculate a list of enriched motifs and the corresponding target genes for all modules
        df = prune2df(dbs, modules, "resources/motifs-v9-nr-mgi.txt")
    regulons = df2regulons(df)  # create regulons from this table of enriched motifs
    
    # save the discovered motifs and regulons
    df.to_csv(motifs_filename)  
    with open(regulons_filename, "wb") as f:
        pickle.dump(regulons, f)
    
    # load the discovered motifs and regulons if saved previously
    df = load_motifs(motifs_filename)
    with open(regulons_filename, "rb") as f:
        regulons = pickle.load(f)
        
    # AUCell process: finds enrichment of each discovered regulon
    auc_matrix = aucell(ex_matrix, regulons, num_workers = 4)

    
    
    
    
