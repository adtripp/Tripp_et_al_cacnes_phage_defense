#!/usr/bin/python

## Created Sept 2023 by ADT

#Thi is a script to find truncations in genome assemblies for genes of interest

#Import packages
import os
import sys
import pandas as pd
import numpy as np
import glob
from Bio import SeqIO
import io
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Conda env requirements
# conda install -c bioconda bioawk
# conda install -c conda-forge biopython
#

# # # FUNCTIONS # # #

def sample_file_guide(clean_samples_list, path_output_dir, prokka_output_dir):
    out_file = os.path.join(path_output_dir, "sample_files_status.csv")
    if not os.path.exists(out_file):
        print("Creating sample file guide.")
        col_names = ["sample", "path", "fna_file", "fsa_file", "gff_file", "tbl_file", "contig_base"]
        status_df = pd.DataFrame()
        contig_id = np.nan
        for idx, sample in enumerate(clean_samples_list):
            if idx % 100 == 0:
                print("On idx {} of {}".format(idx, len(clean_samples_list)))
            sample_path = os.path.join(prokka_output_dir, sample)
            prokka_fna_path = os.path.join(sample_path, "prokka_out.fna")
            prokka_fsa_path = os.path.join(sample_path, "prokka_out.fsa")
            prokka_gff_path = os.path.join(sample_path, "prokka_out.gff")
            prokka_tbl_path = os.path.join(sample_path, "prokka_out.tbl")

            contig_id = np.nan
            if os.path.exists(prokka_gff_path):
                with open (prokka_gff_path, "r") as file:
                    line = file.readline()
                    line = file.readline()
                    contig_id = line.split(" ")[1].split("|")[-1].split("_")[0]

            df_row = [sample, sample_path]
            for path in [prokka_fna_path, prokka_fsa_path, prokka_gff_path, prokka_tbl_path]:
                if os.path.exists(path):
                    df_row.append(True)
                else:
                    df_row.append(False)
            df_row.append(contig_id)
            row = pd.DataFrame(df_row).T
            row.columns = col_names
            status_df = pd.concat([status_df, row])
        status_df.to_csv(out_file)


def makeblastdb(path_output_files, path_clean_samples, input_prokka_files):
    os.system("mkdir -p {}".format(path_output_files))

    # # # Concatenate contigs for database
    blastdb_input = os.path.join(path_output_files, "blast_db_input.fna")
    if os.path.exists(blastdb_input):
        print("The file {} already exists, continuing on ... ".format(blastdb_input))
    if not os.path.exists(blastdb_input):
        print("Need to make blast database... Concatenating prokka nucleotide annotations and storing to {}".format(blastdb_input))
        with open(blastdb_input, 'w+'):
            count = 0
            for idx, path in enumerate(path_clean_samples):
                path_full = os.path.join(path, input_prokka_files)
                if idx % 100 == 0:
                    print("{} / {} done".format(idx, len(path_clean_samples)))
                if os.path.exists(path_full):
                    cat_command = "cat {} >> {}".format(path_full, blastdb_input)
                    os.system(cat_command)
                else:
                    count += 1
            print("{} files did not exist.".format(count))

    # # # Create blast database
    blastdb_output = os.path.join(path_output_files, "blast_db_output/")
    if not os.path.exists(blastdb_output):
        print("Creating BLAST directory...")
        os.system("mkdir -p {}".format(blastdb_output))
        os.system("makeblastdb -in {} -out {} -dbtype nucl -title concat_db_prokka_contigs -parse_seqids".format(blastdb_input, blastdb_output))
    else:
        print("Blast database already exists...")

def run_blastn(input_blastdb, input_query, path_output_files):
    os.system("mkdir -p {}".format(path_output_files))
    os.system("module load engaging/ncbi-blast/2.6.0+")
    # # BLAST query sequences against contig database
    output_blast_hits = os.path.join(path_output_files, "output_blast_hits.csv")
    blast_options = "10 qseqid sseqid length pident mismatch qcovs qcovhsp qstart qend qlen sstart send slen"
    if os.path.exists(output_blast_hits):
        print("Output {} already exists. Not running BLAST...".format(output_blast_hits))
    if not os.path.exists(output_blast_hits):
        print("Running BLAST...")
        os.system("blastn -query {} -db {}/contig_db -outfmt '{}' -out {}".format(input_query, input_blastdb, blast_options, output_blast_hits))
    #

# # # IMPORT DATA # # #

#Import input data
# Path to fasta file of nucleotide sequences for gene of interest
# path_goi = "/home/atripp/mit_lieberman/projects/adt_general/ACERA_PROJECT/DATA/ASSEMBLIES/ISOLATE/6-orthologinfo/3-PhageDefense/pde-hit_ids.fna"
# path_goi = "/home/atripp/mit_lieberman/projects/adt_general/ACERA_PROJECT/DATA/ASSEMBLIES/ISOLATE/6-orthologinfo/3-PhageDefense/defense_finder_hmmer.fna"
path_goi = "/home/atripp/mit_lieberman/projects/rotations/2023_08_rashi/data/cacnes_clade_pde_systems.fna"

#PAth to directory where prokka annotations are stored
# path_prokka = "/home/atripp/mit_lieberman/projects/adt_general/ACERA_PROJECT/DATA/ASSEMBLIES/ISOLATE/4-annotation/Isolate_annotation_filtered/"
path_prokka = "/home/atripp/mit_lieberman/projects/adt_general/ACERA_PROJECT/DATA/ASSEMBLIES/ISOLATE/4-Pure_Isolate_Annotations/C_acnes-Baker2023/"

# Path to file with SampleNames of clean assemblies
path_clean_samples = "/home/atripp/mit_lieberman/projects/adt_general/ACERA_PROJECT/DATA/ASSEMBLIES/ISOLATE/clean_assemblies_stats_meta-c_acnes.txt" #These assemblies were filtered based on bracken, assembly length, number of total contigs

# PAth contig database
path_blast_db = "/home/atripp/mit_lieberman/projects/adt_general/ACERA_PROJECT/DATA/ASSEMBLIES/ISOLATE/9-Truncation_Analysis/C_acnes/blast_db/"

# Output directory
path_output = sys.argv[1]
os.system("mkdir -p {}".format(path_output))

#Read in clean sample names
clean_samples = np.loadtxt(path_clean_samples, dtype = 'str')
# Store path to prokka annotations for clean samples
path_prokka_clean_samples = [os.path.join(path_prokka, sample) for sample in clean_samples]


# # # ANALYSIS # # #

# Check for presence of required prokka annotation files for each samples
# sample_file_guide(clean_samples, path_output, path_prokka)
# makeblastdb(path_output, path_clean_samples, "prokka_out.fna")
# run_blastn(path_blast_db, path_goi, path_output)

#
# ## Dataframe Parsing # # #
# # # Turn BLAST results into df
blast_options = "10 qseqid sseqid length pident mismatch qcovs qcovhsp qstart qend qlen sstart send slen"
blast_hits = pd.read_csv(os.path.join(path_output, "output_blast_hits.csv"), names = blast_options.split(" ")[1:], header=None)

# Create and populate wrangled dataframe
blast_df = blast_hits[['qseqid','sseqid','length', 'qcovs','qcovhsp','pident', 'mismatch','qlen', 'slen']]

start_condition = blast_hits[['sstart','send']].min(axis=1) < 5
end_condition = blast_hits[['sstart','send']].max(axis=1) > blast_hits["slen"] - 5
blast_df["near_end"] = np.where(start_condition | end_condition, True, False)
# blast_df["qrelpos_f"] = blast_hits[['sstart','send']].min(axis=1) / blast_hits["slen"]
# blast_df["qrelpos_r"]  = blast_hits["slen"] - blast_hits[['sstart','send']].max(axis=1) / blast_hits["slen"]
blast_df['blast_minloc'] = blast_hits[['sstart','send']].min(axis=1)
blast_df['blast_maxloc'] = blast_hits[['sstart','send']].max(axis=1)
blast_df['contig_base'] = blast_hits['sseqid'].str.split("|").str[-1].str.split("_") .str[0]

# Filter out things that are at the ends of contigs
blast_df.to_csv(os.path.join(path_output, "output_blast_hits-near_ends.csv"))
# blast_df = blast_df.loc[(blast_df["near_end"] == False)]

# Merge sample info with blast hits
sample_info = pd.read_csv(os.path.join(path_output, "sample_files_status.csv")).reset_index(drop=True)
blast_df = blast_df.merge(sample_info, how="left", on="contig_base")
blast_df = blast_df.reset_index(drop=True)

# Filter out things that don't have an associated fna or gff
blast_df = blast_df.loc[(blast_df["fna_file"] == True) & (blast_df["gff_file"] == True)].sort_values(by="sseqid")
blast_df.drop(columns=['fna_file','fsa_file','gff_file','tbl_file'], inplace=True)

grouped = blast_df.groupby("sample")
chunks = [group for _, group in grouped]

filtered_prokka_hits = pd.DataFrame()

#--------------------------------------------------------------------------------------------------------------------------------
count=0

for idx, chunk in enumerate(chunks):
    if idx % 100 == 0:
        print("On {} of {}".format(idx, len(chunks)))
    prokka_file = "prokka_out.gff"
    sample_path = list(chunk["path"])[0]
    full_path = os.path.join(sample_path, prokka_file)

    # Open and read the file line by line
    lines_to_keep = []
    with open(full_path, 'r') as file:
        keep = True
        for line in file:
            if not line.startswith('>'):
                if not line.startswith('#'):
                    if line.__contains__("Prodigal"):
                        lines_to_keep.append(line.strip())
            else:
                break  # Stop reading after the first line starting with '>'

    # Create a DataFrame from the lines to keep
    gff_fields = ["contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    prokka_labels = ["prokka_" + i for i in gff_fields]

    prokka_annos_df = pd.read_csv(io.StringIO('\n'.join(lines_to_keep)), sep="\t", header=None, names=prokka_labels)
    prokka_annos_df['prokka_minloc'] = prokka_annos_df[['prokka_start','prokka_end']].min(axis=1)
    prokka_annos_df['prokka_maxloc'] = prokka_annos_df[['prokka_start','prokka_end']].max(axis=1)
    prokka_annos_df.drop(columns=["prokka_source", "prokka_score", "prokka_type", "prokka_phase",'prokka_start','prokka_end'], inplace=True)
    prokka_annos_df = prokka_annos_df[(prokka_annos_df["prokka_contig"].isin(set(chunk['sseqid'])))] #grab entries from gff files that are contigs that had a blast hit

    # Extract rows from prokka gff files where there are annotations in the regions matched in the blast search.
    filt_prokka = pd.DataFrame(columns = prokka_annos_df.columns)

    for row in chunk.itertuples(index=False):
        count+=1

        sseqid = row.sseqid
        min_pos = row.blast_minloc
        max_pos = row.blast_maxloc

        prefilt_prokka = prokka_annos_df[(prokka_annos_df["prokka_contig"] == sseqid)]
        tmp_filt_prokka = pd.concat([
            prefilt_prokka.loc[(prefilt_prokka['prokka_minloc'] <= min_pos) & (prefilt_prokka['prokka_maxloc'] >= max_pos)], #extended left and right
            prefilt_prokka.loc[(prefilt_prokka['prokka_minloc'] >= min_pos) & (prefilt_prokka['prokka_maxloc'] <= max_pos)], #full length or truncated from left or right
            prefilt_prokka.loc[(prefilt_prokka['prokka_minloc'] <= min_pos) & (prefilt_prokka['prokka_maxloc'] > min_pos + 100)], #left shift
            prefilt_prokka.loc[(prefilt_prokka['prokka_maxloc'] >= max_pos) & (prefilt_prokka['prokka_minloc'] < max_pos - 100)]], #right shift
            ignore_index=True)
        tmp_filt_prokka.drop_duplicates(inplace=True)
        tmp_filt_prokka[[
        'blast_minloc', 'blast_maxloc', 'slen', 'qseqid', 'qlen', 'near_end', 'mismatch', 'pident', 'qcovhsp', 'sample'
        ]] = pd.DataFrame([[
        row.blast_minloc, row.blast_maxloc, row.slen, row.qseqid, row.qlen, row.near_end, row.mismatch, row.pident, row.qcovhsp, row.sample
        ]], index=tmp_filt_prokka.index)

        if len(tmp_filt_prokka) > 0:
            filt_prokka = pd.concat([filt_prokka, tmp_filt_prokka], ignore_index=True)
        else:
            tmp_filt_prokka = pd.DataFrame(columns=prefilt_prokka.columns)
            tmp_filt_prokka['prokka_contig'] = "no_prokka"
            tmp_filt_prokka[[
            'blast_minloc', 'blast_maxloc', 'slen', 'qseqid', 'qlen', 'near_end', 'mismatch', 'pident', 'qcovhsp', 'sample'
            ]] = pd.DataFrame([[
            row.blast_minloc, row.blast_maxloc, row.slen, row.qseqid, row.qlen, row.near_end, row.mismatch, row.pident, row.qcovhsp, row.sample
            ]], index=tmp_filt_prokka.index)
            filt_prokka = pd.concat([filt_prokka, tmp_filt_prokka], ignore_index=True)

    filtered_prokka_hits = pd.concat([filtered_prokka_hits, filt_prokka], ignore_index=True)

filtered_prokka_hits.reset_index(drop=True)
filtered_prokka_hits.to_csv(os.path.join(path_output, "filtered_prokka_hits.csv"))
#
# truncation_analysis = filtered_prokka_hits
# truncation_analysis[["category", "offset_minloc", "offset_maxloc"]] = pd.DataFrame([[np.nan, np.nan, np.nan]], index=truncation_analysis.index)
#
# truncation_analysis["offset_minloc"]= truncation_analysis.iloc[0]["blast_minloc"] - truncation_analysis.iloc[0]["prokka_minloc"]
# truncation_analysis["offset_maxloc"]= truncation_analysis.iloc[0]["prokka_maxloc"] - truncation_analysis.iloc[0]["blast_maxloc"]
#
# truncation_analysis.sort_values(by="qseqid")
# truncation_analysis.reset_index(drop=True)
# truncation_analysis.to_csv(os.path.join(path_output, "truncation_analysis_output.csv"))
#


#






#
