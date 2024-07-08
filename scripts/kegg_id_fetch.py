import requests
import json
import pandas as pd
from pathlib import Path
from time import sleep
from logger import custom_logger

# Method for logger current states of the program.
logger = custom_logger(__name__)


def find_kegg_entry(data):
    """Extract the first KEGG entry from UniProtKB cross-reference data."""
    # Iterate through results and find the first KEGG database entry
    for result in data.get('results', []):
        for entry in result.get('uniProtKBCrossReferences', []):
            if entry['database'] == 'KEGG':
                return entry
    return None


def fetch_kegg_id(gene_symbol: str) -> str:
    logger.info(f"Getting KEGG ID for {gene_symbol}, with KEGG.")
    """
    Get KEGG ID for a given gene symbol by querying the KEGG API.

    Parameters:
    gene_symbol (str): The gene symbol for which to retrieve the KEGG ID.

    Returns:
    str: The KEGG ID associated with the gene symbol if found, otherwise an empty string.
    """
    url = f"http://rest.kegg.jp/find/genes/{gene_symbol}"
    try:
        response = requests.get(url)
        # Raises an HTTPError for the response returned a 4xx or 5xx status code
        response.raise_for_status()
        lines = response.text.strip().split('\n')
        for line in lines:
            if line.startswith('hsa:'):
                kegg_id, genes = line.split('\t')
                splitted_genes = genes.split(';')
                if len(splitted_genes) == 2:
                    gene_list = [gene.strip()
                                 for gene in splitted_genes[0].split(',')]
                    if gene_symbol in gene_list:
                        logger.result(f"Found KEGG ID: {kegg_id}")
                        return kegg_id
    except requests.HTTPError as e:
        logger.error(f"HTTP error occurred: {e}")
    except requests.RequestException as e:
        logger.error(f"Request error: {e}")
    logger.warning(f"No hit found with KEGG!")
    return "No hit"


def fetch_uniprot_kegg_id(gene_symbol):
    logger.info(f"Getting KEGG ID for {gene_symbol}, with Uniprot.")
    """Fetch UniProt ID for a human gene symbol and print its associated KEGG ID."""
    base_url = "https://rest.uniprot.org"
    search_params = {
        'query': f'gene:"{gene_symbol}" AND organism_id:9606 AND reviewed:true',
        'format': 'json',
        'size': 1
    }
    response = requests.get(
        f"{base_url}/uniprotkb/search", params=search_params)
    if response.status_code == 200:
        data = response.json()
        kegg_entry = find_kegg_entry(data)
        if kegg_entry:
            kegg_id_uniprot = kegg_entry.get("id", "No KEGG ID found")
            logger.result(f"Found KEGG ID: {kegg_id_uniprot}")
            return kegg_id_uniprot
        else:
            logger.warning(f"No hit found with Uniprot!")
            return "No hit"
    else:
        logger.error(f"Failed to retrieve data: HTTP {response.status_code}")
        logger.error("Response details:", response.text)


def fetch_human_kegg_id(gene_symbol):
    uniprot_id = fetch_uniprot_kegg_id(gene_symbol)
    kegg_id = fetch_kegg_id(gene_symbol)
    return uniprot_id, kegg_id


def get_kegg_id_multiple(row):
    """
    Get KEGG IDs for multiple gene symbols.
    """
    gene_symbol = row["Gene_Symbol"]
    filter_gene_symbol = gene_symbol.strip()
    uniprot_id, kegg_id = fetch_human_kegg_id(filter_gene_symbol)
    row["KEGG_ID_UNIPROT"] = uniprot_id
    row["KEGG_ID"] = kegg_id
    sleep(2)
    return row


def clean_and_filter_df(df):
    """
    Cleans and filters the DataFrame by:
    1. Selecting specific columns.
    2. Removing rows with any NaN values.
    3. Removing rows where any cell contains '---'.
    4. Excluding rows where 'Gene_Symbol' or 'Gene_Description' contains RNA-related terms.

    Parameters:
    df (pd.DataFrame): The input DataFrame.

    Returns:
    pd.DataFrame: The cleaned and filtered DataFrame.
    """
    # Select relevant columns and drop any rows with NaN values
    df_cleaned = df[["Gene_Accession",
                     "Gene_Symbol", "Gene_Description"]].dropna()

    # Remove rows containing '---' in any column
    df_cleaned = df_cleaned[~df_cleaned.isin(['---']).any(axis=1)]

    # Regex pattern to match RNA-related terms
    rna_pattern = r'\b(?:RNA|miRNA|mir|RN|Y_RNA)\b(?!\s+polymerase|\s+binding)'

# Filter the DataFrame to exclude RNA-related terms in 'Gene_Symbol' or 'Gene_Description'
    df_cleaned = df_cleaned[
        ~df_cleaned['Gene_Symbol'].str.contains(rna_pattern, case=False, na=False) &
        ~df_cleaned['Gene_Description'].str.contains(
            rna_pattern, case=False, na=False)
    ]

    return df_cleaned


if __name__ == '__main__':
    cwd = Path.cwd()
    files = (cwd / 'data').glob("Significant*")
    # files = [cwd / 'data' / 'Significant_PK_iRBCvsuRBC_spleen_analysis.xlsx']
    for file in files:
        excel_sheets = pd.ExcelFile(file)
        if "Analysis" in excel_sheets.sheet_names:
            logger.info(f"Processing {file}")
            df = pd.read_excel(file, sheet_name='Analysis')
            filtered_df = clean_and_filter_df(df)
            kegg_id_df = filtered_df.apply(get_kegg_id_multiple, axis=1)
            result = pd.merge(df, kegg_id_df, on=[
                "Gene_Accession", "Gene_Symbol", "Gene_Description"], how='outer')
            result.to_excel(
                f'{file.with_suffix("")}_extended.xlsx', index=False)
