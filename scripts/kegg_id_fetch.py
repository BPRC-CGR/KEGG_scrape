import requests
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
from logger import custom_logger
import signal
import sys

# Logger for tracking the state of the program
logger = custom_logger(__name__)


class KeggIdFetcher:
    def __init__(self, max_workers=5):
        self.max_workers = max_workers
        self.keep_running = True
        signal.signal(signal.SIGINT, self.signal_handler)

    def signal_handler(self, sig, frame):
        """Handle Ctrl+C (KeyboardInterrupt) gracefully."""
        logger.warning(
            "KeyboardInterrupt received. Shutting down gracefully...")
        self.keep_running = False
        sys.exit(0)

    def find_kegg_entry(self, data):
        """Extract the first KEGG entry from UniProtKB cross-reference data."""
        for result in data.get('results', []):
            for entry in result.get('uniProtKBCrossReferences', []):
                if entry['database'] == 'KEGG':
                    return entry
        return None

    def fetch_kegg_id(self, gene_symbol: str, species: str) -> str:
        """Get KEGG ID for a given gene symbol by querying the KEGG API."""
        logger.info(f"Getting KEGG ID for {gene_symbol} in {species}.")
        url = f"http://rest.kegg.jp/find/genes/{gene_symbol}"
        try:
            response = requests.get(url)
            response.raise_for_status()
            lines = response.text.strip().split('\n')
            for line in lines:
                if line.startswith(f'{species}:'):
                    kegg_id, genes = line.split('\t')
                    if gene_symbol in genes.split(';')[0]:
                        logger.info(f"Found KEGG ID: {kegg_id}")
                        return kegg_id
        except requests.HTTPError as e:
            logger.error(f"HTTP error occurred: {e}")
        except requests.RequestException as e:
            logger.error(f"Request error: {e}")
        logger.warning(f"No KEGG hit found for {gene_symbol} in {species}.")
        return "No hit"

    def fetch_uniprot_kegg_id(self, gene_symbol, species_code):
        """Fetch UniProt ID for a gene symbol and species code using UniProt API."""
        logger.info(
            f"Fetching UniProt KEGG ID for {gene_symbol} in species {species_code}.")
        base_url = "https://rest.uniprot.org"
        search_params = {
            'query': f'gene:"{gene_symbol}" AND organism_id:{species_code} AND reviewed:true',
            'format': 'json',
            'size': 1
        }
        try:
            response = requests.get(
                f"{base_url}/uniprotkb/search", params=search_params)
            response.raise_for_status()
            data = response.json()
            kegg_entry = self.find_kegg_entry(data)
            if kegg_entry:
                kegg_id = kegg_entry.get("id", "No KEGG ID found")
                logger.info(f"Found UniProt KEGG ID: {kegg_id}")
                return kegg_id
            logger.warning(
                f"No UniProt hit found for {gene_symbol} in species {species_code}.")
        except requests.HTTPError as e:
            logger.error(f"HTTP error occurred: {e}")
        except requests.RequestException as e:
            logger.error(f"Request error: {e}")
        return "No hit"

    def fetch_kegg_ids(self, gene_symbol):
        """Fetch KEGG and UniProt IDs for human and chimp species."""
        human_uniprot_id = self.fetch_uniprot_kegg_id(gene_symbol, '9606')
        human_kegg_id = self.fetch_kegg_id(gene_symbol, 'hsa')
        chimp_uniprot_id = self.fetch_uniprot_kegg_id(gene_symbol, '9598')
        chimp_kegg_id = self.fetch_kegg_id(gene_symbol, 'ptr')
        return (human_uniprot_id, human_kegg_id, chimp_uniprot_id, chimp_kegg_id)

    def get_kegg_id_multiple(self, row):
        """Process a row to retrieve KEGG and UniProt IDs for each gene symbol."""
        gene_symbol = row["Gene_Symbol"].strip()
        human_uniprot_list, human_kegg_list = [], []
        chimp_uniprot_list, chimp_kegg_list = [], []

        if not self.keep_running:
            return row  # Gracefully exit if interrupted

        # Create ThreadPoolExecutor here and ensure it's shut down properly
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = []
            for symbol in gene_symbol.split("/"):
                symbol = symbol.partition(' ')[0].replace('-', '')
                futures.append(executor.submit(self.fetch_kegg_ids, symbol))

            for future in as_completed(futures):
                try:
                    if not self.keep_running:
                        break
                    human_uniprot_id, human_kegg_id, chimp_uniprot_id, chimp_kegg_id = future.result()
                    human_uniprot_list.append(human_uniprot_id)
                    human_kegg_list.append(human_kegg_id)
                    chimp_uniprot_list.append(chimp_uniprot_id)
                    chimp_kegg_list.append(chimp_kegg_id)
                except Exception as e:
                    logger.error(f"Error in thread execution: {e}")

        row["KEGG_ID_UNIPROT_HUMAN"] = "/".join(human_uniprot_list)
        row["KEGG_ID_HUMAN"] = "/".join(human_kegg_list)
        row["KEGG_ID_UNIPROT_CHIMP"] = "/".join(chimp_uniprot_list)
        row["KEGG_ID_CHIMP"] = "/".join(chimp_kegg_list)
        sleep(2)
        return row

    def clean_and_filter_df(self, df):
        """Clean and filter the DataFrame based on RNA and miRNA patterns."""
        df_cleaned = df[["Gene_Accession",
                         "Gene_Symbol", "Gene_Description"]].dropna()
        df_cleaned = df_cleaned[~df_cleaned.isin(['---']).any(axis=1)]

        rna_pattern = r'\b(?:RNA|miRNA|mir|RN|Y_RNA)\b(?!\s+polymerase|\s+binding)'
        df_cleaned = df_cleaned[
            ~df_cleaned['Gene_Symbol'].str.contains(rna_pattern, case=False, na=False) &
            ~df_cleaned['Gene_Description'].str.contains(
                rna_pattern, case=False, na=False)
        ]
        return df_cleaned


def process_files():
    """Main function to process Excel files and enrich data with KEGG/UniProt IDs."""
    cwd = Path.cwd()
    fetcher = KeggIdFetcher()
    files = (cwd / 'data').glob("Significant*")

    for file in files:
        output_file = Path(f'{file.with_suffix("")}_extended.xlsx')
        if not output_file.is_file():
            excel_sheets = pd.ExcelFile(file)
            if "Analysis" in excel_sheets.sheet_names:
                logger.info(f"Processing {file}")

                df = pd.read_excel(file, sheet_name='Analysis')
                filtered_df = fetcher.clean_and_filter_df(df)
                kegg_id_df = filtered_df.apply(
                    fetcher.get_kegg_id_multiple, axis=1)

                result = pd.merge(df, kegg_id_df, on=[
                                  "Gene_Accession", "Gene_Symbol", "Gene_Description"], how='outer')
                result.to_excel(output_file, index=False)
        else:
            logger.info(
                f"Analysis file {output_file} already exists! Skipping!")


if __name__ == '__main__':
    process_files()
