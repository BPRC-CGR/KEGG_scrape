import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from kegg import KEGGExtractor
from human_protein_atlas import HPAExtractor
from logger import custom_logger

logger = custom_logger(__name__)


class DataProcessor:
    def __init__(self):
        self.kegg_extractor = KEGGExtractor()
        self.hpa_extractor = HPAExtractor()
        self.combined_data = {}
        self.convertions = {
            "BM": "bone+marrow",
            "AXILN": "lymph+node",
            "SPLEEN": "spleen",
            "LIVER": "liver",
        }

    def process_kegg_and_hpa_data(self, kegg_id, single_cell_type):
        """
        Fetches KEGG data for a given KEGG ID and then fetches HPA data using the gene symbol.
        """
        try:
            kegg_data = self.kegg_extractor.fetch_kegg_data(kegg_id)
            if not kegg_data:
                logger.error(f"Failed to fetch KEGG data for {kegg_id}")
                return

            gene_symbol = kegg_data.get('SYMBOL', None)
            if gene_symbol:
                gene_symbol = gene_symbol.split(',')[0]
            else:
                logger.warning(f"No gene symbol found for {kegg_id}")
                return

            logger.info(f"Gene symbol for {kegg_id} is {gene_symbol}")

            hpa_data = self.hpa_extractor.fetch_hpa_data(
                gene_symbol, single_cell_type)
            if not hpa_data:
                logger.error(f"Failed to fetch HPA data for {gene_symbol}")
                return

            self.combined_data[kegg_id] = {
                "kegg_data": kegg_data,
                "hpa_data": hpa_data
            }

            logger.info(
                f"Successfully processed KEGG and HPA data for {kegg_id}")

        except Exception as e:
            logger.error(
                f"Error processing KEGG and HPA data for {kegg_id}: {e}")

    def update_data(self, complete_date):
        self.combined_data.update(complete_date)

    def get_combined_data(self):
        """
        Returns the combined KEGG and HPA data.
        """
        return self.combined_data

    def get_conversions(self, key):
        return self.convertions.get(key.upper())


def get_cell_type_and_ids(excel_file, data_processor):
    """
    Extracts the single cell type and KEGG IDs from the provided Excel file.
    """
    key = excel_file.stem.split("_")[-3]
    single_cell_type = data_processor.get_conversions(key)
    df = pd.read_excel(excel_file)[["KEGG_ID_UNIPROT_HUMAN", "KEGG_ID_HUMAN"]]
    ids = pd.unique(
        df.dropna()
        .replace("No hit", pd.NA)
        .dropna()
        .values.ravel()
    )
    return single_cell_type, [x for i in ids for x in i.split("/")]


def get_json(path):
    if path.is_file():
        with open(path, "r") as f:
            return json.load(f)
    return {}


if __name__ == "__main__":
    data_processor = DataProcessor()
    cwd = Path.cwd()
    hsa_entries = get_json(cwd / "scrape.json")
    hsa_entries_keys = list(hsa_entries.keys())
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = []
        for excel_file in (cwd / "KEGG_scrape" / "data").glob("*extended.xlsx"):
            single_cell_type, kegg_ids = get_cell_type_and_ids(
                excel_file, data_processor)
            filtered_hsa_keys = list(set(kegg_ids) - set(hsa_entries_keys))
            if filtered_hsa_keys:
                for kegg_id in filtered_hsa_keys:
                    futures.append(executor.submit(
                        data_processor.process_kegg_and_hpa_data, kegg_id, single_cell_type))
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Error during processing: {e}")
    data_processor.update_data(hsa_entries)
    combined_data = data_processor.get_combined_data()
    with open("scrape.json", "w") as w:
        w.write(json.dumps(combined_data, indent=4))
