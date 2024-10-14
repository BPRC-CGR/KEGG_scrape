import argparse
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from kegg import KEGGExtractor
from human_protein_atlas import HPAExtractor

from logger import custom_logger
from json_to_excel import generate_excel
from kegg_id_fetch import process_files

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

    def process_kegg_and_hpa_data(self, kegg_id, single_cell_type, excel_file_name):
        """
        Fetches KEGG data for a given KEGG ID and then fetches HPA data using the gene symbol.
        Groups data under the corresponding Excel file name.
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

            # Group by the Excel file name in the combined data
            if excel_file_name not in self.combined_data:
                self.combined_data[excel_file_name] = {}

            self.combined_data[excel_file_name][kegg_id] = {
                "kegg_data": kegg_data,
                "hpa_data": hpa_data
            }

            logger.info(
                f"Successfully processed KEGG and HPA data for {kegg_id} from {excel_file_name}")

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
    return single_cell_type, [x for i in ids for x in i.split("/")], excel_file.stem


def get_json(path):
    if path.is_file():
        if path.stat().st_size == 0:  # Check if the file is empty
            logger.warning(f"{path} is empty. Returning an empty dictionary.")
            return {}
        try:
            with open(path, "r") as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            logger.error(f"Error decoding JSON from {path}: {e}")
            return {}
    else:
        logger.warning(f"{path} does not exist. Returning an empty dictionary.")
        return {}



def run_scrape(cwd, input_dir):
    data_processor = DataProcessor()
    hsa_entries = get_json(cwd / "scrape.json")
    hsa_entries_keys = list(hsa_entries.keys())
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = []
        for excel_file in (cwd / input_dir).glob("*extended.xlsx"):
            single_cell_type, kegg_ids, excel_file_name = get_cell_type_and_ids(
                excel_file, data_processor)
            filtered_hsa_keys = list(set(kegg_ids) - set(hsa_entries_keys))
            if filtered_hsa_keys:
                for kegg_id in filtered_hsa_keys:
                    futures.append(executor.submit(
                        data_processor.process_kegg_and_hpa_data, kegg_id, single_cell_type, excel_file_name))
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Error during processing: {e}")
    data_processor.update_data(hsa_entries)
    combined_data = data_processor.get_combined_data()
    return combined_data

    

def main():
    parser = argparse.ArgumentParser(description="Process KEGG and HPA data and generate Excel reports.")
    parser.add_argument("-i", "--input-dir", type=str, help="The directory containing input Excel files")
    parser.add_argument("-o", "--output-dir", type=str, default="output_data", help="The directory to save the output files (default: 'output_data')")

    args = parser.parse_args()
    
    
    cwd = Path.cwd()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    json_file = Path(cwd / "scrape.json")
    # Create the output directory if it doesn't exist
    if not json_file.is_file():
        json_file.touch(exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process KEGG and HPA data
    process_files(input_dir)

    combined_data = run_scrape(cwd, input_dir)
    
    # Save the combined data to the scrape.json file in the output directory
    with open(json_file, "w") as w:
        w.write(json.dumps(combined_data, indent=4))
    
    # Generate the Excel reports using the combined data
    generate_excel(output_dir)


if __name__ == "__main__":
    main()