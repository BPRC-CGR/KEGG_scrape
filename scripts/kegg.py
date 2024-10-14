import requests
from Bio.KEGG.REST import kegg_get
from logger import custom_logger

logger = custom_logger(__name__)


class KEGGExtractor:
    def __init__(self):
        self.kegg_parsed_data = {}

    def fetch_kegg_data(self, kegg_id):
        """
        Fetches KEGG data for a given KEGG ID and parses it.
        """
        try:
            logger.info(f"Fetching KEGG data for {kegg_id}")
            response = kegg_get(kegg_id)
            data = response.read()
            gene_dict = self.parse_kegg_data(data)
            if gene_dict:
                self.kegg_parsed_data[kegg_id] = gene_dict
            return gene_dict
        except Exception as e:
            logger.error(f"Error fetching KEGG data for {kegg_id}: {e}")
            return None

    @staticmethod
    def parse_kegg_data(data):
        """
        Parses the raw KEGG data and returns structured information.
        """
        kegg_dict = {}
        current_key = ""
        for line in data.split("\n"):
            if line.startswith("///"):
                break
            if " " in line:
                key, value = line.split(" ", 1)
                key, value = key.strip(), value.strip()
                if key:
                    current_key = key
                if current_key in kegg_dict:
                    if isinstance(kegg_dict[current_key], list):
                        kegg_dict[current_key].append(value)
                    else:
                        kegg_dict[current_key] = [kegg_dict[current_key], value]
                else:
                    kegg_dict[current_key] = value
        return kegg_dict
