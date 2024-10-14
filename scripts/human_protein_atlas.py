import json
import requests
import shutil
import time
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options
from logger import custom_logger

logger = custom_logger(__name__)


class HPAExtractor:
    def __init__(self):
        self.hpa_data = {}

    def get_geckodriver_path(self):
        geckodriver_path = shutil.which("geckodriver")
        if geckodriver_path is None:
            raise FileNotFoundError("GeckoDriver not found in PATH.")
        return geckodriver_path

    def search_gene_in_hpa(self, gene_id):
        """
        Search for the gene in Human Protein Atlas and find the row with the specific Ensembl ID.
        Returns the gene path if found, else None.
        """
        search_url = f"https://www.proteinatlas.org/search/{gene_id}"
        logger.info(f"Extracting information for gene ID: {gene_id}")

        try:
            response = requests.get(search_url)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            search_result_table = soup.find('table', class_='searchResult')

            if not search_result_table:
                logger.info(f"No search results found for gene ID: {gene_id}")
                return None

            for row in search_result_table.find_all('tr'):
                link = row.find('a', href=True)
                if link and link['href'].endswith(gene_id):
                    path = f"https://www.proteinatlas.org{link['href']}"
                    return path
            logger.info(f"No matching path found for gene ID: {gene_id}")
            return None

        except requests.RequestException as e:
            logger.error(f"Error fetching data for gene ID {gene_id}: {e}")
            return None

    def fetch_immune_info(self, gene_link):
        """
        Fetch immune information from the gene link.
        Extract the 'GENERAL INFORMATION' and 'IMMUNE CELL SECTION SUMMARY' tables.
        """
        immune_url = f"{gene_link}/immune+cell"
        try:
            response = requests.get(immune_url)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')

            general_info = self.extract_table_as_dict(
                soup, "GENERAL INFORMATION")
            immune_cell_summary = self.extract_table_as_dict(
                soup, "IMMUNE CELL SECTION SUMMARY")

            return {
                "general_info": general_info,
                "immune_cell_summary": immune_cell_summary
            }

        except requests.RequestException as e:
            logger.error(f"Error fetching immune data from HPA page: {e}")
            return None

    def extract_table_as_dict(self, soup, header_text):
        """
        Extracts the first table that contains a <th> element with class 'head' and the specified header_text.
        The key is the clean text inside <span> within <th>, and the value is the text inside <td>, excluding sub-elements like <sup>.
        """
        th_elements = soup.find_all('th', class_='head')

        for th in th_elements:
            if header_text in th.get_text():
                table = th.find_parent('table')
                if table:
                    table_data = {}
                    for row in table.find_all('tr'):
                        th = row.find('th')
                        td = row.find('td')
                        if th and td:
                            span = th.find('span')
                            if span:
                                for sub in span.find_all('sup'):
                                    sub.decompose()
                                key = span.get_text(strip=True)
                                value = td.get_text(strip=True)
                                table_data[key] = value
                    return table_data
        return None

    def extract_immune_rna_expression(self, soup):
        immune_cell_data = {}
        svg = soup.find('svg', class_='barchart')
        if svg:
            for g_tag in svg.find_all('g', class_='bar_g'):
                title = g_tag.get('title')
                if title:
                    parts = title.split("<br>")
                    cell_type = parts[0].replace(
                        "<b>", "").replace("</b>", "").strip()
                    nTPM = float(parts[1].replace("nTPM:", "").strip())
                    organ = parts[2].replace("Organ:", "").strip()
                    immune_cell_data[cell_type] = {
                        "cell_type": cell_type, "nTPM": nTPM, "organ": organ}
        return immune_cell_data if immune_cell_data else None

    def extract_single_cell_data(self, soup):
        single_cell_data = {}
        table = soup.find('table', class_='sc_legend_table')
        tbody = table.find('tbody', class_='hover') if table else None

        if tbody:
            for row in tbody.find_all('tr'):
                title = row.get('title')
                if title:
                    parts = title.split("<br>")
                    if len(parts) >= 7:
                        cluster = parts[0].strip()
                        cell_type = parts[1].strip()
                        organ = parts[2].strip()
                        reliability = parts[3].replace(
                            "Reliability:", "").strip()
                        cell_count = parts[4].replace("Cell count:", "").strip()
                        read_count = parts[5].replace("Read count:", "").strip()
                        expression = float(parts[6].replace(
                            "Expression:", "").strip().replace("nTPM", ""))

                        single_cell_data[cluster] = {
                            "cluster": cluster,
                            "cell_type": cell_type,
                            "organ": organ,
                            "reliability": reliability,
                            "cell_count": cell_count,
                            "read_count": read_count,
                            "expression": expression
                        }
        return single_cell_data if single_cell_data else None

    def fetch_data_with_selenium(self, url, data_extraction_function):
        try:
            geckodriver_path = self.get_geckodriver_path()

            options = Options()
            options.add_argument("--headless")

            driver = webdriver.Firefox(service=Service(
                executable_path=geckodriver_path), options=options)
            driver.get(url)
            time.sleep(2)

            soup = BeautifulSoup(driver.page_source, 'html.parser')
            driver.quit()

            return data_extraction_function(soup)

        except Exception as e:
            logger.error(f"An error occurred while processing {url}: {e}")
            return None

    def get_highest_and_lowest(self, data, key='nTPM'):
        if not data:
            return None, None
        lowest = min(data.values(), key=lambda x: x[key])
        highest = max(data.values(), key=lambda x: x[key])
        return lowest, highest

    def fetch_hpa_data(self, gene_id, single_cell_type="bone+marrow"):
        path = self.search_gene_in_hpa(gene_id)
        if not path:
            logger.error(f"Gene {gene_id} not found in HPA.")
            return None

        immune_info = self.fetch_immune_info(path)
        immune_rna_url = f"{path}/immune+cell"
        single_cell_url = f"{path}/single+cell+type/{single_cell_type}"

        immune_data = self.fetch_data_with_selenium(
            immune_rna_url, self.extract_immune_rna_expression)
        single_cell_data = self.fetch_data_with_selenium(
            single_cell_url, self.extract_single_cell_data)

        immune_lowest, immune_highest = self.get_highest_and_lowest(immune_data)
        single_cell_lowest, single_cell_highest = self.get_highest_and_lowest(
            single_cell_data, key='expression')

        structured_data = {
            "urls": {
                "immune_rna_url": immune_rna_url,
                "single_cell_url": single_cell_url
            },
            "immune": {
                "info": immune_info,
                "data": immune_data,
                "lowest": immune_lowest,
                "highest": immune_highest
            },
            "single_cell": {
                "data": single_cell_data,
                "lowest": single_cell_lowest,
                "highest": single_cell_highest
            }
        }

        return structured_data
