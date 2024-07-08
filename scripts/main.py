from pathlib import Path
import pandas as pd
import Bio.KEGG.REST as kegg
import json
from human_protein_atlas_scrape import get_ensamble_id, get_human_protein_atlas_entry


class FileLoader:
    def get_file(self, file_path):
        kegg_file = Path(file_path)
        if not kegg_file.exists():
            raise FileNotFoundError(f"The file {file_path} does not exist.")

        if kegg_file.suffix == ".xlsx":
            df = pd.read_excel(kegg_file)
        else:
            df = pd.read_csv(kegg_file, sep=None, engine="python")
        ids = df.iloc[:, 0].dropna()
        return ids


class KeggDataParser:
    @staticmethod
    def parse_gene_data(data):
        kegg_dict = {}
        current_key = ""
        for line in data.split("\n"):
            if line.startswith("///"):
                break
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

    @staticmethod
    def get_info(sub_list):
        temp_dict = {}
        if isinstance(sub_list, str):
            sub_list = [sub_list]
        for element in sub_list:
            try:
                id, description = element.split(" ", 1)
                temp_dict[id] = description
            except ValueError as e:
                print(f"The following error was thrown {e}")
        return temp_dict


class InfoBuilder:
    def __init__(self):
        self.filtered_kegg = {}
        self.HPA = {}

    def add_info(self, id, kegg_variable, kegg_dict):
        if id in self.filtered_kegg:
            self.filtered_kegg[id].update({kegg_variable: kegg_dict})
        else:
            self.filtered_kegg[id] = {kegg_variable: kegg_dict}

    def get_gene_symbol(self, kegg_info, id):
        gene_symbol_dict = kegg_info.get("SYMBOL", -1)
        if gene_symbol_dict != -1:
            symbol = [i.strip() for i in gene_symbol_dict.split(",")][0]
            self.add_info(id, "symbol", symbol)

    def get_name(self, kegg_info, id):
        name_dict = kegg_info.get("NAME", -1)
        if name_dict != -1:
            description = KeggDataParser.get_info(name_dict)
            description = description.pop('(RefSeq)', None)
            self.add_info(id, "description", name_dict)

    def get_pathway(self, kegg_info, id):
        pathway_dict = kegg_info.get("PATHWAY", -1)
        if pathway_dict != -1:
            # Your logic for 'pathway' can be placed here
            self.add_info(id, "pathway", pathway_dict)

    def get_brite(self, kegg_info, id):
        brite_dict = kegg_info.get("BRITE", -1)
        if brite_dict != -1:
            # Your logic for 'brite' can be placed here
            self.add_info(id, "brite", brite_dict)


class DataWriter:
    @staticmethod
    def write_data_to_file(file_name, data, file_type="json"):
        if file_type == "json":
            with open(f"{file_name}.json", "w") as f:
                json.dump(data, f, indent=4)


class KeggAnalyzer:
    def __init__(self):
        self.file_loader = FileLoader()
        self.kegg_parser = KeggDataParser()
        self.info_builder = InfoBuilder()
        self.data_writer = DataWriter()

    def main(self, file_path):
        ids = self.file_loader.get_file(file_path)
        for id in ids:
            if id.startswith("hsa"):
                print(f"Fetching {id}")
                response = kegg.kegg_get(id)
                data = response.read()
                gene_dict = self.kegg_parser.parse_gene_data(data)

                self.info_builder.get_gene_symbol(gene_dict, id)
                self.info_builder.get_name(gene_dict, id)
                self.info_builder.get_pathway(gene_dict, id)
                self.info_builder.get_brite(gene_dict, id)

        self.data_writer.write_data_to_file(
            "kegg", self.info_builder.filtered_kegg)
        self.data_writer.write_data_to_file(
            "human_protein_atlas", self.info_builder.HPA)


if __name__ == "__main__":
    analyzer = KeggAnalyzer()
    analyzer.main()
