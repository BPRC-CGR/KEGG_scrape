import pandas as pd
import Bio.KEGG.REST as kegg
import json
import os

filtered_kegg = {}


def get_file():
    # Can filter for other files is necessary
    current_cwd = os.getcwd()
    kegg_file_path = os.listdir(os.path.join(current_cwd, "input"))
    kegg_file = os.path.join(current_cwd, "input", kegg_file_path[0])
    if kegg_file.endswith(".xlsx"):
        df = pd.read_excel(kegg_file)
    else:
        df = pd.read_csv(kegg_file, sep=None, engine="python")
    ids = df.iloc[:, 0].dropna()
    return ids


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


def add_info(id, kegg_variable, kegg_dict):
    if id in filtered_kegg:
        filtered_kegg[id].update({kegg_variable: kegg_dict})
    else:
        filtered_kegg[id] = {kegg_variable: kegg_dict}


def get_gene_symbol(kegg_info, id):
    gene_symbol_dict = kegg_info.get("SYMBOL", -1)
    if gene_symbol_dict != -1:
        symbol = [i.strip() for i in gene_symbol_dict.split(",")]
        add_info(id, "symbol", symbol)


def get_name(kegg_info, id):
    name_dict = kegg_info.get("NAME", -1)
    if name_dict != -1:
        description = get_info(name_dict)
        description = description.pop('(RefSeq)', None)
        add_info(id, "description", description)


def get_pathway(kegg_info, id):
    pathway_dict = kegg_info.get("PATHWAY", -1)
    if pathway_dict != -1:
        pathway = get_info(pathway_dict)
        add_info(id, "pathway", pathway)


def get_brite(kegg_info, id):
    brite_dict = kegg_info.get("BRITE", -1)
    header = ""
    brite = {}
    if brite_dict != -1:
        for element in brite_dict:
            if "hsa" in element and element[0].isalpha():
                header = element
                brite[header] = []
            else:
                brite[header].append(element)
            add_info(id, "brite", brite)


def write_kegg_data():
    with open("pathway.json", "w") as f:
        json.dump(filtered_kegg, f, indent=4)


def main():
    ids = get_file()
    for id in ids:
        if id.startswith("hsa"):
            print(f"Fetching {id}")
            response = kegg.kegg_get(id)
            data = response.read()
            gene_dict = parse_gene_data(data)
            get_gene_symbol(gene_dict, id)
            get_name(gene_dict, id)
            get_pathway(gene_dict, id)
            get_brite(gene_dict, id)


if __name__ == "__main__":
    main()
    write_kegg_data()
