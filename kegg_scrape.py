import pandas as pd
import Bio.KEGG.REST as kegg
import requests
import json

test_dictionary = {}


def parse_gene_data(data):
    gene_dict = {}
    current_key = ""
    for line in data.split("\n"):
        if line.startswith("///"):
            break
        key, value = line.split(" ", 1)
        key = key.strip()
        value = value.strip()
        if key:
            current_key = key
        if current_key in gene_dict:
            if isinstance(gene_dict[current_key], list):
                gene_dict[current_key].append(value)
            else:
                gene_dict[current_key] = [gene_dict[current_key], value]
        else:
            gene_dict[current_key] = value
    return gene_dict


def get_info(sub_list):
    pathways = {}
    
    if isinstance(sub_list, str):
        sub_list = [sub_list]
        
    for element in sub_list:
        try:
            id, description = element.split(" ", 1)
            pathways[id] = description
        except ValueError as e:
            print(f"The following error was thrown {e}") 
    return pathways



def get_gene_symbol(kegg_info, id):
    gene_symbol_dict = kegg_info.get("SYMBOL", -1)
    if gene_symbol_dict != -1:
        symbol = [i.strip() for i in gene_symbol_dict.split(",")]
        if id in test_dictionary:
            test_dictionary[id].update({"symbol": symbol})
        else:
            test_dictionary[id] = {"symbol": symbol}

def get_name(kegg_info, id):
    name_dict = kegg_info.get("NAME", -1)
    if name_dict != -1:
        description = get_info(name_dict)
        description = description.pop('(RefSeq)', None)
        if id in test_dictionary:
            test_dictionary[id].update({"description": description})
        else:
            test_dictionary[id] = {"description": description}

def get_pathway(kegg_info, id):
    pathway_dict = kegg_info.get("PATHWAY", -1)
    if pathway_dict != -1:
        pathway = get_info(pathway_dict)
        if id in test_dictionary:
            test_dictionary[id].update({"pathway": pathway})
        else:
            test_dictionary[id] = {"pathway": pathway}


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
            if id in test_dictionary:
                test_dictionary[id].update({"brite": brite})
            else:
                test_dictionary[id] = {"brite": brite} 
                

if __name__ == "__main__":
    df = pd.read_excel("Copy of KEGG_BRITE_AxLN-spleen.xlsx")
    ids = df.iloc[:, 0].dropna()
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
            
        
    with open("pathway.json", "w") as f:
        json.dump(test_dictionary, f, indent=4)