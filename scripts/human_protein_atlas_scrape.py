import requests
import json


def get_ensamble_id(kegg_symbol):
    # Base URL for the Human Protein Atlas API
    api_url = "https://www.proteinatlas.org/api/search_download.php"

    # Specify query parameters
    params = {
        'search': kegg_symbol,
        'format': 'json',
        'columns': 'eg',
        'compress': 'no'
    }
    # Make the API request
    response = requests.get(api_url, params=params)

    # Check for successful request
    if response.status_code == 200:
        data = json.loads(response.text)
        ensamble_id = json.dumps(data[0]["Ensembl"], indent=4)
        ensamble_id = ensamble_id.replace('"', '')
        return ensamble_id
    else:
        print(f"Failed to fetch data: {response.status_code}")


def get_human_protein_atlas_entry(ensamble_id):
    api_url = f"https://www.proteinatlas.org/{ensamble_id}.json"
    response = requests.get(api_url)
    if response.status_code == 200:
        data = json.loads(response.text)
        return json.dumps(data, indent=4)
    else:
        print(f"Failed to fetch data: {response.status_code}")


if __name__ == '__main__':
    A = {}
    ensamble_id = get_ensamble_id("IL10RB")
    entry = get_human_protein_atlas_entry(ensamble_id)
    en_file = json.loads(entry)
    A["IL10RB"] = en_file
    print(json.dumps(A, indent=4))
