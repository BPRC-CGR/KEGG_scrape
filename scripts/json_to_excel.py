import pandas as pd
import json
from pathlib import Path

from logger import custom_logger

logger = custom_logger(__name__)


def flatten_with_full_path(d, parent_key='', sep='_'):
    """
    Recursively flattens a nested dictionary, preserving the full path of keys.
    """
    items = []
    for k, v in d.items():

        k_clean = str(k).replace(' ', '_').replace(':', '_').replace('&', 'and')
        new_key = f"{parent_key}{sep}{k_clean}" if parent_key else k_clean
        if isinstance(v, dict):
            items.extend(flatten_with_full_path(v, new_key, sep=sep).items())
        elif isinstance(v, list):

            list_str = ", ".join([str(item) for item in v])
            items.append((new_key, list_str))
        else:
            items.append((new_key, str(v)))
    return dict(items)


def insert_empty_rows_between_subgroups(df):
    """
    Inserts empty rows between subgroups in the DataFrame.
    Subgroups are determined by the attribute key prefixes at each level.
    """

    df['Attribute_Parts'] = df['Attribute'].str.split('_')

    rows = []
    previous_parts = []
    for idx, row in df.iterrows():
        current_parts = row['Attribute_Parts']

        diff_level = None
        for i, (prev_part, curr_part) in enumerate(zip(previous_parts, current_parts)):
            if prev_part != curr_part:
                diff_level = i
                break

        if diff_level is None:
            if len(current_parts) > len(previous_parts):
                diff_level = len(previous_parts)
            elif len(current_parts) < len(previous_parts):
                diff_level = len(current_parts)

        if diff_level is not None and rows:

            rows.append({'Attribute': '', 'Value': ''})

        rows.append({'Attribute': row['Attribute'], 'Value': row['Value']})

        previous_parts = current_parts

    final_df = pd.DataFrame(rows)
    return final_df


def create_excel_files_from_json(json_data, output_dir="output_excel_files"):
    """
    Generates an Excel file for each Excel entry in the JSON. Each KEGG ID
    within an entry gets its own sheet in the corresponding Excel file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for excel_file_name, gene_data in json_data.items():
        excel_file_path = output_dir / f"{excel_file_name}.xlsx"
        with pd.ExcelWriter(excel_file_path) as writer:
            for kegg_id, kegg_data in gene_data.items():
                combined_flattened_data = flatten_with_full_path(kegg_data)
                flattened_items = list(combined_flattened_data.items())
                flattened_df = pd.DataFrame(
                    flattened_items, columns=['Attribute', 'Value'])
                final_df = insert_empty_rows_between_subgroups(flattened_df)
                sheet_name = kegg_id.replace(":", "_").replace("/", "_")
                sheet_name = sheet_name[:31]
                final_df.to_excel(writer, sheet_name=sheet_name, index=False)
                logger.info(
                    f"Created sheet for KEGG ID: {kegg_id} in {excel_file_name}.xlsx")

    logger.info(f"Excel files created in {output_dir}")

def generate_excel(output_directory):
    with open("scrape.json", "r") as f:
        json_data = json.load(f)
    create_excel_files_from_json(json_data, output_directory)
