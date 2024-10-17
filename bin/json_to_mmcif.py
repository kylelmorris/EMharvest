import json
import argparse
from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from mmcif.io.PdbxWriter import PdbxWriter

prog = "JSON to mmCIF"
usage = """
        Converting JSON to mmCIF or adding more items from JSON to mmCIF
        Example for converting JSON:
        python json_to_mmcif.py -j emharvest/dep/dea57c9e-3b44-4250-8d95-750878c8ecb3_dep.json -f json 
        or f:or adding JSON to an existing mmcif information
        python json_to_mmcif.py -j emharvest/dep/dea57c9e-3b44-4250-8d95-750878c8ecb3_dep.json -f cif -c emharvest/dep_tomo/88415364-c070.cif
        """

parser = argparse.ArgumentParser(description="JSON to mmCIF")
parser.add_argument("-j", "--input_json_file", required=True, help="input json file to convert or add to mmcif")
parser.add_argument("-c", "--input_cif_file", help="input mmCIF file to add the json information to")
parser.add_argument("-f", "--input_format", choices=["json", "cif"], required=True,
                    help="json for converting directly from json, cif for adding the json file to the exisiting cif file")
args = parser.parse_args()

def main():
    convert_input_file(args.input_json_file, args.input_cif_file)

def json_to_dict(input_json_file):
    """Convert a JSON file to a Python dictionary"""
    try:
        with open(input_json_file, 'r') as file:
            data = json.load(file)
    except FileNotFoundError:
        print(f"Error: File '{input_json_file}' not found.")
    except json.JSONDecodeError:
        print(f"Error: File '{input_json_file}' is not a valid JSON file.")

    return data

def mmcif_to_json(input_cif_file):
    """
    Converts an mmCIF file to a JSON format.
    """
    json_data_dict = {}

    with open(input_cif_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('_'):
                key, value = line.split(maxsplit=1)
                if '.' in key:
                    category, sub_key = key[1:].split('.', 1)
                    if category not in json_data_dict:
                        json_data_dict[category] = {}
                    json_data_dict[category][sub_key] = value
    return json_data_dict

def write_mmcif_file(data_list):
    """
    pdbx writer is used to write data stored in self.__dataList
    :return written: a boolean; True when pdf writer is finished
    """
    written = False
    if args.input_json_file:
        mmcif_filename = args.input_json_file.split(".")[0] + '_conv.cif'
        with open(mmcif_filename, "w") as cfile:
            pdbx_writer = PdbxWriter(cfile)
            pdbx_writer.write(data_list)
        written = True
    return written


def add_container(data_list, container_id):
    """
    Adds a container with the specified container_id to the data_list.
    """
    container = DataContainer(container_id)
    data_list.append(container)
    return container


def add_category(container, category_id, items):
    """
    Adds a category with the specified category_id and items to the container.
    """
    category = DataCategory(category_id)
    for item in items:
        category.appendAttribute(item)
    container.append(category)


def insert_data(container, category_id, data_list):
    """
    Inserts data_list into the specified category_id within the container.
    """
    cat_obj = container.getObj(category_id)
    if cat_obj is None:
        return

    if all(isinstance(i, list) for i in data_list):
        list_values = [list(t) for t in zip(*data_list)]
        cat_obj.extend(list_values)
    else:
        cat_obj.append(data_list)

def convert_input_file(input_json_file, input_cif_file):
    container_dict = {}
    if args.input_format == "json":
        with open(input_json_file, 'r') as file:
            container_dict = json.load(file)

    if args.input_format == "cif":
        with open(input_json_file, 'r') as file:
            container_dict = json.load(file)

        cif_dict = mmcif_to_json(input_cif_file)

        for category, values in cif_dict.items():
            if category in container_dict:
                container_dict[category].update(values)
            else:
                container_dict[category] = values

    translate_json_to_cif(container_dict)

def translate_json_to_cif(container_dict):
    """
    Translates input json data into a CIF file.
    """
    cif_data_list = []

    # Add all accumulated categories and values to the CIF container
    container_id = args.input_json_file.split(".")[0]
    container = add_container(cif_data_list, container_id)

    for category_name, category_data in container_dict.items():
        category_list = list(category_data.keys())
        cif_values_list = list(category_data.values())

        add_category(container, category_name, category_list)
        insert_data(container, category_name, cif_values_list)

    # Write the modified CIF data to a file
    return write_mmcif_file(cif_data_list)

main()
