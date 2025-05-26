import json
import xml.etree.ElementTree as ET
from .rdf_parser import parse_sbol
from .clean_json import clean_json

sbol_file = "../../rep_loica_model.xml"
output_data = parse_sbol(sbol_file)

json_output_path = f"{sbol_file}_parsed_new.json"
with open(json_output_path, "w", encoding="utf-8") as json_file:
    json.dump(output_data, json_file, indent=4, ensure_ascii=False)

output_clean_json = f"{sbol_file}_clean.json"

new_data = clean_json(json_output_path)

with open(output_clean_json, "w", encoding="utf-8") as f:
        json.dump(new_data, f, indent=4, ensure_ascii=False)
        
print("File generated: ", output_clean_json)
