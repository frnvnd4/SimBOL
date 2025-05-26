import json
from .rdf_parser import parse_sbol
from .clean_json import clean_json
from google.colab import files

def process_sbol(sbol_file, nombre_base):
  print("Processing SBOL file...")

  output_data = parse_sbol(sbol_file)
  json_output_path = f"{nombre_base}_parsed_new.json"
  with open(json_output_path, "w", encoding="utf-8") as json_file:
    json.dump(output_data, json_file, indent=4, ensure_ascii=False)

  output_clean_json = f"{nombre_base}_clean.json"
  new_data = clean_json(json_output_path)
  with open(output_clean_json, "w", encoding="utf-8") as f:
    json.dump(new_data, f, indent=4, ensure_ascii=False)

  print(f"JSON file generated successfully: {output_clean_json}")
  files.download(output_clean_json)

  return output_clean_json
