import json
import re

def is_uri(value):
    uri_pattern = re.compile(r'^(https?|ftp):\/\/[^\s/$.?#].[^\s]*$')
    return isinstance(value, str) and bool(uri_pattern.match(value))

def clean_json(input_file):
    with open(input_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    filtered_components = []
    for component in data.get("components", {}).values():
        cleaned_component = {
            key: value for key, value in component.items() if not is_uri(value)
        }
        filtered_components.append(cleaned_component)
    

    cleaned_ed = []
    for ed in data.get("ED", []):
        cleaned_ed.append({
            key: value for key, value in ed.items() if key not in {"id", "definition"}
        })
    
    cleaned_interactions = []
    for interaction in data.get("interactions", []):
        cleaned_interactions.append({
            key: value for key, value in interaction.items() if key != "id"
        })
    
    new_data = {
        "components": filtered_components,
        "ED": cleaned_ed, 
        "hierarchy": data.get("hierarchy", {}),  
        "interactions": cleaned_interactions,  
        #"order": data.get("order", {}),
        #"constraints": data.get("constraints", {}),
        "orientation": data.get("orientation", {}),
        "equivalence": data.get("equivalence", {}),
        "replacements": data.get("replacements", {}),
        "spatial_relations": data.get("spatial_relations", {})
        
    }

    return new_data