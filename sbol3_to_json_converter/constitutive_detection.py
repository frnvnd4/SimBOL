from .utils import get_display_id_from_name
from .mappings import NS

def is_constitutive(promoter_id, root, display_id_to_uri):
  promoter_uri = display_id_to_uri.get(promoter_id, None) 
  if promoter_uri is None:
    return False  

  promoter_desc = root.find(f"rdf:Description[@rdf:about='{promoter_uri}']", NS)
  if promoter_desc is not None:
    role = promoter_desc.find("sbol:role", NS)
    if role is not None:
      role_uri = role.attrib.get(f"{{{NS['rdf']}}}resource", "")
      if "SO:0002050" in role_uri:
        return True
  return False

def has_constitutive_regulator(root):
  constitutive_operons = [] 
  for desc in root.findall("rdf:Description", NS):
    rdf_about = desc.attrib.get(f"{{{NS['rdf']}}}about", "")

    type_element = desc.find("rdf:type", NS)
    if type_element is None or \
      type_element.attrib.get(f"{{{NS['rdf']}}}resource", "") != "http://sbols.org/v3#LocalSubComponent":
      continue

    role_element = desc.find("sbol:role", NS)
    if role_element is None or "SO:0002050" not in role_element.attrib.get(f"{{{NS['rdf']}}}resource", ""):
      continue

    for constraint in root.findall("rdf:Description", NS):
      constraint_type = constraint.find("rdf:type", NS)
      if constraint_type is None or constraint_type.attrib.get(f"{{{NS['rdf']}}}resource", "") != "http://sbols.org/v3#Constraint":
        continue

      subject = constraint.find("sbol:subject", NS)
      if subject is not None and subject.attrib.get(f"{{{NS['rdf']}}}resource", "") == rdf_about:

        object_ref = constraint.find("sbol:object", NS)
        if object_ref is not None:
          object_uri = object_ref.attrib.get(f"{{{NS['rdf']}}}resource", "")

          object_desc = root.find(f"rdf:Description[@rdf:about='{object_uri}']", NS)
          if object_desc is not None:
            object_name = object_desc.find("sbol:name", NS)
            if object_name is not None:
              constitutive_operons.append(object_name.text) 
  
  return constitutive_operons  

def detect_constitutive_operons(hierarchy, root, display_id_to_uri):
  
    if isinstance(hierarchy, dict):
        for operon, components in hierarchy.items():
            hierarchy[operon] = {
                "components": components,
                "constitutive": False  
            }

    for operon, data in hierarchy.items():
        components = data["components"]
        constitutive = False

        for component in components:
            component_display_id = get_display_id_from_name(component, root)

            if component_display_id:
                if is_constitutive(component_display_id, root, display_id_to_uri):
                    constitutive = True
                    break
        
        hierarchy[operon]["constitutive"] = constitutive
        
    constitutive_operons = has_constitutive_regulator(root)
    for operon in constitutive_operons:
        if operon in hierarchy:
            hierarchy[operon]["constitutive"] = True

    return hierarchy
