import xml.etree.ElementTree as ET
from .utils import get_label_from_uri, resolve_component_name
from .hierarchy import restructure_hierarchy
from .constitutive_detection import detect_constitutive_operons
from .mappings import NS

def parse_sbol(file):
  tree = ET.parse(file)
  root = tree.getroot()

  components = {}
  externally_defined = []
  interactions = []
  constraints = []
  hierarchy = {} 
  order = []  
  orientation = []
  equivalence = []
  replacements = []
  spatial_relations = []

  
  display_id_to_uri = {}

  extract_components(root, components, display_id_to_uri)
  resolve_component_references(root, components)
  extract_externally_defined(root, externally_defined)
  extract_interactions_and_participations(root, interactions, components, externally_defined)
  extract_constraints(root, constraints, hierarchy, order, orientation, equivalence, replacements, spatial_relations)

  hierarchy_order = restructure_hierarchy(root, hierarchy, order, components)
  hierarchy_order = detect_constitutive_operons(hierarchy_order, root, display_id_to_uri)
  
  output_data = { 
    "components": components,
    "ED": externally_defined,
    "hierarchy": hierarchy_order,
    "interactions": interactions,
    "orientation": orientation,
    "equivalence": equivalence,
    "replacements": replacements,
    "spatial_relations": spatial_relations,
    #"order": order,
    #"constraints": constraints,
  }

  return output_data

#extraer componentes
def extract_components(root, components, display_id_to_uri):
  for desc in root.findall("rdf:Description", NS):
    rdf_about = desc.attrib.get(f"{{{NS['rdf']}}}about", "")
    comp_type = desc.find("rdf:type", NS)
    
    if comp_type is not None and comp_type.attrib.get(f"{{{NS['rdf']}}}resource", "") == "http://sbols.org/v3#Component":
      display_id = desc.find("sbol:displayId", NS)
      role = desc.find("sbol:role", NS)
      sbol_type = desc.find("sbol:type", NS)
      
      display_id_text = display_id.text if display_id is not None else "Unknown"

      components[rdf_about] = {
        "displayId": display_id.text if display_id is not None else "Unknown",
        "role": get_label_from_uri(role.attrib.get(f"{{{NS['rdf']}}}resource", "Unknown") if role is not None else "Unknown"),
        "type": get_label_from_uri(sbol_type.attrib.get(f"{{{NS['rdf']}}}resource", "Unknown") if sbol_type is not None else "Unknown")
      }
      
      display_id_to_uri[display_id_text] = rdf_about 

# resolver referencias de componentes
def resolve_component_references(root, components):
  for desc in root.findall("rdf:Description", NS):
    rdf_about = desc.attrib.get(f"{{{NS['rdf']}}}about", "")
    if desc.find("rdf:type", NS) is not None and \
      desc.find("rdf:type", NS).attrib.get(f"{{{NS['rdf']}}}resource", "") == "http://sbols.org/v3#ComponentReference":

      ref = desc.find("sbol:refersTo", NS)
      if ref is not None:
        ref_target = ref.attrib.get(f"{{{NS['rdf']}}}resource", "")
        if ref_target in components:
          components[rdf_about] = components[ref_target]  

# extraer `ExternallyDefined`
def extract_externally_defined(root, externally_defined):
  for desc in root.findall("rdf:Description", NS):
    rdf_about = desc.attrib.get(f"{{{NS['rdf']}}}about", "")
    if desc.find("rdf:type", NS) is not None and desc.find("rdf:type", NS).attrib.get(f"{{{NS['rdf']}}}resource", "") == "http://sbols.org/v3#ExternallyDefined":
      display_id = desc.findtext("sbol:displayId", "Unknown", NS)
      name = desc.findtext("sbol:name", "Unknown", NS)
      sbol_type = desc.find("sbol:type", NS).attrib.get(f"{{{NS['rdf']}}}resource", "Unknown") if desc.find("sbol:type", NS) is not None else "Unknown"
      definition = desc.find("sbol:definition", NS).attrib.get(f"{{{NS['rdf']}}}resource", "Unknown") if desc.find("sbol:definition", NS) is not None else "Unknown"

      externally_defined.append({
        "id": rdf_about,
        "displayId": display_id,
        "name": name,
        "type": get_label_from_uri(sbol_type),
        "definition": definition
      })

def extract_interactions_and_participations(root, interactions, components, externally_defined):
    participations = {}
    
    for desc in root.findall("rdf:Description", NS):
        rdf_about = desc.attrib.get(f"{{{NS['rdf']}}}about", "")

        if desc.find("rdf:type", NS) is not None and \
            desc.find("rdf:type", NS).attrib.get(f"{{{NS['rdf']}}}resource", "") == "http://sbols.org/v3#Participation":
            
            participant_ref = desc.find("sbol:participant", NS)
            if participant_ref is not None:
                participant_uri = participant_ref.attrib.get(f"{{{NS['rdf']}}}resource", "")
                
                participant_name = resolve_participant_name(root, participant_uri, components, externally_defined)
                participations[rdf_about] = participant_name

    for desc in root.findall("rdf:Description", NS):
        rdf_about = desc.attrib.get(f"{{{NS['rdf']}}}about", "")

        if desc.find("rdf:type", NS) is not None and \
            desc.find("rdf:type", NS).attrib.get(f"{{{NS['rdf']}}}resource", "") == "http://sbols.org/v3#Interaction":
            
            interaction_type = desc.find("sbol:type", NS)
            interaction_type = interaction_type.attrib.get(f"{{{NS['rdf']}}}resource", "Unknown") if interaction_type is not None else "Unknown"
            
            interaction_participations = []
            
            for part in desc.findall("sbol:hasParticipation", NS):
                participation_id = part.attrib.get(f"{{{NS['rdf']}}}resource", "")
                
                participant_name = participations.get(participation_id, "Unknown")
                
                participation = root.find(f"rdf:Description[@rdf:about='{participation_id}']", NS)
                if participation is not None:
                 
                    instance_of = participation.find("sbol:instanceOf", NS)
                    if instance_of is not None:
                        instance_uri = instance_of.attrib.get(f"{{{NS['rdf']}}}resource", "")
                        participant_name = resolve_participant_name(root, instance_uri, components, externally_defined)

                    role = participation.find("sbol:role", NS)
                    role = role.attrib.get(f"{{{NS['rdf']}}}resource", "Unknown") if role is not None else "Unknown"

                    interaction_participations.append({
                        "role": get_label_from_uri(role),
                        "participant": participant_name
                    })

            interactions.append({
                "id": rdf_about,
                "type": get_label_from_uri(interaction_type),
                "participants": interaction_participations
            })

def resolve_participant_name(root, participant_uri, components, externally_defined):
    component_ref = root.find(f"rdf:Description[@rdf:about='{participant_uri}']", NS)

    if component_ref is not None:
        instance_of = component_ref.find("sbol:instanceOf", NS)
        if instance_of is not None:
            instance_uri = instance_of.attrib.get(f"{{{NS['rdf']}}}resource", "")
            return resolve_participant_name(root, instance_uri, components, externally_defined) 
          
    if "ComponentReference" in participant_uri:
  
        refers_to = component_ref.find("sbol:refersTo", NS)
        if refers_to is not None:
            referenced_component_id = refers_to.attrib.get(f"{{{NS['rdf']}}}resource", "")

            referenced_component = root.find(f"rdf:Description[@rdf:about='{referenced_component_id}']", NS)
            if referenced_component is not None:
                instance_of = referenced_component.find("sbol:instanceOf", NS)
                if instance_of is not None:
                    instance_uri = instance_of.attrib.get(f"{{{NS['rdf']}}}resource", "")
                    return components.get(instance_uri, {}).get("displayId", instance_uri)
        
            else:
                return "Unknown"
        else:
            return "Unknown"
    
    elif "ExternallyDefined" in participant_uri:

        for ext in externally_defined:
            if ext["id"] == participant_uri:
                return ext["name"]
        return "Unknown"

    else:
        participant_name = components.get(participant_uri, {}).get("displayId", participant_uri)

        return participant_name

# extraer constraints
def extract_constraints(root, constraints, hierarchy, order, orientation, equivalence, replacements, spatial_relations):
  for desc in root.findall("rdf:Description", NS):
    rdf_about = desc.attrib.get(f"{{{NS['rdf']}}}about", "")
    if desc.find("rdf:type", NS) is not None and \
      desc.find("rdf:type", NS).attrib.get(f"{{{NS['rdf']}}}resource", "") == "http://sbols.org/v3#Constraint":

      subject = desc.find("sbol:subject", NS)
      object_ = desc.find("sbol:object", NS)
      restriction = desc.find("sbol:restriction", NS)

      subject_ref = subject.attrib.get(f"{{{NS['rdf']}}}resource", "") if subject is not None else ""
      object_ref = object_.attrib.get(f"{{{NS['rdf']}}}resource", "") if object_ is not None else ""

      subject_name = resolve_component_name(root, subject_ref)
      object_name = resolve_component_name(root, object_ref)

      restriction_type = restriction.attrib.get(f"{{{NS['rdf']}}}resource", "Unknown") if restriction is not None else "Unknown"
      restriction_clean = restriction_type.lower()

      constraints.append({
        "id": rdf_about,
        "subject": subject_name,
        "object": object_name,
        "restriction": restriction_type
      })

      if any(k in restriction_clean for k in ["meets", "precedes", "strictlyprecedes", "finishes", "starts"]):
        order.append((subject_ref, object_ref))

      elif any(k in restriction_clean for k in ["contains", "strictlycontains", "covers"]):
        if object_ref not in hierarchy:
          hierarchy[object_ref] = []
        hierarchy[object_ref].append(subject_ref)

      elif "sameorientationas" in restriction_clean:
        orientation.append({
          "subject": subject_name,
          "object": object_name,
          "type": "sameOrientation"
        })

      elif "oppositeorientationas" in restriction_clean:
        orientation.append({
          "subject": subject_name,
          "object": object_name,
          "type": "oppositeOrientation"
        })

      elif "verifyidentical" in restriction_clean:
        equivalence.append({
          "subject": subject_name,
          "object": object_name,
          "type": "identical"
        })

      elif "differentfrom" in restriction_clean:
        equivalence.append({
          "subject": subject_name,
          "object": object_name,
          "type": "different"
        })

      elif "replaces" in restriction_clean:
        replacements.append({
          "old": object_name,
          "new": subject_name
        })

      elif "equals" in restriction_clean:
        equivalence.append({
          "subject": subject_name,
          "object": object_name,
          "type": "equals"
        })

      elif "isdisjointfrom" in restriction_clean:
        spatial_relations.append({
          "subject": subject_name,
          "object": object_name,
          "type": "disjoint"
        })

      elif "overlaps" in restriction_clean:
        spatial_relations.append({
          "subject": subject_name,
          "object": object_name,
          "type": "overlaps"
        })

