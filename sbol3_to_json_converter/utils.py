import sbol3
import xml.etree.ElementTree as ET
from collections import defaultdict, deque
from .mappings import NS, inverse_interaction_mapping, inverse_role_mapping, inverse_component_map, inverse_ed_mapping

def get_label_from_uri(uri):
  label = (
    inverse_interaction_mapping.get(uri) or
    inverse_role_mapping.get(uri) or
    inverse_component_map.get(uri) or
    inverse_ed_mapping.get(uri)
  )
  
  return label if label is not None else "Unknown"
    
def get_display_id_from_name(component_name, root):
    for desc in root.findall("rdf:Description", NS):
      name_element = desc.find("sbol:name", NS)
      if name_element is not None and name_element.text == component_name:
        
        type_element = desc.find("rdf:type", NS)
        if type_element is not None and type_element.attrib.get(f"{{{NS['rdf']}}}resource", "") == "http://sbols.org/v3#SubComponent":

          instance_of_element = desc.find("sbol:instanceOf", NS)
          if instance_of_element is not None:
            referenced_uri = instance_of_element.attrib.get(f"{{{NS['rdf']}}}resource", "")
            referenced_desc = root.find(f"rdf:Description[@rdf:about='{referenced_uri}']", NS)
            if referenced_desc is not None:
              ref_display_id = referenced_desc.find("sbol:displayId", NS)
              if ref_display_id is not None:
                return ref_display_id.text  

        display_id_element = desc.find("sbol:displayId", NS)
        if display_id_element is not None:
          return display_id_element.text 

    return None  

def sort_components_by_order(components, order):
  graph = defaultdict(list) 
  in_degree = {comp: 0 for comp in components} 

  for first, second in order:
    if first in components and second in components:
      graph[first].append(second)
      in_degree[second] += 1

  queue = deque([node for node in in_degree if in_degree[node] == 0])  
  sorted_order = []

  while queue:
    node = queue.popleft()
    sorted_order.append(node)

    for neighbor in graph[node]:
      in_degree[neighbor] -= 1
      if in_degree[neighbor] == 0:
        queue.append(neighbor)

  remaining_components = set(components) - set(sorted_order)
  sorted_order.extend(remaining_components)

  return sorted_order

def extract_unique_components(order):
  seen = set()
  unique_components = []

  for first, second in order:
    for component in (first, second):
      if component not in seen:
        seen.add(component)
        unique_components.append(component)

  return unique_components

def resolve_component_name(root, component_uri):
    if "ComponentReference" in component_uri:
        component_ref = root.find(f"rdf:Description[@rdf:about='{component_uri}']", NS)
        refers_to = component_ref.find("sbol:refersTo", NS) if component_ref is not None else None
        if refers_to is not None:
            referenced_component_id = refers_to.attrib.get(f"{{{NS['rdf']}}}resource", "")
            referenced_component = root.find(f"rdf:Description[@rdf:about='{referenced_component_id}']", NS)
            if referenced_component is not None:
                component_name = referenced_component.find("sbol:name", NS)
                display_id = referenced_component.find("sbol:displayId", NS)
                
                if component_name is not None:
                    component_name = get_display_id_from_name(component_name.text, root)
                    
                return component_name if component_name is not None else display_id.text if display_id is not None else "Unknown"

    elif "LocalSubComponent" in component_uri or "SubComponent" in component_uri:
        component = root.find(f"rdf:Description[@rdf:about='{component_uri}']", NS)
        if component is not None:
            component_name = component.find("sbol:name", NS)
            display_id = component.find("sbol:displayId", NS)
            
            if component_name is None:
                instanceof = component.find("sbol:instanceOf", NS) 
                if instanceof is not None:
                    instance_uri = instanceof.attrib.get(f"{{{NS['rdf']}}}resource", "")
                    instance_component = root.find(f"rdf:Description[@rdf:about='{instance_uri}']", NS)
                    if instance_component is not None:
                        display_id = instance_component.find("sbol:displayId", NS)  

            return component_name.text if component_name is not None else display_id.text if display_id is not None else "Unknown"
    
    elif "Component" in component_uri:
        component = root.find(f"rdf:Description[@rdf:about='{component_uri}']", NS)
        if component is not None:
            component_name = component.find("sbol:name", NS)
            display_id = component.find("sbol:displayId", NS)
            if component_name is not None:
                component_name = get_display_id_from_name(component_name.text, root)
                
            return component_name if component_name is not None else display_id.text if display_id is not None else "Unknown"

    return "Unknown"
