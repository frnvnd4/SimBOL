from .utils import sort_components_by_order, extract_unique_components, resolve_component_name
from collections import defaultdict, deque
from lxml import etree
from .mappings import NS

def build_region_hierarchy(root, components):
    regions = {}

    for desc in root.findall("rdf:Description", NS):
        rdf_type = desc.find("rdf:type", NS)

        if rdf_type is not None and rdf_type.attrib.get(f"{{{NS['rdf']}}}resource") == "http://sbols.org/v3#Component":
            display_id_elem = desc.find("sbol:displayId", NS)
            if display_id_elem is None:
                continue

            region_id = display_id_elem.text
            features = []

            for feat in desc.findall("sbol:hasFeature", NS):
                features.append(resolve_component_name(root, feat.attrib[f"{{{NS['rdf']}}}resource"]))

            if features:
                regions[region_id] = features

    return regions

def restructure_regions(regions):
    structured_regions = {}

    for region_id, components in regions.items():

        structured_components = []


        for component in components:

            if component in regions: 
                structured_components.append({component: restructure_subcomponents(regions[component], regions, structured_components)})
            else:  
                structured_components.append(component)

        structured_regions[region_id] = structured_components

    return structured_regions

def remove_duplicate_keys(regions):
    nested_keys = set()

    def find_nested_keys(obj):
        nonlocal nested_keys
        if isinstance(obj, dict):
            keys = obj.keys()
            nested_keys.update(keys)
            for value in obj.values():
                find_nested_keys(value)
        elif isinstance(obj, list):
            for item in obj:
                find_nested_keys(item)

    for main_key, main_value in regions.items():
        if isinstance(main_value, list):
            for item in main_value:
                find_nested_keys(item)

    new_regions = {}

    for main_key, main_value in regions.items():
        if main_key not in nested_keys:
            new_regions[main_key] = main_value

    return new_regions

def restructure_subcomponents(subcomponent_list, regions, parent_components):
    structured_subcomponents = []

    for subcomponent in subcomponent_list:
        if subcomponent in regions:
            subcomponent_structured = restructure_subcomponents(regions[subcomponent], regions, structured_subcomponents)

            if subcomponent_structured not in parent_components:
                structured_subcomponents.append({subcomponent: subcomponent_structured})
        else:
            structured_subcomponents.append(subcomponent)

    return structured_subcomponents

def restructure_hierarchy(root, hierarchy, order, components):
    operon_hierarchy = {}

    if not hierarchy:
        unique_components = extract_unique_components(order)
        sorted_order = sort_components_by_order(unique_components, order)
        resolved_order = [resolve_component_name(root, comp) for comp in sorted_order]

        regions = build_region_hierarchy(root, components)
        
        structured_regions = restructure_regions(regions)
        final = remove_duplicate_keys(structured_regions)

        resolved_hierarchy = {
            "Operon": [final]
        }

    else:
        for component, parents in hierarchy.items():
            for parent in parents:
                if parent not in operon_hierarchy:
                    operon_hierarchy[parent] = []
                operon_hierarchy[parent].append(component)

        resolved_hierarchy = {}

        for operon, components in operon_hierarchy.items():
            sorted_order = sort_components_by_order(components, order)
            resolved_order = [resolve_component_name(root, comp) for comp in sorted_order]
            resolved_operon = resolve_component_name(root, operon)

            resolved_hierarchy[resolved_operon] = resolved_order

    return resolved_hierarchy
