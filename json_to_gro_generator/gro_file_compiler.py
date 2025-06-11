import re
import random

# Counter for generating unique names for auxiliary genes related to Non-Covalent Binding
ncb_aux_gene_counter = 0

def normalize_signal_name(name):
  """
  Normalizes a signal name by removing spaces, hyphens, and underscores
  to create a valid GRO identifier.

  Args:
    name (str): The original signal name.

  Returns:
    str: The normalized signal name.
  """
  return re.sub(r'[\s\-_]', '', name)

def flatten_components(component_list, hierarchy_dict):
  """
  Flattens a hierarchical component list.
  This function recursively traverses a list of components, which can contain
  strings (component IDs) or dictionaries (representing nested structures),
  and returns a single flat list of component IDs.

  Args:
    component_list (list): A list of components, where items can be strings
                           (component IDs) or dictionaries (nested hierarchies).
    hierarchy_dict (dict): The overall hierarchy data, though seemingly not
                           directly used in this function's current logic for
                           flattening, it's passed perhaps for future extensions
                           or context.

  Returns:
    list: A flat list of component string IDs.
  """
  flat_list = []
  for item in component_list:
    if isinstance(item, str):
      flat_list.append(item)
    elif isinstance(item, dict):
      for sub_name, sub_content in item.items():
        components_to_flatten = []
        if isinstance(sub_content, dict) and "components" in sub_content:
          components_to_flatten = sub_content["components"]
        elif isinstance(sub_content, list):
          components_to_flatten = sub_content

        if components_to_flatten:
          flat_list.extend(flatten_components(components_to_flatten, hierarchy_dict))
        else:
          flat_list.append(sub_name)
  return flat_list

def find_cds_in_hierarchy(component_name_to_search, full_hierarchy_data, component_lookup_map):
  """
  Searches for a CDS (Coding Sequence) component under a given component within a hierarchical structure.

  Args:
    component_name_to_search (str): The displayId of the component under which to search for a CDS.
    full_hierarchy_data (dict): The complete hierarchy data for all operons/designs.
    component_lookup_map (dict): A dictionary mapping component displayIds to their component data.

  Returns:
    str or None: The displayId of the first CDS found under the specified component, or None if not found.
  """

  def get_hierarchy_structure_for_component(data_subtree, target_id):
    """
    Recursively searches for a specific component's hierarchical structure.

    Args:
      data_subtree (dict or list): The current part of the hierarchy being searched.
      target_id (str): The component ID whose structure is sought.

    Returns:
      dict or list or None: The hierarchical structure under target_id, or None.
    """
    if isinstance(data_subtree, dict):
      for key, value in data_subtree.items():
        if key == target_id:
          return value
        result = get_hierarchy_structure_for_component(value, target_id)
        if result is not None:
          return result
    elif isinstance(data_subtree, list):
      for list_item in data_subtree:
        result = get_hierarchy_structure_for_component(list_item, target_id)
        if result is not None:
          return result
    return None

  def traverse_hierarchy_for_cds(current_data_structure):
    """
    Traverses a given data structure (part of the hierarchy) to find a CDS.

    Args:
      current_data_structure (dict or list or str): The structure to traverse.

    Returns:
      str or None: The ID of a CDS if found, otherwise None.
    """
    if isinstance(current_data_structure, dict):
      for comp_id_key, children_structure in current_data_structure.items():
        comp_entry_data = component_lookup_map.get(comp_id_key)
        if comp_entry_data and comp_entry_data.get("role") == "CDS":
          return comp_id_key
        result = traverse_hierarchy_for_cds(children_structure)
        if result:
          return result
    elif isinstance(current_data_structure, list):
      for item_in_list in current_data_structure:
        if isinstance(item_in_list, str): 
          comp_id_str = item_in_list
          comp_entry_data = component_lookup_map.get(comp_id_str)
          if comp_entry_data and comp_entry_data.get("role") == "CDS":
            return comp_id_str
        else: # It's a nested structure
          result = traverse_hierarchy_for_cds(item_in_list)
          if result:
            return result
    elif isinstance(current_data_structure, str): # Base case: a single component ID string
      comp_id_str_single = current_data_structure
      comp_entry_data = component_lookup_map.get(comp_id_str_single)
      if comp_entry_data and comp_entry_data.get("role") == "CDS":
        return comp_id_str_single
    return None

  starting_search_structure = None
  for operon_key_name in full_hierarchy_data:
    structure_within_operon = get_hierarchy_structure_for_component(full_hierarchy_data[operon_key_name], component_name_to_search)
    if structure_within_operon is not None:
      starting_search_structure = structure_within_operon
      break

  if starting_search_structure is not None:
    return traverse_hierarchy_for_cds({component_name_to_search: starting_search_structure})
  return None

def get_template_for_participant(participant_id, ed_definitions_list, all_interactions, hierarchy_map, all_components_list):
  """
  Resolves a participant's name/ID to its corresponding GRO template name.
  This is typically the ID of a CDS for proteins, the original name for simple chemicals,
  or the participant_id itself as a fallback.

  Args:
    participant_id (str): The original name/ID of the participant.
    ed_definitions_list (list): List of "ED" (External Definitions) objects.
    all_interactions (list): List of all interaction objects.
    hierarchy_map (dict): The hierarchy data.
    all_components_list (list): List of all component objects.

  Returns:
    str: The resolved GRO template name for the participant.
  """
  component_lookup = {comp["displayId"]: comp for comp in all_components_list}
  ed_lookup = {ed_item.get("name", ""): ed_item for ed_item in ed_definitions_list}

  ed_entry = ed_lookup.get(participant_id)
  if ed_entry and ed_entry.get("type") == 'Simple chemical':
    return participant_id

  component_entry = component_lookup.get(participant_id)
  protein_product_name_to_find = None

  if (component_entry and component_entry.get("type") == "Protein") or \
     (ed_entry and ed_entry.get("type") == "Protein"):
    protein_product_name_to_find = participant_id
  if component_entry and component_entry.get("role") == "CDS":
    return participant_id 

  # Try to find the gene (Template) that produces the protein
  if protein_product_name_to_find:
    genetic_production_interaction = next(
      (inter for inter in all_interactions if inter.get('type') == 'Genetic Production' and
       any(p.get('participant') == protein_product_name_to_find and p.get('role') == 'Product' for p in inter.get('participants', []))),
      None
    )
    if genetic_production_interaction:
      template_id = next((p.get('participant') for p in genetic_production_interaction.get('participants', []) if p.get('role') == 'Template'), None)
      if template_id:
        template_component_data = component_lookup.get(template_id)
        if template_component_data and template_component_data.get("role") == "CDS":
          return template_id # Found the CDS directly

        # Try to find CDS within the hierarchy of the template
        found_cds_id = find_cds_in_hierarchy(template_id, hierarchy_map, component_lookup)
        if found_cds_id:
          return found_cds_id
        #print(f"WARNING [get_template_for_participant]: Could not find CDS for template '{template_id}' of product '{protein_product_name_to_find}'. Returning template name.")
        return template_id 

  #print(f"WARNING [get_template_for_participant]: Could not resolve '{participant_id}' to chemical or protein CDS. Returning original name.")
  return participant_id

def is_promoter_under_control_of_affected_component(affected_component_id, promoter_id_to_check, full_hierarchy_data):
  """
  Checks if a specific promoter (promoter_id_to_check) is hierarchically under
  a given affected_component_id (e.g., an operator or another DNA region).

  Args:
    affected_component_id (str): The ID of the component that is directly affected
                                 in an interaction (e.g., an operator).
    promoter_id_to_check (str): The ID of the promoter whose hierarchical relationship
                                is being checked.
    full_hierarchy_data (dict): The complete hierarchy data for all designs.

  Returns:
    bool: True if promoter_id_to_check is found under affected_component_id, False otherwise.
  """

  def recursive_search_for_promoter_in_structure(structure_of_affected, target_promoter_id):
    """Helper to recursively search for target_promoter_id within structure_of_affected."""
    if isinstance(structure_of_affected, dict):
      for key, child_structure in structure_of_affected.items():
        if key == target_promoter_id: return True
        if recursive_search_for_promoter_in_structure(child_structure, target_promoter_id): return True
    elif isinstance(structure_of_affected, list):
      for child_item in structure_of_affected:
        if isinstance(child_item, str) and child_item == target_promoter_id: return True
        if recursive_search_for_promoter_in_structure(child_item, target_promoter_id): return True
    elif isinstance(structure_of_affected, str):
      if structure_of_affected == target_promoter_id: return True
    return False

  def find_affected_then_search_promoter(current_hierarchy_level, target_affected_id, target_promoter_id):
    """Helper to first find target_affected_id, then search for target_promoter_id within its structure."""
    if isinstance(current_hierarchy_level, dict):
      for key, value_structure in current_hierarchy_level.items():
        if key == target_affected_id:
          if recursive_search_for_promoter_in_structure(value_structure, target_promoter_id): return True
        # Recursively search in deeper levels even if current key is not the target_affected_id
        if find_affected_then_search_promoter(value_structure, target_affected_id, target_promoter_id): return True
    elif isinstance(current_hierarchy_level, list):
      for item_structure in current_hierarchy_level:
        if find_affected_then_search_promoter(item_structure, target_affected_id, target_promoter_id): return True
    return False

  for operon_data_structure in full_hierarchy_data.values():
    if find_affected_then_search_promoter(operon_data_structure, affected_component_id, promoter_id_to_check):
      return True
  return False

def get_reactant_gro_name_for_ncb(reactant_original_name, ed_definitions_list, component_lookup_map,
                                  qs_actions_setup_map_ref, all_interactions_list, full_hierarchy_data):
  """
  Gets the GRO-compatible name for a reactant involved in Non-Covalent Binding (NCB).
  If the reactant is a chemical signal, it may map to a QS_Protein.
  Otherwise, it resolves to a CDS or its original name.
  Modifies qs_actions_setup_map_ref if a new QS_Protein is needed for a chemical.

  Args:
    reactant_original_name (str): The original name of the reactant.
    ed_definitions_list (list): List of "ED" (External Definitions) objects.
    component_lookup_map (dict): Map of component displayIds to component data.
    qs_actions_setup_map_ref (dict): Mutable map for QS action setups.
                                      Format: {normalized_chem_name: (QS_Protein_ID, config_type, original_chem_name)}
    all_interactions_list (list): List of all interaction objects.
    full_hierarchy_data (dict): The complete hierarchy data.


  Returns:
    str: The GRO-compatible name for the reactant.
  """
  ed_entry_data = next((ed_item for ed_item in ed_definitions_list if ed_item.get("name") == reactant_original_name), None)

  if ed_entry_data and ed_entry_data.get("type") == "Simple chemical":
    normalized_reactant_chem_id = normalize_signal_name(reactant_original_name)
    if normalized_reactant_chem_id in qs_actions_setup_map_ref:
      return qs_actions_setup_map_ref[normalized_reactant_chem_id][0] # Return the QS_Protein name
    else:
      # If this chemical is a reactant in NCB and wasn't a primary QS signal -> create a QS_Protein for it to act as a TF.
      qs_protein_for_reactant_chem = f"QS_{normalized_reactant_chem_id}"
      if normalized_reactant_chem_id not in qs_actions_setup_map_ref:
        qs_actions_setup_map_ref[normalized_reactant_chem_id] = (
          qs_protein_for_reactant_chem,
          "SENSING_CONFIG_FROM_UI_NCB_REACTANT", # Indicates origin for UI/parameter lookup
          reactant_original_name
        )
      return qs_protein_for_reactant_chem

  # If not a chemical, try to resolve to CDS or use protein name
  cds_id_name = get_template_for_participant(reactant_original_name, ed_definitions_list, all_interactions_list, full_hierarchy_data, list(component_lookup_map.values()))

  if cds_id_name and cds_id_name != reactant_original_name: 
    return cds_id_name
  elif component_lookup_map.get(reactant_original_name, {}).get('role') == 'CDS': # It was already a CDS ID
    return reactant_original_name
  elif ed_entry_data and ed_entry_data.get("type") == "Protein": # It's a protein from ED without direct CDS in components
    return reactant_original_name # Use the ED protein name

  #print(f"WARNING: Could not resolve reactant '{reactant_original_name}' to QS signal or CDS for NCB.")
  return reactant_original_name 

def extract_ncb_production_genes_and_actions(interactions_list, ed_definitions_list, components_list,
                                             hierarchy_map, qs_actions_map_reference):
  """
  Extracts auxiliary genes and action details for Non-Covalent Binding (NCB)
  interactions that result in the production/emission of a simple chemical.

  Modifies qs_actions_map_reference by adding entries for reactants that are
  chemicals and need to be represented as QS_Proteins (acting as TFs).

  Args:
    interactions_list (list): List of all interaction objects.
    ed_definitions_list (list): List of "ED" (External Definitions) objects.
    components_list (list): List of all component objects.
    hierarchy_map (dict): The hierarchy data.
    qs_actions_map_reference (dict): Mutable map for QS action setups, passed to
                                     get_reactant_gro_name_for_ncb.

  Returns:
    tuple: A tuple containing:
      - created_auxiliary_genes (list): List of dictionaries, each defining an auxiliary gene.
      - ncb_emission_action_details (list): List of dictionaries, each detailing an NCB emission for UI/actions.
  """
  global ncb_aux_gene_counter
  component_lookup = {comp["displayId"]: comp for comp in components_list}
  created_auxiliary_genes = []
  ncb_emission_action_details = []

  for interaction_item in interactions_list:
    if interaction_item.get("type") == "Non-Covalent Binding":
      reactant_original_ids = [p_data["participant"] for p_data in interaction_item["participants"] if p_data["role"] == "Reactant"]
      product_original_id = next((p_data["participant"] for p_data in interaction_item["participants"] if p_data["role"] == "Product"), None)

      if not product_original_id or not reactant_original_ids:
        continue

      # Check if the product is a simple chemical defined in ED
      product_ed_data = next((ed_item for ed_item in ed_definitions_list if ed_item.get("name") == product_original_id and ed_item.get("type") == "Simple chemical"), None)

      if product_ed_data: # If an NCB produces a chemical, create an auxiliary gene system
        ncb_aux_gene_counter += 1
        aux_protein_id = f"P_aux_NCB{ncb_aux_gene_counter}" # Protein produced by the aux gene
        normalized_product_id_for_gene = normalize_signal_name(product_original_id)
        # Gene name includes product and counter for uniqueness
        aux_gene_id = f"Operon_NCB_Emit_{normalized_product_id_for_gene}_{ncb_aux_gene_counter}"

        # Resolve reactant names to their GRO representation (TF names)
        reactant_tf_gro_ids = set()
        for r_id in reactant_original_ids:
          tf_gro_id = get_reactant_gro_name_for_ncb(r_id, ed_definitions_list, component_lookup, qs_actions_map_reference, interactions_list, hierarchy_map)
          if tf_gro_id:
            reactant_tf_gro_ids.add(f'"{tf_gro_id}"')

        if not reactant_tf_gro_ids:
          continue # Skip if no valid TFs resolved for reactants

        # Determine promoter logic based on number of TFs
        promoter_logic_function = '"AND"' if len(reactant_tf_gro_ids) > 1 else '"YES"'

        aux_gene_definition = {
          "name": f'"{aux_gene_id}"',
          "proteins": {f'"{aux_protein_id}"'},
          "promoter": {"function": promoter_logic_function, "transcription_factors": reactant_tf_gro_ids},
          "is_aux_ncb_gene": True, # Flag for special handling in GRO generation
          "fixed_params": { # Predefined, non-user-configurable parameters for these aux genes
            "act_time": 0.0, "act_var": 0.0, "deg_time": 0.0, "deg_var": 0.0,
            "toOn": 0.0, "toOff": 0.0, "noise_time": 0.0
          }
        }
        created_auxiliary_genes.append(aux_gene_definition)

        # Information for UI and subsequent action generation
        ui_and_action_info = {
          "type": "ncb_emission", # Identifies the type of action
          "condition_protein": aux_protein_id, # The auxiliary protein that triggers emission
          "emitted_signal_original_name": product_original_id,
          "emitted_signal_gro_id": normalize_signal_name(product_original_id),
          "reactants_original_names": reactant_original_ids,
          "aux_gene_name": aux_gene_id
        }
        ncb_emission_action_details.append(ui_and_action_info)

  return created_auxiliary_genes, ncb_emission_action_details

def extract_biochemical_reactions(interactions_list, ed_definitions_list, hierarchy_map, components_list,
                                  all_interactions_for_resolution=None):
  """
  Extracts biochemical reaction data, focusing on enzymatic conversions of
  a substrate signal (S1) to a product signal (S2) catalyzed by an enzyme (protein).

  Args:
    interactions_list (list): List of interaction objects to parse.
    ed_definitions_list (list): List of "ED" (External Definitions) objects.
    hierarchy_map (dict): The hierarchy data.
    components_list (list): List of all component objects.
    all_interactions_for_resolution (list, optional): A comprehensive list of interactions
        to use for resolving participant templates (e.g., finding CDS for proteins).
        Defaults to interactions_list if None.

  Returns:
    list: A list of dictionaries, where each dictionary represents a detected
          S1-to-S2 enzymatic conversion, containing details about the
          reactant signal, product signal, and enzyme.
  """
  identified_reactions = []
  component_lookup = {comp["displayId"]: comp for comp in components_list}
  ed_lookup_by_name = {ed_item.get("name", ""): ed_item for ed_item in ed_definitions_list}

  interactions_to_use_for_resolution = all_interactions_for_resolution if all_interactions_for_resolution is not None else interactions_list

  for interaction_item in interactions_list:
    if interaction_item.get("type") == "Biochemical Reaction":
      reactants_details = []
      products_details = []
      modifiers_details = []

      for p_data in interaction_item.get("participants", []):
        participant_original_id = p_data["participant"]
        gro_identifier = participant_original_id 
        is_participant_signal = False
        is_participant_protein_modifier = False

        ed_data = ed_lookup_by_name.get(participant_original_id)

        if ed_data and ed_data.get("type") == "Simple chemical":
          gro_identifier = normalize_signal_name(participant_original_id)
          is_participant_signal = True
        # For Modifiers (enzymes) and potentially other roles if they are proteins
        elif (ed_data and ed_data.get("type") == "Protein") or \
             (component_lookup.get(participant_original_id, {}).get("type") == "Protein"):
          # Resolve protein to its CDS/template if possible
          gro_identifier = get_template_for_participant(participant_original_id, ed_definitions_list, interactions_to_use_for_resolution, hierarchy_map, components_list)
          if p_data["role"] == "Modifier":
            is_participant_protein_modifier = True # Mark if it's a protein acting as a modifier

        participant_entry = {
          "original_name": participant_original_id,
          "gro_id": gro_identifier,
          "is_signal": is_participant_signal,
          "is_protein_modifier": is_participant_protein_modifier
        }

        if p_data["role"] == "Reactant":
          reactants_details.append(participant_entry)
        elif p_data["role"] == "Product":
          products_details.append(participant_entry)
        elif p_data["role"] == "Modifier":
          modifiers_details.append(participant_entry)

      # Specifically look for S1 (signal) + E (protein modifier) -> S2 (signal)
      if len(reactants_details) == 1 and reactants_details[0]["is_signal"] and \
         len(products_details) == 1 and products_details[0]["is_signal"] and \
         len(modifiers_details) == 1 and modifiers_details[0]["is_protein_modifier"]:

        identified_reactions.append({
          "type": "enzymatic_conversion_S1_to_S2", # Custom type for this specific pattern
          "reactant_signal": reactants_details[0],
          "product_signal": products_details[0],
          "enzyme": modifiers_details[0],
          "sbol_interaction_id": interaction_item.get("displayId", f"biochem_{reactants_details[0]['original_name']}_{products_details[0]['original_name']}")
        })
  return identified_reactions

def extract_genes_and_qs_actions(interactions_list, hierarchy_map, components_list, ed_definitions_list):
  """
  (Versión final y robusta)
  Extrae definiciones de genes, manejando correctamente la lógica de reguladores modulados.
  """
  genes_definitions = []
  component_lookup = {comp["displayId"]: comp for comp in components_list}
  ed_lookup = {ed_item.get("name", ""): ed_item for ed_item in ed_definitions_list}
  qs_action_setup_map = {}

  # Paso 1: Identificar qué químicos modulan a qué proteínas
  chemical_effects_on_protein_regulators = {} 
  for interaction in interactions_list:
    interaction_type = interaction.get("type")
    if interaction_type not in ["Inhibition", "Stimulation", "Control"]:
      continue
    
    participants = interaction.get("participants", [])
    actors = [p for p in participants if p["role"] in ["Inhibitor", "Stimulator", "Modifier"]]
    target_id = next((p["participant"] for p in participants if p["role"] in ["Inhibited", "Stimulated", "Modified"]), None)

    if not actors or not target_id or ed_lookup.get(target_id, {}).get("type") != "Protein":
      continue # Solo nos interesan interacciones que afectan a proteínas

    for actor in actors:
      actor_id = actor["participant"]
      if ed_lookup.get(actor_id, {}).get("type") == "Simple chemical":
        logic = "NOT" if interaction_type == "Inhibition" else "YES"
        chemical_effects_on_protein_regulators[target_id] = (actor_id, logic)

  # Paso 2: Construir las definiciones de los genes
  gene_counter = 0
  for operon_data in hierarchy_map.values():
    raw_elements_in_operon = operon_data.get("components", [])
    flattened_elements_list = flatten_components(raw_elements_in_operon, hierarchy_map)
    is_operon_constitutive = operon_data.get("constitutive", False)

    gene_groups_by_promoter = []
    i = 0
    while i < len(flattened_elements_list):
      element_id_current = flattened_elements_list[i]
      if component_lookup.get(element_id_current, {}).get('role') == 'Promoter':
        promoter_id_current = element_id_current
        # cds_list_for_this_promoter = [
        #     elem for elem in flattened_elements_list[i+1:] 
        #     if component_lookup.get(elem, {}).get('role') == 'CDS'
        # ]
        # next_promoter_index = next((idx for idx, elem in enumerate(flattened_elements_list[i+1:]) if component_lookup.get(elem, {}).get('role') == 'Promoter'), -1)
        # if next_promoter_index != -1:
        #     cds_list_for_this_promoter = cds_list_for_this_promoter[:next_promoter_index]
        search_slice = flattened_elements_list[i+1:]
        next_promoter_index_in_slice = next((idx for idx, elem in enumerate(search_slice) if component_lookup.get(elem, {}).get('role') == 'Promoter'), -1)
        if next_promoter_index_in_slice != -1:
            components_in_operon = search_slice[:next_promoter_index_in_slice]
        else:
            components_in_operon = search_slice
        cds_list_for_this_promoter = [
            elem for elem in components_in_operon
            if component_lookup.get(elem, {}).get('role') == 'CDS'
        ]
        if cds_list_for_this_promoter:
          gene_groups_by_promoter.append((promoter_id_current, cds_list_for_this_promoter))
        #i += len(cds_list_for_this_promoter)
      i += 1

    for promoter_id_val, cds_ids_list in gene_groups_by_promoter:
      gene_counter += 1
      gene_name_gro_formatted = f'"Operon_{gene_counter}"'
      gene_definition = {"name": gene_name_gro_formatted, "proteins": {f'"{cds_id}"' for cds_id in cds_ids_list}, "promoter": {"function": '"TRUE"' if is_operon_constitutive else '"UNKNOWN"', "transcription_factors": set()}}
      effective_regulators_for_promoter = {}
      promoter_logic_type = "AND" # Por defecto

      for reg_interaction in interactions_list:
        reg_interaction_type = reg_interaction.get("type")
        if reg_interaction_type not in ["Control", "Stimulation", "Inhibition"]:
          continue
        
        participants = reg_interaction.get("participants", [])
        affected_comp = next((p["participant"] for p in participants if p["role"] in ["Modified", "Stimulated", "Inhibited", "Template"]), None)
        
        if not (affected_comp == promoter_id_val or is_promoter_under_control_of_affected_component(affected_comp, promoter_id_val, hierarchy_map)):
            continue

        regulators = [p for p in participants if p["role"] in ["Modifier", "Stimulator", "Inhibitor"]]

        for regulator in regulators:
          regulator_id = regulator["participant"]
          base_logic = "NOT" if reg_interaction_type == "Inhibition" else "YES"

          if ed_lookup.get(regulator_id, {}).get("type") == "Simple chemical":
            # Caso 1: El regulador es un químico que afecta directamente al promotor
            norm_chem_id = normalize_signal_name(regulator_id)
            qs_id = f"QS_{norm_chem_id}"
            if norm_chem_id not in qs_action_setup_map:
              qs_action_setup_map[norm_chem_id] = (qs_id, "SENSING_CONFIG_FROM_UI", regulator_id)
            effective_regulators_for_promoter[qs_id] = base_logic
          else:
            # Caso 2: El regulador es una proteína
            protein_cds_id = get_template_for_participant(regulator_id, ed_definitions_list, interactions_list, hierarchy_map, components_list)
            effective_regulators_for_promoter[protein_cds_id] = base_logic
            
            # Y esta proteína, a su vez, es modulada por un químico
            if regulator_id in chemical_effects_on_protein_regulators:
              chem_modulator_id, chem_logic = chemical_effects_on_protein_regulators[regulator_id]
              
              # La lógica booleana (NOT A) OR B es la correcta para un represor inducible
              promoter_logic_type = "OR"
              
              # El efecto neto del químico es activación
              net_effect_logic = "YES" if (base_logic == "NOT" and chem_logic == "NOT") else "NOT"

              norm_chem_id = normalize_signal_name(chem_modulator_id)
              qs_id = f"QS_{norm_chem_id}"
              if norm_chem_id not in qs_action_setup_map:
                qs_action_setup_map[norm_chem_id] = (qs_id, "SENSING_CONFIG_FROM_UI", chem_modulator_id)
              effective_regulators_for_promoter[qs_id] = net_effect_logic
      
      # Ensamblaje final de la compuerta
      effective_regulators_list = list(effective_regulators_for_promoter.items())
      current_gene_promoter_tfs = gene_definition["promoter"]["transcription_factors"]
      
      if not effective_regulators_list:
        if not is_operon_constitutive:
          gene_definition["promoter"]["function"] = '"FALSE"'
      elif len(effective_regulators_list) == 1:
        tf_id, logic = effective_regulators_list[0]
        gene_definition["promoter"]["function"] = f'"{logic}"'
        current_gene_promoter_tfs.add(f'"{tf_id}"')
      else:
        # Usa el tipo de lógica determinado (AND por defecto, OR para casos complejos)
        gene_definition["promoter"]["function"] = f'"{promoter_logic_type}"'
        for tf_id, logic in effective_regulators_list:
          tf_entry = f'"-{tf_id}"' if logic == "NOT" else f'"{tf_id}"'
          current_gene_promoter_tfs.add(tf_entry)
          
      genes_definitions.append(gene_definition)

  return genes_definitions, qs_action_setup_map
  
def extract_signal_definitions(ed_definitions_list, user_signal_parameters):
  """
  Extracts signal definitions for GRO from External Definitions (ED) and user-provided parameters.

  Args:
    ed_definitions_list (list): List of "ED" objects, where signals are "Simple chemical".
    user_signal_parameters (dict): Dictionary of parameters for each signal, keyed by original signal name.
                                   Expected keys: "kdiff", "kdeg".

  Returns:
    list: A list of strings, each representing a GRO `s_signal` definition.
  """
  signal_definitions_for_gro = []
  for ed_item in ed_definitions_list:
    if ed_item.get("type") == "Simple chemical":
      original_signal_name = ed_item.get("name")
      normalized_signal_id = normalize_signal_name(original_signal_name)
      if original_signal_name and normalized_signal_id:
        params_for_signal = user_signal_parameters.get(original_signal_name, {"kdiff": 1.0, "kdeg": 0.1}) 
        diffusion_rate = params_for_signal.get("kdiff", 1.0)
        degradation_rate = params_for_signal.get("kdeg", 0.1)
        signal_definitions_for_gro.append(f'{normalized_signal_id} := s_signal([kdiff := {diffusion_rate}, kdeg := {degradation_rate}]);')
  return signal_definitions_for_gro

def generate_protein_paint_actions(protein_to_color_map):
  """
  Generates GRO "paint" actions for proteins based on a color mapping.

  Args:
    protein_to_color_map (dict): A dictionary mapping protein IDs (strings) to color names (strings).

  Returns:
    list: A list of strings, each a GRO "paint" action.
  """
  color_map = {"green": 0, "red": 1, "yellow": 2, "cyan": 3}
  actions = []
  for protein, color in protein_to_color_map.items():
    index = color_map.get(color, 0)
    color_array = ['"0"', '"0"', '"0"', '"0"']
    color_array[index] = f'"{random.randint(100, 200)}"'
    actions.append(f'action({{"{protein}"}}, "paint", {{{", ".join(color_array)}}});')
    color_array[index] = f'"{random.randint(-100, -50)}"'
    actions.append(f'action({{"-{protein}"}}, "paint", {{{", ".join(color_array)}}});')
  return actions

def generate_qs_signal_sensing_actions(qs_action_setup_map, user_signal_parameters):
  """
  Generates GRO `s_get_QS` actions based on the QS setup map and user-defined signal parameters.

  Args:
    qs_action_setup_map (dict): Map from `extract_genes_and_qs_actions`.
                                Format: {normalized_chem_name: (QS_Protein_ID, config_type, original_chem_name)}
    user_signal_parameters (dict): Parameters for signals, keyed by original signal name.
                                   Expected to contain "Symbol_getQS" and "Threshold_getQS".

  Returns:
    list: A list of strings, each representing an `s_get_QS` action for GRO.
  """
  sensing_actions = []
  for _normalized_chem_name, (qs_protein_id, config_type, original_signal_id) in qs_action_setup_map.items():
    if config_type in ["SENSING_CONFIG_FROM_UI", "SENSING_CONFIG_FROM_UI_NCB_REACTANT"]: 
      signal_ui_params = user_signal_parameters.get(original_signal_id, {})
      operator_symbol = signal_ui_params.get("Symbol_getQS", ">") 
      threshold_value = signal_ui_params.get("Threshold_getQS", "0.3") 
      signal_gro_id_tostring = f"tostring({normalize_signal_name(original_signal_id)})"
      sensing_actions.append(f'action({{}}, "s_get_QS", {{{signal_gro_id_tostring}, "{operator_symbol}", "{threshold_value}", "{qs_protein_id}"}});')
  return sensing_actions

def generate_gro_file(simulation_params, gene_definitions_list, qs_actions_map_data, signal_definitions_list, protein_paint_actions_map, biochemical_reactions_data_list, output_file_path):
  """
  Generates the content of a .gro file based on processed simulation parameters and biological constructs.

  Args:
    simulation_params (dict): Dictionary of global and component-specific parameters from the UI.
    gene_definitions_list (list): List of gene definitions (dictionaries) from extract_genes_and_qs_actions.
    qs_actions_map_data (dict): Map for QS actions from extract_genes_and_qs_actions.
    signal_definitions_list (list): List of s_signal definition strings from extract_signal_definitions.
    protein_paint_actions_map (dict): Map of protein IDs to color names for paint actions.
    biochemical_reactions_data_list (list): List of biochemical reaction data from extract_biochemical_reactions.
    output_file_path (str): Path to write the generated .gro file.
  """
  gro_content_lines = ["include gro\n"] 

  # Global Simulation Parameters 
  gro_content_lines.append("// Global Simulation Parameters")
  gro_content_lines.append(f'set ("dt", {simulation_params.get("dt", 0.1)} ); // Timestep in minutes')
  gro_content_lines.append(f'set ("population_max", {simulation_params.get("population_max", 20000)} ); // Maximum cell population\n')

  # Signal Definitions and Setup 
  if signal_definitions_list:
    gro_content_lines.append("// Signal Diffusion Parameters")
    gro_content_lines.append('set("signals", 1.0); // Enable signal module')
    gro_content_lines.append('set("signals_draw", 1.0); // Enable signal drawing')
    gro_content_lines.append('grid("continuous", "gro_original", 10, 10, 8); // Grid (type, diffusion_method, length, cell_size, neighborhood_depth)\n')
    gro_content_lines.append("// Signal Definitions (kdiff = diffusion rate, kdeg = degradation rate)")
    gro_content_lines.extend(signal_definitions_list)
    gro_content_lines.append("") 

  # Gene Definitions
  gro_content_lines.append("// Gene Definitions")
  for gene_data_item in gene_definitions_list:
    sorted_proteins = sorted(list(gene_data_item["proteins"]))
    protein_list_str = ", ".join(sorted_proteins)
    tf_list_str = ", ".join(gene_data_item["promoter"]["transcription_factors"]) 

    gene_gro_name = gene_data_item["name"] 
    promoter_function_str = gene_data_item["promoter"]["function"] 
    gene_name_key_for_params = gene_gro_name.strip('"')

    if gene_data_item.get("is_aux_ncb_gene"):
      gro_content_lines.append(f"// Auxiliary gene for Non-Covalent Binding induced emission: {gene_name_key_for_params}")
      fixed_gene_params = gene_data_item["fixed_params"]
      act_times_str = str(fixed_gene_params["act_time"])
      act_vars_str = str(fixed_gene_params["act_var"])
      deg_times_str = str(fixed_gene_params["deg_time"])
      deg_vars_str = str(fixed_gene_params["deg_var"])
      to_on_noise, to_off_noise, noise_t = fixed_gene_params["toOn"], fixed_gene_params["toOff"], fixed_gene_params["noise_time"]
    else:
      default_timings = {"act_time": 10.0, "act_var": 1.0, "deg_time": 20.0, "deg_var": 1.0,
                         "toOn": 0.0, "toOff": 0.0, "noise_time": 100.0}
      gene_timing_params = simulation_params.get("gene_parameters", {}).get(gene_name_key_for_params, {})

      act_times = gene_timing_params.get("act_times", [default_timings["act_time"]])
      act_vars = gene_timing_params.get("act_vars", [default_timings["act_var"]])
      deg_times = gene_timing_params.get("deg_times", [default_timings["deg_time"]])
      deg_vars = gene_timing_params.get("deg_vars", [default_timings["deg_var"]])

      act_times_str = ", ".join(map(str, act_times))
      act_vars_str = ", ".join(map(str, act_vars))
      deg_times_str = ", ".join(map(str, deg_times))
      deg_vars_str = ", ".join(map(str, deg_vars))

      to_on_noise = float(gene_timing_params.get("toOn", default_timings["toOn"]))
      to_off_noise = float(gene_timing_params.get("toOff", default_timings["toOff"]))
      noise_t = float(gene_timing_params.get("noise_time", default_timings["noise_time"]))

    # Ensamblar el bloque del gen con los nuevos strings
    gro_content_lines.extend([
      f'genes([',
      f'  name := {gene_gro_name},',
      f'  proteins := {{{protein_list_str}}},',
      f'  promoter := [function := {promoter_function_str},',
      f'    transcription_factors := {{{tf_list_str}}},',
      f'    noise := [toOn := {to_on_noise}, toOff := {to_off_noise}, noise_time := {noise_t}]',
      f'  ],',
      f'  prot_act_times := [times := {{{act_times_str}}}, variabilities := {{{act_vars_str}}}],',
      f'  prot_deg_times := [times := {{{deg_times_str}}}, variabilities := {{{deg_vars_str}}}]',
      f']);\n'
    ])

  #  Plasmid Definitions
  plasmid_definitions_map = simulation_params.get("plasmid_configuration", {}).get("defined_plasmids", {})
  if plasmid_definitions_map:
    gro_content_lines.append("// Plasmid Definitions")
    plasmid_entries_list = []
    for plasmid_id_str, genes_on_plasmid_list in plasmid_definitions_map.items():
      formatted_genes_on_plasmid_str = ", ".join(genes_on_plasmid_list)
      plasmid_entries_list.append(f'{plasmid_id_str} := {{{formatted_genes_on_plasmid_str}}}')

    if plasmid_entries_list:
        gro_content_lines.append(f'plasmids_genes([ {", ".join(plasmid_entries_list)} ]);\n')


  # Quorum Sensing Signal Detection Actions (s_get_QS)
  if qs_actions_map_data:
    qs_sensing_action_strings = generate_qs_signal_sensing_actions(qs_actions_map_data, simulation_params.get("signal_parameters", {}))
    if qs_sensing_action_strings:
      gro_content_lines.append("// Quorum Sensing Signal Detection Actions")
      gro_content_lines.extend(qs_sensing_action_strings)
      gro_content_lines.append("")

  # Non-Covalent Binding Product Emission Actions
  ncb_ui_config_params = simulation_params.get("ncb_emission_parameters", {})
  ncb_action_details_from_extraction = simulation_params.get("info_for_ncb_emission_actions", [])

  if ncb_action_details_from_extraction and ncb_ui_config_params:
    ncb_emission_action_strings = []
    for ncb_action_info_item in ncb_action_details_from_extraction:
      if ncb_action_info_item.get("type") != "ncb_emission": continue 

      triggering_protein_id = ncb_action_info_item["condition_protein"] # This is P_aux_NCBx
      original_emitted_signal_id = ncb_action_info_item["emitted_signal_original_name"]
      emitted_signal_gro_normalized_id = ncb_action_info_item["emitted_signal_gro_id"]
      emission_config = ncb_ui_config_params.get(original_emitted_signal_id)
      if emission_config:
        emission_concentration = emission_config.get("concentration", 1.0) 
        emission_behavior_type = emission_config.get("emission_type", "exact") 
        signal_id_for_gro_tostring = f"tostring({emitted_signal_gro_normalized_id})"
        action_str = f'action({{"{triggering_protein_id}"}}, "s_emit_signal", {{{signal_id_for_gro_tostring}, "{emission_concentration}", "{emission_behavior_type}"}});'
        ncb_emission_action_strings.append(action_str)
      else:
        ncb_emission_action_strings.append(f"// WARNING: UI parameters for NCB induced emission of '{original_emitted_signal_id}' (via {triggering_protein_id}) not found.")

    if ncb_emission_action_strings:
      gro_content_lines.append("// Non-Covalent Binding Product Emission Actions")
      gro_content_lines.extend(ncb_emission_action_strings)
      gro_content_lines.append("")

  # Biochemical Signal Conversion Actions (S1 + E -> S2) 
  biochem_ui_config_params = simulation_params.get("biochemical_conversion_parameters", {})
  if biochemical_reactions_data_list and biochem_ui_config_params:
    biochem_action_strings = []
    for reaction_data_item in biochemical_reactions_data_list:
      if reaction_data_item.get("type") == "enzymatic_conversion_S1_to_S2":
        enzyme_gro_id_val = reaction_data_item["enzyme"]["gro_id"] 
        reactant_signal_gro_id_val = reaction_data_item["reactant_signal"]["gro_id"] 
        product_signal_gro_id_val = reaction_data_item["product_signal"]["gro_id"]

        reaction_config_key = reaction_data_item.get("sbol_interaction_id",
                                            f"{reaction_data_item['reactant_signal']['original_name']}_to_{reaction_data_item['product_signal']['original_name']}_by_{reaction_data_item['enzyme']['original_name']}")
        reaction_ui_params = biochem_ui_config_params.get(reaction_config_key)

        if reaction_ui_params:
          conversion_rate = float(reaction_ui_params.get("rate", 0.0))
          absorption_behavior = reaction_ui_params.get("absorption_type", "area")
          emission_behavior = reaction_ui_params.get("emission_type", "area")

          if conversion_rate > 0.0: 
            # Absorb S1
            biochem_action_strings.append(f'action({{"{enzyme_gro_id_val}"}}, "s_absorb_signal", {{tostring({reactant_signal_gro_id_val}), "{conversion_rate}", "{absorption_behavior}"}});')
            # Emit S2
            biochem_action_strings.append(f'action({{"{enzyme_gro_id_val}"}}, "s_emit_signal", {{tostring({product_signal_gro_id_val}), "{conversion_rate}", "{emission_behavior}"}});')
        else:
          biochem_action_strings.append(f"// WARNING: UI parameters for biochemical reaction key '{reaction_config_key}' not found.")


    if biochem_action_strings:
      gro_content_lines.append("// Biochemical Signal Conversion Actions")
      gro_content_lines.extend(biochem_action_strings)
      gro_content_lines.append("")

  # Cell Painting Actions
  all_gene_protein_ids = {p_id.strip('"') for gene in gene_definitions_list for p_id in gene["proteins"]}
  valid_protein_paint_config = {p_id: color for p_id, color in protein_paint_actions_map.items() if p_id in all_gene_protein_ids}

  if valid_protein_paint_config:
    paint_action_strings = generate_protein_paint_actions(valid_protein_paint_config)
    if paint_action_strings:
      gro_content_lines.append("// Cell Painting Actions (based on protein presence)")
      gro_content_lines.extend(paint_action_strings)
      gro_content_lines.append("")

  # Bacterial Conjugation Actions
  conjugation_params = simulation_params.get("conjugation_parameters", {})
  if conjugation_params.get("enabled", False):
    conjugation_settings_map = conjugation_params.get("settings", {})
    if conjugation_settings_map:
      conjugation_action_strings = ["// Bacterial Conjugation Actions"]
      for plasmid_id_conj_str, conj_rate in conjugation_settings_map.items():
        conjugation_action_strings.append(f'action({{}}, "conjugate", {{"{plasmid_id_conj_str.strip()}", "{conj_rate}"}});')
      if len(conjugation_action_strings) > 1: # If any actions were added
          gro_content_lines.extend(conjugation_action_strings)
          gro_content_lines.append("")

  # Program Definition 
  gro_content_lines.append(f'program p() := {{\n  skip();\n}};\n')
  
  # Main Program: Initial Cell and Signal Setup
  main_program_block = ["program main() := {"]

  main_program_block.append(f'  set("ecoli_growth_rate", {simulation_params.get("growth_rate", 0.0346)});')
  if simulation_params.get("signal_parameters", {}) or simulation_params.get("initial_ecoli_populations"):
      main_program_block.append("")
    
  # Initial signal placements (one-time)
  initial_signal_placements_block = []
  # Constant signal emissions (every step)
  constant_signal_emissions_block = []

  user_signal_parameters_map = simulation_params.get("signal_parameters", {})
  for signal_original_id, signal_params_from_ui in user_signal_parameters_map.items():
    normalized_signal_gro_id = normalize_signal_name(signal_original_id)
    initial_points_list = signal_params_from_ui.get("initial_points", [])

    for point_config in initial_points_list:
      concentration_val = point_config.get("conc", 0.0)
      if concentration_val > 0: # Only place if concentration is positive
        coord_x = point_config.get("x", 0.0)
        coord_y = point_config.get("y", 0.0)
        is_constant_emission = point_config.get("constant_emission", False)
        set_signal_call = f's_set_signal({normalized_signal_gro_id}, {concentration_val}, {coord_x}, {coord_y});'

        if is_constant_emission:
          constant_signal_emissions_block.append(f'    {set_signal_call} // Constant emission from this point')
        else:
          initial_signal_placements_block.append(f'  {set_signal_call} // Initial one-time placement')

  if initial_signal_placements_block:
    main_program_block.extend(initial_signal_placements_block)
    if constant_signal_emissions_block or simulation_params.get("initial_ecoli_populations"):
        main_program_block.append("")

  if constant_signal_emissions_block:
    main_program_block.append("  // Constant Signal Emissions (every timestep from specified points)")
    main_program_block.append("  true:") 
    main_program_block.append("  {")
    main_program_block.extend(constant_signal_emissions_block)
    main_program_block.append("  }")
    if simulation_params.get("initial_ecoli_populations"): #
      main_program_block.append("")

  # Initial E. coli cell populations
  initial_ecoli_defs_list = simulation_params.get("initial_ecoli_populations", [])
  if not initial_ecoli_defs_list:
    #print("Warning: No initial E.coli populations defined by user. Creating a default population.")
    default_cell_count = 100; default_pos_x = 0.0; default_pos_y = 0.0; default_pop_radius = 100.0
    default_plasmids_for_ecoli_str = "{}" 
    if plasmid_definitions_map:
      all_defined_plasmid_ids_quoted = [f'"{p_id.strip()}"' for p_id in plasmid_definitions_map.keys()]
      if all_defined_plasmid_ids_quoted:
        default_plasmids_for_ecoli_str = f'{{{", ".join(all_defined_plasmid_ids_quoted)}}}'
    main_program_block.append(f'  c_ecolis({default_cell_count}, {default_pos_x}, {default_pos_y}, {default_pop_radius}, {default_plasmids_for_ecoli_str}, program ecoli_program());')
  else:
    main_program_block.append("  // Initial Cell Population(s) Setup")
    for pop_config_item in initial_ecoli_defs_list:
      num_cells_in_pop = pop_config_item.get("num_ecolis", 100)
      center_x_pos = pop_config_item.get("center_x", 0.0)
      center_y_pos = pop_config_item.get("center_y", 0.0)
      population_radius = pop_config_item.get("radius", 100.0)
      plasmids_for_this_pop_quoted = [f'"{p_id.strip()}"' for p_id in pop_config_item.get("plasmids", [])]
      plasmids_set_str = f'{{{", ".join(plasmids_for_this_pop_quoted)}}}' if plasmids_for_this_pop_quoted else "{}"

      main_program_block.append(f'  c_ecolis({num_cells_in_pop}, {center_x_pos}, {center_y_pos}, {population_radius}, {plasmids_set_str}, program p());')

  main_program_block.append("};") 
  gro_content_lines.extend(main_program_block)
  gro_content_lines.append("")

  # Write to .gro File
  with open(output_file_path, 'w', encoding='utf-8') as gro_file_handle:
    gro_file_handle.write("\n".join(gro_content_lines))

