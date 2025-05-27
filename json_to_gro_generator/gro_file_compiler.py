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
  # Removes spaces and hyphens to create a valid GRO identifier
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
    # Pass a dictionary to traverse_hierarchy_for_cds, with component_name_to_search as the key
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

def extract_degradation_actions(interactions_list, ed_definitions_list, hierarchy_map, components_list):
  """
  Extracts degradation actions from interactions.
  Note: This function is defined but its output is not directly used within this script's
  `generate_gro_file` flow. It's assumed to be called externally if its results are needed.

  Args:
    interactions_list (list): List of all interaction objects.
    ed_definitions_list (list): List of "ED" (External Definitions) objects.
    hierarchy_map (dict): The hierarchy data.
    components_list (list): List of all component objects.

  Returns:
    list: A list of tuples, where each tuple represents a degradation action:
          (degrader_gro_id, target_gro_id, target_type_string).
  """
  degradation_actions_tuples = []
  component_id_to_type_map = {c["displayId"]: c["type"] for c in components_list}
  ed_id_to_type_map = {e["name"]: e["type"] for e in ed_definitions_list}

  for interaction_item in interactions_list:
    if interaction_item["type"] == "Degradation":
      reactant_ids = []
      degrader_ids = []
      for p_data in interaction_item["participants"]:
        participant_role = p_data["role"]
        participant_id = p_data["participant"]
        if participant_role == "Reactant":
          reactant_ids.append(participant_id)
        else:
          degrader_ids.append(participant_id)

      for degrader_id_val in degrader_ids:
        for target_id_val in reactant_ids:
          degrader_gro_id = get_template_for_participant(degrader_id_val, ed_definitions_list, interactions_list, hierarchy_map, components_list)
          target_gro_id = get_template_for_participant(target_id_val, ed_definitions_list, interactions_list, hierarchy_map, components_list)
          target_type_str = component_id_to_type_map.get(target_gro_id) or \
                              ed_id_to_type_map.get(target_gro_id) or \
                              component_id_to_type_map.get(target_id_val) or \
                              ed_id_to_type_map.get(target_id_val) or \
                              "Unknown"

          degradation_actions_tuples.append((degrader_gro_id, target_gro_id, target_type_str))
  return degradation_actions_tuples

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
      else:
        #print(f"DEBUG: Biochemical reaction not matching S1->S2 pattern: {interaction_item.get('displayId')}")
        #print(f"DEBUG: Reactants: {reactants_details}, Products: {products_details}, Modifiers: {modifiers_details}")
  return identified_reactions

def extract_genes_and_qs_actions(interactions_list, hierarchy_map, components_list, ed_definitions_list):
  """
  Extracts gene definitions for GRO and sets up a map for Quorum Sensing (QS) actions.
  This function processes interactions to determine gene regulation (promoters, TFs)
  and identifies chemicals that act as effectors (inducing/repressing), mapping them
  to "QS_Proteins" that represent their activity in the GRO model.

  Args:
    interactions_list (list): List of all interaction objects.
    hierarchy_map (dict): The hierarchy data, showing component organization.
    components_list (list): List of all component objects (DNA parts, proteins etc.).
    ed_definitions_list (list): List of "ED" (External Definitions) for chemicals, etc.

  Returns:
    tuple: A tuple containing:
      - genes_definitions (list): List of dictionaries, each defining a gene for GRO.
      - qs_action_setup_map (dict): Map where keys are normalized chemical names and
                                    values are (QS_Protein_ID, config_type, original_chem_name).
                                    This map is used to generate s_get_QS actions.
  """
  genes_definitions = []
  component_lookup = {comp["displayId"]: comp for comp in components_list}
  ed_lookup = {ed_item.get("name", ""): ed_item for ed_item in ed_definitions_list}

  # To store how chemicals affect regulatory proteins (TFs)
  chemical_effects_on_protein_regulators = {} # protein_id -> (qs_protein_effector, logic_of_effect)
  # To store chemicals that directly induce/repress DNA elements (promoters/operators)
  direct_chemical_inducer_to_qs_protein_map = {} # chemical_original_name -> qs_protein_effector
  # Master map for all chemicals that need a QS_Protein representation
  qs_action_setup_map = {} # normalized_chem_name -> (QS_Protein_ID, config_type, original_chem_name)

  # --- Step 1: Pre-process interactions to identify chemical effectors and their targets ---
  for interaction_item in interactions_list:
    interaction_type = interaction_item.get("type")
    # Focus on regulatory interactions: Inhibition, Stimulation, Control
    if interaction_type not in ["Inhibition", "Stimulation", "Control"]:
      continue

    participants_data = interaction_item.get("participants", [])
    # Actor: Inhibitor, Stimulator, Modifier
    actor_id = next((p_data["participant"] for p_data in participants_data if p_data["role"] in ["Inhibitor", "Stimulator", "Modifier"]), None)
    # Target: Inhibited, Stimulated, Modified, Template (for genetic production context)
    target_id = next((p_data["participant"] for p_data in participants_data if p_data["role"] in ["Inhibited", "Stimulated", "Modified", "Template"]), None)

    if not actor_id or not target_id:
      continue # Skip if actor or target is missing

    # Check if the actor is a simple chemical from ED
    actor_is_ed_chemical = ed_lookup.get(actor_id, {}).get("type") == "Simple chemical"

    if actor_is_ed_chemical:
      normalized_actor_chem_id = normalize_signal_name(actor_id)
      # Create a unique QS_Protein name for this chemical effector
      qs_protein_for_actor_chemical = f"QS_{normalized_actor_chem_id}"

      # Store this mapping for later s_get_QS action generation
      if normalized_actor_chem_id not in qs_action_setup_map:
        qs_action_setup_map[normalized_actor_chem_id] = (qs_protein_for_actor_chemical, "SENSING_CONFIG_FROM_UI", actor_id)

      target_component_data = component_lookup.get(target_id)
      # Case 1: Chemical directly affects a DNA component (promoter, operator, etc.)
      target_is_dna_regulatory_component = target_component_data and \
                                        (target_component_data.get("type") == "DNA" or
                                         target_component_data.get("role") in ["Promoter", "Operator", "Engineered-Region"])

      # Case 2: Chemical affects a protein TF
      target_is_protein_tf = (component_lookup.get(target_id, {}).get("type") == "Protein") or \
                             (ed_lookup.get(target_id) and ed_lookup.get(target_id, {}).get("type") == "Protein")

      if target_is_dna_regulatory_component:
        direct_chemical_inducer_to_qs_protein_map[actor_id] = qs_protein_for_actor_chemical
      elif target_is_protein_tf:
        # Determine the logic of the chemical's effect on the protein TF
        logic_chemical_on_tf = "NOT" if interaction_type == "Inhibition" else "YES"
        chemical_effects_on_protein_regulators[target_id] = (qs_protein_for_actor_chemical, logic_chemical_on_tf)

  # --- Step 2: Create Gene Definitions for GRO based on hierarchy and processed regulatory info ---
  gene_counter = 0
  for operon_data in hierarchy_map.values(): # Iterate through each operon/design in the hierarchy
    raw_elements_in_operon = operon_data.get("components", [])
    # Flatten the component list for this operon for easier processing
    flattened_elements_list = flatten_components(raw_elements_in_operon, hierarchy_map)
    is_operon_constitutive = operon_data.get("constitutive", False)

    # Group CDSs by their preceding promoters
    gene_groups_by_promoter = []
    i = 0
    while i < len(flattened_elements_list):
      element_id_current = flattened_elements_list[i]
      if component_lookup.get(element_id_current, {}).get('role') == 'Promoter':
        promoter_id_current = element_id_current
        cds_list_for_this_promoter = []
        j = i + 1
        # Collect all subsequent CDSs until another promoter is found or list ends
        while j < len(flattened_elements_list):
          next_element_id_in_list = flattened_elements_list[j]
          next_element_role_str = component_lookup.get(next_element_id_in_list, {}).get('role')
          if next_element_role_str == 'CDS':
            cds_list_for_this_promoter.append(next_element_id_in_list)
          elif next_element_role_str == 'Promoter':
            break # Found the next promoter, stop collecting for the current one
          j += 1

        if cds_list_for_this_promoter:
          gene_groups_by_promoter.append((promoter_id_current, cds_list_for_this_promoter))
        i = j - 1 # Continue from the element that broke the inner loop (or end)
      i += 1

    # Process each promoter-CDS group to define a gene
    for promoter_id_val, cds_ids_list in gene_groups_by_promoter:
      gene_counter += 1
      gene_name_gro_formatted = f'"Operon_{gene_counter}"'

      gene_definition = {
        "name": gene_name_gro_formatted,
        "proteins": {f'"{cds_id}"' for cds_id in cds_ids_list},
        "promoter": {
          "function": '"TRUE"' if is_operon_constitutive else '"UNKNOWN"', # Default logic
          "transcription_factors": set()
        }
      }
      
      effective_regulators_for_promoter = {} 

      # Determine TFs and their logic for the current promoter
      for reg_interaction_item in interactions_list:
        reg_interaction_type = reg_interaction_item.get("type")
        if reg_interaction_type not in ["Control", "Stimulation", "Inhibition"]:
          continue

        reg_participants_list = reg_interaction_item.get("participants", [])
        # The component directly targeted by this regulatory interaction
        affected_comp_in_interaction = next((p_data["participant"] for p_data in reg_participants_list if p_data["role"] in ["Modified", "Stimulated", "Inhibited", "Template"]), None)
        regulator_original_id = next((p_data["participant"] for p_data in reg_participants_list if p_data["role"] in ["Modifier", "Stimulator", "Inhibitor"]), None)

        if not regulator_original_id or not affected_comp_in_interaction:
          continue

        # Check if this interaction affects the current promoter (directly or via a controlled element)
        if affected_comp_in_interaction == promoter_id_val or \
           is_promoter_under_control_of_affected_component(affected_comp_in_interaction, promoter_id_val, hierarchy_map):

          original_interaction_logic = "YES"
          if reg_interaction_type == "Inhibition":
            original_interaction_logic = "NOT"

          final_tf_gro_id = None
          final_tf_logic = original_interaction_logic

          # Case A: The original regulator is a protein TF that is itself modulated by a chemical
          if regulator_original_id in chemical_effects_on_protein_regulators:
            qs_protein_modulator_id, logic_chem_on_tf = chemical_effects_on_protein_regulators[regulator_original_id]
            final_tf_gro_id = qs_protein_modulator_id # The chemical's QS_Protein acts as the TF
            # Adjust logic: e.g., if chemical inhibits an activator (YES), final logic is NOT
            if logic_chem_on_tf == "NOT":
              final_tf_logic = "NOT" if original_interaction_logic == "YES" else "YES"

          # Case B: The original regulator is a chemical that directly affects a DNA element (promoter/operator)
          elif regulator_original_id in direct_chemical_inducer_to_qs_protein_map:
            final_tf_gro_id = direct_chemical_inducer_to_qs_protein_map[regulator_original_id]
            # Logic is already set by original_interaction_logic

          # Case C: The original regulator is a protein TF acting directly (not chemically modulated)
          else:
            # Only resolve to CDS if it's not a known chemical (already handled or not relevant here)
            if not ed_lookup.get(regulator_original_id, {}).get("type") == "Simple chemical":
              protein_tf_cds_id = get_template_for_participant(regulator_original_id, ed_definitions_list, interactions_list, hierarchy_map, components_list)
              final_tf_gro_id = protein_tf_cds_id if protein_tf_cds_id else regulator_original_id # Fallback to original
            else:
              # This case (chemical not in direct_chemical_inducer_map) implies it might not be directly regulating this promoter's DNA.
              # Or it's a chemical whose QS_Protein setup was missed. For safety, use original name as fallback.
              final_tf_gro_id = regulator_original_id

          if final_tf_gro_id:
            # Store the effective TF and its logic. Handle potential conflicts if multiple interactions define the same TF differently.
            if final_tf_gro_id not in effective_regulators_for_promoter:
              effective_regulators_for_promoter[final_tf_gro_id] = final_tf_logic
            # else: Conflict detected, current implementation keeps the first logic found. Could be refined.

      # Finalize promoter logic and TFs based on effective regulators
      effective_regulators_list = list(effective_regulators_for_promoter.items())
      current_gene_promoter_tfs = gene_definition["promoter"]["transcription_factors"]

      if not effective_regulators_list: # No regulators found
        if not is_operon_constitutive:
          gene_definition["promoter"]["function"] = '"FALSE"' # Non-constitutive and no activators means OFF
      elif len(effective_regulators_list) == 1: # Single regulator
        tf_id_val, logic_val = effective_regulators_list[0]
        gene_definition["promoter"]["function"] = f'"{logic_val}"'
        current_gene_promoter_tfs.add(f'"{tf_id_val}"')
      else: # Multiple regulators -> default to AND logic
        gene_definition["promoter"]["function"] = '"AND"'
        for tf_id_val, logic_val in effective_regulators_list:
          tf_entry = f'"-{tf_id_val}"' if logic_val == "NOT" else f'"{tf_id_val}"'
          current_gene_promoter_tfs.add(tf_entry)

      # If operon was marked constitutive but has regulators, update function from simple "TRUE"
      if is_operon_constitutive and effective_regulators_list:
        if gene_definition["promoter"]["function"] == '"TRUE"': # Default constitutive "TRUE"
          if len(effective_regulators_list) > 1:
            gene_definition["promoter"]["function"] = '"AND"'
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
        # Get parameters using the original name, as keys in user_signal_parameters
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
  # qs_action_setup_map has: norm_chem_name -> (QS_Protein_ID, config_type, original_chem_name)
  for _normalized_chem_name, (qs_protein_id, config_type, original_signal_id) in qs_action_setup_map.items():
    # Only generate s_get_QS for signals that are meant to be sensed from UI parameters
    # NCB_REACTANT types are handled by NCB aux genes, not direct s_get_QS.
    if config_type == "SENSING_CONFIG_FROM_UI": # This was the key from extract_genes
      # Use original_signal_id to look up UI-defined parameters (Symbol, Threshold)
      signal_ui_params = user_signal_parameters.get(original_signal_id, {})
      operator_symbol = signal_ui_params.get("Symbol_getQS", ">") # Default to ">"
      threshold_value = signal_ui_params.get("Threshold_getQS", "0.3") # Default threshold

      # Use tostring() for the signal name in GRO, with its normalized form
      signal_gro_id_tostring = f"tostring({normalize_signal_name(original_signal_id)})"

      sensing_actions.append(f'action({{}}, "s_get_QS", {{{signal_gro_id_tostring}, "{operator_symbol}", "{threshold_value}", "{qs_protein_id}"}});')
  return sensing_actions

def generate_gro_file(simulation_params, gene_definitions_list, qs_actions_map_data, signal_definitions_list, protein_paint_actions_map, _degradation_actions_list, biochemical_reactions_data_list, output_file_path):
  """
  Generates the content of a .gro file based on processed simulation parameters and biological constructs.

  Args:
    simulation_params (dict): Dictionary of global and component-specific parameters from the UI.
    gene_definitions_list (list): List of gene definitions (dictionaries) from extract_genes_and_qs_actions.
    qs_actions_map_data (dict): Map for QS actions from extract_genes_and_qs_actions.
    signal_definitions_list (list): List of s_signal definition strings from extract_signal_definitions.
    protein_paint_actions_map (dict): Map of protein IDs to color names for paint actions.
    _degradation_actions_list (list): List of degradation action tuples (currently not used in this function).
    biochemical_reactions_data_list (list): List of biochemical reaction data from extract_biochemical_reactions.
    output_file_path (str): Path to write the generated .gro file.
  """
  gro_content_lines = ["include gro\n"] 

  # --- Global Simulation Parameters ---
  gro_content_lines.append("// Global Simulation Parameters")
  gro_content_lines.append(f'set ("dt", {simulation_params.get("dt", 0.1)} ); // Timestep in minutes')
  gro_content_lines.append(f'set ("population_max", {simulation_params.get("population_max", 20000)} ); // Maximum cell population\n')

  # --- Signal Definitions and Setup ---
  if signal_definitions_list:
    gro_content_lines.append("// Signal Diffusion Parameters")
    gro_content_lines.append('set("signals", 1.0); // Enable signal module')
    gro_content_lines.append('set("signals_draw", 1.0); // Enable signal drawing')
    # Example grid setup, consider making this configurable via simulation_params if needed
    gro_content_lines.append('grid("continuous", "gro_original", 10, 10, 8); // Grid (type, diffusion_method, length, cell_size, neighborhood_depth)\n')
    gro_content_lines.append("// Signal Definitions (kdiff = diffusion rate, kdeg = degradation rate)")
    gro_content_lines.extend(signal_definitions_list)
    gro_content_lines.append("") 

  # --- Gene Definitions ---
  gro_content_lines.append("// Gene Definitions")
  for gene_data_item in gene_definitions_list:
    protein_list_str = ", ".join(gene_data_item["proteins"]) 
    tf_list_str = ", ".join(gene_data_item["promoter"]["transcription_factors"]) 

    gene_gro_name = gene_data_item["name"] 
    promoter_function_str = gene_data_item["promoter"]["function"] 
    gene_name_key_for_params = gene_gro_name.strip('"')

    if gene_data_item.get("is_aux_ncb_gene"):
      gro_content_lines.append(f"// Auxiliary gene for Non-Covalent Binding induced emission: {gene_name_key_for_params}")
      fixed_gene_params = gene_data_item["fixed_params"]
      act_time, act_var = fixed_gene_params["act_time"], fixed_gene_params["act_var"]
      deg_time, deg_var = fixed_gene_params["deg_time"], fixed_gene_params["deg_var"]
      to_on_noise, to_off_noise, noise_t = fixed_gene_params["toOn"], fixed_gene_params["toOff"], fixed_gene_params["noise_time"]
    else:
      # Default gene timing parameters if not found in simulation_params
      default_timings = {"act_time": 10.0, "act_var": 1.0, "deg_time": 20.0, "deg_var": 1.0,
                         "toOn": 0.0, "toOff": 0.0, "noise_time": 100.0}
      gene_timing_params = simulation_params.get("gene_parameters", {}).get(gene_name_key_for_params, default_timings)

      act_time = float(gene_timing_params.get("act_time", default_timings["act_time"]))
      act_var = float(gene_timing_params.get("act_var", default_timings["act_var"]))
      deg_time = float(gene_timing_params.get("deg_time", default_timings["deg_time"]))
      deg_var = float(gene_timing_params.get("deg_var", default_timings["deg_var"]))
      to_on_noise = float(gene_timing_params.get("toOn", default_timings["toOn"]))
      to_off_noise = float(gene_timing_params.get("toOff", default_timings["toOff"]))
      noise_t = float(gene_timing_params.get("noise_time", default_timings["noise_time"]))

      # Ensure parameters are within valid ranges
      act_time = max(0.0, act_time); act_var = max(0.0, act_var)
      deg_time = max(0.0, deg_time) 
      deg_var = max(0.0, deg_var)
      to_on_noise = max(0.0, min(1.0, to_on_noise)) # Probability [0,1]
      to_off_noise = max(0.0, min(1.0, to_off_noise)) # Probability [0,1]
      noise_t = max(0.0, noise_t)

    gro_content_lines.extend([
      f'genes([',
      f'  name := {gene_gro_name},',
      f'  proteins := {{{protein_list_str}}},',
      f'  promoter := [function := {promoter_function_str},',
      f'    transcription_factors := {{{tf_list_str}}},',
      f'    noise := [toOn := {to_on_noise}, toOff := {to_off_noise}, noise_time := {noise_t}]',
      f'  ],',
      f'  prot_act_times := [times := {{{act_time}}}, variabilities := {{{act_var}}}],',
      f'  prot_deg_times := [times := {{{deg_time}}}, variabilities := {{{deg_var}}}]',
      f']);\n'
    ])

  # --- Plasmid Definitions ---
  plasmid_definitions_map = simulation_params.get("plasmid_configuration", {}).get("defined_plasmids", {})
  if plasmid_definitions_map:
    gro_content_lines.append("// Plasmid Definitions")
    plasmid_entries_list = []
    for plasmid_id_str, genes_on_plasmid_list in plasmid_definitions_map.items():
      formatted_genes_on_plasmid_str = ", ".join(genes_on_plasmid_list)
      plasmid_entries_list.append(f'{plasmid_id_str} := {{{formatted_genes_on_plasmid_str}}}')

    if plasmid_entries_list:
        gro_content_lines.append(f'plasmids_genes([ {", ".join(plasmid_entries_list)} ]);\n')


  # --- Quorum Sensing Signal Detection Actions (s_get_QS) ---
  if qs_actions_map_data:
    qs_sensing_action_strings = generate_qs_signal_sensing_actions(qs_actions_map_data, simulation_params.get("signal_parameters", {}))
    if qs_sensing_action_strings:
      gro_content_lines.append("// Quorum Sensing Signal Detection Actions")
      gro_content_lines.extend(qs_sensing_action_strings)
      gro_content_lines.append("")

  # --- Non-Covalent Binding Product Emission Actions ---
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

  # --- Biochemical Signal Conversion Actions (S1 + E -> S2) ---
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

          if conversion_rate > 0.0: # Only add actions if rate is positive
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

  # --- Cell Painting Actions ---
  all_gene_protein_ids = {p_id.strip('"') for gene in gene_definitions_list for p_id in gene["proteins"]}
  valid_protein_paint_config = {p_id: color for p_id, color in protein_paint_actions_map.items() if p_id in all_gene_protein_ids}

  if valid_protein_paint_config:
    paint_action_strings = generate_protein_paint_actions(valid_protein_paint_config)
    if paint_action_strings:
      gro_content_lines.append("// Cell Painting Actions (based on protein presence)")
      gro_content_lines.extend(paint_action_strings)
      gro_content_lines.append("")

  # --- Bacterial Conjugation Actions ---
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

  # --- Program Definition ---
  gro_content_lines.append(f'program p() := {{\n  skip();\n}};\n')
  # --- Main Program: Initial Cell and Signal Setup ---
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

  # --- Write to .gro File ---
  try:
    with open(output_file_path, 'w', encoding='utf-8') as gro_file_handle:
      gro_file_handle.write("\n".join(gro_content_lines))
    #print(f".gro file generated successfully: {output_file_path}")
  except IOError as e:
    #print(f"Error writing .gro file: {e}")
