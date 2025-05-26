import sbol3

NS = {
  "sbol": "http://sbols.org/v3#",
  "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
}

inverse_interaction_mapping = {
  sbol3.SBO_INHIBITION: 'Inhibition',
  sbol3.SBO_STIMULATION: 'Stimulation',
  sbol3.SBO_BIOCHEMICAL_REACTION: 'Biochemical Reaction',
  sbol3.SBO_NON_COVALENT_BINDING: 'Non-Covalent Binding',
  sbol3.SBO_DEGRADATION: 'Degradation',
  sbol3.SBO_GENETIC_PRODUCTION: 'Genetic Production',
  sbol3.SBO_CONTROL: 'Control'
}

inverse_ed_mapping = {
  sbol3.SBO_SIMPLE_CHEMICAL: 'Simple chemical',
  sbol3.SBO_PROTEIN: 'Protein',
  # sbol3.SBO_PROTEIN: 'Restriction enzyme'
}

inverse_role_mapping = {
  sbol3.SBO_INHIBITOR: 'Inhibitor',
  sbol3.SBO_INHIBITED: 'Inhibited',
  sbol3.SBO_STIMULATOR: 'Stimulator',
  sbol3.SBO_STIMULATED: 'Stimulated',
  sbol3.SBO_REACTANT: 'Reactant',
  sbol3.SBO_PRODUCT: 'Product',
  sbol3.SBO_MODIFIER: 'Modifier',
  sbol3.SBO_MODIFIED: 'Modified',
  sbol3.SBO_TEMPLATE: 'Template',
  sbol3.SBO_PROMOTER: 'Promoter'
}

inverse_component_map = {
  'https://identifiers.org/SO:0000032': 'Aptamer-DNA',
  'https://identifiers.org/SO:0000033': 'Aptamer-RNA',
  'https://identifiers.org/SO:0002223': 'Inert-DNA-Spacer',
  'https://identifiers.org/SO:0000655': 'Ncrna',
  'https://identifiers.org/SO:0000553': 'PolyA-Site',
  sbol3.SO_PROMOTER: 'Promoter',
  sbol3.SO_RBS: 'RBS',
  sbol3.SO_CDS: 'CDS',
  sbol3.SO_TERMINATOR: 'Terminator',
  sbol3.SO_OPERATOR: 'Operator',
  sbol3.SO_ENGINEERED_REGION: 'Engineered-Region',
  'https://identifiers.org/SO:0005850': 'Primer-Binding-Site',
  'https://identifiers.org/SO:0000724': 'Origin-Of-Transfer',
  'https://identifiers.org/SO:0000296': 'Origin-Of-Replication',
  sbol3.SO_MRNA: 'mRNA',
  sbol3.SBO_DNA: 'DNA',
  sbol3.SBO_RNA: 'RNA',
  sbol3.SBO_FUNCTIONAL_ENTITY: 'Functional Entity'
}
