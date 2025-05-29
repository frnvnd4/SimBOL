# SimBOL - GRO Plugin: User Manual

Welcome to the User Manual for the SimBOL GRO Plugin. This guide provides detailed instructions on how to use SimBOL to convert your SBOL3 genetic circuit designs into simulation-ready `.gro` files, focusing on the requirements for your SBOL3 input and how to effectively use the parameter configuration interface.

## Table of Contents

1.  [Introduction to the SimBOL-GRO Plugin](#1-introduction-to-the-simbol-gro-plugin)
2.  [Core Workflow](#2-core-workflow)
3.  [SBOL3 Input File Requirements (Crucial for GRO Conversion)](#3-sbol3-input-file-requirements-crucial-for-gro-conversion)
    * [3.1 Component Typing (`type` and `role`)](#31-component-typing-type-and-role)
    * [3.2 Operon Definitions (Transcriptional Units)](#32-operon-definitions-transcriptional-units)
    * [3.3 Comprehensive Interaction Specification](#33-comprehensive-interaction-specification)
    * [3.4 Identifiers (`displayId`, `name`)](#34-identifiers-displayid-name)
    * [3.5 Hierarchy (Recommended)](#35-hierarchy-recommended)
4.  [Using the SimBOL-GRO Plugin via Google Colab](#4-using-the-simbol-gro-plugin-via-google-colab)
    * [4.1 Step 1: Environment Setup](#41-step-1-environment-setup)
    * [4.2 Step 2: Upload SBOL3 File](#42-step-2-upload-sbol3-file)
    * [4.3 Step 3: Automatic SBOL3 to Optimized JSON Conversion](#43-step-3-automatic-sbol3-to-optimized-json-conversion)
    * [4.4 Step 4: GRO Simulation Parameter Configuration (UI Detailed)](#44-step-4-gro-simulation-parameter-configuration-ui-detailed)
        * [4.4.1 Global Simulation Parameters](#441-global-simulation-parameters)
        * [4.4.2 Gene Parameters](#442-gene-parameters)
        * [4.4.3 Plasmid Definition](#443-plasmid-definition)
        * [4.4.4 Initial Cell Populations](#444-initial-cell-populations)
        * [4.4.5 Signal Parameters](#445-signal-parameters)
        * [4.4.6 Quorum Sensing (QS) Parameters](#446-quorum-sensing-qs-parameters)
        * [4.4.7 Non-Covalent Binding (NCB) Induced Emission](#447-non-covalent-binding-ncb-induced-emission)
        * [4.4.8 Biochemical Conversion (S1+E->S2)](#448-biochemical-conversion-s1es2)
        * [4.4.9 Bacterial Conjugation](#449-bacterial-conjugation)
        * [4.4.10 Saving Parameters](#4410-saving-parameters)
    * [4.5 Step 5: Generate and Download `.gro` File](#45-step-5-generate-and-download-gro-file)
5.  [Understanding the Generated `.gro` File](#5-understanding-the-generated-gro-file)
6.  [Tips, Warnings, and Best Practices](#6-tips-warnings-and-best-practices)
7.  [Troubleshooting Common Issues](#7-troubleshooting-common-issues)

---

## 1. Introduction to the SimBOL-GRO Plugin

SimBOL is a framework that translates SBOL3 designs into formats for various simulators. This manual specifically covers the **SimBOL-GRO Plugin**, which enables the generation of `.gro` files for the GRO simulator.

The plugin leverages an intermediate, optimized JSON representation derived from your SBOL3 file. This JSON simplifies the complex SBOL3 data, making it easier to map to GRO's simulation constructs and allowing you to fine-tune simulation-specific parameters through an interactive user interface.

## 2. Core Workflow

The SimBOL-GRO plugin follows these main steps (typically within a Google Colab environment):

1.  **SBOL3 Input:** You provide your genetic circuit design as an SBOL3 file.
2.  **JSON Conversion:** SimBOL automatically processes the SBOL3 file into an optimized intermediate JSON format. This JSON provides a summarized, structured view of your circuit and is also available for download.
3.  **Parameter Configuration:** An interactive UI appears, populated with data from the JSON. Here, you specify or adjust parameters required for the GRO simulation.
4.  **`.gro` File Generation:** The plugin uses the JSON data and your configured parameters to create a `.gro` file.
5.  **Download:** You can download the generated `.gro` file for use in the GRO simulator.

## 3. SBOL3 Input File Requirements (Crucial for GRO Conversion)

For the SimBOL-GRO plugin to correctly interpret your design and generate a functional `.gro` file, your input SBOL3 file must adhere to specific guidelines. A well-specified SBOL3 design is key for a successful conversion.

### 3.1 Component Typing (`type` and `role`)

* **Component `type`:** Every `Component` in your SBOL3 design (e.g., DNA, RNA, Protein, Simple Chemical) MUST have its `type` property clearly defined using appropriate ontology terms (e.g., from the Systems Biology Ontology - SBO).
* **Component `role`:** For DNA components that form transcriptional units, it is **essential** that their `role` properties are correctly specified using Sequence Ontology (SO) terms.
    * The GRO plugin **requires at least a Promoter (`SO:0000167`) and one or more Coding Sequences (CDS) (`SO:0000316`)** to define a gene (operon) in the `.gro` file.
    * Other roles like RBS (`SO:0000139`) and Terminator (`SO:0000141`) are good practice for a complete SBOL3 design but are not explicitly parsed into distinct GRO elements by the current plugin beyond their inclusion in the overall gene structure.

### 3.2 Operon Definitions (Transcriptional Units)

The GRO simulator defines genes as transcriptional units. SimBOL attempts to map SBOL3 constructs to this concept.

* **Preferred Method (Clear & Robust):**
    * Define each operon as a distinct SBOL3 `Component` (typically of `type` DNA).
    * Within this "operon" `Component`, instantiate all its constituent genetic parts (e.g., promoter, RBS, CDSs, terminator) as `SubComponent` objects linked via the `hasFeature` property.
    * The order and relationship of these `SubComponent` instances within the parent operon `Component` define the transcriptional unit. This explicit hierarchical structure is the most reliable way for SimBOL to interpret your operons.

* **Automatic Grouping (Fallback Behavior):**
    * If operons are not explicitly structured as described above, the SimBOL-GRO plugin will attempt to infer transcriptional units.
    * It does this by flattening the hierarchy of all DNA components found in the design and then grouping a `Promoter` with all subsequent `CDS`s until another `Promoter` is encountered in sequence. Each such grouping is then treated as a separate "gene" (operon) in the generated `.gro` file.
    * **Warning:** This automatic grouping is a heuristic. While it attempts to build sensible units, it may not always match complex biological realities or designer intent if the SBOL3 hierarchy is ambiguous. **Explicit hierarchical definition of operons in SBOL3 is strongly recommended for accurate conversion.**

### 3.3 Comprehensive Interaction Specification

All molecular interactions relevant to the circuit's behavior MUST be present and correctly defined in the SBOL3 file, using appropriate SBO terms for interaction types and participant roles. This information is crucial for generating corresponding `action` statements in the `.gro` file.

* **Genetic Production (`SBO:0000589`):**
    * **Absolutely essential.** These interactions define which `Component` (acting as a template, typically a CDS) produces which protein `Component` (product).
    * The GRO plugin uses this to:
        * Map proteins back to the CDS that encodes them (this CDS becomes the "protein name" in many GRO contexts).
        * Identify transcription factors and other proteins involved in regulation.
* **Regulatory Interactions (e.g., `Control` SBO:0000168, `Inhibition` SBO:0000169, `Stimulation` SBO:0000170):**
    * These define how transcription factors (proteins, complexes) or chemical signals affect promoters (or associated DNA elements like operators).
    * Clearly specify the `participant` with the `role` of `Modifier`, `Stimulator`, or `Inhibitor`, and the `participant` with the `role` of `Modified`, `Stimulated`, or `Inhibited` (which is typically the promoter `Feature` or an operator `Feature` linked to that promoter).
* **Non-Covalent Binding (`SBO:0000177`):**
    * If the binding of two or more molecules (e.g., Protein + Simple Chemical) forms a complex that subsequently acts as a transcription factor or results in the emission of a new signal, this interaction is vital.
    * The GRO plugin uses this to:
        * Identify reactants that form a complex.
        * If the product is a `Simple chemical`, it generates an auxiliary gene system in GRO to model the conditional emission of this chemical product (see Section 4.4.7).
* **Biochemical Reaction (`SBO:0000176`):**
    * Used for enzymatic conversions, typically S1 (substrate) + E (enzyme) -> S2 (product).
    * The GRO plugin specifically looks for reactions where S1 and S2 are `Simple chemical` type Components, and E is a `Protein` type Component acting as a `Modifier`.
    * This allows the generation of `s_absorb_signal` (for S1) and `s_emit_signal` (for S2) actions in GRO, triggered by the presence of the enzyme.
* **Degradation (`SBO:0000179`):**
    * While SBOL3 defines a "Degradation" interaction type, **GRO does not have a direct `action` for one entity to actively degrade another** (e.g., `action({"ProteinX"}, "degrade", {"ProteinY"})`).
    * In GRO, degradation is primarily an *intrinsic property*:
        * Proteins: Degradation is set by `prot_deg_times` in the `genes` definition.
        * Signals: Degradation is set by `kdeg` in the `s_signal` definition.
    * Your `extract_degradation_actions` function helps identify these SBOL interactions. However, the SimBOL-GRO plugin currently **does not translate these into specific GRO `action` statements**. This information could be used to manually inform the setting of intrinsic `prot_deg_times` or `kdeg` values, or to model complex indirect effects if active, conditional degradation is critical. (See Section 6 for more).

### 3.4 Identifiers (`displayId`, `name`)

* Use consistent and meaningful `displayId` properties for your SBOL3 `Component` and `Feature` objects. These are often used as the basis for names in the generated `.gro` file.
* Be aware that `normalize_signal_name` is used for chemical signals to make them valid GRO identifiers (typically removing spaces and hyphens).

### 3.5 Hierarchy (Recommended)

* While SimBOL's SBOL-to-JSON phase can flatten hierarchies, structuring your SBOL3 designs hierarchically (e.g., an "Operon" `Component` containing "Promoter" and "CDS" `SubComponent`s) greatly aids in unambiguous interpretation and is considered good design practice.

## 4. Using the SimBOL-GRO Plugin via Google Colab

The primary way to use the SimBOL-GRO plugin is through the provided Google Colab notebook (e.g., `SimBOL.ipynb`), typically found in the `notebooks/` directory of the repository.

### 4.1 Step 1: Environment Setup

* **Action:** Run the first code cell in the notebook.
* **Purpose:** This cell installs all necessary Python libraries (`sbol3`, `ipywidgets`, etc.) and clones the latest version of the SimBOL code from GitHub into your Colab environment. It also adds the SimBOL project directory to Python's path so its modules can be imported.
* **Output:** You'll see installation and cloning logs, followed by a confirmation message. If errors occur here, it's usually related to network connectivity or GitHub access.

### 4.2 Step 2: Upload SBOL3 File

* **Action:** Run the second code cell.
* **Purpose:** This cell will prompt you to upload your SBOL3 design file (usually with a `.sbol` or `.xml` extension).
* **UI:** A file upload button will appear. Click it to select your file.
* **Output:** Confirmation that your file has been uploaded.

### 4.3 Step 3: Automatic SBOL3 to Optimized JSON Conversion

* **Action:** Run the third code cell.
* **Purpose:** SimBOL processes the uploaded SBOL3 file and converts it into an intermediate, optimized JSON format. This JSON summarizes your circuit's components, hierarchy, and interactions in a structured way that's easier for the subsequent steps.
* **Output:**
    * A message indicating successful processing.
    * **The generated intermediate JSON file will be automatically downloaded by your browser.** This allows you to inspect how SimBOL has interpreted your design before proceeding.

### 4.4 Step 4: GRO Simulation Parameter Configuration (UI Detailed)

* **Action:** Run the fourth code cell.
* **Purpose:** This displays an interactive user interface (UI) built with `ipywidgets`. The UI is pre-populated with information extracted from the intermediate JSON. Here, you will review and configure parameters specifically for the GRO simulation.
* **UI Sections:**

    * #### 4.4.1 Global Simulation Parameters
        * **Timestep (dt, minutes):** The duration of each simulation timestep.
        * **Population max:** The maximum number of cells allowed in the simulation.
        * **Cell Doubling Time (minutes):** Used to calculate the default growth rate.

    * #### 4.4.2 Gene Parameters
        * For each "gene" (transcriptional unit identified by SimBOL, e.g., "Operon\_1", "Operon\_2"):
            * **Produces:** Lists the protein(s) (CDS IDs from SBOL3) produced by this gene.
            * **Timing (Activation/Degradation, minutes):**
                * `Act. Time`: Average time for protein activation/expression.
                * `Act. Var.`: Variability (standard deviation) for activation time.
                * `Deg. Time`: Average time for protein degradation.
                * `Deg. Var.`: Variability for degradation time.
            * **Noise Parameters:**
                * `Noise Act. Time`: Timespan over which noise probabilities apply.
                * `P(noise ON)` (`toOn`): Probability of the promoter being constitutively ON due to noise.
                * `P(noise OFF)` (`toOff`): Probability of the promoter being constitutively OFF due to noise.

    * #### 4.4.3 Plasmid Definition
        * **Add Plasmid:** Allows you to define new plasmids by name.
        * **Assign Genes:** For each defined plasmid, select (using checkboxes) which of the above "genes" (Operon\_1, Operon\_2, etc.) it contains.
        * **Warning:** As discovered, ensure that plasmid definitions are unique if they have different names (i.e., two differently named plasmids should ideally not contain the exact same set of operons if this causes issues in GRO).

    * #### 4.4.4 Initial Cell Populations
        * **Add Population:** Define one or more groups of E. coli cells.
        * For each group:
            * `Number of Cells`: How many cells in this group.
            * `Center X`, `Center Y`: Coordinates for the center of the initial cell cluster.
            * `Radius`: Radius of the circular area for initial cell placement.
            * `Select Plasmids`: Using checkboxes, assign which of the defined plasmids this cell population group will initially contain.

    * #### 4.4.5 Signal Parameters
        * For each `Simple chemical` identified from your SBOL3 EDs:
            * `Diffusion Rate (kdiff)`: The diffusion coefficient for the signal.
            * `Degradation Rate (kdeg)`: The degradation coefficient for the signal.
            * **Initial Signal Placement:** If the signal is not generated internally by NCB or biochemical conversion:
                * `Add Additional Initial Signal Point`: Add one or more locations for initial signal placement.
                * For each point: `X`, `Y` coordinates, `Conc` (concentration).
                * `Constant Emission`: Checkbox to indicate if this concentration should be emitted from this point at every timestep (`s_set_signal` within a `true: {}` block in `main()`). If unchecked, it's a one-time placement (`s_set_signal` directly in `main()`).

    * #### 4.4.6 Quorum Sensing (QS) Parameters
        * If a signal was identified as an effector in SBOL3 regulatory interactions, you can configure its QS behavior:
            * `QS Operator`: The comparison operator (`>`, `<`) for the `s_get_QS` action.
            * `QS Threshold`: The concentration threshold for the `s_get_QS` action.
            * *(The QS_Protein name is generated automatically by SimBOL based on the signal name).*

    * #### 4.4.7 Non-Covalent Binding (NCB) Induced Emission
        * If an NCB interaction produces a `Simple chemical` (as identified by `extract_ncb_production_genes_and_actions`):
            * The UI will list the **Emitted Signal** and the auxiliary protein (`P_aux_NCBx`) that triggers its emission.
            * `Emission Conc`: Concentration of the signal to be emitted when `P_aux_NCBx` is present.
            * `Emission Type`: How the signal is emitted (e.g., "exact", "area", "average") for the `s_emit_signal` action.

    * #### 4.4.8 Biochemical Conversion (S1+E->S2)
        * If an enzymatic conversion of signal S1 to S2 by enzyme E was identified:
            * The UI will list the **Reaction** (S1 → S2, Catalyzed by: Enzyme).
            * `S1 to S2 Conversion Speed`: (Formerly "Conversion Rate (S1 Abs & S2 Em, units/dt)") The amount of S1 absorbed and S2 emitted per cell per timestep when the enzyme is present.
            * `S1 Abs. Type`: Absorption type for S1 (e.g., "exact", "area", "average") for `s_absorb_signal`.
            * `S2 Em. Type`: Emission type for S2 (e.g., "exact", "area", "average") for `s_emit_signal`.

    * #### 4.4.9 Bacterial Conjugation
        * `Enable Bacterial Conjugation`: Checkbox to enable/disable conjugation features.
        * If enabled:
            * `Select Plasmids for Conjugation`: Checkboxes to select which defined plasmids can be conjugated.
            * For each selected plasmid:
                * `Rate (events/doubling time) for [plasmid_name]`: The average number of conjugation events for this plasmid per cell doubling time.

    * #### 4.4.10 Saving Parameters
        * **Crucial Step:** After configuring all parameters, click the **"Save Parameters"** button at the bottom of the form. This action stores your settings for the next step. You should see a confirmation message "Parameters saved successfully."

### 4.5 Step 5: Generate and Download `.gro` File

* **Action:** Run the final code cell in the notebook.
* **Purpose:** This cell takes the parameters you saved from the UI and the processed JSON data to:
    1.  Prepare final data structures (e.g., signal definitions, actions lists).
    2.  Generate the complete content of the `.gro` file.
    3.  Initiate the download of the `.gro` file to your computer.
* **Output:**
    * A message confirming successful `.gro` file generation.
    * Your browser will start downloading the `.gro` file.

## 5. Understanding the Generated `.gro` File

The generated `.gro` file will contain several sections based on your SBOL3 input and UI configurations:

* `include gro`
* **Global Parameters:** `set("dt", ...)`, `set("population_max", ...)`.
* **Signal Definitions:** `SignalName := s_signal([kdiff := ..., kdeg := ...]);`
* **Gene Definitions:** `genes([... proteins:={}, promoter:=[...], ...]);` blocks for each transcriptional unit.
* **Plasmid Definitions:** `plasmids_genes([ p1 := {"Operon_1"}, ... ]);`
* **Actions:**
    * `action({}, "s_get_QS", ...);` for quorum sensing.
    * `action({"Protein"}, "s_emit_signal", ...);` for NCB-induced emissions.
    * `action({"Enzyme"}, "s_absorb_signal", ...);` and `action({"Enzyme"}, "s_emit_signal", ...);` for biochemical conversions.
    * `action({"Protein"}, "paint", ...);` for cell painting (if configured).
    * `action({}, "conjugate", ...);` for plasmid conjugation.
    * `action({}, "set_growth_rate", ...);` to set a global growth rate.
* **Program Definitions:**
    * `program p() := { skip(); };` (The default program for E. coli, as discussed).
    * `program main() := { ... };` (Contains initial signal placements, constant emissions, and initial cell population setup using `c_ecolis`).

## 6. Tips, Warnings, and Best Practices

* **SBOL3 File Quality is Key:** The more accurately and completely your SBOL3 file describes your circuit (especially component types, roles, and interactions), the better SimBOL can interpret it and generate a meaningful GRO file.
* **Review the Intermediate JSON:** Take advantage of the downloaded intermediate JSON. It shows you how SimBOL "sees" your circuit before you even get to the GRO parameters. This can help you spot issues in your SBOL3 file early.
* **Iterative Refinement:** Start with simpler circuits if you are new to SBOL3 or GRO. Generate a `.gro` file, inspect it, run a short simulation, and then iteratively add complexity.
* **Parameter Sensitivity:** Simulation outcomes in GRO can be sensitive to parameters like `dt`, degradation/diffusion rates, and interaction thresholds. Experiment with these values.
* **Degradation Modeling:** Remember that direct protein-mediated degradation of other proteins is not a standard GRO `action`. Model this through intrinsic degradation rates (`prot_deg_times` in gene definitions, or `kdeg` for signals) or through more complex indirect regulatory pathways.
* **Auxiliary Genes for NCB:** For Non-Covalent Binding interactions that produce a chemical signal, SimBOL automatically creates an "auxiliary gene" (e.g., `Operon_NCB_Emit_SignalX_1`) that produces an auxiliary protein (e.g., `P_aux_NCB1`). This auxiliary protein then triggers the `s_emit_signal` action for the chemical product. You will see these auxiliary genes in the "Gene Parameters" section of the UI, but their kinetic parameters are fixed by SimBOL to be very fast/efficient as they are meant to be direct relays of the binding event.

## 7. Troubleshooting Common Issues

* **`ModuleNotFoundError` in Colab:**
    * Ensure the "Environment Setup" cell (Cell 1) has run successfully and cloned the SimBOL repository.
    * Confirm `sys.path.append('/content/SimBOL')` was executed.
* **Errors during "Convert SBOL3 to JSON" (Cell 3):**
    * Usually indicates an issue with the input SBOL3 file. Validate your SBOL3 file against the SBOL3 specification.
    * Ensure all required information (types, roles, interactions as detailed in Section 3) is present.
* **Parameter UI Not Appearing or Empty (Cell 4):**
    * The SBOL3 file might not have been processed correctly into JSON in the previous step. Check for errors in Cell 3.
    * The SBOL3 file might lack recognizable components or interactions that the UI populates.
* **`.gro` File Not Generating or Incorrect (Cell 5):**
    * **Did you click "Save Parameters" in the UI?** This is a common oversight.
    * Check the Colab cell output for any error messages from the `generate_gro_file` function.
    * Review the parameters you set in the UI for any obvious inconsistencies.
