# SimBOL 🧬➡️⚙️: SBOL3 to GRO Simulation File Converter

SimBOL is a tool designed to convert genetic circuit designs specified in the SBOL3 (Synthetic Biology Open Language v3.1.0) format [cite: 1, 5] into `.gro` simulation files, ready for use with the GRO simulator[cite: 975]. This process streamlines the transition from a biological design's structural and functional specification to a simulatable model by utilizing an intermediate JSON representation for simplified data handling.

## 🌟 Overview

Synthetic biology applies engineering principles to the design of biological systems. The SBOL standard was developed to facilitate the exchange of biological design information, describing sequences, molecular interactions, and expected system behavior. GRO is a complementary language and simulator for modeling the dynamics of these circuits.

SimBOL connects these two worlds:

1.  **SBOL3 Parsing 📄:** Accepts an SBOL3 file as input.
2.  **Intermediate JSON Conversion 🔀:** Transforms complex SBOL3 data into a summarized, more accessible JSON format.
3.  **User-Friendly Parameter Input 💻:** Through an included Google Colab notebook, an interactive UI (using `ipywidgets`) allows users to review extracted data and fine-tune parameters for the GRO simulation.
4.  **`.gro` File Generation 🚀:** Generates a GRO simulation file based on the processed JSON and user-specified parameters.

## ✨ Key Features

* **SBOL3 Conversion:** Interprets SBOL3 files detailing genetic components (genes, promoters, CDS, etc.) [cite: 259, 286] and their interactions (regulation, genetic production, non-covalent binding, biochemical reactions)[cite: 446].
* **Intermediate JSON Representation:** Simplifies SBOL data structures for efficient processing.
* **Interactive User Interface (in Colab):** Allows easy configuration of crucial GRO simulation parameters, including:
    * Global settings (`dt`, `population_max`)[cite: 983].
    * Operon/gene dynamics (activation/degradation times, noise)[cite: 1014, 1016].
    * Plasmid construction and gene assignments[cite: 1042].
    * Signal properties (diffusion, degradation) [cite: 998] and initial/constant emission setups.
    * GRO actions such as Quorum Sensing (`s_get_QS`)[cite: 1105], cell painting (`paint`)[cite: 1063], conjugation (`conjugate`)[cite: 1070], and others.
* **`.gro` File Output:** Produces a complete specification file ready for the GRO simulator.

## 📁 Project Structure

The project is organized for clarity and ease of use:

* `sbol3_to_json_converter/`: Modules for parsing SBOL3 files and converting them to our intermediate JSON format.
    * `sbol_processor.py`: Core SBOL processing logic.
    * `upload_sbol_file.py`: Utility for file uploads in Colab.
    * `rdf_parser.py`, `clean_json.py`: Helper modules for parsing and data cleaning.
* `json_to_gro_generator/`: Modules that take the intermediate JSON and user parameters to generate the final `.gro` file.
    * `ui.py`: Defines the interactive user interface.
    * `builder.py`: Main logic for constructing the `.gro` file.
    * `params.py`: Prepares data and parameters for `.gro` generation.
* `notebooks/` (or `colab_notebooks/`): **Contains the primary Google Colab notebook(s)  नोटबुक for running the entire conversion workflow.**
* `example_data/`: Sample SBOL3 input files, intermediate JSONs, and generated `.gro` files.
* `README.md`: This file.
* `requirements.txt`: Python dependencies needed for the project.

## 🚀 Getting Started with SimBOL (via Google Colab)

The easiest way to use SimBOL is through the provided Google Colab notebook located in the `notebooks/` (or `colab_notebooks/`) directory of this repository.

1.  **Open in Colab ☁️:**
    * Navigate to the `notebooks/` directory in this GitHub repository.
    * Click on the main `.ipynb` notebook file.
    * You should see an "Open in Colab" badge or link at the top of the notebook preview, or you can manually open it by providing the GitHub URL to Google Colab.

2.  **Run the Notebook ▶️:**
    * Once the notebook is open in Colab, simply follow the instructions and run the cells in order.
    * The notebook will guide you through:
        * Cloning this repository into the Colab environment.
        * Installing all necessary dependencies.
        * Uploading your SBOL3 file (or using an example).
        * Interacting with the UI to set parameters for your GRO simulation.
        * Generating and downloading your `.gro` file.

    *No local installation is required if you use the Colab notebook.* The notebook handles the setup within the Colab environment.

## 🔄 Workflow within the Colab Notebook

1.  **Setup:** The initial cells clone the SimBOL repository and install required Python packages.
2.  **Input SBOL3 File:** You'll be prompted to upload your SBOL3 file.
3.  **(Automatic) JSON Conversion:** The tool processes your SBOL3 file into an internal, summarized JSON format.
4.  **Parameter Configuration UI:** An interactive interface will appear, pre-filled with information from your design. Here, you can adjust simulation parameters for GRO.
5.  **`.gro` File Generation:** After confirming parameters, the `.gro` file is generated.
6.  **Download 💾:** You'll be able to download the generated `.gro` file directly from Colab.

## 🌱 Future Enhancements


## 🤝 Contributing



## 📜 License

This project is licensed under the MIT License - see the `LICENSE.md` file for details.