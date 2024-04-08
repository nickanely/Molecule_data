# Molecular Properties Processor
This project is a Python-based tool that processes a CSV file containing information about molecules, computes various molecular properties using the RDKit library, and writes the results to an Excel file.

### Key Features
1. **CSV File Processing:** The program reads a CSV file containing information about molecules, including the SMILES (Simplified Molecular-Input Line-Entry System) strings, pIC50 values, number of atoms, and logP values.


2. **Molecular Property Computation**: The program uses the RDKit library to compute the following molecular properties for each molecule:
   * Molecular weight
   * Topological Polar Surface Area (TPSA)
   * LogP (octanol-water partition coefficient)
   * Number of hydrogen bond acceptors
   * Number of hydrogen bond donors
   * Ring count
   * Lipinski rule of five (a set of criteria for drug-likeness)
   

3. Multiprocessing: To improve performance, the program utilizes multiprocessing to split the data into smaller chunks and process them concurrently. This helps to take advantage of the available CPU cores and speed up the computation.


4. Excel Output: The computed molecular properties are then written to an Excel file named "library.xlsx".


### Key Implementation Details
1. **Class Structure:** The program is organized as a `MolecularPropertiesProcessor` class, which encapsulates the data processing and computation logic.


2. **Column Finder:** The `_column_finder` method is used to identify the column names in the input CSV file, ensuring that the program can work with different file formats.


3. **Chunk-based Processing:** The `_compute_molecule_properties` method splits the data into smaller chunks and processes them concurrently using the multiprocessing.Pool class. This helps to distribute the workload across available CPU cores.


4. **Molecule Property Computation:** The `_compute_molecule_properties_chunk` method is responsible for computing the molecular properties for each chunk of data using the RDKit library.


5. **Error Handling:** The program includes a custom `MoleculeProcessingException` class to handle any exceptions that may occur during the data processing.

   


### Usage
To use this program, follow these steps:

1. Ensure that you have Python 3.6 or higher installed, as well as the required dependencies (pandas, numpy, rdkit, and multiprocessing).


2. Place your CSV file containing the molecular data in the same directory as the Python script.


3. Update the `input_file_path` variable in the `MolecularPropertiesProcessor` class to point to your CSV file.


4. Run the script, and the resulting Excel file "library.xlsx" will be generated in the same directory.