# 3D Protein Structure Comparison Tool

This project provides a comprehensive tool for quantitative and qualitative comparison of three-dimensional protein structures using PDB files, Kabsch algorithm for structural alignment, and RMSD calculation for similarity measurement.

---

## 🚀 How to Run the Project

### 1. Prerequisites
- Python 3.7 or higher
- Required Python libraries (install via requirements.txt)

### 2. Installation
1. Clone or download the project
2. Install required packages:
```bash
pip install -r requirements.txt
```

### 3. Prepare PDB Files
1. Place your PDB files in the `Samples/` directory
2. Update file paths in `pdb_comparator.py`:
```python
pdb_file_1 = r"Samples/1A00.pdb"
pdb_file_2 = r"Samples/1A01.pdb"
```

### 4. Run the Program
```bash
python PDB_Comparator.py
```

---

## 🎯 Project Goal
### This project aims to:
- Compare 3D structures of two proteins quantitatively
- Perform optimal structural alignment using Kabsch algorithm
- Calculate RMSD (Root Mean Square Deviation) as similarity metric
- Generate comprehensive visualizations for structural analysis
- Identify conserved and variable regions in protein structures

---

## 📂 Project Structure

| Folder/File               | Description                                  |
|---------------------------|----------------------------------------------|
| **📄 PDB_Comparator.py**  | Main program code                           |
| **📄 requirements.txt**   | Required Python libraries                   |
| **📄 README.md**          | Project documentation                       |
| **📄 Explanation.md**     | Detailed scientific explanation             |
|                           |                                              |
| **📁 Samples/**          | Sample PDB files                            |
| - `5TIM.pdb`              | Sample protein structure 1                  |
| - `1TIM.pdb`              | Sample protein structure 2                  |
| - `1A01.pdb`              | Sample protein structure 3                  |
| - `1A00.pdb`              | Sample protein structure 4                  |
|                           |                                              |
| **📁 Images/**            | Program output examples                     |
| - `Compare_Similar.jpg`   | Comparison of similar structures            |
| - `Compare_Different.jpg` | Comparison of different structures          |
| - `Console_Results.jpg`   | Console output example                      |

---

## 🔬 Scientific Background

### Protein Data Bank (PDB) Files
- Standard format for 3D structural data of biological macromolecules
- Contains atomic coordinates, sequence information, and structural data

### Alpha Carbon (Cα) Atoms
- Backbone atoms that define protein fold and structure
- Each amino acid has one Cα atom
- Using only Cα atoms provides 20x faster computation while preserving 95% of structural information

### Kabsch Algorithm
- Mathematical method for optimal structural alignment
- Finds best rotation and translation to superimpose two structures
- Uses Singular Value Decomposition (SVD) for optimal rotation matrix

### RMSD (Root Mean Square Deviation)
- Primary metric for structural similarity
- Lower RMSD = higher structural similarity
- Standard interpretation:
  - RMSD < 2.0 Å: High structural similarity ✓
  - RMSD 2.0-4.0 Å: Moderate structural similarity ○
  - RMSD > 4.0 Å: Low structural similarity ✗

---

## 📊 Output Features

### Console Output:
- RMSD value in Angstroms
- Number of compared residues
- Automatic structural similarity interpretation
- Additional statistics (mean, median, standard deviation)

### Visualization Dashboard (9 plots):
1. **3D Overlay**: Both structures superimposed
2. **Protein A Individual**: First protein structure
3. **Protein B Original**: Second protein (original)
4. **Protein B Aligned**: Second protein (aligned)
5. **Transformation**: Before/after alignment comparison
6. **Side View**: Side perspective with distance lines
7. **Distance Histogram**: Distribution of atomic distances
8. **Distance per Residue**: Distance along protein chain
9. **Cumulative Distribution**: Cumulative distance probability

---

## 🖼️ Program Output Examples

### Console Results
![Console Output](images/Console_Results.jpg)
<br>
*Console output showing RMSD calculation and similarity interpretation*

### Similar Structures Comparison
![Similar Structures](images/Compare_Similar.jpg)
<br>
*Visualization dashboard for highly similar protein structures (RMSD < 2.0 Å)*

### Different Structures Comparison  
![Different Structures](images/Compare_Different.jpg)
<br>
*Visualization dashboard for structurally different proteins (RMSD > 4.0 Å)*

---

## 🔍 Key Analysis Features

### Conserved Regions Identification
- **Conserved regions**: Parts of protein that cannot change (like engine in car)
- **Variable regions**: Parts that can vary without affecting function (like car color)
- Program automatically identifies conserved regions through distance analysis

### Structural Similarity Assessment
- Quantitative RMSD measurement
- Qualitative similarity interpretation
- Regional variation analysis
- Functional implication insights

---

## 🛠️ Technologies Used
- Python 3.x
- NumPy for mathematical computations
- Matplotlib for visualization
- PDB file parsing and processing
- Structural bioinformatics algorithms

---

## 📜 License [ [![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) ]
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

Special thanks to:

- The structural bioinformatics community for developing fundamental algorithms like Kabsch alignment
- Protein Data Bank (PDB) for providing standardized structural data format
- Scientific researchers who established RMSD as standard metric for structural comparison
- Python scientific computing community for excellent libraries (NumPy, Matplotlib)

---

## 🔬 Scientific Applications

This tool is particularly useful for:
- Protein structure-function relationship studies
- Homology modeling and comparative modeling
- Protein engineering and design
- Evolutionary biology research
- Drug discovery and target identification
- Bioinformatics education and training

The development of this project follows established principles in structural bioinformatics and implements industry-standard algorithms for protein structure comparison and analysis.
