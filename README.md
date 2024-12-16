# Project SamAnalysis

## Description
This project analyzes SAM (Sequence Alignment/Map) files to extract statistics on alignments. It includes:
- A Bash script to validate the authenticity of SAM files (`veracity_test.sh`).
- A Python script to analyze and extract key biological information (such as the number of mapped reads, MAPQ scores, alignment types) and produce statistics on the alignments as well as distribution graphs.

## Prerequisites
In order to run our script `script_sam_analysis.py` optimally, here are the configurations and packages the script uses:

- **Operating System**: Linux or Windows with Bash.
- **Python version**: 3.8+.
- **Python Libraries**:
  - `matplotlib` for visualization.
  - `sys` for argument handling.

## Features
1. **Validate the structure of a SAM file (Bash script)**:  
    - Verifies if the file is valid and not a directory.  
    - Checks for the presence of required headers starting with `@`.
   
2. **Analyze mapped, unmapped, partially mapped alignments, and read pairs**.

3. **Generate detailed statistics**:
   - Distribution of FLAGS.
   - Distribution of MAPQ scores.
   - Distribution by chromosomes.

4. **Visualization**:  
    - Generates a histogram of the MAPQ score distribution.

## Requirements
- **Python 3.x**  
- Standard Python modules: `collections`, `sys`

## Installation
Clone this Git repository and navigate to the directory:

```bash
git clone https://github.com/Anysdn/ProjetSamAnalysis.git
cd ProjetSAMAnalyzer
```

Ensure Python 3 and the dependencies are installed:
```bash
pip install matplotlib
```

Make the scripts executable:
```bash
chmod +x scripts/veracity_test.sh
chmod +x scripts/script_sam_analysis.py
```
## Usage
Step 1: Run the Bash script to test the authenticity of the SAM file
```bash
veracity_test.sh <file.sam>
```
Step 2: Analyze the SAM file
```bash
python3 script_sam_analysis.py <file.sam>
```
## Example Output

```yaml
=== Results of Analysis ===
Total reads: 351330
Mapped reads: 350015
Unmapped reads: 1315
MAPQ score distribution:
  MAPQ [0-9]: 35068 reads
  MAPQ [60-69]: 301109 reads
CIGAR alignment types:
  simple_match: 350015 reads
```

A histogram file (mapq_distribution.png) will also be generated in the current directory.

## Contribution

Contributions are welcome! Please follow these steps:

1. Fork this repository.
2. Create a new branch:

    ```bash
    git checkout -b feature/new-feature
    ```

3. Make your changes and test them.
4. Submit a Pull Request with a clear description.



## Licence
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**AUTHOR**:
ANIA SAIDANI - Project carried out as part of the UE System (Master's in Bioinformatics)
ania.saidani@etu.umontpellier.fr

