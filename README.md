## **GCMSnote: Simplified GC-MS Data Annotation**

_GCMSnote_ is a specialized Python library designed to streamline the annotation of Gas Chromatography-Mass Spectrometry (GC-MS) data. By leveraging an internal metabolite library, GCMSnote efficiently matches measured metabolites, offering a straightforward solution for researchers and analysts in the field of metabolomics.

### **Features**

- **Automated Metabolite Annotation**: Match your GC-MS data against a curated library with ease.
- **Customizable Matching Criteria**: Tailor the matching process to your specific research needs by adjusting key parameters such as mass-to-charge ratio (m/z) differences, retention time shifts, and etc.
- **Grouping of Similar Annotations:**: Automatically group metabolites based on similarity to enhance the clarity and utility of your analysis.

#### Installation

Install GCMSnote directly from GitHub:
```bash
pip install git+https://github.com/BM-Boris/gcmsnote.git
```

#### Getting Started

Below is a quick example to get you started with GCMSnote:
```python
from gcmsnote import match_meta

# Path to your GC-MS data file
data_path = 'path/to/your/data.csv'
output_path = 'path/for/output.csv'

# Perform metabolite annotation
annotated_data = match_meta(data_path, sep='\t', mz_col='mz', rt_col='rt', shift=16, save=output_path)

# Explore your annotated data
print(annotated_data.head())

```

#### Parameters:
```python
"""
- data (str): Path to the input CSV file containing GC-MS data.
- sep (str): Separator used in the input CSV file. Defaults to '\t'.
- mz_col (str): Column name for mass-to-charge ratio (m/z) in the input data. Defaults to 'mz'.
- rt_col (str): Column name for retention time in the input data. Defaults to 'rt'.
- save (str): Path to save the annotated and grouped data as a CSV file. If None, the data is not saved. Defaults to None.
- shift (str or float): Adjustment to the retention time to align with the library. Defaults to 'auto' - calculates the shift based on 4,4'-DDE.
- time_diff (float): Maximum allowed difference in retention time for a match. Defaults to 0.05.
- mz_diff (float): Maximum allowed m/z difference for a match. Defaults to 5e-6.
- time_range (float): Time range within which metabolites are considered for grouping. Defaults to 2.
- ngroup (int): Minimum number of metabolites required to form a group. Defaults to 3.

"""
```

#### Contact
For questions, suggestions, or feedback, please contact boris.minasenko@emory.edu
