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

Below is a quick example to get you started with GCMSnote. This demonstrates how to load your GC-MS data and perform metabolite annotation:
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

#### Contact
For questions, suggestions, or feedback, please contact boris.minasenko@emory.edu