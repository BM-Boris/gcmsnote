import pandas as pd
import numpy as np
from tqdm import tqdm
import pkg_resources

def load_library_data():
    data_path = pkg_resources.resource_filename(__name__, 'data/gcms_lib.csv')
    return pd.read_csv(data_path)

lib = load_library_data()

def pre_process(feat,lib):
    
    feat.dropna(subset=[feat.columns[0],feat.columns[1]],inplace = True)

    mz_array = feat.iloc[:,0].to_numpy().astype('float64')
    time_array = feat.iloc[:,1].to_numpy().astype('float64')
    lib_mz_array = lib.iloc[:,0].to_numpy().astype('float64')
    lib_time_array = lib.iloc[:,1].to_numpy().astype('float64')
    
    return mz_array,time_array,lib_mz_array,lib_time_array

def match(data, sep='\t', mz_col='mz', rt_col='rt', save=None, shift=0, time_diff=0.05, mz_diff=5e-6, time_range=2, ngroup=3):
    """
    Annotates GC-MS data by matching measured metabolites to a library, and groups annotations by similarity.

    This function processes GC-MS feature data, matches it against an internal library, annotates the metabolites,
    and groups similar annotations. It allows for the adjustment of match criteria, including mass-to-charge ratio (m/z) differences,
    retention time shifts, and grouping thresholds.

    Parameters:
    - data (str): Path to the input CSV file containing GC-MS data.
    - sep (str): Separator used in the input CSV file. Defaults to '\t'.
    - mz_col (str): Column name for mass-to-charge ratio (m/z) in the input data. Defaults to 'mz'.
    - rt_col (str): Column name for retention time in the input data. Defaults to 'rt'.
    - save (str): Path to save the annotated and grouped data as a CSV file. If None, the data is not saved. Defaults to None.
    - shift (float): Adjustment to the retention time to align with the library. Defaults to 0.
    - time_diff (float): Maximum allowed difference in retention time for a match. Defaults to 0.05.
    - mz_diff (float): Maximum allowed m/z difference for a match. Defaults to 5e-6.
    - time_range (float): Time range within which metabolites are considered for grouping. Defaults to 2.
    - ngroup (int): Minimum number of metabolites required to form a group. Defaults to 3.

    Raises:
    - ValueError: If the input data or library cannot be processed.

    Returns:
    - DataFrame: A pandas DataFrame containing the original data with additional columns for annotations, notes,
      secondary annotations, and secondary notes, filtered to include only the best groups of annotations.
    """
    
    feat = pd.read_csv(data, sep=sep)
    feat=feat[[mz_col,rt_col]].astype('float64')

    mz_array,time_array,lib_mz_array,lib_time_array = pre_process(feat,lib)

    annotations = np.empty(len(mz_array), dtype=object)
    second_annotations = np.empty(len(mz_array), dtype=object)
    notes = np.empty(len(mz_array), dtype=object)
    second_notes = np.empty(len(mz_array), dtype=object)

    lib_time_array = lib_time_array*60
    time_array = time_array-shift

    for i in tqdm(range(len(mz_array))):
        dmz = mz_array[i] - lib_mz_array
        dt = time_array[i] - lib_time_array
        distances = dmz**2 + 1e-7*(dt**2)
        closest_indices = np.argpartition(distances, 1)[:5]
        if((np.abs(lib_mz_array[closest_indices[0]] - mz_array[i])/lib_mz_array[closest_indices[0]] < mz_diff) &
       (np.abs(lib_time_array[closest_indices[0]] - time_array[i])/lib_time_array[closest_indices[0]] < time_diff)):
            
            annotations[i] = lib.iloc[closest_indices[0]]['Name']
            second_annotations[i] = lib.iloc[closest_indices[1]]['Name']
            notes[i] = lib.iloc[closest_indices[0]]['note']
            second_notes[i] = lib.iloc[closest_indices[1]]['note']


    feat['annotations'] = annotations
    feat['notes'] = notes
    feat['second_annotations'] = second_annotations
    feat['second_notes'] = second_notes

    feat_ = feat.dropna(subset=['annotations'])
    feat_ = feat_.sort_values(by=['annotations',rt_col])
    
    data = feat_.reset_index().drop('index',axis=1)
    best_group_indexes = []
    
    for chemical in data['annotations'].unique():
        # Filter rows for the current chemical
        chemical_data = data[data['annotations'] == chemical]
    
        groups = []
    
        i = 0
        while i < len(chemical_data):
            group_indexes = [chemical_data.index[i]]
            has_m0 = 'mz0' in chemical_data.iloc[i]['notes']
    
            # Check the next rows to find a group within a range of 2 units
            j = i + 1
            while j < len(chemical_data) and abs(chemical_data.iloc[i][rt_col] - chemical_data.iloc[j][rt_col]) <= time_range:
                group_indexes.append(chemical_data.index[j])
                if 'mz0' in chemical_data.iloc[j]['notes']:
                    has_m0 = True
                j += 1
            group_size = len(group_indexes)
            # If more than two rows are in the group, add their indexes to the unique set
            if len(group_indexes) >= ngroup:
                groups.append((group_indexes, has_m0))
    
            i = j
    
        # Filter groups based on m0 and size
        if(groups!= []):
            if any(has_m0 for _, has_m0 in groups):
                groups = [grp for grp in groups if grp[1]]  # Keep only groups with m0
            max_size = max(len(grp[0]) for grp in groups)
            groups = [grp for grp in groups if len(grp[0]) == max_size]
            
            for grp, _ in groups:
                best_group_indexes.extend(grp)
    
    best_group_indexes = sorted(set(best_group_indexes))

    data = data.loc[best_group_indexes].sort_values(by=['annotations',rt_col]).reset_index().drop('index',axis=1)

    if save:
        data.to_csv(save, index=False)

    return data


















