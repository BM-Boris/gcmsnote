import pandas as pd
import numpy as np
from tqdm import tqdm
import pkg_resources


def load_library_data(library_type: str = "lc") -> pd.DataFrame:
    """
    Load the appropriate library file depending on the library_type.

    Parameters
    ----------
    library_type : {"gc", "lc"}
        - "gc": loads data/gcms_lib.csv (GC-MS library with RT).
        - "lc": loads hmdb.csv (LC-MS library with monoisotopic masses).

    Returns
    -------
    pandas.DataFrame
        Loaded library table.
    """
    if library_type == "gc":
        filename = "data/gcms_lib.csv"
    elif library_type == "lc":
        filename = "data/hmdb.csv"
    else:
        raise ValueError("library_type must be 'gc' or 'lc'")

    data_path = pkg_resources.resource_filename(__name__, filename)
    return pd.read_csv(data_path, low_memory=False)


def pre_process(features: pd.DataFrame, lib_gc: pd.DataFrame):
    """
    Preprocess feature and GC library tables for GC-MS matching.

    Parameters
    ----------
    features : pandas.DataFrame
        Feature table with m/z and RT in the first two columns.
    lib_gc : pandas.DataFrame
        GC-MS library table, first two columns are m/z and RT.

    Returns
    -------
    tuple of numpy.ndarray
        (feature_mz, feature_rt, library_mz, library_rt)
    """
    features.dropna(subset=[features.columns[0], features.columns[1]], inplace=True)

    mz_array = features.iloc[:, 0].to_numpy().astype("float64")
    time_array = features.iloc[:, 1].to_numpy().astype("float64")
    lib_mz_array = lib_gc.iloc[:, 0].to_numpy().astype("float64")
    lib_time_array = lib_gc.iloc[:, 1].to_numpy().astype("float64")

    return mz_array, time_array, lib_mz_array, lib_time_array


def match_meta(
    data: str,
    sep: str = "\t",
    mz_col: str = "mz",
    rt_col: str = "rt",
    save: str | None = None,
    shift: str | float = "auto",
    time_diff: float = 0.05,
    mz_diff: float = 5e-6,
    time_range: float = 2.0,
    ngroup: int = 3,
    type: str = "gc",
    ion_mode: str = "pos",
    adducts: list[str] | None = None,
    mass_tol: float = 1e-6,
) -> pd.DataFrame:
    """
    Annotate GC-MS or LC-MS data by matching features to a reference library.

    Two modes are supported:
      - GC mode (type="gc"): uses a GC-MS library with retention times.
        Matching is performed using both m/z and RT, followed by grouping
        of peaks belonging to the same chemical.
      - LC mode (type="lc"): uses an HMDB-based LC-MS library with
        monoisotopic masses. Matching is performed using m/z only, with
        optional support for multiple adducts (e.g., [M+H]+, [M+Na]+).
        No RT or grouping is used.

    LC behavior:
      - 'annotation': best hit (single metabolite name)
      - 'adduct': best adduct for that hit
      - 'ppm_error': ppm error for the best hit
      - 'candidates':
            * None if there is only a single hit for that mass
            * a list of dicts only if there are multiple candidates
              with the same neutral mass as the best hit:
              [{"name", "adduct", "ppm"}, ...]

    Parameters
    ----------
    data : str
        Path to the input feature table (CSV/TSV).
    sep : str, default "\\t"
        Column separator for the input file.
    mz_col : str, default "mz"
        Name of the m/z column in the input feature table.
    rt_col : str, default "rt"
        Name of the retention time column (GC mode only).
    save : str or None, default None
        If provided, save the annotated table to this path (CSV).
    shift : {"auto", float}, default "auto"
        GC mode only. If "auto", infer RT shift from 4,4'-DDE. If a number,
        use that value (in seconds) as a fixed RT shift.
    time_diff : float, default 0.05
        GC mode only. Maximum allowed relative RT difference for a match.
    mz_diff : float, default 5e-6
        Maximum allowed relative m/z difference for a match
        (used in both GC and LC modes).
    time_range : float, default 2.0
        GC mode only. RT window (in seconds) within which peaks are grouped.
    ngroup : int, default 3
        GC mode only. Minimum group size required to keep a group.
    type : {"gc", "lc"}, default "gc"
        Matching mode:
          - "gc": GC-MS with RT and grouping.
          - "lc": LC-MS with exact masses and adducts.
    ion_mode : {"pos", "neg"}, default "pos"
        LC mode only. Ionization mode:
          - "pos": positive ESI, uses [M+H]+, [M+Na]+, [M+NH4]+, [M+K]+ by default.
          - "neg": negative ESI, uses [M-H]- and [M+Cl]- by default.
    adducts : list of str or None, default None
        LC mode only. Subset of adducts to use. If None, all defaults for the
        given ion_mode are used. Possible values (depending on ion_mode):
          - positive: ["[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+"]
          - negative: ["[M-H]-", "[M+Cl]-"]
    mass_tol : float, default 1e-6
        LC mode only. Absolute tolerance (in Da) used to decide whether
        two candidates share the same neutral mass and should be treated
        as isomeric candidates for the same feature.

    Returns
    -------
    pandas.DataFrame
        - LC mode:
            Columns include:
              * mz_col
              * annotation (best hit)
              * adduct (best adduct)
              * ppm_error (for the best hit)
              * candidates:
                    None if only one unique mass candidate,
                    or a list of dicts if multiple candidates share
                    the same neutral mass as the best hit.
        - GC mode:
            Original m/z and RT plus:
              * annotations
              * notes
              * second_annotations
              * second_notes
            Filtered to keep only the best groups per chemical.
    """
    features = pd.read_csv(data, sep=sep)

    # ========================= LC MODE =========================
    if type == "lc":
        # Keep only m/z column from the feature table
        features = features[[mz_col]].astype("float64")
        mz_array = features[mz_col].to_numpy()

        # Load LC library (HMDB-based)
        lib_lc = load_library_data("lc")

        # Infer mass and name column names
        mass_candidates = [
            "monisotopic_molecular_weight",
            "monoisotopic_mass",
            "exact_mass",
            "mz",
            "mass",
            "MASS",
        ]
        name_candidates = ["name", "Name", "NAME"]

        mass_col = next((c for c in mass_candidates if c in lib_lc.columns), None)
        name_col = next((c for c in name_candidates if c in lib_lc.columns), None)

        if mass_col is None:
            raise ValueError(
                "No suitable mass column found in hmdb.csv. "
                f"Available columns: {list(lib_lc.columns)}"
            )
        if name_col is None:
            raise ValueError(
                "No 'name' column found in hmdb.csv. "
                f"Available columns: {list(lib_lc.columns)}"
            )

        lib_lc[mass_col] = pd.to_numeric(lib_lc[mass_col], errors="coerce")
        lib_lc = lib_lc.dropna(subset=[mass_col])

        lib_mz_array = lib_lc[mass_col].to_numpy()
        lib_name_array = lib_lc[name_col].astype(str).to_numpy()

        # Adduct mass shifts (in Da)
        proton = 1.007276
        mass_Na = 22.989218
        mass_NH4 = 18.033823
        mass_K = 38.963158
        mass_Cl = 34.969402

        if ion_mode == "pos":
            default_adducts = {
                "[M+H]+": (proton, 1),
                "[M+Na]+": (mass_Na, 1),
                "[M+NH4]+": (mass_NH4, 1),
                "[M+K]+": (mass_K, 1),
            }
        elif ion_mode == "neg":
            default_adducts = {
                "[M-H]-": (-proton, 1),
                "[M+Cl]-": (mass_Cl, 1),
            }
        else:
            raise ValueError("ion_mode must be 'pos' or 'neg'")

        # Restrict to a subset of adducts if requested
        if adducts is not None:
            used_adducts = {k: v for k, v in default_adducts.items() if k in adducts}
            if not used_adducts:
                raise ValueError(
                    f"None of the requested adducts {adducts} are valid. "
                    f"Allowed adducts: {list(default_adducts.keys())}"
                )
        else:
            used_adducts = default_adducts

        annotations: list[str | None] = []
        adduct_hits: list[str | None] = []
        ppm_best: list[float | None] = []
        candidates_all: list[list[dict] | None] = []

        # Main LC matching loop
        for mz in tqdm(mz_array):
            best_ppm = None
            best_name = None
            best_adduct = None
            best_mass = None
            candidates_raw: list[dict] = []

            for adduct_name, (mass_shift, charge) in used_adducts.items():
                # Theoretical m/z for each metabolite in the library under this adduct
                theo_mz = (lib_mz_array + mass_shift) / charge
                rel_error = np.abs((theo_mz - mz) / theo_mz)

                matched_indices = np.where(rel_error < mz_diff)[0]
                if matched_indices.size == 0:
                    continue

                for idx in matched_indices:
                    ppm = float(rel_error[idx] * 1e6)
                    name = lib_name_array[idx]
                    neutral_mass = float(lib_mz_array[idx])
                    candidates_raw.append(
                        {
                            "name": name,
                            "adduct": adduct_name,
                            "ppm": ppm,
                            "mass": neutral_mass,
                        }
                    )

                    if best_ppm is None or ppm < best_ppm:
                        best_ppm = ppm
                        best_name = name
                        best_adduct = adduct_name
                        best_mass = neutral_mass

            if best_name is None:
                annotations.append(None)
                adduct_hits.append(None)
                ppm_best.append(None)
                candidates_all.append(None)
            else:
                # Keep only candidates that share the same neutral mass as the best hit
                same_mass_candidates = [
                    c
                    for c in candidates_raw
                    if best_mass is not None and abs(c["mass"] - best_mass) < mass_tol
                ]

                # If there is only one such candidate, we consider it unambiguous -> no candidates list
                if len(same_mass_candidates) <= 1:
                    candidates_all.append(None)
                else:
                    # Drop the internal 'mass' field from output
                    for c in same_mass_candidates:
                        c.pop("mass", None)
                    candidates_all.append(same_mass_candidates)

                annotations.append(best_name)
                adduct_hits.append(best_adduct)
                ppm_best.append(best_ppm)

        features["annotation"] = annotations
        features["adduct"] = adduct_hits
        features["ppm_error"] = ppm_best
        features["candidates"] = candidates_all

        out = features.dropna(subset=["annotation"]).reset_index(drop=True)

        if save is not None:
            out.to_csv(save, index=False)

        return out

    # ========================= GC MODE =========================
    if type != "gc":
        raise ValueError("type must be 'gc' or 'lc'")

    # GC mode: use m/z + RT and the GC library
    features = features[[mz_col, rt_col]].astype("float64")

    lib_gc = load_library_data("gc")
    mz_array, time_array, lib_mz_array, lib_time_array = pre_process(features, lib_gc)

    annotations = np.empty(len(mz_array), dtype=object)
    second_annotations = np.empty(len(mz_array), dtype=object)
    notes = np.empty(len(mz_array), dtype=object)
    second_notes = np.empty(len(mz_array), dtype=object)

    # Convert library RT to seconds
    lib_time_array = lib_time_array * 60

    # Auto RT shift using 4,4'-DDE
    # Auto RT shift using 4,4'-DDE
    if shift == "auto":
        tol = 5e-6  # relative m/z tolerance
        mz = np.asarray(mz_array)
        rt = np.asarray(time_array)
    
        # mz1: 4,4'-DDE
        row1 = lib_gc[(lib_gc.Name == "4,4'-DDE ") & (lib_gc.note == "mz1")].iloc[0]
        mz1, rt1_ref = row1.mz, row1.time * 60
    
        dmz1 = np.abs((mz - mz1) / mz1)
        cand_idx = np.where(dmz1 <= tol)[0]
    
        if cand_idx.size == 0:
            print("No 4,4'-DDE mz1 found in data, using shift = 0")
            shift = 0
        else:
            cand_shifts = rt[cand_idx] - rt1_ref
            base_i = int(np.argmin(dmz1[cand_idx]))
            shift = cand_shifts[base_i]
            confirmed = None
            confirmed_i = base_i
    
            def anchor_shift(mz_a, rt_a):
                dmz = np.abs((mz - mz_a) / mz_a)
                idx = np.where(dmz <= tol)[0]
                if idx.size == 0:
                    return None
                j = idx[int(np.argmin(dmz[idx]))]
                return rt[j] - rt_a
    
            def match_anchor(a_shift):
                if a_shift is None:
                    return None
                diffs = np.abs(cand_shifts - a_shift)
                hits = np.where(diffs <= time_range)[0]
                return int(hits[0]) if hits.size > 0 else None
    
            # 1) C13_pp_DDE: can override base candidate
            row_c13 = lib_gc[lib_gc.Name == "C13_pp_DDE"].iloc[0]
            c13_shift = anchor_shift(row_c13.mz, row_c13.time * 60)
            if c13_shift is None:
                print("C13_pp_DDE not found in data")
            else:
                k = match_anchor(c13_shift)
                if k is not None:
                    shift = cand_shifts[k]
                    confirmed = "C13"
                    confirmed_i = k
                    print("Shift chosen by mz1 candidate confirmed by C13_pp_DDE")
                else:
                    print("C13_pp_DDE found in data, but no mz1 candidate matches its RT shift")
    
            # 2) 4,4'-DDE mz0
            row0 = lib_gc[(lib_gc.Name == "4,4'-DDE ") & (lib_gc.note == "mz0")].iloc[0]
            mz0_shift = anchor_shift(row0.mz, row0.time * 60)
            if mz0_shift is None:
                print("4,4'-DDE mz0 not found in data")
            else:
                if confirmed == "C13":
                    # only check consistency with C13-confirmed candidate
                    if abs(mz0_shift - cand_shifts[confirmed_i]) <= time_range:
                        print("Shift also consistent with 4,4'-DDE mz0")
                    else:
                        print("4,4'-DDE mz0 shift not consistent with C13-based shift")
                else:
                    # C13 did not confirm; mz0 can choose any candidate
                    k = match_anchor(mz0_shift)
                    if k is not None:
                        shift = cand_shifts[k]
                        confirmed = "mz0"
                        confirmed_i = k
                        print("Shift chosen by mz1 candidate confirmed by 4,4'-DDE mz0")
                    else:
                        print("4,4'-DDE mz0 found in data, but no mz1 candidate matches its RT shift")
    
            if confirmed is None:
                if cand_idx.size > 1:
                    print("Multiple mz1 candidates; using best m/z match (no confirmation)")
                else:
                    print("Single mz1 candidate; using best m/z match")
    
            print(f"Shift = {shift} seconds based on 4,4'-DDE mz1")



    # Apply RT shift
    time_array = time_array - float(shift)

    # GC matching loop
    for i in tqdm(range(len(mz_array))):
        dmz = mz_array[i] - lib_mz_array
        dt = time_array[i] - lib_time_array
        distances = dmz**2 + 1e-7 * (dt**2)
        closest_indices = np.argpartition(distances, 1)[:5]

        cond_mz = (
            np.abs(lib_mz_array[closest_indices[0]] - mz_array[i])
            / lib_mz_array[closest_indices[0]]
            < mz_diff
        )
        cond_time = (
            np.abs(lib_time_array[closest_indices[0]] - time_array[i])
            / lib_time_array[closest_indices[0]]
            < time_diff
        )

        if cond_mz and cond_time:
            annotations[i] = lib_gc.iloc[closest_indices[0]]["Name"]
            second_annotations[i] = lib_gc.iloc[closest_indices[1]]["Name"]
            notes[i] = lib_gc.iloc[closest_indices[0]]["note"]
            second_notes[i] = lib_gc.iloc[closest_indices[1]]["note"]

    features["annotations"] = annotations
    features["notes"] = notes
    features["second_annotations"] = second_annotations
    features["second_notes"] = second_notes

    feat_nonnull = features.dropna(subset=["annotations"])
    feat_nonnull = feat_nonnull.sort_values(by=["annotations", rt_col])

    data_gc = feat_nonnull.reset_index().drop("index", axis=1)
    best_group_indices: list[int] = []

    # Group peaks for each chemical and keep best group(s)
    for chemical in data_gc["annotations"].unique():
        chem_data = data_gc[data_gc["annotations"] == chemical]
        groups: list[tuple[list[int], bool]] = []

        i = 0
        while i < len(chem_data):
            group_indices = [chem_data.index[i]]
            has_m0 = "mz0" in chem_data.iloc[i]["notes"]

            j = i + 1
            while (
                j < len(chem_data)
                and abs(chem_data.iloc[i][rt_col] - chem_data.iloc[j][rt_col])
                <= time_range
            ):
                group_indices.append(chem_data.index[j])
                if "mz0" in chem_data.iloc[j]["notes"]:
                    has_m0 = True
                j += 1

            if len(group_indices) >= ngroup:
                groups.append((group_indices, has_m0))

            i = j

        if groups:
            # Prefer groups that contain m0 if any exist
            if any(has_m0 for _, has_m0 in groups):
                groups = [grp for grp in groups if grp[1]]

            # Among remaining groups, keep only those with maximum size
            max_size = max(len(grp[0]) for grp in groups)
            groups = [grp for grp in groups if len(grp[0]) == max_size]

            for grp, _ in groups:
                best_group_indices.extend(grp)

    best_group_indices = sorted(set(best_group_indices))

    data_gc = (
        data_gc.loc[best_group_indices]
        .sort_values(by=["annotations", rt_col])
        .reset_index()
        .drop("index", axis=1)
    )

    if save is not None:
        data_gc.to_csv(save, index=False)

    return data_gc
