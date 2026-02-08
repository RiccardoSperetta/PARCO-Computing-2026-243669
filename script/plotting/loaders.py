from pathlib import Path
import pandas as pd

def read_one_file(path: Path) -> pd.DataFrame:
    """Reads one .txt file and returns it as DataFrame"""
    df = pd.read_csv(
        path,
        sep=r'\s+',           # any amount of spaces
        header=None,
        names=[
            'total_time',
            'comm_time',
            'TEPS',
            'max_over_mean',
            'CV'
        ],
        comment='#',          # in case you have comment lines later
    )
    return df


def load_all_results(results_root: str = "./results") -> pd.DataFrame:
    """
    Finds all graph folders and all basic/hybrid .txt files.
    Adds columns: graph, variant, cores
    """
    root = Path(results_root).resolve()
    all_dfs = []

    # loop over each graph folder (graph_name_A, graph_name_B, ...)
    for graph_folder in root.iterdir():
        if not graph_folder.is_dir():
            continue

        graph_name = graph_folder.name

        # loop over every .txt file in this graph folder
        for txt_file in graph_folder.glob("*.txt"):
            stem = txt_file.stem  # "basic16", "hybrid64", ...

            if stem.startswith("basic"):
                variant = "basic"
                cores_str = stem.replace("basic", "")
            elif stem.startswith("hybrid"):
                variant = "hybrid"
                cores_str = stem.replace("hybrid", "")
            else:
                print(f"Skipping unknown file: {txt_file}")
                continue

            try:
                cores = int(cores_str)
            except ValueError:
                print(f"Cannot parse cores from: {stem}")
                continue

            # read file
            df = read_one_file(txt_file)

            # add metadata columns
            df['graph']     = graph_name
            df['variant']   = variant
            df['cores']     = cores
            df['filename']  = txt_file.name

            all_dfs.append(df)

    if not all_dfs:
        raise FileNotFoundError("No valid .txt files found!")

    combined = pd.concat(all_dfs, ignore_index=True)
    return combined