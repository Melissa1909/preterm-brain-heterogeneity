import os 
import pandas as pd


def read_process_output(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Step 2: Find the start of the table
    start_index = None
    end_index = None

    for i, line in enumerate(lines):
        if "Data for visualizing the conditional effect of the focal predictor:" in line:
            start_index = i + 2  # Table starts two lines below the header
        if "********************" in line and start_index is not None:
            end_index = i
            break

    # Step 3: Extract the table data
    table_lines = lines[start_index:end_index]
    table_data = [line.strip().split() for line in table_lines]

    # Step 4: Convert to a DataFrame
    columns = ["GA", "SES_at_birth", "PC1_CT"]
    df = pd.DataFrame(table_data, columns=columns)

    # Convert columns to numeric types
    df = df.apply(pd.to_numeric, errors='coerce').dropna(how='all')

    return df