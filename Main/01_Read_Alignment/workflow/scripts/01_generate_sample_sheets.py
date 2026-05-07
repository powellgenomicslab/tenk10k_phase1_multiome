# qsub -I -q copyq -l storage=gdata/<HPC_PROJECT>+gdata/<HPC_PROJECT_2>
# mm activate airtable-env

# Script to query the airtable via API and get the data required for automated samplesheet generation

import os
from pyairtable import Api
import json
import pandas as pd


def get_access_token(token_json):
    # api_keys = {}
    with open(token_json, "r") as f:
        api_key = json.loads(f.read())["airtable"]
    return api_key


def update_samplesheet_master_table():
    api = Api(
        get_access_token(
            "/g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/.credentials/token.json"
        )
    )
    sequencing_table = api.table("<AIRTABLE_APP_ID>", "Sequencing")
    # get bioinformatics view data
    sequencing_records = sequencing_table.all(view="Grid view")
    sequencing_df = pd.DataFrame(
        [record["fields"] for record in sequencing_records],
    )
    sample_sheet_info = (
        sequencing_df.loc[
            :,
            [
                "Sequencing ID",
                "Index",
                "FluidX ID (from FluidX tube ID)",
                "Targeted Cells# (from Capture ID)",
            ],
        ]
        .rename(
            columns={
                "FluidX ID (from FluidX tube ID)": "FluidX ID",
                "Targeted Cells# (from Capture ID)": "Targeted Cell Number",
            }
        )
        .dropna()
    )
    sample_sheet_info["Targeted Cell Number"] = [
        i[0] for i in sample_sheet_info["Targeted Cell Number"]
    ]  # unlist
    sample_sheet_info["FluidX ID"] = (
        sample_sheet_info["FluidX ID"]
        .apply(lambda x: "".join(x) if isinstance(x, list) else x)
        .str.strip()
    )
    sample_sheet_info["Pool"] = sample_sheet_info["Sequencing ID"].str.extract(
        r"^(?:[^_]*_){3}([^_]+)"
    )
    # this reflects the order in the airtable spreadsheet, it corresponds to the order in which the samples were processed throughout the week
    sample_sheet_info["Weekly capture number"] = (
        sample_sheet_info["Sequencing ID"]
        .str.extract(r"^(?:[^_]*_){4}(.+)$")[0]
        .str.replace("_R|R_", "", regex=True)
        .astype(int)
    )

    sample_sheet_info = sample_sheet_info.sort_values(
        by=[
            "Targeted Cell Number",
            "FluidX ID",
            "Weekly capture number",
            "Sequencing ID",
        ]
    ).reset_index(drop=True)

    sample_sheet_info["pool_replicate_number"] = (
        sample_sheet_info.groupby("Pool").cumcount() + 1
    )
    # now construct the pool id's to use for cellranger \
    sample_sheet_info["sample_id_cellranger"] = (
        sample_sheet_info["Pool"]
        + "_"
        + sample_sheet_info["pool_replicate_number"].astype(str)
    )
    print(sample_sheet_info)
    sample_sheet_info.to_csv(
        "/g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/sample_sheets/mastertable.csv"
    )


if __name__ == "__main__":
    # get the latest data from airtable
    update_samplesheet_master_table()

    # Example: generate a sample sheet for a given flowcell using the FluidX tube ID (fsid)
    # make_samplesheet(fsid = '<FSID>', output_name = '/g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/sample_sheets/automated/sample_sheet_<DATE>_<FLOWCELL_ID>.csv')
