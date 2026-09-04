# ------------------------------------------------------------
# Create required folders and download data if necessary
# ------------------------------------------------------------
from pathlib import Path
from urllib.request import Request, urlopen
import json


RESULTS_DIR = Path("Results")
DATA_DIR = Path("data")

GITHUB_DATA_API = (
    "https://api.github.com/repos/muharif/"
    "SysBioPhD2026/contents/Transcriptomics/data?ref=main"
)


# Create Results folder if it does not exist
if not RESULTS_DIR.exists():
    RESULTS_DIR.mkdir(parents=True)
    print("Created Results/ folder.")
else:
    print("Results/ folder already exists.")


# Helper function to download a GitHub folder recursively
def download_github_folder(api_url, local_dir):
    local_dir.mkdir(parents=True, exist_ok=True)

    request = Request(
        api_url,
        headers={"User-Agent": "SysBioPhD2026"}
    )

    with urlopen(request) as response:
        items = json.load(response)

    for item in items:
        destination = local_dir / item["name"]

        if item["type"] == "file":
            print(f"Downloading {item['name']}...")

            request = Request(
                item["download_url"],
                headers={"User-Agent": "SysBioPhD2026"}
            )

            with urlopen(request) as response:
                destination.write_bytes(response.read())

        elif item["type"] == "dir":
            download_github_folder(
                item["url"],
                destination
            )


# Download data only if the data folder does not already exist
if not DATA_DIR.exists():
    print("data/ folder not found.")
    print("Downloading exercise data from GitHub...")

    download_github_folder(
        GITHUB_DATA_API,
        DATA_DIR
    )

    print("Data download complete.")

else:
    print("data/ folder already exists. Skipping download.")
