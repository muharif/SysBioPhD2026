from pathlib import Path
from urllib.request import Request, urlopen

# ------------------------------------------------------------
# Check/download co-expression network
# ------------------------------------------------------------

NETWORK_FILE = Path("Network/coexpression_network.txt")

NETWORK_URL = (
    "https://raw.githubusercontent.com/"
    "muharif/SysBioPhD2026/main/Network/coexpression_network.txt"
)

if not NETWORK_FILE.exists():
    print("coexpression_network.txt not found.")
    print("Downloading network file from GitHub...")

    NETWORK_FILE.parent.mkdir(parents=True, exist_ok=True)

    request = Request(
        NETWORK_URL,
        headers={"User-Agent": "SysBioPhD2026"}
    )

    with urlopen(request) as response:
        NETWORK_FILE.write_bytes(response.read())

    print("Network download complete.")

else:
    print("coexpression_network.txt already exists. Skipping download.")