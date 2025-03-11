import requests
from Bio import SeqIO
from io import StringIO

url = "https://mobidb.org/browse?limit=15&projection=level,acc,name,organism,length,prediction-disorder-mobidb_lite.content_fraction,prediction-disorder-alphafold.content_fraction&level=2&sort_asc=_id"



output_file = "microproteins.fasta"

# Clear the output file before starting.
with open(output_file, "w") as f:
    pass

last_id = None
page = 1

while True:
    if last_id:
        params["_id_gt"] = last_id
    else:
        params.pop("_id_gt", None)
    
    print("------------------------------------------------")
    print(f"Fetching page {page} with params: {params}")
    response = requests.get(url, params=params)
    print(f"HTTP status code: {response.status_code}")
    
    try:
        result = response.json()
    except Exception as e:
        print("Error decoding JSON:", e)
        print("Response text:", response.text)
        break
    
    fasta_data = result.get("data", "")
    metadata = result.get("metadata", {})
    last_id = metadata.get("last_id")
    
    print(f"Page {page}: Received {len(fasta_data)} characters of FASTA data.")
    print("Metadata:", metadata)
    
    if not fasta_data:
        print("No FASTA data returned on this page. Ending loop.")
        break
    
    # Parse the FASTA data.
    handle = StringIO(fasta_data)
    try:
        records = list(SeqIO.parse(handle, "fasta"))
    except Exception as e:
        print("Error parsing FASTA data:", e)
        break
    
    print(f"Page {page}: Parsed {len(records)} sequences.")
    
    # Append these records to the output file.
    with open(output_file, "a") as out_handle:
        SeqIO.write(records, out_handle, "fasta")
    print(f"Page {page}: Written {len(records)} sequences to {output_file}.")
    
    page += 1
    
    # If no more pages (last_id is None), break the loop.
    if not last_id:
        print("No last_id in metadata. Ending pagination.")
        break

print("Download complete. Sequences saved to", output_file)
