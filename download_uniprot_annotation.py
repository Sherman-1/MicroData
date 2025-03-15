import requests

import requests
import re

def download_uniprot_paginated(url, output_file):
    session = requests.Session()
    with open(output_file, 'wb') as f:
        while url:
            print(f'Preparing download from: {url}')
            response = session.get(url, stream=True)

            if response.status_code != 200:
                raise Exception(f"Failed to retrieve data: {response.status_code}")

            for chunk in response.iter_content(chunk_size=8192):
                file_written = f.write(chunk)

            link_header = response.headers.get('Link', None)
            url = None  # Reset url
            if link_header:
                match = re.search(r'<([^>]+)>; rel="next"', link_header)
                url = match.group(1) if match else None

url = "https://rest.uniprot.org/uniprotkb/search?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Cgene_orf%2Corganism_id%2Clineage_ids%2Clineage%2Clength%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id%2Cabsorption%2Cannotation_score%2Ccc_caution%2Ckeyword%2Ccomment_count%2Cfeature_count%2Cuniparc_id%2Cft_intramem%2Cft_transmem%2Cft_chain%2Cft_peptide%2Cft_signal%2Cft_transit%2Cft_propep%2Cft_disulfid%2Cft_crosslnk%2Cft_carbohyd%2Cft_init_met%2Cft_lipid%2Cft_mod_res%2Ccc_ptm&format=tsv&query=%28%28length%3A%5B20+TO+100%5D%29%29&size=500"
output_file = "uniprot_annotation.tsv"

download_uniprot_paginated(url, output_file)

print("Download completed!")
