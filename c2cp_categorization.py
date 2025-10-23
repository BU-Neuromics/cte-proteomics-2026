from litellm import completion
import os
import dotenv
import random
import pandas as pd
import re

dotenv.load_dotenv()

import json

with open('/restricted/projectnb/cteseq/projects/somascan/data/c2.cp.v2025.1.Hs.json') as f:
    c2cp = json.load(f)

print(f"{len(c2cp)} gene sets in the database")

system_prompt = {
        "content": ("You are an expect at molecular biology and genetics. "
                    "You can examine lists of genes and other basic gene set "
                    "information and categorize genes into high level categories "
                    "based on the genes function."
                    ),
        "role": "system"
}

user_prompt = """
Consider the following MSigDB gene set information:

Gene set name: {gset_name}

{gset}

Attempt to label the gene set into one of the following categories:

- Inflammation/ Immune System
- Cell Growith/ Survival
- ECM
- Neuronal
- Metabolism/ Mitochondria
- Synaptic
- Proteasome
- Vascular
- Ribosomal
- Other

The 'Other' category should include gene sets that don't fit well into
any of the other categories.

Respond with valid JSON output of the form:

{{
    "category": "your label",
    "rationale": "brief explanation of why you made the choice"
}}

Make sure to provide a rationale based on the content of the gene set
description.
"""

gsets_df = pd.read_excel("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/Unique_fgsea_pathways.xlsx")
gsets = gsets_df.to_dict(orient="records")
#gsets = gsets_df['Pathways'].tolist()
results = []
for Pathways in gsets:

    gset_prompt = user_prompt.format(
            gset_name=Pathways,
            gset=json.dumps(gsets)
    )

    response = completion(
      model="bedrock/anthropic.claude-3-5-sonnet-20240620-v1:0",
      messages=[
          system_prompt,
          {"content": gset_prompt, "role": "user"}
      ]
    )

    resp_json = response['choices'][0]['message']['content']
    
    resp_json_clean = resp_json.strip()

    match = re.search(r'\{.*\}', resp_json_clean, re.DOTALL)
    if match:
       resp_json_clean = match.group(0) 
    try:
        parsed = json.loads(resp_json_clean)
        category = parsed.get('category')
        rationale = parsed.get('rationale')
    except json.JSONDecodeError:
        category = None
        rationale = None
    results.append({
        "pathway": Pathways, 
        "response": resp_json,
        "category": category,
        "rationale": rationale
    })

df = pd.DataFrame(results)

df.to_csv('/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/PMI_imp1/output_fgsea_groups_bedrock.csv', index=False, encoding ='utf-8')
print(f"Saved {len(df)} rows")