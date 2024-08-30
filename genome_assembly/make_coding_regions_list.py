import sys
import pandas as pd

def main():

    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv("/local/home/ssigfstead/GWAS/genome_assembly/ncbi_dataset/data/GCF_000001405.25/genomic.gff", sep="\t", comment="#", header=None, names=columns)
    seqid_to_chr = {}
    for i, row in gff_df[gff_df["type"]=="region"].iterrows():
        if row["type"] == "region":
            all_attributes = row["attributes"].split(";")
            all_attributes_dict = {}
            for attribute in all_attributes:
                key, value = attribute.split("=")
                all_attributes_dict[key] = value
            if all_attributes_dict.get("genome") == "chromosome":
                seqid_to_chr[row["seqid"]] = all_attributes_dict["chromosome"]
    seqid_to_chr.pop('NC_000023.10')
    seqid_to_chr.pop('NC_000024.9')
    #chr_to_seqid = {float(v): k for k, v in seqid_to_chr.items()}
    print(seqid_to_chr)

    # Filter for exons (protein coding regions)
    exons = gff_df[gff_df['type'] == 'exon'].reset_index(drop=True)
    print(exons[exons['seqid'] == 'NC_000010.10'].head())
    exons['chr'] = exons['seqid'].map(seqid_to_chr)
    columns = ["chr", "seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    exons = exons[columns]
    print("after cols")
    print(exons[exons['seqid'] == 'NC_000010.10'].head())

    exons_bp_locations = pd.DataFrame(exons[['start', 'end', 'chr']])
    exons_bp_locations = exons_bp_locations[exons_bp_locations['chr'].notna()]
    print("after basepairs")
    print(exons_bp_locations[exons_bp_locations['chr'] == '10'].head())


    # Save these as csv for use in 1Mb_procedure_v3.py
    exons_bp_locations.to_csv('./genome_assembly/exon_regions_v2.csv', index=False)

    return



if __name__ == "__main__":
    main()
