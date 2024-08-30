import susie_fine_mapping_utils as sfm
import pandas as pd
import os
import subprocess

INPUTS_DIR = "/local/home/ssigfstead/GWAS/PAINTOR_inputs"
LD_SCRIPT =  "/local/home/ssigfstead/tools/PAINTOR_V3.0/PAINTOR_Utilities/CalcLD_1KG_VCF.py"
OVERLAP_SCRIPT = "/local/home/ssigfstead/tools/PAINTOR_V3.0/PAINTOR_Utilities/AnnotateLocus.py"
PAINTOR_SCRIPT = "/local/home/ssigfstead/tools/PAINTOR_V3.0/PAINTOR"

def reformat_locus(rsid, flank_size, save=True):
    locus_df = pd.read_csv(
        os.path.join(sfm.LOCI_DIR, f"locus_{rsid}_{flank_size}.tsv"), sep="\t"
    )

    locus_df["ZSCORE"] = locus_df["BETA"] / locus_df["SE"]
    if save:
        locus_df.to_csv(os.path.join(INPUTS_DIR, f"locus_{rsid}_{flank_size}"), sep=" ", index=None)
    return locus_df


if __name__ == "__main__":
    loci = sfm.get_loci("/local/home/ssigfstead/GWAS/filtered_snps_threshold=0.5.csv", chr_col="chr", pos_col="pos", rsid_col="snp")

    with open("/local/home/ssigfstead/GWAS/PAINTOR_inputs/loci.txt", "w") as f:
        for locus in loci:
            f.write(f"locus_{locus[3]}_{locus[4]}.processed\n")

    subprocess.run(
        [
            PAINTOR_SCRIPT,
            "-in",
            "/local/home/ssigfstead/GWAS/PAINTOR_inputs/",
            "-input",
            "/local/home/ssigfstead/GWAS/PAINTOR_inputs/loci.txt",
            "-Zhead",
            "ZSCORE",
            "-LDname",
            "ld",
            "-annotations",
            "brain_differentially_expressed_enhancers.bed,"
            "GenCode.CDS.hg19,"
            "GenCode.exon.hg19,"
            "GenCode.gene.hg19,"
            "GenCode.Selenocysteine.hg19,"
            "GenCode.start_codon.hg19,"
            "GenCode.stop_codon.hg19",
            "GenCode.transcript.hg19,"
            "GenCode.UTR.hg19,"
            "fBrain-DS11872.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS11877.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS14464.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS14717.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS14718.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS14803.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS14815.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS15453.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS16302.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS20221.hg19.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS20226.hg19.hotspot.twopass.fdr0.05.merge.bed,"
            "fBrain-DS20231.hg19.hotspot.twopass.fdr0.05.merge.bed,"
            "All.TFBS.bed",
            "-out",
            "/local/home/ssigfstead/GWAS/PAINTOR_outputs",
            "-mcmc",
        ],
        text=True,
    )