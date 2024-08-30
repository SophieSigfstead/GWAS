import gwaslab as gl
import os
import pandas as pd
import subprocess
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as numpy2ri
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse
from tqdm import tqdm

SUMSTATS_DIR = "/local/home/ssigfstead/GWAS/gwas_1_and_2_summary_statistics_data"
REF_SEQ = "/local/home/ssigfstead/GWAS/reference_genomes/hg19.fa"
# TSV and VCF files obtained from the GWASLab tutorial at https://cloufield.github.io/GWASTutorial/finemapping_susie/
REF_RSID_TSV = (
    "/local/home/ssigfstead/GWAS/SNP_rsid_references/1kg_dbsnp151_hg19_auto.txt.gz"
)
REF_RSID_VCF = "/local/home/ssigfstead/GWAS/SNP_rsid_references/GCF_000001405.25"
LOCI_DIR = "/local/home/ssigfstead/GWAS/finemapping/finemapping_loci"
PLINK2 = "/local/home/ssigfstead/tools/plink2/plink2"
LD_DIR = "/local/home/ssigfstead/GWAS/linkage_disequilibrium"
REF_GENOTYPES = "/local/home/ssigfstead/GWAS/1KG_genotypes/all_phase3"
RESULTS_DIR = "/local/home/ssigfstead/GWAS/fine_mapping_results"


def load_sumstats(sumstats_name: str):
    path = os.path.join(SUMSTATS_DIR, sumstats_name)
    return gl.Sumstats(path, fmt="ssf")


def harmonize(sumstats: gl.Sumstats):
    # TODO: Check what the hell is going on with the reference files
    # and why we need to run harmonize twice
    sumstats.basic_check()
    sumstats.check_ref(ref_seq=REF_SEQ)
    sumstats.fix_id(fixid=True, forcefixid=True, overwrite=True)
    sumstats.harmonize(
        ref_seq=REF_SEQ,
        ref_rsid_tsv=REF_RSID_TSV,
        ref_rsid_vcf=REF_RSID_VCF,
        n_cores=24,
    )
    sumstats.fix_id(fixid=True, forcefixid=True, overwrite=True)
    sumstats.harmonize(ref_seq="/local/home/ssigfstead/GWAS/reference_genomes/hg19.fa",
                   ref_rsid_tsv=gl.get_path("1kg_dbsnp151_hg19_auto"),
                   ref_rsid_vcf="/local/home/ssigfstead/GWAS/SNP_rsid_references/GCF_000001405.25",
                   n_cores=24)


def get_loci(
    snp_file: str,
    flank_size: int = 5e5,
    chr_col: str = "CHR",
    pos_col: str = "POS",
    rsid_col: str = "rsID",
):
    df = pd.read_csv(snp_file)
    loci = []
    for i, row in df.iterrows():
        loci.append(
            (
                row[chr_col],
                row[pos_col] - flank_size,
                row[pos_col] + flank_size,
                row[rsid_col],
                flank_size,
            )
        )
    return loci


def output_loci_files(sumstats: gl.Sumstats, loci: list):
    for locus in loci:
        chr, start, end, rsid, flank_size = locus
        locus = sumstats.filter_value(f"CHR=={chr} & POS>{start} & POS<{end}")
        locus.fill_data(to_fill=["BETA"])
        locus.data = locus.data.dropna()
        assert locus.data["BETA"].isna().sum() == 0
        assert locus.data["SE"].isna().sum() == 0
        locus.data.where(locus.data["rsID"] != "<NA>").dropna().to_csv(
            os.path.join(LOCI_DIR, f"locus_{rsid}_{flank_size}.tsv"),
            sep="\t",
            index=None,
        )
        locus.data.where(locus.data["rsID"] != "<NA>").dropna()["rsID"].to_csv(
            os.path.join(LOCI_DIR, f"locus_{rsid}_{flank_size}.snplist"),
            sep="\t",
            index=None,
            header=None,
        )


def calculate_LD(locus, phased=True):
    _, _, _, rsid, flank_size = locus
    snplist_file = os.path.join(LOCI_DIR, f"locus_{rsid}_{flank_size}.snplist")
    out_file = os.path.join(LD_DIR, f"locus_{rsid}_{flank_size}_ld")
    phased_suffix = "phased" if phased else "unphased"
    subprocess.run(
        [
            PLINK2,
            f"--r-{phased_suffix}",
            "square",
            "--extract",
            snplist_file,
            "--out",
            out_file,
            "--pfile",
            REF_GENOTYPES,
        ],
        text=True,
    )
    subprocess.run(
        [
            PLINK2,
            f"--r2-{phased_suffix}",
            "square",
            "--extract",
            snplist_file,
            "--out",
            out_file,
            "--pfile",
            REF_GENOTYPES,
        ],
        text=True,
    )


def fine_map_locus(locus, phased=True, L=20, n=337126, min_abs_corr=0.5):
    susieR = importr("susieR")
    numpy2ri.activate()

    _, _, _, rsid, flank_size = locus
    phased_suffix = "phased" if phased else "unphased"

    locus_df = pd.read_csv(
        os.path.join(LOCI_DIR, f"locus_{rsid}_{flank_size}.tsv"), sep="\t"
    )
    assert locus_df["BETA"].isna().sum() == 0
    assert locus_df["SE"].isna().sum() == 0
    vars = pd.read_csv(
        os.path.join(LD_DIR, f"locus_{rsid}_{flank_size}_ld.{phased_suffix}.vcor1.vars"), sep="\t", header=None
    )
    assert set(vars[0]).issubset(set(locus_df["rsID"]))
    locus_df = locus_df[locus_df["rsID"].isin(vars[0])]
    assert len(vars) == len(locus_df)
    assert (list(vars[0]) == list(locus_df["rsID"]))
    locus_df = locus_df.reset_index(drop=True)
    assert locus_df["BETA"].isna().sum() == 0
    assert locus_df["SE"].isna().sum() == 0

    ld = pd.read_csv(
        os.path.join(LD_DIR, f"locus_{rsid}_{flank_size}_ld.{phased_suffix}.vcor1"),
        sep="\t",
        header=None,
    )
    R_df = ld.values
    ld2 = pd.read_csv(
        os.path.join(LD_DIR, f"locus_{rsid}_{flank_size}_ld.{phased_suffix}.vcor2"),
        sep="\t",
        header=None,
    )
    R_df2 = ld2.values

    ro.r("set.seed(123)")
    fit = susieR.susie_rss(
        bhat=locus_df["BETA"].values.reshape((len(R_df), 1)),
        shat=locus_df["SE"].values.reshape((len(R_df), 1)),
        R=R_df,
        L=L,
        n=n,
        check_prior=True,
    )

    locus_df["cs"] = 0
    n_cs = len(
        susieR.susie_get_cs(fit, coverage=0.95, min_abs_corr=min_abs_corr, Xcorr=R_df)[
            0
        ]
    )
    for i in range(n_cs):
        cs_index = susieR.susie_get_cs(
            fit, coverage=0.95, min_abs_corr=min_abs_corr, Xcorr=R_df
        )[0][i]
        locus_df.loc[np.array(cs_index) - 1, "cs"] = i + 1
    locus_df["pip"] = np.array(susieR.susie_get_pip(fit))

    locus_df.to_csv(os.path.join(RESULTS_DIR, f"locus_{rsid}_{flank_size}_mapped.tsv"))
    return locus_df

def plot_mapped_locus(locus, mapped_locus_df: pd.DataFrame, ld):
    fig ,axes = plt.subplots(nrows=2,sharex=True,figsize=(15,7),height_ratios=(4,1))
    df = mapped_locus_df.copy(deep=True)
    df["MLOG10P"] = -np.log10(df["P"])
    col_to_plot = "MLOG10P"
    p=axes[0].scatter(df["POS"],df[col_to_plot],c=ld[df["P"].idxmin()]**2)

    axes[0].scatter(df.loc[df["cs"]>0,"POS"],df.loc[df["cs"]>0,col_to_plot],
            marker='o',s=40,c="None",edgecolors='black',label="Variants in credible sets")

    plt.colorbar( p , label="Rsq with the lead variant")
    axes[0].set_xlabel("position")
    axes[0].set_xlim((locus[1], locus[2]))
    axes[0].set_ylabel(col_to_plot)
    axes[0].legend()

    p=axes[1].scatter(df["POS"],df["pip"],c=ld[df["P"].idxmin()]**2)

    axes[1].scatter(df.loc[df["cs"]>0,"POS"],df.loc[df["cs"]>0,"pip"],
            marker='o',s=40,c="None",edgecolors='black',label="Variants in credible sets")

    plt.colorbar( p , label="Rsq with the lead variant")
    axes[1].set_xlabel("position")
    axes[1].set_xlim((locus[1], locus[2]))
    axes[1].set_ylabel("PIP")
    axes[1].legend()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fine Mapping")
    parser.add_argument("sumstats_name", type=str, help="Name of the sumstats file")
    parser.add_argument("snp_file", type=str, help="Name of file containing snps of interest")
    args = parser.parse_args()

    sumstats = load_sumstats(args.sumstats_name)
    harmonize(sumstats)

    loci = get_loci(args.snp_file, chr_col="chr", pos_col="pos", rsid_col="snp")
    output_loci_files(sumstats, loci)

    print("Calculating LD and fine-mapping loci")
    for locus in tqdm(loci):
        calculate_LD(locus)
        fine_map_locus(locus)