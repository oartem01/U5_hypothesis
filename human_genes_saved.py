#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 11:31:24 2020

@author: olga
"""

def exon_files():
    """List of exon sequences for 132 human genes downloaded from ensembl and
    saved locally."""

    exon_files=['001_DMD_203_exons.txt','002_F5_202_exons.txt',
    '003_CFTR_201_exons.txt','004_F8_202_exons.txt',
    '005_F9_201_exons.txt','006_VWF_201_exons.txt',
    '007_SMN1_202_exons.txt','008_PAH_215_exons.txt',
    '009_APOE_201_exons.txt','010_BRCA1_210_exons.txt',
    '011_HPRT1_201_exons.txt','012_HBB_206_exons.txt',
    '013_FAS_221_exons.txt','014_TP53_223_exons.txt',
    '015_TNF_208_exons.txt','016_UBC_201_exons.txt',
    '017_EGFR_201_exons.txt','018_VEGFA_205_exons.txt',
    '019_IL6_204_exons.txt','020_TGFB1_201_exons.txt',
    '021_ESR1_207_exons.txt','022_HLA_DRB1_201_exons.txt',
    '023_NFKB1_201_exons.txt','024_IL10_202_exons.txt',
    '025_AKT1_208_exons.txt','026_CD4_201_exons.txt',
    '027_GRB2_203_exons.txt','028_ELP1_201_exons.txt',
    '029_NOVA1_206_exons.txt','030_TOP1_201_exons.txt',
    '031_INS_203_exons.txt','032_FHL5_201_exons.txt',
    '033_TPM1_206_exons.txt','034_DAG1_222_exons.txt',
    '035_POLR2A_201_exons.txt','036_PRPF8_201_exons.txt',
    '037_DHX8_201_exons.txt','038_SNRNP200_201_exons.txt',
    '039_MYC_207_exons.txt','040_TBP_202_exons.txt',
    '041_LAMB1_201_exons.txt','042_NOTCH1_205_exons.txt',
    '043_SGCA_201_exons.txt','044_TUBA1A_202_exons.txt',
    '045_MYO7A_202_exons.txt','046_NUP155_201_exons.txt',
    '047_PTBP1_203_exons.txt','048_U2AF2_202_exons.txt',
    '049_GAPDH_201_exons.txt','050_TK1_201_exons.txt',
    '051_ATM_201_exons.txt','052_TERT_201_exons.txt',
    '053_ADAMTS13_206_exons.txt','054_NEU1_206_exons.txt',
    '055_GLB1_201_exons.txt','056_CTSA_204_exons.txt',
    '057_ACE_202_exons.txt','058_ERN1_201_exons.txt',
    '059_BMP2_201_exons.txt','060_CDK6_201_exons.txt',
    '061_CCND1_201_exons.txt','062_EZH1_202_exons.txt',
    '063_HAT1_201_exons.txt','064_ITGB3_201_exons.txt',
    '065_IFNG_201_exons.txt','066_ERBB2_219_exons.txt',
    '067_APP_201_exons.txt','068_EGF_201_exons.txt',
    '069_CTNNB1_201_exons.txt','070_IGF1_203_exons.txt',
    '071_MAPK1_201_exons.txt','072_IGHG1_202_exons.txt',
    '073_TFRC_201_exons.txt','074_CRP_201_exons.txt',
    '075_MET_201_exons.txt','076_RGS6_209_exons.txt',
    '077_BCL2_202_exons.txt','078_XRCC6_202_exons.txt',
    '079_RAD52_202_exons.txt','080_MSTN_201_exons.txt',
    '081_PNMT_201_exons.txt','082_CELF1_214_exons.txt',
    '083_BLM_201_exons.txt','084_PSEN1_201_exons.txt',
    '085_AMY1A_201_exons.txt','086_SPSB1_201_exons.txt',
    '087_G6PD_202_exons.txt','088_USB1_201_exons.txt',
    '089_ERCC6_201_exons.txt','090_ATP2B1_201_exons.txt',
    '091_ETV4_201_exons.txt','092_HNF1A_201_exons.txt',
    '093_NAT2_201_exons.txt','094_LPL_207_exons.txt',
    '095_FADS1_201_exons.txt','096_GDF5_201_exons.txt',
    '097_ACADM_202_exons.txt','098_GATA1_202_exons.txt',
    '099_GSTM1_201_exons.txt','100_CDH1_201_exons.txt',
    '101_FMR1_205_exons.txt','102_POLD1_205_exons.txt',
    '103_RPA1_201_exons.txt','104_WRN_201_exons.txt',
    '105_SP1_201_exons.txt','106_TWIST1_201_exons.txt',
    '107_TAF1_203_exons.txt','108_POLB_201_exons.txt',
    '109_TERF2_201_exons.txt','110_RFC1_201_exons.txt',
    '111_MAPT_204_exons.txt','112_ACTA2_201_exons.txt',
    '113_PCGF2_207_exons.txt','114_SSBP3_204_exons.txt',
    '115_CHEK1_213_exons.txt','116_BRIP1_201_exons.txt',
    '117_NFATC2_201_exons.txt','118_PRKDC_201_exons.txt',
    '119_OPTN_202_exons.txt','120_ZBP1_201_exons.txt',
    '121_LDHA_205_exons.txt','122_RBPJ_205_exons.txt',
    '123_FKTN_201_exons.txt','124_SFPQ_201_exons.txt',
    '125_MAD1L1_201_exons.txt','126_RAD51_201_exons.txt',
    '127_POLA1_201_exons.txt','128_PCNA_201_exons.txt',
    '129_PRIM2_201_exons.txt','130_ITPR1_203_exons.txt',
    '131_EFNA5_201_exons.txt','132_OPA1_205_exons.txt']
    print(len(exon_files))
    return(exon_files)
    
def intron_files():
    """List of intron sequences for 132 human genes downloaded from ensembl and
    saved locally.""" 
    intron_files=['001_DMD_203_introns.txt','002_F5_202_introns.txt',
    '003_CFTR_201_introns.txt','004_F8_202_introns.txt',
    '005_F9_201_introns.txt','006_VWF_201_introns.txt',
    '007_SMN1_202_introns.txt','008_PAH_215_introns.txt',
    '009_APOE_201_introns.txt','010_BRCA1_210_introns.txt',
    '011_HPRT1_201_introns.txt','012_HBB_206_introns.txt',
    '013_FAS_221_introns.txt','014_TP53_223_introns.txt',
    '015_TNF_208_introns.txt','016_UBC_201_introns.txt',
    '017_EGFR_201_introns.txt','018_VEGFA_205_introns.txt',
    '019_IL6_204_introns.txt','020_TGFB1_201_introns.txt',
    '021_ESR1_207_introns.txt','022_HLA_DRB1_201_introns.txt',
    '023_NFKB1_201_introns.txt','024_IL10_202_introns.txt',
    '025_AKT1_208_introns.txt','026_CD4_introns.txt',
    '027_GRB2_203_introns.txt','028_ELP1_201_introns.txt',
    '029_NOVA1_206_introns.txt','030_TOP1_201_introns.txt',
    '031_INS_203_introns.txt','032_FHL5_201_introns.txt',
    '033_TPM1_206_introns.txt','034_DAG1_222_introns.txt',
    '035_POLR2A_201_introns.txt','036_PRPF8_201_introns.txt',
    '037_DHX8_201_introns.txt','038_SNRNP200_201_introns.txt',
    '039_MYC_207_introns.txt','040_TBP_202_introns.txt',
    '041_LAMB1_201_introns.txt','042_NOTCH1_205_introns.txt',
    '043_SGCA_201_introns.txt','044_TUBA1A_202_introns.txt',
    '045_MYO7A_202_introns.txt','046_NUP155_201_introns.txt',
    '047_PTBP1_203_introns.txt','048_U2AF2_202_introns.txt',
    '049_GAPDH_201_introns.txt','050_TK1_201_introns.txt',
    '051_ATM_201_introns.txt','052_TERT_201_introns.txt',
    '053_ADAMTS13_206_introns.txt','054_NEU1_206_introns.txt',
    '055_GLB1_201_introns.txt','056_CTSA_204_introns.txt',
    '057_ACE_202_introns.txt','058_ERN1_201_introns.txt',
    '059_BMP2_201_introns.txt','060_CDK6_201_introns.txt',
    '061_CCND1_201_introns.txt','062_EZH1_202_introns.txt',
    '063_HAT1_201_introns.txt','064_ITGB3_201_intons.txt',
    '065_IFNG_201_introns.txt','066_ERBB2_219_introns.txt',
    '067_APP_201_introns.txt','068_EGF_201_introns.txt',
    '069_CTNNB1_201_introns.txt','070_IGF1_203_introns.txt',
    '071_MAPK1_201_introns.txt','072_IGHG1_202_introns.txt',
    '073_TFRC_201_introns.txt','074_CRP_201_introns.txt',
    '075_MET_201_introns.txt','076_RGS6_209_introns.txt',
    '077_BCL2_202_introns.txt','078_XRCC6_202_introns.txt',
    '079_RAD52_202_introns.txt','080_MSTN_201_introns.txt',
    '081_PNMT_201_introns.txt','082_CELF1_214_introns.txt',
    '083_BLM_201_introns.txt','084_PSEN1_201_introns.txt',
    '085_AMY1A_201_introns.txt','086_SPSB1_201_introns.txt',
    '087_G6PD_202_introns.txt','088_USB1_201_introns.txt',
    '089_ERCC6_201_introns.txt','090_ATP2B1_201_introns.txt',
    '091_ETV4_201_introns.txt','092_HNF1A_201_introns.txt',
    '093_NAT2_201_introns.txt','094_LPL_207_introns.txt',
    '095_FADS1_201_introns.txt','096_GDF5_201_introns.txt',
    '097_ACADM_202_introns.txt','098_GATA1_202_introns.txt',
    '099_GSTM1_201_introns.txt','100_CDH1_201_introns.txt',
    '101_FMR1_205_introns.txt','102_POLD1_205_introns.txt',
    '103_RPA1_201_introns.txt','104_WRN_201_introns.txt',
    '105_SP1_201_introns.txt','106_TWIST1_201_introns.txt',
    '107_TAF1_203_introns.txt','108_POLB_201_introns.txt',
    '109_TERF2_201_introns.txt','110_RFC1_201_introns.txt',
    '111_MAPT_204_introns.txt','112_ACTA2_201_introns.txt',
    '113_PCGF2_207_introns.txt','114_SSBP3_204_introns.txt',
    '115_CHEK1_213_introns.txt','116_BRIP1_201_introns.txt',
    '117_NFATC2_201_introns.txt','118_PRKDC_201_introns.txt',
    '119_OPTN_202_introns.txt','120_ZBP1_201_introns.txt',
    '121_LDHA_205_introns.txt','122_RBPJ_205_introns.txt',
    '123_FKTN_201_introns.txt','124_SFPQ_201_introns.txt',
    '125_MAD1L1_201_introns.txt','126_RAD51_201_introns.txt',
    '127_POLA1_201_introns.txt','128_PCNA_201_introns.txt',
    '129_PRIM2_201_introns.txt','130_ITPR1_203_introns.txt',
    '131_EFNA5_201_introns.txt','132_OPA1_205_introns.txt']
    print(len(intron_files))
    return (intron_files)