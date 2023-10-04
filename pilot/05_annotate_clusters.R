library(tidyverse)
source("dataset.R")


main = function() {
    annotate_whole()
    annotate_glu()
    annotate_gaba()
    annotate_subglu()
}

annotate_whole = function(output_dir = "analysis/whole") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster_name,
        1, "GLU_1 (L2/3 IT)",
        2, "SubCTX_1",
        3, "GABA_1 (?)",
        4, "SubCTX_2",
        5, "GLU_2 (L5 PT)",
        6, "GLU_3 (?)",
        7, "GLU_4 (Piriform)",
        8, "GABA_2 (?)",
        9, "GABA_3 (CTX)",
        10, "GLU_5 (L4/5 IT)",
        11, "GLU_6 (NP)",
        12, "GABA_4 (SubCTX)",
        13, "GLU_7 (hippocampus)",
        14, "GLU_8 (L5 IT)",
        15, "GLU_9 (L6 CT+b)",
        16, "GLU_10 (?)",
        17, "GLU_11 (Car3)",
        18, "SubCTX_3",
        19, "GABA_5 (CTX)",
        20, "SubCTX_4",
        21, "SubCTX_5",
        22, "SubCTX_6",
        23, "GABA_6 (subCTX)",
        24, "SubCTX_7"
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_glu = function(output_dir = "analysis/glu") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster_name,
        1, "IT_12 (TH PG/TRN/MD?+OB)",
        2, "IT_1 (ENT)",
        3, "CT",
        4, "L6b",
        5, "NP",
        6, "PT",
        7, "IT_2 (HIP Sub.)",
        8, "IT_3 (PIRI UL+ENT+AMY CoA/LA)",
        9, "IT_4 (HIP CA+AMY BLA?+PIRI DL+TTd+mPFCa)",
        10, "IT_5 (PPP+RSP UL)",
        11, "IT_6 (L6 IT + HIP DG)",
        12, "IT_7 (L2/3 IT)",
        13, "IT_8 (DL medial, RSP DL?)",
        14, "IT_9 (L5 IT+AMY CoA/BLA?)",
        15, "Car3",
        16, "IT_10 (Low quality)",
        17, "Car3-like (TH Habenula?)",
        18, "IT_11 (L4/5 IT)"
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_gaba = function(output_dir = "analysis/gaba") {
    cluster_annotation = tribble(
        ~cluster_id, ~cluster, ~subclass,
        1, "Pvalb","Pvalb",
        2, "Sst","Sst",
        3, "GABA_1","GABA_1",
        4, "Vip/Sncg","Vip/Sncg",
        5, "GABA_2","GABA_2",
        6, "Lamp5","Lamp5",
        7, "GABA_3","GABA_3",
        8, "GABA_4","GABA_4",
        9, "Sst Chodl-like","Sst Chodl-like",
        10, "GABA_5","GABA_5",
        11, "Meis2-like","Meis2-like"
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

annotate_subglu = function(output_dir = "analysis/subglu") {
    subclass_annotation = read_csv("analysis/glu/cluster_annotation.csv")
    cluster_annotation = tribble(
        ~cluster_id, ~cluster, ~subclass,~notes,
        "Car3_1","Car3 CTX P-L","Car3","CTX P-L",
        "Car3_2","Car3 CLA/EPd","Car3","CLA/EPd",
        "Car3_3","Car3 CTX UL","Car3","CTX UL",
        "Car3_4","Car3 CTX DL","Car3","CTX DL",
        "Car3_5","Car3 EPd/CLA","Car3","EPd/CLA",
        "Car3-like_1","MH D", "TH MH","D",
        "Car3-like_2","MH V-M", "TH MH","V-M",
        "Car3-like_3","MH V-L", "TH MH","V-L",
        "CT_1","CT CTX UL","CT","CTX UL",
        "CT_2","CT CTX P/RSP","CT","CTX P/RSP",
        "CT_3","CT RSP","CT","RSP",
        "CT_4","CT CTX A","CT","CTX A",
        "CT_5","CT Unclear/LQ","CT","Unclear/LQ?",
        "CT_6","CT CTX A/ACA","CT","CTX A/ACA",
        "IT_1_1","ECT P_1?","ENTl L2/3","ECT P?",
        "IT_1_2","ENTl L2/3 D","ENTl L2/3","D",
        "IT_1_3","ENTl L2/3 V","ENTl L2/3","V",
        "IT_1_4","ECT P_2?","ENTl L2/3","ECT P?",
        "IT_10_1","IT LQ","Low Quality","LQ",
        "IT_10_2","IT LQ","Low Quality","LQ",
        "IT_10_3","IT LQ","Low Quality","LQ",
        "IT_10_4","IT LQ","Low Quality","LQ",
        "IT_10_5","IT LQ","Low Quality","LQ",
        "IT_11_1","L4/5 IT M-L","L4/5 IT", "CTX M-L",
        "IT_11_2","L4/5 IT UL","L4/5 IT", "CTX UL",
        "IT_11_3","L4/5 IT ML-A","L4/5 IT", "CTX ML-A",
        "IT_11_4","L4/5 IT ML-P","L4/5 IT", "CTX ML-P",
        "IT_11_5","L4/5 IT DL","L4/5 IT", "CTX DL",
        "IT_11_6","L4/5 IT P-L/LA","L4/5 IT", "CTX P-L/LA",
        "IT_11_7","L4/5 IT M","L4/5 IT", "CTX M",
        "IT_12_1","P V?","Unclear", "Pons?",
        "IT_12_2","P D?","Unclear","Pons?",
        "IT_12_3","TH AD?","TH AD", "AD?",
        "IT_12_4","IT_12_4","Unclear","LQ?",
        "IT_12_5","PIR L3 A?","Unclear", "PIR L3 A?",
        "IT_12_6","IT_12_6","Unclear","LQ?",
        "IT_12_7","MOB UL","MOB UL", "MOB UL",
        "IT_12_8","AON UL","AON UL", "AON UL",
        "IT_2_1","SUBd-m","HPF SUBd", "Molecular layer?",
        "IT_2_2","SUB P?","HPF SUBd", "P",
        "IT_2_3","SUBd-sp","HPF SUBd", "Pyramidal layer?",
        "IT_2_4","SUBd-sr","HPF SUBd", "Stratum radiatum?",
        "IT_3_1","ENTl L1","ENT UL", "L1",
        "IT_3_2","LA D","LA","D",
        "IT_3_3","COAp L2/3","COAp","COAp L2/3/ENTl L1 V/LA V",
        "IT_3_4","COAp L1","COAp", "COAp L1/AIp L1",
        "IT_3_5","PIR L2/3","PIR UL","PIR/ENTl L2/3 V",
        "IT_3_6","PIR L1 DL-A-V","PIR UL","DL-A-V",
        "IT_3_7","PIR L1 DL-P","PIR UL","DL-P",
        "IT_3_8","PIR L1 UL","PIR UL","UL",
        "IT_4_1","BLA/PA/TTd?","BLA/PA/TTd","BLA/PA/TTd?",
        "IT_4_2","SUBv?","HIP SUB","SUBv?",
        "IT_4_3","EPd","EPd","EPd",
        "IT_4_4","CA4","HIP CA","CA3 tip",
        "IT_4_5","CA3","HIP CA","CA3",
        "IT_4_6","AONm","AONm", "AONm/EPv?/ENT?",
        "IT_4_7","CA2","HIP CA","CA2",
        "IT_4_8","CA1","HIP CA","CA1",
        "IT_5_1","ENT?","PPP","ENT?",
        "IT_5_2","PARA","PPP","PARA",
        "IT_5_3","PARA?","PPP","PARA?, slice 3",
        "IT_5_4","POST","PPP","POST",
        "IT_5_5","PRE","PPP","PRE",
        "IT_5_6","RSPd UL","RSP UL","RSPd UL",
        "IT_5_7","RSPv UL","RSP UL","RSPv UL",
        "IT_6_1","L6 IT DL","L6 IT","CTX DL",
        "IT_6_2","DG","HIP DG","DG",
        "IT_6_3","L6 IT ML-P","L6 IT","CTX ML-P",
        "IT_6_4","L6 IT ML-A","L6 IT","CTX ML-A",
        "IT_6_5","L6 IT L","L6 IT","CTX L/BLA?",
        "IT_6_6","L6 IT UL","L6 IT","CTX UL",
        "IT_6_7","PIR/TT A","PIR L6 IT-like","TT/PIR A",
        "IT_6_8","AON DL/PIR A","AON DL","AON DL/PIR A",
        "IT_7_1","L2/3 IT UL-P/ENT","L2/3 IT", "CTX UL-P-L/ENT",
        "IT_7_2","L2/3 IT ML-P","L2/3 IT", "CTX ML-P-M",
        "IT_7_3","L2/3 IT M","L2/3 IT", "CTX M",
        "IT_7_4","L2/3 IT M-L/PIR","L2/3 IT", "CTX M-L/PIR",
        "IT_7_5","L2/3 IT DL","L2/3 IT", "CTX DL",
        "IT_7_6","L2/3 IT ML-A","L2/3 IT", "CTX ML-A",
        "IT_7_7","L2/3 IT Unclear/LQ","L2/3 IT", "Unclear? Very distinct",
        "IT_8_1","PT P RSP/IT-like?","PT", "CTX P, RSP/IT-like?",
        "IT_8_2","RSP L4","RSP DL","Overclustered",
        "IT_8_3","RSP L4","RSP DL","Overclustered",
        "IT_8_4","RSP L4","RSP DL","Overclustered",
        "IT_8_5","RSP L4","RSP DL","Overclustered",
        "IT_8_6","RSP L4","RSP DL","Overclustered",
        "IT_9_1","ENT DL-D/EPv","ENT L5 IT-like","ENT/EPv",
        "IT_9_2","PT AUD","PT", "CTX Unclear AUD?",
        "IT_9_3","PIR L5 IT-like/BLAa","PIR L5 IT-like","PIR L5-IT-like/BLAa",
        "IT_9_4","RSP L5","RSP DL","RSP L5 IT-like",
        "IT_9_5","L5 IT UL","L5 IT","CTX UL",
        "IT_9_6","L5 IT M-L/ENT DL/EPv","L5 IT","CTX/ENT/EPv",
        "IT_9_7","L5 IT M","L5 IT","CTX M",
        "IT_9_8","L5 IT DL","L5 IT","CTX DL",
        "IT_9_9","OLF BA?","OLF","OLF BA?, maybe COA",
        "L6b_1","L6b CTX P","L6b","CTX P",
        "L6b_2","L6b CTX A/EPd?","L6b","CTX A/EPd?",
        "L6b_3","ENT L6b/CT-like","ENT L6b/CT-like", "ENT DL L6b/CT-like",
        "L6b_4","CT CTX A/L-V","CT","CTX A/CTX L-V/CTX around CLA/slice 16",
        "L6b_5","L6b CTX UL","L6b","CTX, slightly UL",
        "L6b_6","PIR L6b-like","PIR L6b-like", "PIR DL",
        "L6b_7","L6b CTX A L-V/EPd","L6b","EPd/CTX A L-V",
        "NP_1","NP RSP","NP","RSP",
        "NP_2","NP PPP","NP","POST3",
        "NP_3","NP CTX L/RSP/SUBv-sp","NP","SUBv pyramidal/RSP/CTX L",
        "NP_4","NP SUBd-sp","NP","SUBd pyramidal",
        "NP_5","NP CTX L5","NP","CTX L5",
        "NP_6","NP CTX L6/M/L","NP","CTX L6/M/L",
        "PT_1","PT CTX P-L/RSPd-agl/TEa","PT","RSP/TEa+ECT+PERI",
        "PT_2","PT CTX P","PT","CTX P",
        "PT_3","PT CTX Unclear","PT","CTX Unclear?",
        "PT_4","PT CTX ML","PT", "CTX ML",
        "PT_5","PT RSPv/PPP","PT","POST2/RSPv",
        "PT_6","PT CTX M/RSPd A","PT","CTX M/RSPd A",
        "PT_7","PT CTX UL","PT","CTX UL",
        "PT_8","PT CTX MOp/s DL","PT", "CTX MOp/s DL",
        "PT_9","PT CTX A-M/L","PT","CTX A-M/L"
    )
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

if (sys.nframe() == 0) {
    main()
}
