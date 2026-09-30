# ------------------------------------------------------ #

read_POSTAR3_mouse = function(
    path = "/home/vfranke/Projects/MKim_MTJ/Documentation/POSTAR3/mouse.txt.gz",
    chain_path = "/data/local/vfranke/Base/UCSC_LiftOver/mm10ToMm39.over.chain"
){
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(rtracklayer)
        library(GenomicRanges)
    })

    # Initialize AnnotationHub
    message("chain ...")
    chain = rtracklayer::import.chain(chain_path)

    message("tab ...")
    gtab = fread(path, header = FALSE, sep = "\t") %>%
        magrittr::set_colnames(c("seqnames","start","end","db","strand","gene_name","experiment","cell","gseid","score")) %>%
        as.data.frame() %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)

    message("liftover ...")
    gtab_mm39 = liftOver(gtab, chain)
    gtab_mm39 = unlist(gtab_mm39[elementNROWS(gtab_mm39) == 1])

    return(gtab_mm39)
}


# ------------------------------------------------------ #
read_Tigd3_RIP = function(
    path = "/home/vfranke/Projects/MKim_MTJ/Documentation/Tigd4_RIP/Tigd4 RIP-Seq all genes.xlsx"
){
    suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(readxl)
    })
    tab = readxl::read_xlsx(path) %>%
        magrittr::set_colnames(tolower(colnames(.))) %>%
        magrittr::set_colnames(str_replace_all(colnames(.)," ","_")) %>%
        dplyr::select(-contains("median")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"_\\(normalize.+","_ln")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\(","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\)","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"-","_"))  %>%
        mutate(log2_fc_ip_vs_input_t4 = as.numeric(log2_fc_ip_vs_input_t4)) %>%
        mutate(adjusted_p_value_ip_vs_input_t4 = as.numeric(adjusted_p_value_ip_vs_input_t4)) 
    return(tab)

}

# ------------------------------------------------------ #
read_RNAseq_DifferentialExpression = function(
    path = "/home/vfranke/Projects/MKim_MTJ/Documentation/RNAseq/CK8 vs Tigd4 only FC and p-values.xlsx"
){
    tab = readxl::read_xlsx(path) %>%
        magrittr::set_colnames(tolower(colnames(.))) %>%
        magrittr::set_colnames(str_replace_all(colnames(.)," ","_")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"-","_"))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\(",""))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\)",""))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"/","")) %>%
        dplyr::select(gene_name, contains("wo_105")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"_wo_105",""))
    return(tab)
}



# ------------------------------------------------------ #
read_Mass_Spectrometry_MTJ = function(
    path = "/home/vfranke/Projects/MKim_MTJ/Documentation/Mass_Spectrometry/250922_MK_RatioPerseus.xlsx"
    
){
    tab = readxl::read_xlsx(path, sheet=2) %>%
        magrittr::set_colnames(tolower(colnames(.))) %>%
        magrittr::set_colnames(str_replace_all(colnames(.)," ","_")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"#","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),":","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\[","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\]","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"-","_"))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\(",""))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\)",""))  %>%
        magrittr::set_colnames(str_replace(colnames(.),"^_",""))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"/","_"))%>%
        dplyr::select(gene_name, accession, logp_value, difference_log2, ratio_ck8_mtj, peptides, unique_peptides, contains("xic"))
    
    return(tab)
}

read_MS_Raw_MTJ = function(
    path = "/home/vfranke/Projects/MKim_MTJ/Documentation/Mass_Spectrometry/250922_MK_Results.xlsx"
){
   
    tab = readxl::read_xlsx(path, sheet=1) %>%
        magrittr::set_colnames(tolower(colnames(.))) %>%
        magrittr::set_colnames(str_replace_all(colnames(.)," ","_")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"#","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),":","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\[","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\]","")) %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"-","_"))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\(",""))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"\\)",""))  %>%
        magrittr::set_colnames(str_replace(colnames(.),"^_",""))  %>%
        magrittr::set_colnames(str_replace_all(colnames(.),"/","_")) %>%
        dplyr::rename(
            CK8_r1 = abundance_f1_sample,
            CK8_r2 = abundance_f2_sample,
            CK8_r3 = abundance_f3_sample,
            MTJ_r1 = abundance_f4_sample,
            MTJ_r2 = abundance_f5_sample,
            MTJ_r3 = abundance_f6_sample,
            input_r1 = abundance_f7_sample,
            input_r2 = abundance_f8_sample,
            input_r3 = abundance_f9_sample
        ) %>%
        mutate(
            mean_CK8 = rowMeans(dplyr::select(., CK8_r1, CK8_r2, CK8_r3), na.rm=TRUE),
            mean_MTJ = rowMeans(dplyr::select(., MTJ_r1, MTJ_r2, MTJ_r3), na.rm=TRUE),
            mean_input = rowMeans(dplyr::select(., input_r1, input_r2, input_r3), na.rm=TRUE),
        ) %>%
        mutate(ratio_mtj_ck8 = mean_MTJ / mean_CK8) %>%
        mutate(lfc_mtj_ck8   = log2(ratio_mtj_ck8)) %>%
        mutate(ratio_mtj_input = mean_MTJ / mean_input) %>%
        mutate(lfc_mtj_input   = log2(ratio_mtj_input)) %>%
        arrange(lfc_mtj_ck8) %>%
        filter(!contaminant) %>%
        filter(unique_peptides >= 2)  %>%
        dplyr::select(accession, contains("r\\d"),contains("mean"), contains("ratio"), contains("lfc"))  
    return(tab)

}