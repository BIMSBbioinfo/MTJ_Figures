# These DE functions build their cdata with stringr (str_detect/str_replace),
# dplyr and tidyr. Attach them when the file is sourced so the functions don't
# depend on whatever the calling report happens to have on the search path.
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(stringr)
})

Riboseq_Differential_Expression = function(
    inpath  = "/data/local/Projects/MKim_MTJ/Results/Processing/pigx-riboseq",
    outpath = Path(file.path(.scripts.path, "Results"),"differential_expression"),
    path_gtf = "/data/local/vfranke/Annotation/mm39/EnsemblGenes/Mus_musculus.GRCm39.113.chr.gtf",
    contlist = list(
        MTJ_HA.Rpl22_CK8_HA.Rpl22 = c("Factor","MTJ_HA.Rpl22","CK8_HA.Rpl22")
    )
){


    path_salmon = file.path(inpath, "salmon_output")
    cdata = list.files(path_salmon) %>%
        data.frame(sample_name = .) %>%
        tidyr::separate(sample_name, into=c("cell","type","replicate"), sep="_", remove=FALSE) %>%
        filter(!str_detect(type, "input")) %>%
        filter(!str_detect(type, "GFP")) %>%
        mutate(condition = str_replace(sample_name, "_br.+","")) %>%
        mutate(Factor = condition) %>%
        mutate(sample = sample_name)



    source.lib("Read_Annotation.R")
    annot = read_Gene_Annotation(path_gtf)


    source.lib("Bacon.R")

    res = summarize_PigX_RNAseq(
    path_rnaseq = file.path(path_salmon, cdata$sample_name),
    contlist    = contlist,
    cdata       = cdata,
    path_gtf    = path_gtf,
    design      = "Factor",
    input_type  = "salmon"
    )

    res$res = res$res %>%
        mutate(diff.MTJ_HA.Rpl22_CK8_HA.Rpl22 = case_when(
            log2FoldChange.MTJ_HA.Rpl22_CK8_HA.Rpl22 > 0.75 & padj.MTJ_HA.Rpl22_CK8_HA.Rpl22 < 0.05 ~ "Up",
            log2FoldChange.MTJ_HA.Rpl22_CK8_HA.Rpl22 < -0.75 & padj.MTJ_HA.Rpl22_CK8_HA.Rpl22 < 0.05 ~ "Down",
            TRUE ~ "No"
        ))
    res$dat = res$dat %>%
        mutate(diff.MTJ_HA.Rpl22_CK8_HA.Rpl22 = case_when(
            log2FoldChange.MTJ_HA.Rpl22_CK8_HA.Rpl22 > 0.75 & padj.MTJ_HA.Rpl22_CK8_HA.Rpl22 < 0.05 ~ "Up",
            log2FoldChange.MTJ_HA.Rpl22_CK8_HA.Rpl22 < -0.75 & padj.MTJ_HA.Rpl22_CK8_HA.Rpl22 < 0.05 ~ "Down",
            TRUE ~ "No"
        ))
    return(res)
}

# ------------------------------------------------------------------- #
Differential_Expression_Exercise = function(
    inpath  = "/data/local/Projects/MKim_MTJ/Results/Riboseq_Exercise/pigx-rnaseq",
    outpath = Path(file.path(.scripts.path, "Results"),"differential_expression_exercise"),
    path_gtf = "/data/local/vfranke/Annotation/mm39/EnsemblGenes/Mus_musculus.GRCm39.113.chr.gtf",
    contlist = list(
        MTJ_HA_CK8_HA = c("Factor","MTJ_HA","CK8_HA"),
        MTJ_HA_MTJ_Input = c("Factor","MTJ_HA","MTJ_Input")
    )
){


    path_salmon = file.path(inpath, "salmon_output")
    cdata = list.files(path_salmon) %>%
        data.frame(sample_name = .) %>%
        tidyr::separate(sample_name, into=c("cell","type","replicate"), sep="_", remove=FALSE) %>%
        filter(!str_detect(type, "GFP")) %>%
        mutate(condition = str_replace(sample_name, "_br.+","")) %>%
        mutate(Factor = condition) %>%
        mutate(sample = sample_name)



    source.lib("Read_Annotation.R")
    annot = read_Gene_Annotation(path_gtf)


    source.lib("Bacon.R")

    res = summarize_PigX_RNAseq(
        path_rnaseq = file.path(path_salmon, cdata$sample_name),
        contlist    = contlist,
        cdata       = cdata,
        path_gtf    = path_gtf,
        design      = "Factor",
        input_type  = "salmon"
    )

    res$res = res$res %>%
        mutate(diff.MTJ_HA_CK8_HA = case_when(
            log2FoldChange.MTJ_HA_CK8_HA > 0.75 & padj.MTJ_HA_CK8_HA < 0.05 ~ "Up",
            log2FoldChange.MTJ_HA_CK8_HA < -0.75 & padj.MTJ_HA_CK8_HA < 0.05 ~ "Down",
            TRUE ~ "No"
        ))
    res$dat = res$dat %>%
        mutate(diff.MTJ_HA_MTJ_Input = case_when(
            log2FoldChange.MTJ_HA_MTJ_Input > 0.75 & padj.MTJ_HA_MTJ_Input < 0.05 ~ "Up",
            log2FoldChange.MTJ_HA_MTJ_Input < -0.75 & padj.MTJ_HA_MTJ_Input < 0.05 ~ "Down",
            TRUE ~ "No"
        ))
    return(res)
}


# ==================================================================== #
Riboseq_Differential_Expression_Combination = function(
    inpaths  = c("/data/local/Projects/MKim_MTJ/Results/Processing/pigx-riboseq","/data/local/Projects/MKim_MTJ/Results/Riboseq_Exercise/pigx-rnaseq"),
    outpath = Path(file.path(.scripts.path, "Results"),"differential_expression"),
    path_gtf = "/data/local/vfranke/Annotation/mm39/EnsemblGenes/Mus_musculus.GRCm39.113.chr.gtf",
    contlist = list(
        MTJ_HA.Unt_CK8_HA.Unt = c("Factor","MTJ_HA.Rpl22","CK8_HA.Rpl22"),
        MTJ_Ex_CK8_Ex         = c("Factor","MTJ_HA","CK8_HA"),
        MTJ_Ex_MTJ_Input      = c("Factor","MTJ_HA","MTJ_Input"),
        MTJ_Ex_MTJ_Unt        = c("Factor","MTJ_HA","MTJ_HA.Rpl22"),
        CK8_Ex_CK8_Unt        = c("Factor","CK8_HA","CK8_HA.Rpl22")
    )
){


    ltab = lapply(inpaths, function(inpath){
        path_salmon = file.path(inpath, "salmon_output")
        if(str_detect(inpath, "Exercise")){
            
            cdata = list.files(path_salmon) %>%
                data.frame(sample_name = .) %>%
                tidyr::separate(sample_name, into=c("cell","type","replicate"), sep="_", remove=FALSE) %>%
                filter(!str_detect(type, "GFP")) %>%
                mutate(condition = str_replace(sample_name, "_br.+","")) %>%
                mutate(Factor = condition) %>%
                mutate(sample = sample_name) %>%
                mutate(experiment = "Untreated")  %>%
                mutate(file_name = file.path(path_salmon, sample_name))

        }else{
                        
            cdata = list.files(path_salmon) %>%
                data.frame(sample_name = .) %>%
                tidyr::separate(sample_name, into=c("cell","type","replicate"), sep="_", remove=FALSE) %>%
                filter(!str_detect(type, "input")) %>%
                filter(!str_detect(type, "GFP")) %>%
                mutate(condition = str_replace(sample_name, "_br.+","")) %>%
                mutate(Factor = condition) %>%
                mutate(sample = sample_name) %>%
                mutate(experiment = "Exercise") %>%
                mutate(file_name = file.path(path_salmon, sample_name))

        }
    })
    cdata = do.call(rbind, ltab)


    source.lib("Read_Annotation.R")
    annot = read_Gene_Annotation(path_gtf)


    source.lib("Bacon.R")

    res = summarize_PigX_RNAseq(
        path_rnaseq = cdata$file_name,
        contlist    = contlist,
        cdata       = cdata %>% dplyr::select(-file_name),
        path_gtf    = path_gtf,
        design      = "Factor",
        input_type  = "salmon"
    )

    return(res)
}


# ------------------------------------------------------------------- #
RNAseq_Differential_Expression = function(
    inpath  = "/data/local/Projects/MKim_MTJ/Results/Processing/RNAseq/pigx-rnaseq/",
    path_gtf = "/data/local/vfranke/Annotation/mm39/EnsemblGenes/Mus_musculus.GRCm39.113.chr.gtf",
    contlist = list(
        Tigd4_CK8 = c("Factor","Tigd4","CK8")
    ),
    lfc = 0.75,
    samples_to_remove = "Tigd4_br2"
){


    path_salmon = file.path(inpath, "salmon_output")
    cdata = list.files(path_salmon) %>%
        data.frame(sample_name = .) %>%
        tidyr::separate(sample_name, into=c("cell","replicate"), sep="_", remove=FALSE) %>%
        mutate(condition = str_replace(sample_name, "_br.+","")) %>%
        mutate(Factor = condition) %>%
        mutate(sample = sample_name) %>%
        filter(!sample %in% samples_to_remove)

    source.lib("Read_Annotation.R")
    annot = read_Gene_Annotation(path_gtf)

    source.lib("Bacon.R")

    res = summarize_PigX_RNAseq(
        path_rnaseq = file.path(path_salmon, cdata$sample_name),
        contlist    = contlist,
        cdata       = cdata,
        path_gtf    = path_gtf,
        design      = "Factor",
        input_type  = "salmon",
        lfc         = lfc
    )

    return(res)
}
