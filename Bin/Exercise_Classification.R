# ------------------------------------------------------ #
# Classify genes by how their MTJ-vs-CK8 translational
# signal changes between sedentary and exercise.
#
# Expects:
#   de_sedentary = Riboseq_Differential_Expression()$res
#     (contrast: MTJ_HA.Rpl22 vs CK8_HA.Rpl22)
#   de_exercise  = Differential_Expression_Exercise()$res
#     (contrast: MTJ_HA vs CK8_HA)
#   de_combined  = Riboseq_Differential_Expression_Combination()$res
#     (contrast: MTJ_Ex vs MTJ_Unt, CK8_Ex vs CK8_Unt)
# ------------------------------------------------------ #

classify_exercise_response = function(
    de_sedentary,
    de_exercise,
    de_combined,
    contrast_sed      = "MTJ_HA.Rpl22_CK8_HA.Rpl22",
    contrast_exe      = "MTJ_HA_CK8_HA",
    contrast_mtj_swap = "MTJ_Ex_MTJ_Unt",
    contrast_ck8_swap = "CK8_Ex_CK8_Unt",
    join_by           = "gene_id",
    lfc_cut           = 0.75,
    padj_cut          = 0.05
){
    suppressPackageStartupMessages({
        library(dplyr)
    })

    pick = function(df, contrast){
        df %>% dplyr::select(
            all_of(join_by),
            !!paste0("log2FoldChange_", contrast) := !!sym(paste0("log2FoldChange.", contrast)),
            !!paste0("padj_",           contrast) := !!sym(paste0("padj.",           contrast))
        )
    }

    sed = pick(de_sedentary, contrast_sed)
    exe = pick(de_exercise,  contrast_exe)
    mtj = pick(de_combined,  contrast_mtj_swap)
    ck8 = pick(de_combined,  contrast_ck8_swap)

    joined = sed %>%
        inner_join(exe, by = join_by) %>%
        inner_join(mtj, by = join_by) %>%
        inner_join(ck8, by = join_by)

    lfc_sed = paste0("log2FoldChange_", contrast_sed)
    lfc_exe = paste0("log2FoldChange_", contrast_exe)
    p_sed   = paste0("padj_",           contrast_sed)
    p_exe   = paste0("padj_",           contrast_exe)

    is_sig = function(lfc, p) abs(lfc) > lfc_cut & !is.na(p) & p < padj_cut

    joined %>%
        mutate(
            sig_sed = is_sig(.data[[lfc_sed]], .data[[p_sed]]),
            sig_exe = is_sig(.data[[lfc_exe]], .data[[p_exe]]),
            dir_sed = sign(.data[[lfc_sed]]),
            dir_exe = sign(.data[[lfc_exe]]),
            response_class = case_when(
                !sig_sed &  sig_exe & dir_exe > 0 ~ "translationally_up_EXE_only",
                !sig_sed &  sig_exe & dir_exe < 0 ~ "translationally_down_EXE_only",
                 sig_sed & !sig_exe & dir_sed > 0 ~ "up_sedentary_only",
                 sig_sed & !sig_exe & dir_sed < 0 ~ "down_sedentary_only",
                 sig_sed &  sig_exe & dir_sed ==  dir_exe & dir_sed > 0 ~ "up_in_both",
                 sig_sed &  sig_exe & dir_sed ==  dir_exe & dir_sed < 0 ~ "down_in_both",
                 sig_sed &  sig_exe & dir_sed != dir_exe              ~ "reversed",
                TRUE                                                    ~ "unchanged"
            )
        )
}
