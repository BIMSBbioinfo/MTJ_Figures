# ------------------------------------------------------ #
# Render every report under Scripts/Reports/NN_*/NN_*.Rmd
# into Results/<NN_*>/<yymmdd>_<tag>/. Run from the project
# root with R45 (or any R that has the renv library
# activated). Skips the legacy top-level Rmds.
#
#   R45 -e 'source("Scripts/Bin/knit_all.R"); knit_all()'
#
# Parameters
# ----------
# tag       : optional suffix appended to today's date for
#             the output folder. Default "" → folder is just
#             <yymmdd>, which is where the reports write their
#             figures too, so HTML + figures co-locate.
# include   : character vector of report dirnames to knit.
#             NULL = all matching `^[0-9]{2}_[A-Za-z0-9_]+$`.
# skip      : same shape; subtracted from `include`.
# stop_on_error : if FALSE, log the failure and continue.
# ------------------------------------------------------ #

knit_all = function(
    reports_dir   = file.path(.scripts.path, "Scripts", "Reports"),
    results_root  = file.path(.scripts.path, "Results"),
    tag           = "",
    include       = NULL,
    skip          = c(),
    stop_on_error = FALSE
){
    suppressPackageStartupMessages({
        library(rmarkdown)
        library(dplyr)
    })

    pat = "^[0-9]{2}_[A-Za-z0-9_]+$"
    candidates = list.dirs(reports_dir, recursive = FALSE)
    candidates = candidates[grepl(pat, basename(candidates))]
    if(!is.null(include))
        candidates = candidates[basename(candidates) %in% include]
    candidates = candidates[!basename(candidates) %in% skip]
    if(length(candidates) == 0)
        stop("no reports matched filters under ", reports_dir)

    today = format(Sys.Date(), "%y%m%d")
    log = data.frame(report = character(), status = character(),
                     output = character(), stringsAsFactors = FALSE)

    for(dir in candidates){
        name = basename(dir)
        rmd  = file.path(dir, paste0(name, ".Rmd"))
        if(!file.exists(rmd)){
            log = rbind(log,
                data.frame(report = name, status = "no_rmd",
                           output = NA_character_,
                           stringsAsFactors = FALSE))
            next
        }
        outdir = file.path(results_root, name,
                           if(nzchar(tag)) paste0(today, "_", tag) else today)
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

        message("[knit_all] ", name, " -> ", outdir)
        ok = tryCatch({
            rmarkdown::render(
                input       = rmd,
                output_dir  = outdir,
                output_file = paste0(name, ".html"),
                envir       = new.env(parent = globalenv()),
                quiet       = TRUE)
            "ok"
        }, error = function(e){
            msg = conditionMessage(e)
            if(stop_on_error) stop(msg)
            message("  FAILED: ", msg)
            paste0("error: ", msg)
        })
        log = rbind(log,
            data.frame(report = name, status = ok, output = outdir,
                       stringsAsFactors = FALSE))
    }
    log
}
