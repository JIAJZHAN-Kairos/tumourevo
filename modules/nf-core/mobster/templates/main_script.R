#!/usr/bin/env Rscript

parse_args = function(x) {
    x = gsub("\\\\[","",x)
    x = gsub("\\\\]","",x)
    # giving errors when we have lists like c(xxx, xxx) since it will separate it
    # args_list = unlist(strsplit(x, ', ')[[1]])
    args_list = unlist(strsplit(x, ", (?=[^)]*(?:\\\\(|\$))", perl=TRUE))
    args_vals = lapply(args_list, function(x) {
        x_splt = strsplit(x, split=":")[[1]]
        c(x_splt[1],  paste(x_splt[2:length(x_splt)], collapse=":"))
    })

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals = lapply(args_vals, function(z){ length(z) = 2; z})

    parsed_args = structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

opt = list(
    prefix = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')
)
args_opt = parse_args('$task.ext.args')
for ( ao in names(args_opt)) opt[[ao]] = args_opt[[ao]]


# Script #####

library(CNAqc)
library(mobster)
library(dplyr)
library(cli)
library(ggplot2)

# BRCA/SASH guard: subset_by_segment_karyotype() errors when a sample has no
# 1:1 segments (common after WGD). Return an empty mutation table instead so the
# downstream guard can skip the sample with a recorded reason.
diploid_mutations = function(obj) {
  tryCatch(obj %>% CNAqc::subset_by_segment_karyotype("1:1") %>% CNAqc::Mutations(),
           error = function(e) {
             cli::cli_alert_warning(paste0("No 1:1 segments: ", conditionMessage(e)))
             CNAqc::Mutations(obj)[0, ]
           })
}

description = "$meta.patient"
samples = strsplit(x="$meta.tumour_sample", ",") %>% unlist()  # list of samples

## read mCNAqc object
if ( grepl(".rds\$", tolower("$rds_join")) ) {
  obj = readRDS("$rds_join")
  if (as.logical(opt[["qc_chr"]])){
    chroms <- intersect(names(obj), paste0("chr", 1:22))
    input_table <- lapply(chroms,function(c){
      original = obj[[c]] %>% get_sample(sample=samples, which_obj="original")
      joint_table_chr = lapply(names(original),
                               function(sample_name) {
                                 purity = original[[sample_name]][["purity"]]
                                 table_s = original[[sample_name]] %>%
                                   # CNAqc::subset_by_segment_karyotype("1:1") %>%
                                   CNAqc::Mutations() %>%
                                   dplyr::filter(karyotype=="1:1") %>%
                                   dplyr::mutate(sample_id=sample_name, purity=purity) %>%
                                   dplyr::filter(blacklisted==FALSE)
                                 if (nrow(table_s) == 0) {
                                   cli::cli_alert_warning("Sample {sample_name} has no diploid mutations!")
                                 }
                                 return(table_s)
                               }) %>% dplyr::bind_rows()
    })  %>% dplyr::bind_rows()
  } else {
    if (class(obj) == "m_cnaqc") {
      original = obj %>% get_sample(sample=samples, which_obj="original")
      input_table = lapply(names(original),
                           function(sample_name) {
                             purity = original[[sample_name]][["purity"]]
                             table_s = diploid_mutations(original[[sample_name]]) %>%
                               dplyr::mutate(sample_id=sample_name, purity=purity) %>%
                               dplyr::filter(blacklisted==FALSE)
                             if (nrow(table_s) == 0) {
                               cli::cli_alert_warning("Sample {sample_name} has no diploid mutations!")
                             }
                             return(table_s)
                           }) %>% dplyr::bind_rows()
    } else {
      cli::cli_abort("Object of class {class($rds_join)} not supported.")
    }
  }
} else {
  input_table = read.csv("$rds_join")
  input_table = input_table %>%
    dplyr::filter(karyotype=="1:1") %>%
    dplyr::filter(chr %in% c(paste0('chr', 1:22), as.character(1:22))) %>%
    dplyr::filter(blacklisted==FALSE)
}

## Function to run a single MOBSTER fit
run_mobster_fit = function(joint_table, descr) {
    purity = unique(joint_table[["purity"]]) %>% as.numeric()

    if(length(purity) > 1)
        cli::cli_abort("Multiple values of purity are found in one sample (sample ID: {descr}).")

    # remove mutations with VAF < min_VAF
    inp_tb = joint_table %>%
        dplyr::filter(VAF < 1) %>%
        dplyr::filter(VAF>=as.numeric(opt[["min_VAF"]])) %>%
        dplyr::filter(karyotype=="1:1")

    mobster_fit(x = inp_tb,
                K = eval(parse(text=opt[["K"]])),
                tail = eval(parse(text=opt[["tail"]])),
                pi_cutoff = as.numeric(opt[["pi_cutoff"]]),
                N_cutoff = as.integer(opt[["n_cutoff"]]),
                parallel = FALSE,
                description = descr)
}

## Run MOBSTER on each sample
lapply(samples, function(sample_name) {

    sample_table = input_table %>% dplyr::filter(sample_id == !!sample_name)
    # BRCA/SASH guard: MOBSTER fits diploid (1:1) mutations only. Polyploid or
    # very low-purity samples can have too few; skip with a recorded reason
    # (exit 3, ignored by the pipeline config) instead of crashing the run.
    skip_mobster = function(reason) {
        msg = paste0("MOBSTER skipped for ", sample_name, ": ", reason)
        cli::cli_alert_warning(msg)
        writeLines(msg, paste0(opt[["prefix"]], "_mobster_skipped.txt"))
        quit(save = "no", status = 3)
    }
    n_diploid = sum(sample_table[["VAF"]] < 1 & sample_table[["VAF"]] >= as.numeric(opt[["min_VAF"]]) &
                    sample_table[["karyotype"]] == "1:1", na.rm = TRUE)
    if (n_diploid < as.integer(opt[["n_cutoff"]]))
        skip_mobster(paste0(n_diploid, " diploid (1:1) mutations, fewer than n_cutoff=", opt[["n_cutoff"]]))
    fit = tryCatch(run_mobster_fit(joint_table=sample_table,
                                   descr=paste(description, sample_name, sep=":")),
                   error = function(e) skip_mobster(paste0("fit failed with ", n_diploid,
                                                           " diploid mutations: ", conditionMessage(e))))

    if (fit[["best"]][["fit.tail"]]) {
        evoparams = evolutionary_parameters(fit)
        fit[["evolutionary_parameters"]] = evoparams
        fit[["best"]][["evolutionary_parameters"]] = evoparams
    }

    best_fit = fit[["best"]]
    plot_fit = plot(best_fit)

    saveRDS(object=fit, file=paste0(opt[["prefix"]], "_mobster_fit.rds"))
    saveRDS(object=best_fit, file=paste0(opt[["prefix"]], "_mobster_best_fit.rds"))
    saveRDS(object=plot_fit, file=paste0(opt[["prefix"]], "_mobster_best_plots.rds"))

    # save report plots
    report_fig = mobster::plot_model_selection(fit)
    saveRDS(report_fig, file=paste0(opt[["prefix"]], "_mobster_report.rds"))
    ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_mobster_report.pdf"), height=210, width=210, units="mm", dpi = 200)
    ggplot2::ggsave(plot=report_fig, filename=paste0(opt[["prefix"]], "_mobster_report.png"), height=210, width=210, units="mm", dpi = 200)
})


# Version export #####
f = file("versions.yml","w")
cnaqc_version = sessionInfo()\$otherPkgs\$CNAqc\$Version
mobster_version = sessionInfo()\$otherPkgs\$mobster\$Version
cli_version = sessionInfo()\$otherPkgs\$cli\$Version
dplyr_version = sessionInfo()\$otherPkgs\$dplyr\$Version
ggplot2_version = sessionInfo()\$otherPkgs\$ggplot2\$Version
writeLines(paste0('"', "$task.process", '"', ":"), f)
writeLines(paste("    CNAqc:", cnaqc_version), f)
writeLines(paste("    mobster:", mobster_version), f)
writeLines(paste("    cli:", cli_version), f)
writeLines(paste("    dplyr:", dplyr_version), f)
writeLines(paste("    ggplot2:", ggplot2_version), f)
close(f)
