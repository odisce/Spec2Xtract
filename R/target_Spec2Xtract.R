#' targets factory to extract MS2 spectra
#'
#' @param filter_irel Filter ions based on relative intensity
#' @inheritParams init_object
#' @inheritParams add_events
#' @inheritParams add_cpd_events
#' @inheritParams add_xics
#' @inheritParams add_best_peaks
#' @inheritParams add_events_to_extract
#' @inheritParams add_spectra
#' @inheritParams add_mspurity
#' @inheritParams add_annot
#' @inheritParams export_tables
#' @return A list of targets objects.
#' @import targets tarchetypes data.table magrittr ggplot2 ggpubr ggrepel
#' @export
#' @examples
#' \dontrun{
#'   targets::tar_dir({
#'     targets::tar_script(
#'       {
#'         library(Spec2Xtract)
#'         library(targets)
#'         list(
#'           target_Spec2Xtract(
#'             files = get_sample_rawfile(),
#'             cpd = Spec2Xtract:::example_cpdlist_realdt,
#'             firstevent = TRUE,
#'             prec_ppm = 10,
#'             minscan = 3,
#'             rt_limit = 2,
#'             ppm = 10,
#'             save_dir = "./report"
#'           )
#'         )
#'       },
#'       ask = FALSE
#'     )
#'     ## Run pipeline
#'     targets::tar_make()
#'   })
#' }
target_Spec2Xtract <- function(
    files,
    cpd,
    save_dir = NULL,
    firstevent = TRUE,
    prec_ppm = 10,
    minscan = 3,
    rt_limit = 1,
    filter_irel = 0,
    ppm = 3,
    resources = targets::tar_option_get("resources")
    ) {
  list(
    ## Inputs
    tarchetypes::tar_files_raw(
      "FILE_IN",
      substitute(files),
      deployment = "worker",
      resources = resources
    ),
    tar_target_raw(
      "CPD_IN",
      {
        substitute(cpd)
      },
      deployment = "main"
    ),
    ## Check inputs
    tar_target_raw(
      "F_INFO_dt",
      quote({
        temp <- load_files(FILE_IN)
        return(temp)
      }),
      deployment = "main"
    ),
    tar_target_raw(
      "CPD_INFO_dt",
      quote({
        temp <- fun_check_cpd(CPD_IN)
        temp <- Spec2Xtract::cpd_add_ionsmass(temp)
        return(temp[])
      }),
      deployment = "main"
    ),
    ## Get Files Index ----
    tar_target_raw(
      "F_INDEX",
      quote({
        if (isTRUE(F_INFO_dt$FileExist)) {
          tryCatch(
            {
              temp_dt <- F_INFO_dt$file_path %>%
                get_rawindex() %>%
                parse_index_dt()
              temp_dt[, FileIndex := F_INFO_dt$FileIndex]
              temp_dt[]
            },
            error = function(e) {
              warning(e)
              return(FALSE)
            }
          )
        } else {
          return(FALSE)
        }
      }),
      pattern = quote(map(F_INFO_dt)),
      iteration = "list",
      deployment = "worker",
      resources = resources,
      packages = c("magrittr", "data.table", "Spec2Xtract")
    ),

    ## Get events
    tar_target_raw(
      "F_EVENTS",
      substitute({
        temp_dt <- get_events_types(
          F_INDEX,
          raw_path = F_INFO_dt$file_path,
          firstevent = firstevent
        )

        temp_dt[, FileIndex := F_INFO_dt$FileIndex]
        ## Merge events with precmz diff < isowin * 0.1
        temp_dt[, EventIndexUn := paste(
          FileIndex,
          msLevel,
          spec_energy,
          msLevel,
          spec_coltype,
          spec_polarity,
          isolation_window,
          sep = "_"
        )]

        return(temp_dt[])
      }),
      pattern = quote(map(F_INFO_dt, F_INDEX)),
      iteration = "list",
      deployment = "worker",
      resources = resources,
      packages = c("magrittr", "data.table", "Spec2Xtract")
    ),

    ## For /CPD get events
    tar_target_raw(
      "CPD_EVENTS",
      quote({
        temp_dt <- copy(F_EVENTS)
        if (temp_dt[, .N] < 1) {
          return(NULL)
        }
        temp_dt[, EventIndex := seq_len(.N)]
        temp_dt[,
          c("CpdEventPrec", "CPDCheck") := {
            prec_check <- FALSE
            if (spec_polarity == 0) {
              mz_target <- CPD_INFO_dt$mz_neg
            } else if (spec_polarity == 1) {
              mz_target <- CPD_INFO_dt$mz_pos
            } else {
              warning(
                paste0(
                  "polarity not recognized for this event: ",
                  scanType,
                  " FileIndex: ",
                  FileIndex
                )
              )
              prec_check <- FALSE
              mz_target <- as.numeric(FALSE)
            }

            if (!isFALSE(mz_target)) {
              if (msLevel == 1) {
                prec_check <- TRUE
              } else if (msLevel == 2) {
                prec_check <- mz_target %between%
                  (
                    spec_prec +
                      c(
                        -isolation_window,
                        +isolation_window
                      )
                  )
              }
            }
            .("CpdEventPrec" = mz_target, "CPDCheck" = prec_check)
          },
          by = .(EventIndex, EventIndexUn)
        ]
        temp_dt[, CpdIndex := CPD_INFO_dt$CpdIndex]
        output <- temp_dt[CPDCheck == TRUE, ][]

        if (nrow(output) > 1) {
          ## Combine similar events
          grped_event <- split(output, output$EventIndexUn) %>% {
            lapply(., function(x) {
              grpi_dic_full <- x[, .(EventIndexUn, EventIndex, spec_prec)]
              if (nrow(x) > 1) {
                comp_grid <- combn(seq_len(x[, .N]), 2) %>%
                  t()

                mz_in_range <- lapply(
                  seq_len(nrow(comp_grid)),
                  function(y) {
                    A <- x[comp_grid[y, 1], spec_prec]
                    B <- x[comp_grid[y, 2], spec_prec]
                    iso_win <- mean(
                      c(
                        x[comp_grid[y, 1], isolation_window],
                        x[comp_grid[y, 2], isolation_window]
                      ), na.rm = TRUE
                    )
                    AB <- abs(A - B)
                    out <- data.table(
                      "A" = A,
                      "B" = B,
                      "AB" = AB
                    )

                    if (AB <= iso_win * 0.01) {
                      out[, diffok := TRUE]
                    } else {
                      out[, diffok := FALSE]
                    }
                    return(out)
                  }
                ) %>%
                  rbindlist()

                if (mz_in_range[diffok == TRUE, .N] >= 1) {
                  grpi_dt <- igraph::graph_from_data_frame(
                    mz_in_range[diffok == TRUE, ]
                  ) %>%
                    igraph::components() %>% {
                      as.data.table(
                        .$membership,
                        keep.rownames = TRUE
                      )
                    }
                  setnames(grpi_dt, c("spec_prec", "grpi"))
                  grpi_dt[, spec_prec := as.numeric(spec_prec)]
                  grpi_dic <- merge(
                    grpi_dt,
                    grpi_dic_full,
                    by = "spec_prec",
                    all.y = TRUE
                  )
                  grp_i_max <- grpi_dic[, max(grpi, na.rm = TRUE)] + 1
                } else {
                  grpi_dic <- grpi_dic_full
                  grp_i_max <- 1
                }
              } else {
                grpi_dic <- grpi_dic_full
                grp_i_max <- 1
              }
              if (!"grpi" %in% names(grpi_dic)) {
                grpi_dic[, grpi := 1]
              }
              grp_i_uk <- grpi_dic[order(spec_prec), ][is.na(grpi), ]
              grp_i_uk[, grpi := grp_i_max:(grp_i_max + .N - 1)]

              output_i <- rbind(
                grp_i_uk,
                grpi_dic[!is.na(grpi), ]
              )

              output_i[, EventIndexGrp := paste0(EventIndexUn, "_", grpi)]

              output_export <- merge(
                x,
                output_i[, .(EventIndexGrp, EventIndex, grpi)],
                by = "EventIndex"
              )
              return(output_export)
            })
          } %>%
            rbindlist()

          return(grped_event)
        }
      }),
      pattern = quote(cross(CPD_INFO_dt, F_EVENTS)),
      iteration = "list",
      deployment = "worker",
      resources = resources,
      packages = c("data.table", "magrittr", "Spec2Xtract")
    ),

    ## Get Prec XICs (MS1 or MS2)
    tar_target_raw(
      "CPD_XICs", substitute({
        if (is.null(CPD_EVENTS)) {
          return(NULL)
        }
        ## For current CPD/File, Extract XICs
        xic_dt <- CPD_EVENTS[,
          {
            FileIndex_i <- FileIndex
            CpdIndex_i <- CpdIndex
            fun_get_xic(
              rawpath = F_INFO_dt[FileIndex == FileIndex_i, file_path],
              mz = CpdEventPrec,
              ppm = prec_ppm,
              filter = scanType,
              type = "xic"
            )
          },
          by = .(EventIndexGrp, EventIndexUn, EventIndex, FileIndex, CpdIndex)
        ]
        return(xic_dt)
      }),
      pattern = quote(map(CPD_EVENTS)),
      iteration = "list",
      deployment = "worker",
      resources = resources,
      packages = c("data.table", "magrittr", "Spec2Xtract")
    ),

    tar_target_raw(
      "CPD_PEAKS", substitute({
        if (is.null(CPD_XICs)) {
          return(NULL)
        }
        CPD_peaks <- CPD_XICs[
          ,
          {
            if (.N > minscan) {
              xic_peaks <- get_peaks_from_xic(
                xic_cpd = .SD[order(rt), ],
                minscan = minscan
              )
              if (nrow(xic_peaks) == 0) {
                NULL
              } else {
                xic_peaks
              }
            } else {
              NULL
            }
          },
          by = .(EventIndexGrp, FileIndex)
        ]
        ## Flag peaks inside rt_limit
        if (nrow(CPD_peaks) > 0) {
          CPD_peaks[, PeaksInRange := FALSE]
          cpd_rtmin <- CPD_INFO_dt[
            CpdIndex == CPD_EVENTS[, unique(CpdIndex)],
            rtmin
          ]
          if (!is.null(rt_limit)) {
            rt_lim_range <- cpd_rtmin + (rt_limit * c(-1, 1))
            CPD_peaks[rt %between% rt_lim_range, PeaksInRange := TRUE]
          }
          ## Get scores
          pk_scores_dt <- get_peak_scores(
            peaks_dt = CPD_peaks,
            ref_rtmin = cpd_rtmin
          )

          ## Flag best peak
          pk_scores_dt[, "BestPeak" := FALSE]
          pk_scores_dt[, PeakIndex := seq_len(.N)]
          beast_peak_index <- pk_scores_dt[
            order(-peak_score),
          ][
            PeaksInRange == TRUE,
          ]

          beast_peak_index <- merge(
            beast_peak_index,
            unique(CPD_EVENTS[, .(EventIndexGrp, msLevel)]),
            by = "EventIndexGrp"
          )

          if (beast_peak_index[msLevel == 1, .N]) {
            beast_peak_index <- beast_peak_index[
              msLevel == 1,
            ][
              1,
              PeakIndex
            ]
          } else {
            beast_peak_index <- beast_peak_index[
              1,
              PeakIndex
            ]
          }

          pk_scores_dt[PeakIndex == beast_peak_index, "BestPeak" := TRUE]
          return(pk_scores_dt[])
        }
      }),
      pattern = quote(map(CPD_XICs, CPD_EVENTS)),
      iteration = "list",
      deployment = "worker",
      resources = resources,
      packages = c("data.table", "magrittr", "Spec2Xtract")
    ),

    ## Check MSEvent in best peak range
    tar_target_raw(
      "MSEVENT2EXTRACT", quote({
        # CPD_PEAKS <- tar_read(CPD_PEAKS, 1)[[1]]
        # CPD_EVENTS <- tar_read(CPD_EVENTS, 1)[[1]]
        # tar_load(c(CPD_INFO_dt, F_INDEX))
        if (is.null(CPD_PEAKS) || nrow(CPD_PEAKS) <= 0 || CPD_PEAKS[BestPeak == TRUE, .N] <= 0) {
          ## No peaks detected, returning null
          return(NULL)
        } else {
          ## Get file_i index
          FileIndex_i <- CPD_EVENTS[, unique(FileIndex)]
          F_INDEX_i <- which(sapply(F_INDEX, function(x) {
            x[, unique(FileIndex)]
          }) == FileIndex_i)
          F_INDEX_dt <- F_INDEX[[F_INDEX_i]]

          ## Get CPD Events in peak range
          F_INDEX_dt_i <- F_INDEX_dt[
            CPD_EVENTS[, .(scanType)],
            on = "scanType"
          ][StartTime %between% CPD_PEAKS[BestPeak == TRUE, c(rtmin, rtmax)], ]

          ## For each event get scan closest to the apex
          rtref <- CPD_INFO_dt[CpdIndex == CPD_EVENTS[, unique(CpdIndex)], rtmin]

          F_INDEX_dt_i <- merge(
            F_INDEX_dt_i,
            CPD_EVENTS[, .(scanType, EventIndexUn)],
            by = "scanType"
          )

          ## If no ref rt, get the best peak rt
          if (is.na(rtref)) {
            rtref <- CPD_PEAKS[BestPeak == TRUE, rt]
          }
          scan_to_get <- F_INDEX_dt_i[,
            {
              .SD[which.min(abs(StartTime - rtref)), ]
            },
            by = .(EventIndexUn)
          ]
          scan_to_get[, CpdIndex := CPD_EVENTS[, unique(CpdIndex)]]

          if (nrow(scan_to_get) <= 0) {
            return(NULL)
          } else {
            return(scan_to_get[])
          }
        }
      }),
      pattern = quote(map(CPD_PEAKS, CPD_EVENTS)),
      iteration = "list",
      deployment = "worker",
      resources = resources,
      packages = c("data.table", "magrittr", "Spec2Xtract")
    ),
    tar_target_raw(
      "MSEVENT2EXTRACT_dt", quote({
        data.table::rbindlist(
          MSEVENT2EXTRACT,
          fill = TRUE
        )
      }),
      deployment = "main",
      packages = c("data.table")
    ),

    ## Extract each spectra (iterate over each file for optimized
    ## read access)
    tar_target_raw(
      "MSSPECTRA_ITER", quote({
        if (MSEVENT2EXTRACT_dt[, .N] < 1) {
          return(FALSE)
        } else {
          MSEVENT2EXTRACT_dt[, unique(FileIndex)]
        }
      }),
      deployment = "main",
      packages = c("data.table")
    ),

    tar_target_raw(
      "MSSPECTRA", quote({
        if (isFALSE(MSSPECTRA_ITER)) {
          return(FALSE)
        }
        # tar_load(c(F_INFO_dt, MSEVENT2EXTRACT_dt)) ; MSSPECTRA_ITER <- tar_read(MSSPECTRA_ITER)[[1]]
        file_info_i <- F_INFO_dt[FileIndex == MSSPECTRA_ITER, ]
        scan_to_get_i <- MSEVENT2EXTRACT_dt[FileIndex == MSSPECTRA_ITER, ]
        spectra_db <- get_spectrum_db(
          ms_scan_to_get = unique(scan_to_get_i),
          rawpath = file_info_i$file_path
        )
        ## Add file and cpd info in db
        return(spectra_db)
      }),
      pattern = quote(map(MSSPECTRA_ITER)),
      iteration = "list",
      deployment = "worker",
      resources = resources,
      packages = c("data.table", "magrittr", "Spec2Xtract")
    ),
    tar_target_raw(
      "MSSPECTRA_COMB", quote({
        if (all(sapply(MSSPECTRA, isFALSE))) {
          return(NULL)
        }
        spec_info <- data.table()
        spec_list <- list()
        spec_index <- 0
        for (i in seq_len(length(MSSPECTRA))) {
          spec_info_i <- MSSPECTRA[[i]]$spectra_info_dt
          spec_info_i[, SpectrumIndex := SpectrumIndex + spec_index]
          spec_info <- rbind(
            spec_info,
            spec_info_i
          )
          spec_list <- c(spec_list, MSSPECTRA[[i]]$spectra_db)
          spec_index <- spec_info_i[, max(SpectrumIndex)]
        }

        spec_info[, iter := seq_len(.N)]
        set_trail_zero_cpd <- spec_info[, max(CpdIndex) %>% nchar()]
        set_trail_zero_file <- spec_info[, max(FileIndex) %>% nchar()]
        set_trail_zero_spec <- spec_info[, max(iter) %>% nchar()]
        spec_info[
          ,
          SpecID := paste0(
            "CPD",
            sprintf(
              paste0("%0", set_trail_zero_cpd, ".0f"),
              CpdIndex
            ),
            "_",
            "File",
            sprintf(
              paste0("%0", set_trail_zero_file, ".0f"),
              FileIndex
            ),
            "_",
            "SPEC",
            sprintf(
              paste0("%0", set_trail_zero_spec, ".0f"),
              SpectrumIndex
            )
          )
        ]

        return(
          list(
            "spectra_info_dt" = spec_info[],
            "spectra_db" = spec_list
          )
        )
      }),
      deployment = "main",
      packages = c("data.table", "magrittr", "Spec2Xtract")
    ),

    ## Calculate isolation purity
    tar_target_raw(
      "ISOPURITY", substitute({
        if (is.null(MSSPECTRA_COMB)) {
          return(NULL)
        }
        ## Calculate purity on MS1
        
        .SD <- MSSPECTRA_COMB$spectra_info_dt[SpecID == "CPD01_File06_SPEC015",]
        CpdIndex <- .SD[, unique(CpdIndex)]
        SpectrumIndex <- .SD[, unique(SpectrumIndex)]

        MSSPECTRA_COMB$spectra_info_dt[
          msLevel == 2,
          "isopurity" := {
            ## Check if MS1 available
            temp_i <- .SD[, .(FileIndex, CpdIndex, spec_polarity)] %>% unique()
            temp_events <- MSSPECTRA_COMB$spectra_info_dt[temp_i, on = c("FileIndex", "CpdIndex", "spec_polarity")][msLevel == 1, ]
            iso_purity_val <- as.numeric(NA)
            if (nrow(temp_events) > 0) {
              ## Select one of the MS1 spectra
              CpdIndex_i <- unique(CpdIndex)
              SpectrumIndex_i <- unique(temp_events$SpectrumIndex)
              iso_purity_val <- getpurity_from_spectrum(
                spectrum_ms = MSSPECTRA_COMB$spectra_db[[SpectrumIndex_i]],
                msn_info = .SD[, .(spec_polarity, spec_precmz, spec_isowin)],
                cpd_info = CPD_INFO_dt[CpdIndex == CpdIndex_i, .(mz_pos, mz_neg)],
                prec_ppm = prec_ppm
              )
            }
            if (is.na(iso_purity_val) || is.infinite(iso_purity_val)) {
              iso_purity_val <- 0
            }
            iso_purity_val
          },
          by = .(FileIndex, CpdIndex, spec_scan, SpectrumIndex)
        ]
        return(MSSPECTRA_COMB)
      }),
      deployment = "main",
      packages = c("data.table", "magrittr", "Spec2Xtract")
    ),

    ## Add annotations layer
    tar_target_raw(
      "ANNOT_ITER", quote({
        if (is.null(ISOPURITY)) {
          return(FALSE)
        }
        ISOPURITY$spectra_info_dt[msLevel > 1, .(iter, SpectrumIndex, CpdIndex, FileIndex)][, unique(iter)]
      }),
      deployment = "main",
      packages = c("data.table", "magrittr")
    ),

    tar_target_raw(
      "ANNOT", substitute({
        ## Add annotation to a different layer if multiple cpd for the same spectrum
        if (isFALSE(ANNOT_ITER)) {
          return(NULL)
        }
        spec_info_i <- ISOPURITY$spectra_info_dt[iter == ANNOT_ITER, ]
        cpd_info_i <- CPD_INFO_dt[CpdIndex == spec_info_i$CpdIndex, ]
        spec_annot_i <- Spec2Annot::annotate_mz(
          input_spectrum = ISOPURITY$spectra_db[[spec_info_i$SpectrumIndex]],
          ppm = ppm,
          polarity = spec_info_i$spec_polarity,
          compo = cpd_info_i$elemcomposition,
          use_golden_ratio = FALSE
        )
        annot_i <- spec_annot_i[!is.na(formula), ]
        if (nrow(annot_i) > 0) {
          output <- cbind(
            spec_info_i[, .(iter, CpdIndex, SpectrumIndex, FileIndex)],
            annot_i
          )
          return(output)
        } else {
          NULL
        }
      }),
      deployment = "worker",
      resources = resources,
      pattern = quote(map(ANNOT_ITER)),
      packages = c("data.table", "magrittr", "Spec2Xtract", "Spec2Annot")
    ),

    ## Make plots by compound
    tar_target_raw(
      "XIC_data", quote({
        temp_dt <- rbindlist(CPD_XICs)
        if (nrow(temp_dt) == 0) {
          stop("XICs extraction failed")
          data.table(group = as.integer())
        } else {
          temp_dt %>%
            dplyr::group_by(., CpdIndex) %>%
            targets::tar_group()
        }
      }),
      deployment = "worker",
      resources = resources,
      iteration = "group",
      packages = c("data.table", "dplyr", "targets", "magrittr")
    ),
    tar_target_raw(
      "PEAK_dt", quote({
        peaks_dt <- lapply(
          CPD_PEAKS,
          function(x) {
            if ("BestPeak" %in% names(x)) {
              x[BestPeak == TRUE, ]
            } else {
              NULL
            }
          }
        ) %>%
          rbindlist() %>% {
            .[, .(CpdIndex, FileIndex, rtmin, rtmax, rt)]
          }
        return(peaks_dt)
      }),
      deployment = "worker",
      resources = resources,
      packages = c("data.table", "magrittr")
    ),

    tar_target_raw(
      "XIC_ggplot", quote({
        if (is.null(XIC_data)) {
          return(NULL)
        }

        xic_dt <- as.data.table(XIC_data)
        peakdt_i <- PEAK_dt[CpdIndex == xic_dt[, unique(CpdIndex)]]
        if (nrow(peakdt_i) <= 0) {
          message("No peak found for CPD: ", xic_dt[, unique(CpdIndex)])
          return(NULL)
        }

        ## Get peaks trace
        peakdt_i[, peakID := paste0("PK", seq_len(.N))]
        peak_xic <- peakdt_i[, {
          temp_i <- .SD[, .(CpdIndex, FileIndex)]
          rt_range <- c(rtmin, rtmax)
          temp_xic <- xic_dt[temp_i, on = c("CpdIndex", "FileIndex")]
          temp_xic[rt %between% rt_range, ]
        }, by = .(peakID, CpdIndex, FileIndex)][, .(CpdIndex, FileIndex, rt, i, peakID, EventIndex, filter)]

        ## Add events
        if (is.null(ISOPURITY)) {
          return(NULL)
        }

        event_dt <- ISOPURITY$spectra_info_dt[
          CpdIndex == xic_dt[, unique(CpdIndex)],
          .(rt = StartTime, CpdIndex, FileIndex, msLevel, scanType)
        ] %>%
          unique()

        ## Create plot
        rt_lim <- peakdt_i[, c(median(rtmin), median(rtmax))] * c(0.8, 1.2)
        rt_lim[rt_lim <= 0] <- 0
        if (nrow(peakdt_i) == 0) {
          return(NULL)
        }

        plot_out <- ggplot(data = xic_dt, aes(rt, i)) +
          geom_hline(yintercept = 0) +
          geom_vline(
            data = peakdt_i,
            aes(xintercept = rt),
            linetype = 2
          ) +
          geom_line(
            data = xic_dt,
            aes(rt, i, group = filter),
            alpha = 0.5
          ) +
          geom_line(
            data = peak_xic,
            aes(rt, i, group = filter),
            color = "red"
          ) +
          geom_vline(
            data = event_dt,
            aes(xintercept = rt, color = scanType)
          ) +
          facet_grid(FileIndex ~ ., scales = "free_x") +
          theme_bw() +
          xlim(rt_lim) +
          labs(
            title = paste0(
              "XICs and peaks of the compound iter: ",
              xic_dt[, unique(CpdIndex)]
            ),
            x = "Retention time (min.)",
            y = "Intensity"
          )

        return(plot_out)
      }),
      deployment = "worker",
      resources = resources,
      pattern = quote(map(XIC_data)),
      iteration = "list",
      packages = c("data.table", "magrittr", "ggplot2")
    ),

    tar_target_raw(
      "SUMMARY_dt", quote({
        peak_dt <- PEAK_dt[, .(PeakDetected = .N), by = .(CpdIndex)]
        CPD_INFO_dt[, "FileNb" := F_INFO_dt[, .N]]
        if (is.null(ISOPURITY) && is.null(ANNOT)) {
          return(NULL)
        }
        if (is.null(ISOPURITY)) {
          isopurity_dt <- data.table(
            msLevel = as.integer(),
            IsoPurity_mean = as.numeric(),
            CpdIndex = peak_dt[, unique(CpdIndex)]
          )
          msspectra_dt <- data.table(
            MS2Spectra = as.integer(),
            CpdIndex = as.integer()
          )
        } else {
          isopurity_dt <- ISOPURITY$spectra_info_dt[
            msLevel > 1,
            .(IsoPurity_mean = mean(isopurity, na.rm = TRUE)),
            by = .(CpdIndex)
          ]

          msspectra_dt <- ISOPURITY$spectra_info_dt[
            msLevel > 1,
            .(MS2Spectra = .N),
            by = .(CpdIndex)
          ]
        }

        annotatedspec_dt <- ANNOT[
          ,
          .N,
          by = .(CpdIndex, SpectrumIndex)
        ][
          ,
          .(AnnotatedSpectra = .N),
          by = .(CpdIndex)
        ]

        output <- Reduce(
          function(x, y) {
            merge(x, y, by = "CpdIndex", all = TRUE)
          },
          list(
            CPD_INFO_dt,
            peak_dt,
            isopurity_dt,
            msspectra_dt,
            annotatedspec_dt
          )
        )
        return(output)
      }),
      deployment = "main",
      packages = c("data.table", "magrittr")
    ),

    tar_target_raw(
      "SUMMARYEVENT_dt", quote({
        if (is.null(ISOPURITY)) {
          return(NULL)
        }
        col_to_keep <- grep(
          "spec_",
          names(ISOPURITY$spectra_info_dt),
          value = TRUE
        ) %>% {
          c("SpecID", "CpdIndex", "FileIndex", ., "isopurity")
        }

        output <- ISOPURITY$spectra_info_dt[msLevel > 1, ..col_to_keep] %>% {
          merge(
            CPD_INFO_dt,
            .,
            by = "CpdIndex",
            all = TRUE
          )
        } %>% {
          merge(
            .,
            F_INFO_dt[, .(FileIndex, file_name)],
            by = "FileIndex",
            all.x = TRUE
          )
        }

        return(output)
      }),
      deployment = "main",
      packages = c("data.table", "magrittr")
    ),

    tar_target_raw(
      "SPECTRA_DB", substitute(
        {
          if (!is.null(ISOPURITY)) {
            spectra_msn_info <- ISOPURITY$spectra_info_dt[msLevel > 1, ]
            output <- lapply(
              seq_len(spectra_msn_info[, .N]),
              function(x) {
                SpectrumIndex_i <- spectra_msn_info[x, SpectrumIndex]
                temp_sp <- ISOPURITY$spectra_db[[SpectrumIndex_i]]
                ## Add annot
                spec_annot <- ANNOT[SpectrumIndex == SpectrumIndex_i, -c("SpectrumIndex", "CpdIndex", "FileIndex")]
                if (nrow(spec_annot) > 0) {
                  mz_to_exclude <- spec_annot[, mz]
                  spec_out <- rbind(
                    spec_annot,
                    temp_sp[!mz %in% mz_to_exclude, ],
                    fill = TRUE
                  )
                } else {
                  spec_out <- temp_sp
                }
                spec_out[, irel := i / max(i)]
                new_colorder <- c("mz", "i", "irel", "formula", "ppm")
                new_colorder <- intersect(new_colorder, names(spec_out))
                setcolorder(spec_out, new_colorder)
                spec_out[, SpectrumIndex := SpectrumIndex_i]
                return(spec_out[irel >= filter_irel, ])
              }
            ) %>%
              rbindlist(., fill = TRUE)
            list(
              "spectra_info_dt" = spectra_msn_info,
              "spectra_db" = split(output, output$SpectrumIndex)
            )
          } else {
            return(NULL)
          }
        }
      ),
      deployment = 'main',
      packages = c("data.table", "magrittr", "Spec2Xtract", "Spec2Annot")
    ),

    ## Export data ----
    ### Tables ----
    tar_target_raw(
      "EXPORT", substitute({
        if (is.null(save_dir)) {
          return(FALSE)
        } else {
          if (!dir.exists(save_dir)) {
            dir.create(save_dir, recursive = TRUE)
          }

          save_l <- list(
            "summary_table" = file.path(save_dir, "Summary.xlsx"),
            "summary_events" = file.path(save_dir, "EventSummary.xlsx"),
            "Spectra_dir" = file.path(save_dir, "spectra"),
            "Spectra_dir_msp" = file.path(save_dir, "spectra", "msp"),
            "Spectra_dir_xlsx" = file.path(save_dir, "spectra", "xlsx")
          )

          save_l[grepl("dir", names(save_l))] %>% {
            sapply(., function(x) {
              if (!dir.exists(x)) {
                dir.create(x, recursive = TRUE)
              }
            })
          }

          if (!is.null(SUMMARY_dt)) {
            openxlsx::write.xlsx(SUMMARY_dt, file = save_l$summary_table)
          }

          if (!is.null(SUMMARYEVENT_dt)) {
            openxlsx::write.xlsx(SUMMARYEVENT_dt, file = save_l$summary_events)
          }

          ## Export individual spectra
          if (!is.null(SPECTRA_DB)) {
            temp <- lapply(
              SPECTRA_DB$spectra_info_dt[, seq_len(.N)],
              function(x) {
                spectra_info_i <- SPECTRA_DB$spectra_info_dt[x, ]
                spec_out <- SPECTRA_DB$spectra_db[[x]]
                ## Save xlsx
                save_path <- file.path(save_l$Spectra_dir_xlsx, paste0(spectra_info_i$SpecID, ".xlsx"))
                temp <- openxlsx::write.xlsx(
                  x = spec_out[order(mz)],
                  file = save_path,
                  overwrite = TRUE
                )
                ## Save msp
                save_path <- file.path(save_l$Spectra_dir_msp, paste0(spectra_info_i$SpecID, ".msp"))
                temp <- save_as_msp(
                  spectra_dt = spec_out,
                  spectra_info = spectra_info_i,
                  cpd_info = CPD_INFO_dt[CpdIndex == spectra_info_i$CpdIndex, ],
                  file_out = save_path
                )
              }
            )
          }
          return(TRUE)
        }
      }),
      deployment = "main",
      packages = c("data.table", "magrittr", "openxlsx")
    ),

    ### Figures: XICs ----
    tar_target_raw(
      "EXPORT_XICs_gg", substitute({
        if (is.null(save_dir) | is.null(XIC_ggplot)) {
          return(FALSE)
        } else {
          ## Create folder
          save_dir_xics <- file.path(save_dir, "figures", "xics")
          if (!dir.exists(save_dir_xics)) {
            dir.create(save_dir_xics, recursive = TRUE)
          }
          ## Add names
          cpd_index_i <- as.data.table(XIC_data)[, unique(CpdIndex)]
          save_path <- file.path(save_dir_xics, paste0("CPD", cpd_index_i, ".png"))
          ggsave(
            filename = save_path,
            plot = XIC_ggplot,
            w = 16,
            h = 9,
            bg = "white"
          )
          return(TRUE)
        }
      }),
      pattern = quote(map(XIC_ggplot, XIC_data)),
      deployment = "main",
      packages = c("data.table", "magrittr", "ggplot2", "ggpubr")
    ),

    ### Figures: Spectra ----
    tar_target_raw(
      "SPECTRA_gg_ITER", quote({
        if (is.null(SPECTRA_DB)) {
          return(FALSE)
        } else {
          seq_len(SPECTRA_DB$spectra_info_dt[, .N])
        }
      }),
      deployment = "main",
      packages = c("data.table")
    ),

    tar_target_raw(
      "EXPORT_SPECTRA_gg", substitute({
        if (isFALSE(SPECTRA_gg_ITER)) {
          return(NULL)
        }
        # tar_load(c(SPECTRA_DB, CPD_INFO_dt, F_INFO_dt))
        save_dir_xics <- file.path(save_dir, "figures", "spectra")
        if (!dir.exists(save_dir_xics)) {
          dir.create(save_dir_xics, recursive = TRUE)
        }
        spectra_info_i <- SPECTRA_DB$spectra_info_dt[SPECTRA_gg_ITER, ]
        spectra_i <- SPECTRA_DB$spectra_db[[SPECTRA_gg_ITER]]
        cpd_info_i <- CPD_INFO_dt[CpdIndex == spectra_info_i$CpdIndex, ]
        file_info_i <- F_INFO_dt[FileIndex == spectra_info_i$FileIndex, ]
        output <- ggplot(spectra_i, aes(mz, irel)) +
          geom_hline(yintercept = 0) +
          geom_vline(
            xintercept = as.numeric(spectra_info_i$spec_precmz),
            linetype = 2,
            color = "red"
          ) +
          geom_text(
            aes(
              x = as.numeric(spectra_info_i$spec_precmz),
              y = +Inf,
              label = "precursor",
              vjust = 1,
              hjust = 1
            ),
            color = "red"
          ) +
          geom_linerange(
            aes(
              ymin = 0, 
              ymax = irel
            )
          ) +
          geom_linerange(
            data = spectra_i[!is.na(formula), ],
            aes(
              ymin = 0,
              ymax = irel
            ),
            color = "#003049"
          ) +
          ggrepel::geom_text_repel(
            data = spectra_i[!is.na(formula), ],
            aes(
              label = paste0(
                formula, "\n",
                sprintf("%0.4f", mz)
              )
            ),
            color = "#669bbc",
            vjust = 1,
            min.segment.length = 0
          ) +
          theme_bw() +
          scale_x_continuous(
            breaks = pretty(spectra_i$mz, n = 10),
            minor_breaks = pretty(spectra_i$mz, n = 100)
          ) +
          labs(
            title = cpd_info_i$compound,
            subtitle = paste0(
              cpd_info_i$elemcomposition, " (precusor m/Z: ",
              sprintf("%0.5f", as.numeric(spectra_info_i$spec_precmz)), ")", "\n",
              "Isolation Purity: ", sprintf("%0.2f", as.numeric(spectra_info_i$isopurity)), "%"
            ),
            x = "m/Z",
            y = "Relative intensity (max)",
            caption = paste0(
              "File: ", file_info_i$file_name, "\n",
              "MSLevel: ", spectra_info_i$msLevel, "\n",
              "Isolation window: ",
              sprintf("%0.4f", as.numeric(spectra_info_i$spec_precmz)), " +- ",
              spectra_info_i$spec_isowin, "\n",
              "Collision Type: ", spectra_info_i$spec_coltype, "\n",
              "          Energy: ", spectra_info_i$spec_energy, "\n",
              "scanType: ", spectra_info_i$scanType, "\n",
              "ID: ", spectra_info_i$SpecID
            )
          )
        
        new_x <- c(
          spectra_i[, mz],
          as.numeric(spectra_info_i$spec_precmz) + c(-2, +2)
        )
        new_xrange <- range(new_x)
        output <- output +
          scale_x_continuous(
            breaks = pretty(new_x, n = 10),
            minor_breaks = pretty(new_x, n = 100)
          ) +
          xlim(new_xrange[1], new_xrange[2])

        save_path <- file.path(save_dir_xics, paste0(spectra_info_i$SpecID, ".png"))
        ggsave(save_path, plot = output, w = 16, h = 9)
      }),
      pattern = quote(map(SPECTRA_gg_ITER)),
      deployment = "main",
      packages = c("data.table", "magrittr", "ggplot2", "ggpubr", "ggrepel")
    )
  )
}

#' Wrapper to run targets pipeline
#'
#' This function run a full targets pipeline from
#' two paths: one to the directory containing .raw
#' files and one to the compound table.
#' @param files_dir path to the directory containing the
#'                  .raw files
#' @param cpd_path path to the table with compound informations
#' @param ncore Number of parallel workers
#' @inheritParams target_Spec2Xtract
#' @import targets data.table magrittr
#' @importFrom crew crew_controller_local
#' @importFrom tools file_ext
#' @importFrom openxlsx read.xlsx
#' @export
run_Spec2Xtract <- function(
  files_dir,
  cpd_path,
  firstevent,
  prec_ppm,
  minscan,
  rt_limit,
  ppm,
  filter_irel = 0,
  save_dir,
  ncore = 1
) {
  dir.create(save_dir)
  if (!file.exists(cpd_path)) {
    stop("cpd_path doesn't exist it must be a file")
  }
  if (!dir.exists(files_dir)) {
    stop("files_dir doesn't exist it must be a directory with .raw files")
  }

  eval(
    substitute(
      {
        targets::tar_script(
          {
            table_ext <- tools::file_ext(cpd_path)
            if (table_ext == "xlsx") {
              cpd_in <- openxlsx::read.xlsx(cpd_path)
              cpd_dt <- data.table::as.data.table(cpd_in)
            } else if (table_ext %in% c("csv", "tsv", "txt")) {
              cpd_dt <- data.table::fread(cpd_path)
            }

            targets::tar_option_set(
              packages = c("data.table", "magrittr", "Spec2Xtract")
            )

            if (ncore > 1) {
              tar_option_set(
                controller = crew::crew_controller_local(workers = ncore)
              )
            } else {
              tar_option_set(
                controller = NULL
              )
            }

            list(
              Spec2Xtract::target_Spec2Xtract(
                files = list.files(
                  files_dir,
                  pattern = "\\.raw$",
                  full.names = TRUE,
                  ignore.case = TRUE
                ),
                cpd = cpd_dt,
                firstevent = firstevent,
                prec_ppm = prec_ppm,
                minscan = minscan,
                rt_limit = rt_limit,
                ppm = ppm,
                save_dir = save_dir,
                filter_irel = filter_irel
              )
            )
          },
          ask = FALSE
        )
      },
      list(
        cpd_path = eval(cpd_path),
        files_dir = eval(files_dir),
        firstevent = eval(firstevent),
        prec_ppm = eval(prec_ppm),
        minscan = eval(minscan),
        rt_limit = eval(rt_limit),
        ppm = eval(ppm),
        save_dir = eval(save_dir),
        ncore = eval(ncore),
        filter_irel = eval(filter_irel)
      )
    )
  )
  ## Run pipeline
  targets::tar_make()
}