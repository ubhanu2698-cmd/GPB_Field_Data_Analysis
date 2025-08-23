# GPB Field Data Analysis App — Raw vs Mean Aware (Compliant Build)

# ---- Packages ----
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(colourpicker)
library(DT)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(FactoMineR)
library(factoextra)
library(agricolae)
library(GGEBiplots)
library(lavaan)
library(semPlot)
library(lme4)
library(lmerTest)
library(emmeans)
library(Matrix)
library(pheatmap)
library(psych)
library(rlang)
library(gridExtra)
library(metan)

# optional fast CSV
suppressWarnings(requireNamespace("data.table", quietly = TRUE))

# =============================== Helpers ===============================
`%||%` <- function(x, y) if (is.null(x)) y else x
lower_names <- function(x) tolower(gsub("\\s+","_", x))
round_df_if_numeric <- function(x, k = 4) {
  if (!is.data.frame(x)) return(x)
  num_cols <- vapply(x, is.numeric, TRUE)
  x[num_cols] <- lapply(x[num_cols], function(v) round(v, k))
  x
}
.add_sig_stars <- function(df) {
  if (is.null(df) || !NROW(df)) return(df)
  
  # Try to find a p-value column across common namings
  pcol_candidates <- c("Pr(>F)", "Pr(>F)", "Pr..F.", "p.value", "pvalue", "Pr", "Prob", "p")
  p_idx <- which(tolower(names(df)) %in% tolower(pcol_candidates))
  
  if (!length(p_idx)) {
    # Fallback: any numeric column with values in [0,1] that looks like p-values
    p_idx <- which(vapply(df, function(x) is.numeric(x) && all(x >= 0 & x <= 1, na.rm = TRUE), TRUE))
  }
  if (!length(p_idx)) return(df)  # give up quietly if no p-values
  
  pcol <- p_idx[1]
  p <- suppressWarnings(as.numeric(df[[pcol]]))
  
  # Add stars:  *** <0.001, ** <0.01, * <0.05, . <0.1, '' otherwise
  stars <- ifelse(is.na(p), "",
                  ifelse(p < 0.001, "***",
                         ifelse(p < 0.01,  "**",
                                ifelse(p < 0.05,  "*",
                                       ifelse(p < 0.10,  ".", "")))))
  
  df$Signif <- stars
  df
}
safe_message <- function(msg, type = "error", duration = 6) {
  try({
    txt <- tryCatch(
      {
        if (inherits(msg, "condition")) conditionMessage(msg)
        else if (is.language(msg)) deparse1(msg)
        else if (is.list(msg)) paste(capture.output(str(msg, give.attr = FALSE)), collapse = " ")
        else paste(msg)
      },
      error = function(e) paste("Message format error:", conditionMessage(e))
    )
    showNotification(txt, type = type, duration = duration)
  }, silent = TRUE)
}
# CSV/Excel reader (single source of truth)
read_data <- function(file_input) {
  req(file_input)
  ext <- tolower(tools::file_ext(file_input$name))
  if (ext == "csv") {
    if (requireNamespace("data.table", quietly = TRUE)) {
      return(as.data.frame(data.table::fread(file_input$datapath, showProgress = FALSE)))
    } else {
      return(readr::read_csv(file_input$datapath, show_col_types = FALSE) |> as.data.frame())
    }
  } else if (ext %in% c("xls","xlsx")) {
    return(readxl::read_excel(file_input$datapath) |> as.data.frame())
  } else {
    validate(need(FALSE, "Please upload a CSV or Excel file."))
  }
}
# numeric-ish detector (characters that can coerce to numeric)
is_numericish <- function(v) {
  is.numeric(v) || (is.character(v) && !any(is.na(suppressWarnings(as.numeric(v[!is.na(v)])))))
}
numish_cols <- function(df) names(df)[vapply(df, is_numericish, TRUE)]

# ---------- Column-role detection ----------
ROLE_PATTERNS <- list(
  genotype    = list(exact = c("genotype","geno","line","variety","cultivar","entry"),
                     fuzzy = c("^gen(?:otype)?$", "geno", "line", "variet", "cultivar", "entry")),
  replication = list(exact = c("rep","replication","replicate","repl"),
                     fuzzy = c("^rep(?:lication|licate|lications|licates|l)?$", "^r$")),
  block       = list(exact = c("block","blk","subblock","sub_block"),
                     fuzzy = c("block|blk|sub[_ ]?block")),
  environment = list(exact = c("environment","env","site","location"),
                     fuzzy = c("^env(?:ironment)?$", "site", "location", "station")),
  year        = list(exact = c("year","yr","season"),
                     fuzzy = c("^y(?:ear)?$", "season", "^yr$")),
  treatment   = list(exact = c("treatment","trt","factor"),
                     fuzzy = c("treat|^trt$|factor")),
  location    = list(exact = c("location","site","place","station"),
                     fuzzy = c("locat|site|place|station"))
)
detect_col <- function(df, exact = character(), fuzzy = character()) {
  nms <- names(df); low <- lower_names(nms)
  for (e in tolower(exact)) {
    hit <- which(low == e)
    if (length(hit)) return(nms[hit[1]])
  }
  for (pat in fuzzy) {
    hit <- grep(pat, low, perl = TRUE)
    if (length(hit)) return(nms[hit[1]])
  }
  NULL
}
detect_roles <- function(df, roles) {
  out <- list()
  for (r in roles) {
    pat <- ROLE_PATTERNS[[r]]
    out[[r]] <- detect_col(df, pat$exact, pat$fuzzy)
  }
  out
}
role_picker_ui <- function(id_prefix, df, roles) {
  req(df)
  det <- detect_roles(df, roles)
  lapply(roles, function(r) {
    selectInput(
      paste0(id_prefix, "_", r), paste0("Column: ", tools::toTitleCase(r)),
      choices = names(df), selected = det[[r]] %||% names(df)[1]
    )
  })
}
canonicalize <- function(df, map) {
  canon <- df
  if (!is.null(map$genotype)    && map$genotype    %in% names(df)) canon$G <- as.factor(df[[map$genotype]])
  if (!is.null(map$replication) && map$replication %in% names(df)) canon$R <- as.factor(df[[map$replication]])
  if (!is.null(map$block)       && map$block       %in% names(df)) canon$B <- as.factor(df[[map$block]])
  if (!is.null(map$environment) && map$environment %in% names(df)) canon$E <- as.factor(df[[map$environment]])
  if (!is.null(map$year)        && map$year        %in% names(df)) canon$Y <- as.factor(df[[map$year]])
  if (!is.null(map$treatment)   && map$treatment   %in% names(df)) canon$T <- as.factor(df[[map$treatment]])
  if (!is.null(map$location)    && map$location    %in% names(df)) canon$L <- as.factor(df[[map$location]])
  canon
}

# ---------- NEW: mean aggregator to enforce "means-only" tabs ----------
mean_by_scope <- function(df, geno, env = NULL, year = NULL, traits = NULL) {
  stopifnot(geno %in% names(df))
  num_cols <- if (is.null(traits)) names(df)[vapply(df, is.numeric, TRUE)] else traits
  grp <- c(geno, env, year)
  grp <- grp[!is.null(grp) & nzchar(grp)]
  if (!length(grp)) grp <- geno
  df |>
    dplyr::group_by(dplyr::across(all_of(grp))) |>
    dplyr::summarise(dplyr::across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

# =============================== UI ===============================
ui <- dashboardPage(
  dashboardHeader(title = "GPB Field Data Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Boxplots", tabName = "boxplot", icon = icon("chart-bar")),
      menuItem("Histograms", tabName = "histogram", icon = icon("chart-area")),
      menuItem("ANOVA", tabName = "anova", icon = icon("table")),
      menuItem("BLUPs and BLUEs", tabName = "blup_blue", icon = icon("table")),
      menuItem("Genetic Parameters", tabName = "genetic_params", icon = icon("dna")),
      menuItem("Correlation", tabName = "correlation", icon = icon("project-diagram")),
      menuItem("PCA", tabName = "pca", icon = icon("chart-line")),
      menuItem("Clustering", tabName = "clustering", icon = icon("sitemap")),
      menuItem("Path Analysis", tabName = "path", icon = icon("network-wired")),
      menuItem("AMMI Biplot", tabName = "ammi", icon = icon("chart-pie")),
      menuItem("GGE Biplot", tabName = "gge", icon = icon("table")),
      menuItem("Sample Data", tabName = "DemoData", icon = icon("th-large")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  dashboardBody(
    tabItems(
      # -------- Boxplot  --------
      tabItem(tabName = "boxplot",
              fluidRow(
                box(title = "Upload CSV/Excel", width = 4, status = "primary", solidHeader = TRUE,
                    fileInput("box_file", "CSV / Excel", accept = c(".csv", ".xlsx")),
                    uiOutput("box_y_ui"),
                    uiOutput("box_x_ui"),
                    colourInput("box_color", "Pick fill color", value = "#56B4E9"),
                    downloadButton("box_download", "Download PNG")
                ),
                box(title = "Boxplot & Data", width = 8, status = "primary", solidHeader = TRUE,
                    plotOutput("box_plot", height = "450px"),
                    DTOutput("box_table")
                )
              )
      ),
      
      # -------- Histogram  --------
      tabItem(tabName = "histogram",
              fluidRow(
                box(title = "Upload CSV/Excel", width = 4, status = "primary", solidHeader = TRUE,
                    fileInput("hist_file", "CSV / Excel", accept = c(".csv", ".xlsx")),
                    uiOutput("hist_var_ui"),
                    colourInput("hist_color", "Pick fill color", value = "#E69F00"),
                    downloadButton("hist_download", "Download PNG")
                ),
                box(title = "Histogram & Data", width = 8, status = "primary", solidHeader = TRUE,
                    plotOutput("hist_plot", height = "450px"),
                    DTOutput("hist_table")
                )
              )
      ),
      
      # -------- ANOVA  --------
      tabItem(
        tabName = "anova",
        fluidRow(
          box(
            title = "Upload & Options (RAW)", width = 4, solidHeader = TRUE, status = "primary",
            fileInput("anova_file", "Upload CSV / Excel", accept = c(".csv", ".xlsx")),
            uiOutput("anova_map_ui"),
            uiOutput("anova_trait_ui"),
            tags$hr(),
            helpText("Character/logical columns are treated as factors."),
            tags$hr(),
            selectInput("exp_design", "Experimental Design",
                        choices = c("CRD", "RBD", "LSD", "Alpha-Lattice", "Split-Plot", "Augmented"),
                        selected = "CRD"),
            conditionalPanel(
              condition = "input.exp_design == 'Alpha-Lattice'",
              radioButtons(
                "alpha_block_type", "Block Type (Alpha-Lattice)",
                choices = c("Incomplete " = "incomplete",
                            "Complete " = "complete"),
                selected = "incomplete"
              ),
              helpText("Incomplete uses a mixed model (lmer). Complete fits a classical aov with replication + genotype + replication:block.")
            ),
            radioButtons("trial_type", "Trial Type",
                         choices = c("Single Environment" = "single_env", "Multi-Environment" = "multi_env"),
                         selected = "single_env"),
            radioButtons("year_type", "Year Type",
                         choices = c("Single Year" = "single_year", "Multi-Year" = "multi_year"),
                         selected = "single_year"),
            tags$hr(),
            uiOutput("factor_ui"),
            checkboxInput("allow_3way", "Allow 3-way interactions (if available)", value = FALSE),
            uiOutput("interaction_ui"),
            tags$hr(),
            strong("Duncan Test"),
            uiOutput("mc_factor_ui"),
            numericInput("mc_alpha", "Significance level (alpha)",
                         value = 0.05, min = 0.0001, max = 0.2, step = 0.01),
            actionButton("mc_run", "Run Duncan Test", icon = icon("check")),
            br(), br(),
            downloadButton("anova_table_download", "Download ANOVA CSV")
          ),
          box(
            title = "ANOVA Results", width = 8, solidHeader = TRUE, status = "primary",
            tabsetPanel(
              tabPanel("Summary", verbatimTextOutput("anova_summary")),
              tabPanel("Residuals", plotOutput("anova_resid", height = "400px")),
              tabPanel("Multiple Comparisons – Groups", DTOutput("mc_groups"))
            )
          )
        )
      ),
      # -------- BLUPs & BLUEs (RAW) --------
      tabItem(tabName = "blup_blue",
              fluidRow(
                box(title = "Upload & Settings (RAW)", width = 4, solidHeader = TRUE, status = "primary",
                    fileInput("bb_file", "Upload CSV / Excel", accept = c(".csv", ".xlsx")),
                    tags$hr(),
                    strong("Column Mapping"),
                    uiOutput("bb_env_ui"),
                    uiOutput("bb_rep_ui"),
                    uiOutput("bb_block_ui"),
                    uiOutput("bb_gen_ui"),
                    tags$hr(),
                    strong("Environment Options"),
                    radioButtons("bb_env_scope", "Analysis scope",
                                 choices = c("Single Environment" = "single", "Multiple Environments" = "multi"),
                                 selected = "single"),
                    uiOutput("bb_env_filter_ui"),
                    radioButtons("bb_env_is_random", "Treat Environment as",
                                 choices = c("Fixed effect" = "fixed", "Random effect" = "random"),
                                 selected = "fixed"),
                    checkboxInput("bb_include_gxe", "Include G×E (random) in BLUP model", value = TRUE),
                    tags$hr(),
                    uiOutput("bb_traits_ui"),
                    actionButton("bb_run", "Run BLUEs & BLUPs", icon = icon("play")),
                    br(), br(),
                    downloadButton("bb_download_both", "Download Combined (BLUE+BLUP)"),
                    tags$hr(),                                                     # NEW
                    checkboxInput("bb_filter_env",
                                  "For BLUE/BLUP CSV: include only selected environments (if multi-env)",
                                  value = FALSE),
                    div(style="display:flex; gap:8px; flex-wrap:wrap;",
                        downloadButton("bb_download_blues", "Download BLUEs (CSV)"),
                        downloadButton("bb_download_blups", "Download BLUPs (CSV)")
                    ) 
                ),
                box(title = "Results", width = 8, solidHeader = TRUE, status = "primary",
                    tabsetPanel(
                      tabPanel("BLUEs", DTOutput("bb_table_blue")),
                      tabPanel("BLUPs", DTOutput("bb_table_blup")),
                      tabPanel("Logs", verbatimTextOutput("bb_log"))
                    )
                )
              )
      ),
      
      # -------- Genetic Parameters  --------
      tabItem(tabName = "genetic_params",
              fluidRow(
                box(title = "Upload & Settings (RAW)", width = 4, solidHeader = TRUE, status = "primary",
                    fileInput("gp_file", "Upload CSV / Excel", accept = c(".csv", ".xlsx")),
                    tags$hr(), strong("Column Mapping"),
                    uiOutput("gp_env_ui"),
                    uiOutput("gp_rep_ui"),
                    uiOutput("gp_gen_ui"),
                    tags$hr(),
                    uiOutput("gp_traits_ui"),
                    tags$hr(),
                    strong("Genetic Advance (GA) options"),
                    numericInput("gp_i", "Selection intensity (i)", value = 2.06, min = 0.1, step = 0.01),
                    helpText("Typical i: 10%→1.755, 5%→2.063, 1%→2.665. GA = i * sqrt(VPmean) * H2_mean; GAM% = 100*GA/Mean."),
                    actionButton("gp_run", "Compute Genetic Parameters", icon = icon("play")),
                    br(), br(),
                    downloadButton("gp_download", "Download Parameters")
                ),
                box(title = "Results", width = 8, solidHeader = TRUE, status = "primary",
                    div(style = "display:flex; gap:16px; align-items:center; margin-bottom:8px;",
                        checkboxInput("gp_transpose", "Show parameters as rows (transpose)", value = FALSE)
                    ),
                    DTOutput("gp_table"), br(), verbatimTextOutput("gp_log")
                )
              )
      ),
      
      # -------- Correlation  --------
      tabItem(tabName = "correlation",
              fluidRow(
                box(title = "Upload (MEAN-based)", width = 4, status = "primary", solidHeader = TRUE,
                    fileInput("cor_file", "CSV / Excel", accept = c(".csv", ".xlsx")),
                    uiOutput("cor_map_ui"),
                    uiOutput("cor_vars_ui"),
                    downloadButton("cor_table_download", "Download CSV"),
                    downloadButton("cor_plot_download", "Download Plot PNG"),
                    tags$hr(),
                    strong("Partial Correlations"),
                    uiOutput("pcor_target_ui"),
                    uiOutput("pcor_ctrl_ui"),
                    actionButton("pcor_run", "Compute Partial Correlation"),
                    br(), br(),
                    downloadButton("pcor_download", "Download Partial Corr CSV")
                ),
                box(title = "Correlation & Table (on means)", width = 8, status = "primary", solidHeader = TRUE,
                    plotOutput("cor_plot", height = "450px"),
                    DTOutput("cor_table"),
                    br(),
                    h5("Latest Partial Correlation"),
                    DTOutput("correlation_partial_out")
                )
              )
      ),
      
      # -------- PCA  --------
      tabItem(tabName = "pca",
              fluidRow(
                box(title = "Upload (MEAN-based)", width = 4, status = "primary", solidHeader = TRUE,
                    fileInput("pca_file", "CSV / Excel", accept = c(".csv", ".xlsx")),
                    uiOutput("pca_map_ui"),
                    uiOutput("pca_trait_select"),
                    downloadButton("pca_biplot_download", "Download Biplot PNG"),
                    downloadButton("pca_scree_download", "Download Scree PNG"),
                    downloadButton("pca_eigen_download", "Download Eigenvalues CSV"),
                    downloadButton("pca_loadings_download", "Download Loadings CSV")
                ),
                box(title = "PCA Plots (means)", width = 8, status = "primary", solidHeader = TRUE,
                    tabsetPanel(
                      tabPanel("Scree Plot", plotOutput("pca_scree_plot", height = "400px")),
                      tabPanel("Biplot", plotOutput("pca_biplot", height = "600px"))
                    )
                )
              ),
              fluidRow(
                box(title = "Eigenvalues & Eigenvectors", width = 12, status = "primary", solidHeader = TRUE,
                    DTOutput("pca_eigen_table")
                )
              ),
              # Add under the PCA "Eigenvalues & Eigenvectors" row or as a new row:
              
              fluidRow(
                box(title = "Trait Contributions (PC1–PC3)", width = 12, status = "primary", solidHeader = TRUE,
                    p("Percent contribution of each trait to the first three principal components (auto-limits if fewer PCs exist)."),
                    DTOutput("pca_contrib_pc123"),
                    br(),
                    downloadButton("pca_contrib_pc123_download", "Download Contributions (CSV)")
                )
              ),
              
      ),
      # -------- Clustering  --------
      tabItem(
        tabName = "clustering",
        fluidRow(
          # Left: Inputs
          box(
            title = "Upload & Settings (Genotype-means clustering)",
            width = 4, status = "primary", solidHeader = TRUE,
            fileInput("clus_file", "CSV / Excel", accept = c(".csv", ".xlsx")),
            uiOutput("clus_gen_ui"),
            uiOutput("clus_trait_ui"),
            checkboxInput("clus_scale", "Scale traits (z-score)", value = TRUE),
            tags$hr(),
            numericInput("clus_k", "Number of clusters (k)", value = 3, min = 2),
            selectInput("clus_dist", "Distance metric", choices = c("euclidean","manhattan")),
            selectInput("clus_method", "Linkage method", choices = c("ward.D2","complete","average","single")),
            tags$hr(),
            selectInput(
              "clus_dendro_type", "Dendrogram Type",
              choices = c("Genotypes only (rows)" = "rows",
                          "Traits only (columns)" = "cols",
                          "Both (two-way heatmap)" = "both"),
              selected = "rows"
            ),
            tags$hr(),
            div(style="display:flex; gap:8px; flex-wrap:wrap;",
                downloadButton("clus_dendro_download", "Download Plot PNG"),
                downloadButton("clus_table_download",  "Download Clusters CSV")
            )
          ),
          
          # Right: Plot + Clustered table
          box(
            title = "Dendrogram / Heatmap",
            width = 8, status = "primary", solidHeader = TRUE,
            plotOutput("clus_dendro", height = "560px")
          )
        ),
        
        fluidRow(
          box(
            title = "Clustered Genotype Means",
            width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("clus_table")
          )
        ),
        
        fluidRow(
          box(
            title = "Cluster Summary",
            width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("clus_summary"),
            br(),
            downloadButton("clus_summary_download", "Download Cluster Summary CSV")
          )
        )
      ),
      # -------- Path --------
      tabItem(
        tabName = "path",
        
        # Row 1: Uploads + selectors
        fluidRow(
          box(
            title = "Upload Field Data (RAW) & Select Variables",
            width = 8, status = "primary", solidHeader = TRUE,
            helpText("Upload CSV/Excel with raw field observations. Only numeric trait columns are used."),
            fileInput("path_file", "Field data (CSV / Excel)",
                      accept = c(".csv", ".xls", ".xlsx")),
            # ID / ENV / REP depend on uploaded columns
            uiOutput("path_id_ui"),
            uiOutput("path_env_ui"),
            uiOutput("path_rep_ui"),
            # Y / X are static; server updates choices via updateSelectInput()
            selectizeInput(
              "path_dep", "Dependent variable (Y)",
              choices = character(0), selected = NULL,
              options = list(placeholder = "Choose the response (Y)")
            ),
            selectizeInput(
              "path_ind", "Independent variables (X)",
              choices = character(0), selected = NULL, multiple = TRUE,
              options = list(placeholder = "Choose one or more predictors")
            )
          ),
          box(
            title = "Model Options",
            width = 4, status = "info", solidHeader = TRUE,
            checkboxInput("path_compact", "Use compact model (top‑K predictors)", value = FALSE),
            conditionalPanel(
              condition = "input.path_compact == true",
              numericInput("path_top_k", "K (top standardized paths to keep)",
                           value = 6, min = 1, step = 1)
            ),
            checkboxInput("path_show_resid", "Show residuals on diagram", value = TRUE),
            tags$small("If traits are highly collinear, enable the compact model or reduce predictors.")
          ),
          
          box(
            title = "Genotypic Summaries — Diagnostics",
            width = 12, status = "info", solidHeader = TRUE,
            DTOutput("geno_diag_table")
          ),
        ),
        
        # Row 2: Phenotypic path (full width)
        fluidRow(
          box(
            title = "Phenotypic Path",
            width = 12, status = "info", solidHeader = TRUE,
            plotOutput("phen_path_plot", height = "650px"),
            div(style = "display:flex; gap:8px; flex-wrap:wrap; margin-top:10px;",
                downloadButton("phen_path_plot_download", "Download PNG"),
                downloadButton("phen_path_table_download", "Download Coefficients CSV")
            ),
            br(),
            DTOutput("phen_path_table")
          )
        ),
        
        # Row 3: Genotypic path (from BLUPs, same raw file)
        fluidRow(
          box(
            title = "Genotypic Path (BLUP‑based)",
            width = 12, status = "warning", solidHeader = TRUE,
            plotOutput("geno_path_plot", height = "650px"),
            div(style = "display:flex; gap:8px; flex-wrap:wrap; margin-top:10px;",
                downloadButton("geno_path_plot_download", "Download PNG"),
                downloadButton("geno_path_table_download", "Download Coefficients CSV")
            ),
            br(),
            DTOutput("geno_path_table")
          )
        )
      ),
      # -------- AMMI  --------
      tabItem(tabName = "ammi",
              fluidRow(
                box(title = "Upload CSV/Excel (long/RAW → mean G×E)", width = 4, status = "primary", solidHeader = TRUE,
                    fileInput("ammi_file", "CSV / Excel", accept = c(".csv", ".xlsx")),
                    uiOutput("ammi_mapper_ui"),
                    selectInput("ammi_plot_type", "Select AMMI Plot Type(s):",
                                choices = c("AMMI1 Biplot (Mean vs IPCA1)" = 1,
                                            "AMMI2 Biplot (IPCA1 vs IPCA2)" = 2,
                                            "Environment Scores Only" = 4,
                                            "AMMI2 with Polygon" = 6),
                                selected = c(1, 2), multiple = TRUE),
                    uiOutput("ammi_trait_ui"),
                    downloadButton("ammi_plot_download", "Download PNG"),
                    downloadButton("ammi_table_download", "Download CSV")
                ),
                box(title = "AMMI Output (means)", width = 8, status = "primary", solidHeader = TRUE,
                    plotOutput("ammi_plot", height = "500px"),
                    DTOutput("ammi_table")
                )
              )
      ),
      
      # -------- GGE  --------
      tabItem(tabName = "gge",
              fluidRow(
                box(title = "Upload CSV/Excel (long: Env, Genotype, Trait)", width = 4, status = "primary", solidHeader = TRUE,
                    fileInput("gge_file", "CSV / Excel", accept = c(".csv", ".xlsx")),
                    uiOutput("gge_mapper_ui"),
                    uiOutput("gge_trait_ui"),
                    uiOutput("gge_plot_type_ui"),
                    downloadButton("gge_plot_download", "Download PNG"),
                    downloadButton("gge_table_download", "Download CSV")
                ),
                box(title = "GGE Biplot (means)", width = 8, status = "primary", solidHeader = TRUE,
                    plotOutput("gge_plot", height = "600px"),
                    DTOutput("gge_table")
                ),
                br(),
                h4("GGE Rank Tables"),
                tabsetPanel(
                  tabPanel("Genotype Ranking",
                           DTOutput("gge_rankgen_tbl"),
                           br(), downloadButton("gge_rankgen_csv", "Download Genotype Ranks")
                  ),
                  tabPanel("Environment Ranking",
                           DTOutput("gge_rankenv_tbl"),
                           br(), downloadButton("gge_rankenv_csv", "Download Environment Ranks")
                  ),
                  tabPanel("Mean vs Stability",
                           DTOutput("gge_meanstab_tbl"),
                           br(), downloadButton("gge_meanstab_csv", "Download Mean–Stability Table")
                  ),
                  tabPanel("Discrimination vs Representativeness",
                           DTOutput("gge_discrep_tbl"),
                           br(), downloadButton("gge_discrep_csv", "Download Disc–Rep Table")
                  )
                )
              )
      ),
      
      # -------- Demo Data --------
      tabItem(tabName = "DemoData",
              fluidRow(
                box(title = "Demo Data Generator", width = 12, status = "info", solidHeader = TRUE,
                    sidebarLayout(
                      sidebarPanel(
                        numericInput("genotypes", "Number of Genotypes:", value = 10, min = 5),
                        numericInput("replications", "Number of Replications:", value = 3, min = 1),
                        numericInput("traits", "Number of Traits:", value = 3, min = 1, max = 10),
                        selectInput("design_type", "Choose Design:",
                                    choices = c("CRD", "RBD", "LSD", "Alpha Lattice", "Split Plot", "Augmented")),
                        conditionalPanel(condition = "input.design_type == 'Alpha Lattice'", uiOutput("block_size_ui")),
                        conditionalPanel(condition = "input.design_type == 'Augmented'", numericInput("checks", "Number of Checks:", value = 2, min = 1)),
                        actionButton("generate_demo", "Generate Demo Data"),
                        downloadButton("download_demo", "Download CSV")
                      ),
                      mainPanel(DTOutput("demo_preview"))
                    )
                )
              )
      ),
      # ------------- Help Tab  ----------
      tabItem(
        tabName = "help",
        fluidPage(
          h2("Help & User Guide"),    box(title = "Overview", width = 12, status = "info", solidHeader = TRUE,
                                          collapsible = TRUE, collapsed = FALSE,
                                          p("The GPB Field Data Analysis App analyzes field trial data from raw plots or genotype means."),
                                          tags$ul(
                                            tags$li("Exploratory plots (Boxplots, Histograms)"),
                                            tags$li("ANOVA (CRD, RBD, LSD, Alpha-lattice, Split-plot, Augmented)"),
                                            tags$li("BLUEs & BLUPs (with optional G×E)"),
                                            tags$li("Genetic parameters (VG, VE, VP, GCV/PCV, H\u00B2)"),
                                            tags$li("Correlation, PCA, Clustering"),
                                            tags$li("AMMI & GGE biplots for G×E interpretation"),
                                            tags$li("Path analysis (phenotypic/genotypic SEM)")
                                          ),
                                          tags$hr(),
                                          strong("Tab families:"),
                                          tags$ul(
                                            tags$li("RAW tabs: ANOVA, BLUE/BLUP, Genetic parameters, Path, Histogram, Boxplot."),
                                            tags$li("MEAN tabs: Correlation, PCA, Clustering, AMMI, GGE.")
                                          )
          ),
          # Workflow ---------------------------------------------------------
          box(title = "Workflow", width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              tags$ol(
                tags$li(
                  strong("Prepare your file: "),
                  "Raw tabs use LONG format (one row = plot/experimental unit). ",
                  "Mean tabs use genotype means (one row = genotype). ",
                  "Accepted: .csv, .xls, .xlsx."
                ),
                tags$li(strong("Pick a tab & upload: "), "Each tab has its own file input."),
                tags$li(
                  strong("Map columns: "),
                  "Genotype / Treatment; Replication; Block (Alpha-lattice); Environment / Location; Year / Season; Row & Column (LSD). ",
                  "Only design-required roles are enforced."
                ),
                tags$li(strong("Choose trait(s): "), "Select numeric column(s)."),
                tags$li(strong("Run & interpret: "), "View tables/plots; use download buttons (CSV/PNG).")
              )
          ),
          # Tab-by-Tab Notes
          box(title = "Tab-by-Tab Notes", width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              h4("RAW Tabs"),
              tags$ul(
                tags$li(
                  strong("ANOVA: "),
                  "Design selector (CRD, RBD, LSD, Alpha-lattice, Split-plot, Augmented). ",
                  "Outputs ANOVA table; residual diagnostics (Residuals vs Fitted, Q–Q); ",
                  "Duncan test with compact letter display."
                ),
                tags$li(
                  strong("BLUE/BLUP: "),
                  "Map Genotype (+ optional Environment, Rep, Block). ",
                  "BLUEs treat G as fixed; BLUPs treat G (and optional G×E) as random. ",
                  "Environment can be fixed or random; G×E optional."
                ),
                tags$li(
                  strong("Genetic Parameters: "),
                  "Estimates VG, VE, VP (on means), GCV, PCV, H\u00B2; plus descriptive stats."
                ),
                tags$li(
                  strong("Path: "),
                  "SEM path model (phenotypic): choose dependent (Y) and independents (X). ",
                  "Two plot styles + coefficients table."
                ),
                tags$li(
                  strong("Histogram / Boxplot: "),
                  "Quick distribution checks overall or by group; useful to spot outliers/skewness."
                )
              ),
              h4("MEAN Tabs"),
              tags$ul(
                tags$li(
                  strong("Correlation: "),
                  "Aggregates to means by Genotype (optionally within Environment/Year). ",
                  "Heatmap with significance marks and downloadable table."
                ),
                tags$li(
                  strong("PCA: "),
                  "Runs PCA on the means; scree plot + biplot (genotype scores vs trait loadings)."
                ),
                tags$li(
                  strong("Clustering: "),
                  "Clusters genotypes on selected mean traits (optional z-score scaling). ",
                  "Dendrogram(s), two-way heatmap, cluster table, and cluster-wise means/SD."
                ),
                tags$li(
                  strong("AMMI: "),
                  "Maps E, G, Rep; fits via metan::performs_ammi. ",
                  "Plots: AMMI1 (Mean vs IPCA1), AMMI2 (IPCA1 vs IPCA2), Environment scores, Polygon."
                ),
                tags$li(
                  strong("GGE: "),
                  "Builds E×G matrix from mean table; ",
                  "Plots: Which-Won-Where, Mean vs Stability, RankGen, RankEnv, Discrimination–Representativeness."
                )
              )
          ),
          
          # Interpretation Tips ----------------------------------------------
          box(title = "Interpretation Tips", width = 12, status = "warning", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              tags$ul(
                tags$li("ANOVA: inspect p-values for main/interaction effects; check residual plots for homoscedasticity and normality."),
                tags$li("Duncan groups: compact letters indicate significantly different levels."),
                tags$li("BLUEs = adjusted (fixed-effect) genotype means; BLUPs = shrinkage predictors of genetic merit."),
                tags$li("High H\u00B2 and GCV suggest stronger genetic control and more efficient selection."),
                tags$li("AMMI1: grand mean vs IPCA1 (stability axis). AMMI2: interaction patterns among E and G."),
                tags$li("GGE Which‑Won‑Where: sectors show environment groups and their “winning” genotypes.")
              )
          ),
          # Design Requirements 
          box(title = "Design Requirements (Quick Reference)", width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              tags$ul(
                tags$li("CRD: Genotype/Treatment."),
                tags$li("RBD: Genotype/Treatment + Replication."),
                tags$li("LSD: Genotype/Treatment + Row + Column."),
                tags$li("Alpha-lattice: Genotype + Replication + Block (incomplete). Env/Year optional."),
                tags$li("Split-plot: \u22653 factors (e.g., Sowing, Fertilizer, Genotype)."),
                tags$li("Augmented: Treatment + Block; model expects both Year and Location with \u2265 2 levels.")
              )
          ),
          
          # Troubleshooting 
          box(title = "Troubleshooting", width = 12, status = "danger", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              tags$ul(
                tags$li("\u201CNo numeric columns found\u201D \u2192 pick the correct trait column or check your file."),
                tags$li("\u201CFactor has < 2 levels\u201D \u2192 levels may have collapsed after NA removal; verify mapping/filters."),
                tags$li("\u201CVariance is zero\u201D \u2192 trait constant after filtering; choose a different trait or check data.")
              )
          ),
          
          # Abbreviations 
          box(title = "Abbreviations", width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              p("G: Genotype; T: Treatment; E: Environment; Y: Year/Season; R: Replication; B: Block"),
              p("BLUE: Best Linear Unbiased Estimate (fixed G)"),
              p("BLUP: Best Linear Unbiased Predictor (random G)"),
              p("G\u00D7E: Genotype \u00D7 Environment interaction"),
              p("AMMI: Additive Main effects and Multiplicative Interaction"),
              p("GGE: Genotype + Genotype \u00D7 Environment"),
              p("VG, VE, VP: Genetic, Environmental, Phenotypic variances"),
              p("GCV/PCV: Genotypic/Phenotypic Coefficients of Variation")
          ),
          
          # System Requirements 
          box(title = "System Requirements", width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              tags$ul(
                tags$li("R version \u2265 4.5.0"),
                tags$li("Key packages: shiny, shinydashboard, agricolae, lme4, emmeans, FactoMineR, GGEBiplots, lavaan, semPlot, corrplot, DT, dplyr, metan."),
                tags$li("Stable internet if running online.")
              )
          ),
          
          # Contact 
          box(title = "Contact", width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              p("Udaya Bhanu A"),
              p("Email: Ubhanu2698+gpb@gmail.com")
          )
        )
      )   # <-- this closes tabItem("help")
    )     # <-- this closes tabItems(
  )       # <-- this closes dashboardBody(
)         # <-- this closes dashboardPage(



# =============================== SERVER ===============================
server <- function(input, output, session) {
  
  # ===== Boxplot  =====
  box_df <- reactive({ req(input$box_file); read_data(input$box_file) })
  output$box_y_ui <- renderUI({
    req(box_df()); nums <- names(box_df())[vapply(box_df(), is.numeric, TRUE)]
    validate(need(length(nums) > 0, "No numeric columns found.")); 
    selectInput("box_y", "Numeric trait (Y)", choices = nums)
  })
  output$box_x_ui <- renderUI({
    req(box_df()); selectInput("box_x", "Grouping variable (X) - optional",
                               choices = c(None = "", names(box_df())), selected = "")
  })
  build_boxplot <- reactive({
    df <- box_df(); req(input$box_y)
    if (is.null(input$box_x) || input$box_x == "") {
      ggplot(df, aes(x = factor(1), y = .data[[input$box_y]])) +
        geom_boxplot(fill = input$box_color) +
        labs(x = NULL) + theme_minimal() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    } else {
      ggplot(df, aes(x = .data[[input$box_x]], y = .data[[input$box_y]])) +
        geom_boxplot(fill = input$box_color) + theme_minimal()
    }
  })
  output$box_plot  <- renderPlot({ build_boxplot() })
  output$box_table <- renderDT({ datatable(box_df()) })
  output$box_download <- downloadHandler(
    filename = function() "boxplot.png",
    content = function(file) { png(file, 900, 600); print(build_boxplot()); dev.off() }
  )
  
  # ===== Histogram =====
  hist_df <- reactive({ req(input$hist_file); read_data(input$hist_file) })
  output$hist_var_ui <- renderUI({
    req(hist_df()); nums <- names(hist_df())[vapply(hist_df(), is.numeric, TRUE)]
    validate(need(length(nums) > 0, "No numeric columns found.")); 
    selectInput("hist_var", "Select numeric variable", choices = nums)
  })
  build_hist <- reactive({
    df <- hist_df(); req(input$hist_var)
    x <- df[[input$hist_var]]
    validate(need(sum(is.finite(x)) >= 8, "Need at least 8 finite values for a sensible histogram."))
    s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
    p <- ggplot(df, aes(x = .data[[input$hist_var]])) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = input$hist_color, color = "black", alpha = 0.7) +
      theme_minimal()
    if (is.finite(s) && s > 0) p <- p + stat_function(fun = dnorm, args = list(mean = m, sd = s))
    p
  })
  output$hist_plot  <- renderPlot({ build_hist() })
  output$hist_table <- renderDT({ datatable(hist_df()) })
  output$hist_download <- downloadHandler(
    filename = function() "histogram.png",
    content = function(file) { png(file, 900, 600); print(build_hist()); dev.off() }
  )
  
  # ===== ANOVA  =====
  
  factor_cols_ge2 <- function(df) {
    nms <- names(df)
    nms[vapply(df, function(x) is.factor(x) && nlevels(droplevels(x)) >= 2, TRUE)]
  }
  
  # Coerce character/logical -> factor; drop empties
  anova_raw <- reactive({
    req(input$anova_file)
    df <- read_data(input$anova_file)
    for (nm in names(df)) {
      if (is.character(df[[nm]]) || is.logical(df[[nm]])) df[[nm]] <- as.factor(df[[nm]])
      if (is.factor(df[[nm]])) df[[nm]] <- droplevels(df[[nm]])
    }
    df
  })
  
  # UI: trait and factor pickers
  output$anova_trait_ui <- renderUI({
    df <- anova_raw(); req(df)
    nums <- names(df)[vapply(df, is.numeric, TRUE)]
    validate(need(length(nums) > 0, "No numeric columns found."))
    selectInput("anova_trait", "Select Numeric Trait", choices = nums, selected = nums[1])
  })
  output$interaction_ui <- renderUI({
    df <- anova_raw(); req(df)
    facs <- input$factors %||% character()
    if (!length(facs) || length(facs) < 2) {
      return(helpText("Select ≥2 factors to enable interactions."))
    }
    two_way <- combn(facs, 2, FUN = function(x) paste0(x[1], ":", x[2]))
    three_way <- if (isTRUE(input$allow_3way) && length(facs) >= 3)
      combn(facs, 3, FUN = function(x) paste0(x[1], ":", x[2], ":", x[3])) else character(0)
    selectizeInput(
      "interactions", "Interactions",
      choices = c("Two-way" = two_way, "Three-way" = three_way),
      multiple = TRUE, options = list(plugins = list("remove_button"))
    )
  })
  output$factor_ui <- renderUI({
    df <- anova_raw(); req(df)
    facs <- factor_cols_ge2(df)
    validate(need(length(facs) > 0, "No factor columns with ≥2 levels detected."))
    selectizeInput("factors", "Select Main Factors (fixed)", choices = facs,
                   selected = facs, multiple = TRUE,
                   options = list(plugins = list("remove_button")))
  })
  output$anova_map_ui <- renderUI({
    df <- anova_raw(); req(df)
    # Prefer 'environment' role; fall back to 'location'
    det_env  <- detect_roles(df, c("environment"))$environment
    det_loc  <- detect_roles(df, c("location"))$location
    det      <- detect_roles(df, c("genotype","replication","block","year","treatment"))
    
    tagList(
      strong("Column Mapping"),
      # Primary factor for treatment/genotype — used across designs
      selectInput("anova_gen",  "Primary factor (Genotype/Treatment)",choices = names(df),selected = det$genotype %||% det$treatment %||% names(df)[1]),
      fluidRow( 
        column(6,selectInput("anova_rep",  "Replication (optional)",
                             choices = c("", names(df)),
                             selected = det$replication %||% "")),
        column(6,selectInput("anova_block","Block (optional)",
                             choices = c("", names(df)),
                             selected = det$block %||% ""))),
      fluidRow(
        column(6,selectInput("anova_env",  "Environment / Location (optional)",
                             choices = c("", names(df)),
                             selected = det_env %||% det_loc %||% "")),
        column(6,selectInput("anova_year", "Year / Season (optional)",
                             choices = c("", names(df)),
                             selected = det$year %||% ""))),
      helpText("Tip: Map these once per dataset. Designs will use only what they need."),
      
      conditionalPanel(condition = "input.exp_design == LSD ",
                       fluidRow(
                         column(6,selectInput("anova_row", "Row (LSD only)",
                                              choices = c("", names(df)),
                                              selected = detect_col(df, exact = c("Row"), fuzzy = c("^row$","row.*")) %||% "")),
                         column(6,selectInput("anova_col", "Column (LSD only)",
                                              choices = c("", names(df)),
                                              selected = detect_col(df, exact = c("Column","Col"), fuzzy = c("^col(umn)?$","column.*")) %||% ""))
                       )
      )
    )
  })
  # Pick first factors by names if present, else NULL
  .pick <- function(df, nm) if (nm %in% names(df)) nm else NULL
  
  # Tighten data to used columns and drop empty levels
  .prepare_anova_data <- function(df, used_vars, trait) {
    used_vars <- unique(c(trait, used_vars))
    used_vars <- intersect(used_vars, names(df))
    dat <- stats::na.omit(df[, used_vars, drop = FALSE])
    # drop empties
    for (nm in names(dat)) if (is.factor(dat[[nm]])) dat[[nm]] <- droplevels(dat[[nm]])
    # all factors must have ≥2 levels
    bad <- names(dat)[vapply(dat, function(x) is.factor(x) && nlevels(x) < 2, TRUE)]
    validate(need(length(bad) == 0,
                  paste("These factors have < 2 levels in the analyzed data:",
                        paste(bad, collapse = ", "))))
    dat
  }
  
  # Decide scope from columns (overrides radio buttons only if columns are missing)
  .scope_from_data <- function(df, Loc, Year) {
    list(
      hasLoc  = !is.null(Loc)  && nlevels(df[[Loc]])  >= 2,
      hasYear = !is.null(Year) && nlevels(df[[Year]]) >= 2
    )
  }
  
  # Build and fit per design
  .fit_anova_by_design <- function(df, trait, design) {
    Treatment <- input$anova_gen   %||% NULL   # primary factor UI can be either "Genotype" or "Treatment" logically
    Genotype  <- input$anova_gen   %||% NULL
    Rep       <- if (nzchar(input$anova_rep))  input$anova_rep  else NULL
    Block     <- if (nzchar(input$anova_block))input$anova_block else NULL
    Loc       <- if (nzchar(input$anova_env))  input$anova_env  else NULL
    Year      <- if (nzchar(input$anova_year)) input$anova_year else NULL
    Row    <- if (nzchar(input$anova_row)) input$anova_row else .pick(df, "Row")
    Column <- if (nzchar(input$anova_col)) input$anova_col else (.pick(df, "Column") %||% .pick(df, "Col"))
    
    
    # If user left things blank, backstop with your detectors
    if (is.null(Rep))   Rep   <- .pick(df, "Replication") %||% .pick(df, "Rep") %||% .pick(df, "R")
    if (is.null(Block)) Block <- .pick(df, "Block") %||% .pick(df, "Blk") %||% .pick(df, "B")
    if (is.null(Loc))   Loc   <- .pick(df, "Environment") %||% .pick(df, "Location") %||% .pick(df, "Loc") %||% .pick(df, "E")
    if (is.null(Year))  Year  <- .pick(df, "Year") %||% .pick(df, "Y")
    # Primary factor fallback if UI is unset (rare)
    if (is.null(Genotype)) {
      Genotype  <- .pick(df, "Genotype") %||% .pick(df, "G") %||% .pick(df, "Treatment") %||% .pick(df, "Trt") %||% .pick(df, "Factor")
      Treatment <- Genotype
    }
    
    # Split-plot typical names
    Sowing    <- .pick(df, "Sowing")
    Fertilizer<- .pick(df, "Fertilizer")
    
    # Let user-selected factors backstop missing ones
    user_facs <- intersect(input$factors %||% character(), names(df))
    
    # What’s available?
    sc <- .scope_from_data(df, Loc, Year)
    
    # Choose model per design
    if (design == "Augmented") {
      # aov: y ~ T + Y + L + Y:L + T:Y + T:L + Block %in% Y:L
      validate(need(!is.null(Treatment), "Augmented design needs a 'Treatment' column."))
      validate(need(!is.null(Block), "Augmented design needs a 'Block' column."))
      validate(need(sc$hasYear && sc$hasLoc, "Augmented model expects both Year and Location with ≥2 levels."))
      used <- c(Treatment, Year, Loc, Block)
      dat <- .prepare_anova_data(df, used, trait)
      fml <- as.formula(paste0(trait, " ~ ", Treatment, " + ", Year, " + ", Loc,
                               " + ", Year, ":", Loc,
                               " + ", Treatment, ":", Year,
                               " + ", Treatment, ":", Loc,
                               " + ", Block, "%in%", Year, ":", Loc))
      fit <- aov(fml, data = dat)
      return(list(kind = "aov", fit = fit))
    }
    
    if (design == "Alpha-Lattice") {
      validate(need(!is.null(Block), "Alpha-lattice needs a Block column."))
      Genotype <- (Genotype %||% Treatment)
      validate(need(!is.null(Genotype), "Alpha-lattice needs a primary factor (map it above)."))
      validate(need(!is.null(Rep), "Alpha-lattice needs a Replication column."))
      
      # If user didn't map explicitly, backstop common names (already done above, but keep safe)
      Rep   <- Rep   %||% .pick(df, "Replication") %||% .pick(df, "Rep") %||% .pick(df, "R")
      Block <- Block %||% .pick(df, "Block")       %||% .pick(df, "Blk") %||% .pick(df, "B")
      validate(need(!is.null(Rep) && !is.null(Block), "Missing Replication or Block after detection."))
      
      # COMPLETE blocks path (user-chosen)
      if (identical(input$alpha_block_type, "complete")) {
        used <- c(Genotype, Rep, Block)
        dat  <- .prepare_anova_data(df, used, trait)
        
        # Build: y ~ Rep + Genotype + Rep:Block
        fml <- as.formula(paste0(trait, " ~ ", Rep, " + ", Genotype, " + ", Rep, ":", Block))
        fit <- aov(fml, data = dat)
        return(list(kind = "aov", fit = fit))
      }
      
      # INCOMPLETE (alpha) — your existing mixed model
      has_ge2 <- function(d, col) !is.null(col) && col %in% names(d) && nlevels(droplevels(as.factor(d[[col]]))) >= 2
      hasLoc  <- has_ge2(df, Loc)
      hasYear <- has_ge2(df, Year)
      
      used <- c(Genotype, Rep, Block, if (hasLoc) Loc, if (hasYear) Year)
      dat  <- .prepare_anova_data(df, used, trait)
      
      y <- dat[[trait]]
      validate(need(is.numeric(y), "Response must be numeric."))
      
      rand_gen  <- paste0("(1|", Genotype, ")")
      rand_loc  <- if (hasLoc)  paste0("(1|", Loc,  ")") else NULL
      rand_year <- if (hasYear) paste0("(1|", Year, ")") else NULL
      rand_gl   <- if (hasLoc)  paste0("(1|", Genotype, ":", Loc,  ")") else NULL
      rand_gy   <- if (hasYear) paste0("(1|", Genotype, ":", Year, ")") else NULL
      rand_ly   <- if (hasLoc && hasYear) paste0("(1|", Loc, ":", Year, ")") else NULL
      rand_gly  <- if (hasLoc && hasYear) paste0("(1|", Genotype, ":", Loc, ":", Year, ")") else NULL
      rand_repblk <- paste0("(1|", Rep, "/", Block, ")")
      
      rhs_terms <- c(rand_repblk, rand_gen, rand_loc, rand_year, rand_gl, rand_gy, rand_ly, rand_gly)
      f_rhs     <- paste(rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)], collapse = " + ")
      fml       <- as.formula(paste(trait, "~", ifelse(nzchar(f_rhs), f_rhs, "1")))
      
      ctrl <- lme4::lmerControl(
        optimizer = "bobyqa", calc.derivs = FALSE,
        check.nobs.vs.nRE = "ignore", check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev = "ignore"
      )
      fit <- lmerTest::lmer(fml, data = dat, REML = TRUE, control = ctrl)
      return(list(kind = "lmer", fit = fit))
    }
    
    
    if (design == "Split-Plot") {
      # Prefer Sowing:Fertilizer:Genotype if present; else the first 3 selected factors
      facs <- na.omit(c(Sowing, Fertilizer, Genotype%||% Treatment))
      if (length(facs) < 3) {
        facs <- unique(na.omit(c(facs, user_facs)))
        validate(need(length(facs) >= 3,
                      "Split-plot needs three factors (e.g., Sowing, Fertilizer, Genotype). Please select factors in the panel."))
        facs <- facs[1:3]
      }
      used <- facs
      dat <- .prepare_anova_data(df, used, trait)
      fml <- as.formula(paste0(trait, " ~ ", paste(facs, collapse = ":")))
      fit <- aov(fml, data = dat)
      return(list(kind = "aov", fit = fit))
    }
    
    if (design == "CRD") {
      # Completely Randomized: y ~ Treatment or Genotype (whichever exists)
      fac <- Genotype %||% Treatment %||% (user_facs[1] %||% stop("CRD needs a primary factor (map it above)."))
      dat <- .prepare_anova_data(df, fac, trait)
      fml <- as.formula(paste(trait, "~", fac))
      fit <- aov(fml, data = dat)
      return(list(kind = "aov", fit = fit))
    }
    
    if (design == "RBD") {
      # Randomized Block: y ~ Factor + Rep
      fac <- Genotype %||% Treatment %||% (user_facs[1] %||% stop("RBD needs a primary factor (map it above)."))
      validate(need(!is.null(Rep), "RBD requires a Replication factor."))
      dat <- .prepare_anova_data(df, c(fac, Rep), trait)
      fml <- as.formula(paste(trait, "~", paste(c(fac, Rep), collapse = " + ")))
      fit <- aov(fml, data = dat)
      return(list(kind = "aov", fit = fit))
    }
    
    if (design == "LSD") {
      # Latin Square: y ~ Row + Column + Factor
      fac <- Treatment %||% Genotype %||% (user_facs[1] %||% stop("LSD needs a primary factor (map it above)."))
      validate(need(!is.null(Row) && !is.null(Column),
                    "LSD expects 'Row' and 'Column' factors in the data."))
      dat <- .prepare_anova_data(df, c(fac, Row, Column), trait)
      fml <- as.formula(paste(trait, "~", paste(c(Row, Column, fac), collapse = " + ")))
      fit <- aov(fml, data = dat)
      return(list(kind = "aov", fit = fit))
    }
    
    # Fallback: generic aov with chosen factors and interactions
    used_rhs <- unique(c(input$factors %||% character(), input$interactions %||% character()))
    dat <- .prepare_anova_data(df, used_rhs, trait)
    rhs <- if (length(used_rhs)) paste(used_rhs, collapse = " + ") else "1"
    fit <- aov(as.formula(paste(trait, "~", rhs)), data = dat)
    list(kind = "aov", fit = fit)
  }
  
  anova_fit <- reactive({
    df <- anova_raw(); req(df, input$anova_trait, input$exp_design)
    res <- .fit_anova_by_design(df, input$anova_trait, input$exp_design)
    # (you can keep helpers here if you need them, but be sure to return `res`)
    res
  })
  
  output$anova_summary <- renderPrint({
    obj <- tryCatch(anova_fit(), error = function(e) e)
    if (inherits(obj, "error")) {
      cat("ANOVA error:\n", conditionMessage(obj))
      return(invisible())
    }
    fit <- obj$fit
    if (identical(obj$kind, "lmer")) {
      # lmer: show type III if desired; here default summary + anova table
      print(summary(fit))
      cat("\n--- ANOVA (Satterthwaite df) ---\n")
      print(anova(fit))
    } else {
      print(summary(fit))
    }
  })
  
  output$anova_resid <- renderPlot({
    obj <- tryCatch(anova_fit(), error = function(e) NULL); req(obj)
    if (identical(obj$kind, "lmer")) {
      par(mfrow=c(1,2))
      plot(fitted(obj$fit), resid(obj$fit), xlab = "Fitted", ylab = "Residuals",
           main = "Residuals vs Fitted"); abline(h=0,lty=2)
      qqnorm(resid(obj$fit), main = "Normal Q-Q (residuals)"); qqline(resid(obj$fit))
      par(mfrow=c(1,1))
    } else {
      par(mfrow=c(1,2))
      plot(obj$fit, which=1, main="Residuals vs Fitted")
      plot(obj$fit, which=2, main="Normal Q-Q")
      par(mfrow=c(1,1))
    }
  })
  
  output$anova_table_download <- downloadHandler(
    filename = function() paste0("ANOVA_", input$anova_trait, ".csv"),
    content  = function(file) {
      obj <- anova_fit()
      if (identical(obj$kind, "lmer")) {
        # export lmer ANOVA (Satterthwaite)
        tab <- as.data.frame(anova(obj$fit))
        write.csv(tab, file, row.names = TRUE)
      } else {
        s <- summary(obj$fit)
        tabs <- if (is.list(s)) lapply(s, function(tt) as.data.frame(tt)) else list(as.data.frame(s[[1]]))
        out <- dplyr::bind_rows(lapply(seq_along(tabs), function(i) dplyr::mutate(tabs[[i]], Stratum = paste0("Stratum_", i))))
        write.csv(out, file, row.names = TRUE)
      }
    }
  )
  
  
  output$mc_factor_ui <- renderUI({
    df <- anova_raw(); req(df)
    facs <- names(df)[vapply(df, function(x) is.factor(x) && nlevels(droplevels(x)) >= 2, TRUE)]
    validate(need(length(facs) > 0, "No factors detected to compare."))
    sel <- if ("Genotype" %in% facs) "Genotype" else facs[1]
    selectInput("mc_factor", "Factor to compare (Duncan)", choices = facs, selected = sel)
  })
  
  mc_fallback_aov <- reactive({
    df <- anova_raw(); req(df, input$anova_trait, input$mc_factor)
    rhs <- c(input$mc_factor)
    add_if <- function(v) if (!is.null(v) && v %in% names(df)) v else NULL
    rhs <- c(rhs,
             add_if(input$anova_rep),   add_if("Replication"), add_if("Rep"), add_if("R"),
             add_if(input$anova_block), add_if("Block"), add_if("Blk"), add_if("B"),
             add_if(input$anova_env),   add_if("Environment"), add_if("Location"), add_if("Loc"), add_if("E"),
             add_if(input$anova_year),  add_if("Year"), add_if("Y"))
    rhs <- unique(rhs[!is.na(rhs) & nzchar(rhs)])
    fmla <- as.formula(paste(input$anova_trait, "~", paste(unique(rhs), collapse = " + ")))
    aov(fmla, data = df)
  })
  
  mc_results <- eventReactive(input$mc_run, {
    df  <- anova_raw(); req(df, input$anova_trait, input$mc_factor, input$mc_alpha)
    fac <- input$mc_factor; alpha <- input$mc_alpha
    validate(need(is.factor(df[[fac]]) && nlevels(droplevels(df[[fac]])) > 1,
                  "Selected factor must be a factor with at least 2 levels."))
    aov_mod <- tryCatch(mc_fallback_aov(), error = identity)
    validate(need(!inherits(aov_mod, "error"), "Could not build aov() for Duncan test."))
    dk <- tryCatch(agricolae::duncan.test(aov_mod, fac, alpha = alpha, group = TRUE, console = FALSE),
                   error = identity)
    validate(need(!inherits(dk, "error"), paste("Duncan test failed:", dk$message)))
    groups_df <- as.data.frame(dk$groups)
    groups_df[[fac]] <- rownames(groups_df); rownames(groups_df) <- NULL
    groups_df <- groups_df[, c(fac, setdiff(names(groups_df), fac)), drop = FALSE]
    names(groups_df)[names(groups_df) == fac]      <- "Level"
    if ("groups" %in% names(groups_df)) names(groups_df)[names(groups_df) == "groups"] <- "Group"
    if ("means"  %in% names(groups_df)) names(groups_df)[names(groups_df) == "means"]  <- "Mean"
    if ("std"    %in% names(groups_df)) names(groups_df)[names(groups_df) == "std"]    <- "SD"
    list(groups = groups_df)
  })
  
  output$mc_groups <- DT::renderDT({
    res <- mc_results(); req(res)
    grp <- res$groups
    validate(need(!is.null(grp) && NROW(grp) > 0, "No grouping table available."))
    DT::datatable(grp, extensions = "Buttons",
                  options = list(dom = "Bfrtip", buttons = c("copy","csv","excel"), pageLength = 25))
  })
  
  
  # ===== BLUPs & BLUEs  =====
  bb_raw <- reactive({
    req(input$bb_file)
    df <- read_data(input$bb_file)
    names(df) <- gsub("\\s+", "_", names(df))
    df
  })
  output$bb_env_ui   <- renderUI({ df <- bb_raw(); req(df); selectInput("bb_env_col", "Environment column", choices = c(None = "", names(df)),
                                                                        selected = if ("Environment" %in% names(df)) "Environment" else "") })
  output$bb_rep_ui   <- renderUI({ df <- bb_raw(); req(df); selectInput("bb_rep_col", "Replication column", choices = c(None = "", names(df)),
                                                                        selected = if ("Rep" %in% names(df)) "Rep" else "") })
  output$bb_block_ui <- renderUI({ df <- bb_raw(); req(df); selectInput("bb_block_col", "Block column (optional)", choices = c(None = "", names(df)),
                                                                        selected = if ("Block" %in% names(df)) "Block" else "") })
  output$bb_gen_ui   <- renderUI({ df <- bb_raw(); req(df); selectInput("bb_gen_col", "Genotype column", choices = names(df),
                                                                        selected = if ("Genotype" %in% names(df)) "Genotype" else names(df)[1]) })
  output$bb_traits_ui <- renderUI({
    df <- bb_raw(); req(df)
    map_cols <- c(input$bb_env_col, input$bb_rep_col, input$bb_block_col, input$bb_gen_col)
    map_cols <- map_cols[map_cols != ""]
    cand <- setdiff(numish_cols(df), map_cols)
    validate(need(length(cand) > 0, "No numeric response variables found."))
    selectizeInput("bb_traits", "Response variable(s)", choices = cand, selected = cand[1],
                   multiple = TRUE, options = list(plugins = list("remove_button")))
  })
  output$bb_env_filter_ui <- renderUI({
    df <- bb_raw(); req(df)
    if (is.null(input$bb_env_col) || input$bb_env_col == "" || input$bb_env_scope == "single") return(NULL)
    envs <- sort(unique(df[[input$bb_env_col]]))
    selectizeInput("bb_env_pick", "Select environments (used in analysis)", choices = envs, selected = envs,
                   multiple = TRUE, options = list(plugins = list("remove_button")))
  })
  bb_data <- reactive({
    df <- bb_raw(); req(df, input$bb_gen_col)
    G <- as.factor(df[[input$bb_gen_col]])
    E <- if (!is.null(input$bb_env_col)   && input$bb_env_col   != "") as.factor(df[[input$bb_env_col]]) else NULL
    R <- if (!is.null(input$bb_rep_col)   && input$bb_rep_col   != "") as.factor(df[[input$bb_rep_col]]) else NULL
    B <- if (!is.null(input$bb_block_col) && input$bb_block_col != "") as.factor(df[[input$bb_block_col]]) else NULL
    out <- df; out$G <- G
    if (!is.null(E)) out$E <- E
    if (!is.null(R)) out$R <- R
    if (!is.null(B)) out$B <- B
    if (!is.null(out$E) && identical(input$bb_env_scope, "multi")) {
      req(input$bb_env_pick)
      out <- out[out$E %in% input$bb_env_pick, , drop = FALSE]
      out$E <- droplevels(out$E)
    }
    out$G <- droplevels(out$G)
    if (!is.null(out$R)) out$R <- droplevels(out$R)
    if (!is.null(out$B)) out$B <- droplevels(out$B)
    out
  })
  has_gt1 <- function(v) {
    if (is.null(v)) return(FALSE)
    if (!is.factor(v)) v <- factor(v)
    nlevels(droplevels(v)) > 1
  }
  build_formulas <- function(trait, E_gt1, envRandom, includeGxE, R_gt1, B_gt1) {
    rhs_fixed  <- c()
    rhs_random <- c()
    if (E_gt1) {
      if (envRandom) rhs_random <- c(rhs_random, "(1|E)") else rhs_fixed <- c(rhs_fixed, "E")
    }
    if (R_gt1) {
      if (E_gt1) rhs_random <- c(rhs_random, "(1|E:R)") else rhs_random <- c(rhs_random, "(1|R)")
    }
    if (B_gt1) {
      if (E_gt1 && R_gt1) {
        rhs_random <- c(rhs_random, "(1|E:R:B)")
      } else if (R_gt1) {
        rhs_random <- c(rhs_random, "(1|R:B)")
      } else {
        rhs_random <- c(rhs_random, "(1|B)")
      }
    }
    blue_rhs <- c(rhs_fixed, "G", rhs_random)
    if (!length(blue_rhs)) blue_rhs <- "1"
    blue_fml <- as.formula(paste(trait, "~", paste(blue_rhs, collapse = " + ")))
    blup_random <- c(rhs_random, "(1|G)")
    if (E_gt1 && isTRUE(includeGxE)) blup_random <- c(blup_random, "(1|G:E)")
    blup_rhs <- c(rhs_fixed, blup_random)
    if (!length(blup_rhs)) blup_rhs <- "1"
    blup_fml <- as.formula(paste(trait, "~", paste(blup_rhs, collapse = " + ")))
    list(blue = blue_fml, blup = blup_fml)
  }
  bb_results <- eventReactive(input$bb_run, {
    df <- bb_data(); req(df, input$bb_traits)
    G_gt1 <- has_gt1(df$G)
    E_gt1 <- has_gt1(df$E)
    R_gt1 <- has_gt1(df$R)
    B_gt1 <- has_gt1(df$B)
    envRandom  <- identical(input$bb_env_is_random, "random")
    includeGxE <- isTRUE(input$bb_include_gxe)
    scope      <- input$bb_env_scope
    validate(need("G" %in% names(df), "Map Genotype column."))
    if (scope == "multi") validate(need(!is.null(df$E) && E_gt1, "Select an Environment column with >1 levels for multi-environment analysis."))
    lm_ctrl <- lme4::lmerControl(
      optimizer = "bobyqa",
      calc.derivs = FALSE,
      check.nobs.vs.nRE   = "ignore",
      check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nlev  = "ignore"
    )
    if (!is.null(df$G)) df$G <- droplevels(df$G)
    if (!is.null(df$E)) df$E <- droplevels(df$E)
    if (!is.null(df$R)) df$R <- droplevels(df$R)
    if (!is.null(df$B)) df$B <- droplevels(df$B)
    blue_rows <- list(); blup_rows <- list(); logs <- character(0)
    withProgress(message = "Fitting models…", value = 0, {
      n <- length(input$bb_traits)
      for (i in seq_along(input$bb_traits)) {
        incProgress(1/n, detail = input$bb_traits[i])
        trait <- input$bb_traits[i]
        y <- df[[trait]]
        if (!is.numeric(y) && is_numericish(y)) y <- suppressWarnings(as.numeric(y))
        if (!is.numeric(y) || all(is.na(y))) { logs <- c(logs, sprintf("Skipping %s: non-numeric or all NA.", trait)); next }
        if (sd(y, na.rm = TRUE) == 0) { logs <- c(logs, sprintf("Skipping %s: zero variance.", trait)); next }
        df[[trait]] <- y
        if (!G_gt1) { logs <- c(logs, sprintf("Skipping %s: Genotype has < 2 levels.", trait)); next }
        fm <- build_formulas(trait, E_gt1, envRandom, includeGxE, R_gt1, B_gt1)
        blue_fit <- tryCatch({
          has_rand <- grepl("\\|", deparse(fm$blue))
          if (has_rand) lme4::lmer(fm$blue, data = df, REML = TRUE, control = lm_ctrl)
          else stats::lm(fm$blue, data = df)
        }, error = function(e) { logs <<- c(logs, sprintf("BLUE model error (%s): %s", trait, conditionMessage(e))); NULL })
        if (!is.null(blue_fit)) {
          blue_overall <- tryCatch({
            em <- emmeans::emmeans(blue_fit, ~ G)
            out <- as.data.frame(em)
            if (!"SE" %in% names(out)) out$SE <- NA_real_
            out |>
              dplyr::transmute(
                Type = "BLUE", Level = "overall",
                Environment = NA_character_, Trait = trait,
                Genotype = as.character(G), Estimate = emmean, SE = SE
              )
          }, error = function(e) { logs <<- c(logs, sprintf("BLUE emmeans overall failed (%s): %s", trait, conditionMessage(e))); NULL })
          blue_by_env <- NULL
          if (E_gt1) {
            blue_by_env <- tryCatch({
              em2 <- emmeans::emmeans(blue_fit, ~ G | E)
              out <- as.data.frame(em2)
              if (!"SE" %in% names(out)) out$SE <- NA_real_
              out |>
                dplyr::transmute(
                  Type = "BLUE", Level = "by_env",
                  Environment = as.character(E), Trait = trait,
                  Genotype = as.character(G), Estimate = emmean, SE = SE
                )
            }, error = function(e) { logs <<- c(logs, sprintf("BLUE emmeans by-env failed (%s): %s", trait, conditionMessage(e))); NULL })
          }
          blue_rows[[trait]] <- dplyr::bind_rows(blue_overall, blue_by_env)
        }
        blup_fit <- tryCatch({ lme4::lmer(fm$blup, data = df, REML = TRUE, control = lm_ctrl) },
                             error = function(e) { logs <<- c(logs, sprintf("BLUP model error (%s): %s", trait, conditionMessage(e))); NULL })
        if (!is.null(blup_fit)) {
          blup_g <- tryCatch({
            re <- lme4::ranef(blup_fit)[["G"]]
            if (is.null(re)) return(NULL)
            data.frame(
              Type = "BLUP", Level = "overall", Environment = NA_character_, Trait = trait,
              Genotype = rownames(re), Estimate = unname(re[, "(Intercept)"]),
              stringsAsFactors = FALSE, check.names = FALSE
            )
          }, error = function(e) { logs <<- c(logs, sprintf("BLUP ranef(G) failed (%s): %s", trait, conditionMessage(e))); NULL })
          blup_ge <- NULL
          if (E_gt1 && includeGxE) {
            blup_ge <- tryCatch({
              re <- lme4::ranef(blup_fit)[["G:E"]]
              if (is.null(re)) return(NULL)
              ge_names <- rownames(re)
              sp <- strsplit(ge_names, ":", fixed = TRUE)
              Genotype    <- vapply(sp, `[`, "", 1)
              Environment <- vapply(sp, `[`, "", 2)
              data.frame(
                Type = "BLUP", Level = "by_env", Environment = Environment, Trait = trait,
                Genotype = Genotype, Estimate = unname(re[, "(Intercept)"]),
                stringsAsFactors = FALSE, check.names = FALSE
              )
            }, error = function(e) { logs <<- c(logs, sprintf("BLUP ranef(G:E) failed (%s): %s", trait, conditionMessage(e))); NULL })
          }
          blup_rows[[trait]] <- dplyr::bind_rows(blup_g, blup_ge)
        }
        logs <- c(logs, sprintf("Done: %s", trait))
      }
    })
    BLUEs <- dplyr::bind_rows(blue_rows)
    BLUPs <- dplyr::bind_rows(blup_rows)
    if (nrow(BLUEs)) {
      BLUEs$Estimate <- round(BLUEs$Estimate, 4)
      if ("SE" %in% names(BLUEs)) BLUEs$SE <- round(BLUEs$SE, 4)
    }
    if (nrow(BLUPs)) {
      BLUPs$Estimate <- round(BLUPs$Estimate, 4)
      if (!"SE" %in% names(BLUPs)) BLUPs$SE <- NA_real_
    }
    combined <- dplyr::bind_rows(
      dplyr::select(BLUEs, Type, Level, Environment, Trait, Genotype, Estimate, SE),
      dplyr::select(BLUPs, Type, Level, Environment, Trait, Genotype, Estimate, SE)
    )
    list(BLUEs = BLUEs, BLUPs = BLUPs, Combined = combined, Log = logs)
  })
  .bb_filter_env_rows <- function(df, want_filter, picked_envs) {
    if (!want_filter) return(df)
    if (!("Environment" %in% names(df))) return(df)
    if (is.null(picked_envs) || !length(picked_envs)) return(df)
    # keep only by_env rows for the chosen environments
    df[df$Level == "by_env" & df$Environment %in% picked_envs, , drop = FALSE]
  }
  
  output$bb_table_blue <- DT::renderDT({
    req(bb_results()); df <- bb_results()$BLUEs
    validate(need(nrow(df) > 0, "No BLUEs produced."))
    DT::datatable(df, extensions = "Buttons",
                  options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")))
  })
  output$bb_table_blup <- DT::renderDT({
    req(bb_results()); df <- bb_results()$BLUPs
    validate(need(nrow(df) > 0, "No BLUPs produced."))
    DT::datatable(df, extensions = "Buttons",
                  options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")))
  })
  output$bb_log <- renderPrint({ req(bb_results()); cat(paste(bb_results()$Log, collapse = "\n")) })
  output$bb_download_both <- downloadHandler(
    filename = function() paste0("BLUE_BLUP_", Sys.Date(), ".csv"),
    content  = function(file) { req(bb_results()); readr::write_csv(bb_results()$Combined, file, na = "") }
  )
  output$bb_download_blues <- downloadHandler(
    filename = function() paste0("BLUEs_", Sys.Date(), ".csv"),
    content = function(file) {
      req(bb_results()); df <- bb_results()$BLUEs
      validate(need(nrow(df) > 0, "No BLUEs produced."))
      picked <- if (!is.null(input$bb_env_pick)) input$bb_env_pick else NULL
      out <- .bb_filter_env_rows(df, isTRUE(input$bb_filter_env), picked)
      readr::write_csv(out, file, na = "")
    }
  )
  
  output$bb_download_blups <- downloadHandler(
    filename = function() paste0("BLUPs_", Sys.Date(), ".csv"),
    content = function(file) {
      req(bb_results()); df <- bb_results()$BLUPs
      validate(need(nrow(df) > 0, "No BLUPs produced."))
      picked <- if (!is.null(input$bb_env_pick)) input$bb_env_pick else NULL
      out <- .bb_filter_env_rows(df, isTRUE(input$bb_filter_env), picked)
      readr::write_csv(out, file, na = "")
    }
  )
  
  # ===== Genetic Parameters  =====
  gp_raw <- reactive({ req(input$gp_file); read_data(input$gp_file) })
  output$gp_env_ui <- renderUI({ df <- gp_raw(); selectInput("gp_env", "Environment column (optional)", choices = c("", names(df))) })
  output$gp_rep_ui <- renderUI({ df <- gp_raw(); selectInput("gp_rep", "Replication column", choices = names(df)) })
  output$gp_gen_ui <- renderUI({ df <- gp_raw(); selectInput("gp_gen", "Genotype column", choices = names(df)) })
  output$gp_traits_ui <- renderUI({
    df <- gp_raw(); nums <- names(df)[vapply(df, is.numeric, TRUE)]
    validate(need(length(nums) > 0, "No numeric traits found."))
    selectizeInput("gp_traits", "Trait(s)", choices = nums, multiple = TRUE, selected = nums[1])
  })
  gp_results <- eventReactive(input$gp_run, {
    df  <- gp_raw(); req(input$gp_traits, input$gp_gen, input$gp_rep)
    gen <- input$gp_gen; rep <- input$gp_rep
    env <- if (nzchar(input$gp_env)) input$gp_env else NULL
    rows <- list();  logs <- character(0)
    lm_ctrl <- lme4::lmerControl(optimizer = "bobyqa", calc.derivs = FALSE,
                                 check.nobs.vs.nRE = "ignore", check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev = "ignore")
    for (tr in input$gp_traits) {
      if (!(tr %in% names(df))) next
      y <- df[[tr]]
      if (!is.numeric(y)) {
        y2 <- suppressWarnings(as.numeric(y))
        if (any(is.finite(y2))) y <- y2
      }
      if (!is.numeric(y)) { logs <- c(logs, paste("Skipping non-numeric:", tr)); next }
      dn <- data.frame(y = y, 
                       G = factor(df[[gen]]), 
                       R = factor(df[[rep]]))
      
      if (!is.null(env)) dn$E <- factor(df[[env]])
      
      # include G:E for across-environment inference when E exists
      if (!is.null(env)) {
        rand_terms <- c("(1|E)", "(1|E:R)", "(1|G)", "(1|G:E)")
      } else {
        rand_terms <- c("(1|R)", "(1|G)")
      }
      f_vc <- as.formula(paste("y ~", paste(rand_terms, collapse = " + ")))
      fit  <- tryCatch(lme4::lmer(f_vc, data = dn, REML = TRUE, control = lm_ctrl), error = function(e) NULL)
      
      vc <- if (!is.null(fit)) as.data.frame(lme4::VarCorr(fit)) else NULL
      vg  <- if (!is.null(vc) && any(vc$grp == "G"))    vc$vcov[vc$grp == "G"]    else NA_real_
      vge <- if (!is.null(vc) && any(vc$grp == "G:E"))  vc$vcov[vc$grp == "G:E"]  else 0
      ve  <- if (!is.null(vc))                         vc$vcov[vc$grp == "Residual"] else NA_real_
      
      # effective replication (r) per G per E and number of environments (e)
      eff_r <- tryCatch({
        if (is.null(env)) {
          # reps per genotype (fallback: #levels of R)
          length(unique(dn$R))
        } else {
          # median count of distinct R per (G,E)
          as.numeric(median(tapply(dn$R, list(dn$G, dn$E), function(z) length(unique(z)))))
        }
      }, error = function(...) NA_real_)
      
      eff_e <- if (is.null(env)) 1 else nlevels(dn$E)
      
      # Phenotypic variance on means and H2 on entry-means:
      vp_means <- if (!is.na(vg) && !is.na(ve) && !is.na(eff_r) && !is.na(eff_e) && eff_r > 0 && eff_e > 0)
        vg + (vge / eff_e) + (ve / (eff_r * eff_e)) else NA_real_
      h2_means <- if (!is.na(vg) && !is.na(vp_means) && vp_means > 0) vg / vp_means else NA_real_
      mean_val <- mean(y, na.rm = TRUE)
      se_val   <- stats::sd(y, na.rm = TRUE) / sqrt(sum(is.finite(y)))
      sd_val   <- stats::sd(y, na.rm = TRUE)
      min_val  <- suppressWarnings(min(y, na.rm = TRUE))
      max_val  <- suppressWarnings(max(y, na.rm = TRUE))
      cv       <- if (mean_val == 0) NA_real_ else (sd_val / abs(mean_val)) * 100
      gcv      <- if (is.finite(vg) && mean_val != 0) sqrt(vg) / abs(mean_val) * 100 else NA_real_
      pcv      <- if (is.finite(vp_means) && mean_val != 0) sqrt(vp_means) / abs(mean_val) * 100 else NA_real_
      i_val  <- as.numeric(input$gp_i %||% 2.06)  # default to 2.06 if missing
      ga      <- if (!is.na(vp_means) && !is.na(h2_means)) i_val * sqrt(max(vp_means, 0)) * h2_means else NA_real_
      gam_pct <- if (!is.na(ga) && is.finite(mean_val) && mean_val != 0) (ga / abs(mean_val)) * 100 else NA_real_
      env_cv  <- if (!is.na(ve) && is.finite(mean_val) && mean_val != 0) (sqrt(max(ve, 0)) / abs(mean_val)) * 100 else NA_real_
      
      rows[[tr]] <- data.frame(
        Trait = tr,
        Mean = round(mean_val, 3),
        `Std. Error of Mean` = round(se_val, 3),
        `Std. Deviation` = round(sd_val, 3),
        Min = round(min_val, 3),
        Max = round(max_val, 3),
        `Coefficient of Variation (%)` = round(cv, 2),
        `Genetic Variance (VG)` = round(vg, 3),
        `GxE Variance (VGE)` = round(vge, 3),
        `Environmental Variance (VE)` = round(ve, 3),
        `Phenotypic Variance on Means (VPmean)` = round(vp_means, 3),
        `GCV (%)` = round(gcv, 2),
        `PCV (%)` = round(pcv, 2),
        `Broad-Sense Heritability on Means (H2_mean)` = round(h2_means, 3),
        `Genetic Advance (GA)` = round(ga, 3),
        `GAM (%)`               = round(gam_pct, 2),
        `Environmental CV (%)`  = round(env_cv, 2),
        check.names = FALSE
      )
      logs <- c(logs, paste("Done:", tr))
      
    }
    main_tbl <- dplyr::bind_rows(rows)
    list(main = main_tbl, log = logs)
  })
  # Helper that returns the displayed table based on the toggle
  gp_display_table <- reactive({
    req(gp_results())
    df <- gp_results()$main
    if (!isTRUE(input$gp_transpose)) return(df)
    
    # Transpose: parameters as rows, traits as columns
    # Long → Wide: Parameter becomes first column, each Trait becomes a column
    long <- tidyr::pivot_longer(df, -Trait, names_to = "Parameter", values_to = "Value")
    wide <- tidyr::pivot_wider(long, names_from = "Trait", values_from = "Value")
    
    # Optional: order parameters nicely (put Trait first if it slipped in; it won’t here)
    wide
  })
  
  output$gp_table <- DT::renderDT({
    df <- gp_display_table()
    DT::datatable(
      df, rownames = FALSE,
      extensions = "Buttons",
      options = list(dom = "Bfrtip", buttons = c("copy","csv","excel"), pageLength = 25)
    )
  })
  
  output$gp_log <- renderPrint({ req(gp_results()); cat(paste(gp_results()$log, collapse = "\n")) })
  output$gp_download <- downloadHandler(
    filename = function() paste0("GeneticParameters_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(gp_results()$main, file, na = "")
  )
  
  # ===== Correlation  =====
  cor_df <- reactive({ req(input$cor_file); read_data(input$cor_file) })
  output$cor_map_ui <- renderUI({
    df <- cor_df(); req(df)
    det <- detect_roles(df, c("genotype","environment","year"))
    tagList(
      selectInput("cor_gen", "Genotype column", choices = names(df), selected = det$genotype %||% names(df)[1]),
      selectInput("cor_env", "Environment column (optional)", choices = c("", names(df)), selected = det$environment %||% ""),
      selectInput("cor_year", "Year column (optional)", choices = c("", names(df)), selected = det$year %||% ""),
      helpText("This tab uses trait means. If Environment/Year chosen, means are by the keys you select.")
    )
  })
  cor_mean_df <- reactive({
    df <- cor_df(); req(df, input$cor_gen)
    env <- if (nzchar(input$cor_env)) input$cor_env else NULL
    yr  <- if (nzchar(input$cor_year)) input$cor_year else NULL
    mean_by_scope(df, geno = input$cor_gen, env = env, year = yr)
  })
  # -- Partial correlation pickers --
  output$pcor_target_ui <- renderUI({
    df <- cor_mean_df(); req(df)
    nums <- names(df)[vapply(df, is.numeric, TRUE)]
    if (length(nums) < 2) return(helpText("Need ≥2 numeric variables for partial correlation."))
    fluidRow(
      column(6, selectInput("pcor_x", "Trait X", choices = nums)),
      column(6, selectInput("pcor_y", "Trait Y", choices = nums))
    )
  })
  
  output$pcor_ctrl_ui <- renderUI({
    df <- cor_mean_df(); req(df)
    nums <- names(df)[vapply(df, is.numeric, TRUE)]
    selectizeInput("pcor_ctrl", "Control trait(s) (Z)", choices = nums, multiple = TRUE,
                   options = list(plugins = list("remove_button")))
  })
  
  # -- Robust partial correlation via residuals --
  pcor_result <- eventReactive(input$pcor_run, {
    df <- cor_mean_df(); req(df, input$pcor_x, input$pcor_y)
    validate(need(!identical(input$pcor_x, input$pcor_y), "X and Y must be different traits."))
    
    keep <- unique(c(input$pcor_x, input$pcor_y, input$pcor_ctrl %||% character()))
    sub  <- df[, keep, drop = FALSE]
    sub  <- sub[stats::complete.cases(sub), , drop = FALSE]
    validate(need(nrow(sub) >= 5, "Not enough complete rows for partial correlation."))
    
    x <- sub[[input$pcor_x]]
    y <- sub[[input$pcor_y]]
    Z <- sub[, setdiff(keep, c(input$pcor_x, input$pcor_y)), drop = FALSE]
    
    if (ncol(Z) == 0) {
      r <- suppressWarnings(stats::cor(x, y))
      dfree <- nrow(sub) - 2
    } else {
      rx <- stats::resid(stats::lm(x ~ ., data = Z))
      ry <- stats::resid(stats::lm(y ~ ., data = Z))
      r <- suppressWarnings(stats::cor(rx, ry))
      dfree <- nrow(sub) - ncol(Z) - 2
    }
    tval <- r * sqrt(dfree / max(1 - r^2, .Machine$double.eps))
    pval <- 2 * stats::pt(-abs(tval), df = dfree)
    
    data.frame(
      Trait_X  = input$pcor_x,
      Trait_Y  = input$pcor_y,
      Controls = paste(input$pcor_ctrl %||% character(), collapse = ", "),
      Partial_r = round(r, 4),
      t = round(tval, 3),
      df = dfree,
      p_value = signif(pval, 4),
      stringsAsFactors = FALSE
    )
  })
  
  output$correlation_partial_out <- DT::renderDT({
    req(pcor_result())
    DT::datatable(pcor_result(), rownames = FALSE, options = list(dom = 't'))
  })
  
  output$pcor_download <- downloadHandler(
    filename = function() paste0("partial_correlation_", Sys.Date(), ".csv"),
    content  = function(file) readr::write_csv(pcor_result(), file, na = "")
  )
  
  output$cor_vars_ui <- renderUI({
    df <- cor_mean_df(); req(df)
    nums <- names(df)[vapply(df, is.numeric, TRUE)]
    validate(need(length(nums) >= 2, "Need at least 2 numeric variables after averaging."))
    selectInput("cor_vars", "Select numeric variables", choices = nums, multiple = TRUE, selected = nums)
  })
  cor_numeric_df <- reactive({
    req(cor_mean_df(), input$cor_vars)
    df <- cor_mean_df()[, input$cor_vars, drop = FALSE]
    df[, vapply(df, function(x) length(unique(x[!is.na(x)])) > 1, TRUE), drop = FALSE]
  })
  output$cor_plot <- renderPlot({
    df <- cor_numeric_df(); req(ncol(df) > 1)
    corr <- Hmisc::rcorr(as.matrix(df), type = "pearson")
    corrplot(corr$r, method = "color", addCoef.col = "black", tl.cex = 0.8,
             col = colorRampPalette(c("red","white","green"))(200),
             sig.level = c(.001,.01,.05), insig = "label_sig", pch.cex = 1.2, pch.col = "black")
  })
  output$cor_table <- renderDT({
    df <- cor_numeric_df(); req(ncol(df) > 1)
    corr <- Hmisc::rcorr(as.matrix(df), type = "pearson")
    p_mat <- corr$P
    stars <- symnum(p_mat, corr = FALSE, na = FALSE,
                    cutpoints = c(0, .001, .01, .05, 1), symbols = c("***","**","*",""))
    stars_mat <- matrix(as.character(stars), nrow=nrow(p_mat), ncol=ncol(p_mat))
    corr_table <- matrix(paste0(round(corr$r, 3), stars_mat), nrow = nrow(corr$r), ncol = ncol(corr$r))
    rownames(corr_table) <- colnames(corr_table) <- colnames(df)
    datatable(as.data.frame(corr_table), rownames = TRUE, options = list(pageLength = 10))
  })
  output$cor_table_download <- downloadHandler(
    filename = function() paste0("correlation_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- cor_numeric_df()
      corr <- Hmisc::rcorr(as.matrix(df), type = "pearson")
      p_mat <- corr$P
      stars <- symnum(p_mat, corr = FALSE, na = FALSE,
                      cutpoints = c(0, .001, .01, .05, 1), symbols = c("***","**","*",""))
      stars_mat <- matrix(as.character(stars), nrow=nrow(p_mat), ncol=ncol(p_mat))
      corr_table <- matrix(paste0(round(corr$r, 3), stars_mat), nrow = nrow(corr$r), ncol = ncol(corr$r))
      rownames(corr_table) <- colnames(corr_table) <- colnames(df)
      write.csv(corr_table, file, row.names = TRUE)
    }
  )
  output$cor_plot_download <- downloadHandler(
    filename = function() paste0("correlation_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      df <- cor_numeric_df(); req(ncol(df) > 1)
      png(file, 1000, 1300, res = 150)
      corr <- Hmisc::rcorr(as.matrix(df), type = "pearson")
      corrplot(corr$r, method = "color", addCoef.col = "black", tl.cex = 0.8,
               col = colorRampPalette(c("red","white","green"))(200), pch.cex = 1.2, pch.col = "black")
      dev.off()
    }
  )
  
  # ===== PCA  =====
  pca_data <- reactive({ req(input$pca_file); read_data(input$pca_file) })
  output$pca_map_ui <- renderUI({
    df <- pca_data(); req(df)
    det <- detect_roles(df, c("genotype","environment","year"))
    tagList(
      selectInput("pca_gen", "Genotype column", choices = names(df), selected = det$genotype %||% names(df)[1]),
      selectInput("pca_env", "Environment column (optional)", choices = c("", names(df)), selected = det$environment %||% ""),
      selectInput("pca_year", "Year column (optional)", choices = c("", names(df)), selected = det$year %||% "")
    )
  })
  pca_mean_df <- reactive({
    df <- pca_data(); req(df, input$pca_gen)
    env <- if (nzchar(input$pca_env)) input$pca_env else NULL
    yr  <- if (nzchar(input$pca_year)) input$pca_year else NULL
    mean_by_scope(df, geno = input$pca_gen, env = env, year = yr)
  })
  output$pca_trait_select <- renderUI({
    req(pca_mean_df()); nums <- names(pca_mean_df())[vapply(pca_mean_df(), is.numeric, TRUE)]
    validate(need(length(nums) >= 2, "Need at least two numeric traits for PCA (after averaging)."))
    selectInput("pca_traits", "Select Traits for PCA", choices = nums, multiple = TRUE, selected = nums)
  })
  pca_result <- reactive({
    req(input$pca_traits); FactoMineR::PCA(pca_mean_df()[, input$pca_traits, drop = FALSE], graph = FALSE)
  })
  # --- Trait contributions for PC1–PC3 (robust to fewer PCs) ---
  pca_contrib_pc123 <- reactive({
    pr <- pca_result()
    contrib_mat <- pr$var$contrib
    validate(need(!is.null(contrib_mat) && NCOL(contrib_mat) >= 1,
                  "Trait contributions are not available for this PCA."))
    
    # Clamp to the first 3 PCs actually available
    k <- min(3, NCOL(contrib_mat))
    
    sub <- as.data.frame(contrib_mat[, seq_len(k), drop = FALSE])
    
    # Trait names
    tr_names <- rownames(sub)
    if (is.null(tr_names)) tr_names <- paste0("Trait_", seq_len(NROW(sub)))
    sub$Trait <- tr_names
    
    # Rename PC columns to friendlier headers
    pc_headers <- paste0("PC", seq_len(k), "_Contribution_%")
    names(sub)[seq_len(k)] <- pc_headers
    
    # Put Trait first, then PCs
    sub <- sub[, c("Trait", pc_headers), drop = FALSE]
    
    # Round nicely
    num_idx <- vapply(sub, is.numeric, TRUE)
    sub[num_idx] <- lapply(sub[num_idx], function(x) round(x, 2))
    
    rownames(sub) <- NULL
    sub
  })
  
  output$pca_contrib_pc123 <- DT::renderDT({
    DT::datatable(
      pca_contrib_pc123(),
      rownames = FALSE,
      extensions = "Buttons",
      options = list(pageLength = 15, dom = "Bfrtip", buttons = c("copy","csv","excel"))
    )
  })
  
  output$pca_contrib_pc123_download <- downloadHandler(
    filename = function() "PCA_Trait_Contributions_PC1_PC3.csv",
    content  = function(file) readr::write_csv(pca_contrib_pc123(), file, na = "")
  )
  
  output$pca_scree_plot <- renderPlot({ factoextra::fviz_eig(pca_result(), addlabels = TRUE, ylim = c(0, 100)) })
  output$pca_biplot     <- renderPlot({ factoextra::fviz_pca_biplot(pca_result(), repel = TRUE, col.var = "red", col.ind = "blue") })
  output$pca_eigen_table <- DT::renderDT({
    eig <- pca_result()$eig; colnames(eig) <- c("Eigenvalue","Variance (%)","Cumulative Variance (%)"); eig
  })
  output$pca_biplot_download <- downloadHandler(
    filename = function() "PCA_Biplot.png",
    content  = function(file) { png(file, 1200, 900); print(factoextra::fviz_pca_biplot(pca_result(), repel = TRUE, col.var = "red", col.ind = "blue")); dev.off() }
  )
  output$pca_scree_download <- downloadHandler(
    filename = function() "PCA_ScreePlot.png",
    content  = function(file) { png(file, 1200, 900); print(factoextra::fviz_eig(pca_result(), addlabels = TRUE, ylim = c(0, 100))); dev.off() }
  )
  output$pca_eigen_download <- downloadHandler(
    filename = function() "PCA_Eigenvalues.csv",
    content  = function(file) { write.csv(factoextra::get_eig(pca_result()), file, row.names = FALSE) }
  )
  output$pca_loadings_download <- downloadHandler(
    filename = function() "PCA_Loadings.csv",
    content  = function(file) {
      loadings <- as.data.frame(pca_result()$var$coord)
      write.csv(loadings, file, row.names = TRUE)
    }
  )
  
  # ===== Clustering  ============================
  clus_raw <- reactive({
    req(input$clus_file)
    read_data(input$clus_file)
  })
  
  output$clus_gen_ui <- renderUI({
    df <- clus_raw(); req(df)
    det <- detect_roles(df, c("genotype"))
    selectInput(
      "clus_gen", "Genotype column",
      choices  = names(df),
      selected = det$genotype %||% names(df)[1]
    )
  })
  
  output$clus_trait_ui <- renderUI({
    df <- clus_raw(); req(df, input$clus_gen)
    numc <- setdiff(names(df)[vapply(df, is.numeric, TRUE)], input$clus_gen)
    validate(need(length(numc) >= 2, "Need at least 2 numeric traits in the file."))
    selectizeInput(
      "clus_traits", "Traits to use in clustering",
      choices = numc, selected = numc, multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  })
  
  
  clus_gen_means <- reactive({
    df <- clus_raw(); req(df, input$clus_gen, input$clus_traits)
    
    gm <- df |>
      dplyr::group_by(.data[[input$clus_gen]]) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(input$clus_traits), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )
    names(gm)[1] <- "Genotype"
    
    keep <- rowSums(!is.na(gm[, -1, drop = FALSE])) > 0
    gm   <- gm[keep, , drop = FALSE]
    validate(need(nrow(gm) >= 2, "Need at least 2 genotypes after averaging."))
    
    informative <- vapply(gm[, -1, drop = FALSE], function(x) length(unique(x[is.finite(x)])) > 1, TRUE)
    gm <- cbind(gm[1], gm[, which(informative) + 1, drop = FALSE])
    validate(need(ncol(gm) >= 3, "All selected traits are constant/NA after averaging; pick different traits."))
    
    gm
  })
  
  # Matrix used for distance/hclust (optionally z-scored)
  clus_matrix <- reactive({
    gm <- clus_gen_means()
    X  <- as.data.frame(gm[, -1, drop = FALSE])
    if (isTRUE(input$clus_scale)) X <- as.data.frame(scale(X))
    rownames(X) <- gm$Genotype
    X
  })
  
  # hclust on rows (genotypes)
  clus_hclust_rows <- reactive({
    X <- clus_matrix()
    validate(need(ncol(X) >= 2, "Need at least two traits for clustering."))
    d  <- dist(X, method = input$clus_dist)
    hclust(d, method = input$clus_method)
  })
  
  # hclust on columns (traits)
  clus_hclust_cols <- reactive({
    X <- clus_matrix()
    Xt <- t(X)
    d  <- dist(Xt, method = input$clus_dist)
    hclust(d, method = input$clus_method)
  })
  
  # Assign clusters from row dendrogram (used for table + summary)
  clus_results <- reactive({
    hc <- clus_hclust_rows()
    X  <- clus_matrix()
    cl <- cutree(hc, k = input$clus_k)
    tbl <- data.frame(
      Genotype = rownames(X),
      Cluster  = factor(cl),
      X,
      row.names = NULL, check.names = FALSE
    )
    list(hc = hc, table = tbl)
  })
  
  
  output$clus_dendro <- renderPlot({
    type <- input$clus_dendro_type
    X    <- clus_matrix(); req(X)
    
    if (type == "rows") {
      hc <- clus_hclust_rows()
      plot(hc, main = "Dendrogram — Genotypes (rows)", xlab = "", sub = "")
      rect.hclust(hc, k = input$clus_k, border = 2:(input$clus_k + 1))
      return(invisible())
    }
    
    if (type == "cols") {
      hc <- clus_hclust_cols()
      plot(hc, main = "Dendrogram — Traits (columns)", xlab = "", sub = "")
      return(invisible())
    }
    
    # type == "both": clustered heatmap with both dendrograms
    # Use pheatmap and explicitly draw its gtable so it appears reliably in Shiny.
    pal <- colorRampPalette(c("blue","white","red"))(200)
    
    # row annotation with cluster labels
    annot <- data.frame(Cluster = clus_results()$table$Cluster)
    rownames(annot) <- clus_results()$table$Genotype
    
    ph <- pheatmap::pheatmap(
      as.matrix(X),
      cluster_rows = TRUE, cluster_cols = TRUE,
      annotation_row = annot,
      color = pal,
      main  = "Two-way Heatmap (Genotypes + Traits)",
      silent = TRUE, fontsize = 10, border_color = NA
    )
    grid::grid.newpage(); grid::grid.draw(ph$gtable)
  })
  
  output$clus_dendro_download <- downloadHandler(
    filename = function() {
      switch(input$clus_dendro_type,
             rows = "dendrogram_rows.png",
             cols = "dendrogram_columns.png",
             both = "heatmap_twoway.png")
    },
    content = function(file) {
      type <- input$clus_dendro_type
      X    <- clus_matrix()
      
      png(file, width = 1400, height = 1000, res = 120)
      
      if (type == "rows") {
        hc <- clus_hclust_rows()
        plot(hc, main = "Dendrogram — Genotypes (rows)", xlab = "", sub = "")
        rect.hclust(hc, k = input$clus_k, border = 2:(input$clus_k + 1))
      } else if (type == "cols") {
        hc <- clus_hclust_cols()
        plot(hc, main = "Dendrogram — Traits (columns)", xlab = "", sub = "")
      } else {
        pal <- colorRampPalette(c("blue","white","red"))(200)
        annot <- data.frame(Cluster = clus_results()$table$Cluster)
        rownames(annot) <- clus_results()$table$Genotype
        ph <- pheatmap::pheatmap(
          as.matrix(X),
          cluster_rows = TRUE, cluster_cols = TRUE,
          annotation_row = annot,
          color = pal,
          main  = "Two-way Heatmap (Genotypes + Traits)",
          silent = TRUE, fontsize = 10, border_color = NA
        )
        grid::grid.draw(ph$gtable)
      }
      
      dev.off()
    }
  )
  
  
  output$clus_table <- renderDT({
    DT::datatable(
      clus_results()$table,
      options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")),
      extensions = "Buttons"
    )
  })
  
  output$clus_table_download <- downloadHandler(
    filename = function() "clusters_genotype_means.csv",
    content  = function(file) {
      write.csv(clus_results()$table, file, row.names = FALSE)
    }
  )
  
  
  clus_summary <- reactive({
    res <- clus_results()
    gm  <- clus_gen_means()
    
    merged <- dplyr::left_join(
      res$table[, c("Genotype", "Cluster")],
      gm,
      by = "Genotype"
    )
    
    trait_cols <- setdiff(names(merged), c("Genotype", "Cluster"))
    
    summ <- merged |>
      dplyr::group_by(Cluster) |>
      dplyr::summarise(
        `Cluster Genotypes` = paste(sort(unique(Genotype)), collapse = ", "),
        `No. of Genotypes`  = dplyr::n_distinct(Genotype),
        dplyr::across(
          dplyr::all_of(trait_cols),
          list(Mean = ~ mean(.x, na.rm = TRUE),
               SD   = ~ stats::sd(.x, na.rm = TRUE)),
          .names = "{.col}_{.fn}"
        ),
        .groups = "drop"
      )
    
    # Order clusters numerically if possible
    if (is.factor(summ$Cluster) || is.character(summ$Cluster)) {
      suppressWarnings({
        num <- as.numeric(as.character(summ$Cluster))
        if (all(is.finite(num))) summ <- summ[order(num), , drop = FALSE]
      })
    }
    
    num_idx <- vapply(summ, is.numeric, TRUE)
    summ[num_idx] <- lapply(summ[num_idx], function(x) round(x, 4))
    summ
  })
  
  output$clus_summary <- DT::renderDT({
    DT::datatable(
      clus_summary(),
      options = list(pageLength = 10, dom = "Bfrtip", buttons = c("copy","csv","excel")),
      extensions = "Buttons",
      rownames = FALSE
    )
  })
  
  output$clus_summary_download <- downloadHandler(
    filename = function() "cluster_summary_genotype_means.csv",
    content  = function(file) readr::write_csv(clus_summary(), file, na = "")
  )
  
  
  # ---------- Path ----------
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  if (!exists("safe_message")) {
    safe_message <- function(msg, type = c("message","warning","error")[1]) {
      try(shiny::showNotification(msg, type = switch(type, message="message", warning="warning", error="error")), silent = TRUE)
      message(msg)
    }
  }
  
  .path_sig_stars <- function(p) {
    ifelse(is.na(p), "",
           ifelse(p < 0.001, "***",
                  ifelse(p < 0.01,  "**",
                         ifelse(p < 0.05, "*",
                                ifelse(p < 0.10, ".", "")))))
  }
  # Helper to pick sane defaults
  .pick_defaults <- function(cols) {
    if (length(cols) < 2) return(list(dep = NULL, ind = character(0)))
    dep <- cols[1]
    ind <- setdiff(cols, dep)
    ind_sel <- ind[seq_len(min(3, length(ind)))]
    list(dep = dep, ind = ind_sel)
  }
  
  # When numeric columns are known, fill Y/X choices
  observeEvent(path_numeric_df(), {
    df_num <- path_numeric_df()
    cols   <- names(df_num)
    
    validate(need(length(cols) >= 2,
                  "No (or only one) numeric column found. Convert traits to numeric and re-upload.")
    )
    
    defaults <- .pick_defaults(cols)
    
    # Y (dep)
    updateSelectizeInput(
      session, "path_dep",
      choices = cols,
      selected = if (!is.null(input$path_dep) && input$path_dep %in% cols)
        input$path_dep else defaults$dep,
      server = TRUE
    )
    
    # X (ind) — exclude current/selected Y
    dep_now <- input$path_dep %||% defaults$dep
    ind_choices <- setdiff(cols, dep_now)
    keep_sel <- intersect(input$path_ind %||% character(0), ind_choices)
    sel <- if (length(keep_sel)) keep_sel else head(ind_choices, 3)
    
    updateSelectizeInput(
      session, "path_ind",
      choices  = ind_choices,
      selected = sel,
      server = TRUE
    )
  }, ignoreInit = FALSE)
  
  # If user changes Y, refresh X to exclude it (and keep valid selections)
  observeEvent(input$path_dep, {
    req(path_numeric_df())
    cols <- names(path_numeric_df())
    ind_choices <- setdiff(cols, input$path_dep)
    keep_sel <- intersect(input$path_ind %||% character(0), ind_choices)
    sel <- if (length(keep_sel)) keep_sel else head(ind_choices, 3)
    
    updateSelectizeInput(
      session, "path_ind",
      choices  = ind_choices,
      selected = sel,
      server = TRUE
    )
  }, ignoreInit = TRUE)
  
  
  raw_df <- reactive({
    req(input$path_file)
    read_data(input$path_file)   # your existing utility
  })
  
  # Build UI pickers for id/env/rep once data is loaded
  output$path_id_ui <- renderUI({
    df <- raw_df()
    selectInput("id_col", "Genotype column (ID)", choices = names(df))
  })
  output$path_env_ui <- renderUI({
    df <- raw_df()
    selectInput("env_col", "Environment column (optional)", choices = c("", names(df)))
  })
  output$path_rep_ui <- renderUI({
    df <- raw_df()
    selectInput("rep_col", "Replication/Block column (optional)", choices = c("", names(df)))
  })
  
  # Numeric-only view for phenotypic SEM (traits)
  path_numeric_df <- reactive({
    df <- raw_df()
    for (nm in names(df)) {
      if (is.character(df[[nm]])) {
        x <- gsub(",", "", df[[nm]])
        if (all(grepl("^\\s*-?\\d*\\.?\\d*\\s*$", x) | x == "")) {
          suppressWarnings(df[[nm]] <- as.numeric(x))
        }
      }
    }
    num_df <- dplyr::select(df, where(is.numeric))
    validate(need(ncol(num_df) >= 2, "Need at least two numeric columns for path analysis."))
    num_df
  })
  
  # Path model text
  path_model_text <- reactive({
    req(input$path_dep, input$path_ind)
    paste0(input$path_dep, " ~ ", paste(input$path_ind, collapse = " + "))
  })
  
  
  path_fit <- reactive({
    req(path_model_text(), path_numeric_df())
    df <- path_numeric_df()
    
    # zero-variance guard
    if (any(vapply(df[, input$path_ind, drop = FALSE], function(x) var(x, na.rm = TRUE) == 0, TRUE))) {
      safe_message("One or more predictors have zero variance.", "warning")
      return(NULL)
    }
    
    tryCatch(
      lavaan::sem(path_model_text(), data = df, missing = "fiml", estimator = "MLR"),
      error = function(e) { safe_message(paste("Path model error:", e$message), "error"); NULL }
    )
  })
  
  path_std_effects <- reactive({
    fit <- path_fit(); req(fit)
    std <- lavaan::standardizedSolution(fit)
    std[std$op == "~", c("lhs","rhs","est.std","pvalue")]
  })
  
  path_compact_fit <- reactive({
    req(isTRUE(input$path_compact))
    df  <- path_numeric_df(); req(df, input$path_dep, input$path_ind)
    std <- path_std_effects(); req(nrow(std) > 0)
    std_y <- subset(std, lhs == input$path_dep & rhs %in% input$path_ind)
    validate(need(nrow(std_y) > 0, "No standardized paths found for selected variables."))
    k <- min(input$path_top_k %||% 6, nrow(std_y))
    keep <- head(std_y[order(-abs(std_y$est.std)), "rhs", drop = TRUE], k)
    fml <- paste0(input$path_dep, " ~ ", paste(keep, collapse = " + "))
    tryCatch(lavaan::sem(fml, data = df, missing = "fiml", estimator = "MLR"),
             error = function(e) { safe_message(paste("Compact path model error:", e$message), "error"); NULL })
  })
  
  .phen_fit_active <- reactive({
    if (isTRUE(input$path_compact)) path_compact_fit() else path_fit()
  })
  
  .build_path_effect_table <- function(fit, df, y, X) {
    std <- lavaan::standardizedSolution(fit)
    dir_df <- subset(std, op == "~" & lhs == y & rhs %in% X,
                     select = c("rhs","est.std","pvalue"))
    if (is.null(dir_df) || nrow(dir_df) == 0) {
      return(data.frame(
        Predictor = character(), Direct = numeric(), Indirect = numeric(),
        Total = numeric(), SE = numeric(), z = numeric(),
        p_direct = numeric(), Signif = character(), stringsAsFactors = FALSE
      ))
    }
    names(dir_df) <- c("Predictor","Direct","p_direct")
    
    # SE and z
    pe <- lavaan::parameterEstimates(fit)
    se_z <- subset(pe, op == "~" & lhs == y & rhs %in% X,
                   select = c("rhs","se","z","pvalue"))
    if (is.null(se_z) || nrow(se_z) == 0) {
      se_z <- data.frame(Predictor = dir_df$Predictor,
                         SE = NA_real_, z = NA_real_, p_raw = NA_real_)
    } else {
      names(se_z) <- c("Predictor","SE","z","p_raw")
    }
    
    # --- FIX: always define R safely ---
    X_present <- intersect(X, colnames(df))
    if (length(X_present) <= 1) {
      R <- matrix(1, nrow = length(X_present), ncol = length(X_present),
                  dimnames = list(X_present, X_present))
    } else {
      R <- stats::cor(df[, X_present, drop = FALSE],
                      use = "pairwise.complete.obs", method = "pearson")
      R <- as.matrix(R)
      if (is.null(rownames(R))) rownames(R) <- X_present
      if (is.null(colnames(R))) colnames(R) <- X_present
    }
    
    # indirect effects
    dir_vec <- dir_df$Direct; names(dir_vec) <- dir_df$Predictor
    indir <- sapply(dir_df$Predictor, function(i) {
      if (!i %in% rownames(R)) return(0)
      others <- setdiff(colnames(R), i)
      if (length(others) == 0) return(0)
      sum(R[i, others, drop = TRUE] * dir_vec[others], na.rm = TRUE)
    })
    
    out <- dplyr::left_join(dir_df, se_z, by = "Predictor")
    out$Indirect <- indir[match(out$Predictor, names(indir))]
    out$Total    <- out$Direct + out$Indirect
    out$Signif   <- .path_sig_stars(out$p_direct)
    
    out <- dplyr::arrange(out, dplyr::desc(abs(Total)))
    num <- vapply(out, is.numeric, TRUE)
    out[num] <- lapply(out[num], function(x) round(x, 4))
    as.data.frame(out, stringsAsFactors = FALSE)
  }
  
  
  path_effects <- reactive({
    fit <- .phen_fit_active(); req(fit)
    df  <- path_numeric_df(); req(df, input$path_dep, input$path_ind)
    .build_path_effect_table(fit, df, y = input$path_dep, X = input$path_ind)
  })
  
  # ensure lme4 is available
  observe({
    if (!requireNamespace("lme4", quietly = TRUE)) {
      safe_message("Package 'lme4' not installed. Install it to enable genotypic path.", "error")
    }
  })
  
  # Build a BLUP matrix: rows = genotypes, cols = traits (Y and Xs)
  # --- Diagnostics + robust BLUP/mean extractor ---
  geno_blup_matrix <- reactive({
    req(raw_df(), input$id_col, input$path_dep, input$path_ind)
    validate(need(requireNamespace("lme4", quietly = TRUE),
                  "lme4 is required for BLUP estimation"))
    
    df  <- raw_df()
    idc <- input$id_col
    enc <- if (!is.null(input$path_env_ui) && !is.null(input$env_col) && nzchar(input$env_col)) input$env_col else NULL
    rpc <- if (!is.null(input$path_rep_ui) && !is.null(input$rep_col) && nzchar(input$rep_col)) input$rep_col else NULL
    
    # Coerce ID/ENV/REP
    if (!(idc %in% names(df))) validate(need(FALSE, "Genotype ID column not found in data."))
    df[[idc]] <- as.factor(df[[idc]])
    if (!is.null(enc) && enc %in% names(df)) df[[enc]] <- as.factor(df[[enc]]) else enc <- NULL
    if (!is.null(rpc) && rpc %in% names(df)) df[[rpc]] <- as.factor(df[[rpc]]) else rpc <- NULL
    
    # Traits = Y + X
    traits <- unique(c(input$path_dep, input$path_ind))
    miss   <- setdiff(traits, names(df))
    validate(need(length(miss) == 0,
                  paste("Missing trait columns:", paste(miss, collapse = ", "))))
    
    # Keep only numeric traits; silently coerce numeric-looking text
    for (nm in traits) {
      if (!is.numeric(df[[nm]])) {
        x <- gsub(",", "", as.character(df[[nm]]))
        if (all(grepl("^\\s*-?\\d*\\.?\\d*\\s*$", x) | x == "")) {
          suppressWarnings(df[[nm]] <- as.numeric(x))
        }
      }
    }
    not_num <- traits[!vapply(df[traits], is.numeric, TRUE)]
    validate(need(length(not_num) == 0,
                  paste("Non-numeric trait(s):", paste(not_num, collapse = ", "))))
    
    # Random structure: include terms only if >1 level present
    rand_terms <- c()
    if (nlevels(df[[idc]]) > 1) rand_terms <- c(rand_terms, paste0("(1|", idc, ")"))
    if (!is.null(enc) && nlevels(df[[enc]]) > 1) rand_terms <- c(rand_terms, paste0("(1|", enc, ")"))
    if (!is.null(rpc) && nlevels(df[[rpc]]) > 1) rand_terms <- c(rand_terms, paste0("(1|", rpc, ")"))
    validate(need(length(rand_terms) >= 1, "Random-effect structure invalid: need at least Genotype with >1 level."))
    rand_rhs <- paste(rand_terms, collapse = " + ")
    
    G_levels <- levels(df[[idc]])
    B_blup   <- matrix(NA_real_, nrow = length(G_levels), ncol = length(traits),
                       dimnames = list(G_levels, traits))
    B_mean   <- matrix(NA_real_, nrow = length(G_levels), ncol = length(traits),
                       dimnames = list(G_levels, traits))
    
    diag_rows <- list()
    
    for (tr in traits) {
      tr_vec <- df[[tr]]
      # Quick checks
      n_obs   <- sum(!is.na(tr_vec))
      n_genos <- nlevels(df[[idc]])
      # Genotype means fallback
      gm <- tapply(tr_vec, df[[idc]], function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
      B_mean[names(gm), tr] <- as.numeric(gm)
      
      # If near-zero genotype variance (by means), BLUPs won’t help
      if (stats::sd(B_mean[, tr], na.rm = TRUE) < .Machine$double.eps) {
        diag_rows[[tr]] <- data.frame(Trait = tr, Nobs = n_obs, Ngen = n_genos,
                                      Used = "MEAN", Reason = "Zero genotypic variance (means)")
        next
      }
      
      # Try BLUP
      fml_txt <- paste0(tr, " ~ 1 + ", rand_rhs)
      fit <- try(lme4::lmer(as.formula(fml_txt), data = df, REML = TRUE,
                            control = lme4::lmerControl(check.rankX = "ignore",
                                                        check.conv.singular = "ignore",
                                                        check.conv.hess = "ignore")), silent = TRUE)
      if (inherits(fit, "try-error")) {
        diag_rows[[tr]] <- data.frame(Trait = tr, Nobs = n_obs, Ngen = n_genos,
                                      Used = "MEAN", Reason = "lmer failed")
        next
      }
      
      re <- try(lme4::ranef(fit, condVar = FALSE)[[idc]][, 1, drop = TRUE], silent = TRUE)
      if (inherits(re, "try-error") || is.null(re)) {
        diag_rows[[tr]] <- data.frame(Trait = tr, Nobs = n_obs, Ngen = n_genos,
                                      Used = "MEAN", Reason = "BLUP extract failed")
        next
      }
      
      # Align to all genotypes
      B_blup[names(re), tr] <- as.numeric(re)
      diag_rows[[tr]] <- data.frame(Trait = tr, Nobs = n_obs, Ngen = n_genos,
                                    Used = "BLUP", Reason = "")
    }
    
    # Prefer BLUP; fall back to mean when BLUP is NA
    B <- B_blup
    na_pos <- is.na(B)
    if (any(na_pos)) B[na_pos] <- B_mean[na_pos]
    
    # Drop traits with all-NA columns
    keep_cols <- colSums(!is.na(B)) > 1
    if (any(!keep_cols)) {
      dropped <- names(keep_cols)[!keep_cols]
      if (length(dropped))
        safe_message(paste("Dropping traits with no usable genotype summaries:", paste(dropped, collapse = ", ")), "warning")
      B <- B[, keep_cols, drop = FALSE]
    }
    
    validate(need(ncol(B) >= 2, "Not enough traits with BLUPs/means to form genotypic correlations."))
    
    # Save diagnostics table for UI
    blup_diag <- do.call(rbind, diag_rows)
    blup_diag <- blup_diag[match(colnames(B_blup), blup_diag$Trait, nomatch = 0), , drop = FALSE]
    attr(B, "diag") <- blup_diag
    
    B
  })
  
  # Genotypic covariance/correlation from BLUPs
  geno_cov_cor <- reactive({
    B <- geno_blup_matrix(); req(B)
    # By default: correlations for SEM; keep cov too for downloads if needed
    covB <- stats::cov(B, use = "pairwise.complete.obs")
    corB <- stats::cov2cor(covB)
    list(cov = covB, cor = corB, n_geno = nrow(B),diag = attr(B, "diag"))
  })
  
  # Fit genotypic SEM using the BLUP-based correlation matrix
  geno_fit <- reactive({
    req(path_model_text(), geno_cov_cor(), input$path_dep, input$path_ind)
    G <- geno_cov_cor()
    C <- G$cor
    
    vars <- unique(c(input$path_dep, input$path_ind))
    present <- intersect(vars, intersect(rownames(C), colnames(C)))
    validate(need(length(present) == length(vars),
                  paste0("BLUP correlation missing variables: ",
                         paste(setdiff(vars, present), collapse = ", "))))
    C2 <- C[vars, vars, drop = FALSE]
    
    # near-PD fix (rare with BLUP-based cor, but safe)
    .make_pd_cor <- function(M) {
      M[is.na(M)] <- 0
      M <- (M + t(M)) / 2
      diag(M) <- 1
      ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values, error=function(e) NA)
      if (length(ev) && min(ev, na.rm = TRUE) > 1e-8) return(M)
      if (requireNamespace("Matrix", quietly = TRUE)) {
        return(as.matrix(Matrix::nearPD(M, corr = TRUE)$mat))
      }
      M
    }
    C2 <- .make_pd_cor(C2)
    
    model_text <- path_model_text()
    
    # Optional compact model on genotypic std effects
    if (isTRUE(input$path_compact)) {
      fit_rank <- tryCatch(
        lavaan::sem(model = model_text,
                    sample.cov = C2,
                    sample.nobs = G$n_geno,
                    std.lv = TRUE,
                    sample.cov.rescale = FALSE),
        error = function(e) { safe_message(paste("Genotypic rank fit error:", e$message), "error"); NULL }
      )
      req(fit_rank)
      std <- lavaan::standardizedSolution(fit_rank)
      std_y <- subset(std, op == "~" & lhs == input$path_dep & rhs %in% input$path_ind)
      validate(need(nrow(std_y) > 0, "No standardized genotypic paths found for selected variables."))
      k <- min(input$path_top_k %||% 6, nrow(std_y))
      keep <- head(std_y[order(-abs(std_y$est.std)), "rhs", drop = TRUE], k)
      model_text <- paste0(input$path_dep, " ~ ", paste(keep, collapse = " + "))
    }
    
    tryCatch(
      lavaan::sem(model = model_text,
                  sample.cov = C2,
                  sample.nobs = G$n_geno,      # effective N = #genotypes with BLUPs
                  std.lv = TRUE,
                  sample.cov.rescale = FALSE),
      error = function(e) { safe_message(paste("Genotypic SEM error:", e$message), "error"); NULL }
    )
  })
  # Ensure a matrix with given row/col names; keep drop=FALSE semantics
  .ensure_named_matrix <- function(M, rows, cols) {
    # coerce to matrix
    M <- as.matrix(M)
    # fix dims if needed
    if (nrow(M) != length(rows) || ncol(M) != length(cols)) {
      # try to subset/reorder if names exist
      rnames <- rownames(M); cnames <- colnames(M)
      if (!is.null(rnames) && !is.null(cnames)) {
        M <- M[rows, cols, drop = FALSE]
      } else {
        # rebuild a matrix of proper size
        M2 <- matrix(NA_real_, nrow = length(rows), ncol = length(cols))
        rownames(M2) <- rows; colnames(M2) <- cols
        # best effort copy if shapes match
        if (length(M) == length(M2)) M2[] <- as.numeric(M)
        M <- M2
      }
    }
    # ensure dimnames
    if (is.null(rownames(M))) rownames(M) <- rows
    if (is.null(colnames(M))) colnames(M) <- cols
    M
  }
  
  # Safe 1x1 correlation from covariance matrix
  .cov2cor_safe <- function(S) {
    S <- as.matrix(S)
    d <- sqrt(diag(S))
    if (length(d) == 1L) {
      return(matrix(1, nrow = 1, ncol = 1,
                    dimnames = list(rownames(S), colnames(S))))
    }
    Dinv <- diag(1/d)
    R <- Dinv %*% S %*% Dinv
    # carry names
    rownames(R) <- rownames(S); colnames(R) <- colnames(S)
    R
  }
  
  
  # helper to ensure named submatrix even for 1x1
  .ensure_named_matrix <- function(M, rows, cols) {
    M <- as.matrix(M)
    if (is.null(rownames(M))) rownames(M) <- seq_len(nrow(M))
    if (is.null(colnames(M))) colnames(M) <- seq_len(ncol(M))
    M2 <- M[rows, cols, drop = FALSE]
    if (is.null(rownames(M2))) rownames(M2) <- rows
    if (is.null(colnames(M2))) colnames(M2) <- cols
    M2
  }
  
  geno_path_effects <- reactive({
    fit <- geno_fit(); req(fit)
    
    # DIRECT (standardized)
    std <- lavaan::standardizedSolution(fit)
    dir_df <- subset(std, op == "~" & lhs == input$path_dep & rhs %in% input$path_ind,
                     select = c("rhs","est.std","pvalue"))
    if (is.null(dir_df) || nrow(dir_df) == 0) {
      return(data.frame(
        Predictor = character(), Direct = numeric(), Indirect = numeric(),
        Total = numeric(), SE = numeric(), z = numeric(),
        p_direct = numeric(), Signif = character(), stringsAsFactors = FALSE
      ))
    }
    names(dir_df) <- c("Predictor","Direct","p_direct")
    
    # SE and z
    pe <- lavaan::parameterEstimates(fit)
    se_z <- subset(pe, op == "~" & lhs == input$path_dep & rhs %in% input$path_ind,
                   select = c("rhs","se","z","pvalue"))
    if (is.null(se_z) || nrow(se_z) == 0) {
      se_z <- data.frame(Predictor = dir_df$Predictor,
                         SE = NA_real_, z = NA_real_, p_raw = NA_real_)
    } else {
      names(se_z) <- c("Predictor","SE","z","p_raw")
    }
    
    # Implied Sigma → correlation among variables (with names)
    Sigma_hat <- lavaan::inspect(fit, "sigma")
    vars  <- unique(c(input$path_dep, input$path_ind))
    # ensure named submatrix even for 1x1
    Sigma_sub <- .ensure_named_matrix(Sigma_hat, rows = vars, cols = vars)
    
    D <- diag(1 / sqrt(diag(Sigma_sub)))
    D[!is.finite(D)] <- 0
    Rall <- as.matrix(D %*% Sigma_sub %*% D)
    rownames(Rall) <- rownames(Sigma_sub); colnames(Rall) <- colnames(Sigma_sub)
    
    preds <- input$path_ind
    if (length(preds) == 0) {
      return(data.frame(
        Predictor = character(), Direct = numeric(), Indirect = numeric(),
        Total = numeric(), SE = numeric(), z = numeric(),
        p_direct = numeric(), Signif = character(), stringsAsFactors = FALSE
      ))
    }
    Rmat <- .ensure_named_matrix(Rall, rows = preds, cols = preds)
    
    dir_vec <- dir_df$Direct; names(dir_vec) <- dir_df$Predictor
    indir <- sapply(dir_df$Predictor, function(i) {
      if (!i %in% colnames(Rmat)) return(0)
      others <- setdiff(colnames(Rmat), i)
      if (length(others) == 0) return(0)
      sum(Rmat[i, others, drop = TRUE] * dir_vec[others], na.rm = TRUE)
    })
    
    out <- dplyr::left_join(dir_df, se_z, by = "Predictor")
    out$Indirect <- indir[match(out$Predictor, names(indir))]
    out$Total    <- out$Direct + out$Indirect
    out$Signif   <- .path_sig_stars(out$p_direct)
    
    out <- dplyr::arrange(out, dplyr::desc(abs(Total)))
    num <- vapply(out, is.numeric, TRUE)
    out[num] <- lapply(out[num], function(x) round(x, 4))
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    rownames(out) <- NULL
    out
  })
  
  # ---------- outputs ----------
  output$phen_path_plot <- renderPlot({
    fit <- .phen_fit_active(); req(fit)
    semPlot::semPaths(
      fit, whatLabels="est", layout="tree", style="ram", rotation=2,
      residuals=isTRUE(input$path_show_resid), residScale=4, covAtResiduals=FALSE,
      edge.label.cex=0.7, thresholds=FALSE, sizeMan=5, sizeLat=5, curve=2.2, curvature=2.0,
      reorder=TRUE, optimizeLatRes=TRUE, mar=c(5,20,18,5), intercepts=FALSE, edge.color="black"
    )
  })
  
  output$geno_path_plot <- renderPlot({
    fit <- geno_fit(); req(fit)
    semPlot::semPaths(
      fit, whatLabels="est", layout="tree", style="ram", rotation=2,
      residuals=isTRUE(input$path_show_resid), residScale=4, covAtResiduals=FALSE,
      edge.label.cex=0.9, thresholds=FALSE, sizeMan=5, sizeLat=5, curve=2.2, curvature=2.0,
      reorder=TRUE, optimizeLatRes=TRUE, mar=c(5,20,18,5), intercepts=FALSE, edge.color="red"
    )
  })
  
  output$phen_path_table <- DT::renderDT({
    tbl <- path_effects(); req(nrow(tbl) > 0)
    DT::datatable(
      tbl[, c("Predictor","Direct","Indirect","Total","SE","z","p_direct","Signif")],
      options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")),
      extensions = "Buttons", rownames = FALSE
    )
  })
  output$geno_path_table <- DT::renderDT({
    tbl <- geno_path_effects(); req(!is.null(tbl))
    # Ensure data.frame
    tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
    
    wanted <- c("Predictor","Direct","Indirect","Total","SE","z","p_direct","Signif")
    have   <- intersect(wanted, colnames(tbl))
    
    # If nothing to show, render an empty but named table
    if (length(have) == 0 || nrow(tbl) == 0) {
      empty <- as.data.frame(setNames(replicate(length(wanted), logical(0), simplify = FALSE), wanted))
      return(DT::datatable(
        empty,
        options = list(pageLength = 10, searching = FALSE),
        rownames = FALSE
      ))
    }
    
    DT::datatable(
      tbl[, have, drop = FALSE],
      options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")),
      extensions = "Buttons", rownames = FALSE
    )
  })
  
  
  
  output$phen_path_plot_download <- downloadHandler(
    filename = function() paste0("phen_path_", Sys.Date(), ".png"),
    content  = function(file) {
      fit <- .phen_fit_active(); req(fit)
      png(file, width = 1800, height = 1100, res = 150)
      semPlot::semPaths(
        fit, whatLabels="est", layout="tree", style="ram", rotation=2,
        residuals=isTRUE(input$path_show_resid), residScale=4, covAtResiduals=FALSE,
        edge.label.cex=0.9, thresholds=FALSE, sizeMan=5, sizeLat=5, curve=2.2, curvature=2.0,
        reorder=TRUE, optimizeLatRes=TRUE, mar=c(5,18,12,5), intercepts=FALSE, edge.color="black"
      )
      dev.off()
    }
  )
  
  output$geno_path_plot_download <- downloadHandler(
    filename = function() paste0("geno_path_", Sys.Date(), ".png"),
    content  = function(file) {
      fit <- geno_fit(); req(fit)
      png(file, width = 1800, height = 1100, res = 150)
      semPlot::semPaths(
        fit, whatLabels="est", layout="tree", style="ram", rotation=2,
        residuals=isTRUE(input$path_show_resid), residScale=4, covAtResiduals=FALSE,
        edge.label.cex=0.9, thresholds=FALSE, sizeMan=5, sizeLat=5, curve=2.2, curvature=2.0,
        reorder=TRUE, optimizeLatRes=TRUE, mar=c(5,18,12,5), intercepts=FALSE, edge.color="red"
      )
      dev.off()
    }
  )
  
  output$phen_path_table_download <- downloadHandler(
    filename = function() paste0("phen_path_coefficients_", Sys.Date(), ".csv"),
    content  = function(file) {
      write.csv(path_effects(), file, row.names = FALSE, na = "")
    }
  )
  
  output$geno_path_table_download <- downloadHandler(
    filename = function() paste0("geno_path_coefficients_", Sys.Date(), ".csv"),
    content  = function(file) {
      write.csv(geno_path_effects(), file, row.names = FALSE, na = "")
    }
  )
  
  # ===== AMMI  =====
  ammi_raw <- reactive({
    req(input$ammi_file)
    df <- read_data(input$ammi_file)
    names(df) <- iconv(names(df), from = "UTF-8", to = "ASCII//TRANSLIT")
    df
  })
  output$ammi_mapper_ui <- renderUI({
    req(ammi_raw()); tagList(tags$strong("Column Mapping"),
                             role_picker_ui("ammi", ammi_raw(), c("environment","genotype","replication")))
  })
  ammi_map <- reactive({
    df <- ammi_raw(); req(df)
    det <- detect_roles(df, c("environment","genotype","replication"))
    list(
      environment = input$ammi_environment %||% det$environment,
      genotype    = input$ammi_genotype    %||% det$genotype,
      replication = input$ammi_replication %||% det$replication
    )
  })
  ammi_df <- reactive({ canonicalize(ammi_raw(), ammi_map()) })
  output$ammi_trait_ui <- renderUI({
    req(ammi_df()); nums <- names(ammi_df())[vapply(ammi_df(), is.numeric, TRUE)]
    validate(need(length(nums) > 0, "No numeric trait found."))
    selectInput("ammi_trait", "Trait (numeric)", choices = nums)
  })
  # Helper: drop levels that would cause contrasts errors
  .ammi_prepare_raw <- function(df, trait) {
    validate(need(all(c("E","G") %in% names(df)), "Map Environment and Genotype."))
    
    # Ensure factors + drop empties
    df$E <- droplevels(factor(df$E))
    df$G <- droplevels(factor(df$G))
    if (!"R" %in% names(df)) df$R <- 1L
    df$R <- droplevels(factor(df$R))
    
    # Keep rows with finite response
    df <- df[is.finite(df[[trait]]), , drop = FALSE]
    
    # Remove environments with <2 genotypes and genotypes with <2 environments
    # (these cause contrasts errors inside the internal ANOVA)
    keepE <- names(which(tapply(df$G, df$E, function(x) length(unique(x)) >= 2)))
    df <- df[df$E %in% keepE, , drop = FALSE]
    df$E <- droplevels(df$E)
    
    keepG <- names(which(tapply(df$E, df$G, function(x) length(unique(x)) >= 2)))
    df <- df[df$G %in% keepG, , drop = FALSE]
    df$G <- droplevels(df$G)
    
    validate(need(nlevels(df$E) >= 2,
                  "AMMI needs ≥2 environments after filtering (removed envs with <2 genotypes)."))
    validate(need(nlevels(df$G) >= 2,
                  "AMMI needs ≥2 genotypes after filtering (removed gens with <2 environments)."))
    
    # Also require ≥2 rep levels overall (not strictly necessary, but stabilizes SS partition)
    if (nlevels(df$R) < 2) {
      # If there is a replication column but it collapsed, synthesize one by numbering within E
      df <- df |>
        dplyr::group_by(E) |>
        dplyr::mutate(R = factor(dplyr::row_number())) |>
        dplyr::ungroup()
    }
    
    # Final variance guard
    validate(need(stats::var(df[[trait]], na.rm = TRUE) > 0,
                  "Trait variance is zero after filtering; AMMI cannot be fitted."))
    
    df
  }
  ammi_model <- reactive({
    df <- ammi_df(); req(df, input$ammi_trait)
    dat <- .ammi_prepare_raw(df, input$ammi_trait)
    
    # Fit AMMI on RAW with replication
    metan::performs_ammi(
      .data = dat,
      env   = !!rlang::sym("E"),
      gen   = !!rlang::sym("G"),
      rep   = !!rlang::sym("R"),
      resp  = !!rlang::sym(input$ammi_trait),
      verbose = FALSE
    )
  })
  get_trait_obj <- reactive({
    mod <- ammi_model()
    tr  <- input$ammi_trait
    if (tr %in% names(mod)) mod[[tr]] else mod
  })
  # How many PCs are conceivable from design?
  .ammi_np_from_design <- function(df) {
    # min(#G-1, #E-1)
    g <- nlevels(df$G); e <- nlevels(df$E)
    max(0, min(g - 1, e - 1))
  }
  output$ammi_table <- DT::renderDT({
    res <- get_trait_obj()
    anova_tbl <- tryCatch(res$ANOVA, error = function(...) NULL)
    if (is.list(anova_tbl) && !is.data.frame(anova_tbl)) anova_tbl <- anova_tbl[[1]]
    validate(need(!is.null(anova_tbl) && nrow(anova_tbl) > 0, "ANOVA table not available"))
    
    tbl <- .add_sig_stars(anova_tbl)
    
    DT::datatable(
      tbl,
      options = list(pageLength = 12, dom = "Bfrtip", buttons = c("copy","csv","excel")),
      extensions = "Buttons",
      rownames = FALSE
    )
  })
  
  output$ammi_plot <- renderPlot({
    mod <- ammi_model(); req(mod)
    trait_label <- input$ammi_trait
    df_fitted <- ammi_df()
    df_fitted$E <- droplevels(factor(df_fitted$E))
    df_fitted$G <- droplevels(factor(df_fitted$G))
    npc <- .ammi_np_from_design(df_fitted)
    plot_types <- as.numeric(input$ammi_plot_type)
    if (npc < 2) plot_types <- intersect(plot_types, c(1,4))
    validate(need(length(plot_types) > 0, "Only IPCA1 is supported by this design; choose AMMI1 or Environment Scores."))
    plots <- list()
    for (pt in plot_types) {
      if (pt == 1) plots[[length(plots)+1]] <- metan::plot_scores(mod, type = 1, main = paste0("AMMI1 - ", trait_label), max.overlaps = 100)
      if (pt == 2 && npc >= 2) plots[[length(plots)+1]] <- metan::plot_scores(mod, type = 2, main = paste0("AMMI2 - ", trait_label), max.overlaps = 100)
      if (pt == 4) plots[[length(plots)+1]] <- metan::plot_scores(mod, type = 4, main = paste0("Env Scores - ", trait_label), max.overlaps = 100)
      if (pt == 6 && npc >= 2) plots[[length(plots)+1]] <- metan::plot_scores(mod, type = 2, polygon = TRUE, axis.expand = 1.8, main = paste0("AMMI2 Polygon - ", trait_label), max.overlaps = 100)
    }
    do.call(gridExtra::grid.arrange, c(plots, ncol = min(2, length(plots))))
  }, height = 600)
  output$ammi_plot_download <- downloadHandler(
    filename = function() paste0("AMMI_biplots_", input$ammi_trait, ".png"),
    content  = function(file) {
      png(file, 1600, 800, res = 120)
      
      mod <- ammi_model()
      if (is.null(mod)) { plot.new(); text(.5, .5, "AMMI model unavailable."); dev.off(); return() }
      
      trait_label <- input$ammi_trait
      df_fitted   <- ammi_df()
      df_fitted$E <- droplevels(factor(df_fitted$E))
      df_fitted$G <- droplevels(factor(df_fitted$G))
      
      # max PCs supported by design: min(#G-1, #E-1)
      g   <- nlevels(df_fitted$G)
      e   <- nlevels(df_fitted$E)
      npc <- max(0, min(g - 1, e - 1))
      
      plot_types <- as.numeric(input$ammi_plot_type)
      if (npc < 2) plot_types <- intersect(plot_types, c(1, 4))
      if (!length(plot_types)) { plot.new(); text(.5, .5, "Only IPCA1 supported by this design"); dev.off(); return() }
      
      plots <- list()
      for (pt in plot_types) {
        if (pt == 1) plots[[length(plots)+1]] <- metan::plot_scores(mod, type = 1, main = paste0("AMMI1 - ", trait_label), max.overlaps = 100)
        if (pt == 2 && npc >= 2) plots[[length(plots)+1]] <- metan::plot_scores(mod, type = 2, main = paste0("AMMI2 - ", trait_label), max.overlaps = 100)
        if (pt == 4) plots[[length(plots)+1]] <- metan::plot_scores(mod, type = 4, main = paste0("Env Scores - ", trait_label), max.overlaps = 100)
        if (pt == 6 && npc >= 2) plots[[length(plots)+1]] <- metan::plot_scores(mod, type = 2, polygon = TRUE, axis.expand = 1.8, main = paste0("AMMI2 Polygon - ", trait_label), max.overlaps = 100)
      }
      gridExtra::grid.arrange(grobs = plots, ncol = min(2, length(plots)))
      dev.off()
    }
  )
  
  output$ammi_table_download <- downloadHandler(
    filename = function() paste0("AMMI_ANOVA_", input$ammi_trait, ".csv"),
    content  = function(file) {
      res <- get_trait_obj()
      anova_tbl <- tryCatch(res$ANOVA, error = function(...) NULL)
      if (is.list(anova_tbl) && !is.data.frame(anova_tbl)) anova_tbl <- anova_tbl[[1]]
      tbl <- .add_sig_stars(anova_tbl)
      write.csv(tbl, file, row.names = FALSE, na = "")
    }
  )
  output$geno_diag_table <- DT::renderDT({
    G <- geno_cov_cor(); req(G$diag)
    DT::datatable(
      G$diag,
      options = list(pageLength = 10, searching = FALSE),
      rownames = FALSE
    )
  })
  
  
  # ===== GGE  ==============================
  .gge_extract_xy <- function(p) {
    out <- NULL
    try({
      gb <- ggplot2::ggplot_build(p)
      parts <- lapply(gb$data, function(d) {
        keep <- intersect(c("label", "x", "y"), names(d))
        if (length(keep) >= 2) as.data.frame(d[keep]) else NULL
      })
      out <- do.call(rbind, Filter(Negate(is.null), parts))
      if (!is.null(out) && !("label" %in% names(out)) && ("x" %in% names(out))) {
        if (!is.null(p$data) && "label" %in% names(p$data)) {
          out$label <- p$data$label[seq_len(nrow(out))]
        }
      }
      if (!is.null(out)) {
        out <- unique(out)
        rownames(out) <- NULL
      }
    }, silent = TRUE)
    out
  }
  
  # Safe ranking helpers
  .rank_desc <- function(x) rank(-x, ties.method = "min")
  .rank_asc  <- function(x) rank( x, ties.method = "min")
  
  gge_raw <- reactive({
    req(input$gge_file)
    df <- read_data(input$gge_file)
    df <- df[!is.na(df[[1]]) & df[[1]] != "", ]
    if ("Enviorment" %in% names(df)) names(df)[names(df) == "Enviorment"] <- "Environment"
    df
  })
  
  output$gge_mapper_ui <- renderUI({
    req(gge_raw())
    tagList(
      tags$strong("Column Mapping"),
      role_picker_ui("gge", gge_raw(), c("environment", "genotype"))
    )
  })
  
  gge_map <- reactive({
    df  <- gge_raw(); req(df)
    det <- detect_roles(df, c("environment", "genotype"))
    list(
      environment = input$gge_environment %||% det$environment,
      genotype    = input$gge_genotype    %||% det$genotype
    )
  })
  
  gge_df <- reactive({
    canonicalize(gge_raw(), gge_map())
  })
  
  output$gge_trait_ui <- renderUI({
    df <- gge_df(); req(df)
    nums <- names(df)[vapply(df, is.numeric, TRUE)]
    validate(need(length(nums) > 0, "No numeric trait found."))
    selectInput("gge_trait", "Trait column (numeric)", choices = nums)
  })
  
  output$gge_plot_type_ui <- renderUI({
    selectInput(
      "gge_plot_type", "GGE Plot Type",
      choices = c("Which-Won-Where" = "ww",
                  "Mean vs Stability" = "ms",
                  "Genotype Ranking"  = "Gen",
                  "Environment Ranking" = "Env",
                  "Discrimination vs. representativeness" = "DiscRep")
    )
  })
  
  output$gge_table <- renderDT({
    datatable(gge_df())
  })
  
  
  build_gge <- function() {
    df <- gge_df(); req(df, input$gge_trait)
    validate(need(all(c("E", "G") %in% names(df)), "Please map Environment and Genotype."))
    
    wide_prep <- df %>%
      dplyr::select(E, G, !!rlang::sym(input$gge_trait)) %>%
      dplyr::filter(!is.na(E) & !is.na(G) & !is.na(.data[[input$gge_trait]])) %>%
      dplyr::group_by(E, G) %>%
      dplyr::summarise(resp = mean(.data[[input$gge_trait]], na.rm = TRUE), .groups = "drop")
    
    validate(need(nrow(wide_prep) > 0, "No valid data after removing missing values."))
    
    wide <- tidyr::pivot_wider(wide_prep, names_from = E, values_from = resp, values_fill = NA)
    mat  <- as.data.frame(wide); rownames(mat) <- mat$G; mat$G <- NULL; mat <- as.matrix(mat)
    
    # keep only rows/cols with some data
    mat <- mat[rowSums(!is.na(mat)) > 0, colSums(!is.na(mat)) > 0, drop = FALSE]
    validate(need(nrow(mat) >= 2 && ncol(mat) >= 2,
                  "Need ≥2 genotypes and ≥2 environments after filtering."))
    
    GGEBiplots::GGEModel(mat)
  }
  
  output$gge_plot <- renderPlot({
    gge_mod <- tryCatch(build_gge(), error = function(e) { safe_message(paste("GGE error:", e$message)); NULL })
    req(gge_mod)
    switch(
      input$gge_plot_type,
      "ww"      = print(GGEBiplots::WhichWon(gge_mod)),
      "ms"      = print(GGEBiplots::MeanStability(gge_mod)),
      "Gen"     = print(GGEBiplots::RankGen(gge_mod)),
      "Env"     = print(GGEBiplots::RankEnv(gge_mod)),
      "DiscRep" = print(GGEBiplots::DiscRep(gge_mod))
    )
  })
  
  output$gge_plot_download <- downloadHandler(
    filename = function() "gge_plot.png",
    content  = function(file) {
      png(file, 1000, 800)
      gge_mod <- tryCatch(build_gge(), error = function(e) NULL)
      if (is.null(gge_mod)) {
        plot.new(); text(0.5, 0.5, "GGE failed - check data")
      } else {
        switch(
          input$gge_plot_type,
          "ww"      = print(GGEBiplots::WhichWon(gge_mod)),
          "ms"      = print(GGEBiplots::MeanStability(gge_mod)),
          "Gen"     = print(GGEBiplots::RankGen(gge_mod)),
          "Env"     = print(GGEBiplots::RankEnv(gge_mod)),
          "DiscRep" = print(GGEBiplots::DiscRep(gge_mod))
        )
      }
      dev.off()
    }
  )
  
  output$gge_table_download <- downloadHandler(
    filename = function() "gge_table.csv",
    content  = function(file) write.csv(gge_df(), file, row.names = FALSE)
  )
  
  
  gge_tables <- reactive({
    gge_mod <- tryCatch(build_gge(), error = function(e) NULL); req(gge_mod)
    
    # Try to extract directly from the plotted ggplots
    p_rg <- tryCatch(GGEBiplots::RankGen(gge_mod),        error = function(e) NULL)
    p_re <- tryCatch(GGEBiplots::RankEnv(gge_mod),        error = function(e) NULL)
    p_ms <- tryCatch(GGEBiplots::MeanStability(gge_mod),  error = function(e) NULL)
    p_dr <- tryCatch(GGEBiplots::DiscRep(gge_mod),        error = function(e) NULL)
    
    ext_rg <- if (!is.null(p_rg)) .gge_extract_xy(p_rg) else NULL
    ext_re <- if (!is.null(p_re)) .gge_extract_xy(p_re) else NULL
    ext_ms <- if (!is.null(p_ms)) .gge_extract_xy(p_ms) else NULL
    ext_dr <- if (!is.null(p_dr)) .gge_extract_xy(p_dr) else NULL
    
    # Fallback computation if extraction fails
    df <- gge_df(); req(df, input$gge_trait)
    wide_prep <- df %>%
      dplyr::select(E, G, !!rlang::sym(input$gge_trait)) %>%
      dplyr::filter(!is.na(E) & !is.na(G) & !is.na(.data[[input$gge_trait]])) %>%
      dplyr::group_by(E, G) %>%
      dplyr::summarise(resp = mean(.data[[input$gge_trait]], na.rm = TRUE), .groups = "drop")
    wide <- tidyr::pivot_wider(wide_prep, names_from = E, values_from = resp, values_fill = NA)
    mat  <- as.data.frame(wide); rownames(mat) <- mat$G; mat$G <- NULL
    M    <- as.matrix(mat)
    
    approx_scores <- function(M) {
      M_c <- scale(M, center = TRUE, scale = FALSE)   # center by env means (GGE convention)
      pc  <- prcomp(M_c, center = FALSE, scale. = FALSE)
      list(
        gen    = data.frame(Genotype = rownames(M_c), PC1 = pc$x[,1], PC2 = pc$x[,2], row.names = NULL),
        env    = data.frame(Environment = colnames(M_c),
                            PC1 = pc$rotation[,1] * pc$sdev[1],
                            PC2 = pc$rotation[,2] * pc$sdev[2], row.names = NULL),
        mean_g = rowMeans(M, na.rm = TRUE)
      )
    }
    sc <- approx_scores(M)
    
    # RankGen
    rankgen_tbl <- if (!is.null(ext_rg) && "label" %in% names(ext_rg)) {
      d <- ext_rg; names(d)[names(d) == "label"] <- "Genotype"
      d$Rank <- .rank_desc(d$x)
      out <- dplyr::arrange(d, Rank) %>% dplyr::select(Rank, Genotype, x, y)
      names(out)[3:4] <- c("Axis1", "Axis2"); out
    } else {
      tmp <- data.frame(
        Genotype = sc$gen$Genotype,
        Mean     = as.numeric(sc$mean_g[sc$gen$Genotype]),
        PC1      = sc$gen$PC1,
        PC2      = sc$gen$PC2
      )
      tmp$Rank <- .rank_desc(tmp$Mean)
      tmp[order(tmp$Rank), c("Rank","Genotype","Mean","PC1","PC2")]
    }
    
    # RankEnv
    rankenv_tbl <- if (!is.null(ext_re) && "label" %in% names(ext_re)) {
      d <- ext_re; names(d)[names(d) == "label"] <- "Environment"
      d$Discriminativeness <- sqrt(d$x^2 + d$y^2)
      d$Rank <- .rank_desc(d$Discriminativeness)
      out <- dplyr::arrange(d, Rank) %>% dplyr::select(Rank, Environment, Discriminativeness, x, y)
      names(out)[4:5] <- c("Axis1","Axis2"); out
    } else {
      tmp <- data.frame(Environment = sc$env$Environment, PC1 = sc$env$PC1, PC2 = sc$env$PC2)
      tmp$Discriminativeness <- sqrt(tmp$PC1^2 + tmp$PC2^2)
      tmp$Rank <- .rank_desc(tmp$Discriminativeness)
      tmp[order(tmp$Rank), c("Rank","Environment","Discriminativeness","PC1","PC2")]
    }
    
    # Mean–Stability
    meanstab_tbl <- if (!is.null(ext_ms) && "label" %in% names(ext_ms)) {
      d <- ext_ms; names(d)[names(d) == "label"] <- "Genotype"
      d$Stability <- abs(d$y)
      d$RankMean  <- .rank_desc(d$x)
      d$RankStab  <- .rank_asc(d$Stability)
      out <- dplyr::arrange(d, RankMean) %>% dplyr::select(Genotype, x, Stability, RankMean, RankStab)
      names(out)[2] <- "MeanAxis"; out
    } else {
      tmp <- data.frame(Genotype = sc$gen$Genotype,
                        Mean     = as.numeric(sc$mean_g[sc$gen$Genotype]),
                        Stability = abs(sc$gen$PC2))
      tmp$RankMean <- .rank_desc(tmp$Mean)
      tmp$RankStab <- .rank_asc(tmp$Stability)
      tmp[order(tmp$RankMean), c("Genotype","Mean","Stability","RankMean","RankStab")]
    }
    
    # Discrimination–Representativeness
    discrep_tbl <- if (!is.null(ext_dr) && "label" %in% names(ext_dr)) {
      d <- ext_dr; names(d)[names(d) == "label"] <- "Environment"
      d$Discriminativeness <- sqrt(d$x^2 + d$y^2)
      d$Angle <- atan2(d$y, d$x)
      d$Representativeness <- 1 - abs(d$Angle) / pi
      out <- dplyr::arrange(d, dplyr::desc(Discriminativeness)) %>%
        dplyr::select(Environment, Discriminativeness, Representativeness, x, y)
      names(out)[4:5] <- c("Axis1","Axis2"); out
    } else {
      tmp <- data.frame(Environment = sc$env$Environment, PC1 = sc$env$PC1, PC2 = sc$env$PC2)
      tmp$Discriminativeness  <- sqrt(tmp$PC1^2 + tmp$PC2^2)
      tmp$Representativeness  <- abs(tmp$PC1) / sqrt(tmp$PC1^2 + tmp$PC2^2)
      out <- tmp[order(-tmp$Discriminativeness), c("Environment","Discriminativeness","Representativeness","PC1","PC2")]
      names(out)[4:5] <- c("Axis1","Axis2"); out
    }
    
    # Round numerics for display
    round_cols <- function(d) {
      if (is.null(d)) return(d)
      num_idx <- vapply(d, is.numeric, TRUE)
      d[num_idx] <- lapply(d[num_idx], function(x) round(x, 4))
      d
    }
    
    list(
      rank_gen  = round_cols(rankgen_tbl),
      rank_env  = round_cols(rankenv_tbl),
      mean_stab = round_cols(meanstab_tbl),
      disc_rep  = round_cols(discrep_tbl)
    )
  })
  
  
  output$gge_rankgen_tbl  <- DT::renderDT({
    DT::datatable(gge_tables()$rank_gen,
                  options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")),
                  extensions = "Buttons", rownames = FALSE)
  })
  output$gge_rankenv_tbl  <- DT::renderDT({
    DT::datatable(gge_tables()$rank_env,
                  options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")),
                  extensions = "Buttons", rownames = FALSE)
  })
  output$gge_meanstab_tbl <- DT::renderDT({
    DT::datatable(gge_tables()$mean_stab,
                  options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")),
                  extensions = "Buttons", rownames = FALSE)
  })
  output$gge_discrep_tbl  <- DT::renderDT({
    DT::datatable(gge_tables()$disc_rep,
                  options = list(pageLength = 25, dom = "Bfrtip", buttons = c("copy","csv","excel")),
                  extensions = "Buttons", rownames = FALSE)
  })
  
  output$gge_rankgen_csv  <- downloadHandler(
    filename = function() "GGE_RankGen.csv",
    content  = function(file) readr::write_csv(gge_tables()$rank_gen, file, na = "")
  )
  output$gge_rankenv_csv  <- downloadHandler(
    filename = function() "GGE_RankEnv.csv",
    content  = function(file) readr::write_csv(gge_tables()$rank_env, file, na = "")
  )
  output$gge_meanstab_csv <- downloadHandler(
    filename = function() "GGE_MeanStability.csv",
    content  = function(file) readr::write_csv(gge_tables()$mean_stab, file, na = "")
  )
  output$gge_discrep_csv  <- downloadHandler(
    filename = function() "GGE_DiscRep.csv",
    content  = function(file) readr::write_csv(gge_tables()$disc_rep, file, na = "")
  )
  
  # ===== Demo Data =====
  demo_data <- reactiveVal(NULL)
  output$block_size_ui <- renderUI({
    req(input$genotypes); t <- as.numeric(input$genotypes); divisors <- which(t %% 1:t == 0)
    valid_k  <- divisors[divisors > 1 & divisors < t]
    if (!length(valid_k)) return(helpText("No valid block sizes: choose a different number of genotypes."))
    selectInput("block_size", "Block size (k):", choices = valid_k, selected = valid_k[1])
  })
  observeEvent(input$generate_demo, {
    req(input$genotypes, input$replications, input$traits)
    genos <- paste0("G", 1:input$genotypes); data <- NULL
    if (input$design_type == "Augmented") { req(input$checks); checks <- paste0("C", 1:input$checks); test_genos <- genos }
    switch(input$design_type,
           "CRD" = { data <- design.crd(trt = genos, r = input$replications, seed = 123)$book },
           "RBD" = { data <- design.rcbd(trt = genos, r = input$replications, seed = 123)$book },
           "LSD" = {
             t <- as.numeric(input$genotypes)
             s <- floor(sqrt(t))
             if (s * s != t) { safe_message("LSD requires number of genotypes to be a perfect square (e.g., 9, 16, 25)."); return(NULL) }
             data <- design.lsd(trt = genos, seed = 123)$book
           },
           "Alpha Lattice" = {
             t  <- as.numeric(input$genotypes); r <- as.numeric(input$replications); k <- as.numeric(input$block_size)
             if (is.na(t) || is.na(r) || is.na(k)) { safe_message("Alpha lattice inputs must be numeric."); return(NULL) }
             if (r < 2) { safe_message("Alpha lattice requires at least 2 replications."); return(NULL) }
             if (t %% k != 0) { safe_message("Genotypes (t) must be divisible by block size (k)."); return(NULL) }
             alpha_res <- tryCatch(agricolae::design.alpha(trt = genos, r = r, k = k, seed = 123), error = function(e) e)
             if (inherits(alpha_res, "error")) { safe_message(paste0("Alpha lattice failed: ", alpha_res$message)); return(NULL) }
             data <- alpha_res[["book"]]
           },
           "Split Plot" = {
             if (input$genotypes < 4) { safe_message("Split Plot requires at least 4 genotypes."); return(NULL) }
             main_plots <- genos[1:ceiling(input$genotypes/2)]; sub_plots  <- genos
             data <- design.split(trt1 = main_plots, trt2 = sub_plots, r = input$replications, design = "rcbd", seed = 123)$book
           },
           "Augmented" = { data <- design.dau(trt1 = checks, trt2 = test_genos, r = input$replications, seed = 123)$book }
    )
    for (i in seq_len(input$traits)) data[[paste0("Trait", i)]] <- round(rnorm(nrow(data), mean = 50 + i*5, sd = 10), 2)
    demo_data(data)
  })
  output$demo_preview <- DT::renderDT({ req(demo_data()); demo_data() })
  output$download_demo <- downloadHandler(
    filename = function() paste0(input$design_type, "_DemoData.csv"),
    content  = function(file) write.csv(demo_data(), file, row.names = FALSE)
  )
}
 
# ---- Run ----
shinyApp(ui = ui, server = server)
