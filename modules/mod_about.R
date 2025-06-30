aboutUI <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    title = "About",
    shiny::withMathJax(),
    
    shiny::tags$div(
      class = "container",
      
      # ------------------------------ 1. Overview ---------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Overview")),
        shiny::tags$p("Analysis of T‑cell activation‑induced marker (AIM) assay data requires normalization of AIM\u2009+ cell frequencies to ‘background’ AIM\u2009+ frequencies in an unstimulated control. Subtracting or dividing by the unstimulated control each have specific disadvantages and can amplify technical variability in the assay."),
        shiny::tags$p("The Box‑Cox correction, based on the transformation proposed by Box & Cox in 1964\u00b9, combines features of division and linear subtraction to better match the mathematical properties of AIM datasets, reducing technical variability and improving signal detection.")
      ),
      
      # ------------------ 2. Application of the Box‑Cox method -------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Application of the Box‑Cox method to AIM assay data")),
        shiny::tags$ol(
          shiny::tags$li(htmltools::HTML("Apply the transformation: $$\\text{BoxCox}(x,\\;\\lambda \\in [0,1]) = \\begin{cases} \\dfrac{x^{\\lambda}-1}{\\lambda}, & \\lambda \\neq 0 \\; ; \\; \\ln x, & \\lambda = 0 \\end{cases}$$")),
          shiny::tags$li(htmltools::HTML("Net stimulated value: $$\\beta = \\text{BoxCox}(\\%\\text{AIM}^{+}_{\\text{stim}}) - \\text{BoxCox}(\\%\\text{AIM}^{+}_{\\text{unstim}})$$")),
          shiny::tags$li(htmltools::HTML("Return to original scale: $$\\text{SI} = \\text{invBoxCox}(\\beta,\\;\\lambda) = \\begin{cases} (\\lambda\\beta + 1)^{1/\\lambda}, & \\lambda \\neq 0 ;\\; e^{\\beta}, & \\lambda = 0 \\end{cases}$$"))
        ),
        shiny::tags$p("Choose λ (0–1) carefully, considering analysis goals and dataset properties.")
      ),
      
      # ------------------------ 3. Advanced use -----------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Advanced use of the Box‑Cox method")),
        shiny::tags$p("When λ approaches 1, the simple SI introduces a scale‑changing translation. Adding a correction parameter θ yields the modified stimulation index (mSI):"),
        shiny::tags$p(htmltools::HTML("$$m\\text{SI}_{i}=g_{\\lambda,\\theta}(s_{i},u_{i}) = \\begin{cases} \\dfrac{s_{i}-u_{i}+1-\\lambda\\theta}{1-\\lambda}, & 0 < \\lambda \\le 1 ;\\; \\dfrac{s_{i}}{u_{i}}, & \\lambda = 0 \\end{cases}$$")),
        shiny::tags$p("Set θ = 0 to recover the simple SI; the default dynamic θ(λ) is recommended for λ > 0.7.")
      ),
      
      # -------------------- 4. Selecting a value of λ -----------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Selecting a value of λ")),
        shiny::tags$ul(
          shiny::tags$li("Multiplicative unstim–stim relationship → λ ≈ 0 (e.g., 0.1)."),
          shiny::tags$li("Pure background noise → λ ≈ 1."),
          shiny::tags$li("Uncertain relationship → λ = 0.5 (middle ground).")
        )
      ),
      
      # --------------- 5. Mathematical estimation of λ ----------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Mathematical estimation of optimal λ values")),
        shiny::tags$p("The app offers four complementary methods, each aiming for S (stimulated) to be uncorrelated with U (unstimulated):"),
        shiny::tags$ol(
          shiny::tags$li(shiny::strong("Method 1 – β‑coefficient:"), " regression slope S ~ U equals 1."),
          shiny::tags$li(shiny::strong("Method 2 – Spearman correlation:"), " residual correlation 0."),
          shiny::tags$li(shiny::strong("Method 3 – Posterior probability:"), " Bayesian approach maximises posterior density at β = 1."),
          shiny::tags$li(shiny::strong("Method 4 – Likelihood function:"), " maximises Box–Cox log‑likelihood with β fixed at 1.")
        )
      ),
      
      # ------------------ 6. Step‑by‑step guide -----------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Step‑by‑step guide to the Box‑Cox App")),
        shiny::tags$ol(
          shiny::tags$li("Upload your AIM .csv or use the example."),
          shiny::tags$li("(Optional) Estimate λ in the Lambda Estimation tab."),
          shiny::tags$li("Compute SI or mSI."),
          shiny::tags$li("Download transformed data.")
        )
      ),
      
      # --------------------- 7. File requirements ---------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("File requirements for AIM datasets")),
        shiny::tags$ul(
          shiny::tags$li("CSV format; one row per sample."),
          shiny::tags$li("Include grouping columns, a stimulant column, and AIM‑frequency columns."),
          shiny::tags$li("Column names: letters/numbers/underscores, no leading digits."),
          shiny::tags$li("Numeric values in decimal notation (e.g., 0.023)."),
          shiny::tags$li("See the example dataset for reference.")
        )
      ),
      
      # ------------------- 8. Example dataset -------------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Download example dataset")),
        shiny::tags$p("SARS‑CoV‑2 AIM data from 33 solid‑organ transplant recipients after three vaccine doses⁵."),
        shiny::tags$a(href = "20250328_PREVENT_RawData_V3.csv", target = "_blank", download = NA, "Download 20250328_PREVENT_RawData_V3.csv")
      ),
      
      # ----------------------- 9. References --------------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("References")),
        shiny::tags$ol(
          shiny::tags$li("Box GEP, Cox DR. An Analysis of Transformations. J Roy Stat Soc B. 1964."),
          shiny::tags$li("Congdon P. Bayesian Statistical Modeling. Wiley; 2006."),
          shiny::tags$li("Clyde M et al. Introduction to Bayesian Regression. 2025."),
          shiny::tags$li("White KJ. Estimation of the Liquidity Trap with a Generalized Functional Form. Econometrica. 1972."),
          shiny::tags$li("Halvorson T et al. Transplantation. 2024.")
        )
      )
    )
  )
}

aboutServer <- function(id) {
  moduleServer(id, function(input, output, session) {})
}
