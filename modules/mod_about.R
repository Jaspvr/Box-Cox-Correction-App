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
        shiny::tags$p("Analysis of T cell activation-induced marker (AIM) assay data requires normalization of AIM+ cell frequencies to ‘background’ AIM+ frequencies in an unstimulated control. Subtracting or dividing by the unstimulated control each have specific disadvantages in different situations and can amplify technical variability in the assay."),
        shiny::tags$p("The Box-Cox correction, based on the Box-Cox transformation initially proposed by Box and Cox in 1964\u00b9, is an innovative approach to AIM data analysis with features of both division and linear subtraction, allowing a more sophisticated correction for unstimulated AIM+ cell frequencies that more closely aligns with the mathematical properties of AIM datasets, reducing technical variability and confounding effects of variable ‘background’ AIM+ cell frequencies in the unstimulated condition, thereby enhancing signal detection.")
      ),
      
      # ------------------ 2. Application of the Box‑Cox method -------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Application of the Box-Cox method to AIM assay data")),
        shiny::tags$p("The simplest method of Box-Cox transformation for AIM assay data takes the following three steps, incorporating the stimulated and unstimulated AIM values and an arbitrary parameter \u03BB:"),
        
        shiny::tags$ol(
          # --- step 1 ----------------------------------------------------
          shiny::tags$li(
            "Apply the transformation:",
            htmltools::HTML("$$\\operatorname{BoxCox}(x,\\; \\lambda \\in [0,1]) = \\begin{cases} \\dfrac{x^{\\lambda}-1}{\\lambda}, & \\lambda \\neq 0 \\;\\;\\; \\\\ \\ln x, & \\lambda = 0 \\end{cases}$$")
          ),
          # --- step 2 ----------------------------------------------------
          shiny::tags$li(
            "Obtain the net stimulated value (\u03B2) using subtraction:",
            htmltools::HTML("$$\\beta = \\operatorname{BoxCox}(\\%\\text{AIM}^{+}_{\\text{stim}}) - \\operatorname{BoxCox}(\\%\\text{AIM}^{+}_{\\text{unstim}})$$")
          ),
          # --- step 3 ----------------------------------------------------
          shiny::tags$li(
            "Return the transformed data to the original scale to obtain the stimulation index (SI):",
            htmltools::HTML("$$\\text{SI} = \\operatorname{invBoxCox}(\\beta,\\; \\lambda \\in [0,1]) = \\begin{cases} (\\lambda\\beta + 1)^{1/\\lambda}, & \\lambda \\neq 0 \\;\\;\\; \\\\ e^{\\beta}, & \\lambda = 0 \\end{cases}$$")
          )
        ),
        
        shiny::tags$p("Note that the parameter \u03BB \u2208 [0,1] is chosen by the user and should be selected carefully in consideration of the goals of the analysis and the unique properties of the dataset.")
      ),
      
      # ------------------------ 3. Advanced use -----------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Advanced use of the Box‑Cox method")),
        shiny::tags$p("A limitation of the simple Box‑Cox method above is the presence of a linear translation when λ approaches 1, which changes the scale of the data. To overcome this, the method can be modified to introduce a correction parameter θ, generating a modified stimulation index (mSI):"),
        shiny::tags$p(
          htmltools::HTML("$$m\\text{SI}_{i}=g_{\\lambda,\\theta}(s_{i},u_{i}) = \\begin{cases} \\left(s_{i}^{\\lambda}-u_{i}^{\\lambda}+1-\\lambda\\theta\\right)^{\\tfrac{1}{\\lambda}}, & 0 < \\lambda \\le 1 \\\\ \\dfrac{s_{i}}{u_{i}}, & \\lambda = 0 \\end{cases}$$")
        ),
        shiny::tags$p(
          htmltools::HTML("where&nbsp;\\(s_i\\) and&nbsp;\\(u_i\\) are the stimulated and unstimulated AIM<sup>+</sup> frequencies, respectively.&nbsp; The parameter&nbsp;\\(\\theta\\)&nbsp;can be set as a constant&nbsp;(0&nbsp;≤&nbsp;\\(\\theta\\)&nbsp;≤&nbsp;1), or preferably, as a function&nbsp;of&nbsp;λ:")
        ),
        shiny::tags$p(
          htmltools::HTML("$$\\theta(\\lambda) = \\begin{cases} 1 - (1-\\lambda)^{\\tfrac{\\lambda}{1-\\lambda}}, & \\lambda < 0.5 \\\\ \\lambda^{\\tfrac{\\lambda}{1-\\lambda}}, & \\lambda \\ge 0.5 \\end{cases}$$")
        ),
        shiny::tags$p("This mSI is strongly recommended when using larger values of λ (e.g., > 0.7) to avoid the issues of translation and scaling. Note that when θ = 0, the mSI is equivalent to the simple SI described above.")
      ),
      
      # -------------------- 4. Selecting a value of λ -----------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Selecting a value of λ")),
        
        shiny::tags$p(
          "The value of λ should be determined based on the expected relationship ",
          "between stimulated and unstimulated AIM⁺ frequencies in a given dataset. ",
          "In general, if a multiplicative relationship between the unstimulated and ",
          "stimulated AIM⁺ frequencies is observed or expected, a choice of λ close ",
          "to 0 is recommended. However, use of a (small) non-zero value for λ, such ",
          "as 0.1, will reduce the impact of very small observations for unstimulated ",
          "values, compared with simple division."
        ),
        
        shiny::tags$p(
          "If we observe or expect the unstimulated values to be random background ",
          "‘noise’, with no relationship to the stimulated AIM⁺ frequencies, then a ",
          "choice of λ close to 1 is recommended."
        ),
        
        shiny::tags$p(
          "In reality, different factors can result in different relationships ",
          "between unstimulated and stimulated values. Batch effects, unexpected ",
          "activation of cells, differences in reactivity or cell survival after ",
          "freeze/thaw, etc., can have unexpected effects on the experiment. When ",
          "the relationship between unstimulated and stimulated is not clear, a ",
          "choice of λ = 0.5 is suggested as a mathematical intermediate between ",
          "the extremes of division and subtraction."
        )
      ),
      
      # --------------- 5. Mathematical estimation of λ ----------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Mathematical estimation of optimal λ values")),
        
        shiny::tags$p(
          "Although a λ value of 0.5 is adequate for many datasets, in some cases it may ",
          "be desirable to further explore the properties of an AIM dataset to empirically ",
          "estimate an optimal λ-value. We developed four complementary methods for this, ",
          "all based on the assumption that the transformed value of the stimulated AIM⁺ ",
          "frequency (S) should be non-correlated with the unstimulated data x₀. In other ",
          "words, the objective is that the correction made to the stimulated data x₁ to ",
          "obtain the effect of stimulation is effective at removing artifacts or background ",
          "effects."
        ),
        
        shiny::tags$ol(
          
          ## ------------------------------------------------------------------ ##
          shiny::tags$li(
            shiny::strong("Method 1 – β coefficient."), " ",
            "Under a linear regression framework, we can consider:",
            htmltools::HTML("$$f_S = \\alpha + \\beta f_U + \\varepsilon$$"),
            "where β = 1 for the optimal λ. This method applies linear regression, using ",
            "the Box-Cox transformation with a set of values for λ to search for a λ value ",
            "such that the estimated slope β = 1."
          ),
          
          ## ------------------------------------------------------------------ ##
          tags$li(
            # already inside withMathJax() at the top of the tabPanel
            htmltools::HTML(paste0(
              "<strong>Method&nbsp;2 – Spearman correlation.</strong><br/>",
              "Under the assumption above, there should be zero correlation between&nbsp;S ",
              "and the unstimulated AIM⁺ frequencies:<br/>",
              "\\[\\operatorname{cor}\\bigl(f_S - f_U,\\; f_U\\bigr)=0\\]",      # display
              "<p>",
              "where \\( S = (s_1, s_2, \\dots, s_n) \\) and ",
              "\\( U = (u_1, u_2, \\dots, u_n) \\) are the \\( n \\) IID ",
              "(independent, identically distributed) measurements for the stimulated ",
              "and unstimulated conditions, respectively. Therefore, this method finds ",
              "the value of \\( \\lambda \\) (if any) that satisfies this zero-correlation ",
              "relationship for the SI. A positive residual correlation indicates an ",
              "under-correction from the unstimulated data, while a negative correlation ",
              "indicates an over-correction.",
              "</p>"
            ))
          ),
          
          
          
          ## ------------------------------------------------------------------ ##
          ## ---------------- Method 3 – Posterior probability distribution -----------
          tags$li(
            htmltools::HTML(paste0(
              "<strong>Method&nbsp;3 – Posterior probability distribution method.</strong><br/>",
              "Under a Bayesian framework for the linear model described above, with a ",
              "non-informative prior and normality assumptions, the posterior distribution ",
              "of \\(\\beta \\mid S\\) follows a \\(t\\)-distribution with estimated mean ",
              "\\(\\hat{\\beta}\\) (the regression slope) and scale parameter ",
              "\\(SE^{2}\\) (its squared standard error). Hence, this method searches for the ",
              "value of \\(\\lambda\\) that maximises the posterior density at ",
              "\\(\\beta = 1 \\mid S\\)."
            ))
          ),
          
          ## ---------------- Method 4 – Likelihood function --------------------------
          tags$li(
            htmltools::HTML(paste0(
              "<strong>Method&nbsp;4 – Likelihood function method.</strong><br/>",
              "Assuming normally distributed AIM values, the likelihood of ",
              "\\(\\lambda\\) is the probability of the observed dependent variable ",
              "\\(S\\) under the transformation with parameter \\(\\lambda\\) and the linear ",
              "model above with fixed slope \\(\\beta = 1\\). The likelihood equals the ",
              "probability of the transformed data \\(f_S\\) multiplied by the Jacobian. ",
              "It is more convenient to use the maximised log-likelihood, which (Box&nbsp;& Cox 1964) is<br/>",
              "\\[ \\mathscr{L}(\\lambda) = -\\tfrac{1}{2}\\,n\\,\\log\\sigma^{2} + (\\lambda - 1) ",
              "\\sum_{i=1}^{n} \\log s_i \\]<br/>",
              "where \\(\\sigma^{2}\\) is the residual sum of squares divided by \\(n\\) and the ",
              "second term is the log-Jacobian. The optimal \\(\\lambda\\) maximises ",
              "\\(\\mathscr{L}(\\lambda)\\); a 95&nbsp;% confidence interval is given by<br/>",
              "\\[ \\mathscr{L}_{\\max} - \\mathscr{L}(\\lambda) < \\tfrac{1}{2}\\,\\chi^{2}_{1;0.05} \\]"
            ))
          )
          
          
        )  # end <ol>
      ), 
      # ------------------ 6. Step‑by‑step guide -----------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Step-by-step guide to the Box-Cox App")),
        shiny::tags$div(
          
          ## ---------------------------------------------------------------------
          shiny::tags$p(
            "You may either upload your own AIM assay dataset in .csv format, ",
            "or alternatively, start with an example dataset to familiarize ",
            "yourself with the functions of the Box-Cox App. The example dataset ",
            "can be downloaded from the following link:",
            shiny::tags$a(
              href     = "20250328_PREVENT_RawData_V3.csv",
              target   = "_blank",
              download = NA,
              "20250328_PREVENT_RawData_V3.csv"
            )
          ),
          shiny::tags$p(
            "These data are previously published SARS-CoV-2 AIM assay data from 33 ",
            "solid organ transplant recipients measured after 3 doses of a COVID-19 ",
            "mRNA vaccine", shiny::tags$sup("5"), "."
          ),
          
          ## ---------------- File requirements ----------------------------------
          shiny::tags$h5("File requirements for AIM assay datasets:"),
          shiny::tags$ul(
            shiny::tags$li("Data must be in .csv format."),
            shiny::tags$li("Each row should represent one unique sample."),
            shiny::tags$li(
              "Each column should represent a variable that identifies or groups ",
              "samples, or contains the raw AIM+ frequencies for each sample. ",
              "Columns can contain categorical or numeric data."
            ),
            shiny::tags$ul(
              shiny::tags$li(
                shiny::HTML(
                  "<strong>Grouping variables:</strong> These columns contain variables ",
                  "that, together, can uniquely identify each sample (e.g., ‘Donor’, ",
                  "‘Timepoint’ or ‘Treatment Group’ in clinical studies)."
                )
              ),
              shiny::tags$li(
                shiny::HTML(
                  "<strong>Stimulant variable:</strong> One column should contain the ",
                  "names of the stimulants or antigens used to stimulate each sample ",
                  "(e.g., ‘SARS-CoV-2’, ‘CMV’, ‘Unstimulated’ or ‘PHA’)."
                )
              ),
              shiny::tags$li(
                shiny::HTML(
                  "<strong>AIM variables:</strong> These columns should contain numerical ",
                  "data representing frequencies of AIM+ cells measured in each sample. ",
                  "There should be one column for each AIM or AIM pair of interest in ",
                  "CD4+ or CD8+ T cells or other subsets (e.g., ‘CD134pCD25p_CD4p’ could ",
                  "represent frequencies of CD4+ T cells that are CD134+/CD25+)."
                )
              ),
              shiny::tags$li(
                shiny::HTML(
                  "<strong>Other columns:</strong> Other columns may be present in the ",
                  "data that are not useful for calculating SI values or identifying ",
                  "unique samples (e.g., columns containing frequencies of total CD4+ ",
                  "or CD8+ T cells). These columns will be ignored by default."
                )
              )
            ),
            shiny::tags$li(
              "Numerical data should use decimal format (e.g., 0.023 rather than 0,023)."
            ),
            shiny::tags$li(
              "Column names should contain letters, numbers and underscores (‘_’) only ",
              "(no spaces), but should not begin with a number."
            ),
            shiny::tags$li(
              "Please see the example dataset above as a reference for the correct ",
              "file format."
            )
          ),
          
          ## ---------------- Lambda estimation ----------------------------------
          shiny::tags$h5("Estimating an optimal value of λ (optional)"),
          shiny::tags$ol(
            shiny::tags$li("Navigate to the Lambda (λ) Estimation Tools tab."),
            shiny::tags$li(
              "Select ‘Browse’ in the ‘Upload CSV’ box and upload your AIM dataset ",
              "in .csv format. Ensure the .csv file complies with the file ",
              "requirements detailed above. A message reading ‘upload complete’ will appear."
            ),
            shiny::tags$li(
              "Select the column containing the stimulant/antigen values from the ",
              "dropdown menu in the ‘Stimulant Column’ box."
            ),
            shiny::tags$li(
              "Enter the names of the stimulant/antigen condition for which λ will be ",
              "estimated into the ‘Stimulant’ box. Ensure the name entered is identical ",
              "to the spelling in the data table."
            ),
            shiny::tags$li(
              "Enter the name of the unstimulated parameter. Ensure the name entered ",
              "is identical to the spelling in the data table."
            ),
            shiny::tags$li(
              "Enter the desired value of Epsilon (ϵ). This should be a small non-zero ",
              "value that will be added to each data point to avoid zero values, which ",
              "is necessary for the λ-estimating functions to operate. The default is ",
              "set to 0.001."
            ),
            shiny::tags$li(
              "Select the AIM variable for which λ will be estimated from the dropdown ",
              "menu in the ‘AIM Variable’ box."
            ),
            shiny::tags$li(
              shiny::tags$strong("Optional:"), " Filter for a desired subset of the data using the ‘Filter columns (optional)’ tab. ",
              "This can be helpful when one or more of the λ-estimation methods fails to estimate an optimal value of lambda on the entire dataset. ",
              "Analyses can then be run individually on specific subsets of the data.",
              shiny::tags$ul(
                shiny::tags$li("Select one or more columns for which to filter for desired values."),
                shiny::tags$li("A checklist will appear containing values present in that column."),
                shiny::tags$li("De-select the values you wish to exclude from that column, leaving the desired values selected.")
              )
            ),
            shiny::tags$li("Select ‘Run’.")
          ),
          
          ## ---------------- Interpreting the data ------------------------------
          shiny::tags$h5("Interpreting the data:"),
          shiny::tags$p(
            "Four graphs will be displayed, corresponding to the outputs of the four ",
            "methods of estimating λ (described above)."
          ),
          shiny::tags$p(
            "The optimal value of λ determined by each method will be displayed on the ",
            "plot. These values can then be used in calculation of SI and mSI (see below)."
          ),
          
          ## ---------------- Simple SI ------------------------------------------
          shiny::tags$h5("Simple Box-Cox Stimulation Index (SI)"),
          shiny::tags$ol(
            shiny::tags$li(
              "Enter the desired value of λ in the ‘Input lambda value’ box."
            ),
            shiny::tags$li(
              "Enter the name of the unstimulated parameter in the ‘Unstimulated ",
              "Parameter’ box. Ensure the name entered is identical to the spelling ",
              "in the data table."
            ),
            shiny::tags$li(
              "Enter the names of the stimulants/antigens used in the dataset into the ",
              "‘Stimulants’ box. Stimulant names should be separated by commas. Ensure ",
              "the names entered are identical to the spelling in the data table."
            ),
            shiny::tags$li(
              "Upload the AIM assay dataset of interest in .csv format by selecting ",
              "‘browse’ in the ‘Input Patient Data (CSV)’ box. Ensure the .csv file ",
              "complies with the file requirements detailed above. A message will ",
              "appear when the upload is complete, and a preview of the data table ",
              "will be visible in the right panel."
            ),
            shiny::tags$li(
              "Click on the ‘AIM Variables’ box and use the dropdown menu to select ",
              "the column names that correspond to the AIM+ frequencies of each sample. ",
              "These are the columns for which an SI will be calculated. Do not select ",
              "columns that do not represent AIMs (such as total frequencies of CD4+ ",
              "events among CD3+ T cells), as these will give uninformative results."
            ),
            shiny::tags$li(
              "Click on the ‘Grouping Columns’ box and use the dropdown menu to select ",
              "the column names that correspond to the grouping variables. These are ",
              "the variables that together, will uniquely identify each sample."
            ),
            shiny::tags$li(
              "Click on the ‘Stimulant Column’ box and use the dropdown menu to select ",
              "the column name corresponding to the antigen or stimulant variable. ",
              "This column should contain the names of the antigen/stimulant used to ",
              "stimulate each sample and is critical for matching stimulated and ",
              "unstimulated AIM+ frequencies."
            ),
            shiny::tags$li(
              "Select ‘Download Transformed Data’. A .csv file containing the ",
              "Box-Cox-corrected SI values for the AIM dataset will be downloaded in ",
              "the browser. These SI data are now ready for downstream visualization, ",
              "analysis and interpretation. Note that Box-Cox-corrected SI values for ",
              "all samples in the unstimulated control condition will equal 1."
            )
          ),
          
          ## ---------------- Advanced mSI ---------------------------------------
          shiny::tags$h5("Advanced Scale-Independent Modified Stimulation Index (mSI)"),
          shiny::tags$ol(
            shiny::tags$li(
              "Enter the desired value of λ in the ‘Lambda (L)’ box."
            ),
            shiny::tags$li(
              "Check the box labelled ‘Apply F1(L) correction’. This applies the ",
              "scale-independent modified stimulation index (mSI) formula. Leaving this ",
              "box unchecked will instead perform the calculation for the simple ",
              "Box-Cox stimulation index (SI) described above."
            ),
            shiny::tags$li(
              "Optional: Specify a value of theta in the box labelled ‘Theta (H)’. ",
              "Leaving this box blank will apply the default formula for θ, which ",
              "defines θ as a function of λ. We recommend using this default."
            ),
            shiny::tags$li(
              "Specify a value of Epsilon in the ‘Epsilon’ box. This adds a small value ",
              "(e.g., 0.001) to all AIM values prior to calculating the stimulation index, ",
              "to avoid the presence of zeros in the dataset. This is recommended if ",
              "your dataset contains any values of zero for AIM+ frequency values."
            ),
            shiny::tags$li(
              "In the box labelled ‘Replacement value for undefined SI’, enter the ",
              "desired value. This value will serve as a default to replace SI values ",
              "that cannot be calculated, in cases where the AIM+ frequency of the ",
              "unstimulated condition is much greater than the corresponding AIM+ ",
              "frequency in the stimulated condition."
            ),
            shiny::tags$li(
              "Enter the name of the unstimulated parameter in the ‘Unstimulated ",
              "Parameter’ box. Ensure the name entered is identical to the spelling ",
              "in the data table."
            ),
            shiny::tags$li(
              "Enter the names of the stimulants/antigens used in the dataset into the ",
              "‘Stimulants’ box. Stimulant names should be separated by commas. Ensure ",
              "the names entered are identical to the spelling in the data table."
            ),
            shiny::tags$li(
              "Upload the AIM assay dataset of interest in .csv format by selecting ",
              "‘browse’ in the ‘Input Patient Data (CSV)’ box. Ensure the .csv file ",
              "complies with the file requirements detailed above. A message will ",
              "appear when the upload is complete, and a preview of the data table ",
              "will be visible in the right panel."
            ),
            shiny::tags$li(
              "Click on the ‘AIM Variables’ box and use the dropdown menu to select ",
              "the column names that correspond to the AIM+ frequencies of each sample. ",
              "These are the columns for which a mSI will be calculated. Do not select ",
              "columns that do not represent AIMs (such as total frequencies of CD4+ ",
              "events among CD3+ T cells), as these will give uninformative results."
            ),
            shiny::tags$li(
              "Click on the ‘Grouping Columns’ box and use the dropdown menu to select ",
              "the column names that correspond to the grouping variables. These are ",
              "the variables that together, will uniquely identify each sample."
            ),
            shiny::tags$li(
              "Click on the ‘Stimulant Column’ box and use the dropdown menu to select ",
              "the column name corresponding to the antigen or stimulant variable. ",
              "This column should contain the names of the antigen/stimulant used to ",
              "stimulate each sample and is critical for matching stimulated and ",
              "unstimulated AIM+ frequencies."
            ),
            shiny::tags$li(
              "Select ‘Download Transformed Data’. A .csv file containing the ",
              "Box-Cox-corrected mSI values for the AIM dataset will be downloaded in ",
              "the browser. These SI data are now ready for downstream visualization, ",
              "analysis and interpretation. Note that Box-Cox-corrected mSI values for ",
              "all samples in the unstimulated control condition will equal 0.25 when ",
              "default values are used."
            )
          )
        )
      ),
      
      # ----------------------- 9. References --------------------------------
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("References")),
        shiny::tags$ol(
          shiny::tags$li("\u00A0Box GEP, Cox DR. An Analysis of Transformations. Journal of the Royal Statistical Society: Series B (Methodological). 1964;26(2):211-243. doi:10.1111/j.2517-6161.1964.tb00553.x"),
          shiny::tags$li("\u00A0Congdon P. In: Bayesian Statistical Modeling. Second. Wiley; 2006:113."),
          shiny::tags$li("\u00A0Clyde M, Çetinkaya-Rundel M, Rundel C, Banks D, Chai C, Huang L. Chapter 6 Introduction to Bayesian Regression | An Introduction to Bayesian Thinking. Accessed June 7, 2025. https://statswithr.github.io/book/introduction-to-bayesian-regression.html"),
          shiny::tags$li("\u00A0White KJ. Estimation of the Liquidity Trap with a Generalized Functional Form. Econometrica. 1972;40(1):193-199. doi:10.2307/1909730"),
          shiny::tags$li("\u00A0Halvorson T, Ivison S, Huang Q, et al. SARS-CoV-2 Variants Omicron BA.4/5 and XBB.1.5 Significantly Escape T Cell Recognition in Solid-organ Transplant Recipients Vaccinated Against the Ancestral Strain. Transplantation. 2024;108(4):e49-e62. doi:10.1097/TP.0000000000004873")
        )
      )
    )
  )
}

aboutServer <- function(id) {
  moduleServer(id, function(input, output, session) {})
}
