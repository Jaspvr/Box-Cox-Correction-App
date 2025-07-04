aboutUI <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    title = "About",
    shiny::withMathJax(),
    
    shiny::tags$div(
      class = "container",
      
      # ------------------------------ 1. Overview ---------------------------
      # shiny::tags$details(
      #   shiny::tags$summary(shiny::tags$h4("Overview")),
      #   shiny::tags$p("Analysis of T‑cell activation‑induced marker (AIM) assay data requires normalization of AIM\u2009+ cell frequencies to ‘background’ AIM\u2009+ frequencies in an unstimulated control. Subtracting or dividing by the unstimulated control each have specific disadvantages and can amplify technical variability in the assay."),
      #   shiny::tags$p("The Box‑Cox correction, based on the transformation proposed by Box & Cox in 1964\u00b9, combines features of division and linear subtraction to better match the mathematical properties of AIM datasets, reducing technical variability and improving signal detection.")
      # ),
      shiny::tags$details(
        shiny::tags$summary(shiny::tags$h4("Overview")),
        shiny::tags$p("Analysis of T cell activation-induced marker (AIM) assay data requires normalization of AIM+ cell frequencies to ‘background’ AIM+ frequencies in an unstimulated control. Subtracting or dividing by the unstimulated control each have specific disadvantages in different situations and can amplify technical variability in the assay."),
        shiny::tags$p("The Box-Cox correction, based on the Box-Cox transformation initially proposed by Box and Cox in 1964\u00b9, is an innovative approach to AIM data analysis with features of both division and linear subtraction, allowing a more sophisticated correction for unstimulated AIM+ cell frequencies that more closely aligns with the mathematical properties of AIM datasets, reducing technical variability and confounding effects of variable ‘background’ AIM+ cell frequencies in the unstimulated condition, thereby enhancing signal detection.")
      ),
      
      # ------------------ 2. Application of the Box‑Cox method -------------
      # shiny::tags$details(
      #   shiny::tags$summary(shiny::tags$h4("Application of the Box‑Cox method to AIM assay data")),
      #   shiny::tags$ol(
      #     shiny::tags$li(htmltools::HTML("Apply the transformation: $$\\text{BoxCox}(x,\\;\\lambda \\in [0,1]) = \\begin{cases} \\dfrac{x^{\\lambda}-1}{\\lambda}, & \\lambda \\neq 0 \\; ; \\; \\ln x, & \\lambda = 0 \\end{cases}$$")),
      #     shiny::tags$li(htmltools::HTML("Net stimulated value: $$\\beta = \\text{BoxCox}(\\%\\text{AIM}^{+}_{\\text{stim}}) - \\text{BoxCox}(\\%\\text{AIM}^{+}_{\\text{unstim}})$$")),
      #     shiny::tags$li(htmltools::HTML("Return to original scale: $$\\text{SI} = \\text{invBoxCox}(\\beta,\\;\\lambda) = \\begin{cases} (\\lambda\\beta + 1)^{1/\\lambda}, & \\lambda \\neq 0 ;\\; e^{\\beta}, & \\lambda = 0 \\end{cases}$$"))
      #   ),
      #   shiny::tags$p("Choose λ (0–1) carefully, considering analysis goals and dataset properties.")
      # ),
      
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
      # shiny::tags$details(
      #   shiny::tags$summary(shiny::tags$h4("Advanced use of the Box‑Cox method")),
      #   shiny::tags$p("When λ approaches 1, the simple SI introduces a scale‑changing translation. Adding a correction parameter θ yields the modified stimulation index (mSI):"),
      #   shiny::tags$p(htmltools::HTML("$$m\\text{SI}_{i}=g_{\\lambda,\\theta}(s_{i},u_{i}) = \\begin{cases} \\dfrac{s_{i}-u_{i}+1-\\lambda\\theta}{1-\\lambda}, & 0 < \\lambda \\le 1 ;\\; \\dfrac{s_{i}}{u_{i}}, & \\lambda = 0 \\end{cases}$$")),
      #   shiny::tags$p("Set θ = 0 to recover the simple SI; the default dynamic θ(λ) is recommended for λ > 0.7.")
      # ),
      
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
