# Box-Cox-Correction-App
Data transformation web application built for a UBC Immunology lab working on developing cell therapies for immune tolerance in transplantation and autoimmunity.
https://www.bcchr.ca/levingslab

# Background
Analysis of T cell activation-induced marker (AIM) assay data requires normalization of AIM+ cell frequencies to background AIM+ frequencies in an unstimulated control. Subtracting or dividing by the unstimulated control each have specific disadvantages and can amplify technical variability in the assay. The Box-Cox correction is an innovative method with features of both division and linear subtraction, allowing a more sophisticated correction for unstimulated AIM+ cell frequencies that better aligns with the mathematical properties of AIM datasets and reduces technical variability.

# Usage
To take advantage of the Box-Cox correction, upload your full AIM dataset and the set of variables to be corrected. The Box-Cox Correction App will immediately return the corrected values which are then ready for data display or statistical analysis.
