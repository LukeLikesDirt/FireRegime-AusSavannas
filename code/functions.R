
### VARIATION INFLATION FACTOR FUNCTION ###

# To cite the VIF function, use:
# Mixed effects models and extensions in ecology with R. (2009).
# Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.

VIF = function(x) {
  x = as.data.frame(x)
  # VIF calculation
  form    = formula(paste("fooy ~ ", paste(strsplit(names(x), " "), collapse = " + ")))
  x       = data.frame(fooy = 1 + rnorm(nrow(x)) ,x)
  lm_mod  = lm(form, x)
  # End VIF calculation
  cat("\n\nVariance inflation factors\n\n")
  print(supVIF(lm_mod))
}
# VIF function dependency
supVIF = function(mod) {
  v = vcov(mod)
  assign = attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v = v[-1, -1]
    assign = assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms = labels(terms(mod))
  n.terms = length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor) = 0
    if (any(tmp_cor == 1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R = cov2cor(v)
  detR = det(R)
  result = matrix(0, n.terms, 3)
  rownames(result) = terms
  colnames(result) = c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs = which(assign == term)
    result[term, 1] = det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] = length(subs)
  }
  if (all(result[, 2] == 1)) {
    result = data.frame(GVIF=result[, 1])
  } else {
    result[, 3] = result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

### STANDARDISE AND UNSTANDARDISE CONTINUOUS COVARIATES ###

# Function to standardise predictors to have a mean of 0 and a standard 
# deviation of 1 
std = function(x) {(x - mean(x)) / sd(x)}
# This function back-transforms standardised predictors for plotting, where 
# 'x.std' represents the standardised predictor and 'x' is the observed 
# predictor value form the original data frame  
unstd = function(x.std, x) {x.std * sd(x) + mean(x)}

### GGPLOT THEME ###
MyTheme = function() {
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  )
}