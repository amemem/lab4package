#' A Reference Class for linear regression
#'
#' @description This linear regression reference class will compute
#' the coefficients, residuals, predicted values, standard errors
#' for both the coefficients and the residuals, t-values and
#' p-values, and the degrees of freedom of a model matrix based
#' on the provided formula and data set.
#'
#' @field formula A formula object
#' @field data A data set
#'
#' @importFrom ggplot2 ggplot aes geom_point stat_summary geom_hline geom_text xlab ylab labs scale_y_continuous theme element_blank element_rect rel unit element_text
#' @export

linreg = setRefClass(
  "linreg",
  fields = list(
    formula = "formula",
    data = "character",
    name = "vector",
    coeff = "vector",
    resi = "vector",
    predi = "vector",
    se_c = "vector",
    se_r = "numeric",
    df = "numeric",
    t = "vector",
    p = "vector",
    i = "numeric"),
  methods = list(
    initialize = function(formula, data) {
      "The constructor will return an object of class linreg."
      X = model.matrix(formula, data)
      y = all.vars(formula)

      X_T = t(X)
      X_I = solve(X_T %*% X)
      y = data[[all.vars(formula)[1]]]

      coef = X_I %*% X_T %*% y
      pred = X %*% coef
      resid = y - pred
      degf = nrow(X) - ncol(X)
      var_resid = as.numeric((t(resid) %*% resid) / degf)
      std_err_resid = sqrt(var_resid)
      var_coef = var_resid * X_I * diag(ncol(X_I))
      std_err_coef = sqrt(var_coef[var_coef > 0])
      t_val = coef / std_err_coef
      p_val = pt(-abs(t_val), degf)

      .self$formula = formula
      .self$data = deparse(substitute(data))
      .self$coeff = as.vector(coef)
      names(.self$coeff) = colnames(X)
      .self$predi = as.vector(pred)
      .self$resi = as.vector(resid)
      .self$se_c = as.vector(std_err_coef)
      .self$se_r = std_err_resid
      .self$df = degf
      .self$t = as.vector(t_val)
      .self$p = as.vector(p_val)
      .self$i = 1
    },
    coef = function() {
      "Returns the coefficients."
      return(.self$coeff)
    },
    pred = function() {
      "Returns the predicted values."
      return(.self$predi)
    },
    resid = function() {
      "Returns the residuals."
      return(.self$resi)
    },
    print = function() {
      "Prints how the function was called and the coefficients."
      cat(paste0("linreg(formula = ", format(.self$formula),
                 ", data = ", .self$data, ")\n"))
      base::print(.self$coeff)
    },
    summary = function() {
      "Prints a summary of the model."
      #for (i in 1:length(.self$coeff)) {
      if (.self$i <= length(.self$coeff)) {
        cat(paste0(names(.self$coeff)[.self$i], " ",
                  round(.self$coeff[.self$i], 2), " ",
                  round(.self$se_c[.self$i], 2), " ",
                  round(.self$t[.self$i], 2), " ",
                  round(.self$p[.self$i], 2), "\n"))
        .self$i = .self$i + 1
      }
      else {
        cat(paste0("Residual standard error: ", round(.self$se_r, 2),
                 " on ", .self$df, " degrees of freedom\n"))
        .self$i = 1
      }
    },
    plot = function() {
      "Plots two graphs; residuals vs fitted values, and scale-location."
      df_resi = as.data.frame(.self$resi)
      df_resi$Fitted = .self$predi
      colnames(df_resi) = c("Residuals", "Fitted")
      df_resi$rows = as.numeric(rownames(df_resi))
      df_resi$absolute = abs(.self$resi)
      df_resi = df_resi[order(-df_resi$absolute),]

      residuals_vs_fitted = ggplot2::ggplot(df_resi, ggplot2::aes(x = Fitted, y = Residuals)) +
        ggplot2::geom_point(shape = 1, alpha = 0.3, colour = "black", fill = "white", size = 2) +
        ggplot2::stat_summary(fun = "median", colour = "red", geom = "line") +
        ggplot2::geom_hline(yintercept = 0, alpha = 0.3, linetype = 3) +
        ggplot2::geom_text(data = as.data.frame(df_resi[1:3,]), ggplot2::aes(x = Fitted, y = Residuals, label = rows, hjust = "inward")) +
        ggplot2::xlab(paste0("Fitted values\n", "linreg(", format(.self$formula), ")")) +
        ggplot2::labs(title = "Residuals vs Fitted") +
        ggplot2::scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5), labels = c("-1.5", "", "-0.5", "", "0.5", "", "1.5")) +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(NA, "black", ggplot2::rel(1)),
          axis.ticks.length = ggplot2::unit(.25, "cm"),
          axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
          axis.text.y = ggplot2::element_text(angle = 90, hjust = .5),
          plot.title = ggplot2::element_text(hjust = .5)
        )
      base::print(residuals_vs_fitted)

      df_resi = df_resi[order(df_resi$rows),]
      df_resi$Residuals = sqrt(abs(.self$resi / .self$se_r))
      df_resi = df_resi[order(-df_resi$absolute),]

      scale_location = ggplot2::ggplot(df_resi, ggplot2::aes(x = Fitted, y = Residuals)) +
        ggplot2::geom_point(shape = 1, alpha = 0.3, colour = "black", fill = "white", size = 2) +
        ggplot2::stat_summary(fun = "mean", colour = "red", geom = "line") +
        ggplot2::geom_text(data = as.data.frame(df_resi[1:3,]), ggplot2::aes(x = Fitted, y = Residuals, label = rows, hjust = "inward")) +
        ggplot2::xlab(paste0("Fitted values\n", "linreg(", format(.self$formula), ")")) +
        ggplot2::ylab(quote(sqrt(abs("Standardized residuals")))) +
        ggplot2::labs(title = "Scale-Location") +
        ggplot2::scale_y_continuous(limits = c(0, 1.8), breaks = seq(0, 1.5, 0.5), labels = c("0.0", "0.5", "1.0", "1.5")) +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(NA, "black", ggplot2::rel(1)),
          axis.ticks.length = ggplot2::unit(.25, "cm"),
          axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
          axis.text.y = ggplot2::element_text(angle = 90, hjust = .5),
          plot.title = ggplot2::element_text(hjust = .5)
        )
      base::print(scale_location)
    }))
