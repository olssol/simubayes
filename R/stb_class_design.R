## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##   DESCRIPTION:
##       DEFINE SIMULATION TOOLBOX CLASSES OF DESIGNS
##
##
##   DESIGNS:
##      1. STB_DESIGN_BAYES_1ARM
##      2. STB_DESIGN_BAYES_2ARM
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##       HELPER FUNCTIONS
## -----------------------------------------------------------------------------

#' Create a Study Design
#'
#'
#'
#' @export
#'
stb_create_design <- function(type = c("bayes_1arm",
                                       "bayes_2arm",
                                       "nb_borrow")) {

    type <- match.arg(type)
    rst  <- switch(type,
                   bayes_1arm     = new("STB_DESIGN_BAYES_1ARM"),
                   bayes_2arm     = new("STB_DESIGN_BAYES_2ARM"),
                   nb_borrow      = new("STB_DESIGN_NB_BORROW"),
                   new("STB_DESIGN"))

    rst
}


## -----------------------------------------------------------------------------
##             Bayesian 1- and 2-arm design
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN_BAYES",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_BAYES",
          function(x, ...) {
              callNextMethod()
              bayes_describe(x, ...)
          })

setMethod("stb_generate_data",
          "STB_DESIGN_BAYES",
          function(x, ...) {
              bayes_gen_data(x@design_para, ...)
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_BAYES",
          function(x, data_ana, ...) {
              callNextMethod()
              rst_post <- bayes_ana_data(
                  data_ana[[1]],
                  x@design_para$prior_by_arm,
                  n_post = x@design_para$n_post,
                  x      = x@design_para$x_post,
                  ...)

              rst_decision <- bayes_ana_decision(
                  rst_post,
                  decision_ref    = x@design_para$decision_ref,
                  decision_h0     = x@design_para$decision_h0,
                  decision_gl     = x@design_para$decision_gl,
                  decision_thresh = x@design_para$decision_thresh)

              list(rst_post    = rst_post,
                   rst_diff    = rst_decision$rst_diff,
                   rst_success = rst_decision$rst_success)
          })

setMethod("stb_simu_gen_raw",
          "STB_DESIGN_BAYES",
          function(x, lst, ...) {
              rst <- list()
              for (i in seq_len(length(lst))) {
                  cur_lst     <- data.frame(lst[[i]]$rst_success)
                  cur_lst$rep <- i
                  rst[[i]]    <- cur_lst
              }

              list(rbindlist(rst))
          })

setMethod("stb_simu_gen_summary",
          "STB_DESIGN_BAYES",
          function(x, lst, ...) {
              rst <- lst[[1]] %>%
                  group_by(arm) %>%
                  summarize(success = mean(success))
              list(rst)
          })

#'
#' @export
#'
setClass("STB_DESIGN_BAYES_1ARM",
         contains = "STB_DESIGN_BAYES")

setMethod("stb_set_default_para",
          "STB_DESIGN_BAYES_1ARM",
          function(x) {
              internal_bayes1arm_dpara()
          })

#'
#' @export
#'
setClass("STB_DESIGN_BAYES_2ARM",
         contains = "STB_DESIGN_BAYES")

setMethod("stb_set_default_para",
          "STB_DESIGN_BAYES_2ARM",
          function(x) {
              internal_bayes2arm_dpara()
          })


## -----------------------------------------------------------------------------
##                      Negative Binomial Borrowing Data
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN_NB_BORROW",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_NB_BORROW",
          function(x, ...) {
              callNextMethod()
              rcurrent_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_NB_BORROW",
          function(x) {
              internal_rcurrent_dpara()
          })

setMethod("stb_generate_data",
          "STB_DESIGN_NB_BORROW",
          function(x, ...) {
              rcurrent_gen_data(x@design_para, ...)
          })

#'
#' @export
#'
setMethod("stb_create_analysis_set",
          "STB_DESIGN_NB_BORROW",
          function(x, data,
                   type          = c("min_fu", "fix_fu"),
                   fu_days       = 12 * 7,
                   min_fu_days   = 12 * 7,
                   pt_proportion = 1,
                   ...) {


              if (is.null(data))
                  return(NULL)

              type     <- match.arg(type)
              dat_full <- switch(
                  type,
                  min_fu = rcurrent_day_eos_1(data,
                                              min_fu_days   = min_fu_days,
                                              pt_proportion = pt_proportion),

                  fix_fu = rcurrent_day_eos_2(data,
                                              fu_days       = fu_days,
                                              pt_proportion = pt_proportion)
              )

              dat_full <- rcurrent_censor(dat_full)
              dat_nb   <- rcurrent_get_nb(dat_full)

              list(data    = dat_full,
                   data_nb = dat_nb)
          })
