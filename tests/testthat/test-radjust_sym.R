context("test radjust_sym()")
test_that(
  "The examples work and the outputs are equal to the saved copies",
  {
    data(mice)

    pv1 <- ifelse(mice$dir_is_left1, mice$twosided_pv1/2, 1-mice$twosided_pv1/2)
    pv2 <- ifelse(mice$dir_is_left2, mice$twosided_pv2/2, 1-mice$twosided_pv2/2)

    mice_rv_adaptive <- radjust_sym(pv1, pv2, input_type = "all",
                                    directional_rep_claim = T, variant = "adaptive")

    mice_rv_non_adpt_sel <- radjust_sym(pv1, pv2, input_type = "all",
                                    directional_rep_claim = T, variant = "non-adaptive-with-alpha-selection")

    data(examples_outputs_for_testing)
    expect_equivalent(mice_rv_adaptive, examples_outputs_for_testing$mice_rv_adaptive)
    expect_equivalent(mice_rv_non_adpt_sel, examples_outputs_for_testing$mice_rv_adjusted)
})

test_that(
  "a basic example works & the output is in the appropriate format",
  {
    pv1 <- c((1:20)/150)
    pv2 <- c((1:20)/150)
    output <- radjust_sym(pv1, pv2)
    expect_s3_class(output,"radjust")
    expect_type(output,"list")
    expect_length(output,9)
  })

test_that(
  "input validations works, so the function output the appropriate errors\\warnings",
  {
    pv1 <- c((1:20)/150)
    pv2 <- c((1:20)/150)

    expect_error(radjust_sym(pv1,pv2 = "a"), regexp = 'is.vector.*is not TRUE')
    expect_error(radjust_sym(pv1), regexp = 'argument "pv2" is missing, with no default')
    expect_error(radjust_sym(pv1,pv2 = 1:10), regexp = 'pv2 <= 1 are not all TRUE')
    expect_error(radjust_sym(pv1,pv2, w1 = 2), regexp = 'w1 < 1 is not TRUE')
    expect_error(radjust_sym(pv1,pv2, general_dependency = 1234), regexp = 'general_dependency')
    expect_error(radjust_sym(pv1,pv2, directional_rep_claim = 1234), regexp = 'directional_rep_claim')

    expect_error(radjust_sym(pv1,pv2, alpha = NULL), regexp = "variant != 'non-adaptive' but alpha value is not specified")
    expect_warning(radjust_sym(pv1,pv2, variant = "non-adaptive", alpha = 0.05), regexp = "variant == 'non-adaptive' so alpha value is ignored")
  })



test_that(
  "All parameters configurations work, without errors\\warnings, except when needed",
  {
    # TODO :complete this test
    # data(mice)
    #
    # pv1 <- ifelse(mice$dir_is_left1, mice$twosided_pv1/2, 1-mice$twosided_pv1/2)
    # pv2 <- ifelse(mice$dir_is_left2, mice$twosided_pv2/2, 1-mice$twosided_pv2/2)
    #
    # input <- expand.grid(
    #   pv1 = list(pv1),
    #   pv2 = list(pv2),
    #   w1 = 0.5,
    #   input_type = c("all","selected"),
    #   general_dependency = c(T,F),
    #   directional_rep_claim = c(T,F),
    #   variant = c("alpha-adjusted", "alpha-adaptive", "none"),
    #   alpha = .05
    #   )
    #
    # result <- expand.grid(
    #   pv1 = T,
    #   pv2 = T,
    #   w1 = T,
    #   input_type = c(T,T),
    #   general_dependency = c(T,T),
    #   directional_rep_claim = c(T,T),
    #   variant = c(T,T,T),
    #   alpha = enquote('T'))
    #
    # res <- do.call(what = radjust_sym(), args = input)
  })


#
# data(crohn)
# crohn_rv1 <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547, l00 = 0.8)
# crohn_rv2 <- radjust_pf(pv1 = crohn$pv1, pv2 = crohn$pv1, m = 635547, l00 = 0.8, variation="use_threshold",tt = 1e-5)
#
# examples_outputs_for_testing <- list(mice_rv_adaptive = mice_rv_adaptive,
#                                      mice_rv_adjusted = mice_rv_adjusted,
#                                      crohn_rv1 = crohn_rv1,
#                                      crohn_rv2 = crohn_rv2)


