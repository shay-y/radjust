# script to generate internal datasets. they are used to verify the outputs in the unittest -----------
rm(list = ls())
example(topic = "radjust_sym", package = "radjust", lib.loc = .libPaths())
example(topic = "radjust_pf", package = "radjust", lib.loc = .libPaths())
examples_outputs_for_testing <- mget(ls())
devtools::use_data(examples_outputs_for_testing, overwrite = T, internal = T)

