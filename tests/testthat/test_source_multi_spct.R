library(photobiology)

test.print <- FALSE

context("multi_spct")

test_that("constructors", {

  my.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct))
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("sun1", "sun2"))
  expect_true(is.source_spct(my.mspct[["sun1"]]))
  expect_true(is.source_spct(my.mspct[["sun2"]]))

  my.mspct <- as.generic_mspct(my.mspct)
  expect_false(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("sun1", "sun2"))
  expect_true(is.source_spct(my.mspct[["sun1"]]))
  expect_true(is.source_spct(my.mspct[["sun2"]]))
  expect_true(is.any_spct(my.mspct[["sun1"]]))
  expect_true(is.any_spct(my.mspct[["sun2"]]))

  my.mspct <- as.generic_mspct(my.mspct, force.spct.class = TRUE)
  expect_true(is.generic_mspct(my.mspct))
  expect_false(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("sun1", "sun2"))
  expect_false(is.source_spct(my.mspct[["sun1"]]))
  expect_false(is.source_spct(my.mspct[["sun2"]]))
  expect_true(is.any_spct(my.mspct[["sun1"]]))
  expect_true(is.any_spct(my.mspct[["sun2"]]))
  expect_true(is.generic_spct(my.mspct[["sun1"]]))
  expect_true(is.generic_spct(my.mspct[["sun2"]]))

  my.mspct <- as.source_mspct(my.mspct)
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("sun1", "sun2"))
  expect_true(is.source_spct(my.mspct[["sun1"]]))
  expect_true(is.source_spct(my.mspct[["sun2"]]))

  my.mspct <- as.source_mspct(sun.spct)
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.source_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.source_mspct(data.frame(w.length = 400:500, s.e.irrad = 0.1))
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.source_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.source_mspct(tibble::tibble(w.length = 400:500, s.e.irrad = 0.1))
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1"))
  expect_true(is.source_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.source_mspct(list(data.frame(w.length = 400:500, s.e.irrad = 0.1),
                                   data.frame(w.length = 400:500, s.e.irrad = 0.2)))
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("spct_1", "spct_2"))
  expect_true(is.source_spct(my.mspct[["spct_1"]]))
  expect_true(is.source_spct(my.mspct[["spct_1"]]))

  my.mspct <- as.source_mspct(list(A = data.frame(w.length = 400:500, s.e.irrad = 0.1),
                                   B = data.frame(w.length = 400:500, s.e.irrad = 0.2)))
  expect_true(is.source_mspct(my.mspct))
  expect_true(is.any_mspct(my.mspct))
  expect_named(my.mspct, c("A", "B"))
  expect_true(is.source_spct(my.mspct[["A"]]))
  expect_true(is.source_spct(my.mspct[["B"]]))

  expect_message(as.source_mspct(1))
  expect_message(as.source_mspct("abc"))
  expect_message(as.source_mspct(TRUE))
  expect_error(as.source_mspct(list(w.length = 400:500,
                                    s.e.irrad = rep(0.3, 101))))

  expect_error(as.filter_spct(my.mspct))
  expect_error(as.reflector_spct(my.mspct))
  expect_error(as.object_spct(my.mspct))
  expect_error(as.response_spct(my.mspct))
  expect_error(as.cps_spct(my.mspct))
  expect_error(as.raw_spct(my.mspct))

  empty.mspct <- source_mspct()
  expect_true(is.source_mspct(empty.mspct))
  expect_true(is.any_mspct(empty.mspct))
  expect_true(is.null(names(empty.mspct)))

  })

test_that("source_mspct", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  my3.spct <- source_spct(w.length = 400:410, s.e.irrad = 3)
  my4.spct <- source_spct(w.length = 400:410, s.e.irrad = 4)
  my5.spct <- source_spct(w.length = 400:410, s.e.irrad = 5)

  spct.l <- list(my1.spct, my2.spct, my3.spct, my4.spct, my5.spct)
  my.mspct <- source_mspct(spct.l)

  expect_equal(paste("spct", seq_len(length(spct.l)), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct,
                       three = my3.spct,
                       four = my4.spct,
                       five = my5.spct)

  my_named.mspct <- source_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 5)
  expect_equal(length(my_named.mspct), 5)

  expect_true(is.data.frame(irrad(my.mspct)))
  expect_true(is.data.frame(irrad(my_named.mspct)))

  expect_equal(irrad(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(irrad(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))

  expect_equal(round(irrad(my.mspct)[["E_Total"]], 3), 1:5 * 10)
  expect_equal(round(irrad(my_named.mspct)[["E_Total"]], 3), 1:5 * 10)

# irrad -------------------------------------------------------------------

  expect_equal(
    round(irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B"))
    )[["E_A"]], 3), 1:5 * 5
  )
  expect_equal(
    round(irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B"))
    )[["E_B"]], 3), 1:5 * 5
  )

  expect_equal(round(irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["E_A"]], 3), 1:5 * 5)
  expect_equal(round(irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["E_B"]], 3), 1:5 * 5)

  expect_equal(round(irrad(my_named.mspct,
                           list(A = waveband(c(400,405)),
                                B = waveband(c(405,410)))
  )[["E_A"]], 3), 1:5 * 5)
  expect_equal(round(irrad(my_named.mspct,
                           list(A = waveband(c(400,405)),
                                B = waveband(c(405,410)))
  )[["E_B"]], 3), 1:5 * 5)

  expect_equal(round(irrad(my_named.mspct,
                           list(waveband(c(400,405)),
                                waveband(c(405,410)))
  )[["E_range.400.405"]], 3), 1:5 * 5)
  expect_equal(round(irrad(my_named.mspct,
                           list(waveband(c(400,405)),
                                waveband(c(405,410)))
  )[["E_range.405.410"]], 3), 1:5 * 5)

  # relative

  expect_equal(
    round(irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["E/Esum_A"]], 3), rep(0.5, 5)
  )
  expect_equal(
    round(irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["E/Esum_B"]], 3), rep(0.5, 5)
  )

  expect_equal(
    round(irrad(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["E/Esum_A"]], 3), rep(0.5, 5)
  )
  expect_equal(
    round(irrad(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["E/Esum_B"]], 3), rep(0.5, 5)
  )

  # e_irrad -----------------------------------------------------------------

  expect_equal(
    round(e_irrad(my.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B"))
    )[["E_A"]], 3), 1:5 * 5
  )
  expect_equal(
    round(e_irrad(my.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B"))
    )[["E_B"]], 3), 1:5 * 5
  )

  expect_equal(round(e_irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["E_A"]], 3), 1:5 * 5)
  expect_equal(round(e_irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
  )[["E_B"]], 3), 1:5 * 5)

  # relative

  expect_equal(
    round(e_irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["E/Esum_A"]], 3), rep(0.5, 5)
  )
  expect_equal(
    round(e_irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["E/Esum_B"]], 3), rep(0.5, 5)
  )

  expect_equal(
    round(e_irrad(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["E/Esum_A"]], 3), rep(0.5, 5)
  )
  expect_equal(
    round(e_irrad(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["E/Esum_B"]], 3), rep(0.5, 5)
  )

  # contribution --------------------------------------------------------------

  expect_equal(
    round(e_irrad(my_named.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B")),
                  quantity = "contribution"
    )[["E/Etot_A"]], 3), rep(0.5, 5)
  )
  expect_equal(
    round(e_irrad(my_named.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B")),
                  quantity = "contribution"
    )[["E/Etot_B"]], 3), rep(0.5, 5)
  )

  # q_irrad -------------------------------------------------------------------

  expect_equal(
    round(q_irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B"))
    )[["Q_A"]], 7), c(1.68e-5, 3.36e-5,5.05e-5, 6.73e-5,8.41e-5)
  )
  expect_equal(
    round(q_irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B"))
    )[["Q_B"]], 7), c(1.70e-5,3.41e-5,5.11e-5,6.81e-5,8.52e-5)
  )

  expect_equal(round(q_irrad(my_named.mspct,
                             list(waveband(c(400,405), wb.name = "A"),
                                  waveband(c(405,410), wb.name = "B"))
    )[["Q_A"]], 7), c(1.68e-5, 3.36e-5,5.05e-5, 6.73e-5,8.41e-5)
  )
  expect_equal(round(q_irrad(my_named.mspct,
                           list(waveband(c(400,405), wb.name = "A"),
                                waveband(c(405,410), wb.name = "B"))
    )[["Q_B"]], 7), c(1.70e-5,3.41e-5,5.11e-5,6.81e-5,8.52e-5)
  )

  # relative

  expect_equal(
    round(q_irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["Q/Qsum_A"]], 3), rep(0.497, 5)
  )
  expect_equal(
    round(q_irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["Q/Qsum_B"]], 3), rep(0.503, 5)
  )

  expect_equal(
    round(q_irrad(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["Q/Qsum_A"]], 3), rep(0.497, 5)
  )
  expect_equal(
    round(q_irrad(my_named.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                quantity = "relative"
    )[["Q/Qsum_B"]], 3), rep(0.503, 5)
  )

  # contribution ------------------------------------------------------------

  expect_equal(
    round(q_irrad(my.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B")),
                  quantity = "contribution"
    )[["Q/Qtot_A"]], 3), rep(0.497, 5)
  )
  expect_equal(
    round(q_irrad(my.mspct,
                  list(waveband(c(400,405), wb.name = "A"),
                       waveband(c(405,410), wb.name = "B")),
                  quantity = "contribution"
    )[["Q/Qtot_B"]], 3), rep(0.503, 5)
  )

  # ratios ------------------------------------------------------------------

  expect_equal(
    round(q_ratio(my.mspct,
                  waveband(c(400,405), wb.name = "A"),
                  waveband(c(405,410), wb.name = "B")
    )[["A:B[q:q]"]], 3), rep(0.988, 5)
  )
  expect_equal(
    round(e_ratio(my.mspct,
                  waveband(c(400,405), wb.name = "A"),
                  waveband(c(405,410), wb.name = "B")
    )[["A:B[e:e]"]], 3), rep(1.000, 5)
  )

  expect_equal(
    round(eq_ratio(my.mspct,
                  waveband(c(400,405), wb.name = "A"
    ))[["A[e:q]"]], -2), rep(297200, 5)
  )
  expect_equal(
    round(qe_ratio(my.mspct,
                  waveband(c(400,405), wb.name = "A"
    ))[["A[q:e]"]], 9), rep(3.365e-06, 5)
  )

  # min ---------------------------------------------------------------------

  expect_equal(min(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:5, sep = "_")))
  expect_equal(levels(max(my_named.mspct)[["spct.idx"]]),
               c("one", "two", "three", "four", "five"))
  expect_equal(min(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(min(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(min(my.mspct), wl_min(my.mspct))
  expect_equal(min(my_named.mspct), wl_min(my_named.mspct))

  # max ---------------------------------------------------------------------

  expect_equal(max(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(max(my_named.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(max(my.mspct), wl_max(my.mspct))
  expect_equal(max(my_named.mspct), wl_max(my_named.mspct))

  # range -------------------------------------------------------------------

  expect_equal(wl_range(my.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_range(my_named.mspct)[["min.wl"]], rep(400, 5))
  expect_equal(wl_range(my.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(wl_range(my_named.mspct)[["max.wl"]], rep(410, 5))
  expect_equal(range(my.mspct), wl_range(my.mspct))
  expect_equal(range(my_named.mspct), wl_range(my_named.mspct))

  # midpoint ----------------------------------------------------------------

  expect_equal(midpoint(my.mspct)[["midpoint.wl"]], rep(405, 5))
  expect_equal(midpoint(my_named.mspct)[["midpoint.wl"]], rep(405, 5))
  expect_equal(midpoint(my.mspct), wl_midpoint(my.mspct))
  expect_equal(midpoint(my_named.mspct), wl_midpoint(my_named.mspct))

  # expanse -----------------------------------------------------------------

  expect_equal(expanse(my.mspct)[["expanse_V1"]], rep(10, 5))
  expect_equal(expanse(my_named.mspct)[["expanse_V1"]], rep(10, 5))
  expect_equal(expanse(my.mspct), wl_expanse(my.mspct))
  expect_equal(expanse(my_named.mspct), wl_expanse(my_named.mspct))
  expect_message(spread(my.mspct[[1]]))

  # constructor methods for 'wide' data frames ------------------------------

  my_wide.df <- data.frame(w.length = 300:400, A = 1, B = 2, C = 3)
  my1_df.mspct <- split2source_mspct(my_wide.df)
  expect_equal(c("A", "B", "C"), names(my1_df.mspct))
  expect_equal(nrow(my_wide.df), nrow(my1_df.mspct[[1]]))
  expect_equal(class(my1_df.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(class(my1_df.mspct[[1]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[2]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my1_df.mspct[[3]])[1:2], c("source_spct", "generic_spct") )

# constructor methods for 'long' data frames --------------------------------

  my_long.df <- data.frame(w.length = rep(300:310, 3),
                           s.e.irrad = c(rep(1, 11), rep(2, 11), rep(3, 11)),
                           spct.idx = c(rep("A", 11), rep("B", 11), rep("C", 11)) )
  my2_df.mspct <- subset2mspct(my_long.df, member.class = "source_spct")

  expect_equal(c("A", "B", "C"), names(my2_df.mspct))
  expect_equal(levels(factor(my_long.df[["spct.idx"]])), names(my2_df.mspct))
  expect_equal(class(my2_df.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(class(my2_df.mspct[[1]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[2]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my2_df.mspct[[3]])[1:2], c("source_spct", "generic_spct") )

  # constructor methods for 'long' spct objects -----------------------------

  my_long.spct <- rbindspct(spct.l)
  my3_df.mspct <- subset2mspct(my_long.spct)

  expect_equal(paste("spct", 1:5, sep = "_"), names(my3_df.mspct))
  expect_equal(class(my3_df.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(class(my3_df.mspct[[1]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[2]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[3]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[4]])[1:2], c("source_spct", "generic_spct") )
  expect_equal(class(my3_df.mspct[[5]])[1:2], c("source_spct", "generic_spct") )

  # print -------------------------------------------------------------------

  if (test.print) expect_equal(print(my.mspct), my.mspct)

  # clean -------------------------------------------------------------------

  expect_equal(clean(my.mspct), my.mspct)
})

test_that("source_mspct_attr", {

  my1.spct <- source_spct(w.length = 400:410, s.e.irrad = 1)
  setWhatMeasured(my1.spct, "first spectrum")
  setWhenMeasured(my1.spct, lubridate::ymd("2018-03-03", tz = "UTC"))
  setWhereMeasured(my1.spct, data.frame(lat = -30, lon = +80))
  my2.spct <- source_spct(w.length = 400:410, s.e.irrad = 2)
  setWhatMeasured(my2.spct, "second spectrum")
  setWhenMeasured(my2.spct, lubridate::ymd_hm("2018-03-03 12:30", tz = "UTC"))
  setWhereMeasured(my2.spct, data.frame(lat = 5, lon = 20))


  spct.l <- list(my1.spct, my2.spct)
  my.mspct <- source_mspct(spct.l)

  expect_equal(paste("spct", seq_len(length(spct.l)), sep = "_"), names(my.mspct))

  expect_equal(class(my.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(attr(my.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my.mspct), 1)
  expect_equal(nrow(my.mspct), length(spct.l))
  expect_equal(attr(my.mspct, "mspct.byrow", exact = TRUE), FALSE)

  named_spct.l <- list(one = my1.spct,
                       two = my2.spct)

  my_named.mspct <- source_mspct(named_spct.l)

  expect_equal(names(named_spct.l), names(my_named.mspct))
  expect_equal(class(my_named.mspct)[1:2], c("source_mspct", "generic_mspct") )
  expect_equal(attr(my_named.mspct, "mspct.version", exact = TRUE), 2)
  expect_equal(ncol(my_named.mspct), 1)
  expect_equal(nrow(my.mspct), length(named_spct.l))
  expect_equal(attr(my_named.mspct, "mspct.byrow", exact = TRUE), FALSE)

  expect_equal(length(my.mspct), 2)
  expect_equal(length(my_named.mspct), 2)

  expect_true(is.data.frame(irrad(my.mspct)))
  expect_true(is.data.frame(irrad(my_named.mspct)))

  expect_equal(irrad(my.mspct)[["spct.idx"]],
               factor(paste("spct", 1:2, sep = "_")))
  expect_equal(levels(irrad(my_named.mspct)[["spct.idx"]]),
               c("one", "two"))

  # irrad -------------------------------------------------------------------

  expect_equal(
    names(
      irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )), c("spct.idx", "E_A", "E_B",
         "when", "what", "latitude", "longitude")
  )

  expect_equal(
    names(
      irrad(my.mspct,
            list(waveband(c(400,405), wb.name = "A"),
                 waveband(c(405,410), wb.name = "B")),
            attr2tb = c(when.measured = "when",
                        what.measured = "what",
                        lat = "latitude",
                        lon = "longitude"),
            idx = "spectrum"
      )), c("spectrum", "E_A", "E_B",
            "when", "what", "latitude", "longitude")
  )

  expect_equal(
    irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
                attr2tb = c(when.measured = "when",
                            what.measured = "what",
                            lat = "latitude",
                            lon = "longitude")
    )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
  )
  expect_equal(
    irrad(my.mspct,
                list(waveband(c(400,405), wb.name = "A"),
                     waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["what"]], c("first spectrum", "second spectrum")
  )
  expect_equal(
    irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["latitude"]], c(-30, 5)
  )
  expect_equal(
    irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["longitude"]], c(80, 20)
  )

  expect_equal(
    irrad(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
  )
  expect_equal(
    irrad(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["what"]], c("first spectrum", "second spectrum")
  )
  expect_equal(
    irrad(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["latitude"]], c(-30, 5)
  )
  expect_equal(
    irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude")
    )[["longitude"]], c(80, 20)
  )

  # relative

  expect_equal(
    irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
  )
  expect_equal(
    irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["what"]], c("first spectrum", "second spectrum")
  )
  expect_equal(
    irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["latitude"]], c(-30, 5)
  )
  expect_equal(
    irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["longitude"]], c(80, 20)
  )

  expect_equal(
    irrad(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["when"]], lubridate::ydm_hm(c("2018-03-03 00:00", "2018-03-03 12:30"))
  )
  expect_equal(
    irrad(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["what"]], c("first spectrum", "second spectrum")
  )
  expect_equal(
    irrad(my_named.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["latitude"]], c(-30, 5)
  )
  expect_equal(
    irrad(my.mspct,
          list(waveband(c(400,405), wb.name = "A"),
               waveband(c(405,410), wb.name = "B")),
          attr2tb = c(when.measured = "when",
                      what.measured = "what",
                      lat = "latitude",
                      lon = "longitude"),
          quantity = "relative"
    )[["longitude"]], c(80, 20)
  )

})



