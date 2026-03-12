library("photobiology")

context("trim_wl")

test_that("trimming works with wavebands", {

  my.spct <- source_spct(w.length = 400:450, s.e.irrad = 0.5, time.unit = "second")
  my.wb <- waveband(c(300,500))
  expect_equal(expanse(trim_wl(my.wb, my.spct)), expanse(my.spct))
  expect_equal(wl_range(trim_wl(my.wb, my.spct)), wl_range(my.spct))

  # NA limits in waveband!
  my.wb <- waveband(c(100, 150))
  expect_equal(trim_wl(my.wb, my.spct), waveband())

  # nested trimming
  my.wb <- waveband(c(100, 1000))
  expect_equal(trim_wl(my.wb, c(400, 700)),
               trim_wl(trim_wl(my.wb, c(300, 800)), c(400, 700)))

  # waveband as range argument
  expect_equal(trim_wl(my.wb, c(400, 700)),
               trim_wl(my.wb, range = waveband(c(400, 700))))

  trimmed.wb <- trim_waveband(my.wb, c(NA, 700))
  expect_equal(wl_range(trimmed.wb), c(100, 700))
  expect_equal(labels(trimmed.wb)$label, "range.100.1000[")
  expect_equal(labels(trimmed.wb)$name, "range.100.1000[")

  trimmed.wb <- trim_waveband(my.wb, high.limit = 700)
  expect_equal(wl_range(trimmed.wb), c(100, 700))
  expect_equal(labels(trimmed.wb)$label, "range.100.1000[")
  expect_equal(labels(trimmed.wb)$name, "range.100.1000[")

  trimmed.wb <- trim_waveband(my.wb, c(400, NA))
  expect_equal(wl_range(trimmed.wb), c(400, 1000))
  expect_equal(labels(trimmed.wb)$label, "]range.100.1000")
  expect_equal(labels(trimmed.wb)$name, "]range.100.1000")

  trimmed.wb <- trim_waveband(my.wb, low.limit = 400)
  expect_equal(wl_range(trimmed.wb), c(400, 1000))
  expect_equal(labels(trimmed.wb)$label, "]range.100.1000")
  expect_equal(labels(trimmed.wb)$name, "]range.100.1000")

  trimmed.wb <- trim_waveband(my.wb, c(400, 700))
  expect_equal(wl_range(trimmed.wb), c(400, 700))
  expect_equal(labels(trimmed.wb)$label, "]range.100.1000[")
  expect_equal(labels(trimmed.wb)$name, "]range.100.1000[")

  trimmed.wb <- trim_waveband(my.wb, c(400, 700), trunc.labels = c(">", "<"))
  expect_equal(wl_range(trimmed.wb), c(400, 700))
  expect_equal(labels(trimmed.wb)$label, ">range.100.1000<")
  expect_equal(labels(trimmed.wb)$name, ">range.100.1000<")

  trimmed.wb <- trim_waveband(my.wb, c(400, 700), trunc.labels = "!")
  expect_equal(wl_range(trimmed.wb), c(400, 700))
  expect_equal(labels(trimmed.wb)$label, "!range.100.1000!")
  expect_equal(labels(trimmed.wb)$name, "!range.100.1000!")

})

test_that("trimming works with list of wavebands", {

  my.spct <- source_spct(w.length = 400:450, s.e.irrad = 0.5, time.unit = "second")
  my.partnamed.wbs.ls <- my.named.wbs.ls <- my.wbs.ls <-
    list(waveband(c(300,500)), waveband(c(400,600)), waveband(c(300,1000)))
  names(my.named.wbs.ls) <- c("first", "second", "third")
  names(my.partnamed.wbs.ls)[c(1,3)] <- c("first", "third")

  expect_equal(unique(sapply(X = trim_wl(my.wbs.ls, my.spct), FUN = expanse)),
               expanse(my.spct))
  expect_equal(range(sapply(X = trim_wl(my.wbs.ls, my.spct), FUN = wl_range)),
               wl_range(my.spct))

  expect_equal(unique(sapply(X = trim_wl(my.named.wbs.ls, my.spct), FUN = expanse)),
               expanse(my.spct))
  expect_equal(range(sapply(X = trim_wl(my.named.wbs.ls, my.spct), FUN = wl_range)),
               wl_range(my.spct))

  expect_equal(unique(sapply(X = trim_wl(my.partnamed.wbs.ls, my.spct), FUN = expanse)),
               expanse(my.spct))
  expect_equal(range(sapply(X = trim_wl(my.partnamed.wbs.ls, my.spct), FUN = wl_range)),
               wl_range(my.spct))

  expect_named(trim_wl(my.wbs.ls, my.spct), NULL)
  expect_named(trim_wl(my.named.wbs.ls, my.spct), names(my.named.wbs.ls))
  expect_named(trim_wl(my.partnamed.wbs.ls, my.spct), names(my.partnamed.wbs.ls))
})
