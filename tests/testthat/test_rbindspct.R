library(photobiology)

context("binding")

test_that("source_mspct", {

  my.spct <- sun.spct
  my.mspct <- source_mspct(list(sun1 = my.spct, sun2 = my.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.irrad", "spct.idx", "s.q.irrad"))
  expect_equal(is_normalized(my.spct), is_normalized(long.spct))
  expect_equal(is_scaled(my.spct), is_scaled(long.spct))
  expect_equal(getTimeUnit(my.spct), getTimeUnit(long.spct))
  expect_equal(is_effective(my.spct), is_effective(long.spct))
  expect_equal(is_tagged(my.spct), is_tagged(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTimeUnit(spct), getTimeUnit(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(is_tagged(spct), is_tagged(long.spct))
  }

  my.spct <- normalize(sun.spct)
  my.mspct <- source_mspct(list(sun1 = my.spct, sun2 = my.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.irrad", "spct.idx"))
  expect_equal(is_normalized(my.spct), is_normalized(long.spct))
  expect_equal(is_scaled(my.spct), is_scaled(long.spct))
  expect_equal(getTimeUnit(my.spct), getTimeUnit(long.spct))
  expect_equal(is_effective(my.spct), is_effective(long.spct))
  expect_equal(is_tagged(my.spct), is_tagged(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTimeUnit(spct), getTimeUnit(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(is_tagged(spct), is_tagged(long.spct))
  }

  my.spct <- fscale(sun.spct)
  my.mspct <- source_mspct(list(sun1 = my.spct, sun2 = my.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.irrad", "spct.idx"))
  expect_equal(is_normalized(my.spct), is_normalized(long.spct))
  expect_equal(is_scaled(my.spct), is_scaled(long.spct))
  expect_equal(getTimeUnit(my.spct), getTimeUnit(long.spct))
  expect_equal(is_effective(my.spct), is_effective(long.spct))
  expect_equal(is_tagged(my.spct), is_tagged(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTimeUnit(spct), getTimeUnit(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(is_tagged(spct), is_tagged(long.spct))
  }

  my.spct <- tag(sun.spct)
  my.mspct <- source_mspct(list(sun1 = my.spct, sun2 = my.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.irrad", "spct.idx", "s.q.irrad"))
  expect_equal(is_normalized(my.spct), is_normalized(long.spct))
  expect_equal(is_scaled(my.spct), is_scaled(long.spct))
  expect_equal(getTimeUnit(my.spct), getTimeUnit(long.spct))
  expect_equal(is_effective(my.spct), is_effective(long.spct))
  expect_true(is_tagged(my.spct))
  expect_true(!is_tagged(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTimeUnit(spct), getTimeUnit(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_true(!is_tagged(spct))
  }

})


test_that("response_mspct", {

  my.spct <- photodiode.spct
  my.mspct <- response_mspct(list(resp1 = my.spct, resp2 = my.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length",  "s.e.response", "spct.idx"))
  expect_equal(is_normalized(my.spct), is_normalized(long.spct))
  expect_equal(is_scaled(my.spct), is_scaled(long.spct))
  expect_equal(getTimeUnit(my.spct), getTimeUnit(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTimeUnit(spct), getTimeUnit(long.spct))
  }

  my.spct <- normalize(photodiode.spct)
  my.mspct <- response_mspct(list(resp1 = my.spct, resp2 = my.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.response", "spct.idx"))
  expect_equal(is_normalized(my.spct), is_normalized(long.spct))
  expect_equal(is_scaled(my.spct), is_scaled(long.spct))
  expect_equal(getTimeUnit(my.spct), getTimeUnit(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTimeUnit(spct), getTimeUnit(long.spct))
  }

  my.spct <- fscale(photodiode.spct)
  my.mspct <- response_mspct(list(resp1 = my.spct, resp2 = my.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.response", "spct.idx"))
  expect_equal(is_normalized(my.spct), is_normalized(long.spct))
  expect_equal(is_scaled(my.spct), is_scaled(long.spct))
  expect_equal(getTimeUnit(my.spct), getTimeUnit(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTimeUnit(spct), getTimeUnit(long.spct))
  }

  my.spct <- tag(photodiode.spct)
  my.mspct <- response_mspct(list(resp1 = my.spct, resp2 = my.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.response", "spct.idx"))
  expect_equal(is_normalized(my.spct), is_normalized(long.spct))
  expect_equal(is_scaled(my.spct), is_scaled(long.spct))
  expect_equal(getTimeUnit(my.spct), getTimeUnit(long.spct))
  expect_equal(is_effective(my.spct), is_effective(long.spct))
  expect_true(is_tagged(my.spct))
  expect_true(!is_tagged(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTimeUnit(spct), getTimeUnit(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_true(!is_tagged(spct))
  }
})

test_that("filter_mspct", {

  flt.spct <- polyester.spct
  my.mspct <- filter_mspct(list(flt1 = flt.spct, flt2 = flt.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Tfr", "spct.idx"))
  expect_equal(is_normalized(flt.spct), is_normalized(long.spct))
  expect_equal(is_scaled(flt.spct), is_scaled(long.spct))
  expect_equal(getTfrType(flt.spct), getTfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTfrType(spct), getTfrType(long.spct))
  }

  flt.spct <- normalize(polyester.spct)
  my.mspct <- filter_mspct(list(flt1 = flt.spct, flt2 = flt.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Tfr", "spct.idx"))
  expect_equal(is_normalized(flt.spct), is_normalized(long.spct))
  expect_equal(is_scaled(flt.spct), is_scaled(long.spct))
  expect_equal(getTfrType(flt.spct), getTfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTfrType(spct), getTfrType(long.spct))
  }

  flt.spct <- fscale(polyester.spct, f = max)
  my.mspct <- filter_mspct(list(flt1 = flt.spct, flt2 = flt.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Tfr", "spct.idx"))
  expect_equal(is_normalized(flt.spct), is_normalized(long.spct))
  expect_equal(is_scaled(flt.spct), is_scaled(long.spct))
  expect_equal(getTfrType(flt.spct), getTfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getTfrType(spct), getTfrType(long.spct))
  }

  flt.spct <- tag(photodiode.spct)
  my.mspct <- response_mspct(list(flt1 = flt.spct, flt2 = flt.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.response", "spct.idx"))
  expect_equal(is_normalized(flt.spct), is_normalized(long.spct))
  expect_equal(is_scaled(flt.spct), is_scaled(long.spct))
  expect_equal(getTfrType(flt.spct), getTfrType(long.spct))
  expect_true(is_tagged(flt.spct))
  expect_true(!is_tagged(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(getTfrType(spct), getTfrType(long.spct))
    expect_true(!is_tagged(spct))
  }
})

test_that("reflector_mspct", {

  # "conversion" retains all data columns
  rfl.spct <- as.reflector_spct(white_body.spct)
  my.mspct <- reflector_mspct(list(flt1 = rfl.spct, flt2 = rfl.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Rfr", "Tfr", "spct.idx"))
  expect_equal(is_normalized(rfl.spct), is_normalized(long.spct))
  expect_equal(is_scaled(rfl.spct), is_scaled(long.spct))
  expect_equal(getRfrType(rfl.spct), getRfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getRfrType(spct), getRfrType(long.spct))
  }

  # "conversion" retains all data columns
  rfl.spct <- normalize(as.reflector_spct(white_body.spct))
  my.mspct <- reflector_mspct(list(flt1 = rfl.spct, flt2 = rfl.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Rfr", "Tfr", "spct.idx"))
  expect_equal(is_normalized(rfl.spct), is_normalized(long.spct))
  expect_equal(is_scaled(rfl.spct), is_scaled(long.spct))
  expect_equal(getRfrType(rfl.spct), getRfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getRfrType(spct), getRfrType(long.spct))
  }

  rfl.spct <- fscale(as.reflector_spct(white_body.spct), f = max)
  my.mspct <- reflector_mspct(list(flt1 = rfl.spct, flt2 = rfl.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Rfr", "Tfr", "spct.idx"))
  expect_equal(is_normalized(rfl.spct), is_normalized(long.spct))
  expect_equal(is_scaled(rfl.spct), is_scaled(long.spct))
  expect_equal(getRfrType(rfl.spct), getRfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(getRfrType(spct), getRfrType(long.spct))
  }

  rfl.spct <- tag(photodiode.spct)
  my.mspct <- response_mspct(list(flt1 = rfl.spct, flt2 = rfl.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "s.e.response", "spct.idx"))
  expect_equal(is_normalized(rfl.spct), is_normalized(long.spct))
  expect_equal(is_scaled(rfl.spct), is_scaled(long.spct))
  expect_equal(getTfrType(rfl.spct), getTfrType(long.spct))
  expect_true(is_tagged(rfl.spct))
  expect_true(!is_tagged(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(getRfrType(spct), getRfrType(long.spct))
    expect_true(!is_tagged(spct))
  }

  })

test_that("object_mspct", {

  obj.spct <- white_body.spct
  my.mspct <- object_mspct(list(flt1 = obj.spct, flt2 = obj.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Rfr", "Tfr", "spct.idx"))
  expect_equal(is_normalized(obj.spct), is_normalized(long.spct))
  expect_equal(is_scaled(obj.spct), is_scaled(long.spct))
  expect_equal(getRfrType(obj.spct), getRfrType(long.spct))
  expect_equal(getTfrType(obj.spct), getTfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(getRfrType(spct), getRfrType(long.spct))
    expect_equal(getTfrType(spct), getTfrType(long.spct))
    expect_true(!is_tagged(spct))
  }

  # there is no normalize method for object_spct
  expect_error(normalize(obj.spct))
  my.mspct <- object_mspct(list(flt1 = obj.spct, flt2 = obj.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Rfr", "Tfr", "spct.idx"))
  expect_equal(is_normalized(obj.spct), is_normalized(long.spct))
  expect_equal(is_scaled(obj.spct), is_scaled(long.spct))
  expect_equal(getRfrType(obj.spct), getRfrType(long.spct))
  expect_equal(getTfrType(obj.spct), getTfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(getRfrType(spct), getRfrType(long.spct))
    expect_equal(getTfrType(spct), getTfrType(long.spct))
    expect_true(!is_tagged(spct))
  }

  # there is no fscale method for object_spct
  expect_error(fscale(white_body.spct, f = max))
  my.mspct <- object_mspct(list(flt1 = obj.spct, flt2 = obj.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Rfr", "Tfr", "spct.idx"))
  expect_equal(is_normalized(obj.spct), is_normalized(long.spct))
  expect_equal(is_scaled(obj.spct), is_scaled(long.spct))
  expect_equal(getRfrType(obj.spct), getRfrType(long.spct))
  expect_equal(getTfrType(obj.spct), getTfrType(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(getRfrType(spct), getRfrType(long.spct))
    expect_equal(getTfrType(spct), getTfrType(long.spct))
    expect_true(!is_tagged(spct))
  }

  objt.spct <- tag(white_body.spct)
  my.mspct <- object_mspct(list(flt1 = obj.spct, flt2 = obj.spct))
  long.spct <- rbindspct(my.mspct)
  expect_named(long.spct, c("w.length", "Rfr", "Tfr", "spct.idx"))
  expect_equal(is_normalized(obj.spct), is_normalized(long.spct))
  expect_equal(is_scaled(obj.spct), is_scaled(long.spct))
  expect_equal(getRfrType(obj.spct), getRfrType(long.spct))
  expect_equal(getTfrType(obj.spct), getTfrType(long.spct))
  expect_true(is_tagged(objt.spct))
  expect_true(!is_tagged(long.spct))

  my_new.mspct <- subset2mspct(long.spct)
  for (spct in my_new.mspct) {
    expect_equal(is_normalized(spct), is_normalized(long.spct))
    expect_equal(is_scaled(spct), is_scaled(long.spct))
    expect_equal(is_effective(spct), is_effective(long.spct))
    expect_equal(getRfrType(spct), getRfrType(long.spct))
    expect_equal(getTfrType(spct), getTfrType(long.spct))
    expect_true(!is_tagged(spct))
  }

})## need to add similar tests for other classes
## need to add additional methods
