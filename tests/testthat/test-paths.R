test_that("replication paths resolve packaged files", {
  expect_true(file.exists(replication_file("DESCRIPTION")))
  expect_true(file.exists(replication_file(
    "data", "empirics", "ff", "ff_2x3_sorts_daily.rds"
  )))
})
