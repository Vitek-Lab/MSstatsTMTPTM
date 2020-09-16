
test_that("groupComparision works", {

  output<-groupComparisonTMTPTM(data.ptm = MSstatsTMTPTM::quant.msstats.ptm,
                                data.protein = MSstatsTMTPTM::quant.msstats.protein)
  expect_equal(test.pairwise, MSstatsTMTPTM::test.pairwise, tolerance=1e-5)

})

test_that("groupComparision handle missing column", {

  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-1]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-2]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-3]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-4]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-5]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-6]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-7]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-8]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.ptm[,-c(1,2)]))

  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-1]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-2]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-3]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-4]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-5]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-6]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-7]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-8]))
  expect_error(groupComparisonTMTPTM(data = MSstatsTMTPTM::quant.msstats.protein[,-c(1,2)]))

})

test_that("groupComparision handle wrong input", {

  expect_error(groupComparisonTMTPTM())
  expect_error(groupComparisonTMTPTM(data.ptm = MSstatsTMTPTM::quant.msstats.ptm,
                                  data.protein = MSstatsTMTPTM::quant.msstats.protein,
                                  contrast.matrix = matrix(c(1,2,3,4),nrow = 2)))
  expect_error(groupComparisonTMTPTM(data.ptm = MSstatsTMTPTM::quant.msstats.ptm,
                                  data.protein = MSstatsTMTPTM::quant.msstats.protein,
                                  remove_norm_channel = "abc"))
  expect_error(groupComparisonTMTPTM(data.ptm = MSstatsTMTPTM::quant.msstats.ptm,
                                  data.protein = MSstatsTMTPTM::quant.msstats.protein,
                                  moderated = "abc"))

})
