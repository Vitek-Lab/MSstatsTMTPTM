
# test_that("groupComparisonTMTPTM works", {
#
#   output<-groupComparisonTMTPTM(data.ptm = MSstatsTMTPTM::quant.msstats.ptm,
#                                 data.protein =
#                                   MSstatsTMTPTM::quant.msstats.protein)
#   expect_equal(output[[2]], MSstatsTMTPTM::model.results.pairwise[[2]], tolerance=1e-3)
#   expect_equal(output[[3]], MSstatsTMTPTM::model.results.pairwise[[3]], tolerance=1e-3)
#
# })

test_that("groupComparisonTMTPTM handle missing column", {

  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-1],
                                     MSstatsTMTPTM::quant.msstats.protein))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-2],
                                     MSstatsTMTPTM::quant.msstats.protein))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-3],
                                     MSstatsTMTPTM::quant.msstats.protein))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-4],
                                     MSstatsTMTPTM::quant.msstats.protein))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-5],
                                     MSstatsTMTPTM::quant.msstats.protein))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-6],
                                     MSstatsTMTPTM::quant.msstats.protein))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-7],
                                     MSstatsTMTPTM::quant.msstats.protein))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-8],
                                     MSstatsTMTPTM::quant.msstats.protein))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm[,-c(1,2)],
                                     MSstatsTMTPTM::quant.msstats.protein))

  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[,-1]))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[,-2]))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[,-3]))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[,-4]))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[,-5]))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[,-6]))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[,-7]))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[,-8]))
  expect_error(groupComparisonTMTPTM(MSstatsTMTPTM::quant.msstats.ptm,
                                     MSstatsTMTPTM::quant.msstats.protein[
                                       ,-c(1,2)]))

})

test_that("groupComparisonTMTPTM handle wrong input", {

  expect_error(groupComparisonTMTPTM())
  expect_error(groupComparisonTMTPTM(
    data.ptm = MSstatsTMTPTM::quant.msstats.ptm,
    data.protein = MSstatsTMTPTM::quant.msstats.protein,
    contrast.matrix = matrix(c(1,2,3,4),nrow = 2)))
  expect_error(groupComparisonTMTPTM(
    data.ptm = MSstatsTMTPTM::quant.msstats.ptm,
    data.protein = MSstatsTMTPTM::quant.msstats.protein,
    remove_norm_channel = "abc"))
  expect_error(groupComparisonTMTPTM(
    data.ptm = MSstatsTMTPTM::quant.msstats.ptm,
    data.protein = MSstatsTMTPTM::quant.msstats.protein,
    moderated = "abc"))

})
