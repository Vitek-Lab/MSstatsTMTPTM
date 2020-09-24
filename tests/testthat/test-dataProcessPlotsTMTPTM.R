test_that("dataProcessPlotsTMTPTM handle missing column", {

  expect_error(dataProcessPlotsTMTPTM(MSstatsTMTPTM::raw.ptm[,-1],
                                      MSstatsTMTPTM::raw.protein,
                                      MSstatsTMTPTM::quant.msstats.ptm,
                                      MSstatsTMTPTM::quant.msstats.protein,
                                      type = 'ProfilePlot'))
  expect_error(dataProcessPlotsTMTPTM(MSstatsTMTPTM::raw.ptm,
                                      MSstatsTMTPTM::raw.protein[,-1],
                                      MSstatsTMTPTM::quant.msstats.ptm,
                                      MSstatsTMTPTM::quant.msstats.protein,
                                      type = 'ProfilePlot'))
  expect_error(dataProcessPlotsTMTPTM(MSstatsTMTPTM::raw.ptm,
                                      MSstatsTMTPTM::raw.protein,
                                      MSstatsTMTPTM::quant.msstats.ptm[,-1],
                                      MSstatsTMTPTM::quant.msstats.protein,
                                      type = 'ProfilePlot'))
  expect_error(dataProcessPlotsTMTPTM(MSstatsTMTPTM::raw.ptm,
                                      MSstatsTMTPTM::raw.protein,
                                      MSstatsTMTPTM::quant.msstats.ptm,
                                      MSstatsTMTPTM::quant.msstats.protein[,-1],
                                      type = 'ProfilePlot'))

})

test_that("dataProcessPlotsTMTPTM handle wrong input", {

  expect_error(dataProcessPlotsTMTPTM())
  expect_error(dataProcessPlotsTMTPTM(MSstatsTMTPTM::raw.ptm,
                                      MSstatsTMTPTM::raw.protein,
                                      MSstatsTMTPTM::quant.msstats.ptm,
                                      MSstatsTMTPTM::quant.msstats.protein,
                                      type = 'abc'))
  expect_error(dataProcessPlotsTMTPTM(MSstatsTMTPTM::raw.ptm,
                                      MSstatsTMTPTM::raw.protein,
                                      MSstatsTMTPTM::quant.msstats.ptm,
                                      MSstatsTMTPTM::quant.msstats.protein,
                                      type = 'ProfilePlot',
                                      which.Protein = 'abc'))
  expect_error(dataProcessPlotsTMTPTM(MSstatsTMTPTM::raw.ptm,
                                      MSstatsTMTPTM::raw.protein,
                                      MSstatsTMTPTM::quant.msstats.ptm,
                                      MSstatsTMTPTM::quant.msstats.protein,
                                      type = 'ProfilePlot',
                                      originalPlot = 'abc'))

})
