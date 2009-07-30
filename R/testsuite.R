runTests <- function() {
  require(RUnit)
  testsuite.GEOquery <- defineTestSuite("GEOquery",
                                        dir=system.file('test',package='GEOquery'))
  tr <- runTestSuite(testsuite.GEOquery)
  printTextProtocol(tr)
}
