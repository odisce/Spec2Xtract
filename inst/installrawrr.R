## Check that rawrr is installed
library(rawrr)
if (!rawrr::.checkDllInMonoPath()) {
  rawrr::installRawFileReaderDLLs()
  rawrr::installRawrrExe()
}
