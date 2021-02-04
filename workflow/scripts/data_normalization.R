library(modules)
library(ArrayExpress)
b = import("../../util/ebits/base/")
io = import('../../util/ebits/base/')
ma = import('../../util/ebits/process/microarray')
idmap = import('../../util/ebits/process/idmap')

ACCESSION = commandArgs(TRUE)[1] %or% 'E-GEOD-40266'
OUTFILE = commandArgs(TRUE)[2] %or% paste0(ACCESSION, ".RData")

# read raw data, normalize, qc
expr = ArrayExpress::ArrayExpress(ACCESSION, drop=FALSE) %>%
    ma$qc() %>%
    ma$normalize() %>%
    ma$annotate(summarize="hgnc_symbol")

# save result
save(expr, file=OUTFILE)
