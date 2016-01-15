source("refund_lib.R")
# Read methylation data
meth.27K.QC <- as.matrix(read.table("meth_27K_QC.txt", check.names = F))
meth.450K.QC <- as.matrix(read.table("meth_450K_QC.txt", check.names = F))

# Tranform the methylation value from beta value to M value
meth.27K.QC <- log2(meth.27K.QC / (1 - meth.27K.QC))
meth.450K.QC <- log2(meth.450K.QC / (1 - meth.450K.QC))

####################
#### Annotation ####
####################

# Read annotation file
annotation.27K <- read.csv("annotation/GPL8490-65.csv.gz", comment.char = "#")
annotation.450K <- read.csv("annotation/GPL13534-10305.csv.gz", comment.char = "#")

annotation.27K <- annotation.27K[(annotation.27K[, 1] %in% rownames(meth.27K.QC)) & (annotation.27K[, 1] %in% rownames(meth.450K.QC)),]
annotation.450K <- annotation.450K[annotation.450K[, 1] %in% rownames(meth.450K.QC),]

# Sort the annotation and methylation data according to probe location
chr <- paste(annotation.27K$Chr)
chr[chr == "X"] <- 23
chr[chr == "Y"] <- 24
chr <- as.integer(chr)
pos <- annotation.27K$MapInfo
order <- order(chr, pos)
annotation.27K <- annotation.27K[order,]
meth.27K.QC <- meth.27K.QC[paste(annotation.27K[, 1]),]

chr <- paste(annotation.450K$CHR)
chr[chr == "X"] <- 23
chr[chr == "Y"] <- 24
chr <- as.integer(chr)
pos <- annotation.450K$MAPINFO
order <- order(chr, pos)
annotation.450K <- annotation.450K[order,]
meth.450K.QC <- meth.450K.QC[paste(annotation.450K[, 1]),]
rownames(annotation.450K) <- annotation.450K[, 1]

annotation.27K$Relation_to_UCSC_CpG_Island <- annotation.450K[rownames(meth.27K.QC), "Relation_to_UCSC_CpG_Island"]

# Choose a set of probes showing large variation for imputation
nprobe <- 1000
meth.var <- apply(2 ^ meth.450K.QC / (2 ^ meth.450K.QC + 1), 1, var)
meth.var <- sort(meth.var, decreasing = T)
meth.var <- meth.var[!(names(meth.var) %in% rownames(meth.27K.QC))]
probe <- names(meth.var)[1:nprobe]
# probe <- names(meth.var)

isl_group <- c("", "Island", "N_Shelf", "N_Shore", "S_Shore", "S_Shelf")
probe.na <- (annotation.450K[probe, "CHR"] == "") | (is.na(annotation.450K[probe, "MAPINFO"])) | (!(annotation.450K[probe, "Relation_to_UCSC_CpG_Island"] %in% isl_group))
probe <- probe[!probe.na]
nprobe <- length(probe)

cov.num <- 10
meth.impute <- matrix(NA, nprobe, ncol(meth.27K.QC))
rownames(meth.impute) <- probe
colnames(meth.impute) <- colnames(meth.27K.QC)
dispersion <- rep(NA, nprobe)
names(dispersion) <- probe

# Build the functional predictor for each group
train.dens <- list()
train.funcs <- list()
test.dens <- list()
test.funcs <- list()

for (i in isl_group) {
  if (i == "") {
    j <- "NA"
  } else {
    j <- i
  }
  train.dens[[j]] <- meth.450K.QC[(annotation.450K[, "Relation_to_UCSC_CpG_Island"] == i) & (rownames(meth.450K.QC) %in% rownames(meth.27K.QC)),]
  train.funcs[[j]] <- apply(train.dens[[j]], 2, function(x) {density(x)$y})
  test.dens[[j]] <- meth.27K.QC[annotation.27K[, "Relation_to_UCSC_CpG_Island"] == i,]
  test.funcs[[j]] <- apply(test.dens[[j]], 2, function(x) {density(x)$y})
}

# Parallelize this for loop can accelerate the imputation
for (i in 1:nprobe) {
  if (i %% 100 == 0) {
    print(paste("Processed ", i, " probes...", sep = ""))
  }
  Y.train <- meth.450K.QC[probe[i],]

  ## locate probe
  chr <- paste(annotation.450K[probe[i], "CHR"])
  pos <- annotation.450K[probe[i], "MAPINFO"]

  # Select local probes
  upstream <- which((annotation.27K$Chr == chr) & (annotation.27K$MapInfo < pos))
  downstream <- which((annotation.27K$Chr == chr) & (annotation.27K$MapInfo > pos))
  if (length(upstream) < cov.num / 2) {
    up.r <- min(min(which(annotation.27K$Chr == chr)) + cov.num / 2 - 1, nrow(annotation.27K))
  } else {
    up.r <- max(upstream)
  }
  if (length(downstream) < cov.num / 2) {
    down.l <- max(1, max(which(annotation.27K$Chr == chr)) - cov.num / 2 + 1)
  } else {
    down.l <- min(downstream)
  }

  cov.probe <- paste(annotation.27K[unique(c(max(1, up.r - cov.num / 2 + 1):up.r, down.l:min(nrow(annotation.27K), down.l + cov.num / 2 - 1))), "ID"])
  cov.train <- meth.450K.QC[cov.probe,]
  cov.test <- meth.27K.QC[cov.probe,]

  # Model: Functional predictor + Local covariates
  j <- paste(annotation.450K[probe[i], "Relation_to_UCSC_CpG_Island"])
  if (j == "") {
    j <- "NA"
  }
  fit = new_pfr(Y.train, t(cov.train), t(train.funcs[[j]]))
  Y.fit <- fit$fitted.vals
  Y.fit <- 2 ^ Y.fit / (2 ^ Y.fit + 1)

  X.test <- create_predictors(t(cov.test), t(test.funcs[[j]]))
  Y.test <- X.test %*% fit$coefs
  Y.test <- 2 ^ Y.test / (2 ^ Y.test + 1)
  meth.impute[probe[i],] <- Y.test
  dispersion[i] <- var(Y.fit) / var(2 ^ Y.train / (2 ^ Y.train + 1))
}

write.table(meth.impute, "impute.txt", quote = F, sep = "\t")
write.table(dispersion, "dispersion.txt", quote = F, col.names = F, sep = "\t")