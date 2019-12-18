# Functions for HGIS data reduction


# Load data

# Calculate per-block fields
calcRaw <- function(data) {
	data %>%
	  mutate(
      ltcorr = CntTotS - CntTotH,
	    cor1412 = doCor1412(he1412, ltcorr, he1312),
	    sig1412 = ,
	    d13c = calcd13c(he1312)
	  )
}
            
# Trim bad blocks

# Add 13C and deadtime correction
doCor1412 <- function(he1412, ltcorr, he1312) {
	he1412 / ltcorr / he1312 ^ 2
}

# Calc internal error for a measurement
calcSig1412 <- function(CntTotH, CntTotS, CntTotGT, cor1412) {
RelErrSq <- (CntTotH - CntTotS) * CntTotH ^ 2 / CntTotS ^ 4 + 
             CntTotH ^ 2 / CntTotGT / CntTotS ^ 2
cor1412 * RelErrSq ^ 0.5
}

# Calculate d13C
calcd13c <- function(he1312) {
	1000 * (he1312 / 1.12372 -1)
}

## Normalize
# Find standards


# Find mean of stds
normStds <- function(cor1412std, defstd) {
  # both vectors, same length or 1 element vector if all standards are same
  mean(cor1412std/defstd) 
}

# Normalize to mean of standards
norm1412 <- function(cor1412, meanstd) {
  cor1412/meanstd
}

# Propagate normalization error

## Blank correction
# Apply large blank
doLBC <- function(fmmeas, fmblank, fmstd) {
	fmmeas - fmblank * (fmstd - fmmeas) / fmstd
}

# Propagate large blank error
doLBCerror <- function(fmmeas, fmblank, fmstd, fmmeaserr, fmblankerr) {
	# why no fmstderr?
	fmmeaserr ^ 2 * (1 + fmblank / fmstd) ^ 2 +
	fmblankerr ^ 2 * ((fmmeas - fmstd) / fmstd) ^ 2
}

# Produce normalized data for wheel per target
