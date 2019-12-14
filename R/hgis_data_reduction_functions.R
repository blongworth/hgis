# Functions for HGIS data reduction


# Load data

# Calculate per-block fields

calcRaw <- function(data) {
	data %>%
	  mutate(
            ltcorr = CntTotS - CntTotH,
	    cor1412 = doCorr(he1412, ltcorr, he1312),
	    sig1412 = ,
	    d13c = calcd13c(he1312)
	  )
}
            
# Trim bad blocks

# Add 13C and deadtime correction
doCorr <- function(he1412, ltcorr, he1312) {
	he1412 / ltcorr / he1312 ^ 2
}

# Calculate d13C
calcd13c <- function(he1312) {
	1000 * (he1312 / 1.12372 -1)
}

## Normalize
# Find standards
# Normalize to mean of standards

# Propagate error

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
