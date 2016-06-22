#!/usr/bin/Rscript --silent

#library("reshape")
#library("FME")
options(digits=22)
MAX_DOUBLE <- 1e6

read_params <- function(parmfilename, u_digits = 16, reformat = FALSE) {
	nlines <- 11
	l_digits <- rep(u_digits, length.out = nlines)
	text <- readLines(parmfilename, encoding="UTF-8")
	i <- 2
	i <- i + 1
	Rm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[1])
	ngenes <- length(Rm)
	i <- i + 1
	i <- i + 1
	Tm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[2])
	for ( k in 2:ngenes ) {
		i <- i + 1
		tm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[2])
		Tm <- rbind(Tm, tm)
	}
	i <- i + 1
	i <- i + 1
	Em <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[3])
	for ( k in 2:ngenes ) {
		i <- i + 1
		em <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[3])
		Em <- rbind(Em, em)
	}
	i <- i + 1
	i <- i + 1
	Mm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[4])
	for ( k in 2:ngenes ) {
		i <- i + 1
		mm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[4])
		Mm <- rbind(Mm, mm)
	}
	i <- i + 1
	i <- i + 1
	Pm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[5])
	for ( k in 2:ngenes ) {
		i <- i + 1
		pm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[5])
		Pm <- rbind(Pm, pm)
	}
	i <- i + 1
	i <- i + 1
	hm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[6])
	i <- i + 1
	i <- i + 1
	Dm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[7])
	i <- i + 1
	i <- i + 1
	lambda <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[8])
	i <- i + 1
	i <- i + 1
	tau <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[9])
	i <- i + 1
	i <- i + 1
	ranget <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[10])
	i <- i + 1
	i <- i + 1
	mmm <- round(as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))), digits = l_digits[11])
	if ( reformat ) {
		eqparams <- list(T = t(Tm[1:4,1:4]), E = t(Em[1:4,1:4]), P = t(Pm), h = hm[1:4], tau = tau[1:4], ranget = ranget[1:4], mm = mmm[1:4], lambda = lambda)
	} else {
		eqparams <- list(R = Rm, T = Tm, E = Em, M = Mm, P = Pm, h = hm, D = Dm, tau = tau, m = ranget, mm = mmm, lambda = lambda)
	}
	return(eqparams)
}

params2vec <- function (pars) {
	lparts <- pars$parts
	paramvec <- c()
	k <- 1
	ldparms <- pars$dparms$R
	for ( i in 1:lparts$R ) {
		paramvec <- c(paramvec, ldparms[i])
		k <- k + 1
	}
	ldparms <- pars$dparms$T
	for ( i in 1:lparts$R ) {
		for ( j in 1:lparts$R ) {
			paramvec <- c(paramvec, ldparms[i,j])
			k <- k + 1
		}
	}
	ldparms <- pars$dparms$E
	for ( i in 1:lparts$R ) {
		for ( j in 1:(lparts$E/lparts$R) ) {
			paramvec <- c(paramvec, ldparms[i,j])
			k <- k + 1
		}
	}
	ldparms <- pars$dparms$M
	for ( i in 1:lparts$R ) {
		paramvec <- c(paramvec, ldparms[i,1])
		k <- k + 1
	}
	ldparms <- pars$dparms$P
	for ( i in 1:lparts$R ) {
		for ( j in 1:(lparts$P/lparts$R) ) {
			paramvec <- c(paramvec, ldparms[i,j])
			k <- k + 1
		}
	}
	ldparms <- pars$dparms$h
	for ( i in 1:lparts$h ) {
		paramvec <- c(paramvec, ldparms[i])
		k <- k + 1
	}
	ldparms <- pars$dparms$D
	for ( i in 1:lparts$D ) {
		paramvec <- c(paramvec, ldparms[i])
		k <- k + 1
	}
	ldparms <- pars$dparms$lambda
	for ( i in 1:lparts$lambda ) {
		paramvec <- c(paramvec, ldparms[i])
		k <- k + 1
	}
	ldparms <- pars$dparms$tau
	for ( i in 1:lparts$tau ) {
		paramvec <- c(paramvec, ldparms[i])
		k <- k + 1
	}
	ldparms <- pars$dparms$m
	for ( i in 1:lparts$m ) {
		paramvec <- c(paramvec, ldparms[i])
		k <- k + 1
	}
	ldparms <- pars$dparms$mm
	for ( i in 1:lparts$mm ) {
		paramvec <- c(paramvec, ldparms[i])
		k <- k + 1
	}
	return(paramvec)
}

vec2params <- function (pars, paramvec) {
	lparts <- pars$parts
	ltweak <- unlist(pars$tweak)
	k <- 1
	n <- 1
	ldparms <- pars$dparms
	for ( i in 1:lparts$R ) {
		if ( ltweak[n] == 1 ) {
			ldparms$R[i] <- paramvec[k]
			k <- k + 1
		}
		n <- n + 1
	}
	for ( i in 1:lparts$R ) {
		for ( j in 1:lparts$R ) {
			if ( ltweak[n] == 1 ) {
				ldparms$T[i,j] <- paramvec[k]
				k <- k + 1
			}
			n <- n + 1
		}
	}
	for ( i in 1:lparts$R ) {
		for ( j in 1:(lparts$E/lparts$R) ) {
			if ( ltweak[n] == 1 ) {
				ldparms$E[i,j] <- paramvec[k]
				k <- k + 1
			}
			n <- n + 1
		}
	}
	for ( i in 1:lparts$R ) {
		if ( ltweak[n] == 1 ) {
			ldparms$M[i,1] <- paramvec[k]
			k <- k + 1
		}
		n <- n + 1
	}
	for ( i in 1:lparts$R ) {
		for ( j in 1:(lparts$P/lparts$R) ) {
			if ( ltweak[n] == 1 ) {
				ldparms$P[i,j] <- paramvec[k]
				k <- k + 1
			}
			n <- n + 1
		}
	}
	for ( i in 1:lparts$h ) {
		if ( ltweak[n] == 1 ) {
			ldparms$h[i] <- paramvec[k]
			k <- k + 1
		}
		n <- n + 1
	}
	for ( i in 1:lparts$R ) {
		if ( ltweak[n] == 1 ) {
			ldparms$D[i] <- paramvec[k]
			k <- k + 1
		}
		n <- n + 1
	}
	for ( i in 1:lparts$D ) {
		if ( ltweak[n] == 1 ) {
			ldparms$lambda[i] <- paramvec[k]
			k <- k + 1
		}
		n <- n + 1
	}
	for ( i in 1:lparts$tau ) {
		if ( ltweak[n] == 1 ) {
			ldparms$tau[i] <- paramvec[k]
			k <- k + 1
		}
		n <- n + 1
	}
	for ( i in 1:lparts$m ) {
		if ( ltweak[n] == 1 ) {
			ldparms$m[i] <- paramvec[k]
			k <- k + 1
		}
		n <- n + 1
	}
	for ( i in 1:lparts$mm ) {
		if ( ltweak[n] == 1 ) {
			ldparms$mm[i] <- paramvec[k]
			k <- k + 1
		}
		n <- n + 1
	}
	return(ldparms)
}

sample_params <- function (logfile, optval, eqparms, lbound, hbound, nsteps, pars, datafile, keys, coeffs) {
	fileConn <- file(logfile, open = "a")
	cat("ipar", "istep", "distP", "distF", rep('x', length(eqparms)), 'F', rep('r', length(keys)), "\n", file = fileConn, sep = ",")
	close(fileConn)
	for ( ipar in 1:length(eqparms) ) {
		step <- ( hbound[ipar] - lbound[ipar] ) / nsteps
		param <- lbound[ipar]
		testparms <- eqparms
		currpars <- pars
		for ( istep in 1:( nsteps + 1 ) ) {
			testparms[ipar] <- param
			f <- objfunc(testparms, currpars, datafile, keys, coeffs, buf = T)
			distF <- 100*abs(optval - f$func)/f$func
			distP <- abs(param - eqparms[ipar])
			fileConn <- file(logfile, open = "a")
			cat(ipar, istep, distP, distF, param, f$func, f$res, "\n", file = fileConn, sep = ",")
			close(fileConn)
			param <- param + step
		}
	}
}

write_params <- function (pars, paramfile) {
# write parameters to file
	lparts <- pars$parts
	unlink(paramfile)
	fileConn <- file(paramfile, open = "a")
	k <- 1
	cat("$input", file = fileConn, sep = "\n", append = T)
	cat("R:", file = fileConn, sep = "\n", append = T)
	s1 <- paste("", sep = "")
	ldparms <- pars$dparms$R
	for ( i in 1:lparts$R ) {
		s1 <- paste(s1, ldparms[i], sep = " ")
		k <- k + 1
	}
	cat(s1, file = fileConn, sep = "\n", append = T)
	cat("T:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$T
	for ( i in 1:lparts$R ) {
		s1 <- paste("", sep = "")
		for ( j in 1:lparts$R ) {
			s1 <- paste(s1, ldparms[i,j], sep = " ")
			k <- k + 1
		}
		cat(s1, file = fileConn, sep = "\n", append = T)
	}
	cat("E:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$E
	for ( i in 1:lparts$R ) {
		s1 <- paste("", sep = "")
		for ( j in 1:(lparts$E/lparts$R) ) {
			s1 <- paste(s1, ldparms[i,j], sep = " ")
			k <- k + 1
		}
		cat(s1, file = fileConn, sep = "\n", append = T)
	}
	cat("M:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$M
	for ( i in 1:lparts$R ) {
		s1 <- paste("", sep = "")
		s1 <- paste(s1, ldparms[i,1], sep = " ")
		k <- k + 1
		cat(s1, file = fileConn, sep = "\n", append = T)
	}
	cat("P:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$P
	for ( i in 1:lparts$R ) {
		s1 <- paste("", sep = "")
		for ( j in 1:(lparts$P/lparts$R) ) {
			s1 <- paste(s1, ldparms[i,j], sep = " ")
			k <- k + 1
		}
		cat(s1, file = fileConn, sep = "\n", append = T)
	}
	cat("h:", file = fileConn, sep = "\n", append = T)
	s1 <- paste("", sep = "")
	ldparms <- pars$dparms$h
	for ( i in 1:lparts$R ) {
		s1 <- paste(s1, ldparms[i], sep = " ")
		k <- k + 1
	}
	cat(s1, file = fileConn, sep = "\n", append = T)
	cat("D:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$D
	s1 <- paste("", sep = "")
	for ( i in 1:lparts$R ) {
		s1 <- paste(s1, ldparms[i], sep = " ")
		k <- k + 1
	}
	cat(s1, file = fileConn, sep = "\n", append = T)
	cat("lambda:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$lambda
	s1 <- paste("", sep = "")
	for ( i in 1:lparts$R ) {
		s1 <- paste(s1, ldparms[i], sep = " ")
		k <- k + 1
	}
	cat(s1, file = fileConn, sep = "\n", append = T)
	cat("tau:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$tau
	s1 <- paste("", sep = "")
	for ( i in 1:lparts$R ) {
		s1 <- paste(s1, ldparms[i], sep = " ")
		k <- k + 1
	}
	cat(s1, file = fileConn, sep = "\n", append = T)
	cat("m:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$m
	s1 <- paste("", sep = "")
	for ( i in 1:lparts$m ) {
		s1 <- paste(s1, ldparms[i], sep = " ")
		k <- k + 1
	}
	cat(s1, file = fileConn, sep = "\n", append = T)
	cat("mm:", file = fileConn, sep = "\n", append = T)
	ldparms <- pars$dparms$mm
	s1 <- paste("", sep = "")
	for ( i in 1:lparts$mm ) {
		s1 <- paste(s1, ldparms[i], sep = " ")
		k <- k + 1
	}
	cat(s1, file = fileConn, sep = "\n", append = T)
	cat("$$", file = fileConn, sep = "\n", append = T)
	close(fileConn)
	return(paramfile)
}

score <- function (datafile, paramfile, keys) {
	out <- system(paste("gcdm_printscore -O hes -X ddmrna -H -i 1 -s dde21d -a 1e-2 -g n -Z -x input -p -A ari", datafile, paramfile, sep = " "), intern = T, ignore.stderr = T)
	res <- unlist(strsplit(out, "\\s+"))
	res = as.numeric(res[keys])
	return(res)
}

objfunc <- function(x, pars, datafile, keys, coeffs, buf = F) {
	paramfile <- tempfile(pattern = "rgcdm", tmpdir = tempdir(), fileext = "")
#	paramfile <- "hhh"
	pars$dparms <- vec2params(pars, x)
	write_params(pars, paramfile)
	res <- score(datafile, paramfile, keys)
#	cat(res[1], " ", res[2], " ", res[3], " ", res[4], " ", res[5], "\n")
	unlink(paramfile)
#	cat(paramfile)
	func <- sum(coeffs*res)
	if (buf == T) {
		return(list(func = func, res = coeffs*res))
	}
	return(func)
}

pobjfunc <- function(x) {
	paramfile <- tempfile(pattern = "rgcdm", tmpdir = tempdir(), fileext = "")
#	paramfile <- "hhh"
	pars$dparms <- vec2params(pars, x)
	write_params(pars, paramfile)

	out <- system(paste("gcdm_printscore -O hes -X ddmrna -H -i 1 -s dde21d -a 1e-2 -g n -Z -x input -p -A ari", datafile, paramfile, sep = " "), intern = T, ignore.stderr = T)
	res <- unlist(strsplit(out, "\\s+"))
	res = as.numeric(res[keys])

#	cat(res[1], " ", res[2], " ", res[3], " ", res[4], " ", res[5], "\n")
	unlink(paramfile)
#	cat(paramfile)

	func <- sum(coeffs*res)
	return(list(func = func, res = coeffs*res))
}

read_deep_log <- function(fn, nkeys) {
	content <- readLines(fn)
	content <- content[length(content)]
	token <- unlist(strsplit(gsub("\\s+", ",", gsub(":", ",", content)), ",", fixed = T))
	inds <- seq(2, 2 * (3 + nkeys), 2)
	token <- token[inds]
	wtime <- token[1]
	tau <- token[2]
	cost <- token[4]
	nf <- token[8]
	return(list(wtime = as.numeric(wtime), cost = as.numeric(cost), tau = as.numeric(tau), nf = as.numeric(nf)))
}

read_deep_log2 <- function(fn, nkeys, nparms) {
	content <- readLines(fn)
	content <- content[length(content)]
	token <- unlist(strsplit(gsub("\\s+", ",", gsub(":", ",", content)), ",", fixed = T))
	inds <- seq(2, 2 * (3 + nkeys + nparms), 2)
	token <- token[inds]
	wtime <- token[1]
	tau <- token[2]
	cost <- token[3]
	x <- token[(4 + nkeys):(3 + nkeys + nparms)]
	return(list(wtime = as.numeric(wtime), cost = as.numeric(cost), tau = as.numeric(tau), x = as.numeric(x)))
}

paramfile <- ""
datafile <- ""
logfile <- ""
template <- ""
nsteps <- 2
action <- ""

args <- commandArgs(trailingOnly = TRUE)
if ( length(args) > 0 ) {
	paramfile <- args[1]
} else {
	cat("Usage: script paramfile datafile logfile nsteps template action\n")
}
if ( length(args) > 1 ) {
	datafile <- args[2]
}
if ( length(args) > 2 ) {
	logfile <- args[3]
}
if ( length(args) > 3 ) {
	nsteps <- as.numeric(args[4])
}
if ( length(args) > 4 ) {
	template <- args[5]
}
if ( length(args) > 5 ) {
	action <- args[6]
}

# Full tweak
text_tweak="1;1;1;1;1;1;1;1; 1;1;1;1;0;0;0;0;1;1;1;1;0;0;0;0;1;1;1;1;0;0;0;0;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0; 1;1;1;1;1;1;1;1;0;0;0;0;1;1;1;1;1;1;1;1;0;0;0;0;1;1;1;1;1;1;1;1;0;0;0;0;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0; 1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1; 1;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0; 1;1;1;1;1;1;1;1; 1;1;1;1;0;0;0;0; 1; 1;"

text_lbound="0.01;0.01;0.01;0.01;1;1;1;1; -700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-700;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1; 0;0;0;-700;-70;-70;-70;0;-1;-1;-1;-1;0;0;-700;-700;-70;-70;-70;0;-1;-1;-1;-1;0;0;-700;-700;-70;-70;-70;0;-1;-1;-1;-1;0;0;-700;-700;-70;-70;-70;0;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1; -1;-1;-1;-1;-1;-1;-1;-1; 0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;1;1;1;1;1;1;1;1;1;1;1;1; 0.000001;-1;-1;-1;-1;-1;-1;-1; 0;0;0;0;0;0;0;0; 2;2;2;2;2;2;2;2; 2;2;2;2;-5;-5;-5;-5; 110; 1;"

text_hbound="1;1;1;1;15;15;15;15; 700;0;0;0;700;700;700;700;0;700;0;0;700;700;700;700;0;0;700;0;700;700;700;700;0;0;0;700;700;700;700;700; 1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1; 700;700;700;0;0;0;0;70;1;1;1;1;700;700;0;0;0;0;0;70;1;1;1;1;700;700;0;0;0;0;0;70;1;1;1;1;700;700;0;0;0;0;0;70;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;5;5;5;5;5;5;5;5;5;5;5;5; 0.00005;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 25;25;25;25;5;5;5;5; 8;8;8;8;0;0;0;0; 170; 4;"

text_mask="0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0; 0; 0;"

text_limited="1;1;1;1;1;1;1;1; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0; 0;0;0;0;0;0;0;0; 1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1; 1;"

text_scale="1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1;1;1;1;1;1;1;1; 1; 1;"

parts <- list(R=8, T=64, E=96, M=8, P=24, h=8, D=8, lambda=8, tau=8, m=1, mm=1)
tweak <- as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text_tweak), ";"))))
lbound <- as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text_lbound), ";"))))
hbound <- as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text_hbound), ";"))))
limited <- as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text_limited), ";"))))
mask <- as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text_mask), ";"))))

dparms <- read_params(paramfile)
pars <- list(parts=parts, dparms=dparms, tweak=tweak, lbound=lbound, hbound=hbound)
keys <- c(148, 5, 9, 14, 18, 23, 27, 32, 36, 41, 45, 50, 54, 59, 63, 68, 72, 77, 81, 95, 99, 104, 108, 152)
coeffs <- c(0, 0.0000001, 1, 0.0000001, 1, 0.0000001, 1, 0.0000001, 1, 0.0000001, 1, 0.0000001, 1, 0.0000001, 1, 0.0000001, 1, 0.00000001, 0.05, 0.00000001, 0.05, 0.00000001, 0.05, 0.0001)
eqparms <- params2vec(pars)
if (action == "params2vec") {
	cat('dparms =', paste0(eqparms, sep = ';'), '\n', sep = '')
	offset <- as.numeric(template)
	if (offset > 0) {
		hbound[eqparms > 0] <- eqparms[eqparms > 0] * (1 + offset)
		hbound[eqparms <= 0] <- eqparms[eqparms <= 0] * (1 - offset)
		lbound[eqparms > 0] <- eqparms[eqparms > 0] * (1 - offset)
		lbound[eqparms <= 0] <- eqparms[eqparms <= 0] * (1 + offset)
	}
	cat('lbound =', paste0(lbound, sep = ';'), '\n', sep = '')
	cat('hbound =', paste0(hbound, sep = ';'), '\n', sep = '')
	cat(length(eqparms), '\n')
	cat(length(hbound), '\n')
	cat(length(lbound), '\n')
	q()
}
eqparms <- eqparms[tweak == 1]
limited <- limited[tweak == 1]
lbound <- lbound[tweak == 1]
hbound <- hbound[tweak == 1]

res <- objfunc(eqparms, pars, datafile, keys, coeffs)
cat(res, '\n')

if (action == "testfiles") {
	lbound[limited == 0] <- -MAX_DOUBLE
	hbound[limited == 0] <- MAX_DOUBLE
	#cat(res[1], " ", res[2], " ", res[3], " ", res[4], " ", res[5], "\n")
	#sample_params(pars, mask, nsteps, datafile, testfile)
	#template <- "toy_a.*.ini.hopt_log_0$"
	testfiles <- list.files(pattern = paste0(template, ".log$"), all.files = FALSE,
	           full.names = FALSE, recursive = FALSE,
	           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
#	cat(length(testfiles), '\n')
	pfiles <- list.files(pattern = paste0(template, "-deep-output$"), all.files = FALSE,
	           full.names = FALSE, recursive = FALSE,
	           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
	nfiles <- length(testfiles)
#	cat(length(pfiles), '\n')
	fileConn <- file(logfile, open = "a")
	cat("N", "fbest", "distPar", "numVal", "cpuTime", seq(1:length(eqparms)), "\n", file = fileConn, append = T, sep = ',')
	close(fileConn)
	for (i in 1:nfiles) {
#		cat(testfiles[[i]], '\n')
		curr <- read_deep_log(testfiles[[i]], length(keys))
#		cat(pfiles[[i]], '\n')
		x <- read_params(pfiles[[i]])
		cpars <- list(parts=parts, dparms=x, tweak=tweak, lbound=lbound, hbound=hbound)
		curr$x <- params2vec(cpars)
		curr$x <- curr$x[tweak == 1]
	#	dist <- 100*sqrt(sum((curr$x - eqparms)^2)/sum(eqparms^2))
	#	dist <- 100*sum(abs(curr$x - eqparms)/abs(eqparms))/length(eqparms)
		dist <- sqrt(sum(w*(curr$x - eqparms)^2))/sqrt(sum(w*eqparms^2))
		fileConn <- file(logfile, open = "a")
		cat(testfiles[[i]], curr$cost, dist, curr$nf, curr$wtime, curr$x, '\n', file = fileConn, append = T, sep = ',')
		close(fileConn)
	}
}
if (action == "sample") {
	offset <- as.numeric(template)
	if (offset > 0) {
		hbound[eqparms > 0] <- eqparms[eqparms > 0] * (1 + offset)
		hbound[eqparms <= 0] <- eqparms[eqparms <= 0] * (1 - offset)
		lbound[eqparms > 0] <- eqparms[eqparms > 0] * (1 - offset)
		lbound[eqparms <= 0] <- eqparms[eqparms <= 0] * (1 + offset)
	}
	sample_params(logfile, res, eqparms, lbound, hbound, nsteps, pars, datafile, keys, coeffs)
	samdat <- read.csv(logfile)
	samdat <- samdat[!is.nan(samdat$distF),]
	gdat <- samdat[samdat$distF < 5,]
	q <- lapply(unique(gdat$ipar), FUN=function(a) { return(max(gdat[gdat$ipar == a, ]$distP)) })
	w <- as.numeric(as.vector(unlist(strsplit(sub("^\\s+", "", text_w), ";"))))
	w[tweak == 1] <- unlist(q)
	cat(w, '\n', sep = ';')
}
if (action == "deoptim") {
	library(DEoptim)
	controlDE <- list(reltol = 1e-8, steptol = 10, itermax = 10000, trace = 5, strategy = 5, p = 0.3, NP = 100, initialpop = NULL)
	fileConn <- file(logfile, open = "a")
	cat("N", "fbest", "distPar", "numVal", "cpuTime", seq(1:length(eqparms)), "\n", file = fileConn, append = T, sep = ',')
	close(fileConn)
	for (i in 1:nsteps) {
		set.seed(i)
		ptm <- proc.time()
		Results <- DEoptim(fn = objfunc, pars = pars, datafile = datafile, keys = keys, coeffs = coeffs, lower = lbound, upper = hbound, control = controlDE)
		ptm <- proc.time() - ptm
#		dist <- 100*sum((Results$optim$bestmem-eqparms)^2)/sum(eqparms^2)
		dist <- sqrt(sum(w*(Results$optim$bestmem - eqparms)^2))/sqrt(sum(w*eqparms^2))
		fileConn <- file(logfile, open = "a")
		cat(i, Results$optim$bestval, dist, Results$optim$nfeval, ptm[3], Results$optim$bestmem, '\n', file = fileConn, append = T, sep = ',')
		close(fileConn)
	}
}
if (action == "meigor") {
	library(MEIGOR)
#	lbound[limited == 0] <- -MAX_DOUBLE
#	hbound[limited == 0] <- MAX_DOUBLE
	# An initial point is specified.
	# The number of integer variables is specified (mandatory).
	# No local solver is available for mixed-integer problems at the moment.
	# Stop criterion determined by the CPU time (2 seconds).
	#========================= PROBLEM SPECIFICATIONS ===========================
	fileConn <- file(logfile, open = "a")
	cat("N", "fbest", "numVal", "cpuTime", seq(1:length(eqparms)), "\n", file = fileConn, append = T, sep = ',')
	close(fileConn)

	if (template == "parallel") {
		problem <- list (f = "pobjfunc", x_L = lbound, x_U = hbound, int_var = 2)
	#Set 1 nodes and 2 cpu ' s per node
		n_nodes <- 1
		n_cpus_per_node <- 4
	#Set different values for dim_refset, bal and n2
	#for each of the 10 cpu's to be used
		counter = 0
		opts = list()
		hosts = c()
		for (i in 1:n_nodes) {
			for (j in 1:n_cpus_per_node) {
				counter = counter + 1
	#Set the name of every thread
				hosts = c(hosts, 'localhost')
				opts[[counter]] = list()
	#Set common options for each thread
#				opts[[counter]]$maxeval = 10000;
				opts[[counter]]$tolc = 10^-8
#				opts[[counter]]$ndiverse = 5
				opts[[counter]]$maxtime = 500
	#Options not set will take default values for every thread
			}
		}
	# Set the address of each machine, defined inside the ' for ' loop
		opts$hosts = hosts
	# Do not define the additional options for cooperative methods (e.g., ce_maxtime, ce_isparallel, et
	# They will take their default values
#		opts$ce_niter = 500
#		opts$ce_maxeval = 1200000
		opts$ce_maxtime = 120000
		opts$ce_type = "SOCKS"
		opts$ce_isparallel = TRUE
#		opts$global_save_list = list("objfunc", "vec2params", "write_params", "score", "pars", "datafile", "keys", "coeffs")
		opts$global_save_list = list("pobjfunc", "vec2params", "write_params", "score", "pars", "datafile", "keys", "coeffs")
	# Call the solver
#		Results <- MEIGO(problem, opts, algorithm = "CeSSR", pars = pars, datafile = datafile, keys = keys, coeffs = coeffs, buf = T)
#		Results <- MEIGO(problem, opts, algorithm = "CeSSR", pars, datafile, keys, coeffs)
		Results <- MEIGO(problem, opts, algorithm = "CeSSR")
	} else {
		problem <- list (f = "objfunc", x_L = lbound, x_U = hbound, int_var = 2)
		for (i in 1:nsteps) {
			opts <- list (maxtime = 120000, tolc = 10^-8)
	# Call optimizer
			Results <- MEIGO (problem, opts, algorithm = "ESS", pars = pars, datafile = datafile, keys = keys, coeffs = coeffs, buf = T)
	#
	#		dist <- 100*sum((Results$xbest-eqparms)^2)/sum(eqparms^2)
	#		dist <- sqrt(sum(w*(Results$optim$xbest - eqparms)^2))/sqrt(sum(w*eqparms^2))
	# Write log
		}
	}
	fileConn <- file(logfile, open = "a")
#		cat(i, Results$fbest, dist, Results$numeval, Results$cpu_time, Results$xbest, '\n', file = fileConn, append = T, sep = ',')
	cat(i, Results$fbest, Results$numeval, Results$cpu_time, Results$xbest, '\n', file = fileConn, append = T, sep = ',')
	close(fileConn)
}
