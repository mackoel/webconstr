#!/usr/bin/Rscript

require("optparse")

# Setup
genenames <- c("hb", "Kr", "gt", "kni")
tfnames <- c("hb", "Kr", "gt", "kni", "constr", "bcd", "cad", "tll", "hkb")
fulldbfile <- "sites_full_db.csv"
parmfilename <- "_con.024k.02.ini-deep-output"
constrgdfile <- "dual_mRNA58_constr.gd"
constrdbfile <- "constr_db.csv"
term <- 'png'
# -1 - don't check
#  1 - accept if TRUE
#  0 - accept if FALSE
check_constr <- -1
check_dnase <- -1
check_codereg <- 0
pwm_low <- -1
pwm_high <- -1
# 1 - logical AND
# 0 - logical OR
and_constr <- 1
and_dnase <- 1
and_pwm <- 1
# Construct description:
# target, global start, global end, local start, local end
constr_target <- "NA"
constr_global_start <- -1
constr_global_end <- -1
constr_local_start <- -1
constr_local_end <- -1
constr_name <- "NA"
outputdir <- "."
outputann <- "annotated_sites.dat"
constrparmfile <- "constr.par"
constrgd <- "constr.gd"

runmodel <- function(constrgd, constrparmfile, GG, cname, plotname, term = "png") {
	# Colors
	Ucolor="#000000"
	Hcolor="#ff0000"
	Kcolor="#00ff00"
	Gcolor="#0000ff"
	Ncolor="#ff00ff"
	Ccolor="#9a9a64"
	Tcolor="#bc1c8d"
	Bcolor="#b693d0"
	Jcolor="#00ff00"
	HKcolor="#007f00"
	HGcolor="#00007f"
	HBcolor="#7f0000"
	KKcolor="#ffa700"
	KNcolor="#edbbff"
	# Options
	options="-O hes -x input -X ddual -H -i 1.0 -s dde21d -a 1e-2 -g n -Z -j times -A ari "
	rdatfile="gt2.dat"
	pdatfile=paste0(constrgd, ".uof")
# Run model
	command <- paste0("gcdm_printscore ", options, constrgd, " ", constrparmfile)
#	cat(command, "\n")
	system(command, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE)
# make plot
	gplt=paste0(cname, ".gplt")
	gff="/usr/share/fonts/dejavu/DejaVuSans.ttf"
# Gnuplot script template
	textgplt <- "
	set macros
	# Placement of the a,b,c,d labels in the graphs
	POS = \"at graph 0.1,0.9 font ',8'\"

	# x- and ytics for each row resp. column
	#NOXTICS = \"set xtics ('' 8192, '' 8211, '' 8227, '' 8240, '' 8249, '' 8258, '' 8267, '' 8271, '' 8275, '' 8284, '' 8291); unset xlabel\"
	#XTICS = \"set xtics ('0' 8192, '20' 8211, '35' 8227, '56' 8249, '74' 8267, '92' 8284); unset xlabel\"

	NOXTICS = \"set xtics ('' 8192, '' 8212, '' 8232, '' 8252, '' 8272); unset xlabel\"
	XTICS = \"set xtics ('0' 8192, '20' 8212, '40' 8232, '60' 8252, '80' 8272); unset xlabel\"

	NOYTICS = \"set ytics ('' 0, '' 100, '' 200); unset ylabel\"
	YTICS = \"set ytics ('0' 0, '100' 100, '200' 200); unset ylabel\"

	set lmargin at screen 0.20
	set rmargin at screen 0.95

	TOP = 0.95
	DY = 0.85
	RIGHT = 0.95
	DX = 0.3

	set yrange [0:400]
	set xrange [8227:8284]
	set xtics nomirror
	set ytics nomirror

	set xlabel font \"@gff@,14\"
	set ylabel font \"@gff@,14\"
	set key font \"@gff@,10\"
	set xtics font \"@gff@,12\"
	set ytics font \"@gff@,12\"
	set title font \"@gff@,10\"

	set term @term@ size 1200,400
	set output \"@me@\"

	set multiplot

	## prot

	# --- GRAPH a
	set label 1 'T7' @POS
	@XTICS; @YTICS

	set tmargin at screen TOP-0*DY
	set bmargin at screen TOP-1*DY

	set lmargin at screen RIGHT-1*DX
	set rmargin at screen RIGHT-0*DX

	set y2tics
	set autoscale y2

	plot '@pdatfile@' using 1:7 index 7 axes x1y2 with lines lc rgbcolor \"@colorGene@\" lw 2 title '@constr@ model prot','@rdatfile@' index 7 using 1:@colGene@ with points lc rgbcolor \"@colorGene@\" pt 1  ps 0.5 title '@strGene@ data'

	# --- GRAPH d
	unset label
	@XTICS; @NOYTICS

	set tmargin at screen TOP-0*DY
	set bmargin at screen TOP-1*DY

	set lmargin at screen RIGHT-3*DX
	set rmargin at screen RIGHT-2*DX

	set key
	set ylabel 'Protein Concentration'
	plot '@pdatfile@' using 1:12 index 7 axes x1y2 with lines lc rgbcolor \"@colorGene@\" lw 2 title '@constr@ model mRNA', '@rdatfile@' index 7 using 1:@colMrnk@ with points lc rgbcolor \"@colorGene@\" pt 1  ps 0.5 title 'mrnk @strGene@ data'

	unset multiplot
	"
# Set up variables
	if (GG == "hb") {
		colGene=3
		colorGene=Hcolor
		strGene="Hb"
		colMrnk=7
		colMrnkM=8
	}
	if (GG == "Kr") {
		colGene=4
		colorGene=Kcolor
		strGene="Kr"
		colMrnk=8
		colMrnkM=9
	}
	if (GG == "gt") {
		colGene=5
		colorGene=Gcolor
		strGene="Gt"
		colMrnk=9
		colMrnkM=10
	}
	if (GG == "kni") {
		colGene=6
		colorGene=Ncolor
		strGene="Kni"
		colMrnk=10
		colMrnkM=11
	}
# Substitute variables
	textgplt <- gsub("@colorGene@", colorGene, textgplt)
	textgplt <- gsub("@colGene@", colGene, textgplt)
	textgplt <- gsub("@colMrnk@", colMrnk, textgplt)
	textgplt <- gsub("@strGene@", strGene, textgplt)
	textgplt <- gsub("@me@", plotname, textgplt)
	textgplt <- gsub("@term@", term, textgplt)
	textgplt <- gsub("@constr@", cname, textgplt)
	textgplt <- gsub("@gff@", gff, textgplt)
	textgplt <- gsub("@rdatfile@", rdatfile, textgplt)
	textgplt <- gsub("@pdatfile@", pdatfile, textgplt)
	cat(textgplt, file = gplt, sep = '\n')
# Run gnuplot
	command <- paste0("gnuplot ", gplt)
#	cat(command, "\n")
	system(command, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE)
}

mkgd <- function(afile, gfile, target, out) {
	text <- readLines(gfile)
	text <- gsub("@ANNFILE@", afile, text)
	text <- gsub("@target@", target, text)
	cat(text, file = out, sep = '\n')
}

mkann <- function(tab, dbfn = 'NA') {
#	dim(tab[tab$target == "hb",])
#	dim(tab[tab$target == "Kr",])
#	dim(tab[tab$target == "gt",])
#	dim(tab[tab$target == "kni",])
#	summary(tab)
#	write.table(cbind(tab$adjstart, tab$loc..strand, tab$tf), file=dbfn, row.names=F, quote=F, sep=" ")
#	unlink(outputann)
	textann <- capture.output({
		for ( i in seq(1, 4)) {
			cat(file = "", append=TRUE, sprintf(">%s %d\n", genenames[i], dim(tab[tab$target == genenames[i],])[1] ))
			write.table(file = "", append=TRUE, cbind(tab[tab$target == genenames[i],]$adjstart, as.character(tab[tab$target == genenames[i],]$loc..strand), tab[tab$target == genenames[i],]$tfidx, -log(tab[tab$target == genenames[i],]$energy), tab[tab$target == genenames[i],]$length, tab[tab$target == genenames[i],]$energy), quote=F, row.names=F, col.names=F)
		}
	})
	return(textann)
}

mkparfile <- function(parmfilename, ctarget) {
	ngenes <- 8
	nmrna <- 4
	negenes <- 4
	mgenes <- 1
	ncouples <- 2
	II <- which(genenames == ctarget)
	text <- readLines(parmfilename, encoding="UTF-8")
	i <- 2
	i <- i + 1
		toks <- unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+"))
		R <- as.vector(toks)
		i <- i + 1
		i <- i + 1
		T <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		for ( k in 2:ngenes ) {
			i <- i + 1
			t <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
			T <- rbind(T, t)
		}
		i <- i + 1
		i <- i + 1
		E <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		for ( k in 2:ngenes ) {
			i <- i + 1
			e <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
			E <- rbind(E, e)
		}
		i <- i + 1
		i <- i + 1
		M <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		for ( k in 2:ngenes ) {
			i <- i + 1
			m <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
			M <- rbind(M, m)
		}
		i <- i + 1
		i <- i + 1
		P <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		for ( k in 2:ngenes ) {
			i <- i + 1
			p <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
			P <- rbind(P, p)
		}
		i <- i + 1
		i <- i + 1
		h <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		i <- i + 1
		i <- i + 1
		D <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		i <- i + 1
		i <- i + 1
		lambda <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		i <- i + 1
		i <- i + 1
		tau <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		i <- i + 1
		i <- i + 1
		ranget <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
		i <- i + 1
		i <- i + 1
		mm <- as.vector(unlist(strsplit(sub("^\\s+", "", text[i]), "\\s+")))
	textpar <- capture.output({
		cat("$input\n")
		cat("R:\n")
		for ( k in 1:nmrna ) {
			cat(" ", R[k])
		}
			cat(" ", R[II])
		for ( k in (nmrna + 1):ngenes ) {
			cat(" ", R[k])
		}
			cat(" ", R[II + nmrna])
		cat("\n")
		cat("T:\n")
		for ( i in 1:nmrna ) {
			for ( k in 1:nmrna ) {
				cat(" ", T[i,k])
			}
				cat(" ", 0.00000000)
			for ( k in (nmrna + 1):ngenes ) {
				cat(" ", T[i,k])
			}
			cat(" ", 0.00000000)
			cat("\n")
		}
			for ( k in 1:nmrna ) {
				cat(" ", T[II,k])
			}
				cat(" ", 0.00000000)
			for ( k in (nmrna + 1):ngenes ) {
				cat(" ", T[II,k])
			}
			cat(" ", 0.00000000)
			cat("\n")
		for ( i in (nmrna + 1):ngenes ) {
			for ( k in 1:nmrna ) {
				cat(" ", T[i,k])
			}
				cat(" ", 0.00000000)
			for ( k in (nmrna + 1):ngenes ) {
				cat(" ", T[i,k])
			}
			cat(" ", 0.00000000)
			cat("\n")
		}
			for ( k in 1:nmrna ) {
				cat(" ", T[II + nmrna,k])
			}
				cat(" ", 0.00000000)
			for ( k in (nmrna + 1):ngenes ) {
				cat(" ", T[II + nmrna,k])
			}
			cat(" ", 0.00000000)
			cat("\n")

		cat("E:\n")
		for ( i in 1:nmrna ) {
			for ( k in 1:(negenes + nmrna) ) {
				cat(" ", E[i,k])
			}
			cat(" ", 0.00000000)
			cat("\n")
		}
			for ( k in 1:(negenes + nmrna) ) {
				cat(" ", E[II,k])
			}
			cat(" ", 0.00000000)
			cat("\n")
		for ( i in (nmrna + 1):ngenes ) {
			for ( k in 1:(negenes + nmrna) ) {
				cat(" ", E[i,k])
			}
			cat(" ", 0.00000000)
			cat("\n")
		}
			for ( k in 1:(negenes + nmrna) ) {
				cat(" ", E[II + nmrna,k])
			}
			cat(" ", 0.00000000)
			cat("\n")

		cat("M:\n")
		for ( i in 1:nmrna ) {
			for ( k in 1:mgenes ) {
				cat(" ", M[i,k])
			}
			cat("\n")
		}
			for ( k in 1:mgenes ) {
				cat(" ", M[II,k])
			}
			cat("\n")
		for ( i in (nmrna + 1):ngenes ) {
			for ( k in 1:mgenes ) {
				cat(" ", M[i,k])
			}
			cat("\n")
		}
			for ( k in 1:mgenes ) {
				cat(" ", M[II + nmrna,k])
			}
			cat("\n")

		PP <- c()
		for ( i in 1:ngenes ) {
			for ( k in 1:ncouples ) {
				PP <- c(PP, P[i,k])
			}
		}

		PP <- c(PP[1:nmrna], PP[II], PP[(nmrna + 1):(nmrna + negenes + nmrna)], PP[(II + negenes + nmrna)], PP[(nmrna + negenes + nmrna + 1):(nmrna + negenes + nmrna + negenes)], 0, 0)

		cat("P:\n")
		j <- 1
		for ( i in 1:(ngenes + 2) ) {
			for ( k in 1:ncouples ) {
				cat(" ", PP[j])
				j <- j + 1
			}
			cat("\n")
		}

		cat("h:\n")
		for ( k in 1:nmrna ) {
			cat(" ", h[k])
		}
			cat(" ", h[II])
		for ( k in (nmrna + 1):ngenes ) {
			cat(" ", h[k])
		}
			cat(" ", h[II + nmrna])
		cat("\n")

		cat("D:\n")
		for ( k in 1:nmrna ) {
			cat(" ", D[k])
		}
			cat(" ", D[II])
		for ( k in (nmrna + 1):ngenes ) {
			cat(" ", D[k])
		}
			cat(" ", D[II + nmrna])
		cat("\n")

		cat("lambda:\n")
		for ( k in 1:nmrna ) {
			cat(" ", lambda[k])
		}
			cat(" ", lambda[II])
		for ( k in (nmrna + 1):ngenes ) {
			cat(" ", lambda[k])
		}
			cat(" ", lambda[II + nmrna])
		cat("\n")

		cat("tau:\n")
		for ( k in 1:nmrna ) {
			cat(" ", tau[k])
		}
			cat(" ", tau[II])
		for ( k in (nmrna + 1):ngenes ) {
			cat(" ", tau[k])
		}
			cat(" ", tau[II + nmrna])
		cat("\n")

		cat("ranget:\n")
		cat(" ", ranget[1])
		cat("\n")

		cat("mm:\n")
		cat(" ", mm[1])
		cat("\n")
		cat("$$\n")
	})
	return(textpar)
}

getsitesdb <- function(fulldbfile) {
	# Read sites db
	ds <- read.csv(fulldbfile)
	# Correct names
	ds$target <- gsub("kr", "Kr", ds$target)
	ds$tf <- gsub("kr", "Kr", ds$tf)
	# Get maximal score
	LLRMAX <- lapply(tfnames, FUN = function(u) {
			if (dim(ds[ds$tf == u,])[1] > 0) {
				result = max(ds[ds$tf==u,]$PWM.score)
			} else {
				result = 1
			}
			return(result)
		})
	LLRMAX <- as.numeric(unlist(LLRMAX))
	# cat(LLRMAX, '\n')
	# Compute normalized energy
	energy <- apply(ds, 1, FUN=function(tab) {
					return(
						as.numeric(tab[ which(colnames(ds) == 'PWM.score') ]) /
						LLRMAX[ which( tfnames == tab[ which(colnames(ds) == 'tf') ] ) ]
					)
				}
			)
	# Compute tf index instead of name
	tfidx <- apply(ds, 1, FUN=function(tab) {
					return(
						which( tfnames == tab[ which(colnames(ds) == 'tf') ] ) - 1
					)
				}
			)
	# get motif length
	lw <- apply(ds, 1, FUN=function(tab) {
					return(
						nchar(tab[ which(colnames(ds) == 'loc..word') ] )
					)
				}
			)
	# add to table
	ds <- cbind(ds, energy, tfidx, length = lw)
	return(ds)
}

applyfilters <- function(db) {
	# Limit to DNAse accebility region?
	if (check_dnase == 1) {
		ddb <- db[db$overlaps.DNAse.acc..regions == "yes",]
	} else {
		ddb <- db
	}
	# Check overlap with constructs?
	if (check_constr == 1) {
		if (and_constr == 1) {
			ddb <- ddb[ddb$overlap_arina == "Yes",]
		}
		if (and_constr == 0) {
			ddb <- db[db$overlap_arina == "Yes" | db$overlaps.DNAse.acc..regions == "yes",]
		}
	}
	# Exclude coding region?
	if (check_codereg == 0) {
		ddb <- ddb[ddb$codereg != T,]
	}
	ddb <- ddb[ddb$PWM.score > 0,]
	return(ddb)
}

constrtext <- function(cre, modeldata) {
	tagsites <- modeldata[modeldata$target == cre[6],]
	tagsites <- tagsites[(as.numeric(as.character(cre[4])) <= tagsites$adjstart & tagsites$adjstart <= as.numeric(as.character(cre[5]))), ]
	tagsites <- rbind(tagsites, tagsites[(as.numeric(as.character(cre[4])) <= tagsites$adjend & tagsites$adjend <= as.numeric(as.character(cre[5]))), ])
	tagsites <- unique(tagsites)
#	if (dim(tagsites)[1] > 0) {
		textconstr <- capture.output({
			cat(">", cre[7], " ", dim(tagsites)[1], "\n", sep = "");
			if (dim(tagsites)[1] > 0) {
			write.table(file = "", cbind(tagsites$adjstart, as.character(tagsites$loc..strand), tagsites$tfidx, -log(tagsites$energy), tagsites$length, tagsites$energy),  row.names = F, col.names = F, sep = " ", quote = F, append = TRUE)
		}
		})
#	}
	return(textconstr)
}

# Main
# Parse options
option_list <- list(
	make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="Print extra output [default=false]"),
	make_option(c("-a", "--and_constr"), type = "integer", default = 1,
	help = "and_constr [default %default]", metavar = "number"),
	make_option(c("-A", "--and_pwm"), type = "integer", default = 1,
	help = "and_pwm [default %default]", metavar = "number"),
	make_option(c("-b", "--and_dnase"), type = "integer", default = 1,
	help = "and_dnase [default %default]", metavar = "number"),
	make_option(c("-c", "--check_constr"), type = "integer", default = -1,
	help = "check_constr [default %default]", metavar = "number"),
	make_option(c("-d", "--check_dnase"), type = "integer", default = -1,
	help = "check_dnase [default %default]", metavar = "number"),
	make_option(c("-e", "--check_codereg"), type = "integer", default = 0,
	help = "check_codereg [default %default]", metavar = "number"),
	make_option(c("-s", "--constr_global_start"), type = "integer", default = -1,
	help = "constr_global_start [default %default]", metavar = "number"),
	make_option(c("-f", "--constr_global_end"), type = "integer", default = -1,
	help = "constr_global_end [default %default]", metavar = "number"),
	make_option(c("-S", "--constr_local_start"), type = "integer", default = -1,
	help = "constr_local_start [default %default]", metavar = "number"),
	make_option(c("-F", "--constr_local_end"), type = "integer", default = -1,
	help = "constr_local_end [default %default]", metavar = "number"),
	make_option(c("-w", "--pwm_low"), type = "double", default = -1,
	help = "pwm_low [default %default]"),
	make_option(c("-g", "--pwm_high"), type = "double", default = -1,
	help = "pwm_high [default %default]"),
	make_option(c("-t", "--constr_target"), type="character", default="NA",
	help="constr_target [default %default]", metavar="string"),
	make_option(c("", "--fulldbfile"), type="character", default=fulldbfile,
	help="fulldbfile [default %default]", metavar="string"),
	make_option(c("", "--constrdbfile"), type="character", default=constrdbfile,
	help="constrdbfile [default %default]", metavar="string"),
	make_option(c("", "--parmfilename"), type="character", default=parmfilename,
	help="parmfilename [default %default]", metavar="string"),
	make_option(c("", "--constrgdfile"), type="character", default=constrgdfile,
	help="constrgdfile [default %default]", metavar="string"),
	make_option(c("-o", "--outputdir"), type="character", default="NA",
	help="outputdir [default %default]", metavar="string"),
	make_option(c("-n", "--constr_name"), type="character", default="NA",
	help="constr_target [default %default]", metavar="string")
)
opt_parser <- OptionParser(usage = "usage: %prog [options] outputann constrparmfile constrgd", option_list = option_list, description = "checkconstr", epilogue = "Send feedback to mackoel@gmail.com")
flaggs <- parse_args(opt_parser, args = commandArgs(trailingOnly = TRUE), print_help_and_exit = TRUE, positional_arguments = TRUE)
opts <- flaggs$options
args <- flaggs$args
if (opts$outputdir == "NA") {
	if (length(args) > 0) {
		outputann <- args[1]
	}
	if (length(args) > 1) {
		constrparmfile <- args[2]
	}
	if (length(args) > 2) {
		constrgd <- args[3]
	}
} else {
	outputdir <- opts$outputdir
	if(!file.exists(outputdir)) dir.create(outputdir, rec=TRUE)
	outputann <- paste0(outputdir, '/', outputann)
	constrparmfile <- paste0(outputdir, '/', constrparmfile)
	constrgd <- paste0(outputdir, '/', constrgd)
}
check_constr <- opts$check_constr
check_dnase <- opts$check_dnase
check_codereg <- opts$check_codereg
pwm_low <- opts$pwm_low
pwm_high <- opts$pwm_high
and_constr <- opts$and_constr
and_dnase <- opts$and_dnase
and_pwm <- opts$and_pwm
constr_global_start <- opts$constr_global_start
constr_global_end <- opts$constr_global_end
constr_local_start <- opts$constr_local_start
constr_local_end <- opts$constr_local_end
constr_name <- opts$constr_name
constr_target <- opts$constr_target
fulldbfile <- opts$fulldbfile
parmfilename <- opts$parmfilename
constrgdfile <- opts$constrgdfile
constrdbfile <- opts$constrdbfile
# Read sites db
ds <- getsitesdb(fulldbfile)
#summary(ds)
#dim(ds[ds$target == "hb",])
#dim(ds[ds$target == "Kr",])
#dim(ds[ds$target == "gt",])
#dim(ds[ds$target == "kni",])
#summary(ds)
dds <- applyfilters(ds)
outann <- mkann(dds)
cat(file = outputann, outann, sep = '\n')
if (constr_name != "NA" && constr_target == "NA") {
	constrdb <- read.csv(constrdbfile, header = T, stringsAsFactors = F)
	constr <- constrdb[constrdb$name == constr_name,]
	consi <- constrtext(unlist(constr), dds)
	cat(file = constrparmfile, sep = '\n', mkparfile(parmfilename, constr$target))
	cat(file = outputann, consi, append = T, sep = '\n')
	mkgd(outputann, constrgdfile, constr$target, constrgd)
	plotname = paste0(outputdir, '/', constr_name, ".", term)
	runmodel(constrgd, constrparmfile, constr$target, constr_name, plotname)
} else if (constr_target != "NA") {
	if (constr_target == "hb") {
		if (constr_local_start == -1) {
			constr_local_start <- 4523540 - constr_global_end + 12001
		}
		if (constr_local_end == -1) {
			constr_local_end <- 4523540 - constr_global_start + 12001
		}
	}
	if (constr_target == "Kr") {
		if (constr_local_start == -1) {
			constr_local_start <- constr_global_start - 21102137 + 1
		}
		if (constr_local_end == -1) {
			constr_local_end <- constr_global_end - 21102137 + 1
		}
	}
	if (constr_target == "gt") {
		if (constr_local_start == -1) {
			constr_local_start <- 2322996 - constr_global_end + 12001
		}
		if (constr_local_end == -1) {
			constr_local_end <- 2322996 - constr_global_start + 12001
		}
	}
	if (constr_target == "kni") {
		if (constr_local_start == -1) {
			constr_local_start <- 20688463 - constr_global_end + 12001
		}
		if (constr_local_end == -1) {
			constr_local_end <- 20688463 - constr_global_start + 12001
		}
	}
	constr_name <- paste0(constr_target, "_", constr_local_start, "_", constr_local_end)
	constr <- list(index=1, start=constr_global_start, end=constr_global_end, adjstart=constr_local_start, adjend=constr_local_end, target=as.character(constr_target), name=as.character(constr_name), numsites=0)
	consi <- constrtext(unlist(constr), dds)
	cat(file = constrparmfile, sep = '\n', mkparfile(parmfilename, constr$target))
	cat(file = outputann, consi, append = T, sep = '\n')
	mkgd(outputann, constrgdfile, constr_target, constrgd)
	plotname = paste0(outputdir, '/', constr_name, ".", term)
	runmodel(constrgd, constrparmfile, constr_target, constr_name, plotname)
} else {

}
cat(plotname, '\n')
