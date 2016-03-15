## LSKM.R
##
##--------------------------------------------------------------------------
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##--------------------------------------------------------------------------

##--------------------------------------------------------------------------
##
##  R code to implement LSKM method of Kwee et al., AJHG 82:386-397, 2008.
##
## version 1.0 4/15/2008
##--------------------------------------------------------------------------


##--------------------------------------------------------------------------
## define the names of the input trait, genotype, and covariate files
## data format details are in the README.txt file
##--------------------------------------------------------------------------
trait_file     <- "trait.dat"
genotype_file  <- "genotype.dat"
covariate_file <- "covariate.dat"


##--------------------------------------------------------------------------
## print information about this script:
##--------------------------------------------------------------------------
version <- "1.0"
output_file <- "LSKM.out"
sink(output_file)

cat("------------------------------------------------------\n")
cat("LSKM.R, version",version,"\n")
cat("Reference: Kwee et al., AJHG 82:386-397, 2008 \n")
analysis_date <- paste("Date & time of analysis: ", date(), sep="")
cat(analysis_date, "\n")
cat("------------------------------------------------------\n")

err = 0  # parameter to keep track of error in input files

##--------------------------------------------------------------------------
## populate genotype, covariate, and trait data arrays from input files:
##--------------------------------------------------------------------------
Y <- data.matrix(read.table(trait_file, header=TRUE))
n <- dim(Y)[1]                 # number of subjects in dataset

gen <- data.matrix(read.table(genotype_file,header=TRUE))
typedloci <- dim(gen)[2]       # number of typed SNPs
if (dim(gen)[1] != n) {
    cat("Error:  different # of individuals in genotype.txt and trait.txt \n")
    err = 1
}


X <- matrix(1,n)
if (file.access(covariate_file) == 0) {
  	Xe <- data.matrix(read.table(covariate_file, header=TRUE))
    n_cov <- dim(Xe)[2]         # number of measured environmental covariates
    if (dim(Xe)[1] != n) {
        cat("Error:  different # of individuals in covariate.txt and trait.txt \n")
        err = 1
    }
    else {
        X=cbind(X,Xe)
    }
} else {
    n_cov=0  	
}
colnames(X)[1] = "Intercept"


if (err == 0) {
    ## calculate XtX_inverse:
    XtX_inv <- solve (t(X) %*% X)

    ## calculate parameter estimates for intercept and environmental covariates:
    betahat <- XtX_inv %*% t(X) %*% Y
    colnames(betahat) <- "Estimate"
    sigsqhat <- (t(Y - X %*% betahat) %*% (Y - X %*% betahat)) / (n - n_cov - 1) 
    beta_var = sigsqhat[1,1] * XtX_inv


    ## calculate K matrix
    K <- array(0, c(n,n))

    ## MAF weights:
	wt <- function(loc) {
  		sum <- sum(gen[,loc])/(2*n)
    	if (sum == 0) sum = 0.01
    	return(1/sqrt(sum))
    }
	weight <- sapply(1:typedloci, wt)
  	sum_weights <- sum(weight)

    ## measure of IBS sharing:
    ibs <- function(sub1, sub2, loc) {
        diff = abs(gen[sub1, loc] - gen[sub2, loc])
        ibs1 = 2 - diff
        return(ibs1)
	}

    ## calculate the MAF-weighted freq of the rare allele at each locus
    for (j in 1:n) {
        for (k in j:n) {
            sum = 0
            for (ell in 1:typedloci) {
                sum = sum + ibs(j, k, ell) * weight[ell]
            } # for ell
            K[j,k] = sum/(sum_weights)
        } # for k
    } # for j

    ## populate remaining components of the matrix by symmetrizing:
    for (j in 2:n) {
        for (k in 1:(j - 1)) {
            K[j,k] = K[k,j]
        }
    }

    ## create identity matrix I
    I <- diag(n)

    ## calculate necessary quantities for the approximate distribution:
	proj <- I - X %*% XtX_inv %*% t(X)
	
	pk   <- proj %*% K
	pkp  <- pk %*% proj
	psq  <- proj %*% proj
	pksq <- pk %*% pk

	i_tt      <- sum(diag(pksq))/2
	e         <- sum(diag(pk))/2
	i_ts      <- sum(diag(pkp))/2
	i_ss      <- sum(diag(psq))/2
    i_tt_reml <- i_tt - (i_ts**2)/i_ss    
	kappa     <- i_tt_reml/(2*e)
	nu        <- 2*(e**2)/i_tt_reml

    ## calculate score statistic:
	XB    <- X %*% betahat
	delta <- Y - XB
	score <- (t(delta) %*% K %*% delta) / (2*sigsqhat*kappa)
	
    ## calculate p-value for score statistic:
	p <- pchisq(score, nu, lower.tail=FALSE)

    ## print results:
    cat("\n")
    cat("Number of subjects:",n,"\n")
    cat("Number of SNPs:",typedloci,"\n")
    cat("Number of environmental covariates:",n_cov,"\n")
    cat("------------------------------------------------------\n \n")

    cat("Association test between SNPs and",colnames(Y),"\n")
    cat("Score statistic:",score,"\n")
    cat("Approximate degrees of freedom:",nu," \n")
    cat("P-value:", p,"\n")
    cat("------------------------------------------------------\n \n")

    cat("\t        Analysis of Fixed Effects \n")
    cat("\t ----------------------------------------- \n")
    cat("\t \t \t Parameter \t  Standard \n")
    cat("\t Variable \t Estimate \t  Error \n")
    cat("------------------------------------------------------ \n")
    for (j in 0:(n_cov + 1)) {
        if (j == 1) {
            cat("\t",rownames(betahat)[j]," \t",betahat[j]," \t",sqrt(diag(beta_var)[j]),"\n")
     	}
        if (j > 1) {
            cat("\t",rownames(betahat)[j]," \t \t",betahat[j]," \t",sqrt(diag(beta_var)[j]),"\n")
        }
    }
    cat("------------------------------------------------------ \n \n")
    cat("Subject-specific variance estimate:", sigsqhat, "\n\n")
}
sink()
