
#' A function for ...
#'
#' This function allows you to ...
#' @param query.p ...
#' @param reference.p ...
#' @param query.gL ...
#' @param reference.gL ...
#' @param na.rm ...
#' @param estimateNullDistr ...
#' @param nullDistrQuantiles ...
#' @param nullDistrN ...
#' @param tolerateWarnings ...
#' @param pAdjust.method.query ...
#' @param pAdjust.method.reference ...
#' @param lambda ...
#' @param rescale ...
#' @param plotRescaled ...
#' @param multicoreN ...
#' @param authors M. Fazel-Najafabadi
#' @keywords test
#' @export 
#' @examples
#' res <- batchTregGRS(....)

batchTregGRS <- function (query.p, reference.p, query.gL = NULL, reference.gL = NULL, 
    na.rm = TRUE, estimateNullDistr = TRUE, nullDistrQuantiles = c(0.9, 
        0.95, 0.99), nullDistrN = 100, tolerateWarnings = TRUE, 
    pAdjust.method.query = NULL, pAdjust.method.reference = NULL, 
    lambda = 0.005,rescale=T,plotRescaled=F,multicoreN=3) 
{
    postP <- function(p, scale = 1) {
        B <- -exp(1) * p * log(p)
        B[p > exp(-1)] <- 1
        1 - 1/(1 + 1/(max(1, scale) * B))
    }
    estimateH0 <- function(p, lambda) {
        pi0 <- 1
        tryCatch(pi0 <- qvalue(p, lambda = lambda)$pi0, error = function(e) {
            pi0 <- 1
# print("#2")
        })
# print("#2")
        pi0 <- min(pi0, 0.999)
# print("#2")
        pi0/(1 - pi0)
# print("#2")
    }
    RSCM <- function(p1, p2, s1, s2) {
        E <- 1/2 * (sum(p1 * s2)/sum(p1) + sum(p2 * s1)/sum(p2))
        mu.p1 <- mean(p1)
        mu.p2 <- mean(p2)
        mu.s1 <- mean(s1)
        mu.s2 <- mean(s2)
        sigma.p1 <- sd(p1)
        sigma.p2 <- sd(p2)
        sigma.s1 <- sd(s1)
        sigma.s2 <- sd(s2)
        gamma1 <- sum((p1 - mu.p1) * (s1 - mu.s1))/(length(s1) - 
            1)
        gamma2 <- sum((p2 - mu.p2) * (s2 - mu.s2))/(length(s2) - 
            1)
        mu.E <- (mu.s1 + mu.s2)/2
        delta <- matrix(c(1/mu.p1, -(mu.s2/mu.p1), 1/mu.p2, -(mu.s1/mu.p2)), 
            4, 1)
        Sigma.11 <- sigma.p1^2 * sigma.s2^2 + sigma.p1^2 * mu.s2^2 + 
            sigma.s2^2 * mu.p1^2
        Sigma.12 <- mu.s2 * sigma.p1^2
        Sigma.13 <- gamma1 * gamma2 + gamma1 * mu.p2 * mu.s2 + 
            gamma2 * mu.p1 * mu.s1
        Sigma.14 <- gamma2 * mu.p1
        Sigma.22 <- sigma.p1^2
        Sigma.23 <- gamma1 * mu.p2
        Sigma.24 <- 0
        Sigma.33 <- sigma.p2^2 * sigma.s1^2 + sigma.p2^2 * mu.s1^2 + 
            sigma.s1^2 * mu.p2^2
        Sigma.34 <- mu.s1 * sigma.p2^2
        Sigma.44 <- sigma.p2^2
        Sigma <- matrix(c(Sigma.11, Sigma.12, Sigma.13, Sigma.14, 
            Sigma.12, Sigma.22, Sigma.23, Sigma.24, Sigma.13, 
            Sigma.23, Sigma.33, Sigma.34, Sigma.14, Sigma.24, 
            Sigma.34, Sigma.44), 4, 4)
        sigma.E <- 0.5 * sqrt(((t(delta) %*% Sigma) %*% delta)/length(p1))
        list(p.value = (pnorm((E - mu.E)/sigma.E, lower.tail = FALSE)),
	    z.score = (E - mu.E)/sigma.E, E.gene = length(p1)/2 * 
            ((p1 * s2)/sum(p1) + (p2 * s1)/sum(p2)))
    }

#     if (pValue.query) 
#          s1 <- -log(query.p, 10)
	  s1 <- -log10(query.p[, 2])
#     if (pValue.reference)
# 	reference.p[reference.p<1e-10]<-1e-10
      s2 <- reference.p[, "Score"]
#        s2 <- -log(reference.p, 10)
#     else s2 <- reference.p[, "Score"]
#     if (is.null(pAdjust.method.query)) {
#         if (pValue.query) {
# print("#1")
             query.m0 <- estimateH0(query.p[,2], lambda)
# print("#1")
             p1 <- postP(query.p[,2], scale = query.m0)
# print("#1")
             query.m0 <- query.m0/(query.m0 + 1)
# print("#1")
#         }
#         else {
#              p1 <- query.p[, 2]
#              query.m0 <- estimateH0(p1, lambda)
#              query.m0 <- query.m0/(query.m0 + 1)
#         }
#     }
#     else {
#         p1 <- 1 - p.adjust(query.p, method = pAdjust.method.query)
#         query.m0 <- NULL
#     }
#     if (is.null(pAdjust.method.reference)) {
#         if (pValue.reference) {
#             reference.m0 <- estimateH0(reference.p, lambda)
#             p2 <- postP(reference.p, scale = reference.m0)
#             reference.m0 <- reference.m0/(reference.m0 + 1)
#         }
#         else {
             p2 <- reference.p[, "Prob"]
             reference.m0 <- estimateH0(1 - p2, lambda)
             reference.m0 <- reference.m0/(reference.m0 + 1)
#         }
#     }
#     else {
#         p2 <- 1 - p.adjust(reference.p, method = pAdjust.method.reference)
#         reference.m0 <- NULL
#     }
   if(rescale){
##	require(preprocessCore)
    	sNorm <- preprocessCore::normalize.quantiles(cbind(s1,s2))
	s1 <- sNorm[,1]
	s2 <- sNorm[,2]
	pNorm <- preprocessCore::normalize.quantiles(cbind(p1,p2))
	p1<-pNorm[,1]
	p2<-pNorm[,2]
	if(plotRescaled){
		X11()
 		par(mfrow=c(1,2))
		plot(s1,s2)
		plot(p1,p2)
	}
   }
# print("#1")
    res <- RSCM(p1, p2, s1, s2)
#    geneTable <- data.frame(geneID = commonGenes, E.value = res$E.gene)
#     ret <- list(p.value = res$p.value, z.score = res$z.score, query.m0 = query.m0, reference.m0 = reference.m0,junk=data.frame(p1,p2,s1,s2))
    ret <- list(p.value = res$p.value, z.score = res$z.score, query.m0 = query.m0, reference.m0 = reference.m0)
     ret
}
