
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
#' @keywords ...
#' @export 
#' @examples
#' res <- batchTregGRS(....)

GRSlocal <- function (query.p, reference.p, query.gL = NULL, reference.gL = NULL, estimateNullDistr = TRUE, nullDistrQuantiles = c(0.9, 0.95,0.99), nullDistrN = 100)
 {
    pAdjust <- function(p, method = "top", pAdjustN = 100, sigLevel = 0.1) {
        if (method == "top" & !is.null(pAdjustN)) {
            index <- rank(p, ties.method = "random") > pAdjustN
            p[index] <- 1
            p[!index] <- sigLevel
            p
        }
        else {
            p.adjust(p, method)
        }
    }
    postP <- function(p, scale = 1) {
        B <- -exp(1) * p * log(p)
        B[p > exp(-1)] <- 1
        1 - 1/(1 + 1/(max(1, scale) * B))
    }
    estimateH0 <- function(p) {
        pi0 <- 1
        tryCatch(pi0 <- qvalue(p)$pi0, error = function(e) {
            pi0 <- 1
        })
        pi0 <- min(pi0, 0.999)
        pi0/(1 - pi0)
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
        list(p.value = (1 - pnorm(abs(E - mu.E)/sigma.E)) * 2,
            z.score = abs(E - mu.E)/sigma.E, E.gene = length(p1)/2 * ((p1 *
                s2)/sum(p1) + (p2 * s1)/sum(p2)))
    }
    if (!is.null(dim(query.p)) | !is.null(query.gL)) {
        if (!is.null(dim(query.p))) {
            query.gL <- query.p[, 1]
            query.p <- as.numeric(query.p[, 2])
            if (is.null(dim(reference.p))) {
                warning("parameter reference.p should be a matrix or data.frame")
            }
            else {
                reference.gL <- reference.p[, 1]
                reference.p <- as.numeric(reference.p[, 2])
            }
            uniqueGenes <- unique(query.gL)
            if (length(uniqueGenes) != length(query.gL)) {
                averagePvalues <- split(query.p, query.gL)
                names(averagePvalues) <- as.character(names(averagePvalues))
                averagePvalues <- unlist(lapply(averagePvalues,
                  function(x) exp(mean(log(x), na.rm = TRUE))))
                query.gL <- uniqueGenes
                query.p <- averagePvalues[match(query.gL, names(averagePvalues))]
            }
        }
        if (is.null(reference.gL)) {
            warning("parameter reference.gL should not be NULL")
        }
        else {
            uniqueGenes <- unique(reference.gL)
            if (length(uniqueGenes) != length(reference.gL)) {
                averagePvalues <- split(reference.p, reference.gL)
                names(averagePvalues) <- as.character(names(averagePvalues))
                averagePvalues <- unlist(lapply(averagePvalues,
                  function(x) exp(mean(log(x), na.rm = TRUE))))
                reference.gL <- uniqueGenes
                reference.p <- averagePvalues[match(reference.gL,
                  names(averagePvalues))]
            }
        }
        commonGenes <- intersect(query.gL, reference.gL)
        if (length(commonGenes) < 2)
            #warning("number of common gene IDs is < 2")
            return("number of common gene IDs is < 2")
        query.p <- query.p[match(commonGenes, query.gL)]
        reference.p <- reference.p[match(commonGenes, reference.gL)]
    }
    s1 <- -log(query.p, 10)
    s2 <- -log(reference.p, 10)
    p1 <- postP(query.p, scale = estimateH0(query.p))
    p2 <- postP(reference.p, scale = estimateH0(reference.p))
    res <- RSCM(p1, p2, s1, s2)
    geneTable <- data.frame(geneID = commonGenes, E.value = res$E.gene)
    ret <- list(p.value = res$p.value, z.score = res$z.score,
        geneTable = geneTable[order(res$E.gene, decreasing = TRUE),
            ])
    if (estimateNullDistr) {
        f <- function(p1, p2, s1, s2) {
            index1 <- sample(length(p1))
            index2 <- sample(length(p2))
            RSCM(p1[index1], p2[index2], s1[index1], s2[index2])$E.gene
        }
        distr <- sapply(1:nullDistrN, function(i) quantile(f(p1,
            p2, s1, s2), probs = nullDistrQuantiles))
        if (length(nullDistrQuantiles) == 1)
            distr <- t(distr)
        distr <- rowMeans(distr)
        ret <- c(ret, list(EvalueNullDistrQ = distr))
    }
    ret
}
