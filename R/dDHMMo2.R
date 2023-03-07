library(nimble)

dDHMMo2 <- nimbleFunction(
    run = function(x = double(1),
                   init = double(1),
                   f = double(),
                   probTrans1  = double(2),
                   probTrans2  = double(2),
                   probTrans3  = double(2),
                   probTrans4  = double(2),
                   probTrans5  = double(2),
                   probTrans6  = double(2),
                   probTrans7  = double(2),
                   probTrans8  = double(2),
                   probTrans9  = double(2),
                   probTrans10 = double(2),
                   probTrans11 = double(2),
                   probTrans12 = double(2),
                   probObs1 = double(2),
                   probObs = double(3),
                   len = double(),
                   checkRowSums = double(0, default = 1),
                   log = integer(0, default = 0)) {
        ##
        ##if (length(init) != dim(probObs)[1]) stop("In dDHMMo: Length of init does not match ncol of probObs in dDHMMo.")
        ##if (length(init) != dim(probTrans)[1]) stop("In dDHMMo: Length of init does not match dim(probTrans)[1] in dDHMMo.")
        ##if (length(init) != dim(probTrans)[2]) stop("In dDHMMo: Length of init does not match dim(probTrans)[2] in dDHMMo.")
        ##if (length(x) != len) stop("In dDHMMo: Length of x does not match len in dDHMM.")
        ##if (len - 1 > dim(probTrans)[3]) stop("In dDHMMo: dim(probTrans)[3] does not match len - 1 in dDHMMo.")
        ##if (len != dim(probObs)[3]) stop("In dDHMMo: dim(probObs)[3] does not match len in dDHMMo.")
        ##if (abs(sum(init) - 1) > 1e-6) stop("In dDHMMo: Initial probabilities must sum to 1.")
        ##
        ##if (checkRowSums) {
        ##    transCheckPasses <- TRUE
        ##    for (i in 1:dim(probTrans)[1]) {
        ##        for (k in 1:dim(probTrans)[3]) {
        ##            thisCheckSum <- sum(probTrans[i,,k])
        ##            if (abs(thisCheckSum - 1) > 1e-6) {
        ##                ## Compilation doesn't support more than a simple string for stop()
        ##                ## so we provide more detail using a print().
        ##                print("In dDHMMo: Problem with sum(probTrans[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
        ##                transCheckPasses <- FALSE
        ##            }
        ##        }
        ##    }
        ##    obsCheckPasses <- TRUE
        ##    for (i in 1:dim(probObs)[1]) {
        ##        for (k in 2:dim(probObs)[3]) {  ## THIS IS CHANGED !! "1:" is now "2:"
        ##            thisCheckSum <- sum(probObs[i,,k])
        ##            if (abs(thisCheckSum - 1) > 1e-6) {
        ##                print("In dDHMMo: Problem with sum(probObs[i,,k]) with i = ", i, " k = ", k, ". The sum should be 1 but is ", thisCheckSum)
        ##                obsCheckPasses <- FALSE
        ##            }
        ##        }
        ##    }
        ##    if(!(transCheckPasses | obsCheckPasses))
        ##        stop("In dDHMMo: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
        ##    if(!transCheckPasses)
        ##        stop("In dDHMMo: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
        ##    if(!obsCheckPasses)
        ##        stop("In dDHMMo: probObs was not specified correctly. Probabilities in each row must sum to 1.")
        ##}
        ##
        pi <- init
        logL <- 0
        nObsClasses <- dim(probObs)[2]
        lengthX <- length(x)
        for (t in 1:lengthX) {
            if (x[t] > nObsClasses | x[t] < 1) stop("In dDHMMo: Invalid value of x[t].")
            if(t==1) Zpi <- probObs1[, x[t]]    * pi
            if(t >1) Zpi <- probObs [, x[t], t] * pi
            sumZpi <- sum(Zpi)
            logL <- logL + log(sumZpi)
            ####original: pi <- ((Zpi %*% probTrans[,,t])/sumZpi)[1, ]
            if(t != lengthX) {
                whichPT <- f + t - 1
                if(whichPT == 1)    pi <- ((Zpi %*% probTrans1[,])/sumZpi)[1, ]
                if(whichPT == 2)    pi <- ((Zpi %*% probTrans2[,])/sumZpi)[1, ]
                if(whichPT == 3)    pi <- ((Zpi %*% probTrans3[,])/sumZpi)[1, ]
                if(whichPT == 4)    pi <- ((Zpi %*% probTrans4[,])/sumZpi)[1, ]
                if(whichPT == 5)    pi <- ((Zpi %*% probTrans5[,])/sumZpi)[1, ]
                if(whichPT == 6)    pi <- ((Zpi %*% probTrans6[,])/sumZpi)[1, ]
                if(whichPT == 7)    pi <- ((Zpi %*% probTrans7[,])/sumZpi)[1, ]
                if(whichPT == 8)    pi <- ((Zpi %*% probTrans8[,])/sumZpi)[1, ]
                if(whichPT == 9)    pi <- ((Zpi %*% probTrans9[,])/sumZpi)[1, ]
                if(whichPT == 10)   pi <- ((Zpi %*% probTrans10[,])/sumZpi)[1, ]
                if(whichPT == 11)   pi <- ((Zpi %*% probTrans11[,])/sumZpi)[1, ]
                if(whichPT == 12)   pi <- ((Zpi %*% probTrans12[,])/sumZpi)[1, ]
            }
        }
        returnType(double())
        if (log) return(logL)
        return(exp(logL))
    }
)






rDHMMo2 <- nimbleFunction(
    run = function(n = integer(),
                   init = double(1),
                   f = double(),
                   probTrans1  = double(2),
                   probTrans2  = double(2),
                   probTrans3  = double(2),
                   probTrans4  = double(2),
                   probTrans5  = double(2),
                   probTrans6  = double(2),
                   probTrans7  = double(2),
                   probTrans8  = double(2),
                   probTrans9  = double(2),
                   probTrans10 = double(2),
                   probTrans11 = double(2),
                   probTrans12 = double(2),
                   probObs1 = double(2),
                   probObs = double(3),
                   len = double(),
                   checkRowSums = double(0, default = 1)) {
        ans <- rep(1,2)
        returnType(double(1))
        return(ans)
    })


registerDistributions(list(
    dDHMMo2 = list(
        BUGSdist = "dDHMMo2(init, f, probTrans1, probTrans2, probTrans3, probTrans4, probTrans5, probTrans6, probTrans7, probTrans8, probTrans9, probTrans10, probTrans11, probTrans12, probObs1, probObs, len, checkRowSums)",
        Rdist    = "dDHMMo2(init, f, probTrans1, probTrans2, probTrans3, probTrans4, probTrans5, probTrans6, probTrans7, probTrans8, probTrans9, probTrans10, probTrans11, probTrans12, probObs1, probObs, len, checkRowSums)",
        discrete = TRUE,
        types = c('value = double(1)',
                  'init = double(1)',
                  'probTrans1  = double(2)',
                  'probTrans2  = double(2)',
                  'probTrans3  = double(2)',
                  'probTrans4  = double(2)',
                  'probTrans5  = double(2)',
                  'probTrans6  = double(2)',
                  'probTrans7  = double(2)',
                  'probTrans8  = double(2)',
                  'probTrans9  = double(2)',
                  'probTrans10 = double(2)',
                  'probTrans11 = double(2)',
                  'probTrans12 = double(2)',
                  'probObs1 = double(2)',
                  'probObs = double(3)',
                  'len = double()',
                  'checkRowSums = double(0)'),
        mixedSizes = TRUE,
        pqAvail = FALSE))
    )
























