
#' Ekbohm MLE based t test
#'
#' @param x numerical vector, should be same length as y
#' @param y numerical vector, should be same length as x
#' @param alternative 'Greater', 'Less', or 'Two.Sided' based on direction of your test
#'
#' @returnreturns warnings if input criteria is not met, or if results
#' may be ineffective, returns test statistic and p value if all is correct
#' @export
#'
#' @examples x <- c(6,9,12,13,15,13,17,NA,13,18,15) y <- c(10,NA,16,15,14,12,12,17,NA,11,12)
#' ekbohmMLE(x,y,alternative = 'two.sided)
ekbohmMLE<- function(x,y, alternative = 'two.sided'){
  # paired samples
  n1.x <- x[(!is.na(x)) & (!is.na(y))]
  n1.y <- y[(!is.na(x)) & (!is.na(y))]
  #tumor and normal samples
  n2 <- x[!is.na(x) & is.na(y)]
  n3 <- y[!is.na(y) & is.na(x)]

  n1len <- length(n1.x)
  n2len <- length(n2)
  n3len <- length(n3)

  Tbar <- mean(n2)
  Nbar <- mean(n3)
  ST <- sd(n2)
  SN <- sd(n3)
  T1bar <- mean(n1.x)
  N1bar <- mean(n1.y)
  ST1 <- sd(n1.x)
  SN1 <- sd(n1.y)
  STN1 <- cov(n1.x, n1.y)

  r <- STN1 / (SN1 * ST1)
  fstar <-n1len * (n1len + n3len + n2len * r) * (((n1len + n2len) * (n1len + n3len) - n2len * n3len * r^2)^-1)
  gstar <-n1len * (n1len + n2len + n3len * r) * (((n1len + n2len) * (n1len + n3len) - n2len * n3len * r^2)^-1)
  sigmahat2 <- (ST1^2 * (n1len  - 1) + SN1^2 * (n1len - 1) + (1 + r^2) * (ST^2 *(n2len - 1) + SN^2 * (n3len - 1))) / (2 * (n1len - 1) + (1 + r^2) * (n2len + n3len - 2))
  V1star <-  sigmahat2 * ((2 * n1len * (1 - r) + (n2len + n3len) * (1 - r^2))) / ((n1len + n2len) * (n1len + n3len) - n2len *n3len * r^2)

  ZE <- (fstar * (T1bar - Tbar) - gstar * (N1bar - Nbar) + Tbar -Nbar) / sqrt(V1star)
  # test variances
  testvarX <- x[!is.na(x)]
  testvarY <- y[!is.na(y)]

  if(var.test(testvarX,testvarY)$p.value <= 0.05){
    cat("warning: the sample variances are significantly different \n",
        "Ekbohm test is not recommended \n \n")
  }
  {if(length(n2) == 0 & length(n3) ==0){
    cat(' all samples have a match \n using paired t test \n',
        'p value is ',t.test(x,y, paired = T, alternative = 'two.sided')$p.value)
  }
    else if ((length(x) != length(y)) | (length(n1.x) == 0) | ((length(n3)<=1 | length(n2)<=1))) {
      cat("warning: length of vectors are unequal\n or no paired samples",
          "\n or length of n2 or n3 is 1 or less", "(", length(n2),",",length(n3),")")
    }}

  {if (alternative == "two.sided") {
    pvalue = 2 * pt(abs(ZE), n1len, lower.tail = F)
    cat("\n Ekbohm MLE based t test \n two sided, p value is  ",pvalue,
        "\n ZE is ", ZE)
  }
    else if (alternative == "greater") {
      pvalue = pt(ZE, n1len, lower.tail = F)
      cat("\n Ekbohm MLE based t test \n upper tail, p value is  ",pvalue,
          "\n ZE is ", ZE)
    }
    else if (alternative == "less") {
      pvalue = pt(ZE, n1len, lower.tail = T)
      cat("\n Ekbohm MLE based t test \n lower tail, p value is  ",pvalue,
          "\n ZE is ", ZE)
    }
  }
}

#' Kim et al Modified T-statistic
#'
#' @param x numerical vector, should be same length as y
#' @param y numerical vector, should be same length as x
#' @param alternative 'Greater', 'Less', or 'Two.Sided' based on direction of your test
#'
#' @return returns warnings if input criteria is not met, or if results
#' may be ineffective, returns test statistic and p value if all is correct
#' @export
#'
#' @examples x <- c(6,13,12,NA,15,13,17,NA,13,18,15) y <- c(10,NA,16,15,14,12,12,17,NA,11,12)
#' modifiedTstatistic(x,y,alternative = 'two.sided)
modifiedTstatistic <- function(x, y, alternative = 'two.sided'){
  # edit vectors to fit n1, n2, n3
  n1.x <- x[(!is.na(x)) & (!is.na(y))]
  n1.y <- y[(!is.na(x)) & (!is.na(y))]
  n2 <- x[!is.na(x) & is.na(y)]
  n3 <- y[!is.na(y) & is.na(x)]

  D <- n1.x - n1.y
  Dbar <- mean(n1.x-n1.y)
  Tbar <- mean(n2)
  Nbar <- mean(n3)
  nh <- 2 /((1/length(n2))+((1/length(n3))))
  VD <- var(D)
  VT <- var(n2)
  VN <- var(n3)

  t3 <-(length(n1.x) * Dbar + nh * (Tbar - Nbar)) / sqrt(length(n1.x)* VD + ((nh ^ 2) * (VN /length(n3) + VT / length(n2))))

  if(length(n2) == 0 & length(n3) ==0){
    cat(' all samples have a match \n using paired t test \n',
        'p value is ',t.test(x,y, paired = T, alternative = 'two.sided')$p.value)
  }
  else if ((length(x) != length(y)) | (length(n1.x) == 0) | ((length(n3)<=1 | length(n2)<=1))) {
    cat("warning: length of vectors are unequal\n or no paired samples",
        "\n or length of n2 or n3 is 1 or less", "(", length(n2),",",length(n3),")")
  }else {
    if (alternative == "two.sided") {
      pvalue = 2 * pnorm(abs(t3), lower.tail = F)
      cat("Modified t test \n two sided, p value is  \n",pvalue,
          "\n t3 value is ", t3)
    }
    if (alternative == "greater") {
      pvalue = pnorm(t3, lower.tail = F)
      cat("Modified t test \n upper tail, p value is  ",pvalue,
          "\n t3 value is ", t3)
    }
    if (alternative == "less") {
      pvalue = pnorm(t3, lower.tail = T)
      cat("Modified t test \n lower tail, p value is  ",pvalue,
          "\n t3 value is ", t3)
    }
  }
}

#' Corrected Z Test by Looney and Jones
#'
#' @param x numerical vector, should be same length as y
#' @param y numerical vector, should be same length as x
#' @param alternative 'Greater', 'Less', or 'Two.Sided' based on direction of your test
#'
#' @return returns warnings if input criteria is not met, or if results
#' may be ineffective, returns test statistic and p value if all is correct
#' @export
#'
#' @examples x <- c(6,9,12,15,NA,13,17,NA,25,18,15) y <- c(10,NA,16,15,14,12,12,17,NA,11,12)
#' correctedZtest(x,y,alternative = 'two.sided)
correctedZtest <- function(x,y,alternative = 'two.sided'){

  n1.x <- x[(!is.na(x)) & (!is.na(y))]
  n1.y <- y[(!is.na(x)) & (!is.na(y))]
  n2 <- x[!is.na(x) & is.na(y)]
  n3 <- y[!is.na(y) & is.na(x)]

  Tbarstar <- mean(c(n1.x,n2))
  Nbarstar <- mean(c(n1.y,n3))
  TVarStar <- var(c(n1.x,n2))
  NVarStar <- var(c(n1.y,n3))
  n1cov <- cov(n1.x,n1.y)

  zcorr<-(Tbarstar-Nbarstar)/sqrt(TVarStar/(length(n1)+length(n2))+NVarStar/(length(n1)+length(n3))-2*length(n1)*n1cov/((length(n1)+length(n2))*(length(n1)+length(n3))))

  cat("Zcorr is ", zcorr, "\n\n")

  if(length(n2) == 0 & length(n3) ==0){
    cat(' all samples have a match \n using paired t test \n',
        'p value is ',t.test(x,y, paired = T, alternative = 'two.sided')$p.value)
  }
  else if ((length(x) != length(y)) | (length(n1.x) == 0) | ((length(n3)<=1 | length(n2)<=1))) {
    cat("warning: length of vectors are unequal\n or no paired samples",
        "\n or length of n2 or n3 is 1 or less", "(", length(n2),",",length(n3),")",
        "\n use t.test")
  }else {

    if (alternative == "two.sided") {
      pvalue = 2 * (pnorm(abs(zcorr), lower.tail = F))
      cat("Corrected Z test \n two sided, p value is  ",pvalue)
    }
    else if (alternative == "greater") {
      pvalue = pnorm(zcorr, lower.tail = F)
      cat("Corrected Z test \n upper tail, p value is  ",pvalue)
    }
    else if (alternative == "less") {
      pvalue = pnorm(zcorr, lower.tail = T)
      cat("Corrected Z test \n lower tail, p value is  ",pvalue)
    }
  }
}

#' weighted Z test
#'
#' @param x numerical vector, should be same length as y
#' @param y numerical vector, should be same length as x
#' @param alternative 'Greater', 'Less', or 'Two.Sided' based on direction of your test
#'
#' @return returns warnings if input criteria is not met, or if results
#' may be ineffective, returns test statistic and p value if all is correct
#' @export
#'
#' @examples x <- c(6,9,12,15,NA,13,17,NA,25,18,15) y <- c(10,NA,16,15,14,12,12,17,NA,11,12)
#' weightedZtest(x,y,alternative = 'two.sided)
weightedZtest <- function(x, y, alternative = 'two.sided'){
  # n1x and n1y represent paired samples
  n1.x <- x[(!is.na(x)) & (!is.na(y))]
  n1.y <- y[(!is.na(x)) & (!is.na(y))]
  # n2 and n3 represent tumor and normal
  n2 <- x[!is.na(x) & is.na(y)]
  n3 <- y[!is.na(y) & is.na(x)]
  if(length(n2) == 0 & length(n3) ==0){
    cat(' all samples have a match \n use paired t test',
        'p value is ',t.test(x,y, paired = T, alternative = 'two.sided')$p.value)
  }
  else if ((length(x) != length(y)) | (length(n1.x) == 0) | ((length(n3)<=1 | length(n2)<=1))) {
    cat("warning: length of vectors are unequal\n or no paired samples",
        "\n or length of tumor or normal samples is 1 or less \n",
        "use t.test")
  }else {
    # weights
    w1 <- sqrt(2*length(n1.x))
    w2 <- sqrt(length(n2)+length(n3))

    # must use one sided test to avoid bias
    P1i <- t.test(n1.x,n1.y, paired = T, alternative = 'greater')$p.value
    P2i <- t.test(n2,n3, alternative = 'greater')$p.value

    Z1i <- qnorm(1-P1i)
    Z2i <- qnorm(1-P2i)

    PCi <- 1 - pnorm((w1*Z1i + w2*Z2i)/sqrt(w1^2 + w2^2))


    if(alternative == 'two.sided'){
      if(PCi < 0.5){
        pvalue <- 2*PCi
        cat("using weighted Z method\n",
            'p value is ', pvalue)
      }else{
        pvalue <- 2*(1-PCi)
        cat("using weighted Z method\n",
            'p value is ', pvalue)
      }
    }else if(alternative == 'greater'){
      pvalue <- PCi
      cat("using weighted Z method\n",
          'p value is ', pvalue)
    }else if(alternative == 'less'){
      pvalue <- 1 - PCi
      cat("using weighted Z method\n",
          'p value is ', pvalue)
    }
  }
  print(cat('\n \n length of (tumor, normal)',"(", length(n2),",",length(n3),")","\n"))
}

#' Lin And Stivers MLE based t test
#'
#' @param x numerical vector, should be same length as y
#' @param y numerical vector, should be same length as x
#' @param alternative 'Greater', 'Less', or 'Two.Sided', based on the direction of your test
#'
#' @return returns warnings if input criteria is not met, or if results
#' may be ineffective, returns test statistic and p value if all is correct
#' @export
#'
#' @examples x <- c(6,9,12,15,NA,13,17,NA,13,18,NA) y <- c(10,NA,16,15,14,12,12,17,NA,11,12)
#' linstiversMLE(x,y,alternative = 'two.sided')

linstiversMLE <- function(x,y,alternative = 'two.sided'){
  # paired samples
  n1.x <- x[(!is.na(x)) & (!is.na(y))]
  n1.y <- y[(!is.na(x)) & (!is.na(y))]
  #n2 represents tumor sample, n3 represents normal sample
  n2 <- x[!is.na(x) & is.na(y)]
  n3 <- y[!is.na(y) & is.na(x)]

  n1len <- length(n1.x)
  n2len <- length(n2)
  n3len <- length(n3)

  Tbar <- mean(n2)
  Nbar <- mean(n3)
  ST <- sd(n2)
  SN <- sd(n3)
  T1bar <- mean(n1.x)
  N1bar <- mean(n1.y)
  ST1 <- sd(n1.x)
  SN1 <- sd(n1.y)
  STN1 <- cov(n1.x, n1.y)

  r <- STN1 / (SN1 * ST1)
  f <- (n1len*(n1len + n3len + (n2len*STN1)/ST1))*(((n1len+n2len)*(n1len+n3len)-(n2len*n3len*r^2))^(-1))
  g <- (n1len*(n1len + n2len + (n3len*STN1)/ST1))*(((n1len+n2len)*(n1len+n3len)-(n2len*n3len*r^2))^(-1))

  Vp1 <- (((f^2)/n1len)+(((1-f)^(2))/n2len))*ST1*(n1len-1)+(((g^2)/n1len)+((1-g)^(2)/n3len))
  Vp2 <- SN1*(n1len - 1) - 2*f*g*STN1*(n1len - 1)/n1len
  V1 <- (Vp1 + Vp2)/(n1len - 1)

  ZLS <- ((f*(T1bar - Tbar)) - (g*(N1bar - Nbar)) + Tbar - Nbar) / sqrt(V1)
  # test variances
  testvarX <- x[!is.na(x)]
  testvarY <- y[!is.na(y)]

  if(var.test(testvarX,testvarY)$p.value >= 0.05){
    cat("warning: the sample variances are similar \n",
        "Ekbohm may be recommended instead of Lin and Stivers \n \n")
  }

  {if(length(n2) == 0 & length(n3) ==0){
    cat(' all samples have a match \n using paired t test \n',
        'p value is ',t.test(x,y, paired = T, alternative = 'two.sided')$p.value)
  }
    else if ((length(x) != length(y)) | (length(n1.x) == 0) | ((length(n3)<=1 | length(n2)<=1))) {
      cat("warning: length of vectors are unequal\n or no paired samples",
          "\n or length of tumor or normal is 1 or less", "(", length(n2),",",length(n3),")",
          "\n use t.test")
    }}

  {if (alternative == "two.sided") {
    pvalue = 2 * pt(abs(ZLS), df = n1len, lower.tail = F)
    cat("\n lin stivers MLE based t test \n two sided, p value is  ",pvalue,
        "\n ZLS is ", ZLS)
  }
    else if (alternative == "greater") {
      pvalue = pt(ZLS, df = n1len, lower.tail = F)
      cat("\n lin stivers MLE based t test \n upper tail, p value is  ",pvalue,
          "\n ZLS is ", ZLS)
    }
    else if (alternative == "less") {
      pvalue = pt(ZLS,df = n1len, lower.tail = T)
      cat("\n lin stivers MLE based t test \n lower tail, p value is  ",pvalue,
          "\n ZLS is ", ZLS)
    }
  }
}

