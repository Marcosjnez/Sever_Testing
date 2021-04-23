library(shiny)
library(shinydashboard)

shinyServer(function(input, output, session) {
  Sever.Z.one <- function(n, dif, sigma, alpha, plot=TRUE, axis=FALSE) {
    validate(
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample sizes cannot be lower than 2"))
    se <- sigma/sqrt(n)
    z <- dif/se
    pvalue <- 1 - pnorm(z)
    crit.v <- qnorm(1-alpha)
    me <- crit.v*se
    ci <- dif - me
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(pvalue <= alpha) {
        find <- function(discrepancy, Severity) (pnorm(z, discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (pnorm(z, -discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- suppressWarnings(optimize(f=find, Severity=0.999, interval=interval, maximum=F, tol = 1e-9))[[1]]
        Discrepancy2 <- -suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval2, maximum=F, tol = 1e-9))[[1]]
        symbol <- bquote(.("SEV(") ~ mu > gamma~.(")"))
        curve(pnorm(z, x/se), xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, xaxt='n', yaxt='n', ylab=NA, xlab=NA, cex.main=1.5)
      }
      else{
        find <- function(discrepancy, Severity) (1-pnorm(z, -discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (1-pnorm(z, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        if(z>=0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        if(z<0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        symbol <- bquote(.("SEV(") ~ mu <= gamma~.(")"))
        curve(1 - pnorm(z, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, ylab=NA, xlab=NA, xaxt='n', yaxt='n')     
      }
      title(ylab="Severity / Power", cex.lab=1.6, line = 3.7)
      title(xlab=bquote(gamma == mu[1] - mu[0]), cex.lab=1.5, line = 3.5)
      Discrepancies=c(Discrepancy, Discrepancy2)
      Discrepancy=Discrepancies[which.min(Discrepancies)]; Discrepancy2=Discrepancies[which.max(Discrepancies)]
      curve(1 - pnorm(crit.v, x/se), add=T, col="#F36565", lwd=4, xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), xaxt='n', yaxt='n', xlab=NA, ylab=NA)
      text(x=Discrepancy2-(Discrepancy2-Discrepancy)*0.15, y=0.7, labels = bquote(hat(mu) == .(round(dif, 2))), cex=1.3)
      axis(2, seq(0, 1, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from, input$to, input$by), las=1, cex.axis=1.3)
      }
      else {
        if(sigma>=1) {
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/8.5, 1))
          axis(1, seq, las=1, cex.axis=1.3)   
        }
        else {
          x <- nchar(strsplit(as.character(sigma), "\\.")[[1]][2])
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/10, x+1))
          axis(1, seq, las=1, cex.axis=1.3)
        }
      }}
    data <- data.frame(Lower.Limit=ci, Upper.Limit=Inf, Z.critic=crit.v, Z.value=z, p.value=pvalue)
    return(data)
  }
  Sever.t.one <- function(n, dif, sigma, alpha, plot=TRUE, axis=FALSE) {
    validate(
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample sizes cannot be lower than 2"))
    se <- sigma/sqrt(n)
    t <- dif/se
    pvalue <- 1 - pt(t, df = n-1)
    crit.v <- qt(1-alpha, df = n-1)
    me <- crit.v*se
    ci <- dif - me   
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(pvalue <= alpha) {
        find <- function(discrepancy, Severity) (pt(t, n-1, discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (pt(t, n-1, -discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- -suppressWarnings(optimize(f=find2, Severity=0.999, interval=interval, maximum=F, tol = 1e-9))[[1]]
        Discrepancy2 <- -suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
        symbol <- bquote(.("SEV(") ~ mu > gamma~.(")"))
        curve(pt(t, n-1, x/se), xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)
      }
      else if(round(pvalue, 5)==1 & t<0) {
        symbol <- bquote(.("SEV(") ~ mu <= gamma~.(")"))
        find2 <- function(discrepancy, Severity) (1-pt(t, n-1, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- suppressWarnings(optimize(f=find2, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
        Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        curve(1 - pt(t, n-1, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)     
      }
      else{
        symbol <- bquote(.("SEV(") ~ mu <= gamma~.(")"))
        find <- function(discrepancy, Severity) (1-pt(t, n-1, -discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (1-pt(t, n-1, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        if(t>=0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        if(t<0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        curve(1 - pt(t, n-1, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)     
      }
      title(ylab="Severity / Power", cex.lab=1.6, line = 3.7)
      title(xlab=bquote(gamma == mu[1] - mu[0]), cex.lab=1.5, line = 3.5)
      Discrepancies=c(Discrepancy, Discrepancy2)
      Discrepancy=Discrepancies[which.min(Discrepancies)]; Discrepancy2=Discrepancies[which.max(Discrepancies)]
      curve(1 - pt(crit.v, n-1, x/se), add=T, col="#F36565", lwd=4, xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), xaxt='n', yaxt='n', xlab=NA, ylab=NA)
      text(x=Discrepancy2-(Discrepancy2-Discrepancy)*0.15, y=0.7, labels = bquote(hat(mu) == .(round(dif, 2))), cex=1.3)
      axis(2, seq(0, 1, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from, input$to, input$by), las=1, cex.axis=1.3)
      }
      else {
        if(sigma>=1) {
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/8.5, 1))
          axis(1, seq, las=1, cex.axis=1.3)   
        }
        else {
          x <- nchar(strsplit(as.character(sigma), "\\.")[[1]][2])
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/10,x+1))
          axis(1, seq, las=1, cex.axis=1.3)
        }   
      }
    }
    data <- data.frame(Lower.Limit=ci, Upper.Limit=Inf, t.critic=crit.v, t.value=t, p.value=pvalue)
    return(data)
  }
  Sever.Z.p <- function(n, dif, sigma, sigma.a, sigma.b, cor, alpha, plot=TRUE, axis=FALSE) {
    if(input$select=='cor'| input$select2=='cor2') sigma <- sqrt(sigma.a^2+sigma.b^2-2*cor*sigma.a*sigma.b)
    validate(
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample size cannot be lower than 2"))
    se <- sqrt(sigma^2/n)
    z <- dif/se
    pvalue <- 1 - pnorm(z)
    crit.v <- qnorm(1-alpha)
    me <- crit.v*se
    ci <- dif - me
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(pvalue <= alpha) {
        find <- function(discrepancy, Severity) (pnorm(z, discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (pnorm(z, -discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- suppressWarnings(optimize(f=find, Severity=0.999, interval=interval, maximum=F, tol = 1e-9))[[1]]
        Discrepancy2 <- -suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval2, maximum=F, tol = 1e-9))[[1]]
        symbol <- bquote(.("SEV(") ~ mu > gamma~.(")"))
        curve(pnorm(z, x/se), xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)
      }
      else{
        find <- function(discrepancy, Severity) (1-pnorm(z, -discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (1-pnorm(z, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        if(z>=0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        if(z<0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }     
        symbol <- bquote(.("SEV(") ~ mu <= gamma~.(")"))
        curve(1 - pnorm(z, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)     
      }
      title(ylab="Severity / Power", cex.lab=1.6, line = 3.7)
      title(xlab=bquote(gamma == mu[1] - mu[0]), cex.lab=1.5, line = 3.5)
      Discrepancies=c(Discrepancy, Discrepancy2)
      Discrepancy=Discrepancies[which.min(Discrepancies)]; Discrepancy2=Discrepancies[which.max(Discrepancies)]
      curve(1 - pnorm(crit.v, x/se), add=T, col="#F36565", lwd=4, xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), xaxt='n', yaxt='n', xlab=NA, ylab=NA)
      text(x=Discrepancy2-(Discrepancy2-Discrepancy)*0.15, y=0.7, labels = bquote(hat(mu) == .(round(dif, 2))), cex=1.3)
      axis(2, seq(0, 1, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from2, input$to2, input$by2), las=1, cex.axis=1.3)
      }
      else {
        if(sigma>=1) {
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/8.5, 1))
          axis(1, seq, las=1, cex.axis=1.3)   
        }
        else {
          x <- nchar(strsplit(as.character(sigma), "\\.")[[1]][2])
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/10,x+1))
          axis(1, seq, las=1, cex.axis=1.3)
        }   
      }
    }
    data <- data.frame(Lower.Limit=ci, Upper.Limit=Inf, Z.critic=crit.v, Z.value=z, p.value=pvalue)
    return(data)
  }
  Sever.t.p <- function(n, dif, sigma, sigma.a, sigma.b, cor, alpha, plot=TRUE, axis=FALSE) {
    if(input$select=='cor'| input$select2=='cor2') sigma <- sqrt(sigma.a^2+sigma.b^2-2*cor*sigma.a*sigma.b)
    validate(
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample size cannot be lower than 2"))
    se <- sqrt(sigma^2/n)
    t <- dif/se 
    pvalue <- 1 - pt(t, df = n-1)
    crit.v <- qt(1-alpha, df = n-1)
    me <- crit.v*se
    ci <- dif - me
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(pvalue <= alpha) {
        find <- function(discrepancy, Severity) (pt(t, n-1, discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (pt(t, n-1, -discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- -suppressWarnings(optimize(f=find2, Severity=0.999, interval=interval, maximum=F, tol = 1e-9))[[1]]
        Discrepancy2 <- -suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
        symbol <- bquote(.("SEV(") ~ mu > gamma~.(")"))
        curve(pt(t, n-1, x/se), xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)
      }
      else if(round(pvalue, 5)==1 & t<0) {
        symbol <- bquote(.("SEV(") ~ mu <= gamma~.(")"))
        find2 <- function(discrepancy, Severity) (1-pt(t, n-1, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- suppressWarnings(optimize(f=find2, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
        Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        curve(1 - pt(t, n-1, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)     
      }
      else{
        symbol <- bquote(.("SEV(") ~ mu <= gamma~.(")"))
        find <- function(discrepancy, Severity) (1-pt(t, n-1, -discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (1-pt(t, n-1, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        if(t>=0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        if(t<0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        } 
        curve(1 - pt(t, n-1, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA, cex.main=1.5)     
      }
      title(ylab="Severity / Power", cex.lab=1.6, line = 3.7)
      title(xlab=bquote(gamma == mu[1] - mu[0]), cex.lab=1.5, line = 3.5)
      Discrepancies=c(Discrepancy, Discrepancy2)
      Discrepancy=Discrepancies[which.min(Discrepancies)]; Discrepancy2=Discrepancies[which.max(Discrepancies)]
      curve(1 - pt(crit.v, n-1, x/se), add=T, col="#F36565", lwd=4, xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), xaxt='n', yaxt='n', xlab=NA, ylab=NA)
      text(x=Discrepancy2-(Discrepancy2-Discrepancy)*0.15, y=0.7, labels = bquote(hat(mu) == .(round(dif, 2))), cex=1.3)
      axis(2, seq(0, 1, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from2, input$to2, input$by2), las=1, cex.axis=1.3)
      }
      else {
        if(sigma>=1) {
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/8.5, 1))
          axis(1, seq, las=1, cex.axis=1.3)   
        }
        else {
          x <- nchar(strsplit(as.character(sigma), "\\.")[[1]][2])
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/10,x+1))
          axis(1, seq, las=1, cex.axis=1.3)
        }   
      }}
    data <- data.frame(Lower.Limit=ci, Upper.Limit=Inf, t.critic=crit.v, t.value=t, p.value=pvalue)
    return(data)
  }
  Sever.Z.two <- function(n1, n2, dif, sigma1, sigma2, alpha, plot=TRUE, axis=FALSE) {
    validate(
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma1>0 & sigma2>0, "both sigmas must be a greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n1>=2 & n2>=2, "the sample sizes cannot be lower than 2"))
    var1 <- sigma1^2
    var2 <- sigma2^2
    sigma <- sqrt(((n1-1)*var1 + (n2-1)*var2)/(n1+n2-2))
    se <- sqrt(var1/n1+var2/n2)
    z <- dif/se
    pvalue <- 1 - pnorm(z)
    crit.v <- qnorm(1-alpha)
    me <- crit.v*se
    ci <- dif - me
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(pvalue <= alpha) {
        find <- function(discrepancy, Severity) (pnorm(z, discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (pnorm(z, -discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- suppressWarnings(optimize(f=find, Severity=0.999, interval=interval, maximum=F, tol = 1e-9))[[1]]
        Discrepancy2 <- -suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval2, maximum=F, tol = 1e-9))[[1]]
        symbol <- bquote(.("SEV(") ~ mu > gamma~.(")"))
        curve(pnorm(z, x/se), xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)
      }
      else{
        find <- function(discrepancy, Severity) (1-pnorm(z, -discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (1-pnorm(z, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        if(z>=0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        if(z<0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }     
        symbol <- bquote(.("SEV(") ~ mu <= gamma~.(")"))
        curve(1 - pnorm(z, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)     
      }
      title(ylab="Severity / Power", cex.lab=1.6, line = 3.7)
      title(xlab=bquote(gamma == mu[1] - mu[0]), cex.lab=1.5, line = 3.5)
      Discrepancies=c(Discrepancy, Discrepancy2)
      Discrepancy=Discrepancies[which.min(Discrepancies)]; Discrepancy2=Discrepancies[which.max(Discrepancies)]
      curve(1 - pnorm(crit.v, x/se), add=T, col="#F36565", lwd=4, xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), xaxt='n', yaxt='n', xlab=NA, ylab=NA)
      text(x=Discrepancy2-(Discrepancy2-Discrepancy)*0.15, y=0.7, labels = bquote(hat(mu) == .(round(dif, 2))), cex=1.3)
      axis(2, seq(0, 1, 0.1), las=2, cex.axis=1.3)
      sigmas <- c(sigma1, sigma2)
      x1 <- nchar(strsplit(as.character(sigma1), "\\.")[[1]][2])
      x2 <- nchar(strsplit(as.character(sigma2), "\\.")[[1]][2])
      x <- c(x1, x2)
      dec <- sigmas[which.max(x)]
      if(axis) {
        axis(1, seq(input$from3, input$to3, input$by3), las=1, cex.axis=1.3)
      }
      else {
        if(sigma>=1) {
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/5, 1))
          axis(1, seq, las=1, cex.axis=1.3)   
        }
        else {
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/5, dec+1))
          axis(1, seq, las=1, cex.axis=1.3)
        }
      }
    }
    data <- data.frame(Lower.Limit=ci, Upper.Limit=Inf, Z.critic=crit.v, Z.value=z, p.value=pvalue)
    return(data)
  }
  Sever.t.two <- function(n1, n2, dif, sigma1, sigma2, alpha, plot=TRUE, axis=FALSE) {
    validate(
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma1>0 & sigma2>0, "both sigmas must be a greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n1>=2 & n2>=2, "the sample sizes cannot be lower than 2"))
    var1 <- sigma1^2
    var2 <- sigma2^2
    sigma <- sqrt(((n1-1)*var1 + (n2-1)*var2)/(n1+n2-2))
    se <- sqrt(var1/n1+var2/n2)
    t <- dif/se
    dfST <- ((var1/n1 + var2/n2)^2) / (var1^2/((n1-1)*n1^2) + var2^2/((n2-1)*n2^2))
    pvalue <- 1 - pt(t, df = dfST)
    crit.v <- qt(1-alpha, df = dfST)
    me <- crit.v*se
    ci <- dif - me
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(pvalue <= alpha) {
        find <- function(discrepancy, Severity) (pt(t, dfST, discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (pt(t, dfST, -discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- -suppressWarnings(optimize(f=find2, Severity=0.999, interval=interval, maximum=F, tol = 1e-9))[[1]]
        Discrepancy2 <- -suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
        symbol <- bquote(.("SEV(") ~ mu > gamma~.(")"))
        curve(pt(t, dfST, x/se), xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)
      }
      else if(round(pvalue, 5)==1 & t<0) {
        symbol <- bquote(.("SEV(") ~ mu > gamma~.(")"))
        find2 <- function(discrepancy, Severity) (1-pt(t, dfST, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        Discrepancy <- suppressWarnings(optimize(f=find2, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
        Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        curve(1 - pt(t, dfST, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)     
      }
      else{
        symbol <- bquote(.("SEV(") ~ mu <= gamma~.(")"))
        find <- function(discrepancy, Severity) (1-pt(t, dfST, -discrepancy/se)-Severity)^2
        find2 <- function(discrepancy, Severity) (1-pt(t, dfST, discrepancy/se)-Severity)^2
        interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
        if(t>=0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        if(t<0) {
          Discrepancy <- -suppressWarnings(optimize(f=find, Severity=0.999, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          Discrepancy2 <- suppressWarnings(optimize(f=find2, Severity=0.001, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }     
        curve(1 - pt(t, dfST, x/se), xlim=c(Discrepancy2, Discrepancy), ylim=c(0,1), type="l", col="#6582F3", lwd=4, main=symbol, cex.main=1.5, xaxt='n', yaxt='n', ylab=NA, xlab=NA)     
      }
      title(ylab="Severity / Power", cex.lab=1.6, line = 3.7)
      title(xlab=bquote(gamma == mu[1] - mu[0]), cex.lab=1.5, line = 3.5)
      Discrepancies=c(Discrepancy, Discrepancy2)
      Discrepancy=Discrepancies[which.min(Discrepancies)]; Discrepancy2=Discrepancies[which.max(Discrepancies)]
      curve(1 - pt(crit.v, dfST, x/se), add=T, col="#F36565", lwd=4, xlim=c(Discrepancy, Discrepancy2), ylim=c(0,1), xaxt='n', yaxt='n', xlab=NA, ylab=NA)
      text(x=Discrepancy2-(Discrepancy2-Discrepancy)*0.15, y=0.7, labels = bquote(hat(mu) == .(round(dif, 2))), cex=1.3)
      axis(2, seq(0, 1, 0.1), las=2, cex.axis=1.3)
      sigmas <- c(sigma1, sigma2)
      x1 <- nchar(strsplit(as.character(sigma1), "\\.")[[1]][2])
      x2 <- nchar(strsplit(as.character(sigma2), "\\.")[[1]][2])
      x <- c(x1, x2)
      dec <- sigmas[which.max(x)]
      if(axis) {
        axis(1, seq(input$from3, input$to3, input$by3), las=1, cex.axis=1.3)
      }
      else {
        if(sigma>=1) {
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/5, 1))
          axis(1, seq, las=1, cex.axis=1.3)   
        }
        else {
          seq <- seq(floor(Discrepancy), ceiling(Discrepancy2), round(sigma/5, dec+1))
          axis(1, seq, las=1, cex.axis=1.3)
        }   
      }
    }
    data <- data.frame(Lower.Limit=ci, Upper.Limit=Inf, t.critic=crit.v, t.value=t, p.value=pvalue)
    return(data)
  }
  output$plot1 <- renderPlot({
    if(input$dist=="Z") {
      if(input$adjust) {
        Sever.Z.one(n=input$obs, dif=input$dif, sigma=input$sd, alpha=input$alpha, plot=TRUE, axis=TRUE)
      }
      else {     
        Sever.Z.one(n=input$obs, dif=input$dif, sigma=input$sd, alpha=input$alpha, plot=TRUE, axis=FALSE)
      }
    }
    else{
      if(input$adjust) {
        Sever.t.one(n=input$obs, dif=input$dif, sigma=input$sd, alpha=input$alpha, plot=TRUE, axis=TRUE)
      }
      else {
        Sever.t.one(n=input$obs, dif=input$dif, sigma=input$sd, alpha=input$alpha, plot=TRUE, axis=FALSE)
      }
    }
  })
  output$plot2 <- renderPlot({
    if(input$dist2=="Z") {
      if(input$adjust2) {
        Sever.Z.p(n=input$obs2, dif=input$dif2, sigma=input$sd2, sigma.a=input$sda, sigma.b=input$sdb, cor=input$correlation, alpha=input$alpha2, plot=TRUE, axis=TRUE)
      }
      else {
        Sever.Z.p(n=input$obs2, dif=input$dif2, sigma=input$sd2, sigma.a=input$sda, sigma.b=input$sdb, cor=input$correlation, alpha=input$alpha2, plot=TRUE, axis=FALSE)
      }
    }
    else{
      if(input$adjust2) {
        Sever.t.p(n=input$obs2, dif=input$dif2, sigma=input$sd2, sigma.a=input$sda, sigma.b=input$sdb, cor=input$correlation, alpha=input$alpha2, plot=TRUE, axis=TRUE)
      }
      else{
        Sever.t.p(n=input$obs2, dif=input$dif2, sigma=input$sd2, sigma.a=input$sda, sigma.b=input$sdb, cor=input$correlation, alpha=input$alpha2, plot=TRUE, axis=FALSE)
      }
    }
  })
  output$plot3 <- renderPlot({
    if(input$dist3=="Z") {
      if(input$adjust3) {
        Sever.Z.two(n1=input$obs3, n2=input$obs4, dif=input$dif3, sigma1=input$sd3, sigma2=input$sd4, alpha=input$alpha3, plot=TRUE, axis=TRUE)
      }
      else{
        Sever.Z.two(n1=input$obs3, n2=input$obs4, dif=input$dif3, sigma1=input$sd3, sigma2=input$sd4, alpha=input$alpha3, plot=TRUE, axis=FALSE)
      }
    }
    else{
      if(input$adjust3) {
        Sever.t.two(n1=input$obs3, n2=input$obs4, dif=input$dif3, sigma1=input$sd3, sigma2=input$sd4, alpha=input$alpha3, plot=TRUE, axis=TRUE)
      }
      else{
        Sever.t.two(n1=input$obs3, n2=input$obs4, dif=input$dif3, sigma1=input$sd3, sigma2=input$sd4, alpha=input$alpha3, plot=TRUE, axis=FALSE)
      }}
  })
  observe({
    if(input$dist=="Z") {
      mytable <- reactive({
        Sever.Z.one(n=input$obs, dif=input$dif, sigma=input$sd, alpha=input$alpha, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.Z.one(n=input$obs, dif=input$dif, sigma=input$sd, alpha=input$alpha, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat <- renderTable(mytable(), digits=4)
      output$dat1 <- renderTable(mytable2(), digits=4)
    }
    else{
      mytable <- reactive({
        Sever.t.one(n=input$obs, dif=input$dif, sigma=input$sd, alpha=input$alpha, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.t.one(n=input$obs, dif=input$dif, sigma=input$sd, alpha=input$alpha, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat <- renderTable(mytable(), digits=4)
      output$dat1 <- renderTable(mytable2(), digits=4)
    }
  })
  observe({
    if(input$dist2=="Z") {
      mytable <- reactive({
        Sever.Z.p(n=input$obs2, dif=input$dif2, sigma=input$sd2, sigma.a=input$sda, sigma.b=input$sdb, cor=input$correlation, alpha=input$alpha2, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.Z.p(n=input$obs2, dif=input$dif2, sigma=input$sd2, sigma.a=input$sda, sigma.b=input$sdb, cor=input$correlation, alpha=input$alpha2, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat2 <- renderTable(mytable(), digits=4)
      output$dat21 <- renderTable(mytable2(), digits=4)
    }
    else{
      mytable <- reactive({
        Sever.t.p(n=input$obs2, dif=input$dif2, sigma=input$sd2, sigma.a=input$sda, sigma.b=input$sdb, cor=input$correlation, alpha=input$alpha2, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.t.p(n=input$obs2, dif=input$dif2, sigma=input$sd2, sigma.a=input$sda, sigma.b=input$sdb, cor=input$correlation, alpha=input$alpha2, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat2 <- renderTable(mytable(), digits=4)
      output$dat21 <- renderTable(mytable2(), digits=4)
    }
  })
  observe({
    if(input$dist3=="Z") {
      mytable <- reactive({
        Sever.Z.two(n1=input$obs3, n2=input$obs4, dif=input$dif3, sigma1=input$sd3, sigma2=input$sd4, alpha=input$alpha3, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.Z.two(n1=input$obs3, n2=input$obs4, dif=input$dif3, sigma1=input$sd3, sigma2=input$sd4, alpha=input$alpha3, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat3 <- renderTable(mytable(), digits=4)
      output$dat31 <- renderTable(mytable2(), digits=4)
    }
    else{
      mytable <- reactive({
        Sever.t.two(n1=input$obs3, n2=input$obs4, dif=input$dif3, sigma1=input$sd3, sigma2=input$sd4, alpha=input$alpha3, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.t.two(n1=input$obs3, n2=input$obs4, dif=input$dif3, sigma1=input$sd3, sigma2=input$sd4, alpha=input$alpha3, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat3 <- renderTable(mytable(), digits=4)
      output$dat31 <- renderTable(mytable2(), digits=4)
    }
  })
  observe({
    dif <- input$dif
    alpha <- input$alpha
    sigma <- input$sd
    n <- input$obs
    severity <- input$severity
    se <- sigma/sqrt(n)
    validate(
      need(exists("severity", mode="integer"), "Set a severity value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample sizes cannot be lower than 2"))
    observeEvent(c(input$dist, severity, dif, n, se, alpha), {
      if(input$dist=="Z") {
        z <- dif/se
        crit.v <- qnorm(1-alpha)
        pvalue <- 1 - pnorm(z)
        if(pvalue <= alpha) {
          if(severity>=0.5) {
            find <- function(discrepancy, Severity) (pnorm(z, discrepancy/se)-Severity)^2
            interval <- c(-15000, dif)
            Discrepancy <- suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
          if(severity<0.5) {
            find <- function(discrepancy, Severity) (pnorm(z, -discrepancy/se)-Severity)^2
            interval <- c(-15000, -dif)
            Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
        }
        else{
          find <- function(discrepancy, Severity) (1-pnorm(z, -discrepancy/se)-Severity)^2
          find2 <- function(discrepancy, Severity) (1-pnorm(z, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(z>=0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
          if(z<0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
        }
        power <- 1 - pnorm(crit.v, Discrepancy/se)
        updateNumericInput(session, inputId = "discrepancy", value = round(Discrepancy, 6))
        updateNumericInput(session, inputId = "power", value = round(power, 6))
      }
      if(input$dist=="t")  {
        t <- dif/se
        crit.v <- qt(1-alpha, n-1)
        pvalue <- 1 - pt(t, n-1)
        if(pvalue <= alpha) {
          if(severity>=0.5) {
            find <- function(discrepancy, Severity) (pt(t, n-1, -discrepancy/se)-Severity)^2
            interval <- c(-15000, dif)
            Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
          if(severity<0.5) {
            find <- function(discrepancy, Severity) (pt(t, n-1, -discrepancy/se)-Severity)^2
            interval <- c(-15000, -dif)
            Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
        }
        else if(round(pvalue, 5)==1 & t<0) {
          find2 <- function(discrepancy, Severity) (1-pt(t, n-1, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(severity>=0.5) Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          if(severity<0.5) Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        else{
          find <- function(discrepancy, Severity) (1-pt(t, n-1, -discrepancy/se)-Severity)^2
          find2 <- function(discrepancy, Severity) (1-pt(t, n-1, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(t>=0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
          if(t<0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
        }
        power <- 1 - pt(crit.v, n-1, Discrepancy/se)
        updateNumericInput(session, inputId = "discrepancy", value = round(Discrepancy, 6))
        updateNumericInput(session, inputId = "power", value = round(power, 6))
      }})
  })
  observe({
    dif <- input$dif2
    alpha <- input$alpha2
    if(input$select=='sigmascorediff') sigma <- input$sd2 
    if(input$select=='cor') sigma <- sqrt(input$sda^2+input$sdb^2-2*input$correlation*input$sda*input$sdb)
    n <- input$obs2
    severity <- input$severity2
    se <- sqrt(sigma^2/(n))
    validate(
      need(exists("severity", mode="integer"), "Set a severity value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample size cannot be lower than 2"))
    observeEvent(c(input$dist2, severity, dif, n, se, alpha), {
      if(input$dist2=="Z") {
        z <- dif/se
        crit.v <- qnorm(1-alpha)
        pvalue <- 1 - pnorm(z)
        if(pvalue <= alpha) {
          if(severity>=0.5) {
            find <- function(discrepancy, Severity) (pnorm(z, discrepancy/se)-Severity)^2
            interval <- c(-15000, dif)
            Discrepancy <- suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
          if(severity<0.5) {
            find <- function(discrepancy, Severity) (pnorm(z, -discrepancy/se)-Severity)^2
            interval <- c(-15000, -dif)
            Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
        }
        else{
          find <- function(discrepancy, Severity) (1-pnorm(z, -discrepancy/se)-Severity)^2
          find2 <- function(discrepancy, Severity) (1-pnorm(z, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(z>=0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
          if(z<0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
        }
        power <- 1 - pnorm(crit.v, Discrepancy/se)
        updateNumericInput(session, inputId = "discrepancy2", value = round(Discrepancy, 6))
        updateNumericInput(session, inputId = "power2", value = round(power, 6))
      }
      if(input$dist2=="t")  {
        t <- dif/se
        crit.v <- qt(1-alpha, n-1)
        pvalue <- 1 - pt(t, n-1)
        if(pvalue <= alpha) {
          if(severity>=0.5) {
            find <- function(discrepancy, Severity) (pt(t, n-1, discrepancy/se)-Severity)^2
            interval <- c(-15000, dif)
            Discrepancy <- suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
          if(severity<0.5) {
            find <- function(discrepancy, Severity) (pt(t, n-1, -discrepancy/se)-Severity)^2
            interval <- c(-15000, -dif)
            Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
        }
        else if(round(pvalue, 5)==1 & t<0) {
          find2 <- function(discrepancy, Severity) (1-pt(t, n-1, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(severity>=0.5) Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          if(severity<0.5) Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        else{
          find <- function(discrepancy, Severity) (1-pt(t, n-1, -discrepancy/se)-Severity)^2
          find2 <- function(discrepancy, Severity) (1-pt(t, n-1, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(t>=0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
          if(t<0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
        }
        power <- 1 - pt(crit.v, n-1, Discrepancy/se)
        updateNumericInput(session, inputId = "discrepancy2", value = round(Discrepancy, 6))
        updateNumericInput(session, inputId = "power2", value = round(power, 6))
      }
    })   
  })
  observe({
    dif <- input$dif3
    alpha <- input$alpha3
    sigma1 <- input$sd3
    sigma2 <- input$sd4
    n1 <- input$obs3
    n2 <- input$obs4
    severity <- input$severity3
    var1 <- sigma1^2
    var2 <- sigma2^2
    se <- sqrt(var1/n1+var2/n2)
    dfST <- ((var1/n1 + var2/n2)^2) / (var1^2/((n1-1)*n1^2) + var2^2/((n2-1)*n2^2))
    validate(
      need(exists("severity", mode="integer"), "Set a severity value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma1>0 & sigma2>0, "both sigmas must be a greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n1>=2 & n2>=2, "the sample sizes cannot be lower than 2"))
    observeEvent(c(input$dist3, severity, dif, se, dfST, alpha), {
      if(input$dist3=="Z") {
        z <- dif/se
        pvalue <- 1 - pnorm(z)
        crit.v <- qnorm(1-alpha)
        if(pvalue <= alpha) {
          if(severity>=0.5) {
            find <- function(discrepancy, Severity) (pnorm(z, discrepancy/se)-Severity)^2
            interval <- c(-15000, dif)
            Discrepancy <- suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
          if(severity<0.5) {
            find <- function(discrepancy, Severity) (pnorm(z, -discrepancy/se)-Severity)^2
            interval <- c(-15000, -dif)
            Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
        }
        else{
          find <- function(discrepancy, Severity) (1-pnorm(z, -discrepancy/se)-Severity)^2
          find2 <- function(discrepancy, Severity) (1-pnorm(z, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(z>=0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
          if(z<0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
        }
        power <- 1 - pnorm(crit.v, Discrepancy/se)
        updateNumericInput(session, inputId = "discrepancy3", value = round(Discrepancy, 6))
        updateNumericInput(session, inputId = "power3", value = round(power, 6))
      }
      if(input$dist3=="t") {
        t <- dif/se
        pvalue <- 1 - pt(t, dfST)
        crit.v <- qt(1-alpha, dfST)
        if(pvalue <= alpha) {
          if(severity>=0.5) {
            find <- function(discrepancy, Severity) (pt(t, dfST, discrepancy/se)-Severity)^2
            interval <- c(-15000, dif)
            Discrepancy <- suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
          if(severity<0.5) {
            find <- function(discrepancy, Severity) (pt(t, dfST, -discrepancy/se)-Severity)^2
            interval <- c(-15000, -dif)
            Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=F, tol = 1e-9))[[1]]
          }
        }
        else if(round(pvalue, 5)==1 & t<0) {
          find2 <- function(discrepancy, Severity) (1-pt(t, dfST, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(severity>=0.5) Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
          if(severity<0.5) Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
        }
        else{
          find <- function(discrepancy, Severity) (1-pt(t, dfST, -discrepancy/se)-Severity)^2
          find2 <- function(discrepancy, Severity) (1-pt(t, dfST, discrepancy/se)-Severity)^2
          interval <- c(-15000, dif); interval2 <- c(-15000, -dif)
          if(t>=0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
          if(t<0) {
            if(severity>=0.5) Discrepancy <- -suppressWarnings(optimize(f=find, Severity=severity, interval=interval2, maximum=FALSE, tol = 1e-9))[[1]]
            else { Discrepancy <- suppressWarnings(optimize(f=find2, Severity=severity, interval=interval, maximum=FALSE, tol = 1e-9))[[1]] }
          }
        }
        power <- 1 - pt(crit.v, dfST, Discrepancy/se)
        updateNumericInput(session, inputId = "discrepancy3", value = round(Discrepancy, 6))
        updateNumericInput(session, inputId = "power3", value = round(power, 6))
      }
    })
  })
  output$test <- renderText("One-sided test:")
  output$test2 <- renderText("One-sided test:")
  output$test3 <- renderText("One-sided test:")
  Sever.Z.one.plots <- function(n, dif, sigma, alpha, discrepancy, plot=TRUE, type, axis=FALSE) {
    validate(
      need(exists("discrepancy", mode="integer"), "Set a discrepancy value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample sizes cannot be lower than 2"))
    se <- sigma/sqrt(n)
    z <- dif/se
    pvalue <- 1 - pnorm(z)
    crit.v <- qnorm(1-alpha)
    delta <- discrepancy/se
    if(pvalue<=alpha) Severity <- pnorm(z, delta)
    else {Severity <- 1 - pnorm(z, delta)}
    power <- 1 - pnorm(crit.v, delta)
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(delta>0) xlim <- c(qnorm(0.0001), qnorm(0.9999, delta))
      else { xlim <- c(qnorm(0.0001, delta), qnorm(0.9999)) }
      mynull <- curve(dnorm(x), n=1e3, xlim=xlim, ylim=c(0,0.4), main="Z statistic distributions", font.main = 1, cex.main=1.5, xlab=NA, ylab=NA, type='n', yaxt='n', xaxt='n')
      mycurve <- curve(dnorm(x, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), type='n', add=T, yaxt='n', xaxt='n')
      if(type=="severity") {
        if(pvalue <= alpha) {
          x1 <- min(which(mycurve$x <= z))
          x2 <- max(which(mycurve$x < z))
        }
        if(pvalue > alpha) {
          x1 <- max(which(mycurve$x >= z))
          x2 <- min(which(mycurve$x > z))
        }
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFFF9B", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="power") {
        q <- qnorm(1-power, delta)
        x1 <- min(which(mycurve$x >= q))
        x2 <- max(which(mycurve$x >= q))
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFC9C9", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="pvalue") {
        x1 <- min(which(mycurve$x >= z))
        x2 <- max(which(mycurve$x >= z))
        with(mynull, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#B7D7FF", border=NA), yaxt='n', xaxt='n')
      }
      curve(dnorm(x), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#6F6BF6", add=T, lwd=4, yaxt='n', xaxt='n')
      curve(dnorm(x, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#F5134F", lwd=4, add=T, yaxt='n', xaxt='n')
      title(ylab="Density", cex.lab=1.6, line = 3.7)
      text(x=xlim[1]-(xlim[1]-xlim[2])*0.9, y=0.35, labels = bquote(Delta == .(round(delta, 2))), cex=1.3)
      axis(2, seq(0, 0.4, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from4, input$to4, input$by4), las=1, cex.axis=1.3)
      }
      else {
        seq <- seq(round(xlim[1]), round(xlim[2]), 1)
        axis(1, seq, las=1, cex.axis=1.3)
      }
    }
    list(Severity=Severity, Power=power)
  }
  Sever.t.one.plots <- function(n, dif, sigma, alpha, discrepancy, plot=TRUE, type, axis=FALSE) {
    validate(
      need(exists("discrepancy", mode="integer"), "Set a discrepancy value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample sizes cannot be lower than 2"))
    se <- sigma/sqrt(n)
    t <- dif/se
    pvalue <- 1 - pt(t, n-1)
    crit.v <- qt(1-alpha, n-1)
    delta <- discrepancy/se
    if(pvalue<=alpha) Severity <- pt(t, n-1, delta)
    else {Severity <- 1 - pt(t, n-1, delta)}
    power <- 1 - pt(crit.v, n-1, delta)
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(delta>0) xlim <- c(qt(0.0001, n-1), qt(0.9999, n-1, delta))
      else { xlim <- c(qt(0.0001, n-1, delta), qt(0.9999, n-1)) }
      mynull <- curve(dt(x, n-1), n=1e3, xlim=xlim, ylim=c(0,0.4), main="t statistic distributions", font.main = 1, cex.main=1.5, xlab=NA,  ylab=NA, type='n', yaxt='n', xaxt='n')
      mycurve <- curve(dt(x, n-1, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), type='n', add=T, yaxt='n', xaxt='n')
      if(type=="severity") {
        if(pvalue <= alpha) {
          x1 <- min(which(mycurve$x <= t))
          x2 <- max(which(mycurve$x < t))
        }
        if(pvalue > alpha) {
          x1 <- max(which(mycurve$x >= t))
          x2 <- min(which(mycurve$x > t))
        }
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFFF9B", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="power") {
        q <- qt(1-power, n-1, delta)
        x1 <- min(which(mycurve$x >= q))
        x2 <- max(which(mycurve$x >= q))
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFC9C9", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="pvalue") {
        x1 <- min(which(mycurve$x >= t))
        x2 <- max(which(mycurve$x >= t))
        with(mynull, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#B7D7FF", border=NA, yaxt='n', xaxt='n'))
      }
      curve(dt(x, n-1), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#6F6BF6", add=T, lwd=4, yaxt='n', xaxt='n')
      curve(dt(x, n-1, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#F5134F", lwd=4, add=T, yaxt='n', xaxt='n')
      title(ylab="Density", cex.lab=1.6, line = 3.7)
      text(x=xlim[1]-(xlim[1]-xlim[2])*0.9, y=0.35, labels = bquote(Delta == .(round(delta, 2))), cex=1.3)
      axis(2, seq(0, 0.4, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from4, input$to4, input$by4), las=1, cex.axis=1.3)
      }
      else {
        seq <- seq(round(xlim[1]), round(xlim[2]), 1)
        axis(1, seq, las=1, cex.axis=1.3)
      }
    }
    list(Severity=Severity, Power=power)
  }
  Sever.Z.p.plots <- function(n, dif, sigma, sigma.a, sigma.b, cor, alpha, discrepancy, plot=TRUE, type, axis=FALSE) {
    if(input$select2=='sigmascorediff2') sigma <- input$sd20 
    if(input$select2=='cor2') sigma <- sqrt(sigma.a^2+sigma.b^2-2*cor*sigma.a*sigma.b)
    validate(
      need(exists("discrepancy", mode="integer"), "Set a discrepancy value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample sizes cannot be lower than 2"))
    se <- sqrt(sigma^2/n)
    z <- dif/se
    pvalue <- 1 - pnorm(z)
    crit.v <- qnorm(1-alpha)
    delta <- discrepancy/se
    if(pvalue<=alpha) Severity <- pnorm(z, delta)
    else {Severity <- 1 - pnorm(z, delta)}
    power <- 1 - pnorm(crit.v, delta)
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(delta>0) xlim <- c(qnorm(0.0001), qnorm(0.9999, delta))
      else { xlim <- c(qnorm(0.0001, delta), qnorm(0.9999)) }
      mynull <- curve(dnorm(x), n=1e3, xlim=xlim, ylim=c(0,0.4), main="Z statistic distributions", font.main = 1, cex.main=1.5, xlab=NA, ylab=NA, cex.lab=1.5, type='n', yaxt='n', xaxt='n')
      mycurve <- curve(dnorm(x, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), type='n', add=T, yaxt='n', xaxt='n')
      if(type=="severity") {
        if(pvalue <= alpha) {
          x1 <- min(which(mycurve$x <= z))
          x2 <- max(which(mycurve$x < z))
        }
        if(pvalue > alpha) {
          x1 <- max(which(mycurve$x >= z))
          x2 <- min(which(mycurve$x > z))
        }
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFFF9B", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="power") {
        q <- qnorm(1-power, delta)
        x1 <- min(which(mycurve$x >= q))
        x2 <- max(which(mycurve$x >= q))
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFC9C9", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="pvalue") {
        x1 <- min(which(mycurve$x >= z))
        x2 <- max(which(mycurve$x >= z))
        with(mynull, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#B7D7FF", border=NA, yaxt='n', xaxt='n'))
      }
      curve(dnorm(x), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#6F6BF6", add=T, lwd=4, yaxt='n', xaxt='n')
      curve(dnorm(x, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#F5134F", lwd=4, add=T, yaxt='n', xaxt='n')
      title(ylab="Density", cex.lab=1.6, line = 3.7)
      text(x=xlim[1]-(xlim[1]-xlim[2])*0.9, y=0.35, labels = bquote(Delta == .(round(delta, 2))), cex=1.3)
      axis(2, seq(0, 0.4, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from5, input$to5, input$by5), las=1, cex.axis=1.3)
      }
      else {
        seq <- seq(round(xlim[1]), round(xlim[2]), 1)
        axis(1, seq, las=1, cex.axis=1.3)
      }
    }
    list(Severity=Severity, Power=power)
  }
  Sever.t.p.plots <- function(n, dif, sigma, sigma.a, sigma.b, cor, alpha, discrepancy, plot=TRUE, type, axis=FALSE) {
    if(input$select2=='sigmascorediff2') sigma <- input$sd20 
    if(input$select2=='cor2') sigma <- sqrt(sigma.a^2+sigma.b^2-2*cor*sigma.a*sigma.b)
    validate(
      need(exists("discrepancy", mode="integer"), "Set a discrepancy value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma>0, "sigma must be greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n>=2, "the sample sizes cannot be lower than 2"))
    se <- sqrt(sigma^2/n)
    t <- dif/se
    pvalue <- 1 - pt(t, n-1)
    crit.v <- qt(1-alpha, n-1)
    delta <- discrepancy/se
    if(pvalue<=alpha) Severity <- pt(t, n-1, delta)
    else {Severity <- 1 - pt(t, n-1, delta)}
    power <- 1 - pt(crit.v, n-1, delta)
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(delta>0) xlim <- c(qt(0.0001, n-1), qt(0.9999, n-1, delta))
      else { xlim <- c(qt(0.0001, n-1, delta), qt(0.9999, n-1)) }
      mynull <- curve(dt(x, n-1), n=1e3, xlim=xlim, ylim=c(0,0.4), main="t statistic distributions", font.main = 1, cex.main=1.5, xlab=NA, ylab=NA, cex.lab=1.5, type='n', yaxt='n', xaxt='n')
      mycurve <- curve(dt(x, n-1, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), type='n', add=T, yaxt='n', xaxt='n')
      if(type=="severity") {
        if(pvalue <= alpha) {
          x1 <- min(which(mycurve$x <= t))
          x2 <- max(which(mycurve$x < t))
        }
        if(pvalue > alpha) {
          x1 <- max(which(mycurve$x >= t))
          x2 <- min(which(mycurve$x > t))
        }
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFFF9B", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="power") {
        q <- qt(1-power, n-1, delta)
        x1 <- min(which(mycurve$x >= q))
        x2 <- max(which(mycurve$x >= q))
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFC9C9", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="pvalue") {
        x1 <- min(which(mycurve$x >= t))
        x2 <- max(which(mycurve$x >= t))
        with(mynull, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#B7D7FF", border=NA, yaxt='n', xaxt='n'))
      }
      curve(dt(x, n-1), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#6F6BF6", add=T, lwd=4, yaxt='n', xaxt='n')
      curve(dt(x, n-1, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#F5134F", lwd=4, add=T, yaxt='n', xaxt='n')
      title(ylab="Density", cex.lab=1.6, line = 3.7)
      text(x=xlim[1]-(xlim[1]-xlim[2])*0.9, y=0.35, labels = bquote(Delta == .(round(delta, 2))), cex=1.3)
      axis(2, seq(0, 0.4, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from5, input$to5, input$by5), las=1, cex.axis=1.3)
      }
      else {
        seq <- seq(round(xlim[1]), round(xlim[2]), 1)
        axis(1, seq, las=1, cex.axis=1.3)
      }
    }
    list(Severity=Severity, Power=power)
  }
  Sever.Z.two.plots <- function(n1, n2, dif, sigma1, sigma2, alpha, discrepancy, plot=TRUE, type, axis=FALSE) {
    validate(
      need(exists("discrepancy", mode="integer"), "Set a discrepancy value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma1>0 & sigma2>0, "both sigmas must be a greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n1>=2 & n2>=2, "the sample sizes cannot be lower than 2"))
    var1 <- sigma1^2
    var2 <- sigma2^2
    sigma <- sqrt(((n1-1)*var1 + (n2-1)*var2)/(n1+n2-2))
    se <- sqrt(var1/n1+var2/n2)
    z <- dif/se
    pvalue <- 1 - pnorm(z)
    crit.v <- qnorm(1-alpha)
    delta <- discrepancy/se
    if(pvalue<=alpha) Severity <- pnorm(z, delta)
    else {Severity <- 1 - pnorm(z, delta)}
    power <- 1 - pnorm(crit.v, delta)
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(delta>0) xlim <- c(qnorm(0.0001), qnorm(0.9999, delta))
      else { xlim <- c(qnorm(0.0001, delta), qnorm(0.9999)) }
      mynull <- curve(dnorm(x), n=1e3, xlim=xlim, ylim=c(0,0.4), main="Z statistic distributions", font.main = 1, cex.main=1.5, xlab=NA, ylab=NA, cex.lab=1.5, type='n', yaxt='n', xaxt='n')
      mycurve <- curve(dnorm(x, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), type='n', add=T, yaxt='n', xaxt='n')
      if(type=="severity") {
        if(pvalue <= alpha) {
          x1 <- min(which(mycurve$x <= z))
          x2 <- max(which(mycurve$x < z))
        }
        if(pvalue > alpha) {
          x1 <- max(which(mycurve$x >= z))
          x2 <- min(which(mycurve$x > z))
        }
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFFF9B", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="power") {
        q <- qnorm(1-power, delta)
        x1 <- min(which(mycurve$x >= q))
        x2 <- max(which(mycurve$x >= q))
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFC9C9", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="pvalue") {
        x1 <- min(which(mycurve$x >= z))
        x2 <- max(which(mycurve$x >= z))
        with(mynull, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#B7D7FF", border=NA, yaxt='n', xaxt='n'))
      }
      curve(dnorm(x), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#6F6BF6", add=T, lwd=4, yaxt='n', xaxt='n')
      curve(dnorm(x, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#F5134F", lwd=4, add=T, yaxt='n', xaxt='n')
      text(x=xlim[1]-(xlim[1]-xlim[2])*0.9, y=0.35, labels = bquote(Delta == .(round(delta, 2))), cex=1.3)
      title(ylab="Density", cex.lab=1.6, line = 3.7)
      axis(2, seq(0, 0.4, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from6, input$to6, input$by6), las=1, cex.axis=1.3)
      }
      else {
        seq <- seq(round(xlim[1]), round(xlim[2]), 1)
        axis(1, seq, las=1, cex.axis=1.3)
      }
    }
    list(Severity=Severity, Power=power)
  }
  Sever.t.two.plots <- function(n1, n2, dif, sigma1, sigma2, alpha, discrepancy, plot=TRUE, type, axis=FALSE) {
    validate(
      need(exists("discrepancy", mode="integer"), "Set a discrepancy value"),
      need(exists("dif", mode="integer"), "Set the difference between means"),
      need(sigma1>0 & sigma2>0, "both sigmas must be a greater than 0"),
      need(alpha>=0 & alpha<=0.5, "alpha must lie between 0 and 0.5"),
      need(n1>=2 & n2>=2, "the sample sizes cannot be lower than 2"))
    var1 <- sigma1^2
    var2 <- sigma2^2
    sigma <- sqrt(((n1-1)*var1 + (n2-1)*var2)/(n1+n2-2))
    se <- sqrt(var1/n1+var2/n2)
    t <- dif/se
    dfST <- ((var1/n1 + var2/n2)^2) / (var1^2/((n1-1)*n1^2) + var2^2/((n2-1)*n2^2))
    pvalue <- 1 - pt(t, dfST)
    crit.v <- qt(1-alpha, dfST)
    delta <- discrepancy/se
    if(pvalue<=alpha) Severity <- pt(t, dfST, delta)
    else {Severity <- 1 - pt(t, dfST, delta)}
    power <- 1 - pt(crit.v, dfST, delta)
    if(plot) {
      par(mar=c(5, 5, 4, 1) + 0.1)
      if(delta>0) xlim <- c(qt(0.0001, dfST), qt(0.9999, dfST, delta))
      else { xlim <- c(qt(0.0001, dfST, delta), qt(0.9999, dfST)) }
      mynull <- curve(dt(x, dfST), n=1e3, xlim=xlim, ylim=c(0,0.4), main="t statistic distributions", font.main = 1, cex.main=1.5, xlab=NA, ylab=NA, cex.lab=1.5, type='n', yaxt='n', xaxt='n')
      mycurve <- curve(dt(x, dfST, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), type='n', add=T, yaxt='n', xaxt='n')
      if(type=="severity") {
        if(pvalue <= alpha) {
          x1 <- min(which(mycurve$x <= t))
          x2 <- max(which(mycurve$x < t))
        }
        if(pvalue > alpha) {
          x1 <- max(which(mycurve$x >= t))
          x2 <- min(which(mycurve$x > t))
        }
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFFF9B", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="power") {
        q <- qt(1-power, dfST, delta)
        x1 <- min(which(mycurve$x >= q))
        x2 <- max(which(mycurve$x >= q))
        with(mycurve, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#FFC9C9", border=NA, yaxt='n', xaxt='n'))
      }
      if(type=="pvalue") {
        x1 <- min(which(mycurve$x >= t))
        x2 <- max(which(mycurve$x >= t))
        with(mynull, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="#B7D7FF", border=NA, yaxt='n', xaxt='n'))
      }
      curve(dt(x, dfST), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#6F6BF6", add=T, lwd=4, yaxt='n', xaxt='n')
      curve(dt(x, dfST, delta), n=1e3, xlim=xlim, ylim=c(0,0.4), col="#F5134F", lwd=4, add=T, yaxt='n', xaxt='n')
      title(ylab="Density", cex.lab=1.6, line = 3.7)
      text(x=xlim[1]-(xlim[1]-xlim[2])*0.9, y=0.35, labels = bquote(Delta == .(round(delta, 2))), cex=1.3)
      axis(2, seq(0, 0.4, 0.1), las=2, cex.axis=1.3)
      if(axis) {
        axis(1, seq(input$from6, input$to6, input$by6), las=1, cex.axis=1.3)
      }
      else {
        seq <- seq(round(xlim[1]), round(xlim[2]), 1)
        axis(1, seq, las=1, cex.axis=1.3)
      }
    }
    list(Severity=Severity, Power=power)
  }
  
  output$plot4 <- renderPlot({
    if(input$dist10=="Z") {
      if(input$adjust4) {
        if(input$display=="S") Sever.Z.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="severity", axis=T)
        if(input$display=="P") Sever.Z.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="power", axis=T)
        if(input$display=="P-val") Sever.Z.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="pvalue", axis=T)
      }
      else {
        if(input$display=="S") Sever.Z.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="severity")
        if(input$display=="P") Sever.Z.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="power")
        if(input$display=="P-val") Sever.Z.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="pvalue")
      }
    }
    else{
      if(input$adjust4) {
        if(input$display=="S") Sever.t.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="severity", axis=T)
        if(input$display=="P") Sever.t.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="power", axis=T)
        if(input$display=="P-val") Sever.t.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="pvalue", axis=T)
      }
      else {
        if(input$display=="S") Sever.t.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="severity")
        if(input$display=="P") Sever.t.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="power")
        if(input$display=="P-val") Sever.t.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=TRUE, type="pvalue")
      }}})
  output$plot5 <- renderPlot({
    if(input$dist20=="Z") {
      if(input$adjust5) {
        if(input$display2=="S") Sever.Z.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="severity", axis=T)
        if(input$display2=="P") Sever.Z.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="power", axis=T)
        if(input$display2=="P-val") Sever.Z.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="pvalue", axis=T)
      }
      else {
        if(input$display2=="S") Sever.Z.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="severity")
        if(input$display2=="P") Sever.Z.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="power")
        if(input$display2=="P-val") Sever.Z.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="pvalue")
      }}
    else{
      if(input$adjust5) {
        if(input$display2=="S") Sever.t.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="severity", axis=T)
        if(input$display2=="P") Sever.t.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="power", axis=T)
        if(input$display2=="P-val") Sever.t.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="pvalue", axis=T)
      }
      else {
        if(input$display2=="S") Sever.t.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="severity")
        if(input$display2=="P") Sever.t.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="power")
        if(input$display2=="P-val") Sever.t.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=TRUE, type="pvalue")
      }}})
  output$plot6 <- renderPlot({
    if(input$dist30=="Z") {
      if(input$adjust6) {
        if(input$display3=="S") Sever.Z.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="severity", axis=T)
        if(input$display3=="P") Sever.Z.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="power", axis=T)
        if(input$display3=="P-val") Sever.Z.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="pvalue", axis=T)
      }
      else {
        if(input$display3=="S") Sever.Z.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="severity")
        if(input$display3=="P") Sever.Z.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="power")
        if(input$display3=="P-val") Sever.Z.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="pvalue")
      }}
    else{
      if(input$adjust6) {
        if(input$display3=="S") Sever.t.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="severity", axis=T)
        if(input$display3=="P") Sever.t.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="power", axis=T)
        if(input$display3=="P-val") Sever.t.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="pvalue", axis=T)
      }
      else {
        if(input$display3=="S") Sever.t.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="severity")
        if(input$display3=="P") Sever.t.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="power")
        if(input$display3=="P-val") Sever.t.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=TRUE, type="pvalue")
      }}})
  output$test10 <- renderText("One-sided test:")
  output$test20 <- renderText("One-sided test:")
  output$test30 <- renderText("One-sided test:")
  observe({
    if(input$dist10=="Z") {
      mytable <- reactive({
        Sever.Z.one(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.Z.one(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat10 <- renderTable(mytable(), digits=4)
      output$dat100 <- renderTable(mytable2(), digits=4)
    }
    else{
      mytable <- reactive({
        Sever.t.one(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.t.one(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat10 <- renderTable(mytable(), digits=4)
      output$dat100 <- renderTable(mytable2(), digits=4)
    }
  })
  observe({
    if(input$dist20=="Z") {
      mytable <- reactive({
        Sever.Z.p(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.Z.p(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat20 <- renderTable(mytable(), digits=4)
      output$dat200 <- renderTable(mytable2(), digits=4)
    }
    else{
      mytable <- reactive({
        Sever.t.p(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.t.p(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat20 <- renderTable(mytable(), digits=4)
      output$dat200 <- renderTable(mytable2(), digits=4)
    }
  })
  observe({
    if(input$dist30=="Z") {
      mytable <- reactive({
        Sever.Z.two(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.Z.two(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat30 <- renderTable(mytable(), digits=4)
      output$dat300 <- renderTable(mytable2(), digits=4)
    }
    else{
      mytable <- reactive({
        Sever.t.two(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, plot=FALSE, axis=FALSE)[1:2]
      })
      mytable2 <- reactive({
        Sever.t.two(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, plot=FALSE, axis=FALSE)[3:5]
      })
      output$dat30 <- renderTable(mytable(), digits=4)
      output$dat300 <- renderTable(mytable2(), digits=4)
    }
  })
  observe({
    if(input$dist10=="Z") {
      Severity <- Sever.Z.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=FALSE)$Severity
      Power <- Sever.Z.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=FALSE)$Power
      updateNumericInput(session, inputId = "severity10", value = round(Severity, 6))
      updateNumericInput(session, inputId = "power10", value = round(Power, 6))
    }
    else{
      Severity <- Sever.t.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=FALSE)$Severity
      Power <- Sever.t.one.plots(n=input$obs10, dif=input$dif10, sigma=input$sd10, alpha=input$alpha10, discrepancy=input$discrepancy10, plot=FALSE)$Power
      updateNumericInput(session, inputId = "severity10", value = round(Severity, 6))
      updateNumericInput(session, inputId = "power10", value = round(Power, 6))
    }
  })
  observe({
    if(input$dist20=="Z") {
      Severity <- Sever.Z.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=FALSE)$Severity
      Power <- Sever.Z.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=FALSE)$Power
      updateNumericInput(session, inputId = "severity20", value = round(Severity, 6))
      updateNumericInput(session, inputId = "power20", value = round(Power, 6))
    }
    else{
      Severity <- Sever.t.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=FALSE)$Severity
      Power <- Sever.t.p.plots(n=input$obs20, dif=input$dif20, sigma=input$sd20, sigma.a=input$sda2, sigma.b=input$sdb2, cor=input$correlation2, alpha=input$alpha20, discrepancy=input$discrepancy20, plot=FALSE)$Power
      updateNumericInput(session, inputId = "severity20", value = round(Severity, 6))
      updateNumericInput(session, inputId = "power20", value = round(Power, 6))
    }
  })
  observe({
    if(input$dist30=="Z") {
      Severity <- Sever.Z.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=FALSE)$Severity
      Power <- Sever.Z.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=FALSE)$Power
      updateNumericInput(session, inputId = "severity30", value = round(Severity, 6))
      updateNumericInput(session, inputId = "power30", value = round(Power, 6))
    }
    else{
      Severity <- Sever.t.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=FALSE)$Severity
      Power <- Sever.t.two.plots(n1=input$obs30, n2=input$obs40, dif=input$dif30, sigma1=input$sd30, sigma2=input$sd40, alpha=input$alpha30, discrepancy=input$discrepancy30, plot=FALSE)$Power
      updateNumericInput(session, inputId = "severity30", value = round(Severity, 6))
      updateNumericInput(session, inputId = "power30", value = round(Power, 6))
    }
  })
})