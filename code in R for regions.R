regions <- function(coord, xcoord, ycoord, col = col){
  t <- seq(0, 2*pi, length = 1000)
  pcoord1 <- coord[1] + xcoord*cos(t)
  pcoord2 <- coord[2] + ycoord*sin(t)
  lines(pcoord1, pcoord2, col = col)
  
}


ca.regions.exe <- function (N, a1 = 1, a2 = 2, alpha = 0.05,delta = 1, cols = c(2, 4), M = min(nrow(N), ncol(N)) - 1, region = 2, scaleplot = 1.2)
  
{
  
  
  I <- nrow(N) # Number of rows of table
  J <- ncol(N) # Number of columns of table
  M <- min(I, J)
  
  Inames <- dimnames(N)[1] # Row category names
  Jnames <- dimnames(N)[2] # Column category names
  
  if (delta == 0){
    for (i in 1:I){
      for (j in 1:J){
        if (N[i,j]==0){
          N[i,j] <- nonzero
        }
      }
    }
  }
  
  n <- sum(N) # Total number of classifications in the table
  p <- N * (1/n) # Matrix of joint relative proportions
  pidot <- apply(p, 1, sum)
  pdotj <- apply(p, 2, sum)
  
  Imass <- as.matrix(apply(p, 1, sum))
  Jmass <- as.matrix(apply(p, 2, sum))
  
  dI <- diag(c(Imass), nrow = I, ncol = I)
  dJ <- diag(c(Jmass), nrow = J, ncol = J)
  Ih <- Imass^-0.5
  Jh <- Jmass^-0.5
  dIh <- diag(c(Ih), nrow = I, ncol = I)
  dJh <- diag(c(Jh), nrow = J, ncol = J)
  
  r <- matrix(0, nrow = I, ncol = J)
  
  if (delta !=0){
    for (i in 1:I){
      for (j in 1:J){
        r[i,j] <-(sqrt(pidot[i]*pdotj[j])/delta)*((p[i,j]/(pidot[i]*pdotj[j]))^delta - 1) 
      }
    }
  } else {
    for (i in 1:I){
      for (j in 1:J){
        #			r[i,j] <- sqrt(pidot[i]*pdotj[j])*log((p[i,j]/(pidot[i]*pdotj[j])))
        r[i,j] <- -sqrt(pidot[i]*pdotj[j])*log(((pidot[i]*pdotj[j])/p[i,j]))
        
      }
    }
  }
  
  sva <- svd(r) # SVD of the Pearson residuals
  
  dmu <- diag(sva$d)
  
  ##########################################################################
  # #
  # Principal Coordinates #
  # #
  ##########################################################################
  
  f <- dIh %*% sva$u %*% dmu # Row Principal Coordinates for Classical CA
  g <- dJh %*% sva$v %*% dmu # Column Principal Coordinates for Classical CA
  
  #  f <- dIh %*% sva$u # Row Standard Coordinates for Classical CA
  #  g <- dJh %*% sva$v # Column Standard Coordinates for Classical CA
  
  dimnames(f) <- list(paste(Inames[[1]]), paste(1:M ))
  dimnames(g) <- list(paste(Jnames[[1]]), paste(1:M ))
  
  ##########################################################################
  # #
  # Calculating the total inertia, Pearson chi-squared statistic, its #
  # p-value and the percentage contribution of the axes to the inertia #
  # #
  ##########################################################################
  
  Principal.Inertia <- diag(dmu[1:M, 1:M]^2)
  
  Total.Inertia <- sum(Principal.Inertia)
  
  Percentage.Inertia <- (Principal.Inertia/Total.Inertia) * 100
  Total.Perc.Inertia.M <- sum(Principal.Inertia[1:M])
  
  
  Total.Perc.Inertia.M <- sum(Principal.Inertia[1:M])
  

  chisq.val <- qchisq(1-alpha, df= (I - 1) * (J - 1))
  
         ##########################################################################
         # #
         # Here we construct the two-dimensional correspondence plot #
         # #
         ##########################################################################
         
         par(pty = "plottype")
         
         par(pty = "s")
         plot(0, 0, pch = " ", xlim = scaleplot * range(f[, 1:M], g[, 1:M]),
              ylim = scaleplot * range(f[, 1:M], g[, 1:M]),
              xlab = paste("Principal Axis ", a1, "(", round(Percentage.Inertia[a1], digits = 2), "%)"), ylab = paste("Principal Axis ", a2, "(", round(Percentage.Inertia[a2], digits = 2), "%)"))
         
         text(f[,1], f[,2], labels = Inames, adj = 0, col = cols[1])
         points(f[, a1], f[, a2], pch = "*", col = cols[1])
         
         text(g[,1], g[,2], labels = Jnames, adj = 1, col = cols[2])
         points(g[, a1], g[, a2], pch = "#", col = cols[2])
         
         abline(h = 0, v = 0)
         title(main = paste(100 * (1 - alpha), "% Confidence Regions"))
        
         
         radii <- sqrt(qchisq(1 - alpha, 2)/(n * Imass))
         radij <- sqrt(qchisq(1 - alpha, 2)/(n * Jmass))
         #############################################################
         # #
         # Calculating the semi-axis lengths for the confidence #
         # ellipses #
         # #
         #############################################################
         hlax1.row <- vector(mode = "numeric", length = I)
         hlax2.row <- vector(mode = "numeric", length = I)
         hlax1.col <- vector(mode = "numeric", length = J)
         hlax2.col <- vector(mode = "numeric", length = J)
         if (M > 2){
           # Semi-axis lengths for the row coordinates in an optimal plot
           for (i in 1:I){
             hlax1.row[i] <- dmu[1,1] * sqrt((chisq.val/(n*Total.Inertia))*
                                               (1/Imass[i] - sum(a[i, 3:M]^2)))
             hlax2.row[i] <- dmu[2,2] * sqrt((chisq.val/(n*Total.Inertia))*
                                               (1/Imass[i] - sum(a[i, 3:M]^2)))
           }
           # Semi-axis lengths for the column coordinates in an optimal plot
           for (j in 1:J){
             hlax1.col[j] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                               (1/Jmass[j] - sum(b[j, 3:M]^2)))
             hlax2.col[j] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                               (1/Jmass[j] - sum(b[j, 3:M]^2)))
           }
         } else {
           # Semi-axis lengths for the row coordinates in a two-dimensional plot
           for (i in 1:I){
             hlax1.row[i] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                               (1/Imass[i]))
             hlax2.row[i] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                               (1/Imass[i]))
           }
           # Semi-axis lengths for the column coordinates in a two-dimensional plot
           for (j in 1:J){
             hlax1.col[j] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                               (1/Jmass[j]))
             hlax2.col[j] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                               (1/Jmass[j]))
           }
         }
         {
           
           
           
           eccentricity <- sqrt(1 - (dmu[a2, a2]/dmu[a1, a1] )^2)
           
           list(N = N, r = round(r, digits = 3), f = round(f, digits = 3), g = round(g, digits = 3), Chi.Squared = round(Chi.Squared, digits = 3), Total.Inertia = round(Total.Inertia, digits = 3), P.Value = round(p.value, digits = 3), Inertia = round(Inertia, digits = 3))
           
         }
         
         Asmadata2=matrix(c(1228,0,0,0, 299, 263, 259, 219, 188,737,491,0,4223,0,0,677,792,855,956,943,2408,1815,0,0,18910,0,3033,3431,3793,4230,4423,9087,9823,0,0,0,14535,2052,2442,2754,3225,4062,6632,7903,299,677,3033,2052,6061,0,0,0,0,2861,3200,263,792,3431,2442,0,6928,0,0,0,3342,3586,259,855,3793,2754,0,0,7661,0,0,3806,3855,219,956,4230,3225,0,0,0,8630,0,4221,4409,188,943,4423,4062,0,0,0,0,9616,4634,4982,737,2408,9087,6632,2861,3342,3806,4221,4634, 18864,0,491,1815,9823,7903,3200,3586,3855,4409,4982,0,20032),nrow=11)
         Asmadata2
         dimnames(Asmadata2)=list(paste(c("Strongly disagree ","disagree","agree","Strongly agree","1","2","3","4","5","time1","time2")),paste(c("Strongly disagree ","disagree","agree","Strongly agree","1","2","3","4","5","time1","time2")))
         Asmadata2
         
         data = Asmadata2
         alph = 0.05
         mindelta = 0.05
         maxdelta = 1.00
         int = (maxdelta - mindelta)/200
         
         delta = seq(mindelta, maxdelta, int)
         ds = length(delta)
         
         dim.one = vector(mode = "numeric", length = ds)
         dim.two = vector(mode = "numeric", length = ds)
         twod = vector(mode = "numeric", length = ds)
         
         dim.one.inertia = vector(mode = "numeric", length = ds)
         dim.two.inertia = vector(mode = "numeric", length = ds)
         twod.inertia = vector(mode = "numeric", length = ds)
         
         
         chisq = vector(mode = "numeric", length = ds)
         
         for (i in 1:ds){
           output = crca.exe(data, delta = delta[i])
           chisq[i] = output$Chi.Squared
           dim.one[i] = output$Inertia[1,2]
           dim.two[i] = output$Inertia[2,2]
           twod[i]=output$Inertia[2,3]
           dim.one.inertia[i] = output$Inertia[1,1]
           dim.two.inertia[i] = output$Inertia[2,1]
           twod.inertia[i]=output$Inertia[1,1]+output$Inertia[2,1]    
         }
         
         # Plot of the % principal inertia values for the first two dimensions
         
         plot(delta, twod, type = "l", ylim = c(0, 100), xlab = expression(delta), ylab = "% Contribution")
         lines(delta, dim.one, lty = 2)
         lines(delta, dim.two, lty = 3)
         #abline(v = c(0, 0.5, 2/3, 1), lty = 2)
         
         # Plot of the raw principal inertia values for the first two dimensions
         
         plot(delta, twod.inertia, type = "l", ylim = c(0, max(twod.inertia)), xlab = expression(delta), ylab = "Principal Inertia")
         lines(delta, dim.one.inertia, lty = 2)
         lines(delta, dim.two.inertia, lty = 3)
         
         plot(delta, chisq, type = "l", ylim = c(0, max(chisq)), xlab = expression(delta), ylab = "Divergence Statistic")
         
         abline(h = qchisq(1 - alph, (nrow(data) - 1)*(ncol(data) - 1)), lty = 2)
         
         
}

eccentricity <- sqrt(1 - (dmu[a2, a2]/dmu[a1, a1])^2)

if (region == 1){
  # Superimposing the confidence circles
  symbols(f[,a1], f[,a2], circles = radii, add = T, fg = cols[1])
  symbols(g[,a1], g[,a2], circles = radij, add = T, fg = cols[2])
} else if (region == 2){
  # Superimposing the confidence ellipses
  for (i in 1:I){
    regions(f[i,], xcoord = hlax1.row[i], ycoord = hlax2.row[i],
            col = cols[1])
  }
  for (j in 1:J){
    regions(g[j,], xcoord = hlax1.col[j], ycoord = hlax2.col[j],
            col = cols[2])
  }
}

ca.regions.exe(Asmadata2,region = 2)

for (i in 1:ds){
  output = ca.regions.exe(data, delta = delta[i])
  chisq[i] = output$Chi.Squared
  dim.one[i] = output$Inertia[1,2]
  dim.two[i] = output$Inertia[2,2]
  twod[i]=output$Inertia[2,3]
  dim.one.inertia[i] = output$Inertia[1,1]
  dim.two.inertia[i] = output$Inertia[2,1]
  twod.inertia[i]=output$Inertia[1,1]+output$Inertia[2,1]  
  
}

$hlax1.row
plot(Asmadata2,plottype = "classic",ell=TRUE)

apply(Asmadata2/sum(Asmadata2),1,sum)