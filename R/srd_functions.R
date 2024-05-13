#################################################
##   srd_functions.R functions file for srd.R  ##
#################################################

draw.srd <- function(d,f=matrix(data=1, ncol=1,nrow=nrow(d)),labelcell=0,new=TRUE, seed=NA,
          colour=2,col,crit,labelrectangle=1, Y=matrix(data=0, ncol=1,nrow=nrow(d)),
          project=c(0,0,0),cells=FALSE,w=matrix(data=1, ncol=1,nrow=nrow(d)),thin=0.5,zname=NA,
         title=NULL,whitespace=TRUE, ws_col, spanr, x, y,perm, borders,is.naz,bgtextcol,indep)
  {
  ##------------------------------------------------------------------------------------------------------
  ## d is data matrix of ones and zeroes
  ## number of data rows n and columns for q rectangles
  ## f is frequency weight count vector (default vector of 1's) (this is "weight" in srd.R)
  ## labelcell=0, no labelcells (default), 1= cell frequencies, 2= area frequencies, 3= % error
  ## new  re-fits a new SRD  (default TRUE). Otherwise the existing configuration is used 
  ## crit=1 least squares (default), 2=least absolute difference, 3= minus log-likelihood
  ## colour=1  pale blue/brown colours, 2=rainbow colours (default), 3=monochrome  (if cells=FALSE)
  ## colour=1,2,3,4 for pale yellow, pale blue, black and white intensity shading if cells=TRUE
  ##       = 4 for full non-see-through rectangles
  ## labelrectangle=0 no rectangle labelcell, =1 (default) inside rectangle, 2=arrow pointer
  ## Y vector of  variable for third axis in a 3-D projection, default  zeroes (gives no projection)
  ## w ?? is it?  thought it was second item of z?? is weight associated with  3-D projection axis.
  ## project - vector of viewing projection angles (degrees) for 3D projection. Default no projection (0,0,0)
  ## cells=FALSE. If TRUE gives representation of cells of the configuation and projected to mean cell value
  ##         
  ##---------------------------------------------------------------------------------------------------------
  
  # 
  # lst <- list("lsqs","labs","logl","chisq")
  # xxxx <- criterion
  # crit <- which(lst==criterion)
  #
  ##print(col)
  
  if(is.matrix(d)==FALSE) {d <- as.matrix(d)}
  if(is.matrix(Y)==FALSE) {Y <- as.matrix(Y)}
  if(is.matrix(f)==FALSE) {f <- as.matrix(f)}
  if(is.matrix(w)==FALSE) {w <- as.matrix(w)} 
  
  
  
  fail <- FALSE
 
  ## n_in and n_out refer to number of records: do not account for w weight
  
  n_in <- nrow(d)
  
  freq_in <-  sum(f,na.rm=TRUE) 
  
  n_out <- n_in
  freq_out <- freq_in 
  q <- ncol(d)

 
## Delete rows of d if incomplete
## shouldnt need this as cc check already done. Patch shutdown
rnames <- colnames(d)

notes <- NULL
if(new){
  

if(freq_out <= 0){
  return(message("Error: Negative or zero sample size?"))
}

##  if new=TRUE  create a new srd configuration by invoking FORTRAN
##  subroutine srd.

   if(q>6){
     ## documentation uses "k" not q
      message(paste("Cannot create an SRD for", q, "rectangles. Up to 6 allowed"))
      
      fail <- TRUE 
      output <- list(fail=fail)
      return(output)
   }
##test for strictly binary  numeric variables
## if any d[,i] has been made a factor, entire d is "character"
## test for this, using first column. 
## This seems flakey, but have checked it out
if(is.character(d[,1])){
message("non-numeric data. Must be binary for rectangle attributes, numeric y and weight")
  fail <- TRUE 
  output <- list(fail=fail)
  return(output)
}

    for (i in 1:q){
      
      if(q>1){
    test <- unique(d[,i])
    sumdi <- sum(d[,i])}
    else
    {test <- unique(d[,])
      sumdi <- sum(d[,])}
      if(sumdi==0){
        ## though zero rectangles filtered in srd.R they
        ## may remain due to missing values after reduction to complete cases
        message(paste("warning: empty rectangle",colnames(d)[i] ))
      }
    
    ##browser()

    if(length(test)!=2){
      
      if(length(test)==1){
        
        ## possibity of single 0 or 1, not good but is allowed, seems OK
        if(test[1] !=1 & test[1] !=0){
          v <- colnames(d)[i]
          message(paste("Variable ",v, " not binary"))
          fail <- TRUE 
          output <- list(fail=fail)
          return(output)
        }else{
      ##    browser()
          v <-  colnames(d)[i]
          note <- paste("Variable ",v, "is entirely value",test)
          message(note)
          notes <- append(notes,note)
        }
      }else{

      v <-colnames(d)[i]
      message(paste("Variable ",v, " not binary"))
      fail <- TRUE 
      output <- list(fail=fail)
      return(output)
    }
    }
    
    else
    {if((test[1]!=0 | test[2] !=1)&(test[1]!=1 | test[2] !=0)){
      v <-colnames(d)[i]
      message(paste("Variable ",v," is two-class but not coded 0 1"))
     fail <- TRUE 
      output <- list(fail=fail)
      return(output)
      
    }}  
    
    }

    # 
    # x <- matrix(nrow=4,ncol=7)
    # y <- matrix(nrow=4,ncol=7)
    # 
    # ## temp routine to match granm in srd.f
    # ## testing process to input into srd.f
  # browser()
  #    XXX <- grandm(d,f,q,x,y)
  # 
  #   x <- XXX[[1]]
  #   y <- XXX[[2]]
    cstat <- matrix(nrow=5,ncol=64)

#   required random number seed for FORTRAN 77 srd subroutine to be
#   passed as argument to establish random starting positions     
    if(is.na(seed)){
    seed <- round(100000*runif(1))
    }     


# 
# for(i in 1:q){
#   size <- sum(d[,i]*f)
# 
# if(size ==0 )
# {message(paste("No data for ", colnames(d)[i], " attribute. Remove from rectangle argument"))
#  
# 
#   }
# }
# ##  remove empty rectangles
# ## b try automatic removal. Running into problems
# if(!is.null(remove))
# {
#   return()
#   q <- q - length(remove)
# d <- d[,-remove]
# rnames <- colnames(d)
# ##browser()
# }
# 

doptm <- 0
dactm <- 0
thinx <- 0  ##returned maximal length:breadth
##names(f) <- "freq"
##df <- cbind(d,f)
##colnames(df)[ncol(df)] <- "freq"
##df <- as.data.frame(df)
  ##browser()
  
  d_ <- d
  f_ <- f
  
  ## BIG data ? Collapsing again?? Why?
  ## because previous collapse includes y variable (s)
  ## This collapse excludes the y variable
  ## which is not part of  .Fortran("srd"  to "fit" the srd.
  ##print(nrow(d))
  ## always do? OK
  if(nrow(d)> 10000 | TRUE){
# 
   df <- as.data.frame( cbind(d_,f_) )
# 
   cw <- colnames(df)[ncol(df)]
  dfx <- SumDupRows(df, 1:(ncol(df)-1), cw,cw)
# 
   d_ <- dfx[ ,  1:(ncol(dfx)-1)  ]
   d_ <- as.matrix( d_, nrow=nrow(dfx), ncol= (ncol(dfx)-1) )
   f_ <- dfx[,ncol(dfx)]
   n_out <- nrow(dfx)

  }
##print(thin)
  
  ##testing patch for an "independence srd".


  if(indep){
  
   XXX <- expectedf(d_, f_)
   d_ <- XXX[[1]]
   f_ <- XXX[[2]]
   n_out <-  nrow(d_)

##print(cbind(d_,f_))
}
#   
#   if(any(!is.wholenumber(f_))){
# #    print(f_)
#     message('non-integer frequencies ')
# ##   f_ <- round(f_,2)
# ##   print(f_)
#   }
##print(f)
  
  ## patch for crit= 5  (stress) for which srf.f isn't functional. 
  ## use lsqs instead.
 
    critX <- crit
    if(crit == 5) critX <- 1
    run_srd <- .Fortran("srd",d=as.integer(d_),f=as.single(f_),n_out=as.integer(n_out),
                q=as.integer(q),x=as.double(x),y=as.double(y),e=as.single(0),cstat=as.single(cstat),
                crit=as.integer(critX),seed=as.integer(seed),thin=as.single(thin), 
                doptm=as.double(doptm),dactm=as.double(dactm),thinx=as.single(thinx),NAOK=TRUE )

   ## print(paste("dactm",signif(run_srd$dactm,3),"doptm",signif(run_srd$doptm,3) ))
    ##print(run_srd$x)
 ##print(run_srd$y)
# save data and configuration as global objects, for subsequent new=FALSE
## doptm and dactm are return values of the criterion.
## doptm is "adjusted/penalised value"  for thinness and empty/missing cells
## , dactm is the actual value
    
    ##==============================================================
    ## use fitted .Fortran as start values to tweak?
    # ## are odd numbered elements
    # 
    # testit <- FALSE
    # 
    # 
    # if(testit){
    #   message("tweaking fit with 4 parameter per rectangle")
    # s1 <- seq(from=1,to=4*(q),by=4)
    # s2 <- s1 +1
    # s3 <- s1 +2
    # 
    # Xpairs <- run_srd$x[sort(c(s1,s3))]
    # Ypairs <- run_srd$y[sort(c(s1,s2))]
    # ## need to get theta as long parameters in blocks of 2. 
    # 
    # theta <- c(Xpairs,Ypairs)
    # ##print(theta)
    # 
    # XXX <- fit_srd(d_,f_,theta)
    # run_srd$cstat <- XXX[[1]]
    # run_srd$x <- XXX[[2]]
    # run_srd$y <- XXX[[3]]
    # }
    # ##=============================================================
     ## use fitted .Fortran as start values to tweak?
    ## with forced rectangle area
   
     
 ### testing lines----------------------------   
   sq <- 5*(seq(1:2^q)-1)
     xcentreX <- run_srd$cstat[1+sq]
     ycentreX <- run_srd$cstat[2+sq]
     areaX <- run_srd$cstat[4+sq]
  # ##  print("run_srd area")
  #   print(round(areaX,3))
     freqX <- run_srd$cstat[3+sq]
  #   
     errorX <- abs(freqX/sum(freqX) - areaX/sum(areaX))
     if(!whitespace) errorX <- errorX[-1]
     diagErrorX <- max(errorX)
  #   lsqs <- sum((areaX-freqX)^2,na.rm=TRUE)
  #   print(paste("run_srd lsqs", signif(lsqs,5)) )
##=========================================================================================
## attempts Oct 2020 to tweak fitted .Fortran()  configuration by
## Powell's algorithm, using  R optim() Nelder-Mead using fitted
##  configuration as starting position.     
## upshot is I am unconvinced that it improves and maybe even worsens fit.
##  one option is to tweak to  making the criterion the diagError measure
##  crit=6 in fit_srd2. Or maintain the same fitted criterion in the tweak.       
## has be a lot of wasted effort to develop the architecture to do "within" 
## R fitting     
    testit2 <-  diagErrorX > 0.001
  
    if(T){
    ##print(paste("pretweak diagError =", signif(diagErrorX,4)))
    ##print("tweaking")  
     ## message("tweaking fit with 3 parameter per rectangle")
    
        
      s1 <- seq(from=1,to=4*(q),by=4)
      s2 <- s1 +1
      s3 <- s1 +2
      
      Xpairs <- run_srd$x[sort(c(s1,s3))]
      Ypairs <- run_srd$y[sort(c(s1,s2))]
      ## need to get theta as long parameters in blocks of 2. 
      Yvalues <- Ypairs[2*seq(1:q)-1]
      theta <- c(Xpairs,Yvalues)
##      print(run_srd$dactm)
 ## metric  5  is "stress" 
      ##-----------------------------------------------------
      ###  fit using  optim() function
      ## tweak wrt  diagError  crit=6
       XXX <- fit_srd2(d_,f_,theta,crit=crit)
      ##-----------------------------------------------------
      ##utilise these tweaked values to maintain code structure
      run_srd$cstat <- XXX[[1]]
      run_srd$x <- XXX[[2]]
      run_srd$y <- XXX[[3]]
      run_srd$dactm <- XXX[[4]]   #minimised alue of criterion
      theta <- XXX[[5]]
##      print(XXX[[4]])
      ##run again  with output theta and diagErrot
      # 
      # XXX <- fit_srd2(d_,f_,theta,crit=6)
      # ##-----------------------------------------------------
      # ##utilise these tweaked values to maintain code structure
      # run_srd$cstat <- XXX[[1]]
      # run_srd$x <- XXX[[2]]
      # run_srd$y <- XXX[[3]]
      # run_srd$dactm <- XXX[[4]]   #minimised alue of criterion
      # theta <- XXX[[5]]
      
    }   ##  if(testit2)
##=========================================================================================
 
 dsave <- NULL
 run_srd_save <- NULL
dsave <<- d
run_srd_save <<- run_srd 

}
  else #  for if(new)   
{
##     message("re-creating existing configuration")
  
  ## new argument may be the input new parameter, but if  in
  ## rotation mode it is also set= FALSE
     if(is.na(run_srd_save[1])==TRUE){
      note <- c(paste("new=FALSE. No previous saved configuration in workspace"))  
      message(note)
      
      fail <- TRUE 
      output <- list(fail=fail)
      return(output)
     }
    if(identical(d,dsave)==FALSE){
      note <- c(paste("new=FALSE expects matching configuration. Re-run with new=TRUE"))  
      message(note)
      fail <- TRUE
      output <- list(fail=fail)
      return(output)
     }
    
    run_srd <- run_srd_save  
    new <- FALSE 
} 

##  Configuration of  rectangles on flat plane in run_srd$x and run_srd$y

## now consider  the 3D Y axis 
##browser()
## could coerce to matrix structure of input with
## matrix(run_srd$x,nrow=4,ncol=7)
##browser()
rot_x <- run_srd$x
rot_y <- run_srd$y
rot_x_ws <- run_srd$x
rot_y_ws <- run_srd$y


z<-  matrix( ncol=1,nrow=28)  
   
   ## account for missing Y

##browser()
   okY <- complete.cases(Y)

   Y <- Y[okY]
   nY <- length(Y)
  

   f <- f[okY]

    test <- sum(f)
if(test < freq_out & identical(project,c(0,0,0))) {
notes <- append(notes, paste("Missing y: y projections based on sub-sample of",test))
}

## set area threshold for cell to be "visible"
## less than 1/10000 of whole 
# Or bigger. freq_out is essentially sample size
    
athreshold <- freq_out*1.0e-4

## 0.05% of area? 
athreshold <- freq_out*5.0e-4
## when is  cell "visible"? 
##athreshold <- freq_out*1.0e-4 
w <-  w[okY]

   if(q>1)
  { 
   d <- d[okY,]}
   else{
     d <- d[okY]
   }

##   obtain z values of the q rectangles, as mean of Y
##   same at each of 4 rectangle corners
##   use dsum as indicator as whether inside one of the rectangles
##   used later to establish whitespace mean.
  dsum <- vector(length=length(Y))

dsum <- 0
   for(i in 1:q){
    ftest <- f 
 ##   d should be matrix, but may not be if q=1, or  only one row
  if(q>1 & is.matrix(d)){  

  YFi  <- Y*f*d[,i]  
  wFi  <- w*f*d[,i]
  
  dsum  <-  dsum +d[,i]}
  else
    {   YFi  <- Y*f*d  
        wFi  <- w*f*d
  dsum  <-  dsum + d}
  
  sumn=sum(YFi,na.rm=TRUE)
  sumd=sum(wFi,na.rm=TRUE)
  ##rectangle mean of projection z, or ratio of means if w != 1
  mean<-sumn/sumd
  if(is.naz) mean <- i
  z[1+(i-1)*4]=mean
  z[2+(i-1)*4]=mean
  z[3+(i-1)*4]=mean
  z[4+(i-1)*4]=mean
  


   }  #for(i 1:q)
#----------------------------------
 ## binary indicator of A v !A 
##redherring?  meanA, meanAc never used!

if(q >1) {
A  <- ifelse(rowSums(d)==0,0,1)
}
else
{A <- ifelse(d==0,0,1)
}
## calculate means of Y in A and in !A
meanA <- sum(Y*f*A)/sum(w*f*A)
meanAc <- sum(Y*f*(1-A))/sum(w*f*(1-A))
#----------------------------------


  ## q+1 th is the 4 corners of surrounding unit square
  #overall mean 
  Yf <- Y*f
  wf <- w*f
   omean <- sum(Yf,na.rm=TRUE)/sum(wf, na.rm=TRUE)
#  mean of whitespace,  1-dsum is 0,1 indicator of being in white space

dac <- ifelse((dsum >0) == TRUE,0,1)
YFc <- Y*f*dac  
wFc <- w*f*dac
##browser()
if(sum(wFc,na.rm=TRUE)==0) {
  ws_mean <- omean
  nowhitespace <- TRUE
  
  ## nowhitespace? force whitespece=F?
  if(whitespace){
    notes <- append(notes, paste("Data has no whitespace. Consider using whitespace=FALSE"))
  }
}
else
{nowhitespace <- FALSE
ws_mean <- sum(YFc,na.rm=TRUE)/sum(wFc,na.rm=TRUE)
}

nowhitespace <- !whitespace

#browser()

  z[1+(q)*4]=ws_mean
  z[2+(q)*4]=ws_mean
  z[3+(q)*4]=ws_mean
  z[4+(q)*4]=ws_mean


  zmax<-max(z,na.rm=TRUE)
  zmin<-min(z,na.rm=TRUE)

fudge <- 0.5

#rescale z  when not cells
if(zmax == zmin & !cells){    
  
##  return(message("No variation in y. 3D projection pointless"))
  fudge <- 0
                     zmin_scaled <- zmax
                     zmax_scaled <- zmin
                     ws_scaled <-   zmax   
}

#----------------------------------------------------------------------
## calculate testYf regardless of cells=T or F. KANEON!!    
####  if(cells==TRUE){
## account for frequencies of frequency f weights 
    Yf <- Y*f
    wf <- w*f

    dYf <- as.data.frame( cbind(d,wf,Yf) )
    ## if cells  need to establish mean values in each cell  
    

## patch up if q=1

if(q==1) {colnames(dYf)[1] <- rnames} 
    
    
    colnames(dYf)[ncol(dYf)]   <- c("Yf")
    colnames(dYf)[ncol(dYf)-1] <- c("wf")
    #ddply fussy about ! character in name. Replace with X
    ##savedcolnames <- colnames(dYf[1:q])


    colnames(dYf) <- gsub(pattern="!",replacement="X",x=colnames(dYf))
    colnames(dYf) <- gsub(pattern=" ",replacement="_",x=colnames(dYf))
    colnames(dYf) <- gsub(pattern="-",replacement="m",x=colnames(dYf))
##Billy Wi's sum function
## take duplicate rows and sum values of "Yf" over the duplicates, returning value "V1" for summ. 
    NUM <- SumDupRows(dYf, 1:(ncol(dYf)-2), "Yf", "V1")
    ## take duplicate rows and sum values of "wf" over the duplicates, returning value "V1" for sum.
    ## normally unless z is two-valued, wf is column of ones
    DEN <- SumDupRows(dYf, 1:(ncol(dYf)-2), "wf", "V1")
    ## V1 is mean of z, if  z single valued, or ratio of sums
    ## now bind this into the data of unique cell combinations.
    ## giving the z projection by cell
    V1 <- NUM$V1/DEN$V1
    if(q>1) {
    testYf <- cbind( NUM[,1:(ncol(dYf)-2)],V1 )}
else
    {
      ## need to patch when q=1
      testYf <- cbind( NUM[1:(ncol(dYf)-2)],V1 )}
##browser()
   
sumYf <- NA
if(labelcell == 8){  ##calculate cell sum of y, numerator sumn        
    sumYf <- NUM} 
if(labelcell == 9){  ##calculate cell sum of w  denominator sumd
   sumYf <- DEN }  

### LATER note:
## trying to re-think this.  If z input into srd is single valued
## which is passed as "y" here NUM is sum (y) and DEN is n number in the cell
## but if z is dual valued with passed values "y"  and "w" then
##  DEN is sum(w). & V1 is ratio of the sums (equivalent to ratio of means)
## 

##browser()

    ## creates object testYf$V1 of means for each binary combination
    if(identical(project,c(0,0,0))){
      

     if(length(which(is.infinite(testYf$V1)))>0){
      fail <- TRUE
 ##   message(paste("Warning: Infinite values of y in ", length(which(is.infinite(testYf$V1))),
 ##                 "cells"))
    
    output <- list(fail=fail)
    ##return(output)
    }

    
    if(length(which(is.nan(testYf$V1)))>0){
      fail <- TRUE
##      message(paste("Warning: NaN values of y in ", length(which(is.nan(testYf$V1))),
##                    "cells "))
      
      output <- list(fail=fail)
      ##return(output)
    }
    }
##-------------------------------------------------------------------
## STUFF to WRITE OUT as flat, no projection
##create outtable
outtable <- NA
outrect <- NA

##  do regardless of angle
## probably this is inefficient -  really only need do
## outtable on Esc.   But WTF.
if(identical(project,c(0,0,0))   | TRUE ){
icell <- as.numeric(vector(length=2^q))
icell <- NA
ncell <- length(testYf$V1)
##    if(zmax != zmin){
##gap <- paste(rep("  ",times=q),collapse='')
##{cat("\n", paste("Cell ", gap,"Frequency : area" ))}
nabs <- 0
nemp <- 0
emptyarea <- 0


combo <- matrix(0, ncol=q,nrow=2^q)

freq <- as.numeric(vector(length=2^q))
area <- as.numeric(vector(length=2^q))
diff <- as.numeric(vector(length=2^q))

error <- as.character(vector(length=2^q))
yval <-  as.numeric(vector(length=2^q))
resid <- as.numeric(vector(length=2^q))

nabsj <- as.numeric(vector(length=q))
nempj <- as.numeric(vector(length=q))
nabsj <- rep(0,times=q)
nempj <- rep(0,times=q)
##for outtable, should it include y values?

tst <- unique(testYf$V1)
tst <- which(is.nan(tst)==FALSE)
## possibility of a NaN returned

if(length(tst) >1){ 
## need to fiddle around to get  Y values of  all cells (not just non-empty)
icell <- as.numeric(vector(length=2^q))
icell <- NA
n_nonempty_cells <- nrow(testYf)
for (i in 1:n_nonempty_cells ){
  
  bin <- as.numeric(testYf[i,1:q])
  
  
  ibin <- sum(2^(which(bin==1)-1)) +1 
  ## keep record of  cell that ibin refers to 
  icell[ibin] <- i
}

for (i in 1:2^q){
  if(is.na(icell[i])==FALSE) {yval[i] <- testYf$V1[icell[i]]}
}
}



# ##first set area and frequency vectors, recompute frequency sum
# for (i in 1:(2^q)){
# ## Statistics cstat return from FORTRAN call in long arrays
# area[i] <- run_srd$cstat[4+5*(i-1)]
# 
# freq[i] <- run_srd$cstat[3+5*(i-1)]
# 
# resid[i]=(freq[i]-area[i])/sqrt(max(area[i],0.5))
# }
##print(freq)

sq <- 5*(seq(1:2^q)-1)
area <- run_srd$cstat[4+sq]
freq <- run_srd$cstat[3+sq]

## replace 0 areas by 0.5
resid <- (freq-area)/sqrt(ifelse(area==0,0.5,area))

# 
 total <- sum(freq)
 tarea <- sum(area)
 
## diagERRortest <- max(abs(area/tarea -freq/total))
 ## print(paste("diagERRortest",diagERRortest))
##--------------------------------------------------------------
## experimental steps in R verification
# verifyarea <- total*cellgrid(unique(run_srd$x),unique(run_srd$y),q)
# print(verifyarea)
# 
# ## construct theta parameter from elements of x,y. These are
# ## each blocks of 4 . x length 4*(q+1)  Need to extract  pairs from each. For x these
# ## are odd numbered elements
# s1 <- seq(from=1,to=4*(q+1),by=4)
# s2 <- s1 +1
# s3 <- s1 +2
# 
# Xpairs <- run_srd$x[sort(c(s1,s3))]
# Ypairs <- run_srd$y[sort(c(s1,s2))]
# ## need to get theta as long parameters in blocks of 2. 
# theta <- c(Xpairs,Ypairs)
# 
# print(fun_optim(theta,freq,crit))
# 
# ##chisqcheck(freq,area)
# ##XXX <- nlm(fun_optim,theta,freq)
# ##print(XXX$minimum)
# XXX <- optim(theta,fun_optim,gr=NULL,freq,crit)
# print(XXX$value)
# 


##--------------------------------------------------------------
## check on criterion=chisq
## print(sum(resid^2))

## without whitespace trim off the 0,0,0,.. cell 
 
 if(!whitespace)resid <- resid[-1]
chisqu <- sum(resid^2,na.rm=TRUE)
# 
# print(paste("chisqu",signif(chisqu,4))  )
# print(freq)
# print(area)

den <- ifelse(area<=0.5,0.5,area)
## use chisq or gof??
gof <- sum(  (area-freq)^2/ (den*(1-den/total)) ,na.rm=TRUE)

dfgof <- length(which(area !=0 | freq !=0)) - !whitespace
pgof <- 1 - pchisq(chisqu, df=dfgof)



nabs <- 0
absfreq <- 0
missedcell <- NULL
emptycell <- NULL
for (i in (1+nowhitespace):(2^q)){
   
  bin <- rev(number2binary(i-1,q))

  combo[i,] <- bin
  
    
  diff[i] <- 100*(freq[i]-area[i])/total
 
  ##  error[i] <- paste(round(diff[i],2))
  error[i] <- signif(area[i]/freq[i],2)

  ## missing cell?
 if(freq[i]>0 & area[i]==0) {
  nabs <- nabs + 1
  absfreq <- absfreq + freq[i]
 error[i] <- paste("missing ",error[i]) 
 ## assign greater visual weight to missing cell

 missedcelli <- paste(colnames(testYf)[ which(as.logical(bin)==1)],collapse="&")
 missedcell <- c(missedcell,missedcelli)
 if(i == 1) emptycelli <- "whitespace"
 note <- paste("missed cell:", missedcelli,paste0("(n=", signif(freq[i],2),")") )
 message(note)
 notes <- append(notes,note )

 for(j in 1:q){
  if(bin[j]==1){ nabsj[j]= nabsj[j] +1}
 
 }
 
 }

if(freq[i]==0 ){
##empty area  
## possibility of whitespace showing thro is not empty
 if( area[i]>athreshold ) {
 

 error[i] <- paste("empty",error[i])
 ## consider cells so small 0.05? as to be not visible when freq=0 to be "no data"
 
emptycelli <- paste(colnames(testYf)[which(as.logical(bin)==1)],collapse="&")
if(i == 1) emptycelli <- "whitespace"
 emptycell <- c(emptycell,emptycelli)
 note <- paste("empty area:", emptycelli, paste0("(area=",signif(area[i],2),")") )
 message(note)
 notes <- append(notes,note)
 
 
 for(j in 1:q){
   if(bin[j]==1){ nempj[j]= nempj[j] +1}  
 }
  } ##if( area[i]>athreshold)
  
  else 
    
 {
   
   ## zero frequency and physical area small
   ## effectively no data
  error[i] <- "no data  -"
}
}  ##if(freq[i]==0)
  
  
##output cell Y values?? No. 
##if(labelcell!=0){cat(paste(sep="", "     (", signif(testYf$V1[i],6),")") )}

if(freq[i]==0 & area[i]==0) {error[i] <- "no data  -"
}
## also ignore if whitespace, for which area[i] may be >0
##if(i==1 & freq[i]==0){error[i] <- "no data"}

if(substr(error[i],1,5)=="empty") {nemp <- nemp +1
emptyarea <- emptyarea +area[i]
}

}  ##for(i in 1:2^q)

 ## browser()

## calculate E% from  cell statistic (should agree with e)
if(whitespace){ 
  E_calc <- sum(abs(diff))

  diagError <- max(abs( area/sum(area) - freq/sum(freq) ))
  }
else
{
  E_calc <- sum(abs(diff[-1]))
  diagError <- max(abs(area[-1]/sum(area[-1]) - freq[-1]/sum(freq[-1])))
  
  error[1] <- "whitespace ignored" 
  area[1] <- 0
  total <- total - freq[1]
  tarea <- tarea - area[1]
}
## print(paste("diagError", signif(diagError,3)))
   ## whitespace could poss be label empty area? Need to remove
## WHY?? temp shutdown
#   if(substr(error[1],1,5)=="empty") {nemp <- nemp -1
# emptyarea <- emptyarea - area[1]}
#    if(substr(error[1],1,6)=="missing") {nabs <- nabs -1
# absfreq  <- absfreq - freq[1]}
 


## check on full rectangle areas
# gfreq <- as.numeric(vector(length=q))
# garea <- as.numeric(vector(length=q))
# 
# for(j in 1:q){
# garea[j] <- sum(area[which(combo[,j]==1)])
# gfreq[j] <- sum(freq[which(combo[,j]==1)])
# 
# }

garea <- apply(combo,2, function(x,a=area) {sum(a[which(x == 1)])} )
gfreq <- apply(combo,2, function(x,f=freq) {sum(f[which(x == 1)])} )




Rlabels <- paste("R",seq(1:q),sep="")
if(length(tst)>1) {
 ## marginal means of y in z in blocks of 4  
  rmeans <- z[seq(1, length(z), by = 4)]
  rmeans <- rmeans[1:q]
outrect <-  as.data.frame(cbind( rnames,round(gfreq,1),round(garea,1),nabsj,nempj,signif(rmeans,3)),drop=FALSE)
names(outrect) <- c("Rectangle","Freq", "Area","missing cells", "empty areas",paste(zname,sep=""))


}else{
  outrect <-  as.data.frame(cbind(rnames,round(gfreq,1),round(garea,1),nabsj,nempj))
  names(outrect) <- c("Rectangle","Freq", "Area","missing cells", "empty areas")
}

if(nemp==0)   {outrect <- subset(outrect, select=c(-5))}
if(nabs==0)   {outrect <- subset(outrect, select=c(-4))}

outrect <- cbind(Rlabels,outrect)
names(outrect)[1] <- ""

##print(outrect, row.names=FALSE)
discrep <- max(abs(garea-gfreq))
  sfreq <- round(sum(freq,na.rm=TRUE),1)
  sarea <- round(sum(area,na.rm=TRUE),2)

area <- round(area,2)

combo <- ifelse(combo==1,"+",".")
if(length(tst)>1) {
  yval <- ifelse(freq==0,NA,yval)
  ymean <- signif(mean(yval,na.rm=TRUE),3)
  yval <- signif(yval,3) 
  ymean <-  signif(ymean,3) 
  outtable <- cbind(combo,round(freq,1),area,error,yval)

names(outtable)[q+4]  <- paste(zname,sep="")



outtable <- rbind(outtable,outtable[1,])


nr <- nrow(outtable)
outtable[nr,1:q] <- rep("",times=q)
outtable[nr,(q+1)] <- sfreq
outtable[nr,(q+2)] <- sarea
outtable[nr,(q+3)] <- paste(signif(sarea/sfreq,2))
## ymean z projection cell value
outtable[nr,(q+4)] <- ymean

outtable <- as.data.frame(outtable)

}else{  #if(length(tst)>1) 


outtable <-  cbind(combo, round(freq,1) ,area,error)

outtable <- rbind(outtable,outtable[1,])


nr <- nrow(outtable)
outtable[nr,1:q] <- rep("",times=q)
outtable[nr,(q+1)] <- sfreq
outtable[nr,(q+2)] <- sarea
outtable[nr,(q+3)] <- paste(signif(sarea/sfreq,2))
outtable <- as.data.frame(outtable)


}

## add in y to outtable.First decided whether y varies (it is 0s if y=NA in spanr()


red_names <- rnames 
for(i in 1:q){

  
  if(nchar(rnames[i])>15 ){

    red_names[i] <- paste(substring(rnames[i],1,15),"~ ",sep="")
  }}
##names(outtable)[(1):(q)] <- red_names
## put R1 R2 R3....along to 
names(outtable)[1:q] <- paste("R",seq(1:q),sep="")

names(outtable)[(q+1):(q+3)] <- c("Freq n_i","Area a_i", "Area/Freq")
 if(length(tst)>1)names(outtable)[q+4] <- paste(zname,sep="")
## only output  cells with data  ??
##outtable <- outtable[which(outtable[,(q+3)] !="no data"),]

Clabel <- paste("C",seq(1:nrow(outtable) ),sep="")
Clabel[nrow(outtable)] <- ""
outtable <- cbind(Clabel,outtable)
names(outtable)[1] <- ""
   ##browser()
if(nrow(outtable)<2^q) {
 notes <- append(notes, 
  paste("Of 2^",q,"=",2^q," cell combinations, ",2^q-nrow(outtable), " are empty and correctly missing",sep=""))
} 

##if(freq[1]==0)cat(" Data has no \"whitespace\" ","\n")

##print(outtable, row.names=FALSE)


if(discrep>0) {
 ## cat("\n", paste("Max whole rectangle |area-freq| discrepancy < ",
##                  round(discrep,digits=8)))
}

if(nabs>0){ 
  rabsfreq <- round(100*absfreq/total,2)
  notes <- append(notes, paste(nabs," missing cells comprise ",rabsfreq
                               ,"% of sample (",round(absfreq,1),"/",round(total,1),")",sep=""))
  ##absfreq <- rabsfreq
  }
 
if(nemp>0){ 
  emptyarea <- round(100*emptyarea/tarea,2)
  notes <- append(notes, paste(nemp," empty areas comprise ",emptyarea,"% of plotting area",sep=""))}
  
}
##--end of STUFF to write out-------------------------------

     
## strip out NaNs. Is this safe KANEON

##browser()

##      nonnans <- which(is.nan(testYf$V1)==FALSE & is.infinite(testYf$V1)==FALSE)
##      testYf <- testYf[nonnans,]

#browser()
ncell <- length(testYf$V1)
## hang on to  testYf before it is re-scaled for plotting & pass it in draw.srd.enc
testYf_actual <- testYf$V1
#@#print(testYf_actual)
##browser()
if(labelcell == 8 | labelcell == 9){

  testYf_actual <- sumYf$V1}
#

if(cells==TRUE) { 
      zmax_cells <- max(testYf$V1,na.rm=TRUE)
      zmin_cells <- min(testYf$V1,na.rm=TRUE)
## fudge is factor to ensure   3d axis contained within plotting area      
fudge <- 0.5 *(zmax-zmin)/(zmax_cells-zmin_cells)
##browser()
      if(zmax !=zmin){
      testYf$V1<- (testYf$V1/(zmax-zmin))*fudge

      }
    
} #if(cells)
#------------------------------------------------------------------
z_actual <- z
z <- (z/(zmax-zmin))*fudge 

## z is possibly NA if cells and zmin-zmax
## how to fix? 
if(zmax==zmin) z <- ws_mean*fudge
##browser()
if(cells)
  {

zmax <- zmax_cells
zmin <- zmin_cells 
zmin_scaled <-min(testYf$V1,na.rm=TRUE)
zmax_scaled <-max(testYf$V1,na.rm=TRUE)

}
else
{
   zmin_scaled <- min(z,na.rm=TRUE)
   zmax_scaled <- max(z,na.rm=TRUE)
}


ws_scaled <- 0
if(zmax != zmin) ws_scaled <- (ws_mean/(zmax-zmin))*fudge
#if(zmax==zmin) z <- ws_mean*fudge

# 28 ?   = 7 * 4 corners, 7= up to 6 rectangles + 1 border
##browser()
  for (i in 1:28){ 
    ## patch to ensure  rotating SRD stays "flat"
  if(zmin==zmax) z[i] <- 0
  xyz <- c(run_srd$x[i],run_srd$y[i],z[i])
  test <- rotate(xyz,project)
#  run_srd$x[i]<- test[1]
rot_y[i]<- test[2]
rot_x[i]<- test[1]
##  z[i]<- test[3] 
## now for back projection onto whitespace to make "hole" 
xyz <- c(run_srd$x[i],run_srd$y[i],ws_scaled)
test <- rotate(xyz,project)
#  run_srd$x[i]<- test[1]
rot_y_ws[i]<- test[2]
rot_x_ws[i]<- test[1]

  }


 
 ncell <- 2^q



#------------------------------------------------------------------
## collapse to get cell combination and frequencies

  df <- as.data.frame( cbind(d,f) )
  ## if cells  need to establish mean values in each cell  
colnames(df)[ncol(df)] <- c("f")
colnames(df) <- gsub(pattern="!",replacement="X",x=colnames(df))
colnames(df) <- gsub(pattern=" ",replacement="_",x=colnames(df))
colnames(df) <- gsub(pattern="-",replacement="m",x=colnames(df))


cellf <- SumDupRows(df, 1:(ncol(df)-1), "f", "V1")

##browser()
## and group frequencies
##browser() 

groupf <- vector(length=q) 
for(i in 1:q){ groupf[i] <- sum(cellf[,i]*cellf$V1)}


##--------------------------------------------------------------

  if(new==TRUE){
   #browser() 
   ocrit <- signif(run_srd$dactm,7) 
   criteria <- c("lsqs","labs","logl","chisq")
   criterion <- criteria[crit]
   notes <- append(notes, paste("Optimised ",criterion,"=",ocrit))
   
  
    ##browser()
   
   if(!is.null(E_calc)) {E_calc <- round(E_calc,digits=2)}
   ## need to hang on to the E% to be output in rotating
   ##  do this by adding into saved run_srd object
     run_srd_save$E_calc <<- E_calc
    good <- ""
    ##if(E>10)  {good <- "(not so good)"}
   ## if(E<10)  {good <- "(some incongruence)"}
  ##  if(E<5)   {good <- "(good enough)"}
    if(E_calc<0.01) {good <- "(exact)"}
    
  
  
  notes <- append(notes,paste("Total area-to-frequency abs error =",E_calc*0.01))

 ## paste("gof =",signif(gof,5),"has P=",signif(pgof,5),"df=",dfgof ) )
 # if(run_srd$E_calc >10){
 #    notes <- append(notes,paste("Note E% =",run_srd$E_calc,
 #    " may be unacceptable"))
 #  }
  
  notes <- append(notes,paste("max(abs(error)) =",signif(diagError,3) ) )
  
    notes <- append(notes,paste("Maximum rectangle length:breadth = ",round(run_srd$thinx,digits=2),":1",sep=""))
  # if(run_srd$thinx > 5)
  # {notes <- append(notes,paste("Thin rectangle(s): could increase thin parameter (current= ", run_srd$thin,")"))}
  
    
  }
#-------------------------------------------------------------------------------------
## call the drawing function


## form ordered sub-cell grid vectors as global
#if(identical(project,c(0,0,0))){
  nn <- 2*(q+1)
  
  xxgrid <- vector(length=nn)
  yygrid <- vector(length=nn)
  nn1 <- nn-1 
  for(k in 1:(q+1)){
    k4 <- 4*k
    k2 <- k*2
    xxgrid[k2-1] <- run_srd$x[k4 -3]
    xxgrid[k2]   <- run_srd$x[k4]
    yygrid[k2-1] <- run_srd$y[k4-3]
    yygrid[k2]   <- run_srd$y[k4-2]
  }
 
  
  xxgrid <- sort(xxgrid, method="quick")
  
 
  
  yygrid <- sort(yygrid, method="quick")
## filter out repeated (possibly boundary repeated, and round)
## rounding does away with virtually same values, and problem of
## missing edges in cell=TRUE mode. What degree of rounding? 4 seems ok
## units of xxgrid and yygrid  are in unit square 0-1.
##browser()
xxgrid <- unique(round(xxgrid,4))
yygrid <- unique(round(yygrid,4))


  

  draw.srd.enc(rot_x,rot_y, rot_x_ws, rot_y_ws, E_calc, run_srd$q,rnames,run_srd$cstat,
               labelcell,colour,col,labelrectangle,   run_srd$x,run_srd$y,
               z, z_actual, ws_scaled, ws_mean,project,
               nowhitespace,cells,testYf,zname,new=new,xxgrid=xxgrid,yygrid=yygrid,
               testYf_actual,heading=title,athreshold=athreshold,resid=resid,spanr,
               nabs,absfreq,rabsfreq, nemp,emptyarea,ws_col,perm, borders,bgtextcol)
#--------------------------------------------------------------------------------------
##  add z axis (y axis) in 3D projection

if(!identical(project,c(0,0,0)) & !identical(project,c(180,180,180)) ){
##  if( project[3]!=180 | project[1]!=180 | project[2]!=180) 
##print(project)
 
  len_q <- 4*(q+1)
  ## origin of the z axis? Use rectangle coordinates of run_srd$x
  ## (usually 0,0 except if no whitespace.
  ## Curious: run_srd$ is retained from initail call with new=TRUE?
  ##  why?? 
  if(nowhitespace) len_q <- 4*q
  xmin <- min(run_srd$x[1:len_q], na.rm=TRUE)
 
  ymin <- min(run_srd$y[1:len_q], na.rm=TRUE)
 ##browser()
  xyz <- c(xmin,ymin,zmin_scaled)
  origin <- rotate(xyz,project)
  xyz<- c(xmin,ymin,zmax_scaled)
  
  end <- rotate(xyz,project)
  if((origin[1] != end[1]) & (origin[2] != end[2]) ) {
  arrows(origin[1],origin[2], end[1], end[2],length = 0.1, angle = 25,lwd=2,col=bgtextcol)
  }
  ##browser()
  # label z axis 
  text(origin[1]+0.8*(end[1]-origin[1]),origin[2] +0.8*(end[2]-origin[2]),
       zname,cex=0.75, c(1,-0.5) ,col=bgtextcol)
  text(end[1],end[2],paste(signif(zmax,4)),cex=0.75,c(0.5,0.5),col=bgtextcol)
  
  text(origin[1],origin[2],paste(signif(zmin,4)),cex=0.75,c(0.5,0.5),col=bgtextcol)
  ##browser()
}

fail <- FALSE
output <- list(fail=fail, outtable=outtable,outrect=outrect,notes=notes,seed_used=seed
  ,new=new,  E=E_calc, dactm= run_srd$dactm, chisqu=chisqu,pgof=pgof,
  x=run_srd$x, y=run_srd$y)
return(output)

}
## end of draw.srd
#########################################################################
##  Drawing and labelling function
##  note input x,y: 
##  x,y  rotated  projected coordinates of full rectangles  to mean value of zname
##  x_ws, y_ws rotated coordinates of back-projected rectangles onto whitespace
##  x_raw,y_raw unrotated x and y coordinates for no projection.
##  

draw.srd.enc <- function(x, y, x_ws, y_ws, e, q, rnames,cstat,labelcell,colour, col,labelrectangle,
    x_raw,y_raw,z, z_actual,ws_scaled,ws_mean, project,nowhitespace,cells, testYf,zname,new=TRUE,
    xxgrid=xxgrid,yygrid=yygrid,testYf_actual,heading=NULL,athreshold,resid,spanr,
    nabs,absfreq,rabsfreq,nemp,emptyarea,ws_col,perm, borders,bgtextcol) 
  { 
 
  
  if(borders) thicken <- 0.75
  blues <- col
  ## ask for prompt between each frame? 
  flat <- identical(project,c(0,0,0))
  
  if(identical(project,c(0,0,0)) | identical(project,c(5,5,5))){ask <- TRUE} 
  else 
  {ask <- FALSE}
  # backgcolour <- "#EEEEEE"
  # par(bg = backgcolour, lwd=1, ask=ask)
  # if(colour==3)par(bg = "#EFEFEF", lwd=1, ask=ask) #change background colour for "mono"
  #graphics.off()
  ## pauses required for relatively smooth pseudo-animation
  # any < 0.04 and you cannot see rotation steps. 
   if(ask) {Sys.sleep(0.04)}
   else
     {Sys.sleep(0.04)}
  
  devAskNewPage(ask=FALSE)
  par(ask=FALSE)
  ## set margins.
  par(mar=c(0, 0, 0, 0))
  plot.new() 
  ## new plotting area

  len_q <- 4*(q+1)
  if(nowhitespace) len_q <- 4*q
  xmin <- min(x[1:len_q], na.rm=TRUE)
  xmax <- max(x[1:len_q], na.rm=TRUE)
  ymin <- min(y[1:len_q], na.rm=TRUE)
  ymax <- max(y[1:len_q], na.rm=TRUE)
  
  xmin_raw <- min(x_raw[1:len_q], na.rm=TRUE)
  xmax_raw <- max(x_raw[1:len_q], na.rm=TRUE)
  ymin_raw <- min(y_raw[1:len_q], na.rm=TRUE)
  ymax_raw <- max(y_raw[1:len_q], na.rm=TRUE)
  ##browser()
  margx <- 0.4
  margx <- 0.2
  ##if(labelrectangle == 2){margx <- 0.8}  ## allow a 
  ## little more on right for arrow labels
  #hang on to plot.window  parameters as global vars
  
  xwin <- c(xmin-.3,xmax+margx)
  ywin <- c(ymin-.2,ymax+.2)
  
  xwin <- c(xmin-.2,xmax+margx)
  ywin <- c(ymin-.15,ymax+.15)
  
  
   plot.window(xlim=xwin,ylim=ywin)
 
##    plot.window(xlim=c(-0.15,1.0),ylim=c(0,1))
  subheading <- " "

  if(labelcell==1){
    subheading <- paste("Cell frequencies")
  }
  if(labelcell==2){
    subheading <- "Cell areas"
  }
  if(labelcell==3){
    subheading <- paste("Area/Freq")
  }
   if(labelcell==4){
    subheading <- paste("Cell error")
   }  
  if(labelcell==5){
    subheading <-  paste("Binary cell codes")
  }
  if(labelcell==6){

    y1y2 <- unlist( strsplit(zname,"/"))
    if(length(y1y2)==2){
      subheading <- paste0("mean(", y1y2[1],")/mean(",y1y2[2],")")
    }
    else
    {
    subheading <-paste0("Mean(", zname,")")}
  }  
  
  if(labelcell==7){
    subheading <- paste("Cell residuals")
  }  
 
  if(labelcell==8){
    subheading <- paste("Cell label: numerator sum of y" )
  } 
  if(labelcell==9){
    subheading <- paste("Cell label: denominator sum of y")
  } 
  

  m.list <- as.list(1:q)     ##set up list of matrices with x, y co-ordinates

#-----------------------------------------------------
##establish icell regardless of whether cell=T or F (repeated later!) 
icell <- as.numeric(vector(length=2^q))
icell <- NA
n_nonempty_cells <- nrow(testYf)
for (i in 1:n_nonempty_cells ){
  
  bin <- as.numeric(testYf[i,1:q])
  
  
  ibin <- sum(2^(which(bin==1)-1)) +1 
  ## keep record of  cell that ibin refers to 
  icell[ibin] <- i
} 
#------------------------------------------------------

# 
# blues11 <- c("#ffffd9b3","#eff9bdb3","#d5eeb3b3","#a9ddb7b3","#73c9bdb3","#45b4c2b3"
#  ,"#2897bfb3",  "#2073b2b3","#234ea0b3","#1c3185b3","#081d58b3")

# ## attempt to ensure continuous spectrum of shades for  spanr() partitions
# if(!is.null(spanr)){
#  if(spanr=="A") blues <- blues11[6:11]
# if(spanr=="!A") blues <- blues11[1:6]
# }
# 

##greens <-  c("#ffffccb3","#d9f0a3b3","#addd8eb3","#78c679b3","#31a354b3","#006837b3")
##mono <-    c("#f7f7f7b3","#d9d9d9b3","#bdbdbdb3","#969696b3","#636363b3","#252525b3")


## transparency codes
# 100% - FF
# 95% - F2
# 90% - E6
# 85% - D9
# 80% - CC
# 75% - BF
# 70% - B3
# 65% - A6
# 60% - 99
# 55% - 8C
# 50% - 80
# 45% - 73
# 40% - 66
# 35% - 59
# 30% - 4D
# 25% - 40
# 20% - 33
# 15% - 26
# 10% - 1A
##orangey <- c("#ffffd4b3","#fee391b3","#fec44fb3","#fe9929b3","#d95f0eb3","#993404b3")




max_z <- max( z_actual , na.rm=TRUE)
min_z <- min( z_actual, na.rm=TRUE )
margin <- 0.01 * ( max_z - min_z)
if( max_z == min_z ) {margin=0.001}
ncuts <- 5
width <-  (max_z-min_z - 2*margin)/(ncuts-1)
cuts <- min_z+margin +(seq(1:ncuts)-1)*width
## pretty cuts better??
# cuts <- pretty(cuts, n=4)
# ncuts  <- length(cuts)
# if(ncuts> 6){
#   cuts <- cuts[1:6]
#   ncuts <- 6
# }
## for  ncuts, there are ncuts + 1 bands
nbands <- ncuts + 1

mid_z <- (max_z+min_z)/2


## lays order for drawing rectangles as 1,2, 3, 
## which is in order of decreasing size of rectangle  
#o <- c(1:q)
## these are heights of projection axis.
##o <- unperm(perm,q)
## I think unperm is simply order(perm)!
o <- order(perm)
zm <- z[ seq(1,4*(q) ,by=4)]
if(any(zm != 0) ){
  ##overrides order with order of  size of projection
  o <- order(zm,decreasing=FALSE)
}
##o <- perm

##-------------------------------------------------------------------------------
##  analysis for cells
if(cells == TRUE)
  {
  

    if(colour==6){
    make_legend(blues,cuts,ncuts,zname,bgtextcol)}
    else
    {
## all white  colour? no point in legend
  
      allwhite <-  all(col[1:q]=="#ffffffb3") 
      if(!allwhite)make_col_legend(col[o],rnames[o],borders,thicken,bgtextcol)}

  
  n_nonempty_cells <- nrow(testYf)
  
## NOTE: this is not necessarily 2^q  
  xpos <- vector(length=2^q)
  ypos <- vector(length=2^q)
cellc <- character(length=2^q)
cellc <- NA
## BUG discovered that testYf$V1  is NOT proction mean but rescaled
##  need to replave testYf$V1 by testYf_actual
  max_cells <- max( testYf_actual[which(!is.infinite(testYf_actual))],na.rm=TRUE )
  min_cells <- min( testYf_actual,na.rm=TRUE )
  
  mid_cells <-  (max_cells +min_cells)/2
  margin <- 0.01 * ( max_cells - min_cells )
  if( max_cells == min_cells) {margin=0.001}
  nbands <- 6
  width <-  (max_cells-min_cells + 2*margin)/nbands

cellcol <- ws_col

## first loop over non-empty areas to get colour bands                                     
for (i in 1:n_nonempty_cells ){
  
  bin <- as.numeric(testYf[i,1:q])
  
  ibin <- sum(2^(which(bin==1)-1)) +1 
  ## keep record of  cell that ibin refers to 
  icell[ibin] <- i
 
# do not assign a colour band if Inf or NA
if(!is.na(testYf_actual[i]) & !is.infinite(testYf_actual[i]) & !is.null(zname) ) {

##band <- floor( (testYf_actual[i]-min_cells+margin)/width ) +1
band <- which_band(cuts,testYf_actual[i])

### BUG: testYf$V1 are not correct mean values but testYF_actual seem to be
## what IS testYf? 
## Ans: testYf is rescaled to ensure projection in plotting area. 
## need to use testYf_actual to get band? 


if(max_cells==min_cells) {delta <-  0.8}
else
  { delta <-  (1 - (1-0.8^n_nonempty_cells)*((testYf_actual[i]-min_cells) /
                                  (1.5*(max_cells - min_cells))) )}
   }
else
{
  
##no Y scale, so make white configuration
  band <- 0
  cellcol <- "#ffffff"

}

if(band>0){

##if( colour == 1 ) {cellcol <- greens[band]} ##earthy
# if( colour == 5 ) {cellcol <- orangey[band]}    #red_y
# if( colour == 6 ) {cellcol <- mono[band]}     ##mono_y
  
  ## cell colour=6 is intensity of Y 
if( colour == 6  ) {
##  browser()
##  make_legend(xmin,xmax,ymin,ymax,blues,cuts,ncuts,zname)
  cellcol <- blues[band]}     ##blue_y
}
# need to ensure that  this cell combination bin corresponds to
# ordering of cells in cstat[]. ibin is binary ordering in cstat[]
## colour=4 is "opaque"

if(colour==4 ){
  
  o <- c(1:q)
  
  zm <- z[ seq(1,4*q,by=4)]
  o <- order(zm,decreasing=FALSE)
  
  done <- FALSE
  iq <- 0
  for(k in 1:q){
    if(bin[k]==1 & !done) {
      iq <- o[k]
      done <- TRUE}
  }

cellcol <- blues[iq+1]
}



###use blended colours as for cells=FALSE

if(colour <= 5 ){

cellc[ibin] <- colblend(col,bin,q)
##cellc[ibin] <- col[ibin]
## why colblend? Loses transparency. test overriding

if(all(bin==0)) cellc[ibin] <- ws_col
}
else
{ ##browser()i.e. for blue_y, red_y, mono_y
  cellc[ibin] <- cellcol
}
## elements 1 and 2 of cstat are x,y coordinates of  cell midpoints for labelcells
##  transformed to rotated positions 


   Z<- testYf$V1[i]



xyz <- c(cstat[1+5*(ibin-1)],cstat[2+5*(ibin-1)],Z)
   test <- rotate(xyz,project)
   xpos[ibin] <- test[1]
   ypos[ibin] <- test[2]  

}  #for(i in 1:n_nonempty_cells



## revise plotting with lowest y first to ensure
## 3-D rotation  has correct layering
ot <- order(testYf$V1,decreasing=FALSE)
    testYf$V1 <- testYf$V1[ot] 
l <-length(ot)

##browser()
for(k in 1:l){
  
  ##  headbending code to  get out correct index of the icell!
  ##  vector!
i <- which(ot[k]== icell)
 
  #ignore if no whitespace  or if nowhitespace=TRUE?

  if(!(i==1 & (cstat[3]==0 | nowhitespace) ) ){
    bin <- rev(number2binary(i-1,q))
    if(cstat[4+5*(i-1)]>athreshold){
      testY <-  testYf$V1[k]
      ##browser()
      drawcell(x_raw,y_raw,q,project, testY,bin,colour=cellc[i],xx=xxgrid,yy=yygrid) 
    }
  }

} 
} ##if(cells=TRUE)
##--------------------------------------------------------------
else ##if(cells=TRUE) within draw.srd.enc
##------------------------------------------------------------
##  if rectangles (subgroups)  are to be displayed 

{  
  
  ## binary representation for  whitespace 
  bin <- vector(length=q)
  bin <- rep(0,times=q)
  


  ##  Draw the q rectangles
zip <- -1.0e20
wsdrawn <- FALSE
   for (ii in 1:q){ 
    i <- o[ii]
    m.list[[i]] <- as.matrix(cbind(x[(4*(i)-3):(4*i)], y[(4*(i)-3):(4*i)]))   
    ## extract components of x and y matrices into m.list 
    ##  which is list of matrices for each rectangle
  }


  if(colour == 6 ){
   make_legend(blues, cuts,ncuts, zname,bgtextcol)
  }
  else
  {## allsame colour? no point in legend
    allwhite <-  all(col[1:q]=="#ffffffb3") 
    # print(rnames[1:q])
    # 
    #       print()
    #       print(rnames[o])
    if(!allwhite) make_col_legend(col[o],rnames[o],borders,thicken,bgtextcol)}



  for (ii in 1:q){
   i <- o[ii] 
   zi <- z_actual[1+(i-1)*4]
## new approach to  problem of layering. Slot in whitespace 
##  where it is between levels of 3D axis


   if(nowhitespace==FALSE){
    ## browser()
   if(ws_mean< zi & ws_mean > zip){
     ## slot in whitespace
   
     if(colour <= 5 ){
       
      ##browser()
    drawcell(x_raw,y_raw,q,project, ws_scaled,bin, colour=ws_col,xx=xxgrid,yy=yygrid)
     }
     
     else
       
     {
  ##  fixing bug use actual z projection 
       za <- z_actual[4*(q+1)]
##       band <- floor( (z_actual[4*(q+1)]-min_z+margin)/width ) +1
       band <- which_band(cuts,za)
       if( colour == 6 ) {cellcol <- blues[band]}    #red_y
       
       drawcell(x_raw,y_raw,q,project, ws_scaled,bin, colour=cellcol,xx=xxgrid,yy=yygrid)
   }
     
     wsdrawn <- TRUE}
   }
##   browser()
## polygon to draw rectangles 
   ##experimental switch for ellipses

if(colour <= 5){  
 if(!borders) {
  polygon(m.list[[i]], col=col[i])
 }
else
  
{
  
##  ellipse(m.list[[i]], col=col[i])
  polygon(m.list[[i]], col=col[i], border=bgtextcol, lwd=1+ii*thicken)
  # xy <- m.list[[i]]
  # 
  # 
  # roundedRect(
  #   xleft=min(xy[,1]), ybottom=min(xy[,2]), xright=max(xy[,1]), ytop=max(xy[,2]), 
  #   rounding=0.05, col=col[i] ,lwd=ii*thicken)
    
  }
  
}
else

{    
  
  ##intensity shading. Choose band according to rectangle level.
  zi <- z_actual[1+(i-1)*4]
##   band <- floor( (zi-min_z+margin)/width ) +1
 
  band <- which_band(cuts,zi)
   

    if( colour == 6 ) {
      cellcol <- blues[band]
     }    
   
   polygon(m.list[[i]], col=cellcol)
    
    
  }
    zip <- zi  
   
  }  ## loop for(ii in 1:q)

##browser()
## exit loop and still not whitespace drawn. Ensure it is done
if(wsdrawn == FALSE & nowhitespace==FALSE ){
##if(nowhitespace==FALSE ){
  if(colour <= 5 ){
    ##put in whitespace

    
     drawcell(x_raw,y_raw,q,project, ws_scaled,bin, colour=ws_col,xx=xxgrid,yy=yygrid)
  
    
    if(borders){
 ## patch: need to redraw rectangles with bolden borders
## why?? 
    for (ii in 1:q){
    i <- o[ii]  
    #  ellipse(m.list[[i]], col=col[i])
  ##  polygon(m.list[[i]], col=col[i] ,border=bgtextcol,lwd=ii*thicken)
    }
    
  }
    
    }
  
  else  ##if(colour <-5)
    
  {
##    band <- floor( (z[4*(q+1)]-min_z+margin)/width ) +1
    za <- z_actual[4*(q+1)]
    ##       band <- floor( (z_actual[4*(q+1)]-min_z+margin)/width ) +1
    band <- which_band(cuts ,za)
    if( colour == 6 ) {cellcol <- blues[band]}    
    
    drawcell(x_raw,y_raw,q,project, ws_scaled,bin, colour=cellcol,xx=xxgrid,yy=yygrid)
  } 
  
 
  }

## colour=4 opaque
if(colour==4 ) {  
## put  rectangle boundaries  
for (ii in 1:q){
  i <- o[ii] 
  if(!borders){
  polygon(m.list[[i]], col=NA, lty=c("dashed"),lwd=1) 
  }
  else
  {

 ##   ellipse(m.list[[i]], col=NA, lty=c("dashed")) 
   polygon(m.list[[i]], col=NA, lty=c("dashed"),lwd=ii*thicken,border=bgtextcol) 
  ##  tvscreen(m.list[[i]], col=NA, lty=c("dashed"),lwd=ii*thicken) 
    # 
    # xy <- m.list[[i]]
    # roundedRect(
    #   xleft=min(xy[,1]), ybottom=min(xy[,2]), xright=max(xy[,1]), ytop=max(xy[,2]), 
    #   rounding=0.05, col=NA ,lty=c("dashed"),lwd=ii*thicken)
    
    }
  }
}


ncell <- 2^q

xpos <- vector(length=ncell)
ypos <- vector(length=ncell)

for (i in 1:ncell ){
  
  
  ## elements 1 and 2 of cstat are x,y coordinates of  cell midpoints
  ##  transformed to rotated positions 
  xm <- cstat[1+5*(i-1)]
  ym <- cstat[2+5*(i-1)]
  if(is.na(ym)==FALSE & is.na(xm)==FALSE){
 ##browser()
    Z <- NA
  for(j in 1:q){
    
    j4 <- 4*j
    maxx <- max(c(x_raw[j4-3],x_raw[j4]))
    minx <- min(c(x_raw[j4-3],x_raw[j4]))
    maxy <- max(c(y_raw[j4-3],y_raw[j4-2]))
    miny <- min(c(y_raw[j4-3],y_raw[j4-2]))
      if(xm < maxx & xm > minx ) {
      if(ym < maxy & ym > miny ) {
        
       Z <- zm[j] 
        
      }}
  
  }
  

   if(is.na(Z)) { Z <- ws_scaled}
  ##browser()
  xyz <- c(cstat[1+5*(i-1)],cstat[2+5*(i-1)],Z)
  test <- rotate(xyz,project)
  xpos[i]<- test[1]
  ypos[i]<- test[2] 
  
  }
} 

 
}  #else if(cells)
#------------------------------------------------------
  

  
##output cell labelcells 

#    cstat is 5 x 2^q matrix of cells statistics returned from .Fortran srd
#         1,2  are x,y coordinates of midpoint
#         3  is cell frequency    
#         4  is cell area of fitted configuration
#         5  is cell error   
## coerce cstat back to this  form (above it is a vector, (how confusing!!)
##  cstat<- matrix(cstat,nrow=5,ncol=64)



##first set area and frequency vectors, recompute frequency sum

ncell<- 2^q

##-------------------------------
## getting areas, freqs
sq <- 5*(seq(1:ncell)-1)
xcentre <- cstat[1+sq]
ycentre <- cstat[2+sq]
area <- cstat[4+sq]
##print("cstat area")
##print(round(area,3))
freq <- cstat[3+sq]
##print(freq)
##-------------------------------

tfreq <- sum(freq)
tarea <- sum(area)
if(nowhitespace) {tfreq  <- tfreq - freq[1]
                  tarea  <- tarea - area[1]}

diagerr <- freq/tfreq - area/tarea
residual <- (freq-area)/sqrt(area) 
##print(diagerr)
##  chi<- sum(r2,na.rm=TRUE)

##  pvalue <-  1- pchisq(chi,df-1)
## cat("\n","P value for area - frequency congruence =",pvalue)

txtcol <- bgtextcol
## or black??
txtcol <- "black"
##why red in cell mode? too hard to read
# if(cells & is.null(zname)==FALSE) {txtcol="red"}
# else {txtcol="black"}




  nodatacell <- 0
  diagError <- 0
#--Put on cell labels cell stats ---------------------------------------------------------  
#  but only on flat projection or if cells 
  ## why? Because cells not "seen" in 3D rotation
  ## later patch: preference to always show cell stats ???  
if(flat | cells | TRUE) {

for (i in (1+nowhitespace):2^q)

  {
  
  check <- NA
## if not flat and cells=TRUE, projection has no coordinates
  ## if cell frequency is zero and y undefined. So put condition
 
 ## patch: always show cell stats ???  
   if(TRUE | !(!flat & cells & freq[i]==0)) {

    
#browser()
  ##add binary codes    
    if(labelcell == 5){
##    bin <- as.numeric(testYf[icell[i],1:q])
    
    bin <- number2binary(i-1, q)
    check <- paste(bin[o],collapse='')
     }
    ## set area threshold for cell to be "visible"
   ## why?? better without
    
    if(area[i] > athreshold  | TRUE){
    
##  y and rr options only allowed with cells=T, for which testYf calculated
        

       if(is.na(icell[i])==FALSE){

       if(labelcell == 6) {
         ##  labelcell="y"
         check <- testYf_actual[icell[i]]

        }
   
       if(labelcell == 7) {
         ## labelcell="y/y0"

         check <- testYf_actual[icell[i]]/testYf_actual[icell[1]]
       }
       
 ##      bin <- as.numeric(testYf[icell[i],1:q])
       
       # if(labelcell == 5) {  #binary codes
       #   ## use  + . symbols?  Not so good
       # ##  binx <- ifelse(bin==1,"+",".")
       # 
       #   check <- paste(bin[o],collapse='')}
       #   
       if(labelcell == 8 | labelcell == 9) {check <- testYf_actual[icell[i]]}
       }
      
      
      #previously I had xp=xpos[i], but xpos not set for nodatacells
      #KAEON!! This is messy, but is a way of displaying cell stats of
      ## no data cells on flat projection
      if(identical(project,c(0,0,0))){
        xp <- xcentre[i]
## try mastubating y a little rnorm(1,cstat[2,i],0.01)
        ## unsatisfactory
        yp <-  ycentre[i] 
       
      }
      else
      {
      xp <- xpos[i]
      yp <- ypos[i]
      }
  
 ## omit if 
 ##browser()
 
 ## dont put a label if nowhitespace and i==1 and if are

 if(!((nowhitespace & i==1) | area[i] < athreshold )){

 ##only add cell labels if sufficient big, 
 ## to agree with not a "no data" cell excluded from outtable output


    ## add  Error% labelcells cstat[5
     ##  add cell frequency labelcells 
    if(labelcell == 1){  
      
      ## cstat index 4 is area, 3 is frequency
      text( col=txtcol,xp,yp, round(freq[i],1), adj = c(0.5,0),cex=0.75)
    }
    ##  add cell areas labelcells 
    if(labelcell == 2){  
      text(col=txtcol, xp,yp, round(area[i],1), adj = c(0.5,0), cex=0.75)
    }
    if(labelcell == 3){
  ##    text( col=txtcol,xp,yp, round(cstat[5,i],2), adj = c(0.5,0),cex=0.75)
## test code for ratio n/area
      
rats <- area[i]/freq[i]
text( col=txtcol,xp,yp, signif(rats,2), adj = c(0.5,0),cex=0.75)
## test code for ratio n/area
## text( col=txtcol,xp,yp, paste(round(cstat[3,i],0),"/",round(cstat[4,i],1),sep=""), adj = c(0.5,0),cex=0.75)  
    }
   ## add  residual labelcells 
   
   ## later: what is point of these "residuals". They are no Pearson type 
   ## residuals from indepenendece?
       if(labelcell == 7){
   ## res=(cstat[3,i]-cstat[4,i])/sqrt(max(cstat[4,i],0.5))
   # ## be honest  with 0's and Inf residulas
         
        text( col=txtcol,xp,yp, paste(round(residual[i],2)), adj = c(0.5,0),cex=0.75)
       }

  ## error 
   if(labelcell == 4){
 ##    print(diagerr[i])
     text( col=txtcol,xp,yp, paste(round(diagerr[i],3)), adj = c(0.5,0),cex=0.75)
   }
   
   if(labelcell == 5){
     check <- paste(bin[o],collapse='')
     text( col=txtcol,xp,yp, check, adj = c(0.5,0),cex=0.75)
   }
   
 
    if(labelcell == 6 ){
      text( col=txtcol,xp,yp, signif(check,3), adj = c(0.5,0),cex=0.75)
    }
 # if(labelcell == 7 ){  ##rr
 #   text( col=txtcol,xp,yp, round(check,2), adj = c(0.5,0),cex=0.75)
 # }
if(labelcell == 8 & is.na(check)==FALSE){  ##sumn
  text( col=txtcol,xp,yp, signif(check,2), adj = c(0.5,0),cex=0.75)
}
  if(labelcell == 9 & is.na(check)==FALSE){  ##sumd
    text( col=txtcol,xp,yp, round(check,1), adj = c(0.5,0),cex=0.75)
  }
##print(paste("i",i, "xp",signif(xp,3),"yp",signif(yp,3)) )
}  ##end if(!(nowhitespace & i==1))

 
}  ##cstat[4,i]>athreshold

}  ## end if(!(!flat & cells & cstat[3,i]==0))
}  ## end for (i in 1:ncell) loop ncell=2^q
}  ## end if(flat | cells ) 
#------------------------------------------------------------------

# add title: heading (at top) and subheading (bottom) 
  xmid <- (xmax+xmin)/2
  xmid <- (xwin[1]+xwin[2])/2
  text(xmid,ywin[2]+0.01*(ywin[2]-ywin[1]),heading,font=2,col=bgtextcol)
  diagError <- max(abs(diagerr[(1+nowhitespace):ncell]))
  yposn <- ywin[1]+0.03*(ywin[2]-ywin[1])
  text(xmid,yposn,paste("max(abs(error))=",signif(diagError,3)),col=bgtextcol,cex=0.8)
  yposn <- ywin[1]+0.055*(ywin[2]-ywin[1])
  if(nemp>0){
    if(nemp==1) {text(xmid,yposn,
                      paste0(nemp," empty area =",signif(emptyarea,2),"% of area" ),col=bgtextcol,cex=0.8 )
    }
    else
    {##plural!
      text(xmid,yposn,
           paste0( nemp," empty areas =",signif(emptyarea,2),"% of area"),col=bgtextcol,cex=0.8 )
    }   
    
  }
  if(nabs>0) {
    if(nemp>0) yposn <- ywin[1]+0.08*(ywin[2]-ywin[1]) 

    ## messy cos absfreq is % passed. to output n calculate
    if(nabs==1) {text(xmid,yposn,
                      paste0(nabs," missing cell ",round(absfreq,1), " observation. ",rabsfreq,"% of sample"),col=bgtextcol,cex=0.8 )
    }
    else
    {text(xmid,yposn,
                      paste0(nabs," missing cells ",round(absfreq,1), " observations. ", rabsfreq,"% of sample"),col=bgtextcol,cex=0.8 )
    }
  }

  text(xmid,ywin[1]-0.00*(ywin[2]-ywin[1]),subheading,col=bgtextcol,cex=0.9)

#---------------------------------------------------------------
## output  rectangle labelcells and arrows (pointers)

 
    ## use qtest instead of q to test whitespace

  if(is.null(spanr)){
    qtest <- q
  }
  else
  {
    ## for spanr() models label the whitespace
      qtest <- q  + 1
      if(spanr== "A")rnames[qtest] <- "!A"
      if(spanr=="!A")rnames[qtest] <- "A"
  }
if(labelrectangle > 0){
  
      
      fontsz <- 0.9  
      
      ## smaller font for  long rnames strings
  if( max(nchar(rnames))  > 25 )fontsz <- 0.7
  
  ##browser()
  if(labelrectangle == 1){

    if(!cells | (cells & flat)){ 
    # position label inside  rectangle
    for (i in 1:qtest){
 
      if(cells) {txtcol="blue"}
      else {txtcol="black"}
      
      txtcol <- "#2c7fb8"
     ## txtcol <- "#002D72"
## bright colours (colour==2) make white or (black?)labels? No better   
##   if(colour==2) txtcol <- "white"
 
## KAEON: there is a prob here with positioning of label. 
## sometimes outside rectangle if  ordering of x,y cordinates
## is not cyclical  0,1,2,3  left -up- right - down
## attempts to fix  based on raw coordiantes do not work
##  in the rotation 
## KAEON
      i4 <- 4*i

    ## random corners, 1, 2  3, 4 clockwise from bottom left
      corner <- floor(4*runif(1))+1
      
      rx <- range(x_raw[(i4-3):(i4-1)])
      xright <- rx[2]
      xleft  <- rx[1]
      
      ry <- range(y_raw[(i4-3):(i4-1)])
      ytop <- ry[2]
      ybot <- ry[1]
      
      if(corner == 1 ){Xlabel <- xleft+0.01
      Ylabel <- ybot +0.01
      adj <- c(0,0)}
      if(corner == 2 ){Xlabel <- xleft+0.01
      Ylabel <- ytop - 0.03
      adj <- c(0,0)}
      if(corner == 3 ){Xlabel <- xright-0.01
      Ylabel <- ytop-0.03
      adj <- c(1,0)}
      if(corner ==4  ){Xlabel <- xright-0.01
      Ylabel <- ybot +0.01
      adj <- c(1,0)}
      
  
  ## rotate positions , throWs OUT OF PLACE    
      
       xyz <- c(Xlabel,Ylabel,z[i4])
       test <- rotate(xyz,project)
       Xlabel <- test[1]
       Ylabel <- test[2] 
      #
     ##  if(!flat) adj <- c(0,0) ## k=3 should be left bottom coordinate marker
 text(col=bgtextcol, Xlabel, Ylabel, rnames[i], adj = adj,cex=fontsz,font=2)
    }
    }
    }

## label arrows on right or left 
  if(labelrectangle == 2 | labelrectangle == 3 ){  ##draw an arrow pointer
    
    if(!cells | (cells & flat)){ 
     yp <- vector(length=qtest)
    iq <- floor(qtest/2)
    for (i in 1:qtest){
      k <- 1
      
      if(i- round((i/2))*2==0) {
        k <- 0
       }
##     position of top (k=1) or bottom (k=0) rh corners      
      yp[i] <- y[4*i - k]

#       
    }
#   
    ot <- order(yp,decreasing=FALSE)
# 
#     
     for (ii in 1:qtest){
#  
      i<- ot[ii] 
      
      k <- 1
# ## odd or even number rectangle 

       if(i- round(i/2)*2==0) {
        k <- 0
       }
#      
      ypoit <- (ymin+ymax)/2 +((ii-iq)/10)*(ymax-ymin)
      
      if(labelrectangle == 2 ){
   ##  pointers on right
      arrows(xmax+0.04,ypoit, x[4*i-k], y[i*4-k],length = 0.1, angle = 25,col="#999999")##draw arrow
      text( xmax+0.04, ypoit, rnames[i], adj = c(0,1),cex=fontsz,col=bgtextcol,font=2)
      }
      
      if(labelrectangle == 3){
   ## pointers on left  k=0 and k=1 transform to lh 3 and 2
    ##   i.e  k <-  3 -k 
         k <- 3-k
      arrows(xwin[1]   + nchar(rnames[i])*0.01,ypoit, x[4*i-k], y[i*4-k],length = 0.1, angle = 25,col="#999999")##draw arrow
      text( xwin[1] , ypoit+0.01, rnames[i], adj = c(0,1),cex=fontsz,col=bgtextcol,font=2)
       }
     
    }
      
    
  }  
  } ##end if(labelrectangle==2 | labelrectangle==3)

  } #3if(labelrectangle > 0)
 
plot.window(xlim=xwin,ylim=ywin)

return("Configuration fitted")

}  ## end of function draw.srd.enc

####################################################################
rotate <- function(xyz,angles) { 
  
  if(identical(angles,c(0,0,0)) | identical(angles,c(180,180,180)) ){
    XD <- xyz[1]
    YD <- xyz[2]
    ZD <- xyz[3]
  }  else { 
    ## NB 0.017453292=pi/180 rad to angles
  xa <- angles[1]*0.017453292
  ya <- angles[2]*0.017453292
  za <- angles[3]*0.017453292
  
  X <- xyz[1]
  Y <- xyz[2]
  Z <- xyz[3]

COSya <- cos(ya)
SINya <- sin(ya)
COSxa <- cos(xa)
SINxa <- sin(xa)
COSza <- cos(za)
SINza <- sin(za)

XD <- X*(COSya*COSza) + Y*(-COSya*SINza) + Z*(SINya) 

YD  <- X*(SINxa*SINya*COSza+COSxa*SINza) + 
       Y*(-SINxa*SINya*SINza+COSxa*COSza) + 
       Z*(-SINxa*COSya) 

ZD <-  X*(-COSxa*SINya*COSza+SINxa*SINza) + 
       Y*(COSxa*SINya*SINza+SINxa*COSza) + 
       Z*(COSxa*COSya)
}
output <- c(XD,YD,ZD)


}  ##end rotate function
######################################################################
drawcell <- function(x,y,q,project,ws_scaled,bin,colour="#FFFFFF",xx,yy)
{  
   faceon <- identical(project,c(0,0,0)) | identical(project,c(180,180,180))
   

##  if(faceon | !is.infinite(ws_scaled) ){
  if(TRUE){
  bintest <- vector(length=q)
  nn <- 2*(q+1)
  xrot <- vector(length=4)
  yrot <- vector(length=4)
  nn1 <- nn-1 
  nnx <- length(xx)-1
  nny <- length(yy)-1

##  xx, yy are global variables of positions of the subcell grid

#---------------------------------------------------------------------------
## loop over  all nn1*nn1 subcells of the grid
## e.g q=6, nn=14, nn1=13 so loop over 13*13
for (k in 1:(nnx)){
  
  xxk <- xx[k]
  xxk1 <- xx[k+1]
  xm <- (xxk1 +xxk)*0.5

  
  for (l in 1:(nny)){
    yyl <- yy[l]
    yyl1 <- yy[l+1]  
    ym <- (yyl1 +yyl)*0.5  

 
  bintest <- rep(0,times=q)
  
   for(j in 1:q){
##browser()
     j4 <- 4*j
##     bintest[j] <- 0
     if(x[j4-3]<x[j4]) {
      maxx <- x[j4]
      minx <- x[j4-3]
     }
     else
     { 
       minx <- x[j4]
       maxx <- x[j4-3]   
     }
     if(xm < maxx   & xm > minx ) {
       
       if(y[j4-3 ]< y[j4-2]) {
         maxy <- y[j4-2]
         miny <- y[j4-3]
       }
       else
       { 
         miny <- y[j4-2]
         maxy <- y[j4-3]
       } 
  #     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
  if(ym < maxy   & ym > miny ) {
        bintest[j] <- 1
        
       }}
     
 
   }
  
  if(identical(bin,bintest) == TRUE) {
 
## filling a  subcell. Extract the unrotated subcell rectangle corners
    ## at z value corresponding 
  ## if( (abs(xxk-xxk1) > 0.0) &  (abs(yyl-yyl1) > 0) ){
    
      p1 <- c(xxk1,yyl1,ws_scaled)
      p2 <- c(xxk1,yyl,ws_scaled)
      p3 <- c(xxk,yyl,ws_scaled)
      p4 <- c(xxk,yyl1,ws_scaled)
##  rotate       
      t1 <- rotate(p1,project)
      t2 <- rotate(p2,project)
      t3 <- rotate(p3,project)
      t4 <- rotate(p4,project)
      
      xrot[1] <- t1[1]
      yrot[1] <- t1[2]
      
      xrot[2] <- t2[1]
      yrot[2] <- t2[2]
      
      xrot[3] <- t3[1]
      yrot[3] <- t3[2]
      
      xrot[4] <- t4[1]
      yrot[4] <- t4[2]
eps <- 0.01
 ## if(abs(yrot[1]-yrot[2])>eps & abs(xrot[1]-xrot[3])>eps) ??

polygon(xrot, yrot, border=NA,col=colour)    ## fill subcell colour   


## determine edges, depending on whether adjacent above, below, left, right

#----------------------------------------------------------------
##  on LEFT 
edge_left <- FALSE
##NOTE changed this from 1 to 2, is it right KANEON
if(k <= 1  ){
  # must be subcell on left boundary
  edge_left <- TRUE}

 
  else
  {
  #centre of cell on the left
  xm <- (xxk +xx[k-1])*0.5
  ym <- (yyl1 +yyl)*0.5
  bintest <- rep(0,times=q)
 for(j in 1:q){
   j4 <- 4*j

#  bintest[j] <- 0
if(x[j4-3]<x[j4]) {
  maxx <- x[j4]
  minx <- x[j4-3]
}
else
{ 
  minx <- x[j4]
  maxx <- x[j4-3]   
}
#     if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
#     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
if(xm < maxx   & xm > minx ) {
  
  if(y[j4-3 ]< y[j4-2]) {
    maxy <- y[j4-2]
    miny <- y[j4-3]
  }
  else
  { 
    miny <- y[j4-2]
    maxy <- y[j4-3]
  } 
  #     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
  if(ym < maxy   & ym > miny ) {
    bintest[j] <- 1
    
  }}
 } 
 
if(identical(bin,bintest)==FALSE){edge_left <- TRUE}


} 



#----------------------------------------------------------------
##  on RIGHT 

edge_right <- FALSE

if(k == nnx){
  # must be subcell on right boundary
  edge_right <- TRUE}
else
{xm <- (xx[k+2] +xxk1)*0.5 
 ym <- (yyl1 +yyl)*0.5 
 bintest <- rep(0,times=q)
 for(j in 1:q){
   j4 <- 4*j
   ##browser()
   #bintest[j] <- 0
   if(x[j4-3]<x[j4]) {
     maxx <- x[j4]
     minx <- x[j4-3]
   }
   else
   { 
     minx <- x[j4]
     maxx <- x[j4-3]   
   }
   #     if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
   #     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
   if(xm < maxx   & xm > minx ) {
     
     if(y[j4-3 ]< y[j4-2]) {
       maxy <- y[j4-2]
       miny <- y[j4-3]
     }
     else
     { 
       miny <- y[j4-2]
       maxy <- y[j4-3]
     } 
     #     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
     if(ym < maxy   & ym > miny ) {
       bintest[j] <- 1
       
     }}
 } 
 if(identical(bin,bintest)==FALSE){edge_right <- TRUE}
} 

#----------------------------------------------------------------
 ##  above 

 edge_above <- FALSE
 
 if(l == nny){
   # must be subcell on above boundary
   edge_above <- TRUE}
 else
 {xm <- (xxk1 +xxk)*0.5 
  ym <- (yy[l+2] +yyl1)*0.5 
  bintest <- rep(0,times=q)
  for(j in 1:q){
    j4 <- 4*j

 #   bintest[j] <- 0
 if(x[j4-3]<x[j4]) {
   maxx <- x[j4]
   minx <- x[j4-3]
 }
 else
 { 
   minx <- x[j4]
   maxx <- x[j4-3]   
 }
 #     if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
 #     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
 if(xm < maxx   & xm > minx ) {
   
   if(y[j4-3 ]< y[j4-2]) {
     maxy <- y[j4-2]
     miny <- y[j4-3]
   }
   else
   { 
     miny <- y[j4-2]
     maxy <- y[j4-3]
   } 
   #     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
   if(ym < maxy   & ym > miny ) {
     bintest[j] <- 1
     
   }}
  } 
  if(identical(bin,bintest)==FALSE){edge_above <- TRUE}
 }  
#----------------------------------------------------------------
##  below 
edge_below <- FALSE

if(l == 1){
  # must be subcell on above boundary
  edge_below <- TRUE}
else
{xm <- (xxk1 +xxk)*0.5
 ym <- (yyl +yy[l-1])*0.5
 bintest <- rep(0,times=q)
 for(j in 1:q){
   j4 <- 4*j

  #bintest[j] <- 0
  if(x[j4-3]<x[j4]) {
    maxx <- x[j4]
    minx <- x[j4-3]
  }
  else
  { 
    minx <- x[j4]
    maxx <- x[j4-3]   
  }
  #     if(xm < max(c(x[j4-3],x[j4]))   & xm > min(c(x[j4-3],x[j4]  )) ) {
  #     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
  if(xm < maxx   & xm > minx ) {
    
    if(y[j4-3 ]< y[j4-2]) {
      maxy <- y[j4-2]
      miny <- y[j4-3]
    }
    else
    { 
      miny <- y[j4-2]
      maxy <- y[j4-3]
    } 
    #     if(ym < max(c(y[j4-3],y[j4-2])) & ym > min(c(y[j4-3],y[j4-2])) ) {
    if(ym < maxy   & ym > miny ) {
      bintest[j] <- 1
      
    }}
 } 
 if(identical(bin,bintest)==FALSE){edge_below <- TRUE}
}  

## Need to now draw the edges 
##if(identical(bin,c(0,0,0,0,1,0) )  ) browser()

if(edge_right){
  
  p1 <- c(xxk1,yyl1,ws_scaled)
  t1 <- rotate(p1,project)
  p2 <- c(xxk1,yyl,ws_scaled)
  t2 <- rotate(p2,project)
  #draw edge by connecting line segments
  segments(t1[1],t1[2],t2[1],t2[2],col="black",lwd=1, lty=1)}


if(edge_left){
  
  p1 <- c(xxk,yyl1,ws_scaled)
  t1 <- rotate(p1,project)
  p2 <- c(xxk,yyl,ws_scaled)
  t2 <- rotate(p2,project)
  #draw edge by connecting line segments
  segments(t1[1],t1[2],t2[1],t2[2],col="black",lwd=1, lty=1)}
 


if(edge_above){
  
  p1 <- c(xxk1,yyl1,ws_scaled)
  t1 <- rotate(p1,project)
  p2 <- c(xxk,yyl1,ws_scaled)
  t2 <- rotate(p2,project)
  #draw edge by connecting line segments
  segments(t1[1],t1[2],t2[1],t2[2],col="black",lwd=1, lty=1)}



if(edge_below){
  
  p1 <- c(xxk1,yyl,ws_scaled)
  t1 <- rotate(p1,project)
  p2 <- c(xxk,yyl,ws_scaled)
  t2 <- rotate(p2,project)
  #draw edge by connecting line segments
  segments(t1[1],t1[2],t2[1],t2[2],col="black",lwd=1, lty=1)}

##}  #if(toosmalltodraw==F)


##if(identical(bin,c(0,0,0,0,1,0)) & l<4) browser()

} 

#if(identical(bin,bintest)



}}   #over k,l loop
##browser()
#---------------------------------------------------------------------------
if(identical(bin,rep(0,times=q)) ){
##   is occasional bug drawing left edge, attempt to fix
##   by always putting in unit square
p1 <- c(min(xx),min(yy),ws_scaled)
p2 <- c(min(xx),max(yy),ws_scaled)
p3 <- c(max(xx),max(yy),ws_scaled)
p4 <- c(max(xx),min(yy),ws_scaled)
##  rotate       
t1 <- rotate(p1,project)
t2 <- rotate(p2,project)
t3 <- rotate(p3,project)
t4 <- rotate(p4,project)

xrot[1] <- t1[1]
yrot[1] <- t1[2]

xrot[2] <- t2[1]
yrot[2] <- t2[2]

xrot[3] <- t3[1]
yrot[3] <- t3[2]

xrot[4] <- t4[1]
yrot[4] <- t4[2]
     
polygon(xrot, yrot, border="black", lwd=1, col=NA)    ## fill subcell colour   

}
#browser()
}
} ##end of drawcell
########################################################################

number2binary <- function(number, noBits) {
  ## inToBits  is from R.utils package
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}
########################################################################
binno <- function(n,q)
{
bin  <- vector(length=q)
bin <- as.numeric(bin)
#  to obtain binary representation of number n
#  in q 0/1 digits. Adapted from Fortran binno.
#  probably a better way, but number2binary above
#  give binary combination in reverse ordering
##browser()
pow <- 2^(q-1)
rem <- n
for (jj in 1:q){
  j <- q-jj +1
i <- floor(rem/pow)
if(i ==1) rem <- rem-pow
bin[j] <- i
pow <- floor(pow/2)
}
return(bin)
}
##############################################################################
#roxygen embedded tags  for srd.summary:
 #' To  summarise a n srd object
 #'
 #' @param fit An output object of  \code{srd}
 #' @author Roger Marshall <rj.marshall@@auckland.ac.nz>, The University of Auckland, New Zealand
 #' 
 #' @details A function to summarise the result of a run of 
 #'  \code{\link{srd}}.
 #' 
srd.summary <- function(fit){
  ## make a sumnmary of the fitted srd.  
 ##---------------------------------------------------------- 
  ## invoke generic summary() for other classes
  if(class(fit) != "srd"){
    message(paste(paste(substitute(fit)),"is not an srd object"))

    return(base::summary(fit))
  }
##------------------------------------------------------------
 
  write.table(fit$notes, quote=FALSE,row.names=FALSE,col.names=FALSE)
   
  cat("\n", "Rectangles R:")
  cat("\n")

  print(fit$rectangles,quote=FALSE,row.names=FALSE)
  
  cat("\n", "Cells C:")
  cat("\n")
  print(fit$cells,quote=FALSE,row.names=FALSE)
  cat("\n", "( + cell C within rectangle R)")
  
  cat("\n")
  
  if(fit$fail) {cat("\n","srd failed")}
   
}
##############################################################################
# Function  due to Billy Wu.
## sum duplicate rows
SumDupRows <- function(input, index, target, name){
  

  input[index] <- lapply(input[index], function(x) ifelse(is.na(x), 9999, x))
  
  eval(parse(text=sprintf("output <- aggregate(input$%s, by=lapply(input[index], unlist), sum)", target)))
  
  output[index] <- lapply(output[index], function(x) ifelse(x==9999, NA, x))
  
  names(output)[length(output)] <- name
  
  return(output)
}
##############################################################################
## colour blending function used for cells=TRUE to allow
## colour schemes to be plotted  by blended cell values

colblend <- function(col,bin,q){
## q=1, no blending required  
  if(q == 1){
    colblend <- col[1]
    return(colblend)
  }
 ## single  one colour cell, return unblended col 
  if(sum(bin) ==1 ){
   ione <- which(bin ==1 ) 
   colblend <- col[ione]
   return(colblend)
  }
  ##  whitespace , return input colour
  if(sum(bin) == 0){
    colblend <- col[1]
    return(colblend)
    
  }
##  must be more than one cell to blend. Start with first colour 
## try blending in others  
  ione <- which(bin ==1 ) 
  blend <- col2rgb(col[ione[1]],alpha=TRUE)
  nones <- length(ione)
 ## browser()
  for(j in 2:nones){
    
    rgba <- col2rgb(col[ione[j]],alpha=TRUE)
    a <- rgba[4]/255
    
    blend[1:3] <- blend[1:3]*(1-a) + rgba[1:3]*a
  }
  
  blend <- blend/255
  ## make blend[4] opacity arbitrary
  blend[4] <- 0.7
  colblend <- rgb(blend[1],blend[2],blend[3],blend[4])
  
  return(colblend)
 
}
###############################################################################
make_legend <- function(blues,cuts,ncuts,zname,bgtextcol){
## use par()$usr to get boundaries of plot (obviously!) 
  xmin <- par()$usr[1]
  xmax <- par()$usr[2]
  ymin <- par()$usr[3]
  ymax <- par()$usr[4]
  
 
  
  boxy <- (ymax-ymin)/50
  boxx <- (xmax-xmin)/25
  
  ypos <-  ymin +  seq(from=1, to=ncuts)*boxy*1.5
  xpos <- rep( xmin+(xmax-xmin)*0.02,times=ncuts)
  ##
  text(xpos[1],ypos[ncuts]+5*boxy,zname, adj = c(0),cex=0.9,col=bgtextcol)  
  text(xpos + 2* boxx,ypos+0.5*boxy, paste0("<",signif(cuts,3)),cex=0.9,col=bgtextcol )

  rect(xpos, ypos, xpos +boxx, ypos + boxy,border=bgtextcol,col=blues)
  ## add final band > largest cut
  xpos <- xpos[ncuts]
  ypos <- ypos[ncuts] + 1.5*boxy
  rect(xpos, ypos, xpos +boxx, ypos + boxy,border=bgtextcol,col=blues[ncuts+1])
  text(xpos + 2*boxx,ypos + 0.5*boxy,paste0(" ",signif(cuts[ncuts],3),"+"),cex=0.9,col=bgtextcol )

}  ## make_legend

##############################################################################
make_col_legend <- function(cols, rnames,borders, thicken,bgtextcol){
## insert a colour legend for each rectangle colour
## use par()$usr to get boundaries of plot (obviously!) 
## to position the legend  
  xmin <- par()$usr[1]
  xmax <- par()$usr[2]
  ymin <- par()$usr[3]
  ymax <- par()$usr[4]

   ncuts <- length(cols)
   
   
   boxy <- (ymax-ymin)/50
   boxx <- (xmax-xmin)/25
   
#  browser()
  ypos <-  ymin + seq(from=1, to=ncuts)*(boxy)*1.5
  xpos <- rep( xmin+(xmax-xmin)*0.02,times=ncuts)
  lw <- rep(1,times=ncuts)
  if(borders)lw <- seq(from=1, to=ncuts)*thicken
##   text(xpos[1],ypos[ncuts]+6*boxy,rnames, adj = c(0))  
  mchar <- max(nchar(rnames))
  tsize <- 0.9
  ## smaller font for long expressions
  if(mchar>20) tsize <- 0.8
  if(mchar>40) tsize <- 0.7
   text(xpos + 1.5* boxx,ypos+boxy,  paste(rnames), adj = c(0,1),cex=tsize,,col=bgtextcol )

  rect(xpos, ypos, xpos +boxx, ypos + boxy,border=bgtextcol,col=cols,lwd =lw)

  
}  ## make_col_legend
##############################################################################
which_band <- function(cuts,z){
  eps <- 1.0e-10
  w <- which(cuts < z+eps)
if(length(w)==0){
  band <- 1
}
else
{band <- max(w) + 1
}
 ## if(z>.99)browser()
return(band)
}
##############################################################################

 grandm <- function(d,f,q,x,y){
  ## mimic of fortran grandm() to establish start
  ## positions of random placed squares inside unit

   total <- sum(f)
   for(i in 1:q){
    
    totali <- sum(d[,i]*f)
    sidei <- sqrt(totali/total)
    u2 <- 100
    while(u2 > 1) {
     u1 <- runif(1)
     u2 <- u1 + sidei
   }

    x[1,i] <- u1
    y[1,i] <- u1
    
    x[2,i] <- u1
    y[2,i] <- u2
    
    x[3,i] <- u2
    y[3,i] <- u2
    
    x[4,i] <- u2
    y[4,i] <- u1
    
     
    
   }
   qq <- q +1 
## add boundary coordinates   
   x[1,qq] <- 0
   y[1,qq] <- 0
   
   x[2,qq] <- 0
   y[2,qq] <- 1
   
   x[3,qq] <- 1
   y[3,qq] <- 1
   
   x[4,qq] <- 1
   y[4,qq] <- 0
   
   return( list(x,y) )
 }  ## grandm()

 ################################################################################
 cellgrid <- function(x,y,q){
   ##superimpose grid over x and y 
   # 
   # x[which(x<0)]  <- 0
   # x[which(x>1)]  <- 1 
   # 
   # y[which(y<0)]  <- 0
   # y[which(y>1)]  <- 1 
   # x_grid <- sort(unique(x))
   # y_grid <- sort(unique(y))## form  surrounding box, 0,1 unit square, otherwise
   ##  expanded to edges of rectangle configuration
   # 
   x_grid <- unique( c(min(0,min(x)),  sort(x), max(1,max(x)) ) )
   y_grid <- unique( c(min(0,min(y)),  sort(y), max(1,max(y)) ) )
   
   ## assume minimum is 0, max 1. But what if not?
   ## force it to be so. 
   lxg <- length(x_grid)
   lyg <- length(y_grid)
   centreX <- (x_grid[-1] + x_grid[-lxg])/2
   lX <- (x_grid[-1] - x_grid[-lxg])
   centreY <- (y_grid[-1] + y_grid[-lyg])/2
   lY <- (y_grid[-1] - y_grid[-lyg])
   nxcells <- length(centreX)
   nycells <- length(centreY)
   
   
   area <- rep(0,times=2^q)

   max_sub_cellarea <- area
   cX <-  rep(NA,times=2^q)
   cY <-  rep(NA,times=2^q)
   ##browser()
   
   for(i in 1:nxcells){
     
     for(j in 1:nycells){
       
       XXX <- cell_comb(x,y,centreX[i],centreY[j],lX[i],lY[j],q)
       ## area combination number
       cellno <- XXX[1]
       ## returned cell numbers are,0,1,-, 2^q-1  for 2^q  possible cells
       ##  cell 0 is whitespace
       cellarea <- XXX[2]
       
       ## which of these sub-cells areas is largest?
       ## ( only needed to get positions for putting cell lables
       ##   so maybe inefficient)
       if(cellarea > max_sub_cellarea[cellno]){
 ##        print( paste( "i",i,"j",j, "max_area", signif(cellarea,5)) )
         max_sub_cellarea[cellno] <- cellarea
         cX[cellno] <- centreX[i]
         cY[cellno] <- centreY[j]
         
       }
       area[cellno] <- area[cellno] + XXX[2]
     }
     
   }  # for (i in 1:nxcell)
   # 
   # print(area)
   # print(sum(area))
   out <- list(area,cX,cY)
   return(out)
 }
 
 ##########################################################
 cell_comb <- function(x,y,centreX,centreY,lX,lY,q){
   ## remember x, y structured in blocks of 4 like this
   ## x1,x2  x1,x2, ......
   ## y1,y2, y1,y2,,....  
   
   pos <- 0
   cell <- 0
   for(iq in 1:q){
     e <- (pos + 1):(pos + 2)
     X <- x[e]
 #    x2 <- x[pos + 2]
  #   y1 <- y[pos + 1]
  #   y2 <- y[pos + 2]
     Y <- y[e]
     is.in.iq  <- between(X,centreX ) & between(Y,centreY)
     
     if(is.in.iq) cell <- cell + (2^(iq-1))
     
      pos <- pos +2
   }
   area <- lX*lY
   cellno <- cell + 1 
   out <- c(cellno,area)  
   return(out)
 }
 #####################################################################
 between <- function(X,centreX ){
   
   between <- (centreX > min(X)) & (centreX < max(X))
   return(between)
 }
 
# ################################################################ 

 fun_optim2 <- function(theta, freq,rfreq, crit){
   ## how will this work for whitespace=F?  Not sure
   ## adapting fun_optim, with 4 parameters per rectangle, to
   ## 3 parameters, forcing exact rectangle congruence
   p <- length(theta)
   q <- p/3 
   nobs <- sum(freq)
   ## first 2* are x cordinates, in pairs
   
   x <- theta[1:(2*q)]
   ## length of rectangle in x ddirection
   
   s2 <- 2*seq(1:q)
   
   lx <- x[s2]-x[s2-1]
   
   ##browser()
   zerox <- which(lx==0)
   if(length(zerox) >0 ) lx[zerox] <- 0.0001
   
   ## fixing length,  fix depth ly according to rectangle area in unit
   ## square assuned to be rfreq/nobs
   ly <- (rfreq/nobs) /lx
## get y position  form last q items of theta   
   y1 <- theta[(2*q+1):p]
   ## bind together with ly to get coordinates of rectangles 
   y2 <- y1+ly
   y <- c(rbind(y1, y2))
## first returned item is  area of cells in unit square
   area <- cellgrid(x,y,q)[[1]]*nobs
   ## print(area)
   deno <- area
   zeros <- which(area == 0 )
   if(length(zeros) > 0 ) deno[zeros] <- nobs*0.001
   ##browser()
   if(crit==1) func <- sqrt(sum( (freq-area)^2  ))
   if(crit==2) func <- sum( abs(freq-area)  )
   if(crit==3) func <- - sum( freq*log(deno)   ) + sum(freq)*log(sum(area))
   if(crit==4) {func <- sum( ((freq-area)^2) / deno  )
   }
   ## include  "stress" metric  
   if(crit==5)  func= sqrt(1- sum(freq*area)^2 / (sum(freq^2) * sum(area^2)))
   if(crit==6)  func= max(abs( area/sum(area) - freq/sum(freq) ) )
   
   ## missing cells :
   #browser()
   missed_cells <- which( area == 0 & freq >0)
   freq_missed_cells <- sum(freq[missed_cells])
   empty_cells <-  which(freq == 0 & area > nobs*0.001)
   area_empty_cells <- sum(area[empty_cells]) 
   ## adjust: scale up for missed cells 
   ##func <- func * (1+ length(missed_cells) ) )
   func <- func * (length(missed_cells)+length(empty_cells) +1)
   
  ## func <- func * (freq_missed_cells+area_empty_cells)
   
   return(func)
   
 }
 ################################################################ 
 chisqcheck <- function(freq,area){
   deno <- area
   zeros <- which(area <= 0.5)
   deno[zeros] <- 0.5
   ##browser()
   chisq <- sum( ((freq-area)^2) / deno  )
  ##print(paste("chisqchk=",signif(chisq,4))) 

 }
 ##############################################################
 
#  fit_srd <- function(d_, f_, theta){
#    q <- ncol(d_)
#    ## d_ is binary matrix of combinations and f_ is 
#    ## associated frequecy of each. 
#    ## NOTE: they do not have all  2^q combinations.
#    ncells <- 2^q
#    
#    ## what  cells are present? 
#    freq <- vector(length=ncells)
#    pow <- 2^(seq(1:q) -1)
#    ## get numbers of rows of d_
#    browser()
#    if(q >1){
#    present <- d_[,1:q]%*%pow +1
#    }
#    else
#    {
#      browser()
#      present <- d_  +1
#    }
#    freq[present] <- f_
#    cells <- seq(1:ncells)
#    empty <- setdiff(cells , present)
#    freq[empty] <- 0
# 
#    
#  ## freq[] is vector of frequencies of all 2^q cells
#   ## with index corresponding to  binary combination + 1   
#    
#   ## get marginal rectangle frequencies
#    rfreq <- vector(length=q)
#    
#    for (i in 1:q){
# ##     browser()
#      rfreq[i] <- sum(f_[which(  d_[,i]== 1 )])
#    }
# ## get starting position    
# ##   theta <- start_position(freq, rfreq )
#    crit <- 4
#    
#    ##initial value 
#    print(fun_optim(theta,freq,crit))
#    
#    ##chisqcheck(freq,area)
# ##   XXX <- nlm(fun_optim,theta,freq,crit)
# ##   print(XXX$minimum)
#    XXX <- optim(theta,fun_optim,gr=NULL,freq,crit,method="Nelder")
#    
#  #  browser()
#    print(XXX$value)
#   ## print(XXX$par)
# 
#    theta <- XXX$par
#    
#    p <- length(theta)
#    q <- p/4 
#    phalf <- p/2
#    x <- c(theta[1:phalf])
#    y <- c(theta[(phalf+1):p])
#    XXX <- cellgrid(x,y,q)
#    area <- XXX[[1]]*sum(freq)
#  
#  #  browser()
#    
#    ## try to simulate .Fortran output
#  
#    cstat <- vector(length=5*64)
#    ii <- 5*(seq(1,2^q)-1) 
#    print(XXX[[2]])
#    cstat[1+ii] <- XXX[[2]]
#    cstat[2+ii] <- XXX[[3]]
#    ## freq[i] is cstat[3+5*(i-1)]  for i=1,2^q
#    cstat[3+ ii]  <- freq
#    cstat[4+ ii]  <- area
#    cstat[5+ ii]  <- area-freq
# 
#    ## fashion X,Y as four corner coordnates
#    ## to include q+1 th whitespace
#    X <- vector(length= 4*(q+1))
#    Y <- vector(length= 4*(q+1))
#    ## must be better way than a loop???
#    for(i in 1:q){
#      i4 <- 4*(i-1)
#      X[1 +4*(i-1)] <- x[1+2*(i-1)]
#      X[2 +4*(i-1)] <- x[1+2*(i-1)]
#      X[3 +4*(i-1)] <- x[2+2*(i-1)]
#      X[4 +4*(i-1)] <- x[2+2*(i-1)]
#      
#      Y[1 +4*(i-1)] <- y[1+2*(i-1)]
#      Y[2 +4*(i-1)] <- y[2+2*(i-1)]
#      Y[3 +4*(i-1)] <- y[2+2*(i-1)]
#      Y[4 +4*(i-1)] <- y[1+2*(i-1)]
#      
#    }
#    X[4*q+1]  <- min(0,min(x))
#    X[4*q+2]  <- min(0,min(x))
#    X[4*q+3]  <- max(1,max(x))
#    X[4*q+4]  <- max(1,max(x))
#    
#    Y[4*q+1]  <- min(0,min(y))
#    Y[4*q+2]  <- max(1,max(y))
#    Y[4*q+3]  <- max(1,max(y))
#    Y[4*q+4]  <- min(0,min(y))
#    
#    
#   out <- list(cstat,X,Y)
#   return(out)
#  }  ##end fit_srd
 ###############################################################
 
 fit_srd2 <- function( d_ , f_ , theta , crit ){
   q <- ncol(d_)
   ## d_ is binary matrix of combinations and f_ is 
   ## associated frequency of each. 
   ## NOTE: they do not have all  2^q combinations.
   ## create  frequency "freq" that does cover all  2^q vcombination
   ## input theta is starting position of configuration
   ## with  3*p items. First 2*q are x1,x2 pairs, remaining q are
   ## y1 values. 
   
   ncells <- 2^q
   
   ## what  cells are present? 
   freq <- vector(length=ncells)
   pow <- 2^(seq(1:q) -1)
   ## get numbers of rows of d_
   if(q >1){
      present <- d_[,1:q]%*%pow +1}
   else
   {
      present <- d_ +1 
   }
   freq[present] <- f_
   cells <- seq(1:ncells)
   empty <- setdiff(cells , present)
   
## fill in  empty cells not represented in the data d_   
   freq[empty] <- 0
   ## freq[] is vector of frequencies of all 2^q cells
   ## with index corresponding to  binary combination + 1   
   
   ## get marginal rectangle frequencies
   # rfreq <- vector(length=q)
   # 
   # for (i in 1:q){
   #   ##     browser()
   #   rfreq[i] <- sum(f_[which(  d_[,i]== 1 )])
   # }
   # 
   
   ## can do this with apply()?  index 2 refers to col-wise
   ## yep, neat!
   rfreq <- apply(d_, 2 , function(x,frq=f_) {sum(frq[which(x == 1)])} )

   #-------------------------------------------------- 
   p <- length(theta)
   q <- p/3 
   ##extract  x pairs and y  from theta
   x <- theta[1:(2*q)]
   ## try purturbation ??   
##   x <- x + runif(2*q,min=0.01,max=0.01)
   
   s2 <- 2*seq(1:q)
   ## get lengths of x side of rectangle
   x1 <- x[s2-1]
   
   lx <- x[s2]-x1
   
   ## degenerate possibility  empty area
  ##browser()
   zerox <- which(lx==0)
   if(length(zerox) >0 ) lx[zerox] <- 0.0001

   ## get depths of y side of rectangle, fixing the rectangle area
   ly <- (rfreq/sum(freq)) /lx
   y1 <- theta[(2*q+1):p]

 ##  y1 <- y1 + runif(q,min=0.01,max=0.01)
   
   
   ## create final y2 coordinate, to fix area of rectangle.
   y2 <- y1+ly
   ##interleave y1 and y2 values to give y  pairs
   
   y <- c(rbind(y1, y2))
   ## get statistics of fitted configuation from cellgrid() 
   
   
   # 
   # x_grid <- unique( c(min(0,min(x)),  sort(x), max(1,max(x)) ) )
   # y_grid <- unique( c(min(0,min(y)),  sort(y), max(1,max(y)) ) )
   # 
   # browser()
   # 
   
   XXX <- cellgrid(x,y,q)
   area <- XXX[[1]]*sum(freq)
##   browser()
   # print("initial area")
   # 
   # print(signif(area,3))
   # 
   # 
##   print("initial crit")
 ## print(signif(   fun_optim2(theta,freq,rfreq,crit=6)  ,4 ) ) 
   
   ## chisqcheck(freq,area)
   if(F){
     
     ## using nlm() function 
   XXX <- nlm(fun_optim2,theta,freq,rfreq,crit,iterlim=100)
##   print(XXX$minimum)
   value <- XXX$minimum
   theta <- XXX$estimate
   }
 ##-------------------------------------------------- 

   if(T){
     
## try minqa:: powell algorithsm
## try overididing theta with new start position??     
 ##  theta <- starting_theta(freq,rfreq,crit)
 ##  print(theta)
   startv <- fun_optim2(theta,freq,rfreq,crit )
   print(paste("startv=",startv))
## recommends maxfun >10length(theta)^2, else get a warning message
     maxfun <- 10*length(theta)^2
   tol <- 1.0e-6
     XXX <- bobyqa(par=theta,fn=fun_optim2,
                lower = -10,upper = +10, control = list(maxfun=maxfun,rhoend=tol),
                freq,rfreq,crit)
     theta <- XXX[[1]]
     value <- XXX[[2]]
    print(paste("bobyqa tweak",value))
     feval <- XXX[[3]]
  #   print(paste("#eval =", feval))
 ##    print(theta)
   }
##-----------------------------------------------------   
   if(F){
     
   #  ngoes <- 0 
   #  conv <- F
   # while ( !conv ){
   #   ##suggested adaptive parameters?? 
  cont <- c(maxit=200, alpha=1, beta=1+2/p, gamma=0.75-0.5/p)
  cont <- c(maxit=200, alpha=1, beta=2, gamma=0.5)
  ##browser()
   XXX <- optim(par=theta,fn=fun_optim2,gr=NULL,freq,rfreq,crit,
          method="Nelder",control=cont)
   theta <- XXX$par
   value <- XXX$value
##   print(paste(value))
##   print(paste("convergence=",XXX$convergence))
   
  #  conv <- XXX$convergence  == 0 | ngoes ==1
  #  ngoes <- ngoes + 1
  # ## theta <- theta + runif(p,min=-0.01,max=0.01)
  #  }
   }
    ##-------------------------------------------------
## trying optimr:optimx)    
  #  lower <- rep( -Inf, times=p) 
  #  upper <- rep( Inf , times=p)
  #  XXX <- optimx(theta,fun_optim2,gr=NULL,
  #                method="Nelder-Mead"
  #             ,freq=freq, rfreq=rfreq, crit=crit)
  #  theta <- as.numeric(XXX[1,1:p])
  #  value <- as.numeric(XXX[1,(p+1)])
  # browser()
  ##----------------------------------------------------------
## supress next line  to "see" initial config
   # print("final theta")
   # print(signif(theta,3))
   # print("final crit")
   # print(XXX$value)
   
   
##-----------------------------------------------------
 # ## run again w optim? 
 #   theta <- XXX$par
 #   
  ## XXX <- optim(theta,fun_optim2,gr=NULL,freq,rfreq,crit,method="Nelder")
   #}
 #           
 # 
 #   #  browser()
 #   print(XXX$value)
 #   
 #   ## returned fitted theta  
  ## theta <- XXX$par
   ##--------------------------------------------------------   
   # ## run again?? with nlm()? 
   # XXX <- nlm(fun_optim2,theta,freq,rfreq,crit)
   # nobs <- sum(freq)
   # #  browser()
   # print(XXX$minimum)
   # 
   # ## returned fitted theta  
   # theta <- XXX$estimate
   # 
   
  
##-------------------------------------------------- 
   p <- length(theta)
   q <- p/3 
##extract  x pairs and y  from theta
   x <- theta[1:(2*q)]
   s <- seq(1:q)
   ## grt length of x side of rectangle
   lx <- x[2*s]-x[2*s-1]
   ##browser()
   zerox <- which(lx==0)
   if(length(zerox) >0 ) lx[zerox] <- 0.0001
   
   sum_freq <- sum(freq)
   ly <- (rfreq/sum_freq) /lx
   y1 <- theta[(2*q+1):p]
   ## create final y2 coordinate, to fix area of rectangle.
   y2 <- y1+ly
   ##interleave y1 and y2 values to give y  pairs
   
   y <- c(rbind(y1, y2))
## get statistics of fitted configuation from cellgrid() 
   
   XY <- cellgrid(x,y,q)
   area <- XY[[1]]*sum_freq
   
   #  browser()
   
   ## simulate .Fortran output vector cstat()
   
   cstat <- vector(length=5*64)
   ii <- 5*(seq(1,2^q)-1) 
   ##print(XXX[[2]])
   cstat[1+ ii] <- XY[[2]]
   cstat[2+ ii] <- XY[[3]]
   ## freq[i] is cstat[3+5*(i-1)]  for i=1,2^q
   cstat[3+ ii]  <- freq
   cstat[4+ ii]  <- area
   cstat[5+ ii]  <- area-freq
   ##print(paste("freq",freq))
   ## fashion X,Y as "four corner coordinates"
   ## to include q+1 th whitespace
   X <- vector(length= 4*(q+1))
   Y <- vector(length= 4*(q+1))
   # ## must be better way than a loop???
   # for(i in 1:q){
   #   i4 <- 4*(i-1)
   #   X[1 +4*(i-1)] <- x[1+2*(i-1)]
   #   X[2 +4*(i-1)] <- x[1+2*(i-1)]
   #   X[3 +4*(i-1)] <- x[2+2*(i-1)]
   #   X[4 +4*(i-1)] <- x[2+2*(i-1)]
   #   
   #   Y[1 +4*(i-1)] <- y[1+2*(i-1)]
   #   Y[2 +4*(i-1)] <- y[2+2*(i-1)]
   #   Y[3 +4*(i-1)] <- y[2+2*(i-1)]
   #   Y[4 +4*(i-1)] <- y[1+2*(i-1)]
   #   
   # }
   ## better done with seq()
   sqXY <- 4*(seq(1:q)-1)
   sqxy <- 2*(seq(1:q)-1)
   X[1+sqXY]  <- x[1+sqxy]
   X[2+sqXY]  <- x[1+sqxy]
   X[3+sqXY]  <- x[2+sqxy]
   X[4+sqXY]  <- x[2+sqxy]
   
   Y[1+sqXY]  <- y[1+sqxy]
   Y[2+sqXY]  <- y[2+sqxy]
   Y[3+sqXY]  <- y[2+sqxy]
   Y[4+sqXY]  <- y[1+sqxy]

      ## enclose unit square, or expanded to edges of rectangles
   q4 <- q*4 
   X[q4 +1]  <- min(0,min(x))
   X[q4 +2]  <- min(0,min(x))
   X[q4 +3]  <- max(1,max(x))
   X[q4 +4]  <- max(1,max(x))
   
   Y[q4 +1]  <- min(0,min(y))
   Y[q4 +2]  <- max(1,max(y))
   Y[q4 +3]  <- max(1,max(y))
   Y[q4 +4]  <- min(0,min(y))
     out <- list(cstat,X,Y,value , theta)
   return(out)
 }## fit_srd2
 ###############################################################
 
 start_position_nested <- function(freq, rfreq ){
   
   ## nested 
   n <- sum(freq)
   q <- length(rfreq)
  
   x <- vector(length=2*q)
   y <- vector(length=2*q)
   
   ## recale rfreq as proportions 
   rfreq <- rfreq/n
   
   o <- order(rfreq,decreasing=TRUE)
   
   
   for(j in 1:q){
     i <- o[j]
     w <- sqrt(rfreq[i])
     
     delta <- 0.05*runif(2)
     
     x[2*i-1] <- 0.5 + delta[1]-w/2
     
     x[2*i]  <- 0.5 + delta[1]+ w/2 
     
     
     y[2*i-1] <- 0.5 + delta[2]-w/2
     
     y[2*i]  <- 0.5 + delta[2] + w/2 
     
   }
   
    XY <- list(x,y)
   return(XY)
 }
 
 #######################################################
 
 start_position_random <- function(freq, rfreq ){
   
   ## nested 
   n <- sum(freq)
   q <- length(rfreq)
   
   x <- vector(length=2*q)
   y <- vector(length=2*q)
   
   ## recale rfreq as proportions 
   rfreq <- rfreq/n
   
  
   
   for(i in 1:q){

     w <- sqrt(rfreq[i])
     
     delta <- w*runif(2)
     
     x[2*i-1] <- delta[1]
     
     x[2*i]  <- delta[1]+ w 
     
     
     y[2*i-1] <- delta[2]
     
     y[2*i]  <- delta[2] + w 
     
   }

   XY <- list(x,y)
   return(XY)
 }
 ##########################################################
 start_position_indep <- function(freq, rfreq ){
   
   ## nested 
   n <- sum(freq)
   q <- length(rfreq)
   
   x <- vector(length=2*q)
   y <- vector(length=2*q)
   
   ## recale rfreq as proportions 
   p <- rfreq/n
   
   if(q ==4){
     
     x[1] <- 0
     x[2] <- p[1]
     y[1] <- 0
     y[2] <- 1
     
     x[3] <- p[1]*(1-p[2])
     x[4] <- x[3] + p[2]
     y[3] <- 0
     y[4] <- 1
     
     x[5] <- 0
     x[6] <- 1
     y[5] <- 0
     y[6] <- p[3]
     
     x[7] <- 0
     x[8] <- 1
     y[7] <- p[3]*(1-p[4])
     y[8] <- y[7] + p[4]
   }

   XY <- list(x,y)
   return(XY)
 }
 
 #######################################################
 
 start_position_twoway <- function(freq, rfreq ){
   
   ## nested 
   n <- sum(freq)
   rfreq <- rfreq/n
   freq <- freq/n
   
   q <- length(rfreq)
   
   x <- vector(length=2*q)
   y <- vector(length=2*q)
   o <- order(rfreq,decreasing=TRUE)
 
##  largest square   
   i <- o[1]
   w <- sqrt(rfreq[i])
   
  
   x1 <- 0.5 -w/2
   x[2*i-1] <- x1
    x2 <- 0.5 + w/2 
    x[2*i]  <- x2
   
   y1 <- 0.5 -w/2
   y[2*i-1] <- y1
   
   y2 <- 0.5  + w/2 
   
   y[2*i]  <- y2
   
   ## recale rfreq as proportions 
  ## set up binary matrix?
   binm <- matrix(nrow=2^q,ncol=q)
   for (b in 1:2^q){
    binm[b,] <- number2binary(b-1,q) 
   }
   

   for(j in 2:q){
     ij <- o[j]
     ## size of overlap with square 1?
     pair <- c(i,ij)
     overlap <- sum(freq[which(  rowSums(binm[,pair])==2  )])
     w <- sqrt(rfreq[i])
     d <- overlap/w +0.01*runif(1)
      if(ij==2){

     x[2*ij-1] <- x2 - d
     x[2*ij]  <- x2 - d + w
     y[2*ij-1] <- 0.5 - w/2
     y[2*ij]  <- 0.5  + w/2 
    }
   
     
     if(ij==3){

       y[2*ij-1] <- y1 +d  -w
       y[2*ij]  <- y1 + d 
       x[2*ij-1] <- 0.5 - w/2
       x[2*ij]  <- 0.5  + w/2 
     }
     
     
     if(ij==4){ 

     x[2*ij-1] <- x1 - d
     x[2*ij]  <-  x1 - d  +w
     y[2*ij-1] <- 0.5 - w/2
     y[2*ij]  <- 0.5  + w/2 
     }
    
     
     if(ij==5){

       y[2*ij-1] <- y2 - d
       y[2*ij]  <-  y2 - d + w
       x[2*ij-1] <- 0.5 - w/2
       x[2*ij]  <- 0.5  + w/2 
     }
     
     }
     
   XY <- list(x,y)
   return(XY)
 }
 
 ##########################################################################
 #roxygen embedded tags  for srd.summary:
 #' To  plot an srd drawing object
 #'
 #' @param fit An output object of  \code{srd}
 #' @author Roger Marshall <rj.marshall@@auckland.ac.nz>, The University of Auckland, New Zealand
 #' 
 #' @details A function to plot an  srd diagram after creating with srd() 
 #'  \code{\link{srd}}.
 #' 
 srd.plot <- function(fit){
   
   if(class(fit) != "srd"){
     message(paste(paste(substitute(fit)),"is not an srd object")) 
     return( graphics::plot(fit) )

   }
   
   replayPlot(fit$image)
   
   
 }
 ############################################################################
 ##expected counts under independence
 
 ##  all 2^q combinations
 expectedf <- function(d_,f_){
   n <- sum(f_)
   p <- colSums(d_*f_)/n
   q <- ncol(d_)
   not_p <- 1 - p

   ef <- vector(length=(2^q) ) # initial vector
   of <- vector(length=(2^q) ) # initial vector
   ##browser()

   edata <- matrix(ncol=q,nrow=(2^q) )
   for (i in 1:(2^q) ){
     
    combo <-  rev( number2binary(i-1,q))
    edata[i,] <- combo
    ef[i] <- (   prod( combo*p + (1-combo)*not_p)  ) *n   
    ##observed frequencies
    of[i] <- 0
    names(combo) <- colnames(d_)
    for(j in 1:length(f_)){
     if(identical(combo,d_[j,])) of[i] <- f_[j]
    }
  
   }
  ## print(cbind(edata,of,ef))
   chisq <- sum((of-ef)^2/ef)
  
   

   pind <- 1 - pchisq(chisq, df=1)
   message("complete independence chi_sq=",signif(chisq,3)," P= ", signif(pind,5))
   
  out <- list(edata,ef)
  
##  print(cbind(edata,ef))
 # browser()
  return(out)
 }
 
 ######################################################################
 
 
 is.wholenumber <-
   function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
 
 ##################################################################
 #roxygen embedded tags  for reweight:
 #' To  plot an srd drawing object
 #'
 #' @param x binary/logical variable
 #' @param weight existing weights
 #' @param alpha  forced probability for \eqn{x=1} for new weights
 #' @author Roger Marshall <rj.marshall@@auckland.ac.nz>, The University of Auckland, New Zealand
 #' 
 #' @details A function to create weights to enable srd() to create diagram with forced
 #' proportion of \eqn{x=1}
 #'  \code{\link{srd}}.
 #' 
 
 reweight <- function(x,weight=NULL,alpha=0.5){
   
 ## to form weights to ensure  pseudo counts
 ## such that, for binary x,  P(x=1) = alpha
 ##  f_ is a pre-existing weighting
  
   if(is.null(weight)){
     px <- mean(x,na.rm=TRUE )
   }
   else
     
   {px <- sum(x*weight,na.rm=TRUE)/sum(weight,na.rm=TRUE)}
   
   
   QQQ <- alpha*x + (1-alpha)*(1-x) 
   PPP <- px*x + (1-px)*(1-x) 
   if(is.null(weight)){  
     reweight <- QQQ/PPP 
     }
   else
     {  reweight <- weight*QQQ/PPP   }
   return(reweight)
 }  ## end reweight
 ####################################################
 
 mycol2hex <- function(col){
  ## convert to rgb, then to hex. 
   ## seems ok if col already  hex
   X <- as.character(as.hexmode( col2rgb( col  )))
   
   ## X is long vector of  blocks of 3 R, G, B elements
   hex <- NULL
   q <- length(col)
##   print(X)
   for (i in 1:q){
     i0 <- (i-1)*3

     hexi <- paste0("#", X[i0 +1], X[i0 + 2], X[i0 + 3] )
  ##   print(hexi)
     hex <- c(hex,hexi)
   }
   
   return(hex)
 }
 
 ##################################################################
 
 starting_theta <- function(freq,rfreq,crit){
   ## play with starting position  (overidding input theta)  
   # 
   #   # 
   Dopt <- Inf
   for(config in 2:2){
     #----------------------------------------------------------
  if(config == 1)  { 
  XY <- start_position_nested(freq, rfreq )
  q <- length(rfreq)
   #   ##browser()
       x <- XY[[1]]  ##  x1,x2  pairs
       y <- XY[[2]]
       s <- 2*seq(1:q) -1  # odd numbers
   # ##---------------------------------------------------   
    theta_s <- c(x,y[s])
  D <- fun_optim2(theta_s, freq,rfreq, crit)
  print(paste("nested", D))
  if(D < Dopt){
    Dopt <- D
    theta <- theta_s
  }
  }  ## s==1
  #-----------------------------------------------------------
     if(config == 2)  { 
       XY <- start_position_random(freq, rfreq )
       q <- length(rfreq)
       #   ##browser()
       x <- XY[[1]]  ##  x1,x2  pairs
       y <- XY[[2]]
       s <- 2*seq(1:q) -1  # odd numbers
       # ##---------------------------------------------------   
       theta_s <- c(x,y[s])
       D <- fun_optim2(theta_s, freq,rfreq, crit)
       print(paste("random",D))
       if(D < Dopt){
         Dopt <- D
         theta <- theta_s
       }
     }  ## s==1
     
  #------------------------------------------------------------   
     if(config == 3)  { 
       XY <- start_position_indep(freq, rfreq )
       q <- length(rfreq)
       #   ##browser()
       x <- XY[[1]]  ##  x1,x2  pairs
       y <- XY[[2]]
       s <- 2*seq(1:q) -1  # odd numbers
       # ##---------------------------------------------------   
       theta_s <- c(x,y[s])
       D <- fun_optim2(theta_s, freq,rfreq, crit)
       print(paste("indep", D))
       if(D < Dopt){
         Dopt <- D
         theta <- theta_s
       }
     }
   }  ## for(config)
  
  return(theta)
  }
 