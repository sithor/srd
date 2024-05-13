#roxygen embedded tags to generate srd.Rd manual file
#' Draws a scaled rectangle diagram.
#'
#' @param  data  A data frame or a table
#' @param rectangles  A character vector to represent the \eqn{k} rectangles. 
#' Each element is a valid logical expression that defines a rectangle (see examples). 
#' 
#' @param criterion Criterion to determine fitting optimisation criterion.   
#' Use \code{"lsqs"} for least squares, \code{"labs"} for least absolute difference, \code{"logl"} for minus log-likelihood, 
#' \code{"chisq"} for chi-square, \code{"stress"} for constrained least squares.
#' 
#' @param z The name of variable(s) to form a 3D projection axis. If a single variable
#' the mean value of the variable is projected, if  two variables the ratio of their means
#' is projected. If \code{NA}, the sequence number of rectangles is assigned. 
#'  
#' @param weight A string giving name of a frequency weight of each row of \code{data}. 
#' @param  new If \code{TRUE}  a new configuration is generated. \code{new=FALSE}  assumes 
#' a configuration from a previous matching \code{srd()} call already exists in the workspace.
#' @param seed Random number seed. \code{NA}  for non-repeatable, possibly different configurations. 
#' @param thin A penalty parameter to avoid thin and elongated rectangles. The \code{criterion} is penalized 
#' by \eqn{e^thin} where \eqn{e} maximal length/breadth of the  rectangles.  
#' @param whitespace If FALSE will not display whitespace, that is, 
#' the space that is not in at least one of the rectangles specified in \code{rectangles}.
#' @param title A title written on the srd graphic.
#' @param layer Determines layering of rectangles from back to front. If TRUE smallest rectangle is at front
#' largest at back. If FALSE, layered as the order listed in \code{rectangles} parameter with first at back last at front.
#' layer is ignored if \code{z} is non-NULL, when layering is by value of \code{z} (largest on top). 
#' @param clickable \code{TRUE} if the graphic is active for on-the-fly mouse input (see Details).
#' @param reweight a character string that matches one of \code{rectangles} 
#' to force analysis of pseudo counts in which the \code{reweight} rectangle
#' is  forced relative frequency 0.5. 
#' @param ...  Other control parameters, for example user defined
#'  rectangle colours (see details).   
#'  
# @import survival
# @import grDevices
# @import graphics
# @import stats
# @import randomcoloR
# @import utils 
# @importFrom grDevices devAskNewPage rgb
# @importFrom graphics arrows par plot.new 
# @importFrom plot.window polygon segments text title
# @importFrom stats complete.cases pchisq runif
#'
#' @return A summary of the graphic.  It gives rectangle and cell 
#' statistics and computational notes. Use
#' \code{\link{srd.summary}} to display. 
#' @author Roger Marshall, <rj.marshall@@auckland.ac.nz>, The University of Auckland, New Zealand
#' @seealso  \code{\link{srd.summary}}  for summary statistics of an \code{srd} object
#'
# @export 
#' @details A function to visualise  \eqn{k} attributes as a Venn-like  
#' diagram using \eqn{k}  rectangles
#' with areas scaled, as best as possible, to be proportional to  frequency. 
#' It allows \eqn{k \le 6}.  
#' Fitting is done by optimising the congruence between cell 
#' area \eqn{a_i} and cell frequency  \eqn{ n_i} 
#' by a chosen  criterion.
#' A starting configuration is fixed by intelligent or
#'  random positioning. A graphic is produced that is active for mouse 
#'  input if \code{clickable=TRUE} with a temporary diaglog box.
#'  Dialog options include: \code{Refit} which toggles fitting criterion to 
#'  usually created a different configuration;  \code{colour} which toggles 
#'  different colour shading arrangements, one of which is random colours (for user specified 
#'  colours see \code{...} parameters);  \code{Mode} which toggles "cell" 
#'  and "rectangle" mode (the distinction is that, 
#'  in an oblique projection, either cells or rectangles are 
#'  rendered); \code{Cell stats} superimposes statistics of 
#'  cells; \code{Labels} toggles how labels to rectangles are drawn;  \code{Rotate} 
#'   throws the diagram into  a 3D projection allowing rotation steps,
#'  the third axis being \code{z}.  
#'  Repeated clicks turns the diagram. This is more quickly 
#'  done with \code{Toggle view} which toggles an oblique, 
#'  flat or full side-on view.  \code{Esc}  leaves the mouse-active graphic. 
#'  
#' Additional ellipsis \code{...} parameters may include 
#'  a colour pallette for rectangles  \code{col}. 
#' Each  colour may include a transparency code, otherwise rectangles are opaque.
#' If a non-white colour for whitespace is preferred it can be specified with \code{col0}.
#' Background colour is specified with \code{bg}. 
#' @references 
#' Marshall, RJ. Displaying  categorical data relationships by scaled rectangle diagrams. Statistics in Medicine 20, 1077-1088, 2001.
#' Marshall, RJ. Scaled rectangle diagrams can be used to visualise clinical and epidemiological data. 
#' Journal of Clinical Epidemiology 58, 974-981, 2005.
#' Marshall, RJ. Determining and visualising at-risk groups in case-control data. 
#' Journal of Epidemiology and Biostatistics 6,  343-348, 2001

#' @examples 
#' # Simulate  and show independent binary A, B, C 
#' A <- rbinom(1000,1,0.1)
#' B <- rbinom(1000,1,0.4)
#' C <- rbinom(1000,1,0.3)
#' x <- cbind(A,B,C)
#' srd(data=x,rectangles=c("A","B","C"))
#' 
#' # Titanic data
#' library(datasets)
#'  Titan <- srd(data=Titanic,
#'      rectangles=c("Sex=='Female'",
#'                   "Age=='Child'",
#'                   "Class=='1st'"),
#'        weight="Freq", title="Titanic survival",z="Survived")
#' srd.summary(Titan) 
#' 
#' ## Co-occurrence of Scabies  and Impetigo in Fiji 
#' scab_imp <- data.frame( "Scabies"=c(1,1,0,0),"Impetigo"=c(1,0,1,0),"n"=c(2021,543,112,8211))
#' srd(rectangles=c("Scabies","Impetigo"),weight="n",
#'            clickable=TRUE,scab_imp, criterion="chisq",
#'            col=c("#ff0012b3","#765388b3") )
#'  
srd <- function( data, rectangles=NULL, 
      weight=NULL, new=TRUE, 
      seed=NA, whitespace=TRUE,criterion="labs",thin=0,z=NULL,
      title=NULL,layer=TRUE, clickable=TRUE,reweight=NULL, ...)
{

## start the accumulation of "notes" on the analysis 
notes <- paste("Using srd version 1.1 Run time: ",Sys.time())
##default title is data argument
if(is.null(title)) title <- paste(substitute(data),"data")

titleX <- title
## is data a data.frame? e.g could be "table". Coerce to be
if(class(data)[1] != "data.frame"){
  data <- as.data.frame(data)
}
## defaults
 cells <- FALSE 
 ws_col <- "white" 
 user_ws_col <- "white"
 c0 <- "white" 
 spanr <- NULL
 ## set no fill colours (white)

 borders <- FALSE
 nofillX <- rep("#ffffffb3",times=6)
 defaultcol <- as.vector(c("#58D6B880","#A7C8CC80","#A9D85E80","#6B785B80","#D8AD7580","#36648B80"))
 ## default col,  colour blind pallette cbp1 RBrewerColor
 cbp1 <- c("#99999980", "#E69F0080", "#56B4E980", "#009E7380",
           "#F0E44280", "#0072B280", "#D55E0080", "#CC79A780")
 #palette using grey
 cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  defaultcol <- cbp1
 ##paste transparemcy 
 defaultcol <- paste0(substr(x=defaultcol, start=1,stop=7),"90") 
## print(defaultcol)
 nofillc <- nofillX  ## no fill colour, may be overrifden by user col col= 
 bg <- "#ffffff"
 ##----------------------------------------------------------
 ## deal with any ellipsis ... three dot parameters
 
 if(length(list(...)) !=0 ){
   ellips <- list(...)
  if("col" %in% names(ellips) ){
    ## ensure transparency ??
    ## convert colour names (if any) to hex
    
    ## whichhash <- substr(x=ellips$col, start=1,stop=1)
    ##requires gplots package
    ##ellips$col <- gplots::col2hex(ellips$col)
    
    ellips$col <- mycol2hex(ellips$col)
    ## paste in transparency code
    ellips$col <- paste0(substr(x=ellips$col, start=1,stop=7),"90") 

   }
  ## test for allowed additional ellipsis parameters
  allowed <- c("spanr","col0","col","bg")
  allow <- is.element(names(ellips),allowed)

  if(!all(allow)){
    
 message(paste("Note: there are non-allowed  ...  additional parameter(s):"))
 message(paste( names(ellips)[which(allow==FALSE)], collapse=" " ))
    
  }

         V <- is.element(names(ellips),"spanr")
        if(length(which(V==TRUE))>0 ) {
         spanr <- as.character(ellips[which(V==TRUE)] )
        }

        V <- is.element(names(ellips),"col0")
        if(length(which(V==TRUE))>0 ) {
          c0 <- as.character(ellips[which(V==TRUE)] )
          user_ws_col <- c0

        }
        
        V <- is.element(names(ellips),"col")
        if(length(which(V==TRUE))>0 ) {
          nofillc <- ellips$col
        }
    
        V <- is.element(names(ellips),"bg")
        if(length(which(V==TRUE))>0 ) {
        bg <- as.character(ellips[which(V==TRUE)] )
        }
        
        
        # V <- is.element(names(ellips),"opaque")
        # if(length(which(V==TRUE))>0 ) {
        #   opaque <- as.logical(ellips[which(V==TRUE)] )
        #   print(opaque)
        # }
 }  ##if(length(list(...)) !=0 )
 ##----------------------------------------------------------
 usercols <- !identical(nofillc,nofillX) | !identical(c0,ws_col)
 ## no fill parameters may be overridden by user colours
 # if(opaque){
 #   nofillc <- paste0(substr(x=nofillc, start=1,stop=7),"FF") 
 #   defaultcol <-  paste0(substr(x=defaultcol, start=1,stop=7),"FF") 
 #   print(defaultcol)
 # }
 ## data input may be  a table?  Convert to data frame.
 if(class(data)[1] !="data.frame"){
   data <- data.frame(data)}
 

start_angle <- c(0,0,0)
##if(dev.cur() > 1) dev.off() ## do I need this? clears screen
  ##=================================================================
  ##  pre-analysis to determine binary data matrix
  ##=================================================================

 if(is.null(rectangles)){
   return(message("rectangles parameter must be specified"))
 }
 ## NTM previously I had "w" as weight argument. For consistency 
 ## with other R packages use "weight" as name of weighting varaible
 ##  Here copy to w to proceed  for name of weighting variable
  ## and make the actual weight values  the matrix "weight"
 w <- weight
  
  if(is.null(w) ){weight <- matrix(data=1, ncol=1,nrow=nrow(data))

    } else {

      weight <- data[,w]
      
    make_0 <- which(is.na(weight) | weight <0 )
if(length(make_0) >0 ){
weight[make_0] <- 0
 
message("Negative or NA weight made zero")}
    
## non-integer weight?
    wtest <- weight - round(weight,0)
    if(any(wtest !=0)){
      message("non-integer weights")
  ##    weight <- round(weight,0)
    }
    
 
    }
 

##repetition in values of rectangles?
 repeats <- NULL
  q <- length(rectangles)
## check for repetitions, removing spaces from rectangles
  rectangles <- gsub(rectangles,pattern=" ", replacement="")
  rectangles <- unique(rectangles)
  
  qr <- length(rectangles)
  if(qr < q){
  message("Repetition in rectangles definition. Repeats ignored")
    q <- qr
  }
  
  #  if(q>6){
  # ## documentation uses "k" not q
  #  return( message(paste("Cannot create an srd for k=", q, "rectangles. Up to 6 allowed")))
  #  }
  # ## check on all factors?
##---------------------------------------------------------------  
## coerce data which must be dataframe of factors into 0/1 matrix
   
  datavars <- colnames(data)
  datavars_in <- datavars
  data_in <- data

  datam <- matrix(nrow=nrow(data),ncol=q)
  datam <- as.data.frame(datam)
  
  ## loop over rectangles[] to get variable and value
  ## from string of form variable=value
  ## the last = sign is used as separator. 
  ## eg in  age<60 & sex=Female = 1  the varaiable is 
  ## age<60 & sex=Female and value 1.
  remove <- NULL
  cexpr <- NULL
  
  for(i in 1:q){
  
   XXX <- eval(parse(text=rectangles[i]),envir=data)

   if(class(XXX) != "logical" & class(XXX) != "integer" & class(XXX) != "numeric" ){
     return(message("rectangles item ",paste0("\"",rectangles[i],
                   "\" is not recognised")))
   }
   datam[,i] <- XXX
   

## is it binary?
U <- unique(XXX)
##print(U)
##  could be an NA, if variable(s) missing. Remove
NA_U <- which(is.na(U))
if(length(NA_U ==1)){
U <- U[-which(is.na(U))]}

##  c(0,1) treated same as c(FALSE,TRUE)

if(!setequal(U,c(0,1))){

   if(length(U)==1){
    if(U==1 | U==0){
    if(U ==0 ){
      note <- paste0("No TRUE values of \"", rectangles[i],"\". Is omitted")
      message(note)
      notes <- append(notes,note)
      remove <- c(remove,i) }
    if(U==1 ){
    note <- paste0("All observations have same \"", rectangles[i],"\" value")
      message(note)
      notes <- append(notes,note)
    }
     }}
     else  ##ilength(U)==1
     { note <- paste0( "rectangles item \"",rectangles[i],"\" does not evaluate to binary. Is omitted")
     message(note)
     notes <- append(notes,note)
     remove <- c(remove,i)
     }
   }  ##if(!setequal(U,c(0,1))){
  
  colnames(datam)[i] <- rectangles[i]
  
  }
  
## replace input data with 0/1 matrix data 
 
  if(!is.null(remove))
  {datam <- datam[,-remove]
  rectangles <- rectangles[-remove]
  q <- q -length(remove)

  }
# 
#   if(!is.null(repeats) ){
#     
#    note <- paste("Repetition in rectangles definition. Repeats ",paste(repeats),"ignored")
#    message(note)
#    notes <- append(notes,note)
#    
#   }
  
  ## are there identical data?

  remove <- NULL
  if(q>1){
  for( i in 1:(q-1)){
    for (j in (i+1):q){

      id_ij <- identical(as.logical(datam[,i]),as.logical(datam[,j]) )
      if(id_ij) {
        note <- paste(rectangles[i],"and",rectangles[j], "match perfectly")
        message(note)
        notes <- append(notes,note)
        note <- paste(rectangles[j],"is ignored")
        message(note)
        notes <- append(notes,note)
        
         remove <- c(remove,j) 
      }
    }
    }
  } 
## also remove perfect matches  

  remove <- unique(remove)
    if(!is.null(remove))
    {datam <- datam[,-remove]
    rectangles <- rectangles[-remove]
    q <- q -length(remove)
    
    }
   
  # ## remove any  zero size rectangles 
  # remove <- NULL
  # browser()
  # for(i in 1:q){
  #   
  #   siz <- sum(datam[,i,drop=FALSE]*weight)
  #    if(siz == 0) remove <- c(remove,i)
  #   print(paste("siz",siz))
  # }
  # if(!is.null(remove)){
  #   q <- q - length(remove)
  #   datam <- datam[,-remove]
  # }
  ##browser()
  data <- datam 
  ## possible that q = 0 due to undefined rectangles?
  if(q  <= 0){
    return(message("q = 0, no non-empty rectangles!"))
  }
  if(q>6){
    ## documentation uses "k" not q
    return( message(paste("Cannot create an srd for", q, "rectangles. <= 6 allowed")))
  }
  
  ## possible data not a dataframe if only one rectangle
  ## coerce
  data <- as.data.frame(data)
  datavars <- colnames(data)
  

##------------------------------------------------  

  
 
  zname <- z
  is.naz <- FALSE
 if(!is.null(z) ){
   if(length(z) > 2){
     message("length(z) argument should be <= 2. Truncated to 2")
     z <- z[1:2]
   }

   if(is.na(z[1]) ){
    is.naz <- TRUE 
 ## form z which will interleave in 3D by values of rectangles
 ##browser()
     
# Option 2. Billy's code to pick out first rectangle assign values 0,1,2,
# y<- apply(data, 1, function(x){
#   v <- which(x==1)[1]
#   replace(v, which(is.na(v)), 0)
#     })
## use, for a cell value of z,  whichever is first listed?
    z <- apply(data,1,function(x){ which(x==1)[1]})
    ## whitespace NA ,make 0
    z[which(is.na(z))]  <- 0 
## note: this is just now used to "fill in" z with non- NA 
## now use rectangle numbers as z value. This done via the 
## is.naz switch globally assigned    
    z <- as.data.frame(z,drop=FALSE)

colnames(z) <- "rectangle #"
##browser()
zname <- colnames(z)
note <- "z=NA makes z axis rectangle number "
message(note)
notes <- append(notes,note)

   }
 else
 {
   for(i in 1:length(z)){
   if(!is.element(z[i],datavars_in) ) {
     return(message(paste0(" z= \"",paste(z[i]), "\" not in data")))
   }
   }
   ## messy, but using z as name on input to now be the  data. Requires interchange
   ##  eg.  or z="A", or z=c("A","B")
      zname <- z

  ## aaagh, dangerous change of use of z. Input as string, now is a dataset    
      z <- as.data.frame(data_in[,zname],drop=FALSE)
      colnames(z) <- zname 
## what if  factors z?
## interesting! Even if data_in[,zname] is "character", it is converted to "factor"
## presumably by as.data.frame() so dealing with character variables for z!
## but if character, default levels correspond to alpha-numeric order of categories
      if(is.factor(z[,1])){

        lev <- levels(z[,1])
        lvX <- min(length(lev), 20)
        z[,1]  <- as.numeric(z[,1]) -1
        note <- paste0("z \"",colnames(z)[1],"\" coerced to numeric factor level - 1")
        message(note)
        notes <- append(notes,note)
        note <- paste(lev[1:lvX],"=", seq(1:lvX) - 1, " ")
        message(note)
        
        notes <- append(notes,note)
        if(length(lev) > lvX){
         note <- paste(length(lev)-lvX,"values omitted")
         message(note)
         notes <- append(notes,note)
        }
        
      }
      ## same for second z item. 
      if(ncol(z) >= 2){
      if(is.factor(z[,2])){
        lev <- levels(z[,2])
        lvX <- min(length(lev), 20)
        z[,2]  <- as.numeric(z[,2]) -1 
        note <- paste0("z \"",colnames(z)[2],"\" coerced to numeric factor level - 1")
        message(note)
        notes <- append(notes,note)
        
        note <- paste(lev[1:lvX],"=", seq(1:lvX) - 1, " ")
        message(note)
        notes <- append(notes,note)
        if(length(lev) > lvX){
          note <- paste(length(lev)-lvX,"values omitted")
          message(note)
          notes <- append(notes,note)
        }
        
        
      }}
      
 }
   ## note NA values of z?
   # ##browser()
   # NAzs <- length(which(is.na(z[,1])))
   # if(NAzs >0){
   #   message("There are ", NAzs," missing values of ",colnames(z)[1])
   # }
   # if(ncol(z) >= 2){
   #   NAzs <- length(which(is.na(z[,2])))
   #   if(NAzs >0){
   #     message("There are ", NAzs," missing values of ",colnames(z)[2])
   #   }
   # }
       
   
}  ##if(!is.null(z))
## do complete cases check
  dfz <- data
if(!is.null(z)) {dfz <- cbind(dfz,z)}
  dfz <- cbind(dfz,weight)
##-------------------------------------------------  
## check on complete cases before any collapse
  
  inall <- nrow(dfz)  
  cc <- complete.cases(dfz)
  ccs <- which(cc==TRUE)
  
  if(sum(cc) < inall){
    note <- paste("Missing data: Complete cases ",sum(cc)," <",inall)
    message(note)
    notes <- append(notes,note)
    
    
  }

  
 ##-------------------------------------------------------- 
 ## big data?  try collapsing, or always?
## yes , always. Better. Ensures automatic encode of factor z.

  if(nrow(data) >10000 | TRUE ){
  nbefore <- nrow(data)
  dfz <- SumDupRows(dfz, 1:(ncol(dfz)-1), "weight", "weight")

  if(!is.null(z)){
    if(ncol(z)==2){
## possible ratio for z  
      data <- dfz[,1:(ncol(dfz)-3)]
      z <- dfz[,(ncol(dfz)-2):(ncol(dfz)-1)] 

    }
    else
    {
    data <- dfz[,1:(ncol(dfz)-2)]
    z <- dfz[,(ncol(dfz)-1),drop=FALSE] 
    } 
  }
  else
  {data <- dfz[,1:(ncol(dfz)-1)]}
  
  weight <- dfz[,(ncol(dfz)),drop=FALSE]
###notes <- append(notes,paste("Data collapse from ",nbefore," to ",nrow(data)," records"))
  }  ## if(nrow(data)>10000)

##--------------------------------------------------------
##post collapse  cc testing

   
  inall <- nrow(dfz)  
  cc <- complete.cases(dfz)
  ccs <- which(cc==TRUE)

  if(sum(cc) < inall){
    data <- as.data.frame(data,drop=FALSE)
    coln <- colnames(data)
    ## messy dealing with q=1 retaining names
    ## what are missing items?
    for(i in 1:ncol(data)){

      ISNA <- which(is.na(data[,i]))
      if(length(ISNA)>0){
        ##account for "weight"
        mssdi <- sum(dfz[ISNA,"weight"])
      note <- paste(colnames(data[i]),"has",mssdi,"missing")
      message(note)
      notes <- append(notes,note)
      }
    }

    data <- data[ccs,]
    if(q==1){
      data <- as.data.frame(data,drop=FALSE)
      colnames(data) <- coln
    }

    if(!is.null(z)){
      z <- as.data.frame(z[ccs,],drop=FALSE)
      colnames(z) <- zname}
    
## pick out ccs rows of  weight too
## patch this weight, may or may not be data.frame

    if(class(weight)=="data.frame"){
      weight<- weight[ccs,]}
    else
    {
    weight<- weight[ccs]}
  }

  ## test reweight?
  ##browser()
  if(!is.null(reweight)){
    reweight <- gsub(reweight,pattern=" ", replacement="")
    wvar <- which(colnames(data)==reweight)
    

    if(length(wvar) == 1){
      ## calculate new weighting for pseudo-counts
    weight <- reweight(data[,wvar],weight=weight,alpha=0.5)
    note <- paste("Pseudo-counts based on reweighting of ",
                  reweight, "to 0.5n")
    message(note)
    notes <- append(notes,note)
    }
    else
    {message("reweight argument does not match one of rectangles")}
  }
  if(!is.null(w)){
    notes <- append(notes,paste("Frequency weighted. Total obs =", signif(sum(weight),4))) 
  }
##---------------------------------------------------------  
  ## re-order data so that largest rectangle is first, 
  ## with smaller on top. But the order may be changed if z != NULL
  ## to ensure correct layering in rotation. 
  ## coerce to ensure dataframe data
  

 data <- as.data.frame(data, drop=FALSE)
 rsize <- vector(length=q)
 emptyrectum <- NULL
 for(i in 1:q){
   
   rsize[i] <- sum(data[,i,drop=FALSE]*weight)
   ##    print(paste("rsize[i]=",rsize[i]))
   if(rsize[i] == 0)emptyrectum <- c(emptyrectum,i)
 }
 ##  remove possible empty rectangles

 if(length(emptyrectum) >0) {
   
   q <- q-length(emptyrectum) 
   
   message(paste("removing empty rectangle(s)",colnames(data)[emptyrectum]))
   
   if(q==0){
     return(message("q = 0, no non-empty rectangles!"))
   }
   rsize <- rsize[-emptyrectum]

   colnameshold <- colnames(data)[-emptyrectum]
   data <- as.data.frame(data[,-emptyrectum], drop=FALSE)
   colnames(data) <- colnameshold
 }

 
 ## et ordering of rectangles  
 
 rorder <-  seq(1:q) 

 if(layer){
   rorder <- order(rsize,decreasing=TRUE)
}  #if(layer)

if(usercols){
## if user colours specified, need to also re-order nofillc so that
## colours col1,col2,col3..line up with  rectangles 1, 2, 3..
  nofillc[1:q] <- nofillc[rorder]
}

if(q>1) data <- data[,rorder]
## need this to deal with q=1:  why??
# browser()
# if(q==1){
# data <- as.data.frame(data, drop=FALSE)
# ##colnames(data)[1] <- rectangles[1]
# 
# }
# 
# browser()

 
##------------------------------------------------------- 
 ## check on whther entirely nested?
 indices <- NULL
 for (i in 1:nrow(data)){
 ibin <- sum(2^(which(data[i,]==1)-1)) +1 
 indices <- c(indices,ibin)
 ##print(ibin)
 }
 allin <- 2^(seq(1:q))
  indices <- unique(indices)
  ## neednt worry about whitespace (value 1) 
  indices <- indices[which(indices >1)]
  nested <-  setequal(allin, indices) 
  
if(nested) {thin <- 0.5 
message("rectangles are entirely nested")
}
  ##-------------------------------------------------------
##  check on whether disjoint 
  disjoint  <- all( rowSums(data) < 2)

  if( disjoint ) {
    thin <- 0
    message("rectangles are entirely disjoint")
  }
  
##============================================================
##  end of additional pre-analysis
##=================================================================

 #  default initial  "shading" (colour). (6 =#BFFFD680)
 colour_num <- 5  ## default
 col <- defaultcol
 ws_col <- "white"
## if user colours colour_num=0
if(usercols) {colour_num <- 0
col <- as.vector(nofillc)
ws_col <- user_ws_col}


# lst <- list("none","earthy","bright","mono","opaque","red_y","mono_y","blue_y")

if(!is.null(z)){ 

blues <- c("#ffffccb3","#c7e9b4b3","#7fcdbbb3","#41b6c4b3",
           "#2c7fb8b3","#253494b3","#000076b3")

blues11 <- c("#ffffd9b3","#eff9bdb3","#d5eeb3b3","#a9ddb7b3","#73c9bdb3","#45b4c2b3"
  ,"#2897bfb3",  "#2073b2b3","#234ea0b3","#1c3185b3","#081d58b3")
## attempt (not v successful!) to ensure continuous spectrum of shades for  spanr() partitions
## later:  dont think this is working well, better to have A and !A dealth with same way?
if(!is.null(spanr) & FALSE){
 if(spanr=="A" ) blues <- blues11[6:11]
 if(spanr=="!A") blues <- blues11[1:6]
}
## default to intensity shading? 
##colour_num <- 6 
##col <- blues

}  ##if(!is.null(z))
 ##---------------------------------------------------------------
## check on criterion
 lst <- list("lsqs","labs","logl","chisq","stress")
  xxxx <- criterion
  crit <- which(lst==criterion)
  if(is.na(crit[1])) {
  return(message(paste("Invalid criterion=",paste("\"",xxxx,"\"",sep=""),". Allowed: \"logl\" \"lsqs\" \"labs\" \"chisq\" \"stress\" ")))

  }
 ##---------------------------------------------------------------- 
  ## this is "pointer" toggle in dialog
  labelrectangle <- 1   ##  2 pointers on right
##  labelrectangle <- 0  ## no labels
  # lst <- list("none","inside","arrow")
  ## this is cell stats toggle in dialog
  labelcell <- 1

## can only accommodate up to 6 rectangles 
  if(ncol(data)> 6){
   return(message(paste("Cannot create an SRD for", ncol(data), "rectangles. Maximum 6 allowed")))
 
  }

thinX <- thin
if(new==TRUE){
## thinness parameter if NA is equivalent to negative - used to signal "automatic" penalty
##  in srd.f  I have default setting of this parameter when thin=0 too  

  if(is.na(thin)==TRUE ) {thin <- -1}


}
else
  
  ##new =FALSE 
{
message("Re-creating existing configuration")
note <- paste("Re-creating existing srd for ", ncol(data)," rectangles")

notes <- append(notes,note)

}
##------------------------------------------------------------------------------
## deal with  z data frame. Note at this point
## z has changed role as the input variable names. See  the "aagh" comment above
## it is here a data frame.  input z  is here "zname" 


if(!is.null(z) ){
 flat <- FALSE 
  isz_dual <- ncol(z) >= 2

 
  if(isz_dual) { 
    ## make a name that is "mean(a)/mean(b)" to replace two valued zname
    zname <- paste(zname[1],"/",zname[2],sep="")
    w <- z[,2] }
 
 else
 { 
   ## coerce a vector of 1's for w
   w <- rep(1,times=nrow(data))
  }
  
  
 ## remove spaces from zname  
  zname <- gsub(pattern=" ",replacement="",x=zname)
  
#  extract first column for z, call it z 
#  second item is w (as set above)   
z <- z[,1] 
##}  ## of if(is.data.frame(z)==TRUE)
}

else  #of if(!is.null(z)
## z=NULL, no point in doing 3D if requested
{
  flat <- TRUE
  z <- rep(0, times=nrow(data))
  w <- rep(1, times=nrow(data))
 
}  #if(!is.null(z)
##------------------------------------------------------------------------------



Esc <- !clickable

##start_angle <- c(0,0,0)
i1 <- 1
# if(!flat){
#   i1 <- 7
#   
#   ## start oblique projection ?  No 
# start_angle <- c(60,60,60)
# }


## keep notes for up to  first draw.srd call
basenotes <- notes

## form initial diagram. data is binary attribute dataframe 
##browser()
ngoes <- 1
E <- 10000
## experiment with looping for initial config. 
## not good because get sequence diagrams flashing
##while(E > 5 & ngoes <4){
## xdummy, z dummy are empty matrices of coordinates in unit square
## set here there are transfered to srd.f as Nan.
## if srd.f is input with actual values initial positions
## are established in srd.f. If non-Nan input values are assumed
## to initail positions for a "nudge". Yet to be implemented!
xdummy <- matrix(nrow=4,ncol=7)
ydummy <- matrix(nrow=4,ncol=7)

## these are dummy ()
## temp routine to match grandm in srd.f
## testing process to input into srd.f
#   to experiment with "nudge" 
#  XXX <- grandm(data,weight,q,xdummy,xdummy)
#  xdummy <- XXX[[1]]
#  xdummy <- XXX[[2]]

##  z is here first value of z,  w is second value
## unless z is dual valued w=1,1,1,1,...
## z passed as Y into draw.srd. Just to confuse!
col <- col[1:q]
perm <- seq(1:q)
bgtextcol <- "black"
backgcolour <- bg
if(backgcolour=="#000000" | backgcolour=="black") bgtextcol <- "#ffffff"
indep <- FALSE 
## can we make bgtextcol global like this??
##assign("bgtextcol", bgtextcol, envir = .GlobalEnv)

par(bg = backgcolour, lwd=1)
##if(colour==3)par(bg = "#EFEFEF", lwd=1) #change background colour for "mono"

whites <- which(rowSums(data)==0)

if( length(whites)==0  & whitespace  ){
##  message("No whitespace. No point in whitespace=TRUE. Made FALSE")
##  whitespace <- FALSE
}
##whitespace=F and  there is whitespace. Remomve whitespace from data
## logically yes, but seem better (e.g. feinsteins ld) not to
# if( length( whites)>0  & !whitespace  ){
#   data <- data[-whites,,drop=FALSE]
#   weight <- weight[-whites]
#   z <- z[-whites]
#   w <- w[-whites]
# }



##-------------------------------------------------------------------------------
## Create Initial configuration
##--------------------------------------------------

drawn <- draw.srd(d=data,weight, labelcell=labelcell,new=new,seed=seed,colour=colour_num,col, crit, 
        labelrectangle=labelrectangle,Y=z,project=start_angle,cells=cells, w,thin,
         zname,title,whitespace=whitespace,ws_col,spanr=spanr,xdummy,ydummy,perm,
        borders,is.naz,bgtextcol,indep)

dataX <- data
 ##xdummy <- drawn$x
 ##ydummy <- drawn$y
# print(xdummy)
##-------------------------------------------------------------------------------

## end  x, y configuration

xend <- drawn$x
yend <- drawn$y


# message(paste("Fitted",criterion,
#               " criterion (optimised value =",signif(drawn$dactm,5),")") )
message(paste("Fitted",criterion," criterion " ))
message(paste0("Fit chi-square (Pvalue) ", 
              signif(drawn$chisqu,4)," (",signif(drawn$pgof,4),")"))

E <- drawn$E

if(drawn$fail)  {
  return(message("Drawing failed. Program stopping"))
}


seed_used <- drawn$seed_used
seedp      <- NA

outtable <- drawn$outtable
outrect <- drawn$outrect
notes <- append(basenotes,drawn$notes)
refit_done <- FALSE
pre_seed <- vector(length=100)
pre_crit <- vector(length=100)
pre_perm <- matrix(ncol=q,nrow=100)
pre_perm[1,] <- seq(1:q)
pre_crit[1] <- crit
pre_seed[1] <- seed_used

##=====================================================================

exp_zname <- NULL

if(!is.null(zname)) {
y1y2 <- unlist( strsplit(zname,"/"))
if(length(y1y2)==2){
  exp_zname <- paste0("Intensity: mean(", y1y2[1],")/mean(",y1y2[2],")")
}
else
{
  exp_zname <-paste0("Mean(", zname,")")}
}

lst <- c("no fill","random","bright","mono","opaque","default",exp_zname)
if(usercols) lst[1] <- "user colours"
message(paste("colours:",lst[colour_num+1]))


xwin <- c(-0.2 , 1.4)
ywin <- c(-0.2 , 1.2)

if(clickable) message('Click on graphic menu expected.  Click Esc to quit')
# i1 <- 1
# if(!is.null(y)) i1 <- 7
minE <- Inf
nrefit <- 1
ngoes <- 1
isflat <- 0
angle <- start_angle

## adding "dialog box"  to change plot type and to Exit
plot.window(xlim=xwin,ylim=ywin)
dialogpos <- ywin[2] 

gap <- (ywin[2] -ywin[1])/40


L <- xwin[1]
M <- xwin[2]

nundo <- NULL
while(Esc != TRUE){
 
#  
  ## adding "dialog box"  to change plot type and to Exit
  plot.window(xlim=xwin,ylim=ywin)
  dialogpos <- ywin[2] 
  
  gap <- (ywin[2] -ywin[1])/40
  
  
  L <- xwin[1]
  M <- xwin[2]
  
## menu box  
##browser()
rect( L, dialogpos-10.5*gap, L+0.1*(M-L), dialogpos+0.5*gap,col="white", border="blue")

text( L+0.05*(M-L),dialogpos,  paste("Esc"),cex=0.7,col="blue")
text( L+0.05*(M-L),dialogpos- 1*gap,  paste("Refit"),cex=0.7,col="blue")
text( L+0.05*(M-L),dialogpos-2*gap,  paste("Undo refit"),cex=0.7,col="blue")

blue <- "blue"
if(is.null(zname)) blue <- "gray"
 text( L+0.05*(M-L),dialogpos-3*gap,  paste("Colour"),cex=0.7,col="blue")
# if(!cells){
# text( L+0.05*(M-L),dialogpos-4*gap,  paste("Cells"),cex=0.7,col="blue")}
# else
# {text( L+0.05*(M-L),dialogpos-4*gap,  paste("Rectangles"),cex=0.7,col="blue")}
text( L+0.05*(M-L),dialogpos-4*gap,  paste("Cell stats"),cex=0.7,col="blue")
text( L+0.05*(M-L),dialogpos-5*gap,  paste("Labels"),cex=0.7,col="blue")
blue <- "blue"
if(is.null(zname)) blue <- "gray"

text( L+0.05*(M-L),dialogpos-6*gap,  paste("Rotate"),cex=0.7,col=blue)
text( L+0.05*(M-L),dialogpos-7*gap,  paste("Mode"),cex=0.7,col=blue)
text( L+0.05*(M-L),dialogpos-8*gap,  paste("Toggle view"),cex=0.7,col=blue)
## adding thickened borders
text( L+0.05*(M-L),dialogpos-9*gap,  paste("Borders"),cex=0.7,col="blue")
text( L+0.05*(M-L),dialogpos-10*gap,  paste("E(counts)"),cex=0.7,col="blue")


XY <- locator(n=1)
  
  
  plot_select <- floor( (dialogpos-XY$y)/gap + 1.5)
  insidebox <- (XY$x >=L &  XY$x <= L+0.1*(M-L))
  ##  click on dialog area


if(length(plot_select) == 0  )   {
##exit via press of Esc button (makes length(plot_select) == 0 )
##  coerce to ensure Esc=TRUE below
  plot_select  <- 1
  insidebox <- TRUE
}


##insidebox <- (XY$x >=L &  XY$x <= L+0.1*(M-L))

if(insidebox) {
##------------------------------------------------------------  
## Esc:
##------------------------------------------------------------
  if(plot_select == 1  )   {
        ## Exit by clicking EXIT menu item
        Esc <- TRUE 
      }
  ##------------------------------------------------------------
  ## refit
  ##------------------------------------------------------------
  if(plot_select == 2 & q == 1){ 
    message("no refitting with q=1")}
  
  if(plot_select == 2 & q > 1){ 
    
    ## record number of one step back from about to be refitted
    nundo <- nrefit 
    if(nundo > 100){
      message("Warning: >100 refits. ")
    }
    crit <- crit +1
    # ## increase thin too ??
    # ##     thin <- thin+0.1 
    if(crit==6) {
      crit <- 1
      
    }
    
    ## also permute the order of rectangles
    ## this screws up layering because data  is
    ## ordered to ensure largest rectangle ar back. 
    ## draw.srd expects this to be the case. Once permuted
    ##  no longer correctly layers.
    ## Unsure how to fix! For time being close off permuting
    ##  May have solved by passing  "perm" thro to draw.srd.enc()
    ## and using "order(perm) to get ordering for laying down rectangles. 
    ##perm <- sample(seq(1:q),replace=FALSE)
    ##print(perm)
    ## also screws up  reproducing. Stick with non-permuting??  
    ##perm <- seq(1:q)
    dataX <- data[,perm,drop=FALSE]
    colX <- col[perm]
##    print(perm)
    if(colour_num == 6 ) colX <- blues
    # renew with a new  seed
    seed <- NA
   ## print(xend)
   
    
    drawn <- draw.srd(d=dataX, weight, labelcell=labelcell,new=TRUE,seed=seed,colour=colour_num,colX,crit, 
                      labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
                      w,thin,zname,title,whitespace=whitespace,ws_col,spanr ,xdummy,ydummy,perm,
                      borders,is.naz,bgtextcol,indep)
    
   ## print(drawn$x)
    xend <- drawn$x
    yend <- drawn$y
    
     ##xdummy <- drawn$x
     ##ydummy <- drawn$y

    
    criteria <- c("lsqs", "labs", "logl", "chisq","stress")
    criterion <- criteria[crit]
    
 ##   message( paste("Refitted with",criterion,"criterion (optimised value =",signif(drawn$dactm,5),")") )
    message( paste("Refitted with",criterion,"criterion ") )
    message( paste0("Fit chi-sq (Pvalue) ", 
                    signif(drawn$chisqu,4)," (",signif(drawn$pgof,4),")"))
  ##  message(paste("refitted seed=",seed))
    
    outtable <- drawn$outtable
    outrect <- drawn$outrect
    
    notes <- append(basenotes,drawn$notes)
    seed_used <- drawn$seed_used
    new <- drawn$new
    minE <- min(minE,drawn$E) 
    ngoes <- ngoes +1 
    # 
    # if(ngoes == 5 & minE>10){
    #   note <-  paste("After",ngoes,"refits min(E%)=",signif(minE,3),
    #                  "Achieving good cell area/freq congruence unlikely")
    #   message(note)
    # }
    
    refit_done <- TRUE
    
    nrefit <- nrefit +1 
    
    ## record all the fitted configurations
    ## (first is initial configuation when nrefit=0)
    pre_crit[nrefit] <- crit
    pre_seed[nrefit] <- seed_used
    pre_perm[nrefit,] <- perm
    
  } 
  
  
  ##------------------------------------------------------------
  ## undo refit: 
  ##------------------------------------------------------------
      if(plot_select == 3 &  !is.null(nundo)){ 
   
    new <- TRUE
    
# nundo is number of config just before last refit
##     
 seed <- pre_seed[nundo]
 crit <- pre_crit[nundo]
 perm <- pre_perm[nundo,]
 dataX <- data[,perm]
 
 
 colX <- col[perm]
 if(colour_num == 6) colX <- blues
 

      drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=TRUE,seed=seed,colour=colour_num,colX,crit, 
                labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
                w,thin,zname,title,whitespace=whitespace,ws_col,spanr,xdummy,ydummy ,perm,
                borders,is.naz,bgtextcol,indep)
      xend <- drawn$x
      yend <- drawn$y
      criteria <- c("lsqs", "labs", "logl", "chisq","stress")
      criterion <- criteria[crit]

      message( paste("Undone Refit") )
 
      
 outtable <- drawn$outtable
 outrect <- drawn$outrect

 notes <- append(basenotes,drawn$notes)
 seed_used <- drawn$seed_used
 new <- drawn$new
 minE <- min(minE,drawn$E)  
 refit_done <- FALSE
 ## browser()
 nrefit <- nundo
## step back for next undo 
  nundo <- max(1,nundo -1)

 
## reset rotation start position to 1.
 ##i1 <- 1
      } ##if(plot_select == 3)
  ##-------------------------------------------------------------
  # 
  # 
  # ## nudge  trials: re-arrange with new seed and fitting criterion
  # ## unsuccessful. Wont actually "nudge", may spin off into a new configuration
  # if(FALSE & plot_select == 10){ 
  #   
  #   ## record number of one step back from about to be refitted
  #   nundo <- nrefit 
  # 
  #   
  #   # renew with a new  seed
  #   seed <- NA
  #   print(xend)
  #   ## refit  call   
  #   drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=TRUE,seed=seed,colour=colour_num,colX,crit, 
  #                     labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
  #                     w,thin,zname,title,whitespace=whitespace,ws_col,spanr ,xend,yend)
  #   
  #   print(drawn$x)
  #   xend <- drawn$x
  #   yend <- drawn$y
  #   
  #   criteria <- c("lsqs", "labs", "logl", "chisq")
  #   criterion <- criteria[crit]
  #   
  #   message( paste("Refitted",criterion," criterion with value =",signif(drawn$dactm,5)) )
  #   message( paste0("Fit chi-sq (Pvalue) ", 
  #                   signif(drawn$chisqu,4)," (",signif(drawn$pgof,4),")"))
  #   
  #   
  #   outtable <- drawn$outtable
  #   outrect <- drawn$outrect
  #   
  #   notes <- append(basenotes,drawn$notes)
  #   seed_used <- drawn$seed_used
  #   new <- drawn$new
  #   minE <- min(minE,drawn$E) 
  #   ngoes <- ngoes +1 
  #   
  #   if(ngoes == 5 & minE>10){
  #     note <-  paste("After",ngoes,"refits min(E%)=",signif(minE,3),
  #                    "Achieving good cell area/freq congruence unlikely")
  #     message(note)
  #   }
  #   
  #   refit_done <- TRUE
  #   
  #   nrefit <- nrefit +1 
  #   
  #   ## record all the fitted configurations
  #   ## (first is initial configuation when nrefit=0)
  #   pre_crit[nrefit] <- crit
  #   pre_seed[nrefit] <- seed_used
  #   
  # } 
  # 
  # 
##------------------------------------------------------------    
##  colours
##------------------------------------------------------------
if(plot_select == 4){

 
  
  lst <- c("no fill","random","bright","mono","opaque","default",
           exp_zname)
  if(usercols) lst[1] <- "user colours"
  xxx <- 5
  if(!flat) xxx<- 6 
  
   ## increment colour toggle  
  colour_num <- colour_num +1

  if(colour_num > xxx) colour_num <- 0
  message(paste("colour:",lst[colour_num+1]))
  
  
  ws_col <- "white" 
  
## reset colours if  user defined
  if(colour_num==0){
    if(usercols){col <- as.vector(nofillc)
    ws_col <- user_ws_col}
    else
    {col <- as.vector(nofillX)
    ws_col <- "white"}
  }


 ## experiment with random colours?? 
  ## note: this doesnt work and colours continually being updated.
  ## need to bring col <- outside of  this function draw.srd.enc()
if(colour_num==1){col <- randomColor(count = 6, luminosity="light")
##paste in tranparency codes  b3 is 80%
col <- paste0(col,"b3")
}
## colour_num=2 bright colours
  if(colour_num==2) {col<-as.vector(c("#0000FF80","#FF000080","#00FF0080","#FFFF0080","#0FF0F080","#F0B0E080"))  
  }
##colour_num=3 black and white "mono"
  if(colour_num==3) {col<-as.vector(c("#bdbdbdb3","#969696b3","#636363b3","#252525b3","#111111b3","#101010b3"))
}
#3 default ordering  o  is as listed 

## colour_num =4  user cols ##opaque
  if(colour_num==4)  {
    col <- defaultcol 
    if(usercols){
      col <- as.vector(nofillc)}
    col <- substr(col, start=1,stop=7)


}

   ## colour_num=1 earthy colour_nums 
  if(colour_num==5) {
    col<-as.vector(c("#58D6B880","#A7C8CC80","#A9D85E80","#6B785B80","#D8AD7580","#BFFFD680"))
  col <- defaultcol}
  
  colX <- col[perm]
  if(colour_num==6) { colX <- blues}

  
  
drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=FALSE,seed=seed_used,colour=colour_num,colX, crit, 
         labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
          w,thin,zname,title,whitespace=whitespace,ws_col,spanr,xdummy,ydummy ,perm,
          borders,is.naz,bgtextcol,indep)
outtable <- drawn$outtable
outrect <- drawn$outrect


}
##-----------------------------------------------------------
## mode toggle mode
##----------------------------------------------------------- 
if(plot_select == 8 & !flat){
## cells
  
  if(!cells){
    message("-> cell mode")
  }
  else
  {message("-> rectangle mode")}
  
  
  new <- FALSE
  if( identical(angle,c(0,0,0) ) ) {new <- TRUE}

  cells <- !cells
  seedp <- seed_used
  colX <- col[perm]
  if(colour_num == 6 ) colX <- blues
  
  drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=FALSE,seed=seed_used,colour=colour_num, colX,crit, 
               labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
    w,thin,zname,title,whitespace=whitespace,ws_col,spanr,xdummy,ydummy ,perm,
    borders,is.naz,bgtextcol ,indep)
  outtable <- drawn$outtable
  outrect <- drawn$outrect
  notes <- append(basenotes,drawn$notes)
  seed_used <- drawn$seed_used
}
 
  ##---------------------------------------------------------- 
  ## cell stats. Toggle labelcell
  ##---------------------------------------------------------- 
 if(plot_select == 5){
   
   labelcell <- labelcell+1 
   xxx <- 7
   if(flat ) xxx <- 5
   if(labelcell > xxx ) labelcell <- 0 


## if( !( identical(angle,c(0,0,0))  ) ) {message("Cell stats only shown for non-oblique projection")}
##else
   colX <- col[perm]
   if(colour_num == 6 ) colX <- blues
   
drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=FALSE,seed=seed_used,colour=colour_num, colX,crit, 
         labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
         w,thin,zname,title,whitespace=whitespace,ws_col,spanr ,xdummy,ydummy,perm, 
         borders,is.naz,bgtextcol,indep)
outtable <- drawn$outtable
outrect <- drawn$outrect
notes   <- append(basenotes,drawn$notes)
##}
##seed_used <- drawn$seed_used
 } ##if(plot_select == 4)
##------------------------------------------------------------ 
## labels
##-------------------------------------------------------------  
   if(plot_select == 6){
  
 ## labels toggle labels
   labelrectangle <- labelrectangle + 1 
   
   if(labelrectangle > 3 ) labelrectangle <- 0 

  if( identical(angle,c(0,0,0) ) ) {new <- TRUE}
   colX <- col[perm]
   if(colour_num == 6 ) colX <- blues
   
drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=FALSE,seed=seed_used,colour=colour_num, colX, crit, 
         labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
         w,thin,zname,title,whitespace=whitespace,ws_col,spanr,xdummy,ydummy ,perm ,
         borders,is.naz,bgtextcol,indep)
outtable <- drawn$outtable
outrect <- drawn$outrect


   } ##if(plot_select == 5)
  
##---------------------------------------------------
## toggle view
##---------------------------------------------------
if(plot_select == 9 ){
 isflat <- isflat +1 
  if(isflat==3)     isflat <- 0
  ## sideways on. make not quite 90 to show colour?
  if(isflat==2)     angle <- c(90,90,90) 
  if(isflat==0)     angle <- c(0,0,0) 
  if(isflat==1)     angle <- c(60,60,60) 
## why need to new=TRUE, re-entering  .FORTRAN(srd)?
##  seems not to need. KANEON
  colX <- col[perm]
  if(colour_num == 6 ) colX <- blues
  
drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=FALSE,seed=seed_used,colour=colour_num, colX,crit, 
         labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
         w,thin,zname,title,whitespace=whitespace,ws_col,spanr ,xdummy,ydummy,perm,
         borders,is.naz,bgtextcol,indep)
outtable <- drawn$outtable
outrect <- drawn$outrect

i1 <- max(1,i1-1)
i1 <- 1
start_angle <- angle
}  ##if(plot_select == 9)

##------------------------------------------------------
## rotate:  create series of frames in 3D
##------------------------------------------------------
  nframes <- 0
 if (plot_select == 7 & !flat ) {
   if(drawn$fail){return(message("srd failed"))
      nframes <- nframes+1
     #browser()
     }
   if(flat){
     message("z is not specified.  Rotation rather pointless, but done anyway!")
   }

 ## this is switch to allow  stepping thro rotation with repeated
   ## clicks on "Rotate". Not sure how to make optional? 
   stop_start <- FALSE
   ## if cells, then turn of labelling of rectangles once into 3D projection
   ##if(cells){labelrectangle <- 0}
    #  manually step rotations, pausing with browser()
 
   seedp <- seed_used
   
   message("Click anywhere to rotate, or Quit")
   
   if(stop_start){
     if(i1==1)message("Repeat click rotate to cycle to side-on view")
     i1 <- i1+1
     ## step in 10 degrees?? 
     degs <- 10 
     i2 <- i1}
   else
   {i1 <- 2
   degs <- 5
   i2 <- 19 }
   ##fewer frames than movie, (1:90) by 5 degrees,  bit jerky but is ok
   i <- 2
   quit <- FALSE
   colX <- col[perm]
   if(colour_num == 6 ) colX <- blues
   
    while (!quit){
      xa <-  start_angle[1]+ degs*(i-1)
      ya <-  start_angle[2]+ degs*(i-1)
      za <-  start_angle[3]+ degs*(i-1)
      lr <- labelrectangle
      lc <- labelcell
      angle <- c(xa,ya,za)
##message(paste("angle",angle))
      drawn <- draw.srd(d=dataX,weight, labelcell=lc,seed=seed_used, new=FALSE,colour=colour_num, colX, crit, 
               labelrectangle=lr,Y=z,project=angle ,cells=cells,
               w,thin,zname,title,whitespace=whitespace,ws_col,spanr,xdummy,ydummy ,perm,
               borders,is.naz,bgtextcol,indep )
      plot.window(xlim=xwin,ylim=ywin)
 ## put in a  Quit  menu box      
      rect( L, dialogpos-0.5*gap, L+0.1*(M-L), dialogpos+0.5*gap,col="white", border="red")
      
      text( L+0.05*(M-L),dialogpos,  paste("Quit"),cex=0.7,col="red")
     
 
      XY <-  locator(n=1)
      
      
      step_select <- floor( (dialogpos-XY$y)/gap + 1.5)
      
      
      if(length(step_select) == 0  )   {
        ##quit via press of Esc button
        quit <- TRUE 
      }
      else
      {
 ## click on Quit item?     
      if(step_select==1 & XY$x >  L-0.02*(M-L) & XY$x <  L+0.12*(M-L)) {quit <- TRUE}
      else {i <- i+1 }
    }
      
if(i==1){

   ## Sys.sleep(0.5)
  outtable <- drawn$outtable
  outrect <- drawn$outrect
notes <- append(basenotes,drawn$notes)}  ##if(i==1)
      
seed_used <- drawn$seed_used
nframes <- nframes+1
##angle <- angle + c(5,5,5)
    }  ## for(i in i1:i2)
   start_angle <- c(xa,ya,za)
   
 } ##if(plot_select==7 & !flat)
##-----------------------------------------------------------
   
  if(plot_select==10) {borders <- !borders
  colX <- col[perm]
  drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=FALSE,seed=seed_used,colour=colour_num, colX, crit, 
                    labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
                    w,thin,zname,title,whitespace=whitespace,ws_col,spanr,xdummy,ydummy ,perm , 
                    borders,is.naz,bgtextcol,indep)
  }
  if(plot_select==11){
    indep <- !indep
    colX <- col[perm]
    title <- titleX
    if(indep) title <- "E(counts) under independence"
    drawn <- draw.srd(d=dataX,weight, labelcell=labelcell,new=TRUE,seed=seed_used,colour=colour_num, colX, crit, 
                      labelrectangle=labelrectangle,Y=z,project=angle,cells=cells, 
                      w,thin,zname,title=title,whitespace=whitespace,ws_col,spanr,xdummy,ydummy ,perm , 
                      borders,is.naz,bgtextcol,indep)
    outtable <- drawn$outtable
    outrect <- drawn$outrect
    
    
  }
  ##print(plot_select)
    
} # if(insidebox)


}  ##while(ESC !=TRUE)

## blank out erase the dialog menu box
rect( L-0.01*(M-L), dialogpos-10.7*gap, L+0.11*(M-L), dialogpos+0.55*gap,col=backgcolour,border=NA)  
# if(minE>10 & ngoes >3 ){
#   
#   note <-  paste("After",ngoes,"refits min(E%)=",signif(minE,3),
#                  "Achieving good cell area/freq congruence unlikely")
#   message(note)
#   notes <- append(notes, note ) }
notes <- append(notes, paste0("Fit chi-sq (Pvalue) ", 
                              signif(drawn$chisqu,4)," (",signif(drawn$pgof,4),")") )

rep <- paste0("seed =",seed_used,", criterion= '",criterion,"', thin=",thin)
notes <- append(notes, rep ) 
note <- paste(paste0("'",col,"'"),collapse=", ")
note <- paste0("col = c(",paste(note),")" )
notes <- append(notes,note)
##note <- paste0("col",seq(1:q),"=","\"",col[1:q],"\",",collapse=" ")
#make note of non-white whitespace

if(c0 != "white"){ 
  note <- paste0("col0 = '",c0,"'")
  notes <- append(notes,note)
}
## make note on graphic of seed and criterion to be able to reproduce
text( (xwin[1]+xwin[2])/2,ywin[1]-0.03*(ywin[2]-ywin[1]),
      paste0(rep), font=3 ,cex=0.8,col=bgtextcol)

## leaves last comma? 
note <- substr(note, start=1,stop=nchar(note)-1)

image <- recordPlot()
output <- list(title=title,cells=outtable,rectangles=outrect, notes=notes,image =image ,fail=drawn$fail)

class(output) <- "srd"

return(output) 

}  ##  end srd()
########################################################################
