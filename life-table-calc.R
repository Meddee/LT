

##############################################################################################################################
#STEP 1: Read in and review the population and death data
##############################################################################################################################

females<-read.table(file="http://www.demog.berkeley.edu/~eddieh/AppliedDemographyToolbox/StanfordCourseLifeTable/StanfordCourseMortalityData.csv",header=TRUE,sep=",")
females

##############################################################################################################################
#STEP 2: Read in or create the fundamental pieces of the life table (age groupings, deaths by age, population by age ->death rates by age
##############################################################################################################################

x <- c(0,1,5,10,15,20,25,35,45,55,65,75,85)
#note that R collapses a single column to a vector when it pulls out the result out of a data.frame
nDx <- females$Death.Count   #other syntax which produces the same result: females[[3]], females[,3], 
nKx <- females$Population
nMx <- nDx / nKx

##############################################################################################################################
#STEP 3: Read in the period life table function
##############################################################################################################################

life.table <- function( x, nMx){
  ## simple lifetable using Keyfitz and Flieger separation factors and 
  ## exponential tail of death distribution (to close out life table)
  b0 <- 0.07;   b1<- 1.7;      
  nmax <- length(x)
  #nMx = nDx/nKx   
  n <- c(diff(x),999)            	#width of the intervals
  nax <- n / 2;		            	# default to .5 of interval
  nax[1] <- b0 + b1 *nMx[1]    		# from Keyfitz & Flieger(1968)
  nax[2] <- 1.5  ;              
  nax[nmax] <- 1/nMx[nmax] 	  	# e_x at open age interval
  nqx <- (n*nMx) / (1 + (n-nax)*nMx)
  nqx<-ifelse( nqx > 1, 1, nqx);		# necessary for high nMx
  nqx[nmax] <- 1.0
  lx <- c(1,cumprod(1-nqx)) ;  		# survivorship lx
  lx <- lx[1:length(nMx)]
  ndx <- lx * nqx ;
  nLx <- n*lx - nax*ndx;       		# equivalent to n*l(x+n) + (n-nax)*ndx
  nLx[nmax] <- lx[nmax]*nax[nmax]
  Tx <- rev(cumsum(rev(nLx)))
  ex <- ifelse( lx[1:nmax] > 0, Tx/lx[1:nmax] , NA);
  lt <- data.frame(x=x,nax=nax,nmx=nMx,nqx=nqx,lx=lx,ndx=ndx,nLx=nLx,Tx=Tx,ex=ex)
  return(lt)
}

##############################################################################################################################
#STEP 4: Apply the function to the data, and review the created life table
##############################################################################################################################

females.life.table<-life.table(x,nMx)
females.life.table

#write.table(###, file="G:/###/###.csv", sep=",")







#####################################################

## Question 2

ourdata<-read.table("JapanMort.txt",header=TRUE)
tail(ourdata)
class(ourdata)
as.data.frame(ourdata)

num.Age<-as.numeric(levels(ourdata$Age))[ourdata$Age]# chnge factor to numeric
#class(yr2000.m)
class(ourdata$Age)

# year 2000 population 
Y00.m<-subset(ourdata, select=c(Year, Age, MalePop, MaleDeath), 
                 subset=(Year==2000 & num.Age<=100 ))  # yr 2000 male pop
Y00.f<-subset(ourdata, select=c(Year, Age, FemPop, FemDeath), 
                 subset=(Year==2000 & num.Age<=100 )) # yr 2000 female pop

Y01.m<-subset(ourdata, select=c(Year,Age, MalePop), 
              subset=(Year==2001 & num.Age<=100 ))    # yr 2001 male pop
Y01.f<-subset(ourdata, select=c(Year, FemPop), 
              subset=(Year==2001 & num.Age<=100 ))#     yr 2001 female pop

#check data
head(Y00.m);head(Y01.m);head(Y00.f);head(Y01.f)
tail(Y00.m)

#Get mid-yr pops
Nx.m<-apply(cbind(Y00.m$MalePop,Y01.m$MalePop),1,mean)
Nx.f<-apply(cbind(Y00.f$FemPop,Y01.f$FemPop),1,mean)

japanMales2000<-data.frame(x=Y00.m$Age, midpop=Nx.m, dea=Y00.m$MaleDeath)
japanFemales2000<-data.frame(x=Y00.f$Age, midpop=Nx.f, dea=Y00.f$FemDeath)

japanMales2000; head(japanMales2000); head(japanFemales2000)


# LIFETABLE FUNCTION

lifetable <- function(x, Nx, Dx) {
  # Estimating a period lifetable
  # Input parameters:
  # x = Age x
  # Nx = Mid Year Population at age x
  # Dx = Number of People dying at age x
  
  max<-length(x) # maximum number of years
 
  #mx
  mx<-Dx/Nx
  
  #qx
  qx<-mx/(1+0.5*mx)
  qx[max] <- 1.0 # close table
   
  #px
  px<-1-qx
  
  #lx
  lx <- 100000 * c(1,cumprod(px[1:max-1]))
  
  #OR
  #lx <- 100000 * c(1,cumprod(px))
  #lx <- lx[1:max]  
  
  
  # lx Using a For loop
  #lx<-c(100000,rep(0,max-1))
  #for( i in 1:max-1){  
  #lx[i+1]<-lx[i]*px[i]
  #}
  
  #dx  
  #dx= lx-lx+1= lx*qx
  dx<-lx*qx
  
  #Lx
  #Lx = lx+1 + 0.5*dx = lx - 0.5*dx 
  Lx= lx - 0.5*dx
  Lx[max]<-0.5*dx[max]
  
  #Lx using For Loop
  #Lx<-c(rep(0,max-1),0.5*dx[max]) 
  #for( i in 1:max-1){  
  # Lx[i]<-lx[i+1]+.5*dx[i]
  # }
  
  # Tx
  Tx <- rev(cumsum(rev(Lx)))
  
  #ex
  ex<-Tx/lx
  
  # prepare the outcomes  
  outcome <- data.frame(x=x, Nx=Nx, Dx=Dx, mx=mx, qx=qx, px=px, lx=lx, dx=dx ,
               Lx=Lx, Tx=Tx, ex=ex)
  
  # give the outcomes
  return(outcome)
}

# Call function-  Males
LTjapanMales2000 <- lifetable(x = japanMales2000$x,
                              Nx = japanMales2000$midpop,
                              Dx = japanMales2000$dea)

head(LTjapanMales2000 )

# Call function -  Females

LTjapanFemales2000 <- lifetable(x = japanFemales2000$x,
                              Nx = japanFemales2000$midpop,
                              Dx = japanFemales2000$dea)


head(LTjapanFemales2000)

