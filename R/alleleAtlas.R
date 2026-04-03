
#fasta is assumed to be imported by read.dna function in ape package with the arguments format="fasta" and as.character=T.
#corres: a list of a numeric vector of sequence of each location.
makeInput<-function(lon,lat,fasta,corres)
{
  if(length(lon)!=length(lat))
  {
    return("Error: lon and lat have different length.")
  }
  #characterize sequence
  seq<-character(length(fasta[,1]))
  for(i in 1:length(fasta[,1]))
  {
    seqNumber<-as.numeric(names(fasta[i,1]))
    seq[seqNumber]<-paste(fasta[i,],collapse="")
  }
  result<-character(length(lon))
  for(i in 1:length(lon))
  {
    temp<-seq[corres[[i]]]
    result[i]<-paste0(lon[i],",",lat[i],",",paste(temp,collapse=","))
  }
  return(result)
}

#export neighbor network file
makeNeighbor<-function(path,lat,lon,northPole,southPole,loop)
{
  if(is.null(lat)||is.null(lon))
  {
    return("Please specify lat and lon.")
  }
  if(lat<=0||lon<=0)
  {
    return("lat and lon must be positive numbers.")
  }
  out<-file(path,"w")
  n<-lat*lon
  temp<-as.character(n)
  writeLines(temp,out,sep="\n")
  for(i in 1:n)
  {
    neighbourCounter<-0
    up<-NULL
    #up neighbour
    if(!(i<=lon))
    {
      if(lon==1)
      {
        up<-i-lon
        neighbourCounter<-neighbourCounter+1
      }
      else if(lon==2)
      {
        up<-numeric(2)
        neighbourCounter<-neighbourCounter+2
        up[1]<-i-lon
        #upright
        if(i%%lon==0)
        {
          up[2]<-i-(lon-1)-lon
        }
        else
        {
          up[2]<-i+1-lon
        }
      }
      else
      {
        up<-numeric(3)
        up[2]<-i-lon
        neighbourCounter<-neighbourCounter+1
        #upright
        if(i%%lon==0&&loop)
        {
          up[3]<-i-(lon-1)-lon
          neighbourCounter<-neighbourCounter+1
        }
        else if(i%%lon!=0)
        {
          up[3]<-i+1-lon
          neighbourCounter<-neighbourCounter+1
        }
        #upleft
        if(i%%lon==1&&loop)
        {
          up[1]<-i+(lon-1)-lon
          neighbourCounter<-neighbourCounter+1
        }
        else if(i%%lon!=1)
        {
          up[1]<-i-1-lon
          neighbourCounter<-neighbourCounter+1
        }
      }
    }
    #neighboring to the North pole
    else if(lon>1&&northPole)
    {
      up<-numeric(lon-1)
      counter<-1
      neighbourCounter<-neighbourCounter+lon-1
      for(j in 1:lon)
      {
        if(i!=j)
        {
          up[counter]<-j
          counter<-counter+1
        }
      }
    }
    #down neighbour
    down<-NULL
    if(lat>1)
    {
      if(!((lat-1)*lon+1<=i&&i<=lat*lon))
      {
        if(lon==1)
        {
          down<-i+lon
          neighbourCounter<-neighbourCounter+1
        }
        else if(lon==2)
        {
          down<-numeric(2)
          down[1]<-i+lon
          neighbourCounter<-neighbourCounter+2
          #downright
          if(i%%lon==0)
          {
            down[2]<-i-(lon-1)+lon
          }
          else
          {
            down[2]<-i+1+lon
          }
        }
        else
        {
          down<-numeric(3)
          down[2]<-i+lon
          neighbourCounter<-neighbourCounter+1
          #downright
          if(i%%lon==0&&loop)
          {
            down[3]<-i-(lon-1)+lon
            neighbourCounter<-neighbourCounter+1
          }
          else if(i%%lon!=0)
          {
            down[3]<-i+1+lon
            neighbourCounter<-neighbourCounter+1
          }
          #upleft
          if(i%%lon==1&loop)
          {
            down[1]<-i+(lon-1)+lon
            neighbourCounter<-neighbourCounter+1
          }
          else if(i%%lon!=1)
          {
            down[1]<-i-1+lon
            neighbourCounter<-neighbourCounter+1
          }
        }
      }
      #neighboring to the South pole
      else if(lon>1&&southPole)
      {
        down<-numeric(lon-1)
        counter<-1
        neighbourCounter<-neighbourCounter+lon-1
        for(j in 1:lon)
        {
          if(i!=j+(lat-1)*lon)
          {
            down[counter]<-(lat-1)*lon+j
            counter<-counter+1
          }
        }
      }
    }
    #right and left neighbour
    right<-NULL
    left<-NULL
    if(!((i<=lon)&&northPole)&&!(((lat-1)*lon+1<=i&&i<=lat*lon)&&southPole))
    {
      if(lon>2)
      {
        #right neighbor
        if(i%%lon==0&&loop)
        {
          right<-i-(lon-1)
          neighbourCounter<-neighbourCounter+1
        }
        else if(i%%lon!=0)
        {
          right<-i+1
          neighbourCounter<-neighbourCounter+1
        }
        #left neighbor
        if(i%%lon==1&&loop)
        {
          left<-i+(lon-1)
          neighbourCounter<-neighbourCounter+1
        }
        else if(i%%lon!=1)
        {
          left<-i-1
          neighbourCounter<-neighbourCounter+1
        }
      }
      else if(lon==2)
      {
        if(i%%lon==0)
        {
          if(loop)
          {
            right<-i-(lon-1)
            neighbourCounter<-neighbourCounter+1
          }
          left<-i-1
          neighbourCounter<-neighbourCounter+1
        }
        else if(i%%lon==1)
        {
          if(loop)
          {
            left<-i+(lon-1)
            neighbourCounter<-neighbourCounter+1
          }
          right<-i+1
          neighbourCounter<-neighbourCounter+1
        }
      }
    }
    #delete 0s which are left as the initial values
    up<-up[!up %in% 0]
    down<-down[!down %in% 0]
    right<-right[!right %in% 0]
    left<-left[!left %in% 0]
    line<-as.character(c(i,neighbourCounter,up,down,left,right))
    line<-paste(line,collapse=" ")
    lastword<-substring(line,nchar(line),nchar(line))
    mat<-match(" ",lastword)
    if(!is.na(mat))
    {
      if(mat==1)
      {
        line<-substr(line,1,nchar(line)-1)
      }
    }
    writeLines(line,out,sep="\n")
  }
  close(out)
}

#make latlon window from location data
makeMinWindow<-function(data.loc,northPole,southPole,loop,forceMaxEast=NULL,margin=0.001)
{
  if(loop)
  {
    westLimit<- -180
    eastLimit<- 180
  }
  else
  {
    #check which makes a smaller window that crosses the meridian on 0 degree or 180 degree
    x<-max(data.loc[,1])-min(data.loc[,1])
    if(x<180)
    {
      westLimit<-min(data.loc[,1])-margin
      eastLimit<-max(data.loc[,1])+margin
    }
    else if(x>180)
    {
      westLimit<--1*max(data.loc[,1])-margin
      eastLimit<-min(data.loc[,1])+margin
    }
    else
    {
      if(is.null(forceMaxEast))
      {
        print("The unique minimum east-west window cannot be decided.
			Consider applying loop=TRUE to this function, or manually decide the window using forceMaxEast variable in this function.")
        return(NULL)
      }
      else if(forceMaxEast)
      {
        westLimit<-min(data.loc[,1])-margin
        eastLimit<-max(data.loc[,1])+margin
      }
      else if(!forceMaxEast)
      {
        westLimit<-max(data.loc[,1])-margin
        eastLimit<-min(data.loc[,1])+margin
      }
    }
  }
  if(northPole)
  {
    northLimit<-90
  }
  else
  {
    northLimit<-max(data.loc[,2])+margin
  }
  if(southPole)
  {
    southLimit<- -90
  }
  else
  {
    southLimit<-min(data.loc[,2])-margin
  }
  win<-c(westLimit,eastLimit,southLimit,northLimit)
  return(win)
}

makeSpatialPolygons<-function(autoWindow,lat,lon,minlon=-180,maxlon=180,minlat=-90,maxlat=90,data.loc=NULL,northPole=NULL,southPole=NULL,loop=NULL,forceMaxEast=NULL,margin=0.001)
{
  #generate automatically the minimum mesh window
  if(autoWindow)
  {
    win<-as.owin(makeMinWindow(data.loc,northPole,southPole,loop,forceMaxEast,margin))
  }
  #generate manually dedicated mesh window
  else
  {
    win<-as.owin(c(minlon,maxlon,minlat,maxlat))
  }
  ppp<-as.ppp(data.loc,win)
  data.poly<-as(as.im(ppp$window,dimyx=c(lat,lon)),"SpatialGridDataFrame")
  data.poly<-as(data.poly,"SpatialPolygons")
  #add number of counts per mesh
  d<-data.frame(Ncounts = rep(0, length(data.poly)))
  row.names(d)<-paste0("g",1:length(data.poly))
  data.counts<-SpatialPolygonsDataFrame(data.poly,d)
  #create an ID column for meshes
  data.counts@data$id<-seq(1,lat*lon)
  return(data.counts)
}

#calculate expected heterozygosity from the result
calcExpectedHeterozygosity<-function(freq)
{
  calcResult<-numeric(length(freq[,1]))
  #check each mesh
  for(i in 1:length(freq[,1]))
  {
    calcResult[i]<-1
    #check each allele
    n<-length(freq[1,])
    for(j in 1:n)
    {
      calcResult[i]<-calcResult[i]-freq[i,j]^2
    }
    #correction for unbiasedness
    #but interpolated mesh does not have sample size
    #calcResult[i]<-calcResult[i]*2*n/(2*n-1)
  }
  return(calcResult)
}
#interprete a nucleotide base, especially mixed bases.
#returns probability of being A, T, G and C.
interpreteBase<-function(char)
{
  rateA<-0
  rateT<-0
  rateG<-0
  rateC<-0
  baseChar<-toupper(char)
  if(baseChar=='A')
  {
    rateA<-1
  }
  else if (baseChar=='T')
  {
    rateT<-1
  }
  else if (baseChar=='G')
  {
    rateG<-1
  }
  else if (baseChar=='C')
  {
    rateC<-1
  }
  #N is frequent, so put it former than the other mixed bases
  else if (baseChar == 'N')
  {
    rateC<-1 / 4
    rateA<-1 / 4
    rateT<-1 / 4
    rateG<-1 / 4
  }
  else if (baseChar=='R')
  {
    rateA<-1 / 2
    rateG<-1 / 2
  }
  else if (baseChar=='M')
  {
    rateA<-1 / 2
    rateC<-1 / 2
  }
  else if (baseChar=='W')
  {
    rateA<-1 / 2
    rateT<-1 / 2
  }
  else if (baseChar=='S')
  {
    rateC<-1 / 2
    rateG<-1 / 2
  }
  else if (baseChar=='Y')
  {
    rateC<-1 / 2
    rateT<-1 / 2
  }
  else if (baseChar=='K')
  {
    rateT<-1 / 2
    rateG<-1 / 2
  }
  else if (baseChar=='H')
  {
    rateC<-1 / 3
    rateA<-1 / 3
    rateT<-1 / 3
  }
  else if (baseChar=='B')
  {
    rateC<-1 / 3
    rateG<-1 / 3
    rateT<-1 / 3
  }
  else if (baseChar=='D')
  {
    rateG<-1 / 3
    rateA<-1 / 3
    rateT<-1 / 3
  }
  else if (baseChar=='V')
  {
    rateC<-1 / 3
    rateA<-1 / 3
    rateG<-1 / 3
  }
  else
  {
    print(paste0("Error: Detected an irregular base: ",char))
    return(c(NA,NA,NA,NA))
  }
  ATGCrate<-c(rateA,rateT,rateG,rateC)
  return(ATGCrate)
}

#count substitution number between two alleles
seq.dist<-function(allele1,allele2)
{
  vec1<-strsplit(allele1,"")[[1]]
  vec2<-strsplit(allele2,"")[[1]]
  if(length(vec1)!=length(vec2))
  {
    return("allele1 and allele2 must be the same length.")
  }
  comparableBase<-0
  dist<-0
  for(i in 1:length(vec1))
  {
    #skip gap bases
    if(vec1[i]=="-"||vec2[i]=="-")
    {
      next
    }
    v1<-interpreteBase(vec1[i])
    v2<-interpreteBase(vec2[i])
    if(any(is.na(v1)) || any(is.na(v2)))
    {
      #skip wrong base
      next
    }
    dist<- dist + sum(abs(v1-v2)/2)
    comparableBase<-comparableBase + 1
  }
  return(list(difference=dist,sites=comparableBase))
}

#calculate nucleotide diversity from the result
calcNucleotideDiversity<-function(freq,alleleSequence)
{
  alleleCount<-length(alleleSequence)
  alleleLength<-nchar(alleleSequence[1])
  #construct base substitution matrix
  subMatrix<-matrix(0,nrow=alleleCount,ncol=alleleCount)
  for(i in 1:alleleCount)
  {
    for(j in 1:alleleCount)
    {
      seqCompare<-seq.dist(alleleSequence[i],alleleSequence[j])
      diff<-seqCompare[[1]]
      comparableSite<-seqCompare[[2]]
      if(comparableSite==0)
      {
        subMatrix[i,j]<-NA
      }
      else
      {
        subMatrix[i,j]<-diff/comparableSite
      }
    }
  }
  alternative<-ncol(freq)
  locNum<-nrow(freq)
  calcResult<-numeric(locNum)
  if(alternative<2)
  {
    #when there is only one allele in all locations
    return(calcResult)
  }
  #check each mesh
  for(i in 1:locNum)
  {
    calcResult[i]<-0
    #check each allele
    for(j in 1:(alternative-1))
    {
      for(k in (j+1):alternative)
      {
        calcResult[i]<-calcResult[i] + freq[i,j] * freq[i,k] * subMatrix[j,k]
      }
    }
    #correction for unbiasedness
    #but interpolated mesh does not have sample size
    #calcResult[i]<-calcResult[i]*n/(n-1)
  }
  return(2*calcResult)
}

#functions to process input data
processInput<-function(data)
{
  data.loc<-matrix(0,nrow=length(data),ncol=2)
  data.maxCol<-0
  for(i in 1:length(data))
  {
    split<-strsplit(data[i],",")
    lon<-parse_lon(split[[1]][1])
    lat<-parse_lat(split[[1]][2])
    if(length(split[[1]])<3)
    {
      print(paste0("Location in row ", i , " has no allele."))
    }
    if(data.maxCol<length(split[[1]]))
    {
      data.maxCol <- length(split[[1]])
    }
    if(lat>90 || lat< -90)
    {
      print(paste0("Latitude in the row", i , "is invalid."))
    }
    if(lon>180 || lon< -180)
    {
      print(paste0("Longitude in the row", i , "is invalid."))
    }
    data.loc[i,1]<-lon
    data.loc[i,2]<-lat
  }

  #remove columns for latlon and assume all sequences are different
  data.max<-(data.maxCol-2)*length(data)
  alleleCounter<-1
  data.alleleCount<-numeric(data.max)
  data.alleleSequence<-character(data.max)
  locationAlleleList<-vector("list",length(data))
  for(i in 1:length(data))
  {
    split<-strsplit(data[i],",")
    for(j in 3:length(split[[1]]))
    {
      sequence<-split[[1]][j]
      if(!sequence=="")
      {
        locationAlleleList[[i]]<-c(locationAlleleList[[i]],sequence)
        #when a new sequence is found
        index <-match(sequence,data.alleleSequence)
        if(is.na(index))
        {
          data.alleleCount[alleleCounter] <- 1
          data.alleleSequence[alleleCounter]<-sequence
          alleleCounter <- alleleCounter + 1
        }
        #when an already found sequence is checked
        else
        {
          data.alleleCount[index] <- data.alleleCount[index] + 1
        }
      }
    }
    locationAlleleList[[i]]<-unique(locationAlleleList[[i]])
  }
  #number of all alleles
  alternative<-(alleleCounter-1)
  #sample size of each allele
  alleleCount<-numeric(alternative)
  #sequence of each allele
  alleleSequence<-character(alternative)

  for(i in 1:alleleCounter-1)
  {
    alleleCount[i] <- data.alleleCount[i]
    alleleSequence[i]<- data.alleleSequence[i]
  }
  #detect location-specific (unique) allele
  isLocUniqueAllele<-vector("logical",length(alleleSequence))
  for(j in 1:length(alleleSequence))
  {
    isFoundInALoc<-FALSE
    for(i in 1:length(data))
    {
      index1<- match(alleleSequence[j],locationAlleleList[[i]])
      #when an allele is found in two locations, it is not location-specific, and check the next allele
      if(!is.na(index1)&&isFoundInALoc)
      {
        isFoundInALoc<-FALSE
        break
      }
      #when an allele is found first, record it as unique
      if(!is.na(index1)&&!isFoundInALoc)
      {
        isFoundInALoc<-TRUE
      }
    }
    if(isFoundInALoc)
    {
      isLocUniqueAllele[j]<-TRUE
    }
  }
  #find location-specific (unique) allele in each location
  uniqueAlleleIndex<-vector("list",length(data))
  for(j in 1:length(alleleSequence))
  {
    if(!isLocUniqueAllele[j])
    {
      next
    }
    for(i in 1:length(data))
    {
      index1<- match(alleleSequence[j],locationAlleleList[[i]])
      if(!is.na(index1))
      {
        uniqueAlleleIndex[[i]]<-c(uniqueAlleleIndex[[i]],j)
      }
    }
  }
  return(list(data.loc,alleleCount,alleleSequence,alternative,uniqueAlleleIndex))
}

sf2SPDF<-function(shape)
{
  dataSize<-length(shape[[1]])
  spdf<-as(shape,"Spatial")
  spdf$id<-seq(1,dataSize)
  return(spdf)
}

getIDsamples<-function(data.loc,SPDF,usest=FALSE)
{
  spatialpoint<-SpatialPoints(data.loc,proj4string=CRS(proj4string(SPDF)))
  sfPolygons<-st_as_sf(SPDF)
  points<-st_as_sf(spatialpoint)
  st_agr(points)<-"constant"
  st_agr(sfPolygons)<-"constant"
  if(usest)
  {
    intersect<-st_intersection(sfPolygons,points)
    #remove duplicated points which are on multiple meshes
    return(intersect$id[!duplicated(intersect$geometry)])
  }
  else
  {
    for(i in 1:length(spatialpoint))
    {
      temp<-s2_intersection(sfPolygons,paste0("POINT (",data.loc[i,1]," ",data.loc[i,2],")"),options=s2_options(model="closed"))
      temp<-match(FALSE,s2_is_empty(temp))
      if(i!=1)
      {
        #use only the first element of temp to exclude the case a point crosses multiple meshes
        intersect<-c(intersect,temp[1])
      }else
      {
        intersect<-temp[1]
      }
    }
    return(intersect)
  }
}
#function for modeling
makeDataMatrix<-function(data,alternative,alleleSequence,id.samples,estimatingLocNo)
{
  totalRow<-estimatingLocNo*alternative
  allY<-rep(NA,totalRow)
  for(i in 1:length(data))
  {
    #split string
    splitAll<-strsplit(data[i],",")[[1]]
    #three and later elements are allele sequence
    sequenceBase<-splitAll[3:length(splitAll)]
    #get id of the location
    currentLocID<-id.samples[i]
    #count the allele
    count<-sapply(alleleSequence,function(seq)sum(sequenceBase==seq,na.rm=TRUE))
    startIndex<-(currentLocID-1)*alternative+1
    endIndex<-currentLocID*alternative
    allY[startIndex:endIndex]<-rowSums(cbind(allY[startIndex:endIndex],count),na.rm=TRUE)
  }
  #allele:(1,1,1,2,2,2,...,100,100,100)
  alleleData<-rep(1:estimatingLocNo,each=alternative)
  #phi:(1,1,1,2,2,2,...,100,100,100)
  phiData<-rep(1:estimatingLocNo,each=alternative)
  #alt:(1,2,3,4,5,6,7,1,2,3,...)
  altData<-rep(1:alternative,times=estimatingLocNo)
  dataMatrix<-data.frame(Y=allY,allele=alleleData,phi=phiData,alt=altData)
  return(dataMatrix)
}
calcCAR<-function(alternative,randomModel,adjacency,dataMatrix,method="group",prior="jeffreystdf",verbose=FALSE)
{
  if(method!="group"&&method!="replicate")
  {
    print("method must be 'group' or 'replicate'.")
    return()
  }
  if(method=="group")
  {
    if(randomModel!="besag"&&randomModel!="besag2")
    {
      formula<-Y~-1+f(phi,hyper=list(prec=list(fixed=TRUE,initial=-10)))+f(alt,constr=TRUE)+f(allele,model=randomModel,graph=adjacency,group=alt,control.group=list(model="iid"))
    }
    else
    {
      formula<-Y~-1+f(phi,hyper=list(prec=list(fixed=TRUE,initial=-10)))+f(alt,constr=TRUE)+f(allele,model=randomModel,graph=adjacency,group=alt,control.group=list(model="iid"),scale.model=TRUE)
    }
  }
  else if(method=="replicate")
  {
    if(randomModel!="besag"&&randomModel!="besag2")
    {
      formula<-Y~-1+f(phi,hyper=list(prec=list(fixed=TRUE,initial=-10)))+f(alt,constr=TRUE)+f(allele,model=randomModel,graph=adjacency,replicate=alt)
    }
    else
    {
      formula<-Y~-1+f(phi,hyper=list(prec=list(fixed=TRUE,initial=-10)))+f(alt,constr=TRUE)+f(allele,model=randomModel,graph=adjacency,replicate=alt,scale.model=TRUE)
    }
  }
  model<-inla(formula,family="poisson",data=dataMatrix,control.compute=list(config=TRUE,cpo=TRUE),control.predictor=list(compute=TRUE),control.inla=list(int.strategy="eb"),verbose=verbose)
  return(model)
}
#resample allele frequency from posterior distribution
resampleFreq<-function(model,n)
{
  alleleCount<-length(model$summary.random)-1
  resampleData<-inla.posterior.sample(n,model)
  alleleSampleList<-vector("list",alleleCount)
  index<-grep(paste0("^allele",1,":"),attr(resampleData[[1]][[2]][,1],"names"))
  meshCount<-length(index)
  #denominator to divide raw values by the sum of allele data
  denominator<-matrix(0,nrow=meshCount,ncol=n)
  #rest of rounded-off
  rest<-matrix(0,nrow=meshCount,ncol=1)
  for(k in 1:alleleCount)
  {
    alleleSampleList[[k]]<-matrix(0,nrow=meshCount,ncol=n)
    index<-grep(paste0("^allele",k,":"),attr(resampleData[[1]][[2]][,1],"names"))
    for(i in 1:n)
    {
      alleleSampleList[[k]][,i]<-exp(resampleData[[i]][[2]][index])
      #calculate sum considering loss of trailing digits
      temp<-denominator[,i]+alleleSampleList[[k]][,i]+rest[,1]
      rest[,1]<-(alleleSampleList[[k]][,i]+rest[,1])-(temp-denominator[,i])
      denominator[,i]<-temp
      #naive calculation:
      #denominator[,i]<-denominator[,i]+alleleSampleList[[k]][,i]
    }
  }
  for(k in 1:alleleCount)
  {
    #divide raw values by the sum of allele data
    alleleSampleList[[k]]<-alleleSampleList[[k]]/denominator
  }
  return(alleleSampleList)
}
#obtain mean, median and quantiles from the result of SPDE estimation
getMeanQtSPDE<-function(model,mesh,data.loc,n=1001,q=c(0.5,0.025,0.975))
{
  samples <- inla.posterior.sample(n, model)
  contents <- model$misc$configs$contents
  #spatial.field index
  startIndex <- contents$start[contents$tag=="spatial.field"]
  indexLength <- contents$length[contents$tag=="spatial.field"]
  fieldIndex <- startIndex:(startIndex+indexLength-1)
  #alt index
  startIndexAlt <- contents$start[contents$tag=="alt"]
  indexLengthAlt <- contents$length[contents$tag=="alt"]
  fieldIndexAlt <- startIndexAlt:(startIndexAlt+indexLengthAlt-1)

  alternative <- indexLength / mesh$n
  locNo <- nrow(data.loc)
  meshProj <- fm_evaluator(mesh, loc=data.loc)
  probSamples <- array(NA, dim=c(locNo, alternative, n))
  for(s in 1:n)
  {
    rhoSample <- samples[[s]]$latent[fieldIndex]
    rhoMatrix <- matrix(rhoSample,nrow=mesh$n,ncol=alternative,byrow=FALSE)
    altSample <- samples[[s]]$latent[fieldIndexAlt]
    #eta(per mesh)
    etaMesh <- sweep(rhoMatrix, 2, altSample, "+")
    #project to locations
    etaLoc <- fm_evaluate(meshProj, field=etaMesh)
    #softmax (per location)
    probLoc <- matrix(NA, nrow=locNo, ncol=alternative)
    for(i in 1:locNo)
    {
      maxEta <- max(etaLoc[i,])
      tmp <- exp(etaLoc[i,] - maxEta)
      probLoc[i,] <- tmp / sum(tmp)
    }
    probSamples[,,s] <- probLoc
  }
  meanData <- matrix(NA, locNo, alternative)
  midCI    <- matrix(NA, locNo, alternative)
  lowerCI  <- matrix(NA, locNo, alternative)
  higherCI <- matrix(NA, locNo, alternative)
  for(i in 1:locNo)
  {
    for(j in 1:alternative)
    {
      x <- probSamples[i,j,]
      meanData[i,j] <- mean(x, na.rm=TRUE)
      midCI[i,j]    <- quantile(x, probs=q[1], na.rm=TRUE, names=FALSE)
      lowerCI[i,j]  <- quantile(x, probs=q[2], na.rm=TRUE, names=FALSE)
      higherCI[i,j] <- quantile(x, probs=q[3], na.rm=TRUE, names=FALSE)
    }
  }
  return(list(mean=meanData,midCI=midCI,lowerCI=lowerCI,higherCI=higherCI))
}

#obtain frequency mean and quantiles from the result of (I)CAR estimation
getMeanQt<- function(model, dataMatrix, n=1001, q=c(0.5,0.025,0.975))
{
  sample <- inla.posterior.sample(n, model)
  contents <- model$misc$configs$contents
  startAlt <- contents$start[contents$tag=="alt"]
  lenAlt <- contents$length[contents$tag=="alt"]
  idxAlt <- startAlt:(startAlt + lenAlt - 1)
  startAllele <- contents$start[contents$tag=="allele"]
  lenAllele <- contents$length[contents$tag=="allele"]
  idxAllele <- startAllele:(startAllele + lenAllele - 1)
  alternative <- max(dataMatrix$alt)
  nObs <- nrow(dataMatrix)
  nAllelePerAlt <- lenAllele / alternative
  probMatrix <- matrix(NA, nrow=nObs, ncol=n)
  for(s in 1:n)
  {
    altSample <- sample[[s]]$latent[idxAlt]
    alleleSample <- sample[[s]]$latent[idxAllele]
    eta <- numeric(nObs)
    for(r in 1:nObs)
    {
      altNow <- dataMatrix$alt[r]
      alleleNow <- dataMatrix$allele[r]
      alleleIndex <- (altNow - 1) * nAllelePerAlt + alleleNow
      eta[r] <- altSample[altNow] + alleleSample[alleleIndex]
    }
    maxEta <- ave(eta, dataMatrix$phi, FUN=max)
    lambda <- exp(eta - maxEta)
    phiSum <- ave(lambda, dataMatrix$phi, FUN=sum)
    probMatrix[, s] <- lambda / phiSum
  }
  meanData <- matrix(NA, nrow=nObs / alternative, ncol=alternative)
  midCI <- matrix(NA, nrow=nObs / alternative, ncol=alternative)
  lowerCI <- matrix(NA, nrow=nObs / alternative, ncol=alternative)
  higherCI <- matrix(NA, nrow=nObs / alternative, ncol=alternative)
  meanVec <- rowMeans(probMatrix, na.rm=TRUE)
  midVec <- apply(probMatrix, 1, quantile, probs=q[1], na.rm=TRUE, names=FALSE)
  lowVec <- apply(probMatrix, 1, quantile, probs=q[2], na.rm=TRUE, names=FALSE)
  hiVec <- apply(probMatrix, 1, quantile, probs=q[3], na.rm=TRUE, names=FALSE)
  meanData <- matrix(meanVec, ncol=alternative, byrow=TRUE)
  midCI <- matrix(midVec, ncol=alternative, byrow=TRUE)
  lowerCI <- matrix(lowVec, ncol=alternative, byrow=TRUE)
  higherCI <- matrix(hiVec, ncol=alternative, byrow=TRUE)
  return(list(mean=meanData, midCI=midCI, lowerCI=lowerCI, higherCI=higherCI))
}
#obtain quantiles of expected heterozygosity
getHeQt<-function(model,n=1001,q=c(0.5,0.025,0.975))
{
  alleleCount<-length(model$summary.random)-1
  alleleSampleList<-resampleFreq(model,n)
  meshCount<-length(alleleSampleList[[1]][,1])
  #convert the data format
  freq<-vector("list",length(q))
  heList<-vector("list",meshCount)
  for(i in 1:meshCount)
  {
    freq[[i]]<-matrix(0,nrow=meshCount,ncol=alleleCount)
    for(k in 1:alleleCount)
    {
      freq[[i]][,k]<-alleleSampleList[[k]][,i]
    }
    heList[[i]]<-calcExpectedHeterozygosity(freq[[i]])
  }
  #calculate quantile
  result<-matrix(0,nrow=meshCount,ncol=length(q))
  for(i in 1:meshCount)
  {
    result[i,]<-quantile(heList[[i]],q)
  }
  colnames(result)<-q
  return(result)
}
#obtain quantiles of nucleotide diversity
getPiQt<-function(model,alleleSequence,n=1001,q=c(0.5,0.025,0.975))
{
  alleleCount<-length(model$summary.random)-1
  alleleSampleList<-resampleFreq(model,n)
  meshCount<-length(alleleSampleList[[1]][,1])
  #convert the data format
  freq<-vector("list",length(q))
  piList<-vector("list",meshCount)
  for(i in 1:meshCount)
  {
    freq[[i]]<-matrix(0,nrow=meshCount,ncol=alleleCount)
    for(k in 1:alleleCount)
    {
      freq[[i]][,k]<-alleleSampleList[[k]][,i]
    }
    piList[[i]]<-calcNucleotideDiversity(freq[[i]],alleleSequence)
  }
  #calculate quantile
  result<-matrix(0,nrow=meshCount,ncol=length(q))
  for(i in 1:meshCount)
  {
    result[i,]<-quantile(piList[[i]],q)
  }
  colnames(result)<-q
  return(result)
}

plotDatum<-function(polygon,data,title,markerLoc=NULL,tile="Esri.WorldGrayCanvas",col1="red",col2="yellow")
{
  lf<-leaflet()
  lf<-addProviderTiles(lf,tile)
  if(class(polygon)[1]=="fm_mesh_2d")
  {
    meshProj<-fm_evaluator(polygon)
    col<-colorNumeric(colorRamp(c(col1,col2)),domain=c(min(fm_evaluate(meshProj,field=data),na.rm=TRUE),max(fm_evaluate(meshProj,field=data),na.rm=TRUE)),na.color="#FF000000")
    rect<-meshProj2Rectangles(meshProj)
    lf<-addRectangles(lf,lng1=rect[[1]],lat1=rect[[2]],lng2=rect[[3]],lat2=rect[[4]],color=col(fm_evaluate(meshProj,field=data)),stroke=F,fillOpacity=0.5)
    lf<-addLegend(lf,pal=col,values=fm_evaluate(meshProj,field=data), title=title)
  }
  else if(class(polygon)[1]=="SpatialPolygonsDataFrame")
  {
    col<-colorNumeric(colorRamp(c(col1,col2)),domain=range(data))
    lf<-addLegend(lf,pal=col,values=data, title=title)
    lf<-addPolygons(lf,data=polygon,color=col(data),stroke=F,fillOpacity=0.5)
  }
  else
  {
    print("Error:The polygon argument must be a fm_mesh_2d object or SpatialPolygonsDataFrame object.")
  }
  if(!is.null(markerLoc))
  {
    lf<-addMarkers(lf,data=markerLoc,group="sample")
  }
  lf<-addLayersControl(lf,overlayGroups=c("sample"),options=layersControlOptions(collapsed=FALSE))
  lf<-addScaleBar(lf,"topleft")
  return(lf)
}

plotAllFreq<-function(polygon,data,groupList=NULL,dimx=10,dimy=10,type="pie",size=30,tile="Esri.WorldGrayCanvas",color=NULL)
{
  if(is.null(color))
  {
    color<-paste0(rainbow(length(data[1,])),"C0")
  }
  #group multiple alleles into a group
  if(!is.null(groupList))
  {
    #check groupList
    check<-unlist(groupList)
    if(length(check)!=length(unique(check)))
    {
      print("Error: One or more allele(s) is included in multiple groupList.")
      return()
    }
    if(length(check)<length(data[1,]))
    {
      dif<-length(data[1,])-length(check)
      print(paste0("Warning: ",dif," allele(s) is not included in any groupList."))
    }
    temp<-matrix(0,length(data[,1]),length(groupList))
    for(i in 1:length(groupList))
    {
      if(length(groupList[[i]])>1)
      {
        temp[,i]<-rowSums(data[,groupList[[i]]])
      }
      #when checking a single row, just pass the data
      else
      {
        temp[,i]<-data[,groupList[[i]]]
      }
    }
    data<-temp
  }
  lf<-leaflet()
  lf<-addProviderTiles(lf,tile)
  if(class(polygon)[1]=="fm_mesh_2d")
  {
    meshProj<-fm_evaluator(polygon,dims=c(dimx,dimy))
    rect<-meshProj2Rectangles(meshProj)
    alleleCount<-length(data[1,])

    for(i in 1:length(rect[[1]]))
    {
      lat<-(rect[[2]][i]+rect[[4]][i])/2
      #check whether the rectangle crosses -180=180 long
      if(rect[[1]][i]*rect[[3]][i]<0 && abs(rect[[1]][i])>90)
      {
        lon<-(rect[[1]][i]+rect[[3]][i])/2-180
      }
      else
      {
        lon<-(rect[[1]][i]+rect[[3]][i])/2
      }
      probability<-numeric(alleleCount)
      for(j in 1:alleleCount)
      {
        probability[j]<-fm_evaluate(meshProj,field=data[,j])[i]
      }
      lf<-addMinicharts(lf,lng=lon,lat=lat,type=type,chartdata=probability,width=size,height=size,colorPalette=color)
    }
  }
  else if(class(polygon)[1]=="SpatialPolygonsDataFrame")
  {
    for(i in 1:length(polygon))
    {
      lon<-coordinates(polygon)[i,1]
      lat<-coordinates(polygon)[i,2]
      lf<-addMinicharts(lf,lng=lon,lat=lat,type=type,chartdata=data[i,],width=size,height=size,colorPalette=color)
    }
  }
  else
  {
    print("Error:The polygon argument must be a fm_mesh_2d object or SpatialPolygonsDataFrame object.")
  }
  return(lf)
}

makeSampleFreq<-function(dataMatrix,SPDF)
{
  locSize<-length(SPDF)
  alleleSize<-max(dataMatrix$alt)
  prob<-matrix(0,locSize,alleleSize)
  #count alleles for each location
  for(i in 1:nrow(dataMatrix))
  {
    if(!is.na(dataMatrix[i,1]))
    {
      if(dataMatrix[i,1]!=0)
      {
        Y<-dataMatrix[i,1]
        alleleNo<-dataMatrix$alt[i]
        locNo<-dataMatrix$allele[i]
        prob[locNo,alleleNo]<-prob[locNo,alleleNo]+Y
      }
    }
  }
  #convert the count into probability
  for(i in 1:locSize)
  {
    locSum<-sum(prob[i,])
    if(locSum!=0)
    {
      prob[i,]<-prob[i,]/locSum
    }
    else
    {
      prob[i,]<-NaN
    }
  }
  return(prob)
}

meshProj2Rectangles<-function(meshProjector)
{
  data.x<-meshProjector$x
  data.y<-meshProjector$y
  xDif<-meshProjector$x[2]-meshProjector$x[1]
  yDif<-meshProjector$y[2]-meshProjector$y[1]
  #Note: data.x and data.y are always longer than 1.

  lng1<-numeric(length(data.x)*length(data.y))
  lat1<-numeric(length(data.x)*length(data.y))
  lng2<-numeric(length(data.x)*length(data.y))
  lat2<-numeric(length(data.x)*length(data.y))
  counter<-1
  for(i in 1:length(data.y))
  {
    for(j in 1:length(data.x))
    {
      lng1[counter]<-data.x[j]-xDif/2
      lng2[counter]<-data.x[j]+xDif/2
      lat1[counter]<-data.y[i]-yDif/2
      lat2[counter]<-data.y[i]+yDif/2
      #adjust lng over 180 and -180
      if(lng1[counter]>180)
      {
        lng1[counter]<-180
      }
      if(lng2[counter]>180)
      {
        lng2[counter]<-180
      }
      #space between < and - is necessary. Otherwise, it is interpreted as <- (substitution).
      if(lng1[counter]< -180)
      {
        lng1[counter]<--180
      }
      if(lng2[counter]< -180)
      {
        lng2[counter]<--180
      }
      #adjust lat over 90 and -90
      if(lat1[counter]>90)
      {
        lat1[counter]<-90
      }
      if(lat2[counter]>90)
      {
        lat2[counter]<-90
      }
      if(lat1[counter]< -90)
      {
        lat1[counter]<--90
      }
      if(lat2[counter]< -90)
      {
        lat2[counter]<--90
      }
      counter<-counter+1
    }
  }
  return(list(lng1,lat1,lng2,lat2))
}
calcSPDE<-function(alternative,dataMatrix,mesh,data.loc,verbose=FALSE,method="group",spdeType="pc",sigma0=1,rangeDenom=5,range=c(NA,0.01),sigma=c(3,0.01))
{
  if(method!="group" && method!="replicate")
  {
    stop("method must be 'group' or 'replicate'")
  }
  #SPDE definition
  if(spdeType=="matern")
  {
    size <- min(c(diff(range(mesh$loc[,1])), diff(range(mesh$loc[,2]))))
    range0 <- size / rangeDenom
    kappa0 <- sqrt(8) / range0
    tau0 <- 1 / (sqrt(4*pi) * kappa0 * sigma0)
    spde <- inla.spde2.matern(mesh,B.tau=cbind(log(tau0),-1,1),B.kappa=cbind(log(kappa0),0,-1))
  }
  else if(spdeType=="pc")
  {
    meshSize <- max(c(diff(range(mesh$loc[,1])), diff(range(mesh$loc[,2]))))
    if(is.na(range[1]))
    {
      range[1] <- meshSize * 0.2
    }
    spde <- inla.spde2.pcmatern(mesh,prior.range=range,prior.sigma=sigma)
  }
  else
  {
    stop("spdeType must be 'matern' or 'pc'")
  }
  coords<-cbind(rep(data.loc[,1], each=alternative),rep(data.loc[,2], each=alternative))
  if(method=="group")
  {
    s_index <- inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde,n.group=alternative)
    A <- inla.spde.make.A(mesh=mesh,loc=coords,group=dataMatrix$alt,n.group=alternative)
  }
  else
  {
    s_index <- inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde,n.repl=alternative)
    A <- inla.spde.make.A(mesh=mesh,loc=coords,replicate=dataMatrix$alt,n.repl=alternative)
  }
  stack <- inla.stack(tag='estimation',data=list(Y=dataMatrix$Y),A=list(A,1),effects=list(s_index,list(phi=dataMatrix$phi, alt=dataMatrix$alt)))
  if(method=="group")
  {
    formula <- Y~-1+f(phi,hyper=list(prec=list(fixed=TRUE,initial=-10)))+f(alt,constr=TRUE)+f(spatial.field, model=spde,group=spatial.field.group,control.group=list(model="iid"))
  }
  else
  {
    formula<-Y~-1+f(phi,hyper=list(prec=list(fixed=TRUE,initial=-10)))+f(alt,constr=TRUE)+f(spatial.field, model=spde,replicate=spatial.field.repl)
  }
  model <- inla(formula,data=inla.stack.data(stack),family="poisson",control.predictor=list(A=inla.stack.A(stack), compute=TRUE),control.compute=list(config=TRUE),control.inla=list(int.strategy="eb"),verbose=verbose)
  return(model)
}

#calculate HTAE: Halved Total Absolute Error of allele frequencies. HTAE ranges [0,1).
calcHTAE<-function(estimatedValue,trueValue)
{
  result<-0
  if(length(estimatedValue)!=length(trueValue))
  {
    print("Error: estimatedValue and trueValue must be the same length.")
    return(0)
  }
  for(i in 1:length(estimatedValue))
  {
    result<-result + abs(estimatedValue[i] - trueValue[i])
  }
  return(result/2)
}
#calculate log loss
calcLogLoss<-function(estimatedFreq,trueCount)
{
  result<-0
  if(length(estimatedFreq)!=length(trueCount))
  {
    print("Error: estimatedFreq and trueCount must be the same length.")
    return(0)
  }
  for(i in 1:length(estimatedFreq))
  {
    if(estimatedFreq[i]==0)
    {
      next
    }
    result<-result-trueCount[i]*log(estimatedFreq[i])
  }
  return(result)
}
#make data where one location is omitted for cross validation
leaveOneOut<-function(outIndex,data)
{
  return(data[-1*outIndex])
}
#compare the frequency estimated by leave-one-out data and the sample frequency for cross validation
validateFrequency<-function(locIndex,id.samples,outData,sampleData,method="both")
{
  meshID<-id.samples[locIndex]
  estimated<-outData[meshID,]
  sample<-sampleData[meshID,]
  if(method=="logloss")
  {
    result<-calcLogLoss(estimated,sample)
  }
  else if(method=="HTAE")
  {
    result<-calcHTAE(estimated,sample)
  }
  else
  {
    result<-"Error: invalid method is designated."
    print(result)
  }
  return(result)
}
#conduct leave-one-out cross validation of CAR-model-based interpolation
cvCAR<-function(data,randomModel,SPDF,data.loc,adjacency,allId.samples,uniqueAlleleIndex,sampleProb,method="group",prior="jeffreystdf")
{
  loss<-vector("list",2)
  loss[[1]]<-vector("numeric",length(data))
  loss[[2]]<-vector("numeric",length(data))
  for(i in 1:length(data))
  {
    print(paste0("location ",i,"/",length(data)))
    newData<-leaveOneOut(i,data)
    temp<-processInput(newData)
    data.loc<-temp[[1]]
    alleleCount<-temp[[2]]
    alleleSequence<-temp[[3]]
    alternative<-temp[[4]]
    id.samples<-getIDsamples(data.loc,SPDF,TRUE)
    dataMatrix<-makeDataMatrix(newData,alternative,alleleSequence,id.samples,length(SPDF))
    models<-calcCAR(alleleSequence,randomModel,graph,dataMatrix,method=method,prior=prior)
    qt<-getMeanQt(models,dataMatrix)
    loss[[1]][i]<-validateFrequency(i,allId.samples,qt[[1]],sampleProb,"HTAE")
    #skip if the location contains location-specific allele
    if(!is.null(uniqueAlleleIndex[[i]]))
    {
      loss[[2]][i]<-NA
      next
    }
    loss[[2]][i]<-validateFrequency(i,allId.samples,qt[[1]],sampleProb,"logloss")
  }
  return(loss)
}
#conduct leave-one-out cross validation of SPDE interpolation
cvSPDE<-function(data,SPDF,estimate.loc,mesh,allId.samples,uniqueAlleleIndex,sampleProb,method="group",spdeType="pc",sigma0=1,rangeDenom=5,range=c(NA,0.01),sigma=c(3,0.01))
{
  loss<-vector("list",2)
  loss[[1]]<-vector("numeric",length(data))
  loss[[2]]<-vector("numeric",length(data))
  for(i in 1:length(data))
  {
    print(paste0("location ",i,"/",length(data)))
    newData<-leaveOneOut(i,data)
    temp<-processInput(newData)
    data.loc<-temp[[1]]
    alleleCount<-temp[[2]]
    alleleSequence<-temp[[3]]
    alternative<-temp[[4]]
    id.samples<-getIDsamples(data.loc,SPDF,TRUE)
    dataMatrix<-makeDataMatrix(newData,alternative,alleleSequence,id.samples,length(SPDF))
    models<-calcSPDE(alternative,dataMatrix,mesh=mesh,data.loc=estimate.loc,method=method,spdeType=spdeType,sigma0=sigma0,range=range,sigma=sigma)
    qt<-getMeanQtSPDE(models,mesh,estimate.loc)
    loss[[1]][i]<-validateFrequency(i,allId.samples,qt[[1]],sampleProb,"HTAE")
    #skip if the location contains location-specific allele
    if(!is.null(uniqueAlleleIndex[[i]]))
    {
      loss[[2]][i]<-NA
      next
    }
    loss[[2]][i]<-validateFrequency(i,allId.samples,qt[[1]],sampleProb,"logloss")
  }
  return(loss)
}
#make a SPDF data for IDW
getAlleleSPDF<-function(data)
{
  #assuming there is no duplicated location data.
  data.loc<-matrix(0,nrow=length(data),ncol=2)
  data.maxCol<-0
  for(i in 1:length(data))
  {
    split<-strsplit(data[i],",")
    lon<-parse_lon(split[[1]][1])
    lat<-parse_lat(split[[1]][2])
    if(length(split[[1]])<3)
    {
      print(paste0("Location in row ", i , " has no allele."))
    }
    if(data.maxCol<length(split[[1]]))
    {
      data.maxCol <- length(split[[1]])
    }
    if(lat>90 || lat< -90)
    {
      print(paste0("Latitude in the row", i , "is invalid."))
    }
    if(lon>180 || lon< -180)
    {
      print(paste0("Longitude in the row", i , "is invalid."))
    }
    data.loc[i,1]<-lon
    data.loc[i,2]<-lat
  }

  #remove columns for latlon and assume all sequences are different
  data.max<-(data.maxCol-2)*length(data)
  alleleCounter<-1
  data.alleleCount<-numeric(data.max)
  data.alleleSequence<-character(data.max)

  for(i in 1:length(data))
  {
    split<-strsplit(data[i],",")
    for(j in 3:length(split[[1]]))
    {
      sequence<-split[[1]][j]
      if(!sequence=="")
      {
        #when a new sequence is found
        index <-match(sequence,data.alleleSequence)
        if(is.na(index))
        {
          data.alleleCount[alleleCounter] <- 1
          data.alleleSequence[alleleCounter]<-sequence
          alleleCounter <- alleleCounter + 1
        }
        #when an already found sequence is checked
        else
        {
          data.alleleCount[index] <- data.alleleCount[index] + 1
        }
      }
    }
  }
  #number of all potential alleles
  alternative<-(alleleCounter-1)
  #sample size of each allele
  alleleCount<-numeric(alternative)
  #sequence of each allele
  alleleSequence<-character(alternative)

  for(i in 1:alleleCounter-1)
  {
    alleleCount[i] <- data.alleleCount[i]
    alleleSequence[i]<- data.alleleSequence[i]
  }

  #calculate allele frequency for each point
  alleleFreq<-vector("list",length(data))
  for(i in 1:length(data))
  {
    alleleFreq[[i]]<-numeric(length(alleleCount))
    totalAllele<-0
    split<-strsplit(data[i],",")
    for(j in 3:length(split[[1]]))
    {
      sequence<-split[[1]][j]
      if(!sequence=="")
      {
        index<-match(sequence,data.alleleSequence)
        alleleFreq[[i]][index]<-alleleFreq[[i]][index]+1
        totalAllele<-totalAllele+1
      }
    }
    alleleFreq[[i]]<-alleleFreq[[i]]/totalAllele
  }
  #make spatial points dataframe
  sp<-SpatialPoints(data.loc)
  d<-data.frame(t(data.frame(alleleFreq)),row.names=1:length(data))
  spd<-SpatialPointsDataFrame(sp,data=d)
  return(spd)
}

#make a SPDF data for sample of IDW
getPopulationAlleleSPDF<-function(data,pdata)
{
  #assuming there is no duplicated location data.
  data.loc<-matrix(0,nrow=length(data),ncol=2)
  data.maxCol<-0
  for(i in 1:length(data))
  {
    split<-strsplit(data[i],",")
    lon<-parse_lon(split[[1]][1])
    lat<-parse_lat(split[[1]][2])
    if(length(split[[1]])<3)
    {
      print(paste0("Location in row ", i , " has no allele."))
    }
    if(data.maxCol<length(split[[1]]))
    {
      data.maxCol <- length(split[[1]])
    }
    if(lat>90 || lat< -90)
    {
      print(paste0("Latitude in the row", i , "is invalid."))
    }
    if(lon>180 || lon< -180)
    {
      print(paste0("Longitude in the row", i , "is invalid."))
    }
    data.loc[i,1]<-lon
    data.loc[i,2]<-lat
  }

  #assuming there is no duplicated location data.
  pdata.loc<-matrix(0,nrow=length(pdata),ncol=2)
  pdata.maxCol<-0
  for(i in 1:length(pdata))
  {
    split<-strsplit(pdata[i],",")
    lon<-parse_lon(split[[1]][1])
    lat<-parse_lat(split[[1]][2])
    if(length(split[[1]])<3)
    {
      print(paste0("Location in row ", i , " has no allele."))
    }
    if(pdata.maxCol<length(split[[1]]))
    {
      pdata.maxCol <- length(split[[1]])
    }
    if(lat>90 || lat< -90)
    {
      print(paste0("Latitude in the row", i , "is invalid in population data."))
    }
    if(lon>180 || lon< -180)
    {
      print(paste0("Longitude in the row", i , "is invalid in population data."))
    }
    pdata.loc[i,1]<-lon
    pdata.loc[i,2]<-lat
  }

  #remove columns for latlon and assume all sequences are different
  data.max<-(data.maxCol-2)*length(pdata)
  alleleCounter<-1
  data.alleleCount<-numeric(data.max)
  data.alleleSequence<-character(data.max)
  for(i in 1:length(pdata))
  {
    psplit<-strsplit(pdata[i],",")
    #create a column based on the population data
    for(j in 3:length(psplit[[1]]))
    {
      sequence<-psplit[[1]][j]
      if(!sequence=="")
      {
        #when a new sequence is found
        index <-match(sequence,data.alleleSequence)
        if(is.na(index))
        {
          data.alleleSequence[alleleCounter]<-sequence
          alleleCounter <- alleleCounter + 1
        }
      }
    }
  }
  for(i in 1:length(data))
  {
    split<-strsplit(data[i],",")
    #input data based on sample data
    for(j in 3:length(split[[1]]))
    {
      sequence<-split[[1]][j]
      if(!sequence=="")
      {
        index <-match(sequence,data.alleleSequence)
        data.alleleCount[index] <- data.alleleCount[index] + 1
      }
    }
  }
  #number of all potential alleles
  alternative<-(alleleCounter-1)
  #sample size of each allele
  alleleCount<-numeric(alternative)
  #sequence of each allele
  alleleSequence<-character(alternative)
  for(i in 1:alleleCounter-1)
  {
    alleleCount[i] <- data.alleleCount[i]
    alleleSequence[i]<- data.alleleSequence[i]
  }
  #calculate allele frequency for each point
  alleleFreq<-vector("list",length(data))
  for(i in 1:length(data))
  {
    alleleFreq[[i]]<-numeric(length(alleleCount))
    totalAllele<-0
    split<-strsplit(data[i],",")
    for(j in 3:length(split[[1]]))
    {
      sequence<-split[[1]][j]
      if(!sequence=="")
      {
        index<-match(sequence,data.alleleSequence)
        alleleFreq[[i]][index]<-alleleFreq[[i]][index]+1
        totalAllele<-totalAllele+1
      }
    }
    alleleFreq[[i]]<-alleleFreq[[i]]/totalAllele
  }
  #make spatial points dataframe
  sp<-SpatialPoints(data.loc)
  d<-data.frame(t(data.frame(alleleFreq)),row.names=1:length(data))
  spd<-SpatialPointsDataFrame(sp,data=d)
  return(spd)
}
#convert SPDF into center cooridinates for SPDE and IDW grid
polygons2centerpoints<-function(spdf)
{
  polygonCount<-length(spdf@polygons)
  coordinates<-matrix(0,ncol=2,nrow=polygonCount)
  for(i in 1:polygonCount)
  {
    points<-spdf@polygons[[i]]@Polygons[[1]]@coords
    #when the polygon does not cross 180 degree meridian
    if(points[1,1]<points[3,1])
    {
      coordinates[i,]<-(points[1,]+points[2,]+points[3,]+points[4,])/4
    }
    #when the polygon crosses 180 degree meridian.
    else
    {
      eastDeviation<-180 - abs(points[1,1]) + 180 - abs(points[4,1])
      westDeviation<-180 - abs(points[2,1]) + 180 - abs(points[3,1])
      totalDeviation<-(eastDeviation - westDeviation) / 4
      if(totalDeviation>=0)
      {
        coordinates[i,1]<- 180 - totalDeviation
      }
      else
      {
        coordinates[i,1]<- -180 - totalDeviation
      }
      coordinates[i,2]<-(points[1,2]+points[2,2]+points[3,2]+points[4,2])/4
    }
  }
  return(coordinates)
}
