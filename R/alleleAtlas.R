#export neighbor network file
makeNeighbor<-function (path, lat, lon, northPole, southPole, loop)
{
  if(is.null(lat) || is.null(lon))
  {
    return("Please specify lat and lon.")
  }
  if(lat <= 0 || lon <= 0)
  {
    return("lat and lon must be positive numbers.")
  }
  out<-file(path,"w")
  n <- lat*lon
  temp <- as.character(n)
  writeLines(temp,out,sep="\n")
  for(i in 1:n)
  {
    neighbourCounter <- 0
    up<-NULL
    #up neighbour
    if(!(i <= lon))
    {
      if(lon == 1)
      {
        up<-i-lon
        neighbourCounter <- neighbourCounter + 1
      }
      else if(lon == 2)
      {
        up<-numeric(2)
        neighbourCounter <- neighbourCounter + 2
        up[1] <- i-lon
        #upright
        if(i %% lon == 0)
        {
          up[2] <- i-(lon-1)-lon
        }
        else
        {
          up[2] <- i+1-lon
        }
      }
      else
      {
        up <-numeric(3)
        up[2] <- i-lon
        neighbourCounter <- neighbourCounter + 1
        #upright
        if(i %% lon == 0 && loop)
        {
          up[3] <- i-(lon-1)-lon
          neighbourCounter <- neighbourCounter + 1
        }
        else if(i %% lon != 0)
        {
          up[3] <- i+1-lon
          neighbourCounter <- neighbourCounter + 1
        }
        #upleft
        if(i %% lon == 1 && loop)
        {
          up[1] <- i+(lon-1)-lon
          neighbourCounter <- neighbourCounter + 1
        }
        else if(i %% lon != 1)
        {
          up[1] <- i-1-lon
          neighbourCounter <- neighbourCounter + 1
        }
      }
    }
    #neighboring to the North pole
    else if(lon > 1 && northPole)
    {
      up <-numeric(lon-1)
      counter<-1
      neighbourCounter <- neighbourCounter + lon - 1
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
    if(lat > 1)
    {
      if(!((lat-1)*lon+1 <= i && i <= lat*lon))
      {
        if(lon == 1)
        {
          down<-i+lon
          neighbourCounter <- neighbourCounter + 1
        }
        else if(lon == 2)
        {
          down <- numeric(2)
          down[1] <- i+lon
          neighbourCounter <- neighbourCounter + 2
          #downright
          if(i %% lon == 0)
          {
            down[2] <- i-(lon-1)+lon
          }
          else
          {
            down[2] <- i+1+lon
          }
        }
        else
        {
          down <- numeric(3)
          down[2] <- i+lon
          neighbourCounter <- neighbourCounter + 1
          #downright
          if(i %% lon == 0 && loop)
          {
            down[3] <- i-(lon-1)+lon
            neighbourCounter <- neighbourCounter + 1
          }
          else if(i %% lon != 0)
          {
            down[3] <- i+1+lon
            neighbourCounter <- neighbourCounter + 1
          }
          #upleft
          if(i %% lon == 1 & loop)
          {
            down[1] <- i+(lon-1)+lon
            neighbourCounter <- neighbourCounter + 1
          }
          else if(i %% lon != 1)
          {
            down[1] <- i-1+lon
            neighbourCounter <- neighbourCounter + 1
          }
        }
      }
      #neighboring to the South pole
      else if(lon > 1 && southPole)
      {
        down <-numeric(lon-1)
        counter<-1
        neighbourCounter <- neighbourCounter + lon - 1
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
    right <- NULL
    left <- NULL
    if(!((i <= lon)&&northPole) && !(( (lat-1)*lon+1 <= i && i <= lat*lon) && southPole))
    {
      if(lon > 2)
      {
        #right neighbor
        if(i %% lon == 0 && loop)
        {
          right <- i-(lon-1)
          neighbourCounter <- neighbourCounter + 1
        }
        else if(i %% lon != 0)
        {
          right <- i+1
          neighbourCounter <- neighbourCounter + 1
        }
        #left neighbor
        if(i %% lon == 1 && loop)
        {
          left <- i+(lon-1)
          neighbourCounter <- neighbourCounter + 1
        }
        else if(i %% lon != 1)
        {
          left <- i-1
          neighbourCounter <- neighbourCounter + 1
        }
      }
      else if(lon == 2)
      {
        if(i %% lon == 0)
        {
          if(loop)
          {
            right <- i-(lon-1)
            neighbourCounter <- neighbourCounter + 1
          }
          left <- i-1
          neighbourCounter <- neighbourCounter + 1
        }
        else if(i %% lon == 1)
        {
          if(loop)
          {
            left <- i+(lon-1)
            neighbourCounter <- neighbourCounter + 1
          }
          right <- i+1
          neighbourCounter <- neighbourCounter + 1
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
        line <- substr(line,1,nchar(line)-1)
      }
    }
    writeLines(line,out,sep="\n")
  }
  close(out)
}

#make latlon window from location data
makeMinWindow<-function(data.loc,northPole, southPole, loop, forceMaxEast=NULL)
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
      westLimit<-min(data.loc[,1])
      eastLimit<-max(data.loc[,1])
    }
    else if(x>180)
    {
      westLimit<--1*max(data.loc[,1])
      eastLimit<-min(data.loc[,1])
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
        westLimit<-min(data.loc[,1])
        eastLimit<-max(data.loc[,1])
      }
      else if(!forceMaxEast)
      {
        westLimit<-max(data.loc[,1])
        eastLimit<-min(data.loc[,1])
      }
    }
  }
  if(northPole)
  {
    northLimit<-90
  }
  else
  {
    northLimit<-max(data.loc[,2])
  }
  if(southPole)
  {
    southLimit<- -90
  }
  else
  {
    southLimit<-min(data.loc[,2])
  }
  win<-c(westLimit,eastLimit,southLimit,northLimit)
  return(win)
}

makeSpatialPolygons<-function(autoWindow,lat,lon,minlon=-180,maxlon=180,minlat=-90,maxlat=90,data.loc=NULL,northPole=NULL,southPole=NULL,loop=NULL,forceMaxEast=NULL)
{
  #generate automatically the minimum mesh window
  if(autoWindow)
  {
    win<-as.owin(makeMinWindow(data.loc,northPole,southPole,loop,forceMaxEast))
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
  row.names(d)<-paste0("g", 1:length(data.poly))
  data.counts<-SpatialPolygonsDataFrame(data.poly, d)
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
      calcResult[i]<-calcResult[i] - freq[i,j]^2
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
  dist<-0
  for(i in 1:length(vec1))
  {
    #skip gap bases
    if(vec1[i]=="-"||vec2[i]=="-")
    {
      next
    }
    dist<- dist + sum(abs(interpreteBase(vec1[i])-interpreteBase(vec2[i]))/2)
  }
  return(dist)
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
      subMatrix[i,j]<-seq.dist(alleleSequence[i],alleleSequence[j])/alleleLength
    }
  }
  calcResult<-numeric(length(freq[,1]))
  #check each mesh
  for(i in 1:length(freq[,1]))
  {
    calcResult[i]<-0
    #check each allele
    n<-length(freq[1,])
    if(n==1)
    {
      next
    }
    for(j in 1:n)
    {
      for(k in 1:n)
      {
        calcResult[i]<-calcResult[i] + freq[i,j] * freq[i,k] * subMatrix[j,k]
      }
    }
    #correction for unbiasedness
    #but interpolated mesh does not have sample size
    #calcResult[i]<-calcResult[i]*n/(n-1)
  }
  return(calcResult)
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
  return(list(data.loc,alleleCount,alleleSequence,alternative))
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
makeDataMatrix<-function(data,alternative,alleleSequence,id.samples)
{
  #alleles per location
  Y<-matrix(0,nrow=length(data),ncol=alternative)
  for(i in 1:length(data))
  {
    split<-strsplit(data[i],",")
    for(j in 1:length(alleleSequence))
    {
      for(k in 3:length(split[[1]]))
      {
        sequence<-split[[1]][k]
        if(!sequence=="")
        {
          #when a new sequence is found
          if(alleleSequence[j]==sequence)
          {
            Y[i,j]<-Y[i,j]+1
          }
        }
      }
    }
  }
  dataMatrix<-matrix(NA,ncol=1,nrow = length(data)*alternative)
  for(i in 1:length(data))
  {
    #raw allele count
    dataMatrix[((i-1)*alternative+1):(i*alternative),1]<-c(Y[i,])
  }
  dataMatrix<-data.frame(dataMatrix)
  names(dataMatrix)<-c('Y')
  #create and add mesh index column
  idx<-matrix(NA,ncol=alternative,nrow=length(data)*alternative)
  for(i in 1:alternative)
  {
    idx[seq(i,length(data)*alternative,by=alternative),i]<-id.samples
  }
  colnames(idx)<-paste0("allele",1:length(alleleSequence))
  dataMatrix<-cbind(dataMatrix,idx)
  return(dataMatrix)
}

calcCARModel<-function(alleleSequence,randomModel,adjacency,dataMatrix,prior="jeffreystdf",verbose=FALSE)
{
  models<-vector("list",length(alleleSequence))
  formulette<-vector("list",length(alleleSequence))
  meshSize<-length(dataMatrix[,1])/length(alleleSequence)
  for(i in 1:(length(alleleSequence)))
  {
    varname<-paste0("allele",i)
    if(randomModel!="besag2")
    {
      formulette[[i]]<-paste("f(",varname,", model=randomModel, hyper=list(theta=list(prior=prior,param=numeric())),graph = adjacency)")
    }
    else
    {
      formulette[[i]]<-paste("f(",varname,", model=randomModel, hyper=list(theta=list(prior=prior,param=numeric())),scale.model=TRUE, graph = adjacency)")
    }
    #separate data
    data<-dataMatrix[0:(meshSize-1)*length(alleleSequence)+i,c(1,i+1)]
    formula<-as.formula(paste("Y~-1+", paste(formulette[[i]],collapse="+")))
    models[[i]]<-inla(formula, family = "poisson", data = data, control.compute=list(config=TRUE,cpo=TRUE),control.predictor = list(compute=TRUE),verbose=verbose)
  }
  return(models)
}

#function to depict results
getFreq<-function(alternative,models,randomModel=NULL)
{
  prob<-list(alternative)
  for(i in 1:alternative)
  {
    prob[[i]]<-models[[i]]$summary.random[[1]][5][[1]]
  }
  meshCount<-length(prob[[1]])
  if(is.null(randomModel)==TRUE)
  {
    prob.matrix<-matrix(0,nrow=meshCount,ncol=alternative)
  }
  else if((randomModel == "bym") || (randomModel == "bym2"))
  {
    prob.matrix<-matrix(0,nrow=meshCount,ncol=alternative)
  }
  else
  {
    prob.matrix<-matrix(0,nrow=meshCount,ncol=alternative)
  }

  for(i in 1:alternative)
  {
    for(j in 1:meshCount)
    {
      prob.matrix[j,i]<-exp(prob[[i]][j])
    }
  }
  if(is.null(randomModel)==TRUE)
  {
    prob.result<-matrix(0,nrow=meshCount,ncol=alternative)
  }
  else if((randomModel=="bym")||(randomModel=="bym2"))
  {
    prob.result<-matrix(0,nrow=meshCount,ncol=alternative)
  }
  else
  {
    prob.result<-matrix(0,nrow=meshCount,ncol=alternative)
  }
  for(i in 1:alternative)
  {
    for(j in 1:meshCount)
    {
      prob.result[j,i]<-prob.matrix[j,i]/sum(prob.matrix[j,])
    }
  }
  return(prob.result)
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
      if(rect[[1]][i]*rect[[3]][i]<0)
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
  alleleSize<-length(dataMatrix[1,])-1
  prob<-matrix(0,locSize,alleleSize)
  #count alleles for each location
  for(i in 1:length(dataMatrix[,1]))
  {
    if(dataMatrix[i,1]!=0)
    {
      Y<-dataMatrix[i,1]
      alleleNo<-i%%(alleleSize)
      if(alleleNo==0)
      {
        alleleNo<-alleleSize
      }
      locNo<-dataMatrix[i,1+alleleNo]
      prob[locNo,alleleNo]<-prob[locNo,alleleNo]+Y
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

calcSPDEModel<-function(alleleSequence,alternative,dataMatrix,prior="flat",verbose=FALSE,mesh,data.loc)
{
  #make spde object
  spde<-inla.spde2.matern(mesh)

  #make indices
  indices<-vector("list",alternative)
  for(i in 1:alternative)
  {
    indices[[i]]<-inla.spde.make.index(name=paste0("allele",i),n.spde=spde$n.spde)
  }

  #make stack object
  X<-rep(0,length(dataMatrix[,1]))
  tempFrame<-data.frame(X)
  for(i in 1:alternative)
  {
    X<-rep(i,length(dataMatrix[,1]))
    tempFrame<-data.frame(cbind(tempFrame,X))
  }
  tempFrame<-tempFrame[,colnames(tempFrame)!="X"]
  SPDEdata<-data.frame(cbind(dataMatrix,tempFrame))
  effectList<-list()
  for(i in 1:(alternative+length(SPDEdata)-1))
  {
    if(i<=alternative)
    {
      effectList[[i]]<-indices[[i]]
    }
    else
    {
      effectList[[i]]<-SPDEdata[,1+i-alternative]
    }
  }
  names(effectList)[(alternative+1):length(effectList)]<-colnames(SPDEdata)[(alternative+2):length(SPDEdata)]

  #make index
  s_index<-inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde)

  #make stack
  Atext<-""
  for(i in 1:1)
  {
    Atext<-paste0(Atext,"A")
  }

  #make model formula
  models<-vector("list",length(alleleSequence))
  formulette<-vector("list",length(alleleSequence))
  for(i in 1:(length(alleleSequence)))
  {
    effectList[[i]]<-effectList[[i]][1]
    ydata<-SPDEdata$Y
    j<-i
    if(i==length(alleleSequence))
    {
      j<-0
    }
    index<-(1:length(ydata))%%13==j
    ydata<-ydata[index]
    tempLoc<-data.loc[ydata>0,]
    ydata<-ydata[ydata>0]
    #make A
    if(length(ydata)>1)
    {
      A<-inla.spde.make.A(mesh=mesh,loc=cbind(tempLoc[,1],tempLoc[,2]))
    }
    #when length(ydata)==1, it is no longer a matrix, but a vector
    else
    {
      A<-inla.spde.make.A(mesh=mesh,loc=cbind(tempLoc[1],tempLoc[2]))
    }
    stack<-(eval(parse(text=paste0("inla.stack(data=list(Y=ydata),effects=effectList[[i]],A=list(",Atext,"),tag='tag')"))))

    varname<-paste0("allele",i)
    formulette[[i]]<-paste("f(",varname,",model=spde,group=s_index$spatial.field.group,hyper=list(theta=list(prior=prior,param=numeric())))")
    formula<-as.formula(paste("Y~-1+", paste(formulette[[i]],collapse="+")))
    models[[i]]<-inla(formula,data=inla.stack.data(stack,spde=spde),family="poisson",control.predictor=list(A=inla.stack.A(stack),compute=TRUE),verbose=verbose)
  }
  return(models)
}
