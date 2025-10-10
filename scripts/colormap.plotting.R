plotHaplotypes <- function(filename)
{
  #filename <- "outfile.rs7215219.ihs.out";
  
  # filename.der <- paste(filename,".der.colormap",sep="")
  # filename.anc <- paste(filename,".anc.colormap",sep="")
  # data<-as.matrix(read.table(filename.der))
  rawData<-read.table(filename)


  print("Plotting EHH decay...")
  setEPS()
  postscript(file=paste(filename,".ehh.eps",sep=""))
  plot(rawData$V1/1000000,rawData$V2,pch=" ",xlab="Distance from locus (Mb)",ylab="EHH",ylim=c(0,1))
  lines(rawData$V1/1000000,rawData$V3,col="red")
  lines(rawData$V1/1000000,rawData$V4,col="blue")
  legend("topright",legend=c("Derived","Ancestral"),col=c("red","blue"),lty=1)
  dev.off()
  
  numHaps <- dim(data)[1]
  numSites <- dim(data)[2]
  numCols<-range(data)[2]-range(data)[1]
  
  pos <- rawData$V1
  hap <- seq(1,numHaps)

  numSortCols <- min(numCols,5)
  
  ordering <- as.data.frame(matrix(rep(0,numSortCols*numHaps),nrow=numHaps,ncol=numSortCols))

  print("Sorting derived colors...")
  
  for (h in 1:numSortCols)
    {
      for (i in 1:numHaps)
        {
          for (j in 1:numSites)
            {
              if ( data[i,j] == h-1 )
                {
                  ordering[i,h] <- ordering[i,h]+1
                }
            }
        }
    }
  
  order.string <- "sorted.order<-order("
  for (h in 1:numSortCols)
    {
      order.string <- paste(order.string,"ordering$V",h,",",sep="")
    }
  substr(order.string,nchar(order.string),nchar(order.string)) <- ")"
  eval(parse(text = order.string))
  
  blank<-rep("",length(pos))
  
  
  data2<-as.matrix(read.table(filename.anc))
  numHaps2 <- dim(data2)[1]
  numSites2 <- dim(data2)[2]
  numCols2 <- range(data2)[2]-range(data2)[1]
  
  pos2 <- rawData$V1
  hap2 <- seq(1,numHaps2)

  numSortCols2 <- min(numCols2,5)
  
  ordering2 <- as.data.frame(matrix(rep(0,numSortCols2*numHaps2),nrow=numHaps2,ncol=numSortCols2))

  print("Sorting ancestral colors...")
  
  for (h in 1:numSortCols2)
    {
      for (i in 1:numHaps2)
        {
          for (j in 1:numSites2)
            {
              if ( data2[i,j] == h-1 )
                {
                  ordering2[i,h] <- ordering2[i,h]+1
                }
            }
        }
    }
  
  order.string2 <- "sorted.order2<-order("
  for (h in 1:numSortCols2)
    {
      order.string2 <- paste(order.string2,"ordering2$V",h,",",sep="")
    }
  substr(order.string2,nchar(order.string2),nchar(order.string2)) <- ")"
  eval(parse(text = order.string2))
  
  blank2<-rep("",length(pos2))
  space <- rep(-1,numSites)

  padding <- as.integer((numHaps+numHaps2) * 0.05)
  if(padding < 2)
    {
      padding <- 2
    } else if(padding %% 2 == 1)
    {
      padding <- padding-1
    }
  
  combined.data <- rbind(data[sorted.order,],space);

  for (i in 1:padding)
    {
      combined.data <- rbind(combined.data,space);
    }
  
  combined.data <- rbind(combined.data,data2[sorted.order2,])
  combined.numHaps <- dim(combined.data)[1]
  combined.hap <- seq(1,combined.numHaps)
  combined.numCols <- max(numCols,numCols2)
  
  col.der <- c(rgb(227/255, 26/255, 28/255),rgb(51/255, 160/255, 44/255),rgb(31/255, 120/255, 180/255),rgb(255/255, 127/255, 0/255),rgb(106/255, 61/255, 154/255),rgb(251/255, 154/255, 153/255),rgb(178/255, 223/255, 138/255),rgb(166/255, 206/255, 227/255),rgb(253/255, 191/255, 111/255),rgb(202/255, 178/255, 214/255))
  colors.der <- rep("",combined.numCols)
  
  col.anc <- c(rgb(31/255, 120/255, 180/255),rgb(255/255, 127/255, 0/255),rgb(106/255, 61/255, 154/255),rgb(227/255, 26/255, 28/255),rgb(51/255, 160/255, 44/255),rgb(166/255, 206/255, 227/255),rgb(253/255, 191/255, 111/255),rgb(202/255, 178/255, 214/255),rgb(251/255, 154/255, 153/255),rgb(178/255, 223/255, 138/255))
  colors.anc <- rep("",combined.numCols)
  
  index <- 1;
  for(i in 0:(combined.numCols-1))
    {
      colors.der[i+1] <- col.der[(i %% 10) + 1]
      colors.anc[i+1] <- col.anc[(i %% 10) + 1]
    }

  pos <- pos/1000000
  
  print("Plotting haplotype colors...")
  setEPS()
  postscript(file=paste(filename,".hapcolor.eps",sep=""))
  
  image(pos,combined.hap,t(combined.data),zlim=c(-0.1,-0.01),yaxt="n",xaxt="n",ylab="",xlab="Distance from locus (Mb)",ylim=c(1-(padding/2),combined.numHaps+(padding/2)))
  axis(3,at=pos,lab=blank,tck=(1/(combined.numHaps) * padding/4))
  axis(3)
  axis(1,at=pos,lab=blank,tck=(1/(combined.numHaps) * padding/4))
  axis(1)
  abline(h=numHaps+(padding/2)+1)
  
  image(pos,combined.hap[1:numHaps],t(combined.data[1:numHaps,]),zlim=c(0,max(combined.data)),col=colors.der,add=TRUE)#,yaxt="n",xaxt="n",ylab="",ylim=c(1-(padding/10),combined.numHaps+(padding/10)),xlab="")
  mtext("Derived",2,at=(numHaps/2 + padding/2))  
  image(pos,combined.hap[(numHaps+1):(combined.numHaps)],t(combined.data[(numHaps+1):(combined.numHaps),]),zlim=c(0,max(combined.data)),col=colors.anc,add=TRUE)#,yaxt="n",xaxt="n",ylab="",ylim=c(-1,combined.numHaps+2),xlab="")
  mtext("Ancestral",2,at= (numHaps + (padding+1) + numHaps2/2 ) )

  dev.off()

  return
}

plotHaplotypes("/Users/amatur/code/selscan/src/out/outfile.ihs.out")