#Organize Data
leafNames <- c("B1A","B1B","B1C","B1D","B1E","B1F",
               "B2A","B2B","B2C","B2D","B2E","B2F",
               "B3A","B3B","B3C","B3D","B3E","B3F",
               "O1A","O1B","O1C","O1D","O1E","O1F",
               "O2A","O2B","O2C","O2D","O2E","O2F",
               "O3A","O3B","O3C","O3D","O3E","O3F",
               "O4A","O4B","O4C","O4D","O4E","O4F")
trees <- c("B1","B2","B3","O1","O2","O3","O4")
letterSequence <- c("a","b","c","d","e","f")
dataDirectory <- "Data/CCI_Measurements/"
heightData <- read.csv('Data/leafHeights.csv',stringsAsFactors = FALSE)
heightData$Height <- heightData$Height*2.54 #Convert from inches to cm
leafNames <- heightData$LeafName[order(heightData$Height,decreasing = FALSE)]
output <- matrix(nrow=0,ncol=4)
for(t in 1:length(trees)){
  tree <- trees[t]
  #print(tree)
  files <- dir(path=dataDirectory,pattern=tree)
  for(f in 1:length(files)){
    dat <- read.csv(paste0(dataDirectory,files[f]),stringsAsFactors = FALSE,header=TRUE)
    dte <- strsplit(as.character(dat$Time.Date[2])," ")[[1]][2]
    # if(is.na(dte)){
    #   print(tree)
    #   print(f)
    # }
    sds <- numeric()
    means <- numeric()
    vls <- numeric()
    lfNums <- character()
    for(i in 1:nrow(dat)){
      #for(i in 1:5){
      if(as.character(dat$Units[i])==" CCI"){
        sds <- c(sds,sd(vls))
        means <- c(means,mean(vls,na.rm=TRUE))
        vls <- numeric()
      }else if((is.na(dat$Sample[i]) | dat$Sample[i]=="")){
        vls <- c(vls,as.numeric(dat$Reading[i]))
      }else if(substr(dat$Sample[i],1,1)%in%c("A","B","C","D","E","F","M")){
        lfNums <- c(lfNums,substr(dat$Sample[i],1,1))
      }
    }
      output <- rbind(output,
                      cbind(paste0(tree,lfNums),rep(dte,length(means)),means,sds))
  }
}
output <- as.data.frame(output)
colnames(output) <- c("Leaf","Date","CCI_mean","CCI_sd")
output$Date <- as.Date(output$Date,format="%m/%d/%Y")
colfunc <- colorRampPalette(c("black", "white"))
cols <- c(colfunc(25)[1:18],colfunc(30)[1:24])
cols2 <- c(rep(cols[1],3),rep(cols[4],3),rep(cols[7],3),rep(cols[10],3),
           rep(cols[13],3),rep(cols[16],3),rep(cols[19],3),rep(cols[22],3),
           rep(cols[25],3),rep(cols[28],3),rep(cols[31],3),rep(cols[33],3),
           rep(cols[36],3),rep(cols[39],3))
#test <- output[order(output$Date),]
jpeg("leafCCI_data.jpeg",width=7,height=6.4,units = "in",res=1000)
par(mfrow=c(2,1))
par(mai=c(0.8,1,0.5,0.1))
subDat <- subset(output,Leaf==leafNames[19])
leafHeight <- heightData$Height[heightData$LeafName==leafNames[19]]
plot(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),type = "l",
     ylim=c(0,25),xlim=range(output$Date),main="Oak Leaves (light is higher)",ylab="CCI",xlab="",bty="n",col=cols2[19])
for(i in 20:length(leafNames)){
  subDat <- subset(output,Leaf==leafNames[i])
  leafHeight <- heightData$Height[heightData$LeafName==leafNames[i]]
  lines(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),col=cols2[i])
}

subDat <- subset(output,Leaf==leafNames[1])
leafHeight <- heightData$Height[heightData$LeafName==leafNames[1]]
plot(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),type = "l",
     ylim=c(0,25),xlim=range(output$Date),
     xlab="Date", ylab="CCI",bty="n",main="Beech Leaves (light is higher)",col=cols2[1])
for(i in 2:18){
  subDat <- subset(output,Leaf==leafNames[i])
  leafHeight <- heightData$Height[heightData$LeafName==leafNames[19]]
  lines(subDat$Date,as.numeric(as.character(subDat$CCI_mean)),col=cols2[i])
}

dev.off()







               