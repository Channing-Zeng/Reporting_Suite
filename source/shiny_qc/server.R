
library(shiny)
#library(XML)


#theurl <- "../NextSeqRunQC/results/Reports/html/H23FMAFXX/all/all/all/laneBarcode.html"
#tables <- readHTMLTable(theurl)
#sbatch NextSeqQC.sh $1
#  src=/projects/NGS/projects/DS/fromNextSeq500/150707_NS500580_0024_AH23FMAFXX
#sheet=/projects/NGS/projects/DS/fromNextSeq500/150707_NS500580_0024_AH23FMAFXX/SampleSheet.csv
#res="."
library("XML")
require("stringr")
require("raster")
require("data.table")
require(gplots)

shinyServer(function(input, output) {

  statsHtml = reactive({
    rn=input$txtInRunName
    # rn="150709_NS500580_0025_AH23FWAFXX"  
    #rn="150310_NS500580_0004_AH2CJ7AFXX"
    rnX=str_split_fixed(rn,pattern="_",10)[,4]
    rnX=str_replace(rnX,"^A","")

    fl=paste(rn,"/Reports/html/",rnX,"/all/all/all/laneBarcode.html",sep="")
    #fl=      "/projects/NGS/projects/DS/results/150522_NS500580_0012_AH352YBGXX/Reports/html/H352YBGXX/all/all/all/laneBarcode.html"
    readHTMLTable(fl)    
    #statsHtml=readHTMLTable(fl)    
  })

  output$readPlot <- renderPlot({    
    #x=statsHtml[[3]] 
    x=statsHtml()[[3]] 
    x=x[2:nrow(x),]
    rdsTot = as.integer(str_replace_all( as.matrix(x[,9]), ",", ""))/1E6
    xx=aggregate(rdsTot,by=list(paste(x[,2]),paste(x[,3])), FUN=sum ) 
    totRun = sum(rdsTot)
    idx=xx$Group.2 != "unknown"    
    if (sum(idx) > 0){
      xxx=xx[idx,]
      if (length(unique(xxx$Group.1)) > 1){
        lbs=paste(xxx$Group.1,xxx$Group.2)  
        xxx=xxx[order(lbs, decreasing=T),]
        lbs=paste(xxx$Group.1,xxx$Group.2)  
      }else{
        lbs=xxx$Group.2
        xxx=xxx[order(lbs, decreasing=T),]
        lbs=xxx$Group.2
      }    
      par(mar=c(4,10,2,2))
      xv = xxx$x
      xv = xv - min(xv)
      xv = xv/max(xv)
      xv = as.integer(99*xv)+1
      xcols=colorRampPalette(c("red","yellow","green"))(100)[xv]
      bp=barplot(xxx$x, names.arg=lbs, horiz=T,  las=2, main="Reads per sample", col=xcols )
      #text(x, bp, signif(x,2), pos=4)
      p=paste0(formatC(100 * xxx$x/totRun,  digits = 2), "%")     
      text(x=xxx$x,bp, p   , pos=4)
    }else{
      barplot(1, horiz=T, main="No samples found")
    }
  })
  
  output$BCPlot = renderPlot({
    #dtx=statsHtml[[3]] 
    dtx=statsHtml()[[3]] 
    dtx=dtx[2:nrow(dtx),]
    rdsTot = as.integer(str_replace_all( as.matrix(dtx[,5]), ",", ""))
    xx=aggregate(rdsTot,by=list(paste(dtx[,2]),paste(dtx[,3])), FUN=sum ) 
    totRun = sum(rdsTot)
    
    #dtx=statsHtml[[4]]         
    dtx=statsHtml()[[4]]         
    Ns = t(dtx[1,c(2,5,8,11)]    )
    Ns = as.integer( str_replace_all(Ns,",",""))
    BC =  rbind(
      as.matrix(dtx[3:11,1:2]),
      as.matrix(dtx[3:11,3:4]),
      as.matrix(dtx[3:11,5:6]),
      as.matrix(dtx[3:11,7:8]))
    cts = as.integer( str_replace_all(BC[,1],",",""))    
    regbcs = paste("Registered BCs(",as.integer(100*(totRun - sum(Ns) - sum(cts)) / totRun  ),"%)",sep="")
    cts=rbind(aggregate(cts,by=list(BC[,2]), FUN=sum ), 
              c("NNNNNNNN-NNNNNNNN", sum(Ns)) ,
              c(regbcs, totRun - sum(Ns) - sum(cts)  )
    )
    
    
    x=as.integer(cts[,2])/1E6
    
    xv = x
    xv = xv - min(xv)
    xv = xv/max(xv)
    xv = as.integer(99*xv)+1
    xcols=colorRampPalette(c("red","blue"))(100)[xv]
    
    lbs=cts[,1]
    par(mar=c(4,14,2,2))
    bp=barplot(x, names.arg=lbs, horiz=T,  las=2, main="Reads By Barcode", col=xcols )
    #text(x, bp, signif(x,2), pos=4)
    p=as.integer(cts[,2])/totRun
    p=paste0(formatC(100 * p,  digits = 2), "%")     
    text(x=x,bp, p   , pos=4)
    
  })
  
  
  
  output$FilteredPlot = renderPlot({
    
    #x=statsHtml[[3]] 
    x=statsHtml()[[3]] 
    x=x[2:nrow(x),]
    rdsTot = as.integer(str_replace_all( as.matrix(x[,5]), ",", ""))/1E6
    rdsFilt = as.integer(str_replace_all( as.matrix(x[,9]), ",", ""))/1E6
    xxx=aggregate(cbind(rdsTot,rdsFilt),by=list(paste(x[,2]),paste(x[,3])), FUN=sum )   
    
    
    if (length(unique(xxx$Group.1)) > 1){
      lbs=paste(xxx$Group.1,xxx$Group.2)  
      xxx=xxx[order(lbs, decreasing=T),]
      lbs=paste(xxx$Group.1,xxx$Group.2)  
    }else{
      lbs=xxx$Group.2
      xxx=xxx[order(lbs, decreasing=T),]
      lbs=xxx$Group.2
    }    
    par(mar=c(4,10,2,2))
    xxx$x = xxx$rdsFilt/xxx$rdsTot
    xv = xxx$x
    xv = xv - min(xv)
    xv = xv/max(xv)
    xv = as.integer(99*xv)+1
    xcols=colorRampPalette(c("red","yellow","green"))(100)[xv]
    bp=barplot(1-xxx$x, names.arg=lbs, horiz=T,  las=2, main="Filtered Reads", col=xcols )
    #text(x, bp, signif(x,2), pos=4)
    p=paste0(formatC(100 * (1- xxx$x),  digits = 2), "%")     
    text(x=xxx$x,bp, p   , pos=4)    
  })
  
  output$byGenomePlot = renderPlot({
    
  })
  
  metaMapDtx = reactive({    
    rn=input$txtInRunName
    # rn="150709_NS500580_0025_AH23FWAFXX"  
    #rn="150310_NS500580_0004_AH2CJ7AFXX"
    rnX=str_split_fixed(rn,pattern="_",10)[,4]
    rnX=str_replace(rnX,"^A","")    
    fl=paste("/projects/NGS/projects/DS/PreProcessed/",rn,"/MetaMapping/Mapping1M.tab",sep="")
    dtx=read.table(fl, sep="\t", stringsAsFactors=F)
    dtx$V1= unlist(lapply(paste(dtx$V1),  FUN=basename))    
    dtx$V1 = str_replace_all(dtx$V1  ,"_R1_001.fastq.gz","")    
    dtx=dtx[dtx$V5 == "HIGH",]    
    colnames(dtx) = c("Fastq","genome","samp","chromo","MapQ","Cts")
    dtx=dtx[order(dtx$samp),]
    dtx$Fastq = NULL
    dtx    
  })
  
  
  output$DataTableMapping = renderDataTable({
    metaMapDtx()
  })
  
  
  output$MetaGenomeMappingPlot = renderPlot({
    dtx = metaMapDtx()
    srt=list(hg38=1,mm10=2,rn5=3,canFam3=4,dm6=5,ce10=6,sacCer3=7,K12=8,PhiX174=9)
    
    dtxx = data.table(dtx[dtx$MapQ == "HIGH",])
    dtxx = dtxx[,by=list(genome,samp), sum(Cts) ]
    dtxx$X = as.integer(srt[dtxx$genome])
    dtxx$Y = as.integer(factor(dtxx$samp))
    mtx=matrix(nrow=max(dtxx$Y),ncol=9 )
    mtx[cbind(dtxx$Y,dtxx$X)] = dtxx$V1
    par(mar=c(4,4,4,4))
    
    y=cbind( as.integer(factor(dtxx$samp)), dtxx$samp)        
    y=data.frame(y[!duplicated(y),], stringsAsFactors=F)  
    
    
    rownames(mtx) = y$X2
    colnames(mtx) = labels(srt)
    
    
    
    par(mar=c(3,3,3,3))
    #dev.off()
    heatmap.2(mtx,Rowv=NULL,Colv=NULL, scale="none", trace="none",dendrogram='none',
              col=colorRampPalette(colors=c("white","black","red"))(100) ,               
              lmat =rbind(c(4,3),c(1,2),c(0,0) ),
              lwid = c(4,1),
              lhei = c(1,4,1), 
              cexRow=1.5, cexCol=1.5
    )  
    
    
  })
  
  output$MetaMappingPlot = renderPlot({
    dtx = metaMapDtx()    
    srt=list(hg38=1000,mm10=2000,rn5=3000,canFam3=4000,dm6=5000,ce10=6000,sacCer3=7000,K12=8000,PhiX174=9000)
    
    dtx$Y = as.integer(factor(dtx$samp))
    dtChromos = read.table("/projects/NGS/projects/DS/Scripts/MetaChromos.tab", stringsAsFactors=F)
    dtChromos=dtChromos[order(paste( dtChromos$V1, dtChromos$V2)),]
    colnames(dtChromos) = c("genome","chromo")
    chrI=str_replace_all(string= dtChromos$chromo, pattern="chr","")
    chrI=as.integer(chrI)
    chrI[is.na(chrI)] = 100     
    chrI = order(unlist(srt[dtChromos$genome]) +chrI)
    dtChromos = dtChromos[ chrI ,]
    dtChromos$X = 1:nrow(dtChromos)
    
    dtx = merge(dtx,dtChromos,all.y=T)    
    dtx$samp[is.na(dtx$samp)] = ""
    dtx$Cts[is.na(dtx$Cts)] = 0
    dtx$Y[is.na(dtx$Y)] = max(dtx$Y[!is.na(dtx$Y)])+1
    
    ylabs=cbind(dtx$Y, dtx$samp)
    ylabs=ylabs[!duplicated(ylabs),]
    ylabs=ylabs[order(as.integer(ylabs[,1])),]
    ylabs=ylabs[,2]
    
    xlabs=cbind(dtx$X, dtx$chromo, dtx$genome)
    xlabs=xlabs[!duplicated(xlabs),]
    xlabs=xlabs[order(as.integer(xlabs[,1])),]
    xlabs[duplicated(xlabs[,3]) ,3] = ""
    xlabs=paste(xlabs[,3],xlabs[,2])  
    
    xmax=max(dtx$X)+1
    ymax=max(dtx$Y)+1
    mtx=matrix(nrow=ymax,ncol=xmax )
    mtx[cbind(dtx$Y,dtx$X)] = dtx$Cts
    par(mar=c(3,3,3,3))
    #dev.off()
    heatmap.2(mtx,Rowv=NULL,Colv=NULL, scale="none", trace="none",dendrogram='none',
              col=colorRampPalette(colors=c("white","black","red"))(100) , 
              labCol=xlabs, labRow=ylabs,                            
              lmat =rbind(c(4,3),c(1,2),c(0,0) ),
              lwid = c(4,1),
              lhei = c(1,4,1)
              )  
  
    
    })
  
  
})





plotFlowCell = function(dtx){
  
  fl="/projects/NGS/projects/DS/Scripts/FlowCellStatsTest.tab.gz"
  dtx=read.table(fl,sep="\t", header=T)
  
  
  ymax=max(dtx$Y)
  xmax=max(dtx$X)
  tlsYOff = (as.integer(factor(paste(dtx$Lane, dtx$Tile) )) - 1) * ymax  
  yy=length(unique( paste(dtx$Lane, dtx$Tile)  ))*ymax + ymax
  xx=(ncol(dtx)-5 )*xmax  + xmax
  mtx=matrix(nrow=yy,ncol=xx )
  for (i in 0:50){
    c = paste("V",i+1,sep="")
    Xoff = xmax*i
    y=dtx$Y+tlsYOff
    x=dtx$X+Xoff    
    mtx[cbind(y,x)] = dtx[,get(c)]/dtx$ct
    print(i)    
  }
  par(mar=c(0,0,0,0))
  idx=mtx< 60
  idx[is.na(idx)] = F
  mtx[idx] = 60
  plot(
    raster(mtx ), col=rev(topo.colors(100))
  )
}


