viewR<-function(data=NULL,ret=FALSE){
  if(is.null(data)){
    file<-selectR()
    if(length(file)==0){stop("You must select a file")}
    data<-readNii(file)
    }
  if(is.character(data)){
    file<-data
    if(length(file)==0){stop("You must select a file")}
    data<-readNii(file)
  }
    
  if(!is.array(data)){
    stop("This is not an array with 3 or more dimensions")
  }
  orig<-data
  ##############################
  ##### DESIRED HEIGHT #########
  ##############################
  screenHeight<-as.numeric(strsplit(tclvalue(tkwm.maxsize("."))," ")[[1]])[2]
  desHeight<-round(screenHeight/5*3)
  ##############################
  ##### ASPECT RATIO ###########
  ##############################
  d<-dim(data)
  d1<-d[1]
  d2<-d[2]
  d3<-d[3]
  xx<-round(d1/2)
  yy<-round(d2/2)
  zz<-round(d3/2)
  if(length(d)>3){
    d4<-d[4]
  }else{
    dim(data)<-c(dim(data),1)
    d4<-1
  }
  
  xCoords<-1:d1
  yCoords<-1:d2
  zCoords<-1:d3
  tCoords<-1:d4
  at<-attr(data,"pixdim")
  if(is.null(at)){pixdim<-c(1,1,1)}else{pixdim<-at}
  if(length(pixdim)==8){pixdim<-pixdim[2:4]}
  data<-zeroNa(input = data)
  
  w1<-d1*pixdim[1]
  w2<-d2*pixdim[2]
  w3<-d3*pixdim[3]
  relHeights<-c(w3/(w3+w2),w2/(w3+w2))
  relWidths<-c(w1/(w1+w2),w2/(w1+w2))
  relWidths<-relWidths*relHeights[2]/relWidths[2]
  m<-max(c(relWidths,relHeights))
  
  relHeights<-relHeights/m
  relWidths<-relWidths/m
  #########################
  ### IMAGE VARIABLES #####
  #########################
  r<-range(data)
  
  time<-tclVar("1")
  xl<-tclVar(as.character(xx))
  yl<-tclVar(as.character(yy))
  zl<-tclVar(as.character(zz))
  
  X<-tclVar(tclvalue(xl))
  Y<-tclVar(tclvalue(yl))
  Z<-tclVar(tclvalue(zl))
  
  parPlotSize1<-c()
  usrCoords1<-c()
  parPlotSize2<-c()
  usrCoords2<-c()
  parPlotSize3<-c()
  usrCoords3<-c()
  parPlotSize4<-c()
  usrCoords4<-c()
  
  low<-tclVar(as.character(r[1]))
  high<-tclVar(as.character(r[2]))
  intens<-tclVar(as.character(data[xx,yy,zz,1]))
  crossHairsOn<-TRUE
  #########################
  ####### FRAMES ##########
  #########################
  tclServiceMode(FALSE)
  base<-tktoplevel()
  
  tktitle(base)<-"Display"
  master<-tkframe(parent=base)
  f1<-tkframe(parent = master,borderwidth=0,relief="flat")
  f2<-tkframe(parent = master,borderwidth=0,relief="flat")
  f3<-tkframe(parent = master,borderwidth=0,relief="flat")
  #########################
  ##### FUNCTIONS #########
  #########################
  plotf1<-function(){
    
    par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
    xl<-as.numeric(tclvalue(xl))
    yl<-as.numeric(tclvalue(yl))
    zl<-as.numeric(tclvalue(zl))
    x<-round(xl)
    y<-round(yl)
    z<-round(zl)
    tclvalue(X)<<-x
    tclvalue(Y)<<-y
    tclvalue(Z)<<-z
    
    t<-as.numeric(tclvalue(time))
    if(is.na(t)){t<-1}
    tclvalue(time)<<-t
    
    lim<-c(as.numeric(tclvalue(low)),as.numeric(tclvalue(high)))
    if(any(is.na(lim))){lim<-r;tclvalue(low)<<-r[1];tclvalue(high)<<-r[2]}
    tclvalue(intens)<<-data[x,y,z,t]
    
    im<-data[,y,,t]
    im[im<lim[1]]<-lim[1]
    im[im>lim[2]]<-lim[2]
    
    image(1:d1,1:d3,z=im,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = zl,v = xl,col="green")
    
    parPlotSize1 <<- par("plt")
    usrCoords1   <<- par("usr")
    
  }
  OnLeftClick1 <- function(x,y){
    xClick <- x
    yClick <- y
    
    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img1)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img1)))
    
    xMin <- parPlotSize1[1] * width
    xMax <- parPlotSize1[2] * width
    yMin <- parPlotSize1[3] * height
    yMax <- parPlotSize1[4] * height
    
    rangeX <- usrCoords1[2] - usrCoords1[1]
    rangeY <- usrCoords1[4] - usrCoords1[3]
    
    imgXcoords <- (xCoords-usrCoords1[1])*(xMax-xMin)/rangeX + xMin
    imgYcoords <- (zCoords-usrCoords1[3])*(yMax-yMin)/rangeY + yMin
    
    xClick <- as.numeric(xClick)+1
    yClick <- as.numeric(yClick)+1
    yClick <- height - yClick
    
    xl<- usrCoords1[1]+(xClick-xMin)*rangeX/(xMax-xMin)
    zl<- usrCoords1[3]+(yClick-yMin)*rangeY/(yMax-yMin)
    
    if(xl>d1){xl<-d1}
    if(zl>d3){zl<-d3}
    if(xl<1){xl<-1}
    if(zl<1){zl<-1}
    tclvalue(xl)<<-xl
    tclvalue(zl)<<-zl
    tclvalue(X)<<-round(xl)
    tclvalue(Z)<<-round(zl)
    
    tkrreplot(img1)
    
    tkrreplot(img2)
    
    tkrreplot(img3)
    
    if(d4>1){tkrreplot(img4)}
    
  }
  
  plotf2<-function(){
    
    par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
    xl<-as.numeric(tclvalue(xl))
    yl<-as.numeric(tclvalue(yl))
    zl<-as.numeric(tclvalue(zl))
    x<-round(xl)
    y<-round(yl)
    z<-round(zl)
    tclvalue(X)<<-x
    tclvalue(Y)<<-y
    tclvalue(Z)<<-z
    t<-as.numeric(tclvalue(time))
    if(is.na(t)){t<-1}
    tclvalue(time)<<-t
    
    
    lim<-c(as.numeric(tclvalue(low)),as.numeric(tclvalue(high)))
    if(any(is.na(lim))){lim<-r;tclvalue(low)<<-r[1];tclvalue(high)<<-r[2]}
    im<-data[x,,,t]
    im[im<lim[1]]<-lim[1]
    im[im>lim[2]]<-lim[2]
    
    
    image(1:d2,1:d3,z=im,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = zl,v = yl,col="green")
    
    
    parPlotSize2 <<- par("plt")
    usrCoords2   <<- par("usr")
    
  }
  OnLeftClick2 <- function(x,y){
    xClick <- x
    yClick <- y
    
    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img2)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img2)))
    
    xMin <- parPlotSize2[1] * width
    xMax <- parPlotSize2[2] * width
    yMin <- parPlotSize2[3] * height
    yMax <- parPlotSize2[4] * height
    
    rangeX <- usrCoords2[2] - usrCoords2[1]
    rangeY <- usrCoords2[4] - usrCoords2[3]
    
    imgXcoords <- (yCoords-usrCoords2[1])*(xMax-xMin)/rangeX + xMin
    imgYcoords <- (zCoords-usrCoords2[3])*(yMax-yMin)/rangeY + yMin
    
    xClick <- as.numeric(xClick)+1
    yClick <- as.numeric(yClick)+1
    yClick <- height - yClick
    
    yl <- usrCoords2[1]+(xClick-xMin)*rangeX/(xMax-xMin)
    zl <- usrCoords2[3]+(yClick-yMin)*rangeY/(yMax-yMin)
    
    if(yl>d2){yl<-d2}
    if(zl>d3){zl<-d3}
    if(yl<1){yl<-1}
    if(zl<1){zl<-1}
    tclvalue(yl)<<-yl
    tclvalue(zl)<<-zl
    
    tkrreplot(img1)
    
    tkrreplot(img2)
    
    tkrreplot(img3)
    
    if(d4>1){tkrreplot(img4)}
  }
  
  plotf3<-function(){
    
    par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
    xl<-as.numeric(tclvalue(xl))
    yl<-as.numeric(tclvalue(yl))
    zl<-as.numeric(tclvalue(zl))
    x<-round(xl)
    y<-round(yl)
    z<-round(zl)
    tclvalue(X)<<-x
    tclvalue(Y)<<-y
    tclvalue(Z)<<-z
    t<-as.numeric(tclvalue(time))
    if(is.na(t)){t<-1}
    tclvalue(time)<<-t
    
    lim<-c(as.numeric(tclvalue(low)),as.numeric(tclvalue(high)))
    if(any(is.na(lim))){lim<-r;tclvalue(low)<<-r[1];tclvalue(high)<<-r[2]}
    
    im<-data[,,z,t]
    im[im<lim[1]]<-lim[1]
    im[im>lim[2]]<-lim[2]
    image(1:d1,1:d2,z=im,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = yl,v = xl,col="green")
    
    parPlotSize3 <<- par("plt")
    usrCoords3   <<- par("usr")
    
  }
  OnLeftClick3 <- function(x,y){
    xClick <- x
    yClick <- y
    
    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img3)))
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img3)))
    
    xMin <- parPlotSize3[1] * width
    xMax <- parPlotSize3[2] * width
    yMin <- parPlotSize3[3] * height
    yMax <- parPlotSize3[4] * height
    
    rangeX <- usrCoords3[2] - usrCoords3[1]
    rangeY <- usrCoords3[4] - usrCoords3[3]
    
    imgXcoords <- (xCoords-usrCoords3[1])*(xMax-xMin)/rangeX + xMin
    imgYcoords <- (yCoords-usrCoords3[3])*(yMax-yMin)/rangeY + yMin
    
    xClick <- as.numeric(xClick)+1
    yClick <- as.numeric(yClick)+1
    yClick <- height - yClick
    
    xl <- usrCoords3[1]+(xClick-xMin)*rangeX/(xMax-xMin)
    yl <- usrCoords3[3]+(yClick-yMin)*rangeY/(yMax-yMin)
    
    if(xl>d1){xl<-d1}
    if(yl>d2){yl<-d2}
    if(xl<1){xl<-1}
    if(yl<1){yl<-1}
    tclvalue(xl)<<-xl
    tclvalue(yl)<<-yl
    
    tkrreplot(img1)
    
    tkrreplot(img2)
    
    tkrreplot(img3)
    
    if(d4>1){tkrreplot(img4)}
    
  }
  
  if(d4>1){
    plotf4<-function(){
      par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
      xl<-as.numeric(tclvalue(xl))
      yl<-as.numeric(tclvalue(yl))
      zl<-as.numeric(tclvalue(zl))
      tclvalue(time)<<-tclvalue(time)
      x<-round(xl)
      y<-round(yl)
      z<-round(zl)
      t<-as.numeric(tclvalue(time))
      tseries<-data[x,y,z,]
      plot(1:d4,tseries,type="l",axes=FALSE,col="grey")
      points(x = t,y = tseries[t],col="green")
      
      parPlotSize4 <<- par("plt")
      usrCoords4   <<- par("usr")
    }
    
    OnLeftClick4 <- function(x,y){
      xClick <- x
      yClick <- y
      
      width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img4)))
      height <- as.numeric(tclvalue(tkwinfo("reqheight",img4)))
      
      xMin <- parPlotSize4[1] * width
      xMax <- parPlotSize4[2] * width
      yMin <- parPlotSize4[3] * height
      yMax <- parPlotSize4[4] * height
      
      rangeX <- usrCoords4[2] - usrCoords4[1]
      rangeY <- usrCoords4[4] - usrCoords4[3]
      
      imgXcoords <- (xCoords-usrCoords4[1])*(xMax-xMin)/rangeX + xMin
      imgYcoords <- (yCoords-usrCoords4[3])*(yMax-yMin)/rangeY + yMin
      
      xClick <- as.numeric(xClick)+1
      yClick <- as.numeric(yClick)+1
      yClick <- height - yClick
      
      cox<- usrCoords4[1]+(xClick-xMin)*rangeX/(xMax-xMin)
      coy<- usrCoords4[3]+(yClick-yMin)*rangeY/(yMax-yMin)
      xl<-as.numeric(tclvalue(xl))
      yl<-as.numeric(tclvalue(yl))
      zl<-as.numeric(tclvalue(zl))
      tclvalue(time)<<-tclvalue(time)
      x<-round(xl)
      y<-round(yl)
      z<-round(zl)
      
      tseries<-data[x,y,z,]
      #tclvalue(time)<<-which.min(sqrt(colSums((t(cbind(tseries,1:d4))-c(coy,cox))^2)))
      tclvalue(time)<<-which.min(abs((1:d4)-cox))
      
      
      tkrreplot(img1)
      
      tkrreplot(img2)
      
      tkrreplot(img3)    
      
      if(d4>1){tkrreplot(img4)}
      
    }
  }
  
  changeTemp<-function(){
    
    if(is.na(strtoi(tclvalue(time)))){
      tclvalue(time)<<-as.character(1)
    }
    t<-strtoi(tclvalue(time))
    if(t>d4){tclvalue(time)<<-as.character(d4)}
    if(t<1){tclvalue(time)<<-as.character(1)}
    tkrreplot(img1)
    
    tkrreplot(img2)
    
    tkrreplot(img3)
    
    if(d4>1){tkrreplot(img4)}
  }
  tkspinbox <- function(parent, ...) {tkwidget(parent, "tk::spinbox", ...)}
  onSpin<-function(){
    if(is.na(strtoi(tclvalue(X)))||strtoi(tclvalue(X))==0){
      tclvalue(xl)<<-as.character(xx)
      tclvalue(X)<<-as.character(xx)}else{
        tclvalue(xl)<<-tclvalue(X)
      }
    
    if(is.na(strtoi(tclvalue(Y)))||strtoi(tclvalue(Y))==0){
      tclvalue(yl)<<-as.character(yy)
      tclvalue(Y)<<-as.character(yy)}else{
        tclvalue(yl)<<-tclvalue(Y)
      }
    
    if(is.na(strtoi(tclvalue(Z)))||strtoi(tclvalue(Z))==0){
      tclvalue(zl)<<-as.character(zz)
      tclvalue(Z)<<-as.character(zz)}else{
        tclvalue(zl)<<-tclvalue(Z)
      }
    if(is.na(strtoi(tclvalue(time)))||strtoi(tclvalue(time))==0){
      tclvalue(time)<<-as.character(1)}else{
        tclvalue(time)<<-tclvalue(time)
      }
    x<-round(as.numeric(tclvalue(xl)))
    y<-round(as.numeric(tclvalue(yl)))
    z<-round(as.numeric(tclvalue(zl)))
    
    if(x<1){tclvalue(xl)<<-1}
    if(y<1){tclvalue(yl)<<-1}
    if(z<1){tclvalue(zl)<<-1}
    
    if(x>d1){tclvalue(xl)<<-d1}
    if(y>d2){tclvalue(yl)<<-d2}
    if(z>d3){tclvalue(zl)<<-d3}
    
    
    
    tkrreplot(img1)
    
    tkrreplot(img2)
    
    tkrreplot(img3)
    
    if(d4>1){tkrreplot(img4)}
    
  }
  maxminChange<-function(){
    if(is.na(strtoi(tclvalue(high)))||tclvalue(high)==""){
      tclvalue(high)<<-r[2]
    }
    if(is.na(strtoi(tclvalue(low)))||tclvalue(low)==""){
      tclvalue(low)<<-r[1]
    }
    max<-as.numeric(tclvalue(high))
    min<-as.numeric(tclvalue(low))
    if(max<min){tclvalue(high)<<-r[2];tclvalue(low)<<-r[1]}
    tkrreplot(img1)  
    tkrreplot(img2)
    tkrreplot(img3)
  }
  crossHairs<-function(){
    crossHairsOn<<-!crossHairsOn
    tkrreplot(img1)
    tkrreplot(img2)
    tkrreplot(img3)
    if(d4>1){tkrreplot(img4)}
  }
  
  
  
  
  ##########################
  ###### THE PLOTS #########
  ##########################
  testImg<-tkrplot(parent = f1,fun =plotf1,hscale = 1,vscale = 1)
  testHeight<-as.numeric(tkwinfo("reqheight",testImg))
  scaleFactor<-desHeight/(testHeight*2)
  
  img1<-tkrplot(parent = f1,fun = plotf1,hscale=relWidths[1]*scaleFactor,vscale=relHeights[1]*scaleFactor)
  img2<-tkrplot(parent = f2,fun = plotf2,hscale=relWidths[2]*scaleFactor,vscale=relHeights[1]*scaleFactor)
  img3<-tkrplot(parent = f3,fun = plotf3,hscale=relWidths[1]*scaleFactor,vscale=relHeights[2]*scaleFactor)
  
  f4<-tkframe(parent = master,width=as.numeric(tkwinfo("reqwidth",img2)),height=as.numeric(tkwinfo("reqheight",img3)))
  f5<-tkframe(parent=master,
              width=as.numeric(tkwinfo("reqwidth",img1))+as.numeric(tkwinfo("reqwidth",img2)),
              height=100,
              borderwidth=2,
              relief="groove")
  if(d4>1){
    img4<-tkrplot(parent = f4,fun = plotf4,hscale=relWidths[2]*scaleFactor,vscale=relWidths[2]*scaleFactor)
  }
  ##########################
  ###### GEOMETRY ##########
  ##########################
  tkgrid(master,columnspan=2)
  tkgrid(f1,f2)
  tkgrid(f3,f4)
  tkgrid(img1)
  tkgrid(img2)
  tkgrid(img3)
  if(d4>1){
    tkgrid(img4)
  }
  tkgrid(f5,columnspan=2,rowspan=1)
  ##########################
  ##### Widgets ############
  ##########################
  spin<-tkspinbox(f5,textvariable=time,from=1,to=d4,command=changeTemp,increment=1,repeatdelay=200,width=5,wrap=TRUE)
  coX<-tkspinbox(f5,textvariable=X,from=1,to=d1,command=function()onSpin(),increment=1,repeatdelay=10,width=5)
  coY<-tkspinbox(f5,textvariable=Y,from=1,to=d2,command=function()onSpin(),increment=1,repeatdelay=10,width=5)
  coZ<-tkspinbox(f5,textvariable=Z,from=1,to=d3,command=function()onSpin(),increment=1,repeatdelay=10,width=5)
  MAX<-tkentry(f5,textvariable=high,width=7)
  MIN<-tkentry(f5,textvariable=low,width=7)
  intensity<-tkentry(f5,textvariable=intens,state="readonly",readonlybackground="white",width=7)
  coXLab<-tklabel(parent = f5,text="X")
  coYLab<-tklabel(parent = f5,text="Y")
  coZLab<-tklabel(parent = f5,text="Z")
  coTLab<-tklabel(parent = f5,text="Time")
  maxLab<-tklabel(parent = f5,text="Max")
  minLab<-tklabel(parent = f5,text="Min")
  intensLab<-tklabel(parent=f5,text="Intensity")
  tmovie<-tclVar(FALSE)
  
  after_ID <-""
  someFlag<-TRUE
  
  repeat_call<-function(ms = 200 , f) {
    after_ID <<- tcl( "after" , ms,function(){ 
      if(someFlag){
        f()
        after_ID<<-repeat_call(ms,f)
      }else{
        tcl("after" , "cancel" , after_ID)
      }
    })
  }
  movieT<-function(){
    if(tclvalue(tmovie)==1){
      
      someFlag<<-TRUE
      repeat_call(1,function() {tkinvoke(spin,"buttonup");tkrreplot(img4)})
    }else{someFlag<<-FALSE}
  }
  tmovieBut<-tkcheckbutton(f5,variable=tmovie, command=movieT,text="Movie")
  
  
  tkplace(coXLab,relx=0,rely=1/5,y=-10)
  tkplace(coYLab,relx=0,rely=2/5,y=-10)
  tkplace(coZLab,relx=0,rely=3/5,y=-10)
  tkplace(coTLab,relx=0,rely=4/5,y=-10)
  
  tkplace(coX,'in'=coXLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(coY,'in'=coYLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(coZ,'in'=coZLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(spin,'in'=coTLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  
  if(d4>1){
    tkplace(tmovieBut,"in"=spin,rely=.5,anchor="w",x=as.numeric(tkwinfo("reqwidth",spin)),height=100/4)
  }
  wiBox<-as.numeric(tkwinfo("reqwidth",intensity))
  tkplace(intensity,relx=1,x=-wiBox,rely=0/4,y=10)
  tkplace(MAX,relx=1,x=-wiBox,rely=1/4,y=10)
  tkplace(MIN,relx=1,x=-wiBox,rely=2/4,y=10)
  
  tkplace(intensLab,'in'=intensity,x=-1,rely=.5,anchor="e",height=100/4)
  tkplace(maxLab,'in'=MAX,x=-1,rely=.5,anchor="e",height=100/4)
  tkplace(minLab,'in'=MIN,x=-1,rely=.5,anchor="e",height=100/4)
  
  
  ###########################
  ##### BINDINGS ############
  ###########################
  tkbind(img1, "<Button-1>",OnLeftClick1)
  tkbind(img1, "<B1-Motion>",OnLeftClick1)
  tkbind(img2, "<Button-1>",OnLeftClick2)
  tkbind(img2, "<B1-Motion>",OnLeftClick2)
  tkbind(img3, "<Button-1>",OnLeftClick3)
  tkbind(img3, "<B1-Motion>",OnLeftClick3)
  if(d4>1){
    tkbind(img4, "<Button-1>",OnLeftClick4)
    tkbind(img4, "<B1-Motion>",OnLeftClick4)
    tkconfigure(img4,cursor="crosshair")
  }
  tkbind(MAX, "<Return>",function()maxminChange())
  tkbind(MIN, "<Return>",function()maxminChange())
  tkbind(coX, "<Return>",function()onSpin())
  tkbind(coY, "<Return>",function()onSpin())
  tkbind(coZ, "<Return>",function()onSpin())
  tkbind(spin, "<Return>",function()changeTemp())
  tkbind(img1, "<Button-3>",crossHairs)
  tkbind(img2, "<Button-3>",crossHairs)
  tkbind(img3, "<Button-3>",crossHairs)
  
  tkconfigure(img1,cursor="crosshair")
  tkconfigure(img2,cursor="crosshair")
  tkconfigure(img3,cursor="crosshair")
  
  tclServiceMode(TRUE)
  tkwm.resizable(base,FALSE,FALSE)
  tkfocus(base)
  if(ret){return(orig)}
}
