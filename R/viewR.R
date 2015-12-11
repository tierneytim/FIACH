viewR<-function(data=NULL,overlay=NULL,otherData=NULL,xyz=NULL,ret=FALSE){
  ###############################
  ###### TKRPLOT CHECK ##########
  ###############################
  tkrPossible<-requireNamespace("tkrplot",quietly = TRUE)
  if(!tkrPossible){stop("tkrplot is not available or could not be loaded. Please install to use viewR")}
  ###############################
  ###### DATA CHECK #############
  ###############################
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
  #################################
  ######## OVERLAY CHECK ##########
  #################################
  if(is.character(overlay)){
    file<-overlay
    if(length(file)==0){stop("You must select a file")}
    overlay<-readNii(file)
  }
  if(!is.null(overlay)){
    if(length(dim(overlay))==3){dim(overlay)<-c(dim(overlay),1)}
  if(!is.array(overlay)){
    stop("This is not an array with 3 or more dimensions")
  }
  if(length(dim(overlay))==4){if(dim(overlay)[4]>1){stop("Overlay must have 3 and only three dimensions")}}
    if(any(dim(overlay)[1:3]!=dim(data)[1:3])){stop("Overlay must have same dimensions as data")}
  }
  olay<-!is.null(overlay)
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
  xi<-xyz[1]
  yi<-xyz[2]
  zi<-xyz[3]
  xx <- ifelse(is.null(xi)||xi>d1||xi<1||!is.numeric(xi), round(d1/2), xi)
  yy <- ifelse(is.null(yi)||yi>d2||yi<1||!is.numeric(xi), round(d2/2), yi)
  zz <- ifelse(is.null(zi)||zi>d3||zi<1||!is.numeric(xi), round(d3/2), zi)
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
  null2<-is.null(otherData)
  if(null2){
    otherData<-array(0,dim = c(dim(data)[1:3],1))
  }else{
    dOther<-dim(otherData)
    if(any(dOther[1:3]!=d[1:3])){stop("dimensions of data and otherData must match")}
    if(length(dOther)==3){dim(otherData)<-c(dOther,1)
    }
  }
  origOther<-otherData
  w1<-d1*pixdim[1]
  w2<-d2*pixdim[2]
  w3<-d3*pixdim[3]
  #########################
  ### IMAGE VARIABLES #####
  #########################
  r<-round(range(data),digits = 2)
  r2<-round(range(otherData),digits = 2)
  time<-tclVar("1")
  xl<-tclVar(as.character(xx))
  yl<-tclVar(as.character(yy))
  zl<-tclVar(as.character(zz))
  X<-tclVar(tclvalue(xl))
  Y<-tclVar(tclvalue(yl))
  Z<-tclVar(tclvalue(zl))
  wRange<-round(RNiftyReg::voxelToWorld(matrix(c(1,1,1,d1,d2,d3),byrow = TRUE,ncol=3),orig))
  worldInit<-round(RNiftyReg::voxelToWorld(c(xx,yy,zz),orig))
  Xw<-tclVar(worldInit[1])
  Yw<-tclVar(worldInit[2])
  Zw<-tclVar(worldInit[3])
  parPlotSize1<-c()
  usrCoords1<-c()
  low<-tclVar(as.character(r[1]))
  high<-tclVar(as.character(r[2]))
  low2<-tclVar(as.character(r2[1]))
  high2<-tclVar(as.character(r2[2]))
  intens<-tclVar(as.character(round(data[xx,yy,zz,1],digits = 2)))
  if(olay){ro<-round(range(overlay),digits = 2)
  lowo<-tclVar(as.character(ro[1]))
  higho<-tclVar(as.character(ro[2]))
  intenso<-tclVar(as.character(round(overlay[xx,yy,zz,1],digits = 2)))
  }
  crossHairsOn<-TRUE
  tmovie<-tclVar(FALSE)
  after_ID <-""
  someFlag<-TRUE
  cRegMNI<-tclVar("0")
  cReg<-tclVar("0")
  cRegFileN<-tclVar("0")
  cRegFileO<-tclVar("0")
  ovLay<-tclVar("0")
  #########################
  ####### FRAMES ##########
  #########################
  base<-tktoplevel()
  tktitle(base)<-"Display"
  master<-tkframe(parent=base)
  f1<-tkframe(parent = master)
  #########################
  ##### FUNCTIONS #########
  #########################
  plotf<-function(){
    par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
    mat<-matrix(1:4,nrow = 2,ncol = 2,byrow = TRUE)
    nf<-layout(mat, widths = c(w1,w2),heights = c(w3,w2), respect = FALSE)
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
    olim<-c()
    o1<-c()
    o2<-c()
    o3<-c()
    if(any(is.na(lim))){lim<-r;tclvalue(low)<<-r[1];tclvalue(high)<<-r[2]}
    tclvalue(intens)<<-round(data[x,y,z,t],digits = 2)
    if(tclvalue(ovLay)=="1" && olay){
    olim<-c(as.numeric(tclvalue(lowo)),as.numeric(tclvalue(higho)))
    if(any(is.na(olim))){olim<-ro;tclvalue(lowo)<<-ro[1];tclvalue(higho)<<-ro[2]}
    tclvalue(intenso)<<-round(overlay[x,y,z,1],digits = 2)
    o1<-overlay[,y,,1]
    o2<-overlay[x,,,1]
    o3<-overlay[,,z,1]
    }
    im1<-data[,y,,t]
    im1[im1<lim[1]]<-lim[1]
    im1[im1>lim[2]]<-lim[2]
    o1[o1>olim[2]]<-olim[2]
    image(1:d1,1:d3,z=im1,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(tclvalue(ovLay)=="1" && olay){image(1:d1,1:d3,z=o1,useRaster=FALSE,col=hotMetal(100),axes=FALSE,add=TRUE,zlim = olim)}
    if(crossHairsOn) abline(h = zl,v = xl,col="green")
    im2<-data[x,,,t]
    im2[im2<lim[1]]<-lim[1]
    im2[im2>lim[2]]<-lim[2]
    o2[o2>olim[2]]<-olim[2]
    image(1:d2,1:d3,z=im2,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(tclvalue(ovLay)=="1" && olay){image(1:d2,1:d3,z=o2,useRaster=FALSE,col=hotMetal(100),axes=FALSE,add=TRUE,zlim = olim)}
    if(crossHairsOn) abline(h = zl,v = yl,col="green")
    im3<-data[,,z,t]
    im3[im3<lim[1]]<-lim[1]
    im3[im3>lim[2]]<-lim[2]
    o3[o3>olim[2]]<-olim[2]
    image(1:d1,1:d2,z=im3,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(tclvalue(ovLay)=="1" && olay){image(1:d1,1:d2,z=o3,useRaster=FALSE,col=hotMetal(100),axes=FALSE,add=TRUE,zlim = olim)}
    if(crossHairsOn) abline(h = yl,v = xl,col="green")
    if(d4>1 && tclvalue(ovLay)=="0"){
      tseries<-data[x,y,z,]
      plot(1:d4,tseries,type="l",axes=FALSE,col="grey",xlim=c(1,d4))
      points(x = t,y = tseries[t],col="green")
    }
    if(tclvalue(ovLay)=="1" && olay){
      ran<-round(seq(olim[1],to = olim[2],length.out = 5),digits = 2)
      image(as.matrix(1),useRaster=TRUE,col="black",axes=FALSE)
      rasterImage(image = rev(hotMetal(32)),xleft = -.1,xright = .1,ybottom = -.7,ytop = .7)
      lines(x = c(-.15,-.15),y = c(-.7,.7),col="white")
      lines(x = c(-.2,-.15),y = c(-.7,-.7),col="white")
      lines(x = c(-.2,-.15),y = c(-.35,-.35),col="white")
      lines(x = c(-.2,-.15),y = c(0,0),col="white")
      lines(x = c(-.2,-.15),y = c(.7,.7),col="white")
      lines(x = c(-.2,-.15),y = c(.35,.35),col="white")
      text(x = rep(-.5,5),y = seq(-.7,.7,length.out = 5),labels = as.character(ran),col = "white")
      }
    parPlotSize1 <<- par("plt")
    usrCoords1   <<- par("usr")
    
    wc<-round(RNiftyReg::voxelToWorld(points = c(xl,yl,zl),image = orig))
    if(exists("coXw")){
    tkset(widget = coXw,as.character(wc[1]))
    tkset(widget = coYw,as.character(wc[2]))
    tkset(widget = coZw,as.character(wc[3]))
    }
  }
  plotf2<-function(){
    par(oma = rep(0, 4), mar = rep(0, 4), bg = "black")
    mat<-matrix(1:4,nrow = 2,ncol = 2,byrow = TRUE)
    nf<-layout(mat, widths = c(w1,w2),heights = c(w3,w2), respect = FALSE)
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
    if(t>dim(otherData)[4]){t<-1}
    lim<-c(as.numeric(tclvalue(low2)),as.numeric(tclvalue(high2)))
    if(any(is.na(lim))){lim<-r2;tclvalue(low2)<<-r2[1];tclvalue(high2)<<-r2[2]}
    im1<-otherData[,y,,t]
    im1[im1<lim[1]]<-lim[1]
    im1[im1>lim[2]]<-lim[2]
    image(1:d1,1:d3,z=im1,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = zl,v = xl,col="green")
    im2<-otherData[x,,,t]
    im2[im2<lim[1]]<-lim[1]
    im2[im2>lim[2]]<-lim[2]
    image(1:d2,1:d3,z=im2,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = zl,v = yl,col="green")
    im3<-otherData[,,z,t]
    im3[im3<lim[1]]<-lim[1]
    im3[im3>lim[2]]<-lim[2]
    image(1:d1,1:d2,z=im3,useRaster=TRUE,col=grey(1:100/100),axes=FALSE,zlim=lim)
    if(crossHairsOn) abline(h = yl,v = xl,col="green")
    d4<-dim(otherData)[4]
    if(d4>1){
      tseries<-otherData[x,y,z,]
      plot(1:d4,tseries,type="l",axes=FALSE,col="grey",xlim=c(1,d4))
      points(x = t,y = tseries[t],col="green")
    }
    parPlotSize1 <<- par("plt")
    usrCoords1   <<- par("usr")
  }
  OnLeftClick1 <- function(x,y){
    xClick<-x
    yClick<-y
    width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img1)))-2
    height <- as.numeric(tclvalue(tkwinfo("reqheight",img1)))-2
    xBorder<-w1/(w1+w2)*(width)
    yBorder<-w3/(w3+w2)*(height)
    xClick <- as.numeric(xClick)
    yClick <- as.numeric(yClick)
    if(xClick<=xBorder && xClick>0 && yClick>0 && yClick<=height ){
      xl<-xClick/xBorder*d1+.5
      if(xl>d1){xl<-d1}
      if(xl<1){xl<-1}
      tclvalue(xl)<<-xl
      tclvalue(X)<<-round(xl)
    }
    if(yClick<=yBorder && yClick>0&& xClick>0 &&xClick<width){
      zl<-d3-yClick/yBorder*d3+.5
      if(zl>d3){zl<-d3}
      if(zl<1){zl<-1}
      tclvalue(zl)<<-zl
      tclvalue(Z)<<-round(zl)
    }
    if(xClick>xBorder && xClick<width && yClick<=yBorder && yClick>0){
      yl<-(xClick-xBorder)/(width-xBorder)*d2+.5
      if(yl>d2){yl<-d2}
      if(yl<1){yl<-1}
      tclvalue(yl)<<-yl
      tclvalue(Y)<<-round(yl)
    }
    if(yClick>yBorder && xClick<xBorder && yClick<=height && xClick>0){
      yl<-d2-(yClick-yBorder)/(height-yBorder)*d2+.5
      if(yl>d2){yl<-d2}
      if(yl<1){yl<-1}
      tclvalue(yl)<<-yl
      tclvalue(Y)<<-round(yl)
    }
    if(yClick>yBorder&&xClick>xBorder &&d4>1 &&xClick<width && yClick<=height && tclvalue(ovLay)=="0"){
      first<-1.04-.04*d4
      last<- -.04+1.04*d4
      s<-round(first):round(last)
      rangeUsr<-last-first
      rangePix<-width-xBorder
      click<-xClick-xBorder
      tind<-round(click/rangePix*length(s))
      t<-s[tind]
      if(length(t)<1||t<1){t<-1}
      if(t>d4){t<-d4}
      tclvalue(time)<<-t
    }
    tkrplot::tkrreplot(img1)
    if(tclvalue(cReg)=="1"|tclvalue(cRegMNI)=="1"){tkrplot::tkrreplot(img2)}
  }
  changeTemp<-function(){
    if(is.na(strtoi(tclvalue(time)))){
      tclvalue(time)<<-as.character(1)
    }
    t<-strtoi(tclvalue(time))
    if(t>d4){tclvalue(time)<<-as.character(d4)}
    if(t<1){tclvalue(time)<<-as.character(1)}
    tkrplot::tkrreplot(img1)
    if(tclvalue(cReg)=="1"|tclvalue(cRegMNI)=="1" & dim(otherData)[4]==d4){tkrplot::tkrreplot(img2)}
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
    xl<-as.numeric(tclvalue(xl))
    yl<-as.numeric(tclvalue(yl))
    zl<-as.numeric(tclvalue(zl))
    x<-round(xl)
    y<-round(yl)
    z<-round(zl)
    if(x<1){tclvalue(xl)<<-1}
    if(y<1){tclvalue(yl)<<-1}
    if(z<1){tclvalue(zl)<<-1}
    if(x>d1){tclvalue(xl)<<-d1}
    if(y>d2){tclvalue(yl)<<-d2}
    if(z>d3){tclvalue(zl)<<-d3}
    tkrplot::tkrreplot(img1)
    if(tclvalue(cReg)=="1"|tclvalue(cRegMNI)=="1"){
      if(dim(otherData)[4]>1){
        tkrplot::tkrreplot(img2)
      }
    }
  }
  onSpinW<-function(){
    suppressWarnings(xw<-as.numeric(tclvalue(Xw)))
    suppressWarnings(yw<-as.numeric(tclvalue(Yw)))
    suppressWarnings(zw<-as.numeric(tclvalue(Zw)))
    voxCo<-RNiftyReg::worldToVoxel(points = c(xw,yw,zw),image = orig)
    
  if(is.na(voxCo[1])){
      tclvalue(xl)<<-as.character(xx)
    }else{
      xl<-voxCo[1]
      x<-round(xl)
      if(x<1||x>d1){tclvalue(xl)<<-1}else{tclvalue(xl)<<-xl}
    }
    
    if(is.na(voxCo[2])){
      tclvalue(yl)<<-as.character(yy)
    }else{
      yl<-voxCo[2]
      y<-round(yl)
      if(y<1||y>d2){tclvalue(yl)<<-1}else{tclvalue(yl)<<-yl}
    }
    
    if(is.na(voxCo[3])){
      tclvalue(zl)<<-as.character(zz)
    }else{
      zl<-voxCo[3]
      z<-round(zl)
      if(z<1||z>d3){tclvalue(zl)<<-1}else{tclvalue(zl)<<-zl}
    }
    
    tkrplot::tkrreplot(img1)
    if(tclvalue(cReg)=="1"|tclvalue(cRegMNI)=="1"){
      if(dim(otherData)[4]>1){
        tkrplot::tkrreplot(img2)
      }
    }
  }
  maxminChange<-function(){
    if(suppressWarnings(is.na(as.numeric(tclvalue(high))))||tclvalue(high)==""){
      tclvalue(high)<<-r[2]
    }
    if(suppressWarnings(is.na(as.numeric(tclvalue(low))))||tclvalue(low)==""){
      tclvalue(low)<<-r[1]
    }
    max<-as.numeric(tclvalue(high))
    min<-as.numeric(tclvalue(low))
    if(max<min){tclvalue(high)<<-r[2];tclvalue(low)<<-r[1]}
    tkrplot::tkrreplot(img1)
  }
  if(olay){
  maxminChangeO<-function(){
    if(suppressWarnings(is.na(as.numeric(tclvalue(higho))))||tclvalue(higho)==""){
      tclvalue(higho)<<-ro[2]
    }
    if(suppressWarnings(is.na(as.numeric(tclvalue(lowo))))||tclvalue(lowo)==""){
      tclvalue(lowo)<<-ro[1]
    }
    max<-as.numeric(tclvalue(higho))
    min<-as.numeric(tclvalue(lowo))
    if(max<min){tclvalue(higho)<<-ro[2];tclvalue(lowo)<<-ro[1]}
    tkrplot::tkrreplot(img1)
  }
  }
  crossHairs<-function(){
    crossHairsOn<<-!crossHairsOn
    tkrplot::tkrreplot(img1)
    if(tclvalue(cReg)=="1"|tclvalue(cRegMNI)=="1"){tkrplot::tkrreplot(img2)}
  }
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
      repeat_call(1,function() {tkinvoke(spin,"buttonup")})
    }else{someFlag<<-FALSE}
  }
  checkRegMNI<-function(){
    if(tclvalue(cReg)=="1"){tclvalue(cRegMNI)<<-"0";return(0)}
    if(tclvalue(cRegMNI)=="1"){
      otherData<<-readNii(system.file("extdata","mni.nii.gz",package = "FIACH"))
      dim(otherData)<<-c(dim(otherData),1)
      if(any(dim(data)[1:3]!=dim(otherData)[1:3])){tclvalue(cRegMNI)<<-"0";return(0)}
      r2<<-round(range(otherData),digits = 2)
      low2<<-tclVar(as.character(r2[1]))
      high2<<-tclVar(as.character(r2[2]))
      if(asp>1){
        img2<<-tkrplot::tkrplot(parent = f1,fun = plotf2,hscale=scaleFactor,vscale=scaleFactor/asp)
      }else{
        img2<<-tkrplot::tkrplot(parent = f1,fun = plotf2,hscale=scaleFactor*asp,vscale=scaleFactor)
      }
      tkgrid(f6,column=1,row=1)
      tkgrid(img2,column=1,row=0)
      tkbind(img2, "<Button-3>",crossHairs)
      tkbind(img2, "<Button-1>",OnLeftClick1)
      tkbind(img2, "<B1-Motion>",OnLeftClick1)
      tkconfigure(img2,cursor="crosshair")
    }else{
      tkgrid.forget(img2)
      tkgrid.forget(f1)
      tkgrid.forget(f6)
      tkgrid(f1,columnspan=2,rowspan=1,column=0,row=0)
      tkgrid(img1,column=0,row=0)
    }
  }
  checkReg<-function(){
    if(tclvalue(cRegMNI)=="1"|sum(origOther)==0){tclvalue(cReg)<<-"0";return(0)}
    if(tclvalue(cReg)=="1"){
      dim(otherData)<<-c(dim(otherData),1)
      if(any(dim(data)[1:3]!=dim(otherData)[1:3])){tclvalue(cReg)<<-"0";return(0)}
      otherData<<-origOther
      r2<<-round(range(otherData),digits = 2)
      low2<<-tclVar(as.character(r2[1]))
      high2<<-tclVar(as.character(r2[2]))
      if(asp>1){
        img2<<-tkrplot::tkrplot(parent = f1,fun = plotf2,hscale=scaleFactor,vscale=scaleFactor/asp)
      }else{
        img2<<-tkrplot::tkrplot(parent = f1,fun = plotf2,hscale=scaleFactor*asp,vscale=scaleFactor)
      }
      tkgrid(f6,column=1,row=1)
      tkgrid(img2,column=1,row=0)
      tkbind(img2, "<Button-3>",crossHairs)
      tkbind(img2, "<Button-1>",OnLeftClick1)
      tkbind(img2, "<B1-Motion>",OnLeftClick1)
      tkconfigure(img2,cursor="crosshair")
      
    }else{
      tkgrid.forget(img2)
      tkgrid.forget(f1)
      tkgrid.forget(f6)
      tkgrid(f1,columnspan=2,rowspan=1,column=0,row=0)
      tkgrid(img1,column=0,row=0)
    }
  }
  savePNG<-function(){
    Jwin<-tktoplevel()
    width <- tclVar(tkwinfo("reqwidth",img1))
    height<-tclVar(tkwinfo("reqheight",img1))
    res<-tclVar("300")
    WIDTH<-tkentry(Jwin,width="3",textvariable=width)
    HEIGHT<-tkentry(Jwin,width="3",textvariable=height)
    RES<-tkentry(Jwin,width="3",textvariable=res)
    Wlab<-tklabel(Jwin,text="Width")
    Hlab<-tklabel(Jwin,text="Height")
    Rlab<-tklabel(Jwin,text="DPI")
    Wlab2<-tklabel(Jwin,text="px")
    Hlab2<-tklabel(Jwin,text="px")
    tkgrid(Hlab,HEIGHT,Hlab2)
    tkgrid(Wlab,WIDTH,Wlab2)
    tkgrid(Rlab,RES)
    tkfocus(Jwin)
    save<-function(){
      w<-as.numeric(tclvalue(width))
      h<-as.numeric(tclvalue(height))
      r<-as.numeric(tclvalue(res))
      png(filename = tclvalue(tkgetSaveFile()),width = w,height = h,res =r)
      plotf()
      dev.off()
      tkdestroy(Jwin)
    }
    savBut <-tkbutton(Jwin,text="  SAVE   ",command=save)
    tkgrid(savBut)
  }
  overLay<-function(){
    tkrplot::tkrreplot(img1)
  }
  
  ##########################
  ###### THE PLOTS #########
  ##########################
  testImg<-tkrplot::tkrplot(parent = f1,fun =plotf)
  testHeight<-as.numeric(tkwinfo("reqheight",testImg))
  testWidth<-as.numeric(tkwinfo("reqwidth",testImg))
  scaleFactor<-(desHeight)/(testHeight)
  asp<-(w1+w2)/(w2+w3)
  if(asp>1){
    img1<-tkrplot::tkrplot(parent = f1,fun = plotf,hscale=scaleFactor,vscale=scaleFactor/asp)
    img2<-tkrplot::tkrplot(parent = f1,fun = plotf2,hscale=scaleFactor,vscale=scaleFactor/asp)
  }else{
    img1<-tkrplot::tkrplot(parent = f1,fun = plotf,hscale=scaleFactor*asp,vscale=scaleFactor)
    img2<-tkrplot::tkrplot(parent = f1,fun = plotf2,hscale=scaleFactor*asp,vscale=scaleFactor)
  }
  f5<-tkframe(parent=master,
              width=as.numeric(tkwinfo("reqwidth",img1)),
              height=150,
              borderwidth=0,
              relief="groove")
  f6<-tkframe(parent=master,
              width=as.numeric(tkwinfo("reqwidth",img1)),
              height=150,
              borderwidth=0,
              relief="groove")
  tkgrid(master,columnspan=2,rowspan=2)
  tkgrid(f1,columnspan=2,rowspan=1)
  tkgrid(img1,column=0,row=0)
  ##########################
  ##### Widgets ############
  ##########################
  spin<-tkspinbox(f5,textvariable=time,from=1,to=d4,command=changeTemp,increment=1,repeatdelay=200,width=5,wrap=TRUE)
  coX<-tkspinbox(f5,textvariable=X,from=1,to=d1,command=onSpin,increment=1,repeatdelay=10,width=5)
  coY<-tkspinbox(f5,textvariable=Y,from=1,to=d2,command=onSpin,increment=1,repeatdelay=10,width=5)
  coZ<-tkspinbox(f5,textvariable=Z,from=1,to=d3,command=onSpin,increment=1,repeatdelay=10,width=5)
  coXw<-tkspinbox(f5,textvariable=Xw,values=wRange[1,1]:wRange[2,1],repeatdelay=10,width=5,command=onSpinW)
  coYw<-tkspinbox(f5,textvariable=Yw,values=wRange[1,2]:wRange[2,2],repeatdelay=10,width=5,command=onSpinW)
  coZw<-tkspinbox(f5,textvariable=Zw,values=wRange[1,3]:wRange[2,3],repeatdelay=10,width=5,command=onSpinW)
  MAX<-tkentry(f5,textvariable=high,width=7)
  MIN<-tkentry(f5,textvariable=low,width=7)
  if(olay){
    intensityo<-tkentry(f5,textvariable=intenso,state="readonly",readonlybackground="white",width=7)
    MAXo<-tkentry(f5,textvariable=higho,width=7)
    MINo<-tkentry(f5,textvariable=lowo,width=7)
  }
  intensity<-tkentry(f5,textvariable=intens,state="readonly",readonlybackground="white",width=7)
  coXLab<-tklabel(parent = f5,text="X")
  coYLab<-tklabel(parent = f5,text="Y")
  coZLab<-tklabel(parent = f5,text="Z")
  coTLab<-tklabel(parent = f5,text="Time")
  maxLab<-tklabel(parent = f5,text="Max")
  minLab<-tklabel(parent = f5,text="Min")
  intensLab<-tklabel(parent=f5,text="Intensity")
  voxLab<-tklabel(parent=f5,text="Voxel")
  worldLab<-tklabel(parent=f5,text="World")
  baseLab<-tklabel(parent=f5,text="Base")
  oLab<-tklabel(parent=f5,text="Overlay")
  tmovieBut<-tkcheckbutton(f5,variable=tmovie, command=movieT,text="Movie")
  ##########################
  ###### GEOMETRY ##########
  ##########################
  tkgrid(f5,column=0,row=1)
  tkplace(coXLab,relx=0,rely=1/6)
  tkplace(coYLab,relx=0,rely=2/6)
  tkplace(coZLab,relx=0,rely=3/6)
  tkplace(coTLab,relx=0,rely=4/6)
  tkplace(coX,'in'=coXLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(coY,'in'=coYLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(coZ,'in'=coZLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(spin,'in'=coTLab,x=as.numeric(tkwinfo("reqwidth",coTLab)),rely=.5,anchor="w",height=100/5)
  tkplace(coXw,'in'=coX,x=as.numeric(tkwinfo("reqwidth",coX))*1.1,rely=.5,anchor="w",height=100/5)
  tkplace(coYw,'in'=coY,x=as.numeric(tkwinfo("reqwidth",coY))*1.1,rely=.5,anchor="w",height=100/5)
  tkplace(coZw,'in'=coZ,x=as.numeric(tkwinfo("reqwidth",coZ))*1.1,rely=.5,anchor="w",height=100/5)
  tkplace(voxLab,"in"=coX,x=0,y=-as.numeric(tkwinfo("reqheight",coX))*.75,anchor="w")
  tkplace(worldLab,"in"=coXw,x=0,y=-as.numeric(tkwinfo("reqheight",coXw))*.75,anchor="w")
  if(d4>1){tkplace(tmovieBut,"in"=spin,rely=.5,anchor="w",x=as.numeric(tkwinfo("reqwidth",spin)),height=80/4)}
  wiBox<-as.numeric(tkwinfo("reqwidth",intensity))+10
  tkplace(intensity,relx=1,x=-wiBox*2,rely=1/5)
  tkplace(MAX,relx=1,x=-wiBox*2,rely=2/5)
  tkplace(MIN,relx=1,x=-wiBox*2,rely=3/5)
  if(olay){
  tkplace(intensityo,relx=1,x=-wiBox,rely=1/5)
  tkplace(MAXo,relx=1,x=-wiBox,rely=2/5)
  tkplace(MINo,relx=1,x=-wiBox,rely=3/5)
  tkplace(oLab,"in"=intensityo,x=-5,y=-as.numeric(tkwinfo("reqheight",intensityo))*.75,anchor="w")
  tkbind(MAXo, "<Return>",function()maxminChangeO())
  tkbind(MINo, "<Return>",function()maxminChangeO())
  }
  tkplace(intensLab,'in'=intensity,x=-1,rely=.5,anchor="e",height=80/4)
  tkplace(maxLab,'in'=MAX,x=-1,rely=.5,anchor="e",height=80/4)
  tkplace(minLab,'in'=MIN,x=-1,rely=.5,anchor="e",height=80/4)
  tkplace(baseLab,"in"=intensity,x=0,y=-as.numeric(tkwinfo("reqheight",intensity))*.75,anchor="w")
  ###########################
  ######## BINDINGS #########
  ###########################
  tkbind(MAX, "<Return>",function()maxminChange())
  tkbind(MIN, "<Return>",function()maxminChange())
  tkbind(coX, "<Return>",function()onSpin())
  tkbind(coY, "<Return>",function()onSpin())
  tkbind(coZ, "<Return>",function()onSpin())
  tkbind(coXw, "<Return>",function()onSpinW())
  tkbind(coYw, "<Return>",function()onSpinW())
  tkbind(coZw, "<Return>",function()onSpinW())
  tkbind(spin, "<Return>",function()changeTemp())
  tkbind(img1, "<Button-3>",crossHairs)
  tkbind(img1, "<Button-1>",OnLeftClick1)
  tkbind(img1, "<B1-Motion>",OnLeftClick1)
  ###########################
  ##### Configure ###########
  ###########################
  tkconfigure(img1,cursor="crosshair")
  tkwm.resizable(base,FALSE,FALSE)
  tkset(widget = coXw,as.character(worldInit[1]))
  tkset(widget = coYw,as.character(worldInit[2]))
  tkset(widget = coZw,as.character(worldInit[3]))
  ###########################
  ####### MENUS #############
  ###########################
  origWarn<-getOption(x = "warn")
  options(warn = -1)
  topMenu <- tkmenu(base,tearoff=FALSE)
  tkconfigure(base, menu = topMenu)
  fileMenu <- tkmenu(topMenu, tearoff = FALSE)
  tkadd(topMenu, "cascade", label = "File", menu = fileMenu)
  saveAsMenu <- tkmenu(topMenu, tearoff = FALSE)  # Our cacaded menu
  tkadd(fileMenu, "cascade", label = "Save as", menu = saveAsMenu)
  overlayMenu <- tkmenu(topMenu, tearoff = FALSE)  # Our cascaded menu
  tkadd(topMenu,"cascade", label = "Overlay",menu=overlayMenu)
  checkRegMenu <- tkmenu(topMenu, tearoff = FALSE)  # Our cascaded menu
  tkadd(topMenu,"cascade", label = "Check Reg",menu=checkRegMenu)
  tkadd(saveAsMenu, "command", label = ".png",command=savePNG)
  tkadd(checkRegMenu, "checkbutton",label="MNI 2mm Iso", variable=cRegMNI, onvalue=1 ,offvalue=0,command=checkRegMNI)
  tkadd(overlayMenu, "checkbutton", label = "Image", variable=ovLay, onvalue=1 ,offvalue=0,command=overLay)
  if(!null2){tkadd(checkRegMenu, "checkbutton",label="Image", variable=cReg, onvalue=1 ,offvalue=0,command=checkReg)}
  suppressWarnings(options(warn=1))
  suppressWarnings(options(warn=origWarn))
  if(ret){return(orig)}
}
