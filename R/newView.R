viewNew<-function(data=NULL){
  ###################################
  ######### ARG CHECK ###############
  ###################################
  if(is.null(data)){
    file<-selectR()
    if(length(file)==0){stop("You must select a file")}
    func<-readNii(file)
  } else if(is.character(data)){
    file<-data
    if(length(file)==0){stop("You must select a file")}
    func<-readNii(file)
  }else if(!is.array(data)){
    stop("This is not an array with 3 or more dimensions")
  }else{
    func<-data
  }  
  #olay<-readNifti(file)
  func<-zeroNa(func)
  ##################################
  ######## ASPECT RATIO ############
  ##################################
  
  # pixel dimensions needed to get aspect ratio right
  pix<-RNiftyReg::pixdim(func)[1:3]
  # vector containing the aspect ratios for each image
  asp<-c(pix[2]/pix[3],pix[1]/pix[3],pix[1]/pix[2])
  
  # make all image 4d as it is then easier to handle truly 4d images without writing new code
  d<-dim(func)
  if(length(d)==3){
    dim(func)<-c(dim(func),1)
    #dim(olay)<-c(dim(olay),1)
  }
  d<-dim(func)
  # get xorm for transforming to world coordinates
  xf<-RNiftyReg::xform(func)
  ##################################
  ########## FRAMES ################
  ##################################
  # 3 image windows and a utilities window in a master window 
  top <- tktoplevel()
  tktitle(top)<-"Display"
  master<-tkframe(parent=top)
  img<-tkframe(parent = master)
  f1<-tkframe(parent = img)
  f2<-tkframe(parent = img)
  f3<-tkframe(parent = img)
  f4<-tkframe(parent = img)
  f5<-tkframe(parent = master,borderwidth=5,relief="groove")
  myfont <- tkfont.create(family="Arial",size=9)
  ##################################
  ####### GLOBAL VARIABLES #########
  ##################################
  
  # initial coordinate is middle of image
  xyz<-round(c(d[1]/2,d[2]/2,d[3]/2))
  # image expands to 256 pixels or if naturally larger stays as is
  zoom<-pmax(floor(256/max(d[1:3])),1)
  
  #empty vectors to hold crosshair poistion and length
  xyzL<-xyz
  xyzLineLength<-c()
  
  # time variable initialised to 1
  t<-1
  # no subsampling  of image takes place
  subsamp<-1
  # by default crosshairs are on
  crosshairsOn<-TRUE
  
  # tcl varaibles to hold voxel coordinates and time
  X<-tclVar(xyz[1])
  Y<-tclVar(xyz[2])
  Z<-tclVar(xyz[3])
  time<-tclVar(t)
  
  #get range to initialise max and min widget
  r<-range(func,na.rm = TRUE)
  #or<-oRANGE #commented out support for overlay
  
  #tcl varaibles to hold max,min and intensity
  high<-tclVar(r[2])
  low<-tclVar(r[1])
  #ohigh<-tclVar(or[2]) #commented out support for overlay
  #olow<-tclVar(or[1]) #commented out support for overlay
  intens<-tclVar(as.character(round(func[xyz[1],xyz[2],xyz[3],1],digits = 4)))
  #ointens<-tclVar(as.character(round(olay[xyz[1],xyz[2],xyz[3],1],digits = 4))) #commented out support for overlay
  
  # varaibles to describe which window is being clicked upon
  click1<-FALSE
  click2<-FALSE
  click3<-FALSE
  
  # set initial world coordinates and hold them in tcl variables
  worldInit<-round(RNiftyReg::voxelToWorld(xyz,func))
  Xw<-tclVar(worldInit[1])
  Yw<-tclVar(worldInit[2])
  Zw<-tclVar(worldInit[3])
  
  # initialise image matrices to black
  im3<-matrix(data = "#ffffff",nrow = d[1],ncol = d[2])
  im2<-matrix(data = "#ffffff",nrow = d[1],ncol = d[3])
  im1<-matrix(data = "#ffffff",nrow = d[2],ncol = d[3])
  
  # concatenate along columns.. the ugly hack begins
  # effectively the vector elements are rows of pixels
  ccim3<-.concat1(im3,margin = 2)
  ccim2<-.concat1(im2,margin = 2)
  ccim1<-.concat1(im1,margin = 2)
  
  # collapse the object into one string separated by curly braces
  p3<-paste(get("ccim3"),collapse = " } {",sep="")
  cmdz<-paste("{",p3," }",sep = "")
  
  # repeat for the Y image 
  p2<-paste(get("ccim2"),collapse = " } {",sep="")
  cmdy<-paste("{",p2," }",sep = "")
  
  # repeat for the Z image 
  p1<-paste(get("ccim1"),collapse = " } {",sep="")
  cmdx<-paste("{",p1," }",sep = "")
  
  #set the colour palette to greyscale
  palette<-grey(0:255/255)
  
  tmovie<-tclVar(FALSE)
  after_ID <-""
  someFlag<-TRUE
  ##################################
  ####### FUNCTIONS ################
  ##################################
  # function designed to approximate a real number with a rational number
  ratApprox<-function(asp,maxRat = 10){
    zooms<-round(asp*1:maxRat)
    subsamps<-1:maxRat
    subsamp<-which.min(abs(zooms/subsamps-asp))
    zoom<-zooms[subsamp]
    return(c(zoom,subsamp))
  }
  
  # prepares images to be displayed if clicked
  onLeftClick1<-function(x,y){
    
    coord<-as.numeric(c(x,y))
    height<-as.character(tkwinfo("reqheight",f1))
    width<-as.character(tkwinfo("reqwidth",f1))
    xyz[2]<<-round(coord[1]/as.numeric(width)*d[2])
    xyz[3]<<-d[3]-round(coord[2]/as.numeric(height)*d[3])
    xyzL[2]<<-xyzLineLength[2]-coord[1]+2
    xyzL[3]<<-coord[2]-2
    xyzL[xyzL<1]<<-1
    xyzL[xyzL>xyzLineLength]<<-d[xyzL>xyzLineLength]
    click1<<-TRUE
    reslicer(c(xyz[1],xyz[2],xyz[3]),zoom,subsamp)
    click1<<-FALSE
    
  }
  
  # prepares images to be displayed if clicked
  onLeftClick2<-function(x,y){
    
    coord <- as.numeric(c(x,y))
    height <- as.character(tkwinfo("reqheight",f2))
    width <- as.character(tkwinfo("reqwidth",f2))
    xyz[1] <<- round(coord[1] / as.numeric(width) * d[1])
    xyz[3] <<- d[3] - round(coord[2] / as.numeric(height) * d[3])
    xyzL[1] <<- coord[1] - 2
    xyzL[3] <<- coord[2] - 2
    xyzL[xyzL<1]<<-1
    xyzL[xyzL>xyzLineLength]<<-d[xyzL>xyzLineLength]
    click2<<-TRUE
    reslicer(xyz,zoom,subsamp)
    click2<<-FALSE
  }
  
  # prepares images to be displayed if clicked
  onLeftClick3<-function(x,y){
    coord<-as.numeric(c(x,y))
    height<-as.character(tkwinfo("reqheight",f3))
    width<-as.character(tkwinfo("reqwidth",f3))
    xyz[1]<<-round(coord[1]/as.numeric(width)*d[1])
    xyz[2]<<-d[2]-round(coord[2]/as.numeric(height)*d[2])
    xyzL[1]<<-coord[1]-2
    xyzL[2]<<-coord[2]-2
    xyzL[xyzL<1]<<-1
    xyzL[xyzL>xyzLineLength]<<-d[xyzL>xyzLineLength]
    click3<<-TRUE
    reslicer(c(xyz[1],xyz[2],xyz[3]),zoom,subsamp)
    click3<<-FALSE
  }
  onLeftClick4<-function(x,y){
    coord<-as.numeric(c(x,y))
    cw<-as.numeric(tkwinfo("reqwidth",canvas))
    x<-0:(d[4]-1)
    xfac<-cw/max(x)
    xnew<-x*xfac
    tclvalue(time)<<-which.min(abs(xnew-coord[1]))
    reslicer(c(xyz[1],xyz[2],xyz[3]),zoom,subsamp)
  }
  
  # turns crosshairs on and off
  crossHairs<-function(){
    crosshairsOn<<-!crosshairsOn
    reslicer(coNew = xyz,zoom = zoom,subsamp = subsamp)
    gc()
  }
  
  # displays image in response to coordinate being manually changed
  onSpin<-function(){
    xyz[1]<<-as.numeric(tclvalue(X))
    xyz[2]<<-as.numeric(tclvalue(Y))
    xyz[3]<<-as.numeric(tclvalue(Z))
    def<-round(d[1:3]/2)
    xyz[is.na(xyz)]<<-def[is.na(xyz)]
    xyzL<<-round((c(0,d[2],d[3])-xyz)/d[1:3]*xyzLineLength*c(-1,1,1))
    
    reslicer(coNew = xyz,zoom = zoom,subsamp = subsamp)
  }
  
  # redisplays image with max and min changed
  maxminChange<-function(){
    ma<-as.numeric(tclvalue(high))
    mi<-as.numeric(tclvalue(low))
    if(is.na(ma)){tclvalue(high)<<-r[2];ma<-r[2]}
    if(is.na(mi)){tclvalue(low)<<-r[1];mi<-r[1]}
    reslicer(coNew = xyz,zoom = zoom,subsamp = subsamp)
  }
  
  # redisplays image with overlay max and min changed
  # omaxminChange<-function(){
  #   ma<-as.numeric(tclvalue(ohigh))
  #   mi<-as.numeric(tclvalue(olow))
  #   if(is.na(ma)){tclvalue(ohigh)<<-or[2];ma<-or[2]}
  #   if(is.na(mi)){tclvalue(olow)<<-or[1];mi<-or[1]}
  #   maxminChange()
  #   olay[olay>ma]<-ma
  #   check<-olay>mi
  #   olay[olay<mi]<-mi
  #   scfunc<-scaleRGB(olay,max = ma,min = mi)
  #   RGBarr[check]<<-hextest(input = olay,palette = hotMetal(256),currentmax = ma,currentmin = mi)[check]
  #   
  #   reslicer(coNew = xyz,zoom = zoom,subsamp = subsamp)
  # }
  
  # adds the spinbox widget which can be mysteriously absent...
  tkspinbox <- function(parent, ...) {tkwidget(parent, "tk::spinbox", ...)}
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
      repeat_call(1,function() {tkinvoke(coT,"buttonup")})
    }else{someFlag<<-FALSE}
  }
  ##################################
  ######## INITIAL #################
  ##################################
  xlabel<-sample(0:9,size = 10,replace = TRUE)
  ylabel<-sample(0:9,size = 10,replace = TRUE)
  zlabel<-sample(0:9,size = 10,replace = TRUE)
  hmmx<-paste("h",paste(xlabel,collapse = ""),sep = "")
  hmmy<-paste("h",paste(ylabel,collapse = ""),sep = "")
  hmmz<-paste("h",paste(zlabel,collapse = ""),sep = "")
  shx<-paste("sh",paste(xlabel,collapse = ""),sep = "")
  shy<-paste("sh",paste(ylabel,collapse = ""),sep = "")
  shz<-paste("sh",paste(zlabel,collapse = ""),sep = "")
  rhx<-paste("rh",paste(xlabel,collapse = ""),sep = "")
  rhy<-paste("rh",paste(ylabel,collapse = ""),sep = "")
  rhz<-paste("rh",paste(zlabel,collapse = ""),sep = "")
  
  # Create the base x,y and z photo images
  #.Tcl("image create photo hmmx")
  # .Tcl("image create photo hmmy")
  # .Tcl("image create photo hmmz")
  .Tcl(paste("image create photo",hmmx))
  .Tcl(paste("image create photo",hmmy))
  .Tcl(paste("image create photo",hmmz))
  
  # Create the intermediary x,y and z photo images
  # .Tcl("image create photo shmmx")
  # .Tcl("image create photo shmmy")
  # .Tcl("image create photo shmmz")
  .Tcl(paste("image create photo",shx))
  .Tcl(paste("image create photo",shy))
  .Tcl(paste("image create photo",shz))
                 
  # Create the final x,y and z photo images
  # xImageCallback<-.Tcl("image create photo rhmmx")
  # yImageCallback<-.Tcl("image create photo rhmmy")
  # zImageCallback<-.Tcl("image create photo rhmmz")
  xImageCallback<-.Tcl(paste("image create photo",rhx))
  yImageCallback<-.Tcl(paste("image create photo",rhy))
  zImageCallback<-.Tcl(paste("image create photo",rhz))
                       
  
  # attach the final images to a tk widget so they can be displayed
  lx <- tklabel(f1, image = rhx)
  ly <- tklabel(f2, image = rhy)
  lz <- tklabel(f3, image = rhz)
  
  # create the arguements that will be passed to tclk
  tclArgZ<-.Tcl.args.objv(hmmz,"put",cmdz)
  tclArgY<-.Tcl.args.objv(hmmy,"put",cmdy)
  tclArgX<-.Tcl.args.objv(hmmx,"put",cmdx)
  
  # put the actual data in the image using the arguments above
  zPutCallback<-.Tcl.objv(tclArgZ)
  yPutCallback<-.Tcl.objv(tclArgY)
  xPutCallback<-.Tcl.objv(tclArgX)
  
  # The workhorse function that displays new images
  reslicer<-function(coNew,zoom,subsamp){
    # perform bounds checks on coordinates
    xyz[xyz<1]<<-1
    xyzL[xyzL<1]<-1
    xyz[xyz<1]<<-1
    xyzL[xyzL<1]<-1
    xyz[xyz>d[1:3]]<<-(d[1:3])[xyz>d[1:3]]
    t<<-as.numeric(tclvalue(time))
    t[is.na(t)]<<-1
    t[t>d[4]]<-d[4]
    t[t<1]<-1
    tclvalue(time)<<-t
    
    #convert to world coordinates and update variable
    world<-round(xf%*%c(xyz,1))     
    tclvalue(Xw)<<-world[1]
    tclvalue(Yw)<<-world[2]
    tclvalue(Zw)<<-world[3]
    xyz<<-round(xyz)
    tclvalue(X)<<-xyz[1]
    tclvalue(Y)<<-xyz[2]
    tclvalue(Z)<<-xyz[3]
    
    # get max and min and check if non numeric.. correct if necessary
    ma<-as.numeric(tclvalue(high))
    mi<-as.numeric(tclvalue(low))
    if(is.na(ma)){tclvalue(high)<<-r[2];ma<-r[2]}
    if(is.na(mi)){tclvalue(low)<<-r[1];mi<-r[1]}
    
    
    # get rational approximation for aspect ratio
    indAsp<-which(asp<1)
    asp[indAsp]<-1/asp[indAsp]
    aspAdj1<-ratApprox(asp[1])
    rat1<-ratApprox(asp[1])
    rat2<-ratApprox(asp[2])
    rat3<-ratApprox(asp[3])
    
    # caclualte height and widths of images in pixels
    zoom<<-zoom
    subsamp<<-subsamp
    zoom1<-paste(zoom,zoom*rat1[1])
    zoom2<-paste(zoom,zoom*rat2[1])
    zoom3<-paste(zoom,zoom*rat3[1])
    subsamp1<-paste(subsamp,-subsamp*rat1[2])
    subsamp2<-paste(subsamp,-subsamp*rat2[2])
    subsamp3<-paste(subsamp,-subsamp*rat3[2])
    wz<-round(d[1]*zoom/subsamp)
    hz<-round(d[2]*zoom/subsamp*rat3[1]/rat3[2])
    wy<-round(d[1]*zoom/subsamp)
    hy<-round(d[3]*zoom/subsamp*rat2[1]/rat2[2])
    wx<-round(d[2]*zoom/subsamp)
    hx<-round(d[3]*zoom/subsamp*rat1[1]/rat1[2])
    
    # set height and width  of images
    # .Tcl(paste("rhmmz configure -width",wz,"-height",hz)) #1/2
    # .Tcl(paste("rhmmy configure -width",wy,"-height",hy))#1/2
    # .Tcl(paste("rhmmx configure -width",wx,"-height",hx))#1/2
    .Tcl(paste( rhz, "configure -width",wz,"-height",hz))#1/2
    .Tcl(paste( rhy, "configure -width",wy,"-height",hy))#1/2
    .Tcl(paste( rhx, "configure -width",wx,"-height",hx))#1/2
    
    # create a local copy of image 
    locFunc<-get("func")
    
    # get the intensity at crosshair location
    tclvalue(intens)<<-round(locFunc[xyz[1],xyz[2],xyz[3],t],digits = 3)

    if(!click3){
      # in place modification of im3   
      .hextest(input = locFunc[,,xyz[3],t],palette = palette,currentmax = ma,currentmin = mi,out = im3)
      # in place modification of ugly concatenated string
      .concat2(get("im3"),margin = 2,y=get("ccim3"))
      
      # creat the tcl command in string form
      p3<<-paste(get("ccim3"),collapse = " } {",sep="")
      cmdz<<-paste("{",p3," }",sep = "")
      
      # update tcl argument structure in place with string command
      .tclObject( tclArgZ[[3]],update = cmdz)
      
      # pass the ugly cmd string to tcltk (computational bottleneck)
      zPutCallback<<-.Tcl.objv(tclArgZ)
      
      # zoom by the factor zoom 3
      #.Tcl(paste("shmmz copy hmmz -zoom",zoom3))
       .Tcl(paste(shz,"copy",hmmz, "-zoom",zoom3))
    }
    if(!click2){
      # same as previous but for different image 
      .hextest(input = locFunc[,xyz[2],,t],palette = palette,currentmax = ma,currentmin = mi,out = im2)
      .concat2(get("im2"),margin = 2,y=get("ccim2"))
      p2<<-paste(get("ccim2"),collapse = " } {",sep="")
      cmdy<<-paste("{",p2," }",sep = "")
      .tclObject( tclArgY[[3]],update = cmdy)
      yPutCallback<<-.Tcl.objv(tclArgY)
      #.Tcl(paste("shmmy copy hmmy -zoom",zoom2))
      .Tcl(paste(shy,"copy", hmmy,"-zoom",zoom2))
    }
    if(!click1){
      # same as previous but for different image 
      .hextest(input = locFunc[xyz[1],,,t],palette = palette,currentmax = ma,currentmin = mi,out = im1)
      .concat2(get("im1"),margin = 2,y=get("ccim1"))
      p1<<-paste(get("ccim1"),collapse = " } {",sep="")
      cmdx<<-paste("{",p1," }",sep = "")
      .tclObject( tclArgX[[3]],update = cmdx)
      xPutCallback<<-.Tcl.objv(tclArgX)
      #.Tcl(paste("shmmx copy hmmx -zoom",zoom1))
      #.Tcl(paste("rhmmx copy shmmx  -subsample",subsamp1))
      .Tcl(paste(shx,"copy", hmmx,"-zoom",zoom1))
    }
    
    #final image is  the subsampled  version of the zoomed image 
    # .Tcl(paste("rhmmz copy shmmz -subsample",subsamp3))
    # .Tcl(paste("rhmmy copy shmmy -subsample",subsamp2))
    # .Tcl(paste("rhmmx copy shmmx -subsample",subsamp1))
    .Tcl(paste(rhz,"copy", shz, "-subsample",subsamp3))
    .Tcl(paste(rhy,"copy", shy, "-subsample",subsamp2))
    .Tcl(paste(rhx,"copy", shx, "-subsample",subsamp1))
    
    # places green line at coordinate
     if(crosshairsOn){
      cmd<-paste(rhx,"put  #00FF00 -to",xyzLineLength[2]-xyzL[2],0,xyzLineLength[2]-xyzL[2]+1,xyzLineLength[3])
      .Tcl(cmd)
      cmd<-paste(rhx, "put  #00FF00 -to",0,xyzL[3],xyzLineLength[2],xyzL[3]+1)
      .Tcl(cmd)

      cmd<-paste(rhz,"put  #00FF00 -to",0,xyzL[2],xyzLineLength[1],xyzL[2]-1)
      .Tcl(cmd)
      cmd<-paste(rhz,"put  #00FF00 -to",xyzL[1],0,xyzL[1]+1,xyzLineLength[2])
      .Tcl(cmd)

      cmd<-paste(rhy,"put  #00FF00 -to",0,xyzL[3],xyzLineLength[2],xyzL[3]+1)
      .Tcl(cmd)
      cmd<-paste(rhy,"put  #00FF00 -to",xyzL[1],0,xyzL[1]+1,xyzLineLength[3])
      .Tcl(cmd)}
    # 
    # 
    # update green line lengths and size of bottom frame(potential for resizing)
    xyzLineLength[1]<<-as.numeric(tkwinfo("reqwidth",f3))
    xyzLineLength[2]<<-as.numeric(tkwinfo("reqheight",f3))
    xyzLineLength[3]<<-as.numeric(tkwinfo("reqheight",f1))
    tkconfigure(f5,height=100,width=as.numeric(tkwinfo("width",img)))
    
    if(d[4]>1){
      x<-0:(d[4]-1)
      y<-func[xyz[1],xyz[2],xyz[3],]
      ch<-as.numeric(tkwinfo("reqheight",canvas))
      cw<-as.numeric(tkwinfo("reqwidth",canvas))
      xfac<-cw/max(x)
      xnew<-x*xfac
   
      if(sum(y)!=0){
      yfac<-ch*.8/diff(range(y))
      ynew<-ch-((y-min(y))*yfac+ch*.1)
      py<-ch-((y[t]-min(y))*yfac+ch*.1)
      }else{
        ynew<-rep(ch/2,length(y))
        py<-ch/2
      }
      s<-round(c(rbind(xnew,ynew)))
      .Tcl(paste(canvas$ID,"delete all"))
      .Tcl(paste(canvas$ID,"create line" ,paste(s,collapse=" "), '-fill grey -tags "myline"'))
       px<-(t-1)*xfac
      .Tcl(paste(canvas$ID,"create oval" ,
                 paste(px-3,py-3,px+3,py+3,
                       collapse=" "), '-outline #00ff00 -tags "mypoint"'))
      
    }
    
  }
  ##################################
  ######## WIDGETS #################
  ##################################
  # spinbox widget for voxel coordinates
  coX<-tkspinbox(f5,textvariable=X,from=1,to=d[1],command=onSpin,increment=1,repeatdelay=10)
  coY<-tkspinbox(f5,textvariable=Y,from=1,to=d[2],command=onSpin,increment=1,repeatdelay=10)
  coZ<-tkspinbox(f5,textvariable=Z,from=1,to=d[3],command=onSpin,increment=1,repeatdelay=10)
  
  # add a time spinbox if there is more than 1 image 
  if(dim(func)[4]>1){
    coT<-tkspinbox(f5,textvariable=time,from=1,to=d[4],command=onSpin,increment=1,repeatdelay=10,wrap=TRUE)
  }
  
  # read only widget to hold world coordinates
  coXw<-tkentry(f5,textvariable=Xw,state="readonly",readonlybackground="white")
  coYw<-tkentry(f5,textvariable=Yw,state="readonly",readonlybackground="white")
  coZw<-tkentry(f5,textvariable=Zw,state="readonly",readonlybackground="white")
  
  # labels for the voxel coordinates
  coXLab<-tklabel(parent = f5,text="X",font=myfont)
  coYLab<-tklabel(parent = f5,text="Y",font=myfont)
  coZLab<-tklabel(parent = f5,text="Z",font=myfont)
  if(dim(func)[4]>1){
    coTLab<-tklabel(parent = f5,text="T",font=myfont)
  }
  
  # voxel and world labels
  voxLab<-tklabel(parent=f5,text="Voxel",font=myfont)
  worldLab<-tklabel(parent=f5,text="World",font=myfont)
  
  #  max/min/intensity widgets and labels
  MAX<-tkentry(f5,textvariable=high)
  MIN<-tkentry(f5,textvariable=low)
  intensity<-tkentry(f5,textvariable=intens,state="readonly",readonlybackground="white")
  maxLab<-tklabel(parent = f5,text="Max",font=myfont)
  minLab<-tklabel(parent = f5,text="Min",font=myfont)
  intensLab<-tklabel(parent=f5,text="Intensity",font=myfont)
  baseLab<-tklabel(parent=f5,text="Base  ",font=myfont)
  
  # movie button
  tmovieBut<-tkcheckbutton(f5,variable=tmovie, command=movieT,text="Movie",font=myfont)
  
  #oLab<-tklabel(parent=f5,text="Overlay")
  #oMAX<-tkentry(f5,textvariable=ohigh)
  #oMIN<-tkentry(f5,textvariable=olow)
  #ointensity<-tkentry(f5,textvariable=ointens,state="readonly",readonlybackground="white")
  
  ##################################
  ######### GEOMETRY ###############
  ##################################
  tkgrid(master)
  tkgrid(img)
  tkgrid(f2,f1)
  tkgrid(f3,f4)
  tkgrid(lx)
  tkgrid(ly)
  tkgrid(lz)
  # canvas <- tkcanvas(f4, relief="raised",background="black",
  #                    width=as.numeric(.Tcl("image width rhmmy")),
  #                    height=as.numeric(.Tcl("image height rhmmz")))
  canvas <- tkcanvas(f4, relief="raised",background="black",
                     width=as.numeric(.Tcl(paste("image width" ,rhy))),
                     height=as.numeric(.Tcl(paste("image width" ,rhz))))
  #paste("rh",paste(zlabel,collapse = ""))
  tkgrid(canvas)
  reslicer(coNew = xyz,zoom = zoom,subsamp = subsamp)
  crossHairs()
  crossHairs()
  if(d[4]>1){
    # tkconfigure(canvas,
    #             height=as.numeric(.Tcl("image height rhmmz")),
    #             width=as.numeric(.Tcl("image width rhmmy"))
    # )
    tkconfigure(canvas,
                height=as.numeric(.Tcl(paste("image width" ,rhz))),
                width=as.numeric(.Tcl(paste("image width" ,rhy)))
    )
  }
  tkgrid(f5)
  boxWidth<-round(as.numeric(tkwinfo("reqwidth",f5))/10)
  boxHeight<-round(as.numeric(tkwinfo("reqheight",f5))/6)
  tkplace(coXLab,relx=0,rely=1/6)
  tkplace(coYLab,relx=0,rely=2/6)
  tkplace(coZLab,relx=0,rely=3/6)
  if(dim(func)[4]>1){
    tkplace(coTLab,relx=0,rely=4/6)
    tkplace(tmovieBut,"in"=coT,x=boxWidth)
  }
  tkplace(coX,"in"=coXLab,x=as.numeric(tkwinfo("reqwidth",coXLab)),width=boxWidth,height=boxHeight)
  tkplace(coY,"in"=coYLab,x=as.numeric(tkwinfo("reqwidth",coXLab)),width=boxWidth,height=boxHeight)
  tkplace(coZ,"in"=coZLab,x=as.numeric(tkwinfo("reqwidth",coXLab)),width=boxWidth,height=boxHeight)
  Sys.sleep(.1)
  tkplace(coXw,"in"=coXLab,x=as.numeric(tkwinfo("width",coX))*1.3,width=boxWidth,height=boxHeight)
  tkplace(coYw,"in"=coYLab,x=as.numeric(tkwinfo("width",coX))*1.3,width=boxWidth,height=boxHeight)
  tkplace(coZw,"in"=coZLab,x=as.numeric(tkwinfo("width",coX))*1.3,width=boxWidth,height=boxHeight)
  if(dim(func)[4]>1){tkplace(coT,"in"=coTLab,x=as.numeric(tkwinfo("reqwidth",coXLab)),width=boxWidth,height=boxHeight)}
  tkplace(voxLab,"in"=coX,x=0,y=-as.numeric(tkwinfo("reqheight",coX))*.5,anchor="w",height=boxHeight)
  tkplace(worldLab,"in"=coXw,x=0,y=-as.numeric(tkwinfo("reqheight",coXw))*.5,anchor="w",height=boxHeight)
  height<-tkwinfo("reqheight",master)
  width<-tkwinfo("reqwidth",master)
  tkplace(intensity,"in"=coXLab,x=as.numeric(width)-2*boxWidth-11,width=boxWidth,height=boxHeight)
  tkplace(MAX,"in"=coYLab,x=as.numeric(width)-2*boxWidth-11,width=boxWidth,height=boxHeight)
  tkplace(MIN,"in"=coZLab,x=as.numeric(width)-2*boxWidth-11,width=boxWidth,height=boxHeight)
  tkplace(intensLab,'in'=intensity,x=-1,rely=.5,anchor="e",height=boxHeight)
  tkplace(maxLab,'in'=MAX,x=-1,rely=.5,anchor="e",height=boxHeight)
  tkplace(minLab,'in'=MIN,x=-1,rely=.5,anchor="e",height=boxHeight)
  tkplace(baseLab,"in"=intensity,x=-4,y=-as.numeric(tkwinfo("reqheight",intensity))*.5,anchor="w",height=boxHeight)
  
  #tkplace(oLab,"in"=ointensity,x=-9,y=-as.numeric(tkwinfo("reqheight",ointensity))*.5,anchor="w",height=boxHeight)
  ###################################
  ########## BINDINGS ###############
  ###################################
  tkbind(lx, "<Button-1>",onLeftClick1)
  tkbind(lx, "<B1-Motion>",onLeftClick1)
  tkbind(ly, "<Button-1>",onLeftClick2)
  tkbind(ly, "<B1-Motion>",onLeftClick2)
  tkbind(lz, "<Button-1>",onLeftClick3)
  tkbind(lz, "<B1-Motion>",onLeftClick3)
  tkbind(top, "<Button-3>",crossHairs)
  tkbind(canvas, "<Button-1>",onLeftClick4)
  tkbind(canvas, "<B1-Motion>",onLeftClick4)
  tkbind(coX, "<Return>",onSpin)
  tkbind(coY, "<Return>",onSpin)
  tkbind(coZ, "<Return>",onSpin)
  if(dim(func)[4]>1){tkbind(coT, "<Return>",onSpin)}
  tkbind(MAX, "<Return>",maxminChange)
  tkbind(MIN, "<Return>",maxminChange)
  #tkbind(oMAX, "<Return>",omaxminChange)
  #tkbind(oMIN, "<Return>",omaxminChange)
  ##################################
  ########## CONFIGURES ############
  ##################################
  
  # get crosshair to appear
  tkconfigure(img,cursor="crosshair")
  
  # set a fixed width/height in pixels
  geom<-paste(width,height,sep="x")
  tkwm.geometry(top,geom)
  
  #do not make window resizable
  tkwm.resizable(top,FALSE,FALSE)
}