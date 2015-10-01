#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]
arma::mat badData(Rcpp::NumericMatrix X,NumericVector meds,NumericVector mads,double nMads,double t) {
  int n = X.nrow(), k = X.ncol();
  arma::mat x(X.begin(), n, k, true); 
  arma::rowvec up  = meds+(meds/100*t) + (nMads*mads);
  arma::rowvec low  = meds-(meds/100*t) - (nMads*mads);
  arma::mat upMat  = arma::repmat(up,n,1);
  arma::mat lowMat  = arma::repmat(low,n,1);
  x.elem(arma::find(x<lowMat|| x>upMat)).fill(arma::datum::nan);
  return(x);
}

// [[Rcpp::export]]
Rcpp::NumericVector colMedian(Rcpp::NumericMatrix X) {
  int n = X.nrow(), k = X.ncol();
  arma::mat x(X.begin(), n, k, false); 
  arma::mat meds = median(x,0);
  return(NumericVector(meds.begin(),meds.end()));
}

// [[Rcpp::export]]
Rcpp::NumericVector colMad(Rcpp::NumericMatrix X){
  int n = X.nrow(), k = X.ncol();
  arma::mat x(X.begin(), n, k, true); 
  arma::mat meds = median(x,0);
  x.each_row() -= meds;
  arma::mat mads = median(abs(x),0)*1.4826;
  return(NumericVector(mads.begin(),mads.end()));
}

// [[Rcpp::export]]
Rcpp::NumericVector colsd(Rcpp::NumericMatrix X){
  int n = X.nrow(), k = X.ncol();
  arma::mat x(X.begin(), n, k, true); 
  arma::mat sds = stddev(x,0,0);
  return(NumericVector(sds.begin(),sds.end()));
}

// [[Rcpp::export]]
arma::mat convolve1d(arma::mat x,arma::colvec fir,int Nfft,bool subtractMed = true) {
  int N = x.n_rows;
  int L = fir.n_rows;
  int Nconv = N + L - 1;
  int diff =  Nfft - Nconv;
  int useless = (Nfft - (N + diff))/2;                      
  int st = useless;                                               
  int en = Nfft - (useless + diff)-1;            
  
  arma::cx_vec ffKern(Nfft);
  ffKern = arma::fft(fir,Nfft);
  arma::rowvec meds;
  if(subtractMed){
    meds = median(x,0);
    x.each_row() -= meds;}
  
  arma::cx_mat cxMat;
  double k = x.n_cols;
  
  if(k>1){
    if(std::ceil(k/2)>std::floor(k/2)){
      cxMat = arma::cx_mat(x.cols(0,std::ceil(k/2)-1),arma::join_horiz(x.cols(std::ceil(k/2),k-1),arma::zeros(N,1)));
    }else{cxMat = arma::cx_mat(x.cols(0,std::ceil(k/2)-1),x.cols(std::ceil(k/2),k-1));}
    
    cxMat = fft(cxMat,Nfft);
    cxMat.each_col() %= ffKern;
    cxMat = arma::ifft(cxMat);
    arma::mat Im = arma::imag(cxMat.rows(st,en));
    x = arma::join_horiz(arma::real(cxMat.rows(st,en)),Im.cols(0,std::floor(k/2-1)));
    
    if(subtractMed){
      x.each_row() += meds;}
    
    return(x);
  }else{
    cxMat = fft(x,Nfft);
    cxMat %= ffKern;
    x = arma::real(arma::ifft(cxMat));
    
    if(subtractMed){
      x.each_row() += meds;}
    
    return(x.rows(st,en));
    
  }
}

// [[Rcpp::export]]
arma::mat hampel(arma::mat x, int k, double t0) {
  int n = x.n_rows;
  arma::mat y = x;
  for (int i = k; i<(n-k); i++){
    arma::mat wind = x.rows((i-k),(i+k));
    arma::mat x0 = median(wind,0);
    wind.each_row() -= x0;
    arma::mat S0 = median(abs(wind),0);
    arma::mat dev = abs(x.row(i)-x0);
    arma::mat thresh = (S0*t0*1.4826);
    arma::mat j = y.row(i);
    j.elem(find(dev>thresh)) = x0.elem(find(dev>thresh));
    y.row(i)=j;
  }
  return(y); }

// [[Rcpp::export]]
Rcpp::List gmm(Rcpp::NumericVector x,int k,Rcpp::NumericVector imeans=NumericVector(1),Rcpp::NumericVector isd = NumericVector(1),Rcpp::NumericVector ilambda= NumericVector(1),bool print= false,double tol = 1e-8,int maxit = 1000){
  arma::rowvec data(x.begin(),x.size(),false);
  int km  = imeans.size();
  int ks  = isd.size();
  int kl  = ilambda.size();
  
  arma::mat initM(imeans.begin(),1,km,false);
  arma::mat initS(isd.begin(),1,ks,false);
  arma::mat initL(ilambda.begin(),1,kl,false);
  
  arma::gmm_diag model;
  if(k==km && k==ks && k == kl && k>1){
    model.set_params(initM,arma::square(initS),initL);
    model.learn(data, k, arma::maha_dist, arma::keep_existing, 0, 1, 1e-10, print);
  }else{model.learn(data, k, arma::maha_dist, arma::random_subset, 0, 1, 1e-10, print);}
  
  double lastLog = model.avg_log_p(data);
  double diff = lastLog;
  int i = 1;
  
  while(diff<(-tol) && i<maxit){
    ++i;
    if(any(arma::vectorise(model.dcovs)<=1e-10)){
      model.learn(data, k, arma::maha_dist, arma::random_subset, 0, 1, 1e-10, print);
      lastLog = model.avg_log_p(data);
      diff = lastLog;}else{
        model.learn(data, k, arma::maha_dist, arma::keep_existing, 0, 1, 1e-10, print);
        diff = -1*(model.avg_log_p(data) - lastLog);
        lastLog = model.avg_log_p(data);}
  }
  
  return Rcpp::List::create(Named("mu") = model.means, Named("sigma")= sqrt(model.dcovs), Named("lambda") = model.hefts,Named("LL") = model.avg_log_p(data),Named("niter") = i);
  }

// [[Rcpp::export]]
arma::cx_mat fftN(arma::mat X,int N){
  return(arma::fft(X,N));              
}

// [[Rcpp::export]]
arma::mat pseudo(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y,bool residuals = false,bool keepMean = false, bool includeIntercept = true) {
  
  
  Rcpp::IntegerVector xArrayDims = x.attr("dim");
  Rcpp::IntegerVector yArrayDims = y.attr("dim");
  arma::mat X(x.begin(), xArrayDims[0], xArrayDims[1], true); 
  if(includeIntercept==true){
    X.insert_cols(xArrayDims[1],arma::ones(xArrayDims[0]));
  }
  
  arma::mat Pinv  = pinv(X);    
  if(yArrayDims[0]==1){return Pinv;}
  
  if(residuals && keepMean && includeIntercept){
    arma::mat Y(y.begin(), yArrayDims[0], yArrayDims[1], false);
    arma::mat Ident(xArrayDims[1]+1,xArrayDims[1]+1,arma::fill::eye);
    Ident(xArrayDims[1],xArrayDims[1]) = 0;
    arma::mat coef  = Pinv*Y;
    arma::mat resid = Y - (X*Ident)*coef;
    return resid;
  }
  
  if(residuals){
    arma::mat Y(y.begin(), yArrayDims[0], yArrayDims[1], false); 
    arma::mat coef  = Pinv*Y;
    arma::mat resid = Y - X*coef;
    return resid;
  }
  
  if(yArrayDims[0]==xArrayDims[0]){
    arma::mat Y(y.begin(), yArrayDims[0], yArrayDims[1], false); 
    arma::mat coef  = Pinv*Y;
    return coef;
  }
  return Pinv;
}

// [[Rcpp::export]]
Rcpp::NumericVector rowMad(Rcpp::NumericMatrix X) {
  int n = X.nrow(), k = X.ncol();
  arma::mat x(X.begin(), n, k, true); 
  arma::mat meds = median(x,1);
  x.each_col() -= meds;
  arma::mat mads = median(abs(x),1)*1.4826;
  return(NumericVector(mads.begin(),mads.end()));
}

// [[Rcpp::export]]
Rcpp::NumericVector rowMedian(Rcpp::NumericMatrix X) {
  int n = X.nrow(), k = X.ncol();
  arma::mat x(X.begin(), n, k, false); 
  arma::mat meds = median(x,1);
  return(NumericVector(meds.begin(),meds.end()));
}

// [[Rcpp::export]]
Rcpp::NumericVector rowsd(Rcpp::NumericMatrix X){
  int n = X.nrow(), k = X.ncol();
  arma::mat x(X.begin(), n, k, true); 
  arma::mat sds = stddev(x,0,1);
  return(NumericVector(sds.begin(),sds.end()));
}

// [[Rcpp::export]]
arma::cube sepConvolve3d(NumericVector x,arma::colvec kernX,arma::colvec kernY,arma::colvec kernZ,int Nx,int Ny,int Nz) {
  
  NumericVector arrayDims = x.attr("dim");
  arma::cube cubeArray(x.begin(), arrayDims[0], arrayDims[1], arrayDims[2], true);
  arma::mat N(arrayDims.begin(),1,3,false);
  arma::mat L;
  L << kernX.n_rows << kernY.n_rows << kernZ.n_rows << arma::endr;
  arma::mat Nconv = N + L - 1;
  arma::mat Nffts;
  Nffts << Nx << Ny << Nz << arma::endr;
  
  arma::mat diff =  Nffts - Nconv;
  arma::mat useless = (Nffts - (N + diff))/2;                      
  arma::mat st = useless;                                               
  arma::mat en = Nffts - (useless + diff)-1;            
  
  
  arma::cx_vec ffKernX(Nx);
  ffKernX = arma::fft(kernX,Nx);
  arma::colvec tempMaskX(arrayDims[0]);
  tempMaskX.ones();
  arma::cx_colvec cFacX = arma::ifft(arma::fft(tempMaskX,Nx)%ffKernX);
  
  arma::cx_vec ffKernY(Ny);
  ffKernY = arma::fft(kernY,Ny);
  arma::colvec tempMaskY(arrayDims[1]);
  tempMaskY.ones();
  arma::cx_colvec cFacY = arma::ifft(arma::fft(tempMaskY,Ny)%ffKernY);
  
  arma::cx_vec ffKernZ(Nz);
  ffKernZ = arma::fft(kernZ,Nz);
  arma::vec tempMaskZ(arrayDims[2]);
  tempMaskZ.ones();
  arma::cx_colvec cFacZ = arma::ifft(arma::fft(tempMaskZ,Nz)%ffKernZ);
  
  
  
  
  arma::cx_mat cxSlice;
  double k = arrayDims[1];
  
  for (int i=0; i < arrayDims[2]; i++){ 
    if(std::ceil(k/2)>std::floor(k/2)){
      cxSlice = arma::cx_mat(cubeArray.slice(i).cols(0,std::ceil(k/2)-1),arma::join_horiz(cubeArray.slice(i).cols(std::ceil(k/2),k-1),arma::zeros(arrayDims[0],1)));
    }else{cxSlice = arma::cx_mat(cubeArray.slice(i).cols(0,std::ceil(k/2)-1),cubeArray.slice(i).cols(std::ceil(k/2),k-1));}
    cxSlice = fft(cxSlice,Nx);
    cxSlice.each_col() %= ffKernX;
    cxSlice = arma::ifft(cxSlice);
    cxSlice.each_col() /=cFacX;
    arma::mat Im = arma::imag(cxSlice.rows(st.at(0,0),en.at(0,0)));
    cubeArray.slice(i) = arma::join_horiz(arma::real(cxSlice.rows(st.at(0,0),en.at(0,0))),Im.cols(0,std::floor(k/2-1)));
  }
  
  
  arma::mat sliceTempY(arrayDims[1],arrayDims[0]);
  double ky = arrayDims[0];
  
  
  for (int i=0; i < arrayDims[2]; i++){ 
    sliceTempY =  trans(cubeArray.slice(i));
    if(std::ceil(ky/2)>std::floor(ky/2)){
      cxSlice = arma::cx_mat(sliceTempY.cols(0,std::ceil(ky/2)-1),arma::join_horiz(sliceTempY.cols(std::ceil(ky/2),ky-1),arma::zeros(arrayDims[1],1)));
    }else{cxSlice = arma::cx_mat(sliceTempY.cols(0,std::ceil(ky/2)-1),sliceTempY.cols(std::ceil(ky/2),ky-1));}
    cxSlice = fft(cxSlice,Ny);
    cxSlice.each_col() %= ffKernY;
    cxSlice = arma::ifft(cxSlice);
    cxSlice.each_col() /=cFacY;
    arma::mat Im = arma::imag(cxSlice.rows(st.at(0,1),en.at(0,1)));
    cubeArray.slice(i) = arma::trans(arma::join_horiz(arma::real(cxSlice.rows(st.at(0,1),en.at(0,1))),Im.cols(0,std::floor(ky/2-1))));
  }
  
  
  double kz = arrayDims[1];
  
  for (int i=0; i < arrayDims[0]; i++){ 
    arma::mat sliceInit = cubeArray.subcube(arma::span(i),arma::span(),arma::span());
    arma::mat sliceTemp =  trans(sliceInit);
    if(std::ceil(kz/2)>std::floor(kz/2)){
      cxSlice = arma::cx_mat(sliceTemp.cols(0,std::ceil(kz/2)-1),arma::join_horiz(sliceTemp.cols(std::ceil(kz/2),kz-1),arma::zeros(arrayDims[2],1)));
    }else{cxSlice = arma::cx_mat(sliceTemp.cols(0,std::ceil(kz/2)-1),sliceTemp.cols(std::ceil(kz/2),kz-1));}
    cxSlice = fft(cxSlice,Nz);
    cxSlice.each_col() %= ffKernZ;
    cxSlice = arma::ifft(cxSlice);
    cxSlice.each_col() /=cFacZ;
    arma::mat Im = arma::imag(cxSlice.rows(st.at(0,2),en.at(0,2)));
    cubeArray(arma::span(i),arma::span(),arma::span()) = arma::trans(arma::join_horiz(arma::real(cxSlice.rows(st.at(0,2),en.at(0,2))),Im.cols(0,std::floor(kz/2-1))));
  }
  
  
  return(cubeArray);
  }

// [[Rcpp::export]]
NumericVector zero_na(NumericVector input) {
  arma::vec arrayDims = input.attr("dim"); 
  int ndims = arrayDims.n_elem;
  
  if(ndims ==2 ){
    arrayDims.resize(3);
    arrayDims(2) = 1;
  }
  
  int x=arrayDims(0);
  int y=arrayDims(1);
  int z=arrayDims(2);
  ///////////////// DECLARE CUBE /////////////
  arma::cube X(input.begin(),x,y,z);
  
  X.elem( arma::find_nonfinite(X) ).zeros();
  ////////// WRAP IT ALL UP INTO VECTOR ///////////
  Rcpp::NumericVector out = wrap(X);
  if(ndims==2){out.attr("dim")= IntegerVector::create(arrayDims(0),arrayDims(1));}
  return out;
  /////////////// END //////////////////////////////
}

// [[Rcpp::export]]
Rcpp::NumericVector dilate(Rcpp::NumericVector input,int k) {
  ////// GET DIMENSIONS ////////////////////
  arma::vec arrayDims = input.attr("dim"); 
  int ndims = arrayDims.n_elem;
  
  if(ndims ==2 ){
    arrayDims.resize(3);
    arrayDims(2) = 1;
  }
  int x=arrayDims(0);
  int y=arrayDims(1);
  int z=arrayDims(2);
  ///////////////// DECLARE CUBE /////////////
  arma::cube X(input.begin(),x,y,z);
  arma::cube outCu = X;
  ////////////// LOOP ON SLICES //////////////
  
  for(int j = 0;j<z;j++){
    
    arma::mat sli = X.slice(j);
    arma::mat tempSli = sli;
    arma::mat outSli = sli;
    
    int begin;
    int end;
    //// loop on columns/////////
    
    for (int i = 0; i<x ; i++){
      begin = i-k;
      end=i+k;
      if(begin<=0){begin = 0;}
      if(end>=(x-1)){end = x-1;}
      
      arma::mat wind = sli.rows((begin),(end));
      arma::mat check =max(wind,0);       
      tempSli.row(i) =check;          
    }
    //// loop on rows ///////////
    for (int i = 0; i<y; i++){
      begin = i-k;
      end=i+k;
      if(begin<=0){begin = 0;}
      if(end>=(y-1)){end = y-1;}
      
      arma::mat wind = tempSli.cols((begin),(end));
      arma::mat check =max(wind,1);       
      outSli.col(i) = check;          
    }
    
    ////// update the output /////
    outCu.slice(j) = outSli;
    
  }
  
  for (int j=0; j < x; j++){ 
    arma::mat sliceInit = outCu.subcube(arma::span(j),arma::span(),arma::span());
    arma::mat sliceTemp =  trans(sliceInit);
    arma::mat sliceOut = sliceTemp;        
    int begin;
    int end;
    //// loop on columns/////////
    
    for (int i = 0; i<z ; i++){
      begin = i-k;
      end=i+k;
      if(begin<=0){begin = 0;}
      if(end>=(z-1)){end = z-1;}
      
      arma::mat wind = sliceTemp.rows((begin),(end));
      arma::mat check =max(wind,0);       
      sliceOut.row(i) =check;          
    }
    outCu(arma::span(j),arma::span(),arma::span()) = arma::trans(sliceOut);
  }
  ////////// WRAP IT ALL UP INTO VECTOR ///////////
  Rcpp::NumericVector out = Rcpp::wrap(outCu);
  if(ndims==2){out.attr("dim")= IntegerVector::create(arrayDims(0),arrayDims(1));}
  return out;
  /////////////// END //////////////////////////////
}

// [[Rcpp::export]]
Rcpp::NumericVector erode(Rcpp::NumericVector input,int k) {
  ////// GET DIMENSIONS ////////////////////
  arma::vec arrayDims = input.attr("dim"); 
  int ndims = arrayDims.n_elem;
  
  if(ndims ==2 ){
    arrayDims.resize(3);
    arrayDims(2) = 1;
  }
  int x=arrayDims(0);
  int y=arrayDims(1);
  int z=arrayDims(2);
  ///////////////// DECLARE CUBE /////////////
  arma::cube X(input.begin(),x,y,z);
  arma::cube outCu = X;
  ////////////// LOOP ON SLICES //////////////
  
  for(int j = 0;j<z;j++){
    
    arma::mat sli = X.slice(j);
    arma::mat tempSli = sli;
    arma::mat outSli = sli;
    
    int begin;
    int end;
    //// loop on columns/////////
    
    for (int i = 0; i<x ; i++){
      begin = i-k;
      end=i+k;
      if(begin<=0){begin = 0;}
      if(end>=(x-1)){end = x-1;}
      
      arma::mat wind = sli.rows((begin),(end));
      arma::mat check =min(wind,0);       
      tempSli.row(i) =check;          
    }
    //// loop on rows ///////////
    for (int i = 0; i<y; i++){
      begin = i-k;
      end=i+k;
      if(begin<=0){begin = 0;}
      if(end>=(y-1)){end = y-1;}
      
      arma::mat wind = tempSli.cols((begin),(end));
      arma::mat check =min(wind,1);       
      outSli.col(i) = check;          
    }
    
    ////// update the output /////
    outCu.slice(j) = outSli;
    
  }
  
  for (int j=0; j < x; j++){ 
    arma::mat sliceInit = outCu.subcube(arma::span(j),arma::span(),arma::span());
    arma::mat sliceTemp =  trans(sliceInit);
    arma::mat sliceOut = sliceTemp;        
    int begin;
    int end;
    //// loop on columns/////////
    
    for (int i = 0; i<z ; i++){
      begin = i-k;
      end=i+k;
      if(begin<=0){begin = 0;}
      if(end>=(z-1)){end = z-1;}
      
      arma::mat wind = sliceTemp.rows((begin),(end));
      arma::mat check =min(wind,0);       
      sliceOut.row(i) =check;          
    }
    outCu(arma::span(j),arma::span(),arma::span()) = arma::trans(sliceOut);
  }
  ////////// WRAP IT ALL UP INTO VECTOR ///////////
  Rcpp::NumericVector out = Rcpp::wrap(outCu);
  if(ndims==2){out.attr("dim")= IntegerVector::create(arrayDims(0),arrayDims(1));}
  return out;
  /////////////// END //////////////////////////////
}

// [[Rcpp::export]]
Rcpp::NumericVector dcombine(Rcpp::NumericVector X,Rcpp::IntegerVector dim){
  Rcpp::NumericVector out = X;
  out.attr("dim")= IntegerVector::create(dim(0),dim(1),dim(2),dim(3));
  return(out);
}

// [[Rcpp::export]]
IntegerVector icombine(Rcpp::IntegerVector X,Rcpp::IntegerVector dim){
  Rcpp::IntegerVector out = X;
  out.attr("dim")= IntegerVector::create(dim(0),dim(1),dim(2),dim(3));
  return(out);
}
