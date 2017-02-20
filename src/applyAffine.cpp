#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
Rcpp::NumericVector applyAffine(Rcpp::NumericVector yr,arma::mat aff,arma::ivec outDim){
  //////////////////////////////////////////////////////
  ////// Coefficient Cube Construction /////////////////                   
  //////////////////////////////////////////////////////
  Rcpp::IntegerVector arrayDims = yr.attr("dim");
  int x = arrayDims[0];
  int y = arrayDims[1];
  int z = arrayDims[2];
  
  // need arma cube as NumericVector doesn't allow 3d Access. e.g yr(i,j,k)
  arma::cube temp(yr.begin(),x,y,z,false,false);  
  
  
  /* Coefficient cube has n+4 elements in each dimension 
   * due to zeroth and nth +1 element and zero padding for controlling edge effects
   */
  arma::cube coeffx = arma::zeros(x+4, y+4,z+4);  
  
  
  for(int k = 0;k<(z);k++){ 
    for(int j = 0;j<(y);j++){ 
      for(int i = 0;i<(x);i++){
        /* place the data in the middle of the cube of zeros.
         * This effectively zero pads the data.
         */
        coeffx.at(i+2,j+2,k+2) = temp.at(i,j,k);  
      }
    }
  }
  
  // update lengths to reflect zero padding
  x+=2; 
  y+=2;
  z+=2;
  //////////////////////////////////////////////////////
  ////// Spline Coefficient Calculation/////////////////                   
  //////////////////////////////////////////////////////
  
  /* spline coefficients in each dimension can be represented by a 
   * difference equation applied in one direction followed by the reverse.
   * see Unser, M., Aldroubi, A, & Eden, M. (1993). B-Spline Signal-Processing .1. Theory. Ieee Transactions on Signal Processing. 
   * and Unser, M., Aldroubi, A., & Eden, M. (1991). Fast B-Spline Transforms for Continuous Image Representation and Interpolation. IEEE Transactions on Pattern Analysis and Machine Intelligence, 13(3), 277–285. 
   */
  
  
  /* note the value of gi  converges rapidly to the absolute value 
   * of the pole of the filter which is -2-sqrt(3). 
   */
  Rcpp::NumericVector g1(x);
  Rcpp::NumericVector g2(y);
  Rcpp::NumericVector g3(z);
  
  g1[0]=0;
  g2[0]=0;
  g3[0]=0;
  for(int i=1;i<x;i++){g1[i]=1/(4-g1[i-1]);}
  for(int i=1;i<y;i++){g2[i]=1/(4-g2[i-1]);}
  for(int i=1;i<z;i++){g3[i]=1/(4-g3[i-1]);}
  
  /*The B-spline interpolator is an Infinite Impulse Response Filter.
   *  In  B-spline interpolation the data is  represented by a  convolution(*)
   * of a set of unknown coefficients (c) with B-spline Kernels(B). 
   * 
   * y = c * B
   * 
   * c can be obtained by deconvolution.
   * 
   * the deconvolution can be performed by applying the difference equation defined
   * by the  inverse of the z-transform of B to the data y.
   * The  symmetric cubic B-spline kernel has the z trasform:
   * 
   * z+4+z^-1
   * --------
   *    6
   *
   * 
   * The difference equation is c[i]= (y[i]*6-c[i-1])*pole in the forward direction.
   * However, we initialise c to be y. Therefore the difference equation updates to 
   *  
   *   c[i]= (c[i]*6-c[i-1])*pole 
   * 
   * This prevents the copying of large arrays.
   * Then in the negative direction the difference equation is
   * 
   *   c[i] = c[i]-pole*c[i+1].
   *  
   *  edge coefficients are create using natrual splines.
   *  
   *  c[0]=2c[1]-c[2]
   *  c[n] = 2c[n-1]-c[n-2]
   *  
   * however prior to coefficent calcualtion the coefficient array is zero padded.
   * This allows for accruate values to computed when the signal touches the border.
   */  
  
  
  //Loop forward in X//
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int k = 1;k<(z+1);k++){ 
    for(int j = 1;j<(y+1);j++){ 
      for(int i = 2;i<(x);i++){
        coeffx.at(i,j,k) = (coeffx.at(i,j,k)*6-coeffx.at(i-1,j,k))*g1[i-1];
      }
    }
  }
  
  // Loop Backward in X//
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int k = (z);k>0; k--){ 
    for(int j = (y);j>0; j--){ 
      for(int i = (x-1);i>1;  i--){
        coeffx.at(i,j,k) = coeffx.at(i,j,k)-g1[i-1]*coeffx.at(i+1,j,k);
      }
    }
  }
  arma::span span0 = arma::span(0,0);
  arma::span span1 = arma::span(1,1);
  arma::span span2 = arma::span(2,2);
  
  arma::span spanX = arma::span(0,x+1);
  arma::span spanY = arma::span(0,y+1);
  arma::span spanZ = arma::span(0,z+1);
  
  arma::span spanX1 = arma::span(x+1,x+1);
  arma::span spanX0 = arma::span(x,x);
  arma::span spanXM1 = arma::span(x-1,x-1);
  
  //c edge coefficients(c(0)=2*c(1)-c(2))
  //c edge coefficients(c(n+1)=2*c(n)-c(n-1))
  coeffx(span0,spanY,spanZ) = 2*coeffx(span1,spanY,spanZ)-coeffx(span2,spanY,spanZ);
  coeffx(spanX1,spanY,spanZ) = 2*coeffx(spanX0,spanY,spanZ)-coeffx(spanXM1,spanY,spanZ);
  
  
  //Loop forward in y//
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int k = 0; k<(z+2);k++){ 
    for(int j = 2; j<(y); j++){ 
      for(int i = 0; i<(x+2);   i++){
        coeffx.at(i,j,k) = (coeffx.at(i,j,k)*6-coeffx.at(i,j-1,k))*g2[j-1];
      }
    }
  }
  
  // Loop Backward in y//
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int k = (z+1); k>-1;k--){ 
    for(int j = (y-1); j>1;   j--){ 
      for(int i = (x+1); i>-1;i--){
        coeffx.at(i,j,k) = coeffx.at(i,j,k)-g2[j-1]*coeffx.at(i,j+1,k);
      }
    }
  }
  

  arma::span spanY1 = arma::span(y+1,y+1);
  arma::span spanY0 = arma::span(y,y);
  arma::span spanYM1 = arma::span(y-1,y-1);
  //Y edge coefficients(y(0)=2*y(1)-y(2))
  coeffx(spanX,span0,spanZ) = 2*coeffx(spanX,span1,spanZ)-coeffx(spanX,span2,spanZ);
  coeffx(spanX,spanY1,spanZ) = 2*coeffx(spanX,spanY0,spanZ)-coeffx(spanX,spanYM1,spanZ);
  
  //Loop forward in z//
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int j = 0; j<(y+2);j++){ 
    for(int k = 2; k<(z);k++){ 
      for(int i = 0; i<(x+2);i++){
        coeffx.at(i,j,k) = (coeffx.at(i,j,k)*6-coeffx.at(i,j,k-1))*g3[k-1];
      }
    }
  }
  // Loop Backward in z//
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int j = (y+1); j>-1;j--){ 
    for(int k = (z-1); k>1;k--){ 
      for(int i = (x+1); i>-1;i--){
        coeffx.at(i,j,k) = coeffx.at(i,j,k)-g3[k-1]*coeffx.at(i,j,k+1);
      }
    }
  }
 
  
  arma::span spanZ1 = arma::span(z+1,z+1);
  arma::span spanZ0 = arma::span(z,z);
  arma::span spanZM1 = arma::span(z-1,z-1);
  
  //Z edge coefficients(z(0)=2*z(1)-z(2))
  coeffx(spanX,spanY,span0) = 2*coeffx(spanX,spanY,span1)-coeffx(spanX,spanY,span2);
  coeffx(spanX,spanY,spanZ1) = 2*coeffx(spanX,spanY,spanZ0)-coeffx(spanX,spanY,spanZM1); 
  
  //////////////////////////////////////////////////////
  //////            Spline interpolation           /////                   
  //////////////////////////////////////////////////////
  
  // create the output vector
  Rcpp::Dimension d(outDim[0],outDim[1],outDim[2]);                
  Rcpp::NumericVector out(d); 
  int nout = out.size();
  
  // inversion necessary to apply the transform to the coordinates
  arma::mat44 sa = arma::inv(aff);
  
  /* how many elements in each slice?
   * this is necessary to determine the coordinate 
   * I am at on each iteration as I do not fill the output 
   * usign array indexing out(i,j,k) but with pointers, out[i].
   */
  int outxy = outDim[0]*outDim[1];
  
  // make it run in parallel... for "IMPRESSIVE" speedp              
#ifdef _OPENMP
#pragma omp parallel for
#endif                        
  // loop over requested points and...
  for(int i = 0; i<nout;i++){       
    
    arma::vec ip(4);
    /* this calcualtes the coordinate I'm at
     * by first finding the slice number
     * then the column nad finally the row
     * ip[3] = 1 ensures the coordiantes are homogenous
     */
    
    ip[3] = 1;                       
    ip[2] = i/outxy;                  
    int temp = i % outxy;
    ip[1] = temp/outDim[0];
    ip[0] = (temp % outDim[0]);
    
    /* trasform the coordinate using the inverse affine matrix
     * also the index starts at 1 so 1 is added to each point.
     */
    arma::vec p = sa * ip;
    p+=1;
    
    /* cast to int to achieve faster flooring
     * Negative points are always ignored so 
     * this should work.
     */
    int indx = (int)p[0];             
    int indy = (int)p[1];             
    int indz = (int)p[2]; 
    // Rcpp::Rcout<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
    
    // check range
    if(p[0]<0||p[0]>=(x-1)||p[1]>=(y-1)||p[1]<0||p[2]<0||p[2]>=(z-1)){  
      // set to zero if outside range, I'm interpolating, not extrapolating.
      out[i]=0;                                    
    }else{                                         
      
      /* difference between point and floor(point) is the distance between voxels.
       * This is the linear portion of the polynomial
       */
      double dx = p[0]-indx;                       
      
      // square part of the x polynomial
      double dsqx = dx*dx;                         
      
      // cube part of the polynomial
      double dcux = dsqx*dx;                       
      
      /* multiply square and cubic term by 3. 
       * Therefore I'm not multiplying withing nested for loop.
       * This saves 30 multiplications per iteration.
       */
      double dx3 = dx*3;                           
      double dsqx3 = dsqx*3;                       
      
      // linear, square and cubic terms
      double dy = p[1]-indy;                       
      double dsqy = dy*dy;                         
      double dcuy = dsqy*dy;                       
      
      // y B-spline polynomial at the four points
      double yp[4];                                
      yp[0] = (dcuy-3*dsqy +3*dy-1)/-6.0;          
      yp[1] = (dcuy-2*dsqy)/2.0+2.0/3.0;           
      yp[2] = (dcuy-dsqy-dy+1)/-2.0+2.0/3.0;       
      yp[3] = dcuy/6.0;                            
      
      // linear, square and cubic terms
      double dz = p[2]-indz;                       
      double dsqz = dz*dz;                         
      double dcuz = dsqz*dz;                       
      
      // z B-spline polynomial at the four points
      double zp[4];                                
      zp[0] = (dcuz-3*dsqz +3*dz-1)/-6.0;          
      zp[1] = (dcuz-2*dsqz)/2.0+2.0/3.0;           
      zp[2] = (dcuz-dsqz-dz+1)/-2.0+2.0/3.0;       
      zp[3] = dcuz/6.0;
      
      double zit;
      for(int k = 0; k<4; k++){  
        // set the coordinate amd the polynomial
        int cz = indz+k;         
        zit = zp[k];
        
        for(int j = 0; j<4; j++){                 
          int cy  = indy+j;
          
          // the relevant spline coefficients
          double ci = coeffx.at(indx,cy,cz);     
          double cj = coeffx.at(indx+1,cy,cz);    
          double ck = coeffx.at(indx+2,cy,cz);    
          double cl = coeffx.at(indx+3,cy,cz);    
          
          /* The polynomial coefficeints (a,b,c,d)
           * ax^3 +bx^2+cx+d
           */
          double pi = -ci+3*(cj-ck)+cl;            
          double pj = ci-2*cj+ck;                 
          double pk = -ci+ck;                     
          double pl =  pj+6*cj;                   
          
          //explicitly construct the polynomial and accumalate over all j and k iterations 
          out[i]+=(dcux*pi+dsqx3*pj+dx3*pk+pl)*yp[j]*zit; 
        }
        
      }
    }
    /* divide by six in outer loop. 
     * This prevents unnecessary divisions in inner loop
     * this factor is simply the denominator of the transfer 
     * function of the B-spline kernel
     */
    out[i]/=6;                                     
  }
  return(out);
}

