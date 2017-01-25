#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
Rcpp::NumericVector applyAffine(Rcpp::NumericVector yr,arma::mat aff){
//arma::mat applyAffine(Rcpp::NumericVector yr,arma::mat aff){
    
  //////////////////////////////////////////////////////
  ////// Coefficient Cube Construction /////////////////                   
  //////////////////////////////////////////////////////
  Rcpp::NumericVector arrayDims = yr.attr("dim");
  int x = arrayDims[0];
  int y = arrayDims[1];
  int z = arrayDims[2];
  int xy = x*y;
  arma::cube temp(yr.begin(),x,y,z,false,false);  // need arma cube as NumericVector doesn't allow 3d Access. e.g yr(i,j,k)
  arma::cube coeffx = arma::zeros(x+2, y+2,z+2);  // Coefficient cube has n+2 elements in each dimension due to zeroth and nth +1 element
  
  for(int k = 0;k<(z);k++){ 
    for(int j = 0;j<(y);j++){ 
      for(int i = 0;i<(x);i++){
        coeffx.at(i+1,j+1,k+1) = temp.at(i,j,k);  // fill cube with zero padding for edge coefficients
      }
    }
  }
  //////////////////////////////////////////////////////
  ////// Spline Coefficient Calculation/////////////////                   
  //////////////////////////////////////////////////////
  
  // spline coefficients in each dimension can be represented by a difference equation applied in one direction followed by the reverse.
  Rcpp::NumericVector g1(x);
  Rcpp::NumericVector g2(y);
  Rcpp::NumericVector g3(z);
  
  g1[0]=0;
  g2[0]=0;
  g3[0]=0;
  for(int i=1;i<x;i++){g1[i]=1/(4-g1[i-1]);}
  for(int i=1;i<y;i++){g2[i]=1/(4-g2[i-1]);}
  for(int i=1;i<z;i++){g3[i]=1/(4-g3[i-1]);}
  
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
  
  //x edge coefficients(x(0)=2*x(1)-x(2))
  coeffx(arma::span(0,0),arma::span(0,y+1),arma::span(0,z+1)) 
    = 2*coeffx(arma::span(1,1),arma::span(0,y+1),arma::span(0,z+1))
    -coeffx(arma::span(2,2),arma::span(0,y+1),arma::span(0,z+1));
    coeffx(arma::span(x+1,x+1),arma::span(0,y+1),arma::span(0,z+1)) 
      = 2*coeffx(arma::span(x,x),arma::span(0,y+1),arma::span(0,z+1))
      -coeffx(arma::span(x-1,x-1),arma::span(0,y+1),arma::span(0,z+1));
      
      //Loop forward in y//
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int k = 0; k<(z+2);   k++){ 
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
      //Y edge coefficients(x(0)=2*x(1)-x(2))
      coeffx(arma::span(0,x+1),arma::span(0,0),arma::span(0,z+1)) 
        = 2*coeffx(arma::span(0,x+1),arma::span(1,1),arma::span(0,z+1))
        -coeffx(arma::span(0,x+1),arma::span(2,2),arma::span(0,z+1));
        coeffx(arma::span(0,x+1),arma::span(y+1,y+1),arma::span(0,z+1)) 
          = 2*coeffx(arma::span(0,x+1),arma::span(y,y),arma::span(0,z+1))
          -coeffx(arma::span(0,x+1),arma::span(y-1,y-1),arma::span(0,z+1));
          
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
          
          //Z edge coefficients(x(0)=2*x(1)-x(2))
          coeffx(arma::span(0,x+1),arma::span(0,y+1),arma::span(0,0)) 
            = 2*coeffx(arma::span(0,x+1),arma::span(0,y+1),arma::span(1,1))
            -coeffx(arma::span(0,x+1),arma::span(0,y+1),arma::span(2,2));
            coeffx(arma::span(0,x+1),arma::span(0,y+1),arma::span(z+1,z+1)) 
              = 2*coeffx(arma::span(0,x+1),arma::span(0,y+1),arma::span(z,z))
              -coeffx(arma::span(0,x+1),arma::span(0,y+1),arma::span(z-1,z-1)); 
              //////////////////////////////////////////////////////
              //////            Spline interpolation           /////                   
              //////////////////////////////////////////////////////
              int nout = temp.n_elem;
              Rcpp::Dimension d(arrayDims);                // get the dim object
              Rcpp::NumericVector out(d); 
            //  std::fill(out.begin(), out.end(), 0);
                
              arma::mat sa = arma::inv(aff);
              
#ifdef _OPENMP
#pragma omp parallel for
#endif           // make it run in parallel... for "IMPRESSIVE" speedp
              for(int i = 0; i<nout;i++){       // loop over requested points and...
                arma::vec ip(4);
                ip[3] = 1;
                ip[2] = i/xy+1;
                int temp = i % xy;
                ip[1] = temp/x+1;
                ip[0] = (temp % x)+1;
                

                arma::vec p = sa * ip;
                //std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<"\n";  
                int indx = (int)p[0];             // cast to int to achieve faster flooring (hopefully int is always 32 bit)
                int indy = (int)p[1];             // cast to int to achieve faster flooring (hopefully int is always 32 bit)
                int indz = (int)p[2];             // cast to int to achieve faster flooring (hopefully int is always 32 bit)
                //std::cout<<"not Setting to zero";
               if(indx<1||p[0]>x||p[1]>y||indy<1||indz<1||p[2]>z){  // check range
                 out[i]=0;                                    // set to zero if outside range, I'm interpolating, not extrapolating.
                 //std::cout<<"setting to zero";
               }else{                                         // otherwise...
                 
                // std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<"\n";  
                 
                  double dx = p[0]-indx;                         // difference between point and floor(point) forms the spline polynomial
                  //std::cout<<dx<<"\n";
                  double dsqx = dx*dx;                         // square part of the x polynomial
                  double dcux = dsqx*dx;                       // cube part of the polynomial
                  double dx3 = dx*3;                           // multiply by 3 so not multiplying withing nested for loop(saves 15 multiplications per iteration)
                  double dsqx3 = dsqx*3;                       // multiply by 3 so not multiplying withing nested for loop(saves 15 multiplications per iteration)
                  
                  double dy = p[1]-indy;                         //...same as x
                  double dsqy = dy*dy;                         //...same as x
                  double dcuy = dsqy*dy;                       //...same as x
                  
                  double yp[4];                                // holds  y B-spline polynomial
                  yp[0] = (dcuy-3*dsqy +3*dy-1)/-6.0;          // at index 1
                  yp[1] = (dcuy-2*dsqy)/2.0+2.0/3.0;           // at index 2
                  yp[2] = (dcuy-dsqy-dy+1)/-2.0+2.0/3.0;       // at index 3
                  yp[3] = dcuy/6.0;                            // at index 4... zero outside this range so therefore the B-spline is shifted(that's why points are floored)
                  
                  
                  double dz = p[2]-indz;                         //...same as y
                  double dsqz = dz*dz;                         //...same as y
                  double dcuz = dsqz*dz;                       //...same as y
                  
                  double zp[4];                                //...same as y
                  zp[0] = (dcuz-3*dsqz +3*dz-1)/-6.0;          //...same as y
                  zp[1] = (dcuz-2*dsqz)/2.0+2.0/3.0;           //...same as y
                  zp[2] = (dcuz-dsqz-dz+1)/-2.0+2.0/3.0;       //...same as y
                  zp[3] = dcuz/6.0;
                  
                  
                  for(int k = -1; k<3; k++){                   // loop over z polynomial
                    int cz = indz+k;                           // find the coordinate for the appropriate coefficients
                    if(cz==(z+2)){cz=z+1;}                     // check if coordinate greater than 102 and then set to 101
                    double zit = zp[k+1];                      // pick the appropriate polynomial
                    
                    for(int j = -1; j<3; j++){                 // loop over y polynomial
                      int cy  = indy+j;
                      double ci = coeffx.at(indx-1,cy,cz);     // find the first coefficient 
                      double cj = coeffx.at(indx,cy,cz);       // find the second coefficient
                      double ck = coeffx.at(indx+1,cy,cz);     // find the third coefficient
                      double cl = coeffx.at(indx+2,cy,cz);     // find the fourth coefficient
                      double pi = -ci+3*(cj-ck)+cl;            // Construct the coefficiets of the polynomial in X
                      // pj = 3*ci-6*cj+3*ck;                  // Construct the coefficiets of the polynomial in X
                      // pk = -3*ci+3*ck;                      // Construct the coefficiets of the polynomial in X
                      double pj = ci-2*cj+ck;                  // Construct the coefficiets of the polynomial in X
                      double pk = -ci+ck;                      // Construct the coefficiets of the polynomial in X
                      double pl =  pj+6*cj;                    // Construct the coefficiets of the polynomial in X
                      //outx=(dcux*pi+dsqx3*pj+dx3*pk+pl)*yp[j+1]*zit;
                      
                      out[i]+=(dcux*pi+dsqx3*pj+dx3*pk+pl)*yp[j+1]*zit; //explicitly construct the polynomial and accumalte over all j and k iterations 
                      
                    }
                    
                  }
               }
               
                out[i]/=6;                                     // divide by six in outer loop to prevent unnecessaey divisions
              }
              //SEXP p = Rcpp::wrap(out);
              return(out);                                     // Return NumericVector to R 
}

