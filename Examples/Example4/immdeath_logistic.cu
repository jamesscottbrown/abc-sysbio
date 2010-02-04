
#define NXMAX 1
__constant__ int nx[] = {1, 1};
__constant__ float smatrix[] = {1, -1};
    
__device__ void hazards(float& h0, float* x, int& model, MersenneTwisterState &state, unsigned int threadID){

  int nh=2;
  float h[] = {0,0};

  // Immigration - death
  if(model==0){
    
    // Access the parameters via a 2D texture
    float im = tex2D(param_tex, 1, threadID);
    float d =  tex2D(param_tex, 2, threadID);

    // Fill the hazard array
    h[0] = im;
    if(x[0] > 0){
      h[1] = d*x[0];
    } 

    // Sample the reactions given current hazards
    h0 = h[0] + h[1];
    int i = sample(nh,h0,h,state,threadID);

    // Update the stochiometery
    x[0] = x[0] + smatrix[i];
  }
  
  // Logistic growth
  if(model==1){
    if( x[0] > 0 ){ 
      
      // Access the parameters via a 2D texture
      float b1 = tex2D(param_tex, 1, threadID);
      float d1 = tex2D(param_tex, 2, threadID);
      float d2 = tex2D(param_tex, 3, threadID);

      // Fill the hazard array
      h[0] = b1*x[0];
      h[1] = (d1 + d2*x[0])*x[0];
         
      // Sample the reactions given current hazards
      h0 = h[0] + h[1];
      int i = sample(nh,h0,h,state,threadID);

      // Update the stochiometery
      x[0] = x[0] + smatrix[i];
    }else{
      h0 = 0;
    }
  }

}
