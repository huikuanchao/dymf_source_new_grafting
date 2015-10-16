#include "globals.h"
void init_surface_bonds();
void allocate( void ) ;
void calc_A(void) ;
void fft_init( void ) ;
void initialize_potential( void ) ;
void initialize_configuration( void ) ;
void charge_grid( void ) ;
void read_input( void ) ;

void initialize() {
  int i,j ;

  idum =   -long( time(0) ) ; // 9 ; //

//  read_input() ;

  if ( phiP + phiHA + phiHB + phiHC> 1.0 )
    die("Invalid volume fractions!\n") ;

  n_surf_bonds = 0;

  
  lagrange_weights = 0 ;
  spline_weights = 1 ;

  mem_use = 0. ;

  M = 1 ;
  V = 1.0 ;
  grid_per_partic = 1 ;
  for ( j=0 ; j<Dim ; j++ ) {
    V *= L[j] ;
    M *= Nx[j] ;
    dx[j] = L[j] / double( Nx[j] ) ;
    Lh[j] = 0.5 * L[j] ;
    grid_per_partic *= ( pmeorder + 1 ) ;
    printf("dx: %lf\n" , dx[j] ) ;
  }
  gvol = V / double( M ) ;

  // This is used for density fields //
  ntypes = Nsp+1 ;
  
  cout<<"ntypes "<<ntypes<<endl;

  Rg = sqrt( double( Nda + Ndb ) / 6.0 ) ;
  Rg3 = Rg * Rg * Rg ;
  Range = Rg *0.2;
  cout<<"Range "<<Range<<endl;
  Range2 = Range *Range;
  rho0 = C * double( Nda + Ndb ) / Rg3 ;
  chiAB = chiAB * Rg3/ double( Nda + Ndb ) ;
  chiAC = chiAC / double( Nda + Ndb ) ;
  chiBC = chiBC / double( Nda + Ndb ) ;
  
  kappa = kappa * Rg3 / double( Nda + Ndb ) ;/// double( Nda + Ndb ) ;
  kappa_p = kappa_p* Rg3 / double( Nda + Ndb ) ; //kappa_p / double( Nda + Ndb ) ;

  cout<<"kN "<<kappa*( Nda + Ndb )<<" chiABN "<<( Nda + Ndb )*chiAB<<endl;
  
 double  rdc_V = V * (L[2]-wall_thick*2)/L[2];

  nD = int( ( 1.0 - phiHA - phiHB - phiP -phiHC) * rho0 * rdc_V / ( Nda + Ndb ) * CG_ratio ) ;
  nA = int( phiHA * rho0 * rdc_V / Nha * CG_ratio ) ;
  nB = int( phiHB * rho0 * rdc_V / Nhb * CG_ratio ) ;
  nC = int( phiHC * rho0 * rdc_V / Nhc * CG_ratio ) ; 

  Vp = rho0 ;

  if ( Dim == 2 )
    Vp *= PI * Rp * Rp ;
  else if ( Dim == 3 )
    Vp *= 4.0 * PI * Rp * Rp * Rp / 3.0 ;



  ////////////////////////////////////////
  // Initialize nanoparticle parameters //
  ////////////////////////////////////////
  Diff[2] *= 1.0 / Vp ;

  p_m = Vp;

  Diff_rot *= 1.0 / Vp / Rp / Rp ;


  nP = nsP = 0.0 ;

  if ( phiP > 0.0 ) {

    if ( sigma > 0.0 ) {
      if ( Dim == 2 )
        ng_per_partic = int( sigma * PI * Rp * 2.0 ) ;
      else if ( Dim == 3 )
        ng_per_partic = int ( sigma * 4.0 * PI * Rp * Rp ) ;
    }
    else
      ng_per_partic = 0 ;

    nP = int( phiP * rho0 * rdc_V /(Vp+ Ng * ng_per_partic * CG_ratio )) ;
    //nP = int( phiP * rho0 * V /(Vp));

    printf("nP = %d particles with %d grafted chains per particle\n" , 
        nP , ng_per_partic ) ;

  }

  //nD = int( (( 1.0 - phiHA - phiHB - phiP ) * rho0 * V - nP*Ng*ng_per_partic )/ ( Nda + Ndb ) * CG_ratio ) ;
  //nA = int( phiHA * rho0 * V / Nha * CG_ratio ) ;
  //nB = int( phiHB * rho0 * V / Nhb * CG_ratio ) ;

 
  nsD = nD * (Nda + Ndb) ;
  nsA = nA * Nha ;
  nsB = nB * Nhb ;
  nsC = nC * Nhc ; 
  nsP = nP * ( 1 + ng_per_partic * ( Ng + 1 ) ) ; // + 1 because of the graft site

  printf("Input rho0: %lf , " , rho0 ) ;
 
  rho0 = ( nD * (Nda + Ndb ) + nA * Nha + nB * Nhb + nC*Nhc+nP * (Vp + ng_per_partic * Ng ) ) / rdc_V / CG_ratio ;
  printf("actual rho0: %lf\n" , rho0 ) ;

  printf("\nnD: %d\nnA: %d\nnB: %d\nnC: %d\nnP: %d\n\n" , nD, nA, nB, nC, nP ) ;

  // Derived quantities //
  nstot = nA * Nha + nB * Nhb + nD * ( Nda + Ndb ) + nC*Nhc 
    + nP * ( 1 + ng_per_partic * ( Ng + 1 ) ) ;
  
   kappa *= rho0;
   kappa_p *= rho0;
   chiAB *= rho0;
   for(i =0 ; i<4 ; i++)
    wall_lamb[i] *= kappa;
 

  step = 0 ;
  num_averages = 0.0 ;
  
  buff_ind = 0 ;
  if ( stress_freq > print_freq )
    stress_freq = print_freq ;

  if ( stress_freq > 0 )
    buff_size = print_freq / stress_freq + 1;
  else
    buff_size = 0 ;

  I = complex<double>( 0.0 , 1.0 ) ;

    printf("Total segments: %d\n" , nstot ) ;
    printf("grid_vol: %lf\n" , gvol ) ;
    printf("Particles per grid point: %lf\n" , double(nstot) / double(M) ) ;

  fft_init() ;
    printf("FFTW-MPI Initialized\n") ; fflush( stdout ) ;

  allocate() ;
 

  for(j=0; j<ntypes;  j++){

	if(Diff[j] > 0 ){
	    verlet_a[j] = (1 - delt/2.0/Diff[j]/(j==2 ? p_m: 1.0))/(1 + delt/2.0/Diff[j]/(j==2 ? p_m: 1.0));
	     verlet_b[j] = 1.0 / (1 + delt/2.0/Diff[j]/(j==2 ? p_m: 1.0));
			         
	}
	else{
		            verlet_a[j] = -1;
			                verlet_b[j] = 0;
	}

	printf("verlet_a: %lf and verlet_b: %lf for  type: %d\n" , verlet_a[j],verlet_b[j], j ) ;fflush(stdout) ;

  }//ntypes

  printf("Memory allocated: %lf MB\n" , mem_use / 1.0E6 ) ; 
  fflush(stdout) ;
  
  initialize_configuration() ;
 
  //calc_A();//calc A and dAdang
   if(graft_surface) init_surface_bonds(); 

  printf("Initial config generated\n") ; fflush(stdout) ;

  charge_grid() ;

  printf("grid charged\n"); fflush(stdout); 
  
  initialize_potential() ;

  printf("potentials initialized, written\n") ; fflush( stdout ) ; 
}


void adj_L(){
   int j,i,k ;
   double tmp_L[Dim],tmp_prss=0;
   
   for(i=0;i<Dim ;i++){
	aver_Ptens[i][i] /= aver_Ptens[0][1];

   }

   if(optm_L == 1){
   for(i=0; i<Dim ;i++){

	if(L_dim != i)
		tmp_prss += aver_Ptens[i][i];
	
         tmp_L[i] = L[i];
   }
   
   tmp_prss /= (Dim-1);

    L[L_dim] += lam_L*(aver_Ptens[L_dim][L_dim] - tmp_prss); 

     cout<<"px "<<aver_Ptens[L_dim][L_dim]<<" pv "<<tmp_prss<<endl; 

    for ( j=0 ; j<Dim ; j++ ) {
    if(j != L_dim)
      L[j] = pow(V/L[L_dim],1.0/double(Dim-1));
    
      dx[j] = L[j] / double( Nx[j] ) ;
      Lh[j] = 0.5 * L[j] ;
    }
  }//optm_L=1
  else if(optm_L ==2 ){
    for(i=0; i<Dim ;i++){

	if( i<Dim -1  )
		tmp_prss += aver_Ptens[i][i];
	
         tmp_L[i] = L[i];
     }
   
     tmp_prss /= (Dim-1);

  /*   L[Dim-1] += lam_L*(aver_Ptens[Dim-1][Dim-1] - tmp_prss); 

     L[1] = pow(V/L[Dim-1]/sqrt(3),0.5);
     L[0] = L[1]*sqrt(3); 

     cout<<"px "<<aver_Ptens[Dim-1][Dim-1]<<" pv "<<tmp_prss<<endl; 

     for ( j=0 ; j<Dim ; j++ ) {
    
       dx[j] = L[j] / double( Nx[j] ) ;
       Lh[j] = 0.5 * L[j] ;
     }


*/

  }//optm_l=2

/*
#pragma omp parallel for
  for(i =0 ;i<nstot ; i++){
    for(int m=0;m<Dim ;m++){

	x[i][m] *= L[m]/tmp_L[m];
    	x_bac[i][m] *= L[m]/tmp_L[m];

    }

  }
*/
   for(i=0;i<Dim ;i++){
	aver_Ptens[i][i]  = 0;///= aver_Ptens[0][1];

   }
   aver_Ptens[0][1] = 0;


  //initialize_potential( );


}



void initialize_potential( ) {

  int i, j;

  double ro[Dim], rc[Dim], dr[Dim], mdr2 , pref , mdr , k2, kv[Dim] ;
  pref = V / ( pow( 2.0 * sqrt(PI) *Range, Dim ) ) ; // Note: the factor of V comes from the FFT

  for ( j=0 ; j<Dim ; j++ )
    ro[j] = 0.0 ;

  for ( i=0 ; i<M ; i++ ) {

    get_r( i , rc ) ;

    mdr2 = pbc_mdr2( ro, rc, dr ) ;
    mdr = sqrt( mdr2 ) ;

    uG[i] = exp( -mdr2 / 4.0/Range2 ) * pref ;
    r_dudr[i] = -mdr2 * exp( -mdr2 / 4.0/Range2 ) ;

    tmp[i] = rho0 / 2.0 * ( 1.0 - erf( ( mdr - Rp ) / Xi ) ) * V;
    gammaP[i] = tmp[i] ;
    
  }

  // Set up the particle-particle potential //
  fftw_fwd( tmp , ktmp ) ;
  for ( i=0 ; i<M ; i++ ) 
    ktmp2[i] = ktmp[i] * ktmp[i] ;
  fftw_back( ktmp2 , uP ) ;

  // Set up particle-polymer potential //
 #pragma omp parallel for private(k2,kv)
  for ( i=0 ; i<M ; i++ ) {
    k2 = get_k( i , kv ) ;
    ktmp[i] *= exp( -Range2*k2 /2.0) ;
  }
  fftw_back( ktmp , uPG ) ;


  for ( j=0 ; j<Dim ; j++ ) {
    field_gradient( uG , grad_uG[j] , j ) ;
    field_gradient( uP , grad_uP[j] , j ) ;
    field_gradient( uPG , grad_uPG[j] , j ) ;
  }


  int j2; 
  for ( j=0 ; j<Dim ; j++ ) {

  #pragma omp parallel for private(j2, mdr2,rc, dr )
	for ( i=0 ; i<M ; i++ ) {
	get_r( i , rc ) ;
	mdr2 = pbc_mdr2( rc , ro , dr ) ;
	  for ( j2=0 ; j2<Dim ; j2++ ){
		vir_func[j][j2][i] = dr[j2] * -grad_uG[j][i] ;
		vir_funcpp[j][j2][i] = dr[j2] * -grad_uP[j][i] ;
		vir_funcpg[j][j2][i] = dr[j2] * -grad_uPG[j][i] ;
	  }
 
	}
   }

  for ( j=0 ; j<Dim ; j++ )
	for ( j2=0 ; j2<Dim ; j2++ ){
		fftw_fwd( vir_func[j][j2], vir_func_hat[j][j2] ) ;
		fftw_fwd( vir_funcpp[j][j2], vir_funcpp_hat[j][j2] ) ;
		fftw_fwd( vir_funcpg[j][j2], vir_funcpg_hat[j][j2] ) ;
	}

//  write_grid_data( "ug.dat" , uG ) ;
//  write_grid_data( "up.dat" , uP ) ;
//  write_grid_data( "upg.dat" , uPG ) ;


  for ( j=0 ; j<Dim ; j++ ) {
    char nm[20] ;
    fftw_fwd( grad_uG[j] , grad_uG_hat[j] ) ;
//    sprintf( nm , "grad_ug_%d.dat" , j ) ;
 //   write_grid_data( nm , grad_uG[j] ) ;

    fftw_fwd( grad_uP[j] , grad_uP_hat[j] ) ;
 //   sprintf( nm , "grad_up_%d.dat" , j ) ;
 //   write_grid_data( nm , grad_uP[j] ) ;

    fftw_fwd( grad_uPG[j] , grad_uPG_hat[j] ) ;
   // sprintf( nm , "grad_upg_%d.dat" , j ) ;
   // write_grid_data( nm , grad_uPG[j] ) ;
  }
}

void allocate( ) {

  int i, j,k;


  for(i = 0 ;i<Dim ; i++)
    for(j = 0; j <Dim ; j++)
        aver_Ptens[i][j] = 0 ;

  Stress_bond_t = ( double*** ) calloc( Dim , sizeof( double**) ) ;
 
  W_tsn = (double***) calloc(nstot, sizeof( double** )); 
  
  
  for(i =0; i< nstot; i++){
  	W_tsn[i] = (double ** )calloc( Dim , sizeof( double* ) );


  	for ( j=0 ; j<Dim ; j++ ) {

    		W_tsn[i][j] = ( double* ) calloc( pmeorder+1 , sizeof( double ) );

 	}

 }

  epslon = (double***) calloc(3, sizeof(double**));

  vir_func = ( double*** ) calloc( Dim , sizeof( double**) ) ;
  vir_func_hat = ( complex<double>*** ) calloc( Dim , sizeof( complex<double>**) ) ;
  vir_funcpp = ( double*** ) calloc( Dim , sizeof( double**) ) ;

  vir_funcpp_hat = ( complex<double>*** ) calloc( Dim , sizeof( complex<double>**) ) ;
  
  vir_funcpg = ( double*** ) calloc( Dim , sizeof( double**) ) ;

  vir_funcpg_hat = ( complex<double>*** ) calloc( Dim , sizeof( complex<double>**) ) ;

  for(i=0; i<3; i++){
	Stress_bond_t[i] = (double**) calloc(Dim, sizeof(double*));
	epslon[i] = ( double** ) calloc( 3 , sizeof( double* ) ) ;

	vir_func[i] = (double**) calloc(Dim, sizeof(double*));

	vir_func_hat[i] = (complex<double>**) calloc(Dim, 
sizeof(complex<double>*));
	vir_funcpp[i] = (double**) calloc(Dim, sizeof(double*));
	vir_funcpp_hat[i] = (complex<double>**) calloc(Dim, sizeof(complex<double>*));
	vir_funcpg[i] = (double**) calloc(Dim, sizeof(double*));
	vir_funcpg_hat[i] = (complex<double>**) calloc(Dim, sizeof(complex<double>*));
	
	for(j=0 ;j<3 ;j++){
		Stress_bond_t[i][j] = ( double* ) calloc(nthreads , sizeof( double ) ) ;
		vir_func[i][j] = ( double* ) calloc( M , sizeof( double ) ) ;

		vir_func_hat[i][j] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
		vir_funcpp[i][j] = ( double* ) calloc( M , sizeof( double ) ) ;
		vir_funcpp_hat[i][j] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;		
		vir_funcpg[i][j] = ( double* ) calloc( M , sizeof( double ) ) ;
		vir_funcpg_hat[i][j] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
		
		epslon[i][j] = ( double* ) calloc( 3 , sizeof( double ) ) ;
   	  	for(k=0 ;k<3 ;k++){
			
			if( (i==k) or (k==j) or (i==j) )
				epslon[i][j][k] =0;
			else if ((k>i) and (k >j))
				epslon[i][j][k] = (j >i ? 1 :-1.0);
			else if( (j>i) and (j>k) )
				epslon[i][j][k] = (i > k ? 1 : -1.0);
			else if(i>k and i>j) 
				epslon[i][j][k] = ( j > k ? -1: 1.0); 
			else{}
		}
       }

  }

   q_noise = ( double** ) calloc( nP , sizeof( double* )) ;
   for(i=0 ;i<nP ;i++){

	q_noise[i] = ( double* ) calloc(4  , sizeof( double )) ;

   }


   gn_bac = ( double** ) calloc( nstot , sizeof( double* ) ) ;
   x_bac = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  
  x = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  xc = ( char** ) calloc( nstot , sizeof( char* ) ) ;
  f = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  for ( i=0 ; i<nstot ; i++ ) {
     x_bac[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
     gn_bac[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
    x[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
    f[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;
    xc[i] = ( char* ) calloc( 8 , sizeof( char ) ) ;
  }

  mem_use += nstot * 2 * 2 * Dim * sizeof( double ) ; 
  mem_use += nstot * 8 * sizeof( char ) ; 
  mem_use += nstot*3*(pmeorder+1)*sizeof( double ) ;

  tp = ( int* ) calloc( nstot , sizeof( int ) );
  grid_inds = ( int** ) calloc( nstot , sizeof( int* ) );
  grid_W = ( double** ) calloc( nstot , sizeof( double* ) ) ;
  for ( i=0 ; i<nstot ; i++ ) {
    grid_inds[i] = ( int* ) calloc( grid_per_partic , sizeof( int ) ) ;
    grid_W[i] = ( double* ) calloc( grid_per_partic , sizeof( double ) ) ;
  }

  mem_use += ( nstot + nstot * grid_per_partic ) * sizeof( int ) ;
  mem_use += nstot * grid_per_partic * sizeof( double ) ;

  sts_buf = ( double*** ) calloc( buff_size , sizeof( double** ) ) ;
  sts_buf_pp = ( double*** ) calloc( buff_size , sizeof( double** ) ) ;
  sts_buf_ng = ( double*** ) calloc( buff_size , sizeof( double** ) ) ;
  for ( i=0 ; i<buff_size ; i++ ) {
    sts_buf[i] = ( double** ) calloc( Dim+1 , sizeof( double* ) ) ;
     sts_buf_pp[i] = ( double** ) calloc( Dim+1 , sizeof( double* ) ) ;
     sts_buf_ng[i] = ( double** ) calloc( Dim+1 , sizeof( double* ) ) ;
    for ( j=0 ; j<Dim+1 ; j++ ){
      sts_buf[i][j] = ( double* ) calloc( Dim+nP , sizeof( double ) ) ;
       sts_buf_pp[i][j] = ( double* ) calloc( Dim+nP , sizeof( double ) ) ;
    
      sts_buf_ng[i][j] = ( double* ) calloc( Dim+nP , sizeof( double ) ) ;
 
    }
  }
  mem_use += Dim * Dim * buff_size * sizeof( double ) ;


  for ( i=0 ; i<Dim ; i++ ) {
    grad_uG_hat[i] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
    grad_uP_hat[i] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
    grad_uPG_hat[i] = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
  }
  
  mem_use += 3 * Dim * M * sizeof( double ) ; 
  
  k_wall = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
  ktmp = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
  tmp_PP = ( double* ) calloc( M , sizeof( double ) ) ;
  tmp_Ng =  (complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
  ktmp2 = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
   ktmp3 = ( complex<double>* ) calloc( M , sizeof( complex<double> ) ) ;
  
  tmp = ( double* ) calloc( M , sizeof( double ) ) ;
  tmp2 = ( double* ) calloc( M , sizeof( double ) ) ;
  tmp3 = ( double* ) calloc( M , sizeof( double ) ) ;
  uG = ( double* ) calloc( M , sizeof( double ) ) ;
  uP = ( double* ) calloc( M , sizeof( double ) ) ;
  uPG = ( double* ) calloc( M , sizeof( double ) ) ;
  r_dudr = ( double* ) calloc( M , sizeof( double ) ) ;
  
  mem_use += 7 * M * sizeof( double ) ; 


  // Array that maintains the equilibrium bond lengths for
  // the grafting sites on the grafted nanoparticles. 
  graft_req = ( double** ) calloc( nP , sizeof( double* ) ) ;
  for ( i=0 ; i<nP ; i++ )
    graft_req[i] = ( double* ) calloc( ng_per_partic , sizeof( double ) ) ;

  grf_bf_x = ( double** ) calloc( ng_per_partic*nP, sizeof( double *) ) ;
  for ( i=0 ; i<ng_per_partic*nP ; i++ )
	grf_bf_x[i] = ( double* ) calloc( Dim , sizeof( double ) ) ;

  euler_ang = ( double** ) calloc( nP , sizeof( double* ) ) ;
   euler_q = ( double** ) calloc( nP , sizeof( double* ) ) ;
    euler_adot = ( double** ) calloc( nP , sizeof( double* ) ) ;
   euler_B = ( double*** ) calloc( nP , sizeof( double** ) ) ;
   euler_A = ( double*** ) calloc( nP , sizeof( double** ) ) ;
    euler_Q = ( double*** ) calloc( nP , sizeof( double** ) ) ;
   dAdphi = ( double*** ) calloc( nP , sizeof( double** ) ) ;
  dAdtheta = ( double*** ) calloc( nP , sizeof( double** ) ) ;
  dAdpsi = ( double*** ) calloc( nP , sizeof( double** ) ) ;
   real_trq = ( double** ) calloc( nP , sizeof( double* ) ) ;           
  trq = ( double** ) calloc( nP , sizeof( double* ) ) ;           
  for ( i=0 ; i<nP ; i++ ){
        real_trq [i] = ( double* ) calloc(4 , sizeof( double) ) ;
        trq [i] = ( double* ) calloc(4 , sizeof( double) ) ;
 	euler_ang[i] = ( double* ) calloc( 3 , sizeof( double)); 
   	euler_q[i] = ( double* ) calloc( 4 , sizeof( double)); 
	euler_B[i] = ( double** ) calloc( 4 , sizeof( double*));
    	
	euler_A[i] = ( double** ) calloc( 3 , sizeof( double*));
    	euler_adot[i] = ( double* ) calloc( 4 , sizeof( double)); 
	euler_Q[i] = ( double** ) calloc( 3 , sizeof( double*));
        dAdphi[i] = ( double** ) calloc( 3 , sizeof( double*));
  	dAdtheta[i] = ( double** ) calloc( 3 , sizeof( double*));
	dAdpsi[i] = ( double** ) calloc( 3 , sizeof( double*));
   	for(j=0;j<3;j++){
		dAdphi[i][j] = ( double* ) calloc( 3 , sizeof( double)); 
		dAdpsi[i][j] = ( double* ) calloc( 3 , sizeof( double)); 
		dAdtheta[i][j] = ( double* ) calloc( 3 , sizeof( double)); 
	        euler_A[i][j] = ( double* ) calloc( 3 , sizeof( double));		
	        euler_Q[i][j] = ( double* ) calloc( 3 , sizeof( double));		
	}
	for(j=0;j<4;j++){

		euler_B[i][j] = ( double* ) calloc( 4 , sizeof( double));
	}
   }
 
 mem_use += nP * (ng_per_partic+ 3*27+3*2)* sizeof( double ) ;



  rhot = ( double* ) calloc( M , sizeof( double ) ) ;
  rhoga = ( double* ) calloc( M , sizeof( double ) ) ;
  rhoha = ( double* ) calloc( M , sizeof( double ) ) ;
  rhohb = ( double* ) calloc( M , sizeof( double ) ) ;
  rhohc = ( double* ) calloc( M , sizeof( double ) ) ;
  rhoda = ( double* ) calloc( M , sizeof( double ) ) ;
  rhodb = ( double* ) calloc( M , sizeof( double ) ) ;
  rhop = ( double* ) calloc( M , sizeof( double ) ) ;
  rhow =  ( double* ) calloc( M , sizeof( double ) ) ;

  gammaP = ( double* ) calloc( M , sizeof( double ) ) ;
  smrhop = ( double* ) calloc( M , sizeof( double ) ) ;

  mem_use += 8 * M * sizeof( double ) ; 
  
  rhoga_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhoha_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhohb_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhohc_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhoda_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhodb_t = ( double** ) calloc( nthreads , sizeof( double* ) ) ;
  rhop_t =  ( double** ) calloc( nthreads , sizeof( double* ) ) ;

  for ( i=0 ; i<nthreads ; i++ ) {
    rhoga_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhoha_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhohb_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhohc_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhoda_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhodb_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
    rhop_t[i] = ( double* ) calloc( M , sizeof( double ) ) ;
  }

  mem_use += 5 * M * nthreads * sizeof( double ) ;

  verlet_a = ( double* ) calloc( ntypes , sizeof( double) ) ;
   verlet_b = ( double* ) calloc( ntypes , sizeof( double) ) ;
  rho = ( double** ) calloc( ntypes , sizeof( double* ) ) ;
  rho_hat = ( complex<double>** ) calloc( ntypes , sizeof( complex<double>* ) ) ;
  avg_rho = ( double** ) calloc( ntypes , sizeof( double* ) ) ;
  avg_sk = ( complex<double>** ) calloc( ntypes , sizeof( complex<double>* ) ) ;
  w = ( double** ) calloc( ntypes , sizeof( double* ) ) ;
  for ( i=0 ; i<ntypes ; i++ ) {
    rho_hat[i] = ( complex<double> * ) calloc( M , sizeof( complex<double>) ) ;
    rho[i] = ( double * ) calloc( M , sizeof( double ) ) ;
    avg_sk[i] = ( complex<double> * ) calloc( M , sizeof( complex<double> ) ) ;
    avg_rho[i] = ( double * ) calloc( M , sizeof( double ) ) ;
    w[i] = ( double* ) calloc( M , sizeof( double ) ) ;
  }

  mem_use += 2 * ntypes * M * sizeof( double ) ;
  mem_use += ntypes * M * sizeof( complex<double> ) ;

  for ( j=0 ; j<Dim ; j++ ) {
    grad_uG[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    grad_uP[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    grad_uPG[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    
    gradwA[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    gradwB[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    gradwP[j] = ( double* ) calloc( M , sizeof( double ) ) ;
    gradwC[j] = ( double* ) calloc( M , sizeof( double ) ) ;
  }

  mem_use += Dim * M * 6 * sizeof( double ) ;

}
