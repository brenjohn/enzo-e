// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPmUpdate.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPmUpdate class


#include "charm_simulation.hpp"
#include "enzo.hpp"

// #define DEBUG_UPDATE

#ifdef DEBUG_UPDATE
#  define TRACE_PM(MESSAGE)						\
  CkPrintf ("%s:%d %s\n",						\
	    __FILE__,__LINE__,MESSAGE);				
#else
#  define TRACE_PM(MESSAGE) /* ... */
#endif

//----------------------------------------------------------------------

EnzoMethodPmUpdate::EnzoMethodPmUpdate 
( double max_dt ) 
  : Method(),
    max_dt_(max_dt)
{
  TRACE_PM("EnzoMethodPmUpdate()");
  // Initialize default Refresh object

  Refresh & refresh = new_refresh(ir_post_);
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  
  const int rank = cello::rank();
  
  if (rank >= 1) refresh.add_field("acceleration_x");
  if (rank >= 2) refresh.add_field("acceleration_y");
  if (rank >= 3) refresh.add_field("acceleration_z");

  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();

  int num_mass = particle_groups->size("has_mass");

  for (int ipt = 0; ipt < num_mass; ipt++)
    refresh.add_particle(particle_descr->type_index(
                                particle_groups->item("has_mass",ipt)
                                                        ));

  // PM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPmUpdate::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | max_dt_;
}

//----------------------------------------------------------------------

void EnzoMethodPmUpdate::compute ( Block * block) throw()
{
  TRACE_PM("compute()");

  if (block->is_leaf()) {

#ifdef DEBUG_UPDATE    
    double a3sum[3]={0.0};
    double v3sum[3]={0.0};
    double a3sum2[3]={0.0};
    double v3sum2[3]={0.0};
#endif    
    
    EnzoPhysicsCosmology * cosmology = enzo::cosmology();

    enzo_float cosmo_a=1.0,cosmo_dadt=0.0;
    
    if (cosmology) {
      double time = block->time();
      double dt   = block->dt();
      cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
      //      cosmology-> compute_expansion_factor (&av,&dadtv,time+dt);
    }

    const int rank = cello::rank();

    const double dt = block->dt();

    double dt_shift = 0.5*dt/cosmo_a;

    const double cp = dt/cosmo_a;
    const double coef = 0.25*cosmo_dadt/cosmo_a*dt;
    const double cvv = (1.0 - coef) / (1.0 + coef);
    const double cva = 0.5*dt / (1.0 + coef);


    Particle particle = block->data()->particle();
    ParticleDescr * particle_descr = cello::particle_descr();
    Grouping * particle_groups     = particle_descr->groups();

    int num_mass = particle_groups->size("has_mass");

    for (int ipt = 0; ipt < num_mass; ipt++){

      std::string particle_type = particle_groups->item("has_mass",ipt);
      int it = particle.type_index (particle_type);

      //    double dt_shift = 0.0;
      if (rank >= 1) {
        EnzoComputeCicInterp interp_x ("acceleration_x", particle_type, "ax", dt_shift);
        interp_x.compute(block);
      }

      if (rank >= 2) {
        EnzoComputeCicInterp interp_y ("acceleration_y", particle_type, "ay", dt_shift);
        interp_y.compute(block);
      }

      if (rank >= 3) {
        EnzoComputeCicInterp interp_z ("acceleration_z", particle_type, "az", dt_shift);
        interp_z.compute(block);
      }


      int ia_x  = (rank >= 1) ? particle.attribute_index (it, "x") : -1;
      int ia_y  = (rank >= 2) ? particle.attribute_index (it, "y") : -1;
      int ia_z  = (rank >= 3) ? particle.attribute_index (it, "z") : -1;

      int ia_vx = (rank >= 1) ? particle.attribute_index (it, "vx") : -1;
      int ia_vy = (rank >= 2) ? particle.attribute_index (it, "vy") : -1;
      int ia_vz = (rank >= 3) ? particle.attribute_index (it, "vz") : -1;

      int ia_ax = (rank >= 1) ? particle.attribute_index (it, "ax") : -1;
      int ia_ay = (rank >= 2) ? particle.attribute_index (it, "ay") : -1;
      int ia_az = (rank >= 3) ? particle.attribute_index (it, "az") : -1;

      int dp = particle.stride(it, ia_x);
      int dv = particle.stride(it, ia_vx);
      int da = particle.stride(it, ia_ax);

      int nb = particle.num_batches (it);

      // check precisions match

      int ba = particle.attribute_bytes(it,ia_x); // "bytes (actual)"
      int be = sizeof(enzo_float);                // "bytes (expected)"

      ASSERT4 ("EnzoMethodPmUpdate::compute()",
	       "Particle type %s attribute %s defined as %s but expecting %s",
	       particle.type_name(it).c_str(),
	       particle.attribute_name(it,ia_x).c_str(),
	       ((ba == 4) ? "single" :
	        ((ba == 8) ? "double" : "quadruple")),
	       ((be == 4) ? "single" :
	        ((be == 8) ? "double" : "quadruple")),
	       (ba == be));

      for (int ib=0; ib<nb; ib++) {

        enzo_float *x=0, *y=0, *z=0;
        enzo_float *vx=0, *vy=0, *vz=0;
        enzo_float *ax=0, *ay=0, *az=0;

        if (rank >= 1) {
	  x  = (enzo_float *) particle.attribute_array (it, ia_x,  ib);
	  vx = (enzo_float *) particle.attribute_array (it, ia_vx, ib);
	  ax = (enzo_float *) particle.attribute_array (it, ia_ax, ib);
        }
        if (rank >= 2) {
	  y  = (enzo_float *) particle.attribute_array (it, ia_y,  ib);
	  vy = (enzo_float *) particle.attribute_array (it, ia_vy, ib);
	  ay = (enzo_float *) particle.attribute_array (it, ia_ay, ib);
        }
        if (rank >= 3) {
 	  z  = (enzo_float *) particle.attribute_array (it, ia_z,  ib);
	  vz = (enzo_float *) particle.attribute_array (it, ia_vz, ib);
	  az = (enzo_float *) particle.attribute_array (it, ia_az, ib);
        }

        int np = particle.num_particles(it,ib);

        if (rank >= 1) {

	  for (int ip=0; ip<np; ip++) {

	    int ipdv = ip*dv;
	    int ipdp = ip*dp;
	    int ipda = ip*da;

#ifdef DEBUG_UPDATE    
	    v3sum[0]+=std::abs(vx[ipdv]);
	    a3sum[0]+=std::abs(ax[ipda]);
	    v3sum2[0]+=vx[ipdv]*vx[ipdv];
	    a3sum2[0]+=ax[ipda]*ax[ipda];
	    CkPrintf ("DEBUG_UPDATE x %g v %g a %g\n",x[ipdp],vx[ipdv],ax[ipda]);
#endif	  
	    vx[ipdv] = cvv*vx[ipdv] + cva*ax[ipda];
	    x [ipdp] += cp*vx[ipdv];
	    vx[ipdv] = cvv*vx[ipdv] + cva*ax[ipda];

	  }
        }
        if (rank >= 2) {

	  for (int ip=0; ip<np; ip++) {

	    int ipdv = ip*dv;
	    int ipdp = ip*dp;
	    int ipda = ip*da;

#ifdef DEBUG_UPDATE    
	    v3sum[1]+=std::abs(vy[ipdv]);
	    a3sum[1]+=std::abs(ay[ipda]);
	    v3sum2[1]+=vy[ipdv]*vy[ipdv];
	    a3sum2[1]+=ay[ipda]*ay[ipda];
#endif	  
	    vy[ipdv] = cvv*vy[ipdv] + cva*ay[ipda];
	    y [ipdp] += cp*vy[ipdv];
	    vy[ipdv] = cvv*vy[ipdv] + cva*ay[ipda];

	  }
	
        }
        if (rank >= 3) {

	  for (int ip=0; ip<np; ip++) {

	    int ipdv = ip*dv;
	    int ipdp = ip*dp;
	    int ipda = ip*da;

#ifdef DEBUG_UPDATE    
	    v3sum[2]+=std::abs(vz[ipdv]);
	    a3sum[2]+=std::abs(az[ipda]);
	    v3sum2[2]+=vz[ipdv]*vz[ipdv];
	    a3sum2[2]+=az[ipda]*az[ipda];
#endif	  

	    vz[ipdv] = cvv*vz[ipdv] + cva*az[ipda];
	    z [ipdp] += cp*vz[ipdv];
	    vz[ipdv] = cvv*vz[ipdv] + cva*az[ipda];

	  }
        }
      }
    } // end loop over particle types
    
#ifdef DEBUG_UPDATE    
    CkPrintf ("DEBUG_UPDATE asum %g %g %g vsum %g %g %g\n",
	      a3sum[0],a3sum[1],a3sum[2],
	      v3sum[0],v3sum[1],v3sum[2]);
    CkPrintf ("DEBUG_UPDATE asum2 %g %g %g vsum2 %g %g %g\n",
	      a3sum2[0],a3sum2[1],a3sum2[2],
	      v3sum2[0],v3sum2[1],v3sum2[2]);
#endif    
  }

  block->compute_done(); 
  
}

//----------------------------------------------------------------------

double EnzoMethodPmUpdate::timestep ( Block * block ) const throw()
{
  TRACE_PM("timestep()");

  const int rank = cello::rank();

  double dt = std::numeric_limits<double>::max();

  if (block->is_leaf()) {

    Particle particle = block->data()->particle();
    ParticleDescr * particle_descr = cello::particle_descr();
    Grouping * particle_groups = particle_descr->groups();
    Field    field    = block->data()->field();

    double xm,ym,zm;
    double xp,yp,zp;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);

    int nx,ny,nz;
    field.size(&nx,&ny,&nz);

    double hx = (xp-xm)/nx;
    double hy = (yp-ym)/ny;
    double hz = (zp-zm)/nz;
    
    // Adjust for expansion terms if any
    EnzoPhysicsCosmology * cosmology = enzo::cosmology();
    
    if (cosmology) {
      enzo_float cosmo_a=1.0,cosmo_dadt=0.0;
      double time = block->time();
      double dt   = block->dt();
      cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
      hx *= cosmo_a;
      hy *= cosmo_a;
      hz *= cosmo_a;
    }


    int num_mass = particle_groups->size("has_mass");

    for (int ipt = 0; ipt < num_mass; ipt++){

      std::string particle_type = particle_groups->item("has_mass",ipt);
      int it = particle.type_index (particle_type);
      int nb = particle.num_batches (it);

      int ia_vx = particle.attribute_index (it, "vx");
      int ia_vy = particle.attribute_index (it, "vy");
      int ia_vz = particle.attribute_index (it, "vz");

      int ia_ax = particle.attribute_index (it, "ax");
      int ia_ay = particle.attribute_index (it, "ay");
      int ia_az = particle.attribute_index (it, "az");

      int dv = particle.stride(it, ia_vx);
      int da = particle.stride(it, ia_ax);


      for (int ib=0; ib<nb; ib++) {
        int np = particle.num_particles(it,ib);

        if (rank >= 1) {
	  enzo_float * vx = (enzo_float *) 
	    particle.attribute_array (it, ia_vx, ib); 
	  enzo_float * ax = (enzo_float *) 
	    particle.attribute_array (it, ia_ax, ib); 
	  for (int ip=0; ip<np; ip++) {
	    double v = fabs(vx[ip*dv]);
	    double a = fabs(ax[ip*da]);
	    double dt_v = hx /MAX(v,1e-6);
	    double dt_a = sqrt(2.0*hx/MAX(a,1e-6));
	    dt = MIN(dt,dt_v);
	    dt = MIN(dt,dt_a);
	  }
        }

        if (rank >= 2) {
	  enzo_float * vy = (enzo_float *) 
	    particle.attribute_array (it, ia_vy, ib); 
	  enzo_float * ay = (enzo_float *) 
	    particle.attribute_array (it, ia_ay, ib); 
	  for (int ip=0; ip<np; ip++) {
	    double v = fabs(vy[ip*dv]);
  	    double a = fabs(ay[ip*da]);
	    double dt_v = hy /MAX(v,1e-6);
	    double dt_a = sqrt(2.0*hy/MAX(a,1e-6));
	    dt = MIN(dt,dt_v);
	    dt = MIN(dt,dt_a);
	  }
        }

        if (rank >= 3) {
	  enzo_float * vz = (enzo_float *) 
	    particle.attribute_array (it, ia_vz, ib); 
	  enzo_float * az = (enzo_float *) 
	    particle.attribute_array (it, ia_az, ib); 
	  for (int ip=0; ip<np; ip++) {
	    double v = fabs(vz[ip*dv]);
	    double a = fabs(az[ip*da]);
	    double dt_v = hz/MAX(v,1e-6);
	    double dt_a = sqrt(2.0*hz/MAX(a,1e-6));
	    dt = MIN(dt,dt_v);
	    dt = MIN(dt,dt_a);
	  }
        }

      }
    } // end loop over particle types
  }
  
  dt = MIN(dt,max_dt_);
  return dt;
}
