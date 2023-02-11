// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#include "enzo.hpp"

#include "charm_enzo.hpp"

// #define DEBUG_FACTORY
// #define TRACE_FACTORY

//----------------------------------------------------------------------

void EnzoFactory::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Factory::pup(p);
}

//----------------------------------------------------------------------

IoBlock * EnzoFactory::create_io_block () const throw()
{
  return new IoEnzoBlock;
}

//----------------------------------------------------------------------

CProxy_Block EnzoFactory::new_block_proxy
(
 DataMsg * data_msg,
 int nbx, int nby, int nbz) const throw()
{
  CProxy_EnzoBlock enzo_block_array;

  //  CProxy_MappingTree array_map  = CProxy_MappingTree::ckNew(nbx,nby,nbz);
  CProxy_MappingArray array_map  = CProxy_MappingArray::ckNew(nbx,nby,nbz);

  CkArrayOptions opts;
  opts.setMap(array_map);

  //###################
  // int ep = CkIndex_EnzoBlock::p_say_hello();
  // Index ind(1, 0, 1);
  // CkCallback cb = CkCallback(ep, enzo_block_array[ind]);
  // opts.setInitCallback(cb);
  //###################

  enzo_block_array = CProxy_EnzoBlock::ckNew(opts);

  return enzo_block_array;
  
}

//----------------------------------------------------------------------
void EnzoFactory::create_block_array
(
 DataMsg * data_msg,
 CProxy_Block block_array,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks
 ) const throw()
{
#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY create_block_array %d %d %d\n",nbx,nby,nbz);
#endif  
  CProxy_EnzoBlock enzo_block_array = enzo::block_array();

  int count_adapt;

  int    cycle = 0;
  double time  = 0.0;
  double dt    = 0.0;
  int num_face_level = 0;
  int * face_level = 0;

#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY %s:%d\n",__FILE__,__LINE__); fflush(stdout);
#endif
  for (int ix=0; ix<nbx; ix++) {
    for (int iy=0; iy<nby; iy++) {
      for (int iz=0; iz<nbz; iz++) {

        Index index(ix,iy,iz);
        //#################################################################################################
        // const int nb[3] = {2, 2, 2};
        // std::cout << "Now inserting " << index.bit_string(-1, 3, nb) << " " << nbx << nby << nbz << std::endl;
        //#################################################################################################

        MsgRefine * msg = new MsgRefine 
          (index,
          nx,ny,nz,
          num_field_blocks,
          count_adapt = 0,
          cycle,time,dt,
          refresh_same,
          num_face_level, face_level,
                nullptr);

        msg->set_data_msg(data_msg);

      #ifdef BYPASS_CHARM_MEM_LEAK
        enzo::simulation()->set_msg_refine (index,msg);
        enzo_block_array[index].insert (process_type(CkMyPe()), MsgType::msg_refine);
      #else
        enzo_block_array[index].insert (msg);
      #endif
	// --------------------------------------------------

      }
    }
  }

  //############################################################
  // Index parent_index(0, 0, 0);
  // enzo_block_array[parent_index].p_initialize_children();
  //############################################################
  TRACE1("EnzoFactory::create_block_array = %p",&enzo_block_array);
}

//#######################################################################
void EnzoFactory::add_refined_blocks(int lower[3], int upper[3]) const throw()
{
  std::cout << "spawning " << lower[0] << " to " << upper[0] << std::endl;
  const int nb[3] = {2, 2, 2};
  CProxy_EnzoBlock enzo_block_array = enzo::block_array();

  for (int x = 0; x <= 1; x++){
    for (int y = 0; y <= 1; y++){
      for (int z = 0; z <= 1; z++){
        Index ind(x, y, z);
        enzo_block_array[ind].p_say_hello();
      }
    }
  }

  for (int x = lower[0]; x <= upper[0]; x++){
    for (int y = lower[1]; y <= upper[1]; y++){
      for (int z = lower[2]; z <= upper[2]; z++){
        Index parent_index(x, y, z);
        std::cout << parent_index.bit_string(0, 3, nb) << " is now spawning children" << std::endl;
        int ax, ay, az;
        parent_index.array(&ax, &ay, &az);
        std::cout << ax << " " << ay << " " << az << std::endl;
        enzo_block_array[parent_index].p_initialize_children();
        TRACE1("EnzoFactory::create_block_array = %p",&enzo_block_array);
      }
    }
  }
}
//#######################################################################

//----------------------------------------------------------------------

void EnzoFactory::create_subblock_array
(
 DataMsg * data_msg,
 CProxy_Block block_array,
 int min_level,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks) const throw()
{
#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY create_subblock_array %d %d %d\n",nbx,nby,nbz);
#endif  
  TRACE8("EnzoFactory::create_subblock_array(min_level %d na(%d %d %d) n(%d %d %d) num_field_blocks %d",
	 min_level,nbx,nby,nbz,nx,ny,nz,num_field_blocks);

  if (min_level >= 0) {
    WARNING1("EnzoFactor::create_subblock_array",
	     "Trying to create subblock array with min_level %d >= 0",
	     min_level);
  }

  CProxy_EnzoBlock enzo_block_array = enzo::block_array();

  for (int level = -1; level >= min_level; level--) {

    if (nbx > 1) nbx = ceil(0.5*nbx);
    if (nby > 1) nby = ceil(0.5*nby);
    if (nbz > 1) nbz = ceil(0.5*nbz);

    int count_adapt;

    int    cycle = 0;
    double time  = 0.0;
    double dt    = 0.0;
    int num_face_level = 0;
    int * face_level = 0;

#ifdef TRACE_FACTORY
    CkPrintf ("TRACE_FACTORY %s:%d\n",__FILE__,__LINE__); fflush(stdout);
#endif
    for (int ix=0; ix<nbx; ix++) {
      for (int iy=0; iy<nby; iy++) {
	for (int iz=0; iz<nbz; iz++) {

	  int shift = -level;
 
	  Index index(ix<<shift,iy<<shift,iz<<shift);

	  index.set_level(level);

	  TRACE3 ("inserting %d %d %d",ix,iy,iz);

	  MsgRefine * msg = new MsgRefine 
	    (index,
	     nx,ny,nz,
	     num_field_blocks,
	     count_adapt=0,
	     cycle,time,dt,
	     refresh_same,
	     num_face_level, face_level,
             nullptr);

	  msg->set_data_msg(data_msg);

#ifdef BYPASS_CHARM_MEM_LEAK
	  enzo::simulation()->set_msg_refine (index,msg);
	  enzo_block_array[index].insert (process_type(CkMyPe()), MsgType::msg_refine);
#else
	  enzo_block_array[index].insert (msg);
#endif
	  // --------------------------------------------------

	}
      }
    }

  }

}

//----------------------------------------------------------------------

void EnzoFactory::create_block
(
 DataMsg * data_msg,
 CProxy_Block block_array,
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int count_adapt,
 int cycle, double time, double dt,
 int narray, char * array, int refresh_type,
 int num_face_level,
 int * face_level,
 Adapt * adapt,
 Simulation * simulation,
 int io_reader
 ) const throw()
{
#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY create_block %d %d %d\n");
  index.print("TRACE_FACTORY create_block",2);
#endif  

  TRACE3("EnzoFactory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("EnzoFactory::create_block() num_field_blocks %d  count_adapt %d",
	 num_field_blocks,count_adapt);

  CProxy_EnzoBlock enzo_block_array = (CProxy_EnzoBlock) block_array;

  const int rank = cello::rank();

#ifdef DEBUG_FACTORY
  CkPrintf ("EnzoFactory::create_block: ");
  for (int i=0; i<num_face_level; i++) {
    CkPrintf ("%d ",face_level[i]);
  }
  CkPrintf ("\n");
#endif

   MsgRefine * msg = new MsgRefine 
    (index,
     nx,ny,nz,
     num_field_blocks,
     count_adapt,
     cycle,time,dt,
     refresh_type,
     num_face_level, face_level,
     adapt,
     io_reader);

  msg->set_data_msg(data_msg);

#ifdef BYPASS_CHARM_MEM_LEAK
  enzo::simulation()->set_msg_refine (index,msg);
  enzo_block_array[index].insert ( process_type(CkMyPe()), MsgType::msg_refine );
#else
  enzo_block_array[index].insert (msg);
#endif
}

// //----------------------------------------------------------------------

void EnzoFactory::create_block_check
(
 EnzoMsgCheck * msg_check,
 CProxy_Block block_array,
 Index index
 ) const throw()
{
#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY create_block %d %d %d\n");
  index.print("TRACE_FACTORY create_block",2);
#endif  

  TRACE3("EnzoFactory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("EnzoFactory::create_block() num_field_blocks %d  count_adapt %d",
	 num_field_blocks,count_adapt);

  CProxy_EnzoBlock enzo_block_array = (CProxy_EnzoBlock) block_array;

#ifdef BYPASS_CHARM_MEM_LEAK
  enzo::simulation()->set_msg_check (index,msg_check);
  enzo_block_array[index].insert ( process_type(CkMyPe()), MsgType::msg_check );
#else
  enzo_block_array[index].insert (msg_check);
#endif
}

