// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-12
/// @brief    Implementation of the FieldFace class

#include "cello.hpp"
#include "data.hpp"

#define FORTRAN_STORE

long FieldFace::counter[CONFIG_NODE_SIZE] = {0};

#define FORTRAN_NAME(NAME) NAME##_

extern "C" void FORTRAN_NAME(field_face_store_4)
  (float * field, float * array, int * m3, int * n3, int * accumulate);
extern "C" void FORTRAN_NAME(field_face_store_8)
  (double * field, double * array, int * m3, int * n3, int * accumulate);
extern "C" void FORTRAN_NAME(field_face_store_16)
  (long double * field, long double * array, int * m3, int * n3, int * accumulate);

enum enum_op_type {
  op_unknown,
  op_load,
  op_store
};


#define CHECK_COARSE(FIELD,i_f)                 \
  { \
    cello_float * coarse = (cello_float *)FIELD.coarse_values(i_f);     \
    ASSERT ("CHECK_COARSE","coarse array is null", (coarse != nullptr)); \
  }


//----------------------------------------------------------------------

FieldFace::FieldFace (int rank) throw()
  : rank_(rank),
    refresh_type_(refresh_unknown),
    prolong_(NULL),
    restrict_(NULL),
    refresh_(NULL),
    new_refresh_(false)
{
  ++counter[cello::index_static()]; 

  for (int i=0; i<3; i++) {
    face_[i] = 0;
    child_[i] = 0;
  }

}

//----------------------------------------------------------------------

FieldFace::~FieldFace() throw ()
{
  --counter[cello::index_static()];

  if (new_refresh_) {
    delete refresh_;
    refresh_ = NULL;
  }
}

//----------------------------------------------------------------------

FieldFace::FieldFace(const FieldFace & field_face) throw ()
  :  refresh_type_(refresh_unknown),
     prolong_(NULL),
     restrict_(NULL),
     refresh_(NULL),
     new_refresh_(false)

{
  ++counter[cello::index_static()];

  copy_(field_face);
}

//----------------------------------------------------------------------

FieldFace & FieldFace::operator= (const FieldFace & field_face) throw ()
/// @param     field_face  Source object of the assignment
///
/// @return    The target assigned object
{
  copy_(field_face);
  return *this;
}

//----------------------------------------------------------------------

void FieldFace::copy_(const FieldFace & field_face)
{
  for (int i=0; i<3; i++) {
    ghost_[i] = field_face.ghost_[i];
    face_[i]  = field_face.face_[i];
    child_[i] = field_face.child_[i];
  }
  refresh_type_ = field_face.refresh_type_;
  restrict_     = field_face.restrict_;
  prolong_      = field_face.prolong_;
  refresh_      = field_face.refresh_;
  // new_refresh_ must not be true in more than one FieldFace to avoid
  // multiple deletes
  new_refresh_  = false;
}

//----------------------------------------------------------------------

void FieldFace::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  p | rank_;
  PUParray(p,face_,3);
  PUParray(p,ghost_,3);
  PUParray(p,child_,3);
  p | refresh_type_;
  p | prolong_;
  p | restrict_;
  p | refresh_;
  p | new_refresh_;
}

//======================================================================

void FieldFace::face_to_array ( Field field, int * n, char ** array) throw()
{
  ASSERT("FieldFace::face_to_array()",
	 "field_src.size() must be > 0",
	 refresh_->any_fields());

  *n = num_bytes_array(field);
  *array = new char [*n];

  ASSERT("FieldFace::face_to_array()",
	 "array size must be > 0",
	 *n > 0);

  face_to_array (field, *array);

}

//----------------------------------------------------------------------
void FieldFace::face_to_array ( Field field,char * array) throw()
{
  size_t index_array = 0;

  auto field_list_src = refresh_->field_list_src();
  auto field_list_dst = refresh_->field_list_dst();

  for (size_t i_f=0; i_f < field_list_src.size(); i_f++) {

    CHECK_COARSE(field,i_f);

    const size_t index_field = field_list_src[i_f];
    
    precision_type precision = field.precision(index_field);

    void * field_face = field.values(index_field);

    char * array_face  = &array[index_array];

    int m3[3],g3[3],c3[3];

    field.dimensions(index_field,&m3[0],&m3[1],&m3[2]);
    field.ghost_depth(index_field,&g3[0],&g3[1],&g3[2]);
    field.centering(index_field,&c3[0],&c3[1],&c3[2]);

    int index_src = field_list_src[i_f];
    int index_dst = field_list_dst[i_f];
    
    const bool accumulate = accumulate_(index_src,index_dst);

    int i3[3], n3[3];
    field.size(n3,n3+1,n3+2);
    Box box(rank_,n3,g3);
    set_box_(&box);

    if (refresh_type_ == refresh_fine) {
      const int pad = prolong()->padding();
      if (pad>0) box.set_padding(pad);
    }

    box_adjust_accumulate_(&box,accumulate,g3);

    box.compute_region();
    // limits for Send block 
    box.get_limits(Box::BlockType::send,i3,n3);
    n3[0] -= i3[0];
    n3[1] -= i3[1];
    n3[2] -= i3[2];

    // scale by density if needed to convert to conservative form
    mul_by_density_(field,index_field,i3,n3,m3);

    if (refresh_type_ == refresh_coarse) {

      // Restrict field to array

      int nc3[3] = { (n3[0]+1)/2, (n3[1]+1)/2,(n3[2]+1)/2 };

      int i3_array[3] = {0,0,0};

      index_array += restrict()->apply
	(precision, 
	 array_face,nc3,i3_array,nc3, 
	 field_face,m3,i3, n3);

    } else {

      union { float * a4; double * a8; long double * a16; };
      union { float * f4; double * f8;long double * f16;  };
      a4 = (float *) array_face;
      f4 = (float *) field_face;
      
      // Copy field to array
      
      if (precision == precision_single) {
	index_array += load_ ( a4,  f4,  m3,n3,i3, accumulate);
      } else if (precision == precision_double) {
	index_array += load_ ( a8,  f8,  m3,n3,i3, accumulate);
      } else if (precision == precision_quadruple) {
	index_array += load_ ( a16, f16, m3,n3,i3, accumulate);
      } else {
	ERROR("FieldFace::face_to_array", "Unsupported precision");
      }
    }

    // unscale by density if needed to convert back from conservative form
    div_by_density_(field,index_field,i3,n3,m3);
  }

}

//----------------------------------------------------------------------

void FieldFace::array_to_face (char * array, Field field) throw()
{
  size_t index_array = 0;

  auto field_list_src = refresh_->field_list_src();
  auto field_list_dst = refresh_->field_list_dst();
  
  for (size_t i_f=0; i_f < field_list_dst.size(); i_f++) {

    CHECK_COARSE(field,i_f);
    size_t index_field = field_list_dst[i_f];

    precision_type precision = field.precision(index_field);

    char * field_ghost = field.values(index_field);
    
    char * array_ghost  = array + index_array;

    int m3[3],g3[3],c3[3];

    field.dimensions(index_field,&m3[0],&m3[1],&m3[2]);
    field.ghost_depth(index_field,&g3[0],&g3[1],&g3[2]);
    field.centering(index_field,&c3[0],&c3[1],&c3[2]);

    int index_src = field_list_src[i_f];
    int index_dst = field_list_dst[i_f];
    const bool accumulate = accumulate_(index_src,index_dst);

    int i3[3], n3[3];

    // adjust face relative to sender
    invert_face();

    field.size(n3,n3+1,n3+2);
    Box box (rank_,n3,g3);
    set_box_(&box);
    invert_face();

    if (refresh_type_ == refresh_fine) {
      const int pad = prolong()->padding();
      if (pad>0) box.set_padding(pad);
    }

    box_adjust_accumulate_(&box,accumulate,g3);

    box.get_limits(Box::BlockType::receive,i3,n3,BoxType_receive);
    n3[0] -= i3[0];
    n3[1] -= i3[1];
    n3[2] -= i3[2];

    if (refresh_type_ == refresh_fine) {

      // Prolong array to field

      bool need_padding = (g3[0]%2==1) || (g3[1]%2==1) || (g3[2]%2==1);

      ASSERT("FieldFace::array_to_face()",
	     "Odd ghost zones not implemented yet: prolong needs padding",
	     ! need_padding);

      int nc3[3] = { (n3[0]+1)/2, (n3[1]+1)/2, (n3[2]+1)/2 };

      int i3_array[3] = {0,0,0};

      const int pad = prolong()->padding();

      if (pad == 0) {
        index_array += prolong()->apply
          (precision, 
           field_ghost,m3,i3,       n3,
           array_ghost,nc3,i3_array, nc3,
           accumulate);
      } else {

        const int nf = refresh_->field_list_src().size();
        const int nxc = rank_ >= 1 ? (nc3[0]+2*pad) : 1;
        const int nyc = rank_ >= 2 ? (nc3[1]+2*pad) : 1;
        const int nzc = rank_ >= 3 ? (nc3[2]+2*pad) : 1;
        const int oxc = (rank_ >= 1) ? 1 : 0;
        const int oyc = (rank_ >= 2) ? 1 : 0;
        const int ozc = (rank_ >= 3) ? 1 : 0;
        FieldData * field_data = field.field_data();
        cello_float * padded_face_vector = field_data->padded_array_allocate
          (face_[0],face_[1],face_[2],nf,nxc,nyc,nzc);

        const int na = field_data->padded_array_dimensions
          (face_[0],face_[1],face_[2]);
                                                                      
        ASSERT5("EnzoProlong::apply",
                "Array size mismatch %d * %d * %d  !=  %d / %d\n",
                nxc,nyc,nzc,na, nf,
                (nxc*nyc*nzc) == (na/nf));

        const int nc = nxc*nyc*nzc;
        cello_float * padded_array = &padded_face_vector[i_f*nc];

        cello::copy(padded_array,
                    nxc,nyc,nzc, oxc,oyc,ozc,
                    (cello_float *) array_ghost,
                    nc3[0],nc3[1],nc3[2],0,0,0,
                    nc3[0],nc3[1],nc3[2]);

      }


    } else {

      // Copy array to field
      union { float * as4; double * as8; long double * as16; };
      union { float * fd4; double * fd8; long double * fd16; };
      as4 = (float *) array_ghost;
      fd4 = (float *) field_ghost;
      
      // Copy field to array

      if (precision == precision_single) {
	index_array += store_ ( fd4,  as4,  m3,n3,i3, accumulate);
      } else if (precision == precision_double) {
	index_array += store_ ( fd8,  as8,  m3,n3,i3, accumulate);
      } else if (precision == precision_quadruple) {
	index_array += store_ ( fd16, as16, m3,n3,i3, accumulate);
      } else {
	ERROR("FieldFace::array_to_face()", "Unsupported precision");
      }
    }

    // unscale by density if needed to convert back from conservative form
    div_by_density_(field,index_field,i3,n3,m3);

  }
}

//----------------------------------------------------------------------

void FieldFace::face_to_face (Field field_src, Field field_dst)
{
  
  auto field_list_src = refresh_->field_list_src();
  auto field_list_dst = refresh_->field_list_dst();
  
  
  cello_float * padded_face_vector = nullptr;

  for (size_t i_f=0; i_f < field_list_src.size(); i_f++) {

    CHECK_COARSE(field_src,i_f);
    size_t index_src = field_list_src[i_f];
    size_t index_dst = field_list_dst[i_f];

    int m3[3]={0},g3[3]={0},c3[3]={0};

    field_src.dimensions(index_src,&m3[0],&m3[1],&m3[2]);
    field_src.ghost_depth(index_src,&g3[0],&g3[1],&g3[2]);
    field_src.centering(index_src,&c3[0],&c3[1],&c3[2]);
    
    const bool accumulate = accumulate_(index_src,index_dst);

    int is3[3], ns3[3];

    int n3[3];
    field_src.size(n3,n3+1,n3+2);
    Box box (rank_,n3,g3);
    set_box_(&box);

    if (refresh_type_ == refresh_fine) {
      if (prolong()->padding()) {
        box.set_padding(prolong()->padding());
      }      
    }

    box_adjust_accumulate_(&box,accumulate,g3);

    box.get_limits(Box::BlockType::send,is3,ns3);
    ns3[0] -= is3[0];
    ns3[1] -= is3[1];
    ns3[2] -= is3[2];

    int id3[3], nd3[3];

    box.get_limits(Box::BlockType::receive,id3,nd3,BoxType_receive);
    nd3[0] -= id3[0];
    nd3[1] -= id3[1];
    nd3[2] -= id3[2];

    // Adjust loop limits if accumulating to include ghost zones
    // on neighbor axes

    precision_type precision = field_src.precision(index_src);
    
    char * values_src = field_src.values(index_src);
    char * values_dst = field_dst.values(index_dst);
    
    // scale by density if needed to convert to conservative form
    mul_by_density_(field_src,index_src,is3,ns3,m3);

    if (refresh_type_ == refresh_fine) {

      // Prolong field

      bool need_padding = (g3[0]%2==1) || (g3[1]%2==1) || (g3[2]%2==1);

      ASSERT("FieldFace::face_to_face()",
	     "Odd ghost zones not implemented yet: prolong needs padding",
	     ! need_padding);

      const int pad = prolong()->padding();

      if (pad == 0) {

        prolong()->apply (precision, 
                        values_dst,m3,id3, nd3,
                        values_src,m3,is3, ns3,
                        accumulate);

      } else {
        
        int ia3[3], na3[3];
        box.get_limits(Box::BlockType::receive,ia3,na3,BoxType_receive);
        na3[0] -= ia3[0];
        na3[1] -= ia3[1];
        na3[2] -= ia3[2];

        const int nf = refresh_->field_list_src().size();
        const int nxc = rank_ >= 1 ? (ns3[0]+2*pad) : 1;
        const int nyc = rank_ >= 2 ? (ns3[1]+2*pad) : 1;
        const int nzc = rank_ >= 3 ? (ns3[2]+2*pad) : 1;
        const int oxc = (rank_ >= 1) ? 1 : 0;
        const int oyc = (rank_ >= 2) ? 1 : 0;
        const int ozc = (rank_ >= 3) ? 1 : 0;
        FieldData * field_data = field_dst.field_data();
        if (i_f == 0)  {
          padded_face_vector = field_data->padded_array_allocate
            (face_[0],face_[1],face_[2],nf,nxc,nyc,nzc);
        }
        const int na = field_data->padded_array_dimensions
          (face_[0],face_[1],face_[2]);

        ASSERT5("EnzoProlong::apply",
                "Array size mismatch %d * %d * %d  !=  %d / %d\n",
                nxc,nyc,nzc,na, nf,
                ((nxc*nyc*nzc) == na));

        const int nc = nxc*nyc*nzc;
        cello_float * padded_array = &padded_face_vector[i_f*nc];

        print("face_to_face");
        CkPrintf ("DEBUG_COPY padded_array %p i_f %d nc %d\n",
                  padded_array,i_f,nc);
        CkPrintf ("DEBUG_COPY nxc,nyc,nzc %d %d %d\n", nxc,nyc,nzc);
        CkPrintf ("DEBUG_COPY oxc,oyc,ozc %d %d %d\n", oxc,oyc,ozc);
        CkPrintf ("DEBUG_COPY values_src  %p\n",values_src);
        CkPrintf ("DEBUG_COPY m3[]        %d %d %d\n", m3[0],m3[1],m3[2]);
        CkPrintf ("DEBUG_COPY is3[]       %d %d %d\n", is3[0],is3[1],is3[2]);
        CkPrintf ("DEBUG_COPY ns3[]       %d %d %d\n", ns3[0],ns3[1],ns3[2]);
        fflush(stdout);
        for (int i=0; i<nxc*nyc*nzc; i++) {padded_array[i] = 0.0;}
        cello::copy(padded_array,
                    nxc,nyc,nzc, oxc,oyc,ozc,
                    (cello_float *) values_src,
                    m3[0],m3[1],m3[2],is3[0],is3[1],is3[2],
                    ns3[0],ns3[1],ns3[2]);
      }

    } else if (refresh_type_ == refresh_coarse) {

      // Restrict field

      restrict()->apply (precision,
                         values_dst,m3,id3, nd3,
                         values_src,m3,is3, ns3,
                         accumulate);

    } else {

      // Copy faces to ghosts

      union { float * fs4; double * fs8; long double * fs16; };
      union { float * fd4; double * fd8; long double * fd16; };
      fs4 = (float *) values_src;
      fd4 = (float *) values_dst;
      
      // Copy field to array

      if (precision == precision_single) {
	copy_ ( fd4, m3,nd3,id3,fs4, m3, ns3,is3,accumulate);
      } else if (precision == precision_double) {
	copy_ ( fd8, m3,nd3,id3,fs8, m3, ns3,is3,accumulate);
      } else if (precision == precision_quadruple) {
	copy_ ( fd16,m3,nd3,id3,fs16,m3,ns3,is3,accumulate);
      } else {
	ERROR("FieldFace::face_to_face()", "Unsupported precision");
      }
    }
    // unscale by density if needed to convert back from conservative form
    div_by_density_(field_src,index_src,is3,ns3,m3);
    div_by_density_(field_dst,index_dst,id3,nd3,m3);
  }
}

//----------------------------------------------------------------------

int FieldFace::num_bytes_array(Field field) throw()
{
  int array_size = 0;

  auto field_list_src = refresh_->field_list_src();
  auto field_list_dst = refresh_->field_list_dst();

  for (size_t i_f=0; i_f < field_list_src.size(); i_f++) {

    CHECK_COARSE(field,i_f);
    size_t index_field = field_list_src[i_f];

    precision_type precision = field.precision(index_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    int m3[3],g3[3],c3[3];

    field.dimensions (index_field,&m3[0],&m3[1],&m3[2]);
    field.ghost_depth(index_field,&g3[0],&g3[1],&g3[2]);
    field.centering(index_field,&c3[0],&c3[1],&c3[2]);

    int index_src = field_list_src[i_f];
    int index_dst = field_list_dst[i_f];
    const bool accumulate = accumulate_(index_src,index_dst);

    int i3[3], n3[3];
    field.size(n3,n3+1,n3+2);
    Box box (rank_,n3,g3);
    set_box_(&box);

    if (refresh_type_ == refresh_fine) {
      const int pad = prolong()->padding();
      if (pad>0) box.set_padding(pad);
    }

    box_adjust_accumulate_(&box,accumulate,g3);
    box.get_limits(Box::BlockType::send,i3,n3);
    n3[0] -= i3[0];
    n3[1] -= i3[1];
    n3[2] -= i3[2];

    array_size += n3[0]*n3[1]*n3[2]*bytes_per_element;

  }

  ASSERT("FieldFace::num_bytes_array()",
	 "array_size must be > 0, maybe field_list.size() is 0?",
	 array_size);

  return array_size;

}

//----------------------------------------------------------------------

int FieldFace::data_size () const
{
  int count = 0;

  count += 3*sizeof(int);  // face_[3]
  count += 3*sizeof(bool); // ghost_[3]
  count += 3*sizeof(int);  // child_[3];
  count += 1*sizeof(int);  // refresh_type_ (restrict,prolong,copy)
  count += refresh_->data_size(); // refresh_

  return count;

}

//----------------------------------------------------------------------

char * FieldFace::save_data (char * buffer) const
{
  char * p = buffer;
  int n;

  memcpy(p,face_, n=3*sizeof(int));  
  p+=n;

  memcpy(p,ghost_,n=3*sizeof(bool)); 
  p+=n;

  memcpy(p,child_,n=3*sizeof(int));  
  p+=n;

  memcpy(p,&refresh_type_,n=sizeof(int));  
  p+=n;

  p = refresh_->save_data(p);

  ASSERT2("FieldFace::save_data()",
	  "Buffer has size %ld but expecting size %d",
	  (p-buffer),data_size(),
	  ((p-buffer) == data_size()));
  
  return p;
}

//----------------------------------------------------------------------

char * FieldFace::load_data (char * buffer)
{

  char * p = buffer;
  int n;

  memcpy(face_,p, n=3*sizeof(int));
  p+=n;

  memcpy(ghost_,p,n=3*sizeof(bool));
  p+=n;

  memcpy(child_,p,n=3*sizeof(int));
  p+=n;

  memcpy(&refresh_type_,p,n=sizeof(int));
  p+=n;

  Refresh * refresh = new Refresh;
  set_refresh(refresh,true);

  p = refresh_->load_data(p);

  ASSERT2("FieldFace::save_data()",
	  "Buffer has size %ld but expecting size %d",
	  (p-buffer),data_size(),
	  ((p-buffer) == data_size()));

  return p;
}

//======================================================================

template<class T>
size_t FieldFace::load_
( T * array_face, const T * field_face, 
  int m3[3], int n3[3],int i3[3], bool accumulate ) throw()
{
  // NOTE: don't check accumulate since loading array; accumulate
  // is handled in corresponding store_() at the receiving end
  // add values

  for (int iz=0; iz <n3[2]; iz++)  {
    int kz = iz+i3[2];
    for (int iy=0; iy < n3[1]; iy++) {
      int ky = iy+i3[1];
      for (int ix=0; ix < n3[0]; ix++) {
	int kx = ix+i3[0];
	int index_array = ix +   n3[0]*(iy +   n3[1] * iz);
	int index_field = kx + m3[0]*(ky + m3[1] * kz);
	array_face[index_array] = field_face[index_field];
      }
    }
  }

  return (sizeof(T) * n3[0] * n3[1] * n3[2]);

}

//----------------------------------------------------------------------

template<class T> size_t FieldFace::store_
( T * ghost, const T * array,
  int m3[3], int n3[3],int i3[3], bool accumulate) throw()
{

#ifdef FORTRAN_STORE
  const bool use_fortran_store = true;
#else  
  const bool use_fortran_store = false;
#endif

  // This is to get around a bug on SDSC Comet where this function
  // crashes with -O3 (See bugzilla report #90)

  union {
    float *       ghost_4;
    double *      ghost_8;
    long double * ghost_16;
  };
  union {
    float *       array_4;
    double *      array_8;
    long double * array_16;
  };

  ghost_4 = (float *) ghost;
  array_4 = (float *) array;
  
  int im = i3[0] + m3[0]*(i3[1] + m3[1]*i3[2]);

  int iaccumulate = accumulate ? 1 : 0;

  if (use_fortran_store &&
      (sizeof(T) != sizeof(long double)) ) {

    if (sizeof(T)==sizeof(float)) {
      FORTRAN_NAME(field_face_store_4)(ghost_4 + im,   array_4, m3,n3,
                                       &iaccumulate);
    } else if (sizeof(T)==sizeof(double)) {
      FORTRAN_NAME(field_face_store_8)(ghost_8 + im,   array_8, m3,n3,
                                       &iaccumulate);
      // SUSPECTED ERROR IN COMPILER: seg fault for quad test in
      // test_FieldFace
      //
      //  } else if (sizeof(T)==sizeof(long double)) {
      //    FORTRAN_NAME(field_face_store_16)(ghost_16 + im, array_16, m3,n3,
      //				      &iaccumulate);
    } else {
      ERROR1 ("FieldFace::store_()",
              "unknown float precision sizeof(T) = %lu\n",sizeof(T));
    }
  } else {

    if (accumulate) {
      // add values
      for (int iz=0; iz <n3[2]; iz++)  {
        int kz = iz+i3[2];
        for (int iy=0; iy < n3[1]; iy++) {
          int ky = iy+i3[1];
          for (int ix=0; ix < n3[0]; ix++) {
            int kx = ix+i3[0];
            int index_array = ix + n3[0]*(iy + n3[1] * iz);
            int index_field = kx + m3[0]*(ky + m3[1] * kz);
            ghost[index_field] += array[index_array];
          }
        }
      }
    } else {
      // copy values
      for (int iz=0; iz <n3[2]; iz++)  {
        int kz = iz+i3[2];
        for (int iy=0; iy < n3[1]; iy++) {
          int ky = iy+i3[1];
          for (int ix=0; ix < n3[0]; ix++) {
            int kx = ix+i3[0];
            int index_array = ix + n3[0]*(iy + n3[1] * iz);
            int index_field = kx + m3[0]*(ky + m3[1] * kz);
            ghost[index_field] = array[index_array];
          }
        }
      }
    }
  }

  return (sizeof(T) * n3[0] * n3[1] * n3[2]);

}

//----------------------------------------------------------------------

template<class T> void FieldFace::copy_
( T       * vd, int md3[3],int nd3[3],int id3[3],
  const T * vs, int ms3[3],int ns3[3],int is3[3],
  bool accumulate) throw()
{
  const int is0 = is3[0] + ms3[0]*(is3[1] + ms3[1]*is3[2]);
  const int id0 = id3[0] + md3[0]*(id3[1] + md3[1]*id3[2]);
  T * vd0 = vd + id0;
  const T * vs0 = vs + is0;
  const int msx = ms3[0];
  const int msy = ms3[1];
  const int mdx = md3[0];
  const int mdy = md3[1];
  if (accumulate) {
    for (int iz=0; iz <ns3[2]; iz++)  {
      for (int iy=0; iy < ns3[1]; iy++) {
	for (int ix=0; ix < ns3[0]; ix++) {
	  int i_src = ix + msx*(iy + msy*iz);
	  int i_dst = ix + mdx*(iy + mdy*iz);
	  vd0[i_dst] += vs0[i_src];
	}
      }
    }
  } else {
    for (int iz=0; iz <ns3[2]; iz++)  {
      for (int iy=0; iy < ns3[1]; iy++) {
	for (int ix=0; ix < ns3[0]; ix++) {
	  int i_src = ix + msx*(iy + msy*iz);
	  int i_dst = ix + mdx*(iy + mdy*iz);
	  vd0[i_dst] = vs0[i_src];
	}
      }
    }
  }
}

//----------------------------------------------------------------------

void FieldFace::print(const char * message)
{
  CkPrintf (" FieldFace   %s %p\n",message,(void*)this);
  CkPrintf ("    face_    %d %d %d\n",face_[0],face_[1],face_[2]);
  CkPrintf ("    ghost_   %d %d %d\n",ghost_[0],ghost_[1],ghost_[2]);
  CkPrintf ("    child_   %d %d %d\n",child_[0],child_[1],child_[2]);
  CkPrintf ("    prolong_ %p\n",(void*)prolong_);
  CkPrintf ("    restrict_ %p\n",(void*)restrict_);
  CkPrintf ("    refresh_type_ %d\n",refresh_type_);
  if (refresh_) refresh_->print();
}

//----------------------------------------------------------------------

void FieldFace::set_field_list(std::vector<int> field_list)
{
  refresh_->set_field_list(field_list);
}

//----------------------------------------------------------------------

bool FieldFace::accumulate_(int index_src, int index_dst) const
{
  return ((index_src != index_dst) && refresh_->accumulate());
}

//----------------------------------------------------------------------

void FieldFace::mul_by_density_
(Field field, int index_field,
 const int i3[3], const int n3[3], const int m3[3])
{
  if (field.is_temporary(index_field)) return;
  
  precision_type precision = field.precision(index_field);
  void * field_face = field.values(index_field);

  Grouping * groups = cello::field_groups();

  void * field_density = field.values("density");
  
  const std::string field_name = field.field_name(index_field);

  const bool scale_by_density =
    (refresh_type_ != refresh_same) &&
    groups->is_in (field_name,"make_field_conservative");
  if (scale_by_density) {
    union { float * d4; double * d8; long double * d16; };
    union { float * f4; double * f8;long double * f16;  };
    d4 = (float *) field_density;
    f4 = (float *) field_face;

    if (precision == precision_single) {
      for (int iz=i3[2]; iz<i3[2]+n3[2]; iz++) {
        for (int iy=i3[1]; iy<i3[1]+n3[1]; iy++) {
          for (int ix=i3[0]; ix<i3[0]+n3[0]; ix++) {
            const int i=ix + m3[0]*(iy + m3[1]*iz);
            f4[i] *= d4[i];
          }
        }
      }
    } else if (precision == precision_double) {
      for (int iz=i3[2]; iz<i3[2]+n3[2]; iz++) {
        for (int iy=i3[1]; iy<i3[1]+n3[1]; iy++) {
          for (int ix=i3[0]; ix<i3[0]+n3[0]; ix++) {
            const int i=ix + m3[0]*(iy + m3[1]*iz);
            f8[i] *= d8[i];
          }
        }
      }
    } else if (precision == precision_quadruple) {
      for (int iz=i3[2]; iz<i3[2]+n3[2]; iz++) {
        for (int iy=i3[1]; iy<i3[1]+n3[1]; iy++) {
          for (int ix=i3[0]; ix<i3[0]+n3[0]; ix++) {
            const int i=ix + m3[0]*(iy + m3[1]*iz);
            f16[i] *= d16[i];
          }
        }
      }
    } else {
      ERROR("FieldFace::mul_by_density_()", "Unsupported precision");
    }
  }
}

//----------------------------------------------------------------------

void FieldFace::div_by_density_
(Field field, int index_field,
 const int i3[3], const int n3[3], const int m3[3])
{
      
  if (field.is_temporary(index_field)) return;

  precision_type precision = field.precision(index_field);
  void * field_face = field.values(index_field);

  Grouping * groups = cello::field_groups();

  void * field_density = field.values("density");
  
  const std::string field_name = field.field_name(index_field);

  const bool scale_by_density =
    (refresh_type_ != refresh_same) &&
    groups->is_in (field_name,"make_field_conservative");
  if (scale_by_density) {
    union { float * d4; double * d8; long double * d16; };
    union { float * f4; double * f8;long double * f16;  };
    d4 = (float *) field_density;
    f4 = (float *) field_face;

    if (precision == precision_single) {
      for (int iz=i3[2]; iz<i3[2]+n3[2]; iz++) {
        for (int iy=i3[1]; iy<i3[1]+n3[1]; iy++) {
          for (int ix=i3[0]; ix<i3[0]+n3[0]; ix++) {
            const int i=ix + m3[0]*(iy + m3[1]*iz);
            f4[i] /= d4[i];
          }
        }
      }
    } else if (precision == precision_double) {
      for (int iz=i3[2]; iz<i3[2]+n3[2]; iz++) {
        for (int iy=i3[1]; iy<i3[1]+n3[1]; iy++) {
          for (int ix=i3[0]; ix<i3[0]+n3[0]; ix++) {
            const int i=ix + m3[0]*(iy + m3[1]*iz);
            f8[i] /= d8[i];
          }
        }
      }
    } else if (precision == precision_quadruple) {
      for (int iz=i3[2]; iz<i3[2]+n3[2]; iz++) {
        for (int iy=i3[1]; iy<i3[1]+n3[1]; iy++) {
          for (int ix=i3[0]; ix<i3[0]+n3[0]; ix++) {
            const int i=ix + m3[0]*(iy + m3[1]*iz);
            f16[i] /= d16[i];
          }
        }
      }
    } else {
      ERROR("FieldFace::div_by_density_()", "Unsupported precision");
    }
  }
}

//----------------------------------------------------------------------

void FieldFace::set_box_(Box * box)
{
  const int level =
    (refresh_type_==refresh_coarse) ? -1
    : (refresh_type_==refresh_same) ?  0 : +1;

  box->set_block(BoxType_receive,level,face_,child_);

  if (prolong() && prolong()->padding() > 0)
    box->set_padding(prolong()->padding());
  box->compute_region();
}

//----------------------------------------------------------------------

void FieldFace::box_adjust_accumulate_ (Box * box, int accumulate, int g3[3])
{
  if (accumulate) {
    int gs3[3] = {(face_[0]!=0)?g3[0]:0,
                  (face_[1]!=0)?g3[1]:0,
                  (face_[2]!=0)?g3[2]:0};
    box->set_send_ghosts(gs3);
                           
  } else {
    int gs3[3] = {(ghost_[0]&&face_[0]==0)?g3[0]:0,
                  (ghost_[1]&&face_[1]==0)?g3[1]:0,
                  (ghost_[2]&&face_[2]==0)?g3[2]:0};
    box->set_send_ghosts (gs3);
  }
  box->compute_block_start(BoxType_receive);
  box->compute_region();
}
