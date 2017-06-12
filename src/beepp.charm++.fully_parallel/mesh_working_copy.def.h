/* DEFS: array MeshWorkingCopy: ArrayElement{
MeshWorkingCopy(CkMigrateMessage* impl_msg);
MeshWorkingCopy(void);
void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
void calculate_energy(double kappa, double Dsolvent);
void write_output(const std::string &output_filename);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_MeshWorkingCopy::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: MeshWorkingCopy(CkMigrateMessage* impl_msg);
 */

/* DEFS: MeshWorkingCopy(void);
 */
void CProxyElement_MeshWorkingCopy::insert(int onPE)
{ 
  void *impl_msg = CkAllocSysMsg();
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_MeshWorkingCopy::__idx_MeshWorkingCopy_void,onPE);
}

/* DEFS: void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
 */
void CProxyElement_MeshWorkingCopy::init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(Quaternion &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)impl_noname_2;
    implP|impl_noname_3;
    implP|impl_noname_4;
    implP|impl_noname_5;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(Quaternion &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)impl_noname_2;
    implP|impl_noname_3;
    implP|impl_noname_4;
    implP|impl_noname_5;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_init_marshall2,0);
}

/* DEFS: void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
 */
void CProxyElement_MeshWorkingCopy::calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(vanilla_fmm::CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(vanilla_fmm::CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_calculate_rhs_marshall3,0);
}

/* DEFS: void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
 */
void CProxyElement_MeshWorkingCopy::process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_process_returned_eval_pts_Eval_Message,0);
}

/* DEFS: void calculate_energy(double kappa, double Dsolvent);
 */
void CProxyElement_MeshWorkingCopy::calculate_energy(double kappa, double Dsolvent, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa, double Dsolvent
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    implP|Dsolvent;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
    implP|Dsolvent;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_calculate_energy_marshall5,0);
}

/* DEFS: void write_output(const std::string &output_filename);
 */
void CProxyElement_MeshWorkingCopy::write_output(const std::string &output_filename, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const std::string &output_filename
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)output_filename;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)output_filename;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_write_output_marshall6,0);
}

/* DEFS: MeshWorkingCopy(CkMigrateMessage* impl_msg);
 */
 int CkIndex_MeshWorkingCopy::__idx_MeshWorkingCopy_CkMigrateMessage=0;
void CkIndex_MeshWorkingCopy::_call_MeshWorkingCopy_CkMigrateMessage(void* impl_msg,MeshWorkingCopy * impl_obj)
{
  new (impl_obj) MeshWorkingCopy((CkMigrateMessage*)impl_msg);
}

/* DEFS: MeshWorkingCopy(void);
 */
CkArrayID CProxy_MeshWorkingCopy::ckNew(const CkArrayOptions &opts)
{ 
  void *impl_msg = CkAllocSysMsg();
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_MeshWorkingCopy::__idx_MeshWorkingCopy_void,opts);
}
CkArrayID CProxy_MeshWorkingCopy::ckNew(const int s1)
{ 
  void *impl_msg = CkAllocSysMsg();
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_MeshWorkingCopy::__idx_MeshWorkingCopy_void,CkArrayOptions(s1));
}
 int CkIndex_MeshWorkingCopy::__idx_MeshWorkingCopy_void=0;
void CkIndex_MeshWorkingCopy::_call_MeshWorkingCopy_void(void* impl_msg,MeshWorkingCopy * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  new (impl_obj) MeshWorkingCopy();
}

/* DEFS: void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
 */
void CProxy_MeshWorkingCopy::init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(Quaternion &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)impl_noname_2;
    implP|impl_noname_3;
    implP|impl_noname_4;
    implP|impl_noname_5;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(Quaternion &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)impl_noname_2;
    implP|impl_noname_3;
    implP|impl_noname_4;
    implP|impl_noname_5;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_MeshWorkingCopy::__idx_init_marshall2,0);
}
 int CkIndex_MeshWorkingCopy::__idx_init_marshall2=0;
void CkIndex_MeshWorkingCopy::_call_init_marshall2(void* impl_msg,MeshWorkingCopy * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5*/
  PUP::fromMem implP(impl_buf);
  unsigned int impl_noname_0; implP|impl_noname_0;
  Quaternion impl_noname_1; implP|impl_noname_1;
  Vector impl_noname_2; implP|impl_noname_2;
  unsigned int impl_noname_3; implP|impl_noname_3;
  double impl_noname_4; implP|impl_noname_4;
  double impl_noname_5; implP|impl_noname_5;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->init(impl_noname_0, impl_noname_1, impl_noname_2, impl_noname_3, impl_noname_4, impl_noname_5);
}
int CkIndex_MeshWorkingCopy::_callmarshall_init_marshall2(char* impl_buf,MeshWorkingCopy * impl_obj) {
  /*Unmarshall pup'd fields: unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5*/
  PUP::fromMem implP(impl_buf);
  unsigned int impl_noname_0; implP|impl_noname_0;
  Quaternion impl_noname_1; implP|impl_noname_1;
  Vector impl_noname_2; implP|impl_noname_2;
  unsigned int impl_noname_3; implP|impl_noname_3;
  double impl_noname_4; implP|impl_noname_4;
  double impl_noname_5; implP|impl_noname_5;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->init(impl_noname_0, impl_noname_1, impl_noname_2, impl_noname_3, impl_noname_4, impl_noname_5);
  return implP.size();
}
void CkIndex_MeshWorkingCopy::_marshallmessagepup_init_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5*/
  PUP::fromMem implP(impl_buf);
  unsigned int impl_noname_0; implP|impl_noname_0;
  Quaternion impl_noname_1; implP|impl_noname_1;
  Vector impl_noname_2; implP|impl_noname_2;
  unsigned int impl_noname_3; implP|impl_noname_3;
  double impl_noname_4; implP|impl_noname_4;
  double impl_noname_5; implP|impl_noname_5;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_1");
  implDestP|impl_noname_1;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_2");
  implDestP|impl_noname_2;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_3");
  implDestP|impl_noname_3;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_4");
  implDestP|impl_noname_4;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_5");
  implDestP|impl_noname_5;
}

/* DEFS: void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
 */
void CProxy_MeshWorkingCopy::calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(vanilla_fmm::CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(vanilla_fmm::CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_MeshWorkingCopy::__idx_calculate_rhs_marshall3,0);
}
 int CkIndex_MeshWorkingCopy::__idx_calculate_rhs_marshall3=0;
void CkIndex_MeshWorkingCopy::_call_calculate_rhs_marshall3(void* impl_msg,MeshWorkingCopy * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy; implP|ParallelFMMOctreeProxy;
  vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy; implP|VanillaFMMWorkerProxy;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->calculate_rhs(ParallelFMMOctreeProxy, VanillaFMMWorkerProxy, cb);
}
int CkIndex_MeshWorkingCopy::_callmarshall_calculate_rhs_marshall3(char* impl_buf,MeshWorkingCopy * impl_obj) {
  /*Unmarshall pup'd fields: const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy; implP|ParallelFMMOctreeProxy;
  vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy; implP|VanillaFMMWorkerProxy;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->calculate_rhs(ParallelFMMOctreeProxy, VanillaFMMWorkerProxy, cb);
  return implP.size();
}
void CkIndex_MeshWorkingCopy::_marshallmessagepup_calculate_rhs_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy; implP|ParallelFMMOctreeProxy;
  vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy; implP|VanillaFMMWorkerProxy;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("ParallelFMMOctreeProxy");
  implDestP|ParallelFMMOctreeProxy;
  if (implDestP.hasComments()) implDestP.comment("VanillaFMMWorkerProxy");
  implDestP|VanillaFMMWorkerProxy;
  if (implDestP.hasComments()) implDestP.comment("cb");
  implDestP|cb;
}

/* DEFS: void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
 */
void CProxy_MeshWorkingCopy::process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_MeshWorkingCopy::__idx_process_returned_eval_pts_Eval_Message,0);
}
 int CkIndex_MeshWorkingCopy::__idx_process_returned_eval_pts_Eval_Message=0;
void CkIndex_MeshWorkingCopy::_call_process_returned_eval_pts_Eval_Message(void* impl_msg,MeshWorkingCopy * impl_obj)
{
  impl_obj->process_returned_eval_pts((vanilla_fmm::Eval_Message*)impl_msg);
}

/* DEFS: void calculate_energy(double kappa, double Dsolvent);
 */
void CProxy_MeshWorkingCopy::calculate_energy(double kappa, double Dsolvent, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa, double Dsolvent
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    implP|Dsolvent;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
    implP|Dsolvent;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_MeshWorkingCopy::__idx_calculate_energy_marshall5,0);
}
 int CkIndex_MeshWorkingCopy::__idx_calculate_energy_marshall5=0;
void CkIndex_MeshWorkingCopy::_call_calculate_energy_marshall5(void* impl_msg,MeshWorkingCopy * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double kappa, double Dsolvent*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  double Dsolvent; implP|Dsolvent;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->calculate_energy(kappa, Dsolvent);
}
int CkIndex_MeshWorkingCopy::_callmarshall_calculate_energy_marshall5(char* impl_buf,MeshWorkingCopy * impl_obj) {
  /*Unmarshall pup'd fields: double kappa, double Dsolvent*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  double Dsolvent; implP|Dsolvent;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->calculate_energy(kappa, Dsolvent);
  return implP.size();
}
void CkIndex_MeshWorkingCopy::_marshallmessagepup_calculate_energy_marshall5(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double kappa, double Dsolvent*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  double Dsolvent; implP|Dsolvent;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("kappa");
  implDestP|kappa;
  if (implDestP.hasComments()) implDestP.comment("Dsolvent");
  implDestP|Dsolvent;
}

/* DEFS: void write_output(const std::string &output_filename);
 */
void CProxy_MeshWorkingCopy::write_output(const std::string &output_filename, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const std::string &output_filename
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)output_filename;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)output_filename;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_MeshWorkingCopy::__idx_write_output_marshall6,0);
}
 int CkIndex_MeshWorkingCopy::__idx_write_output_marshall6=0;
void CkIndex_MeshWorkingCopy::_call_write_output_marshall6(void* impl_msg,MeshWorkingCopy * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::string &output_filename*/
  PUP::fromMem implP(impl_buf);
  std::string output_filename; implP|output_filename;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->write_output(output_filename);
}
int CkIndex_MeshWorkingCopy::_callmarshall_write_output_marshall6(char* impl_buf,MeshWorkingCopy * impl_obj) {
  /*Unmarshall pup'd fields: const std::string &output_filename*/
  PUP::fromMem implP(impl_buf);
  std::string output_filename; implP|output_filename;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->write_output(output_filename);
  return implP.size();
}
void CkIndex_MeshWorkingCopy::_marshallmessagepup_write_output_marshall6(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::string &output_filename*/
  PUP::fromMem implP(impl_buf);
  std::string output_filename; implP|output_filename;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("output_filename");
  implDestP|output_filename;
}

/* DEFS: MeshWorkingCopy(CkMigrateMessage* impl_msg);
 */

/* DEFS: MeshWorkingCopy(void);
 */

/* DEFS: void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
 */
void CProxySection_MeshWorkingCopy::init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(Quaternion &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)impl_noname_2;
    implP|impl_noname_3;
    implP|impl_noname_4;
    implP|impl_noname_5;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(Quaternion &)impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)impl_noname_2;
    implP|impl_noname_3;
    implP|impl_noname_4;
    implP|impl_noname_5;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_init_marshall2,0);
}

/* DEFS: void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
 */
void CProxySection_MeshWorkingCopy::calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(vanilla_fmm::CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(vanilla_fmm::CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_calculate_rhs_marshall3,0);
}

/* DEFS: void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
 */
void CProxySection_MeshWorkingCopy::process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_process_returned_eval_pts_Eval_Message,0);
}

/* DEFS: void calculate_energy(double kappa, double Dsolvent);
 */
void CProxySection_MeshWorkingCopy::calculate_energy(double kappa, double Dsolvent, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa, double Dsolvent
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    implP|Dsolvent;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
    implP|Dsolvent;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_calculate_energy_marshall5,0);
}

/* DEFS: void write_output(const std::string &output_filename);
 */
void CProxySection_MeshWorkingCopy::write_output(const std::string &output_filename, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const std::string &output_filename
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)output_filename;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)output_filename;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_MeshWorkingCopy::__idx_write_output_marshall6,0);
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_MeshWorkingCopy::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
// REG: MeshWorkingCopy(CkMigrateMessage* impl_msg);
  __idx_MeshWorkingCopy_CkMigrateMessage = CkRegisterEp("MeshWorkingCopy(CkMigrateMessage* impl_msg)",
     (CkCallFnPtr)_call_MeshWorkingCopy_CkMigrateMessage, 0, __idx, 0);
  CkRegisterMigCtor(__idx, __idx_MeshWorkingCopy_CkMigrateMessage);

// REG: MeshWorkingCopy(void);
  __idx_MeshWorkingCopy_void = CkRegisterEp("MeshWorkingCopy(void)",
     (CkCallFnPtr)_call_MeshWorkingCopy_void, 0, __idx, 0);
  CkRegisterDefaultCtor(__idx, __idx_MeshWorkingCopy_void);

// REG: void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
  __idx_init_marshall2 = CkRegisterEp("init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5)",
     (CkCallFnPtr)_call_init_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_init_marshall2,(CkMarshallUnpackFn)_callmarshall_init_marshall2);
  CkRegisterMessagePupFn(__idx_init_marshall2,(CkMessagePupFn)_marshallmessagepup_init_marshall2);

// REG: void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
  __idx_calculate_rhs_marshall3 = CkRegisterEp("calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb)",
     (CkCallFnPtr)_call_calculate_rhs_marshall3, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_calculate_rhs_marshall3,(CkMarshallUnpackFn)_callmarshall_calculate_rhs_marshall3);
  CkRegisterMessagePupFn(__idx_calculate_rhs_marshall3,(CkMessagePupFn)_marshallmessagepup_calculate_rhs_marshall3);

// REG: void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
  __idx_process_returned_eval_pts_Eval_Message = CkRegisterEp("process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg)",
     (CkCallFnPtr)_call_process_returned_eval_pts_Eval_Message, vanilla_fmm::CMessage_Eval_Message::__idx, __idx, 0);

// REG: void calculate_energy(double kappa, double Dsolvent);
  __idx_calculate_energy_marshall5 = CkRegisterEp("calculate_energy(double kappa, double Dsolvent)",
     (CkCallFnPtr)_call_calculate_energy_marshall5, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_calculate_energy_marshall5,(CkMarshallUnpackFn)_callmarshall_calculate_energy_marshall5);
  CkRegisterMessagePupFn(__idx_calculate_energy_marshall5,(CkMessagePupFn)_marshallmessagepup_calculate_energy_marshall5);

// REG: void write_output(const std::string &output_filename);
  __idx_write_output_marshall6 = CkRegisterEp("write_output(const std::string &output_filename)",
     (CkCallFnPtr)_call_write_output_marshall6, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_write_output_marshall6,(CkMarshallUnpackFn)_callmarshall_write_output_marshall6);
  CkRegisterMessagePupFn(__idx_write_output_marshall6,(CkMessagePupFn)_marshallmessagepup_write_output_marshall6);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registermesh_working_copy(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: array MeshWorkingCopy: ArrayElement{
MeshWorkingCopy(CkMigrateMessage* impl_msg);
MeshWorkingCopy(void);
void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
void calculate_energy(double kappa, double Dsolvent);
void write_output(const std::string &output_filename);
};
*/
  CkIndex_MeshWorkingCopy::__register("MeshWorkingCopy", sizeof(MeshWorkingCopy));

}
#endif
