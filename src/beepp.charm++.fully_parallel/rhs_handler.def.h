/* DEFS: chare RHS_Handler: Chare{
RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length);
void add_charges(const std::vector<Charge > &charges);
threaded void get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb);
void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_RHS_Handler::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length);
 */
CkChareID CProxy_RHS_Handler::ckNew(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
  }
  CkChareID impl_ret;
  CkCreateChare(CkIndex_RHS_Handler::__idx, CkIndex_RHS_Handler::__idx_RHS_Handler_marshall1, impl_msg, &impl_ret, impl_onPE);
  return impl_ret;
}
void CProxy_RHS_Handler::ckNew(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, CkChareID* pcid, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
  }
  CkCreateChare(CkIndex_RHS_Handler::__idx, CkIndex_RHS_Handler::__idx_RHS_Handler_marshall1, impl_msg, pcid, impl_onPE);
}
  CProxy_RHS_Handler::CProxy_RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
  }
  CkChareID impl_ret;
  CkCreateChare(CkIndex_RHS_Handler::__idx, CkIndex_RHS_Handler::__idx_RHS_Handler_marshall1, impl_msg, &impl_ret, impl_onPE);
  ckSetChareID(impl_ret);
}
 int CkIndex_RHS_Handler::__idx_RHS_Handler_marshall1=0;
void CkIndex_RHS_Handler::_call_RHS_Handler_marshall1(void* impl_msg,RHS_Handler * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length*/
  PUP::fromMem implP(impl_buf);
  unsigned int max_items_per_node; implP|max_items_per_node;
  Vector universe_centre; implP|universe_centre;
  double universe_edge_length; implP|universe_edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) RHS_Handler(max_items_per_node, universe_centre, universe_edge_length);
}
int CkIndex_RHS_Handler::_callmarshall_RHS_Handler_marshall1(char* impl_buf,RHS_Handler * impl_obj) {
  /*Unmarshall pup'd fields: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length*/
  PUP::fromMem implP(impl_buf);
  unsigned int max_items_per_node; implP|max_items_per_node;
  Vector universe_centre; implP|universe_centre;
  double universe_edge_length; implP|universe_edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) RHS_Handler(max_items_per_node, universe_centre, universe_edge_length);
  return implP.size();
}
void CkIndex_RHS_Handler::_marshallmessagepup_RHS_Handler_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length*/
  PUP::fromMem implP(impl_buf);
  unsigned int max_items_per_node; implP|max_items_per_node;
  Vector universe_centre; implP|universe_centre;
  double universe_edge_length; implP|universe_edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("max_items_per_node");
  implDestP|max_items_per_node;
  if (implDestP.hasComments()) implDestP.comment("universe_centre");
  implDestP|universe_centre;
  if (implDestP.hasComments()) implDestP.comment("universe_edge_length");
  implDestP|universe_edge_length;
}

/* DEFS: void add_charges(const std::vector<Charge > &charges);
 */
void CProxy_RHS_Handler::add_charges(const std::vector<Charge > &charges, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const std::vector<Charge > &charges
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<Charge > &)charges;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<Charge > &)charges;
  }
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_RHS_Handler::__idx_add_charges_marshall2, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_RHS_Handler::__idx_add_charges_marshall2, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_RHS_Handler::__idx_add_charges_marshall2, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_RHS_Handler::__idx_add_charges_marshall2=0;
void CkIndex_RHS_Handler::_call_add_charges_marshall2(void* impl_msg,RHS_Handler * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::vector<Charge > &charges*/
  PUP::fromMem implP(impl_buf);
  std::vector<Charge > charges; implP|charges;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->add_charges(charges);
}
int CkIndex_RHS_Handler::_callmarshall_add_charges_marshall2(char* impl_buf,RHS_Handler * impl_obj) {
  /*Unmarshall pup'd fields: const std::vector<Charge > &charges*/
  PUP::fromMem implP(impl_buf);
  std::vector<Charge > charges; implP|charges;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->add_charges(charges);
  return implP.size();
}
void CkIndex_RHS_Handler::_marshallmessagepup_add_charges_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::vector<Charge > &charges*/
  PUP::fromMem implP(impl_buf);
  std::vector<Charge > charges; implP|charges;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("charges");
  implDestP|charges;
}

/* DEFS: threaded void get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb);
 */
void CProxy_RHS_Handler::get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|Dsolvent;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|Dsolvent;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_RHS_Handler::__idx_get_rhs_marshall3, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_RHS_Handler::__idx_get_rhs_marshall3, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_RHS_Handler::__idx_get_rhs_marshall3, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_RHS_Handler::__idx_get_rhs_marshall3=0;
void CkIndex_RHS_Handler::_call_get_rhs_marshall3(void* impl_msg,RHS_Handler * impl_obj)
{
  CthThread tid = CthCreate((CthVoidFn)_callthr_get_rhs_marshall3, new CkThrCallArg(impl_msg,impl_obj), 0);
  ((Chare *)impl_obj)->CkAddThreadListeners(tid,impl_msg);
  CthAwaken(tid);
}
void CkIndex_RHS_Handler::_callthr_get_rhs_marshall3(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  RHS_Handler *impl_obj = (RHS_Handler *) impl_arg->obj;
  delete impl_arg;
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  double Dsolvent; implP|Dsolvent;
  unsigned int impl_noname_0; implP|impl_noname_0;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->get_rhs(Dsolvent, impl_noname_0, cb);
  delete impl_msg_typed;
}
void CkIndex_RHS_Handler::_marshallmessagepup_get_rhs_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  double Dsolvent; implP|Dsolvent;
  unsigned int impl_noname_0; implP|impl_noname_0;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("Dsolvent");
  implDestP|Dsolvent;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
  if (implDestP.hasComments()) implDestP.comment("cb");
  implDestP|cb;
}

/* DEFS: void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg);
 */
void CProxy_RHS_Handler::process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_RHS_Handler::__idx_process_rhs_results_from_working_mesh_Eval_Message, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_RHS_Handler::__idx_process_rhs_results_from_working_mesh_Eval_Message, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_RHS_Handler::__idx_process_rhs_results_from_working_mesh_Eval_Message, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_RHS_Handler::__idx_process_rhs_results_from_working_mesh_Eval_Message=0;
void CkIndex_RHS_Handler::_call_process_rhs_results_from_working_mesh_Eval_Message(void* impl_msg,RHS_Handler * impl_obj)
{
  impl_obj->process_rhs_results_from_working_mesh((vanilla_fmm::Eval_Message*)impl_msg);
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_RHS_Handler::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeChare);
  CkRegisterBase(__idx, CkIndex_Chare::__idx);
// REG: RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length);
  __idx_RHS_Handler_marshall1 = CkRegisterEp("RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length)",
     (CkCallFnPtr)_call_RHS_Handler_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_RHS_Handler_marshall1,(CkMarshallUnpackFn)_callmarshall_RHS_Handler_marshall1);
  CkRegisterMessagePupFn(__idx_RHS_Handler_marshall1,(CkMessagePupFn)_marshallmessagepup_RHS_Handler_marshall1);

// REG: void add_charges(const std::vector<Charge > &charges);
  __idx_add_charges_marshall2 = CkRegisterEp("add_charges(const std::vector<Charge > &charges)",
     (CkCallFnPtr)_call_add_charges_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_add_charges_marshall2,(CkMarshallUnpackFn)_callmarshall_add_charges_marshall2);
  CkRegisterMessagePupFn(__idx_add_charges_marshall2,(CkMessagePupFn)_marshallmessagepup_add_charges_marshall2);

// REG: threaded void get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb);
  __idx_get_rhs_marshall3 = CkRegisterEp("get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb)",
     (CkCallFnPtr)_call_get_rhs_marshall3, CkMarshallMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(__idx_get_rhs_marshall3,(CkMessagePupFn)_marshallmessagepup_get_rhs_marshall3);

// REG: void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg);
  __idx_process_rhs_results_from_working_mesh_Eval_Message = CkRegisterEp("process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg)",
     (CkCallFnPtr)_call_process_rhs_results_from_working_mesh_Eval_Message, vanilla_fmm::CMessage_Eval_Message::__idx, __idx, 0);

}
#endif

/* DEFS: message RHS_Message{
double data[];
}
;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_RHS_Message::operator new(size_t s){
  return RHS_Message::alloc(__idx, s, 0, 0);
}
void *CMessage_RHS_Message::operator new(size_t s, int* sz){
  return RHS_Message::alloc(__idx, s, sz, 0);
}
void *CMessage_RHS_Message::operator new(size_t s, int* sz,const int pb){
  return RHS_Message::alloc(__idx, s, sz, pb);
}
void *CMessage_RHS_Message::operator new(size_t s, int sz0) {
  int sizes[1];
  sizes[0] = sz0;
  return RHS_Message::alloc(__idx, s, sizes, 0);
}
void *CMessage_RHS_Message::operator new(size_t s, int sz0, const int p) {
  int sizes[1];
  sizes[0] = sz0;
  return RHS_Message::alloc(__idx, s, sizes, p);
}
void* CMessage_RHS_Message::alloc(int msgnum, size_t sz, int *sizes, int pb) {
  CkpvAccess(_offsets)[0] = ALIGN8(sz);
  if(sizes==0)
    CkpvAccess(_offsets)[1] = CkpvAccess(_offsets)[0];
  else
    CkpvAccess(_offsets)[1] = CkpvAccess(_offsets)[0] + ALIGN8(sizeof(double)*sizes[0]);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[1], pb);
}
CMessage_RHS_Message::CMessage_RHS_Message() {
RHS_Message *newmsg = (RHS_Message *)this;
  newmsg->data = (double *) ((char *)newmsg + CkpvAccess(_offsets)[0]);
}
void CMessage_RHS_Message::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_RHS_Message::pack(RHS_Message *msg) {
  msg->data = (double *) ((char *)msg->data - (char *)msg);
  return (void *) msg;
}
RHS_Message* CMessage_RHS_Message::unpack(void* buf) {
  RHS_Message *msg = (RHS_Message *) buf;
  msg->data = (double *) ((size_t)msg->data + (char *)msg);
  return msg;
}
int CMessage_RHS_Message::__idx=0;
#endif

#ifndef CK_TEMPLATES_ONLY
void _registerrhs_handler(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: chare RHS_Handler: Chare{
RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length);
void add_charges(const std::vector<Charge > &charges);
threaded void get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb);
void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg);
};
*/
  CkIndex_RHS_Handler::__register("RHS_Handler", sizeof(RHS_Handler));

/* REG: message RHS_Message{
double data[];
}
;
*/
CMessage_RHS_Message::__register("RHS_Message", sizeof(RHS_Message),(CkPackFnPtr) RHS_Message::pack,(CkUnpackFnPtr) RHS_Message::unpack);

}
#endif
