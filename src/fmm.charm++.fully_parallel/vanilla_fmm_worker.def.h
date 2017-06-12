
namespace vanilla_fmm {
/* DEFS: template < class HType > message FMM_Msg;
 */
#ifdef CK_TEMPLATES_ONLY
template < class HType >  void *CMessage_FMM_Msg < HType > ::operator new(size_t s){
  return FMM_Msg < HType > ::alloc(__idx, s, 0, 0);
}
template < class HType >  void *CMessage_FMM_Msg < HType > ::operator new(size_t s, int* sz){
  return FMM_Msg < HType > ::alloc(__idx, s, sz, 0);
}
template < class HType >  void *CMessage_FMM_Msg < HType > ::operator new(size_t s, int* sz,const int pb){
  return FMM_Msg < HType > ::alloc(__idx, s, sz, pb);
}
template < class HType >  void *CMessage_FMM_Msg < HType > ::operator new(size_t s, const int p) {
  return FMM_Msg < HType > ::alloc(__idx, s, 0, p);
}
template < class HType >  void* CMessage_FMM_Msg < HType > ::alloc(int msgnum, size_t sz, int *sizes, int pb) {
  CkpvAccess(_offsets)[0] = ALIGN8(sz);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[0], pb);
}
template < class HType >  CMessage_FMM_Msg < HType > ::CMessage_FMM_Msg() {
FMM_Msg < HType >  *newmsg = (FMM_Msg < HType >  *)this;
}
template < class HType >  void CMessage_FMM_Msg < HType > ::dealloc(void *p) {
  CkFreeMsg(p);
}
template < class HType >  void* CMessage_FMM_Msg < HType > ::pack(FMM_Msg < HType >  *msg) {
  return (void *) msg;
}
template < class HType >  FMM_Msg < HType > * CMessage_FMM_Msg < HType > ::unpack(void* buf) {
  FMM_Msg < HType >  *msg = (FMM_Msg < HType >  *) buf;
  return msg;
}
template < class HType >  int CMessage_FMM_Msg < HType > ::__idx=0;
#endif

/* DEFS: template < int NTERMS, int NLAMBS, int NWAVES > array VanillaFMMWorkerT: ArrayElement{
VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
void evaluate(Eval_Message* impl_msg);
threaded void asynchronous_check(int impl_noname_0);
void check_waves_complete(int impl_noname_1);
void request_collected_waves(void);
void receive_incoming_wave(CkMessage* impl_msg);
void solve(double impl_noname_2, const CkCallback &impl_noname_3);
void make_plane_waves(int dummy);
};
 */
#ifdef CK_TEMPLATES_ONLY
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx=0;
#endif
#ifdef CK_TEMPLATES_ONLY
/* DEFS: VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DEFS: VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::insert(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length, int onPE, const CkEntryOptions *impl_e_opts)
{ 
  //Marshall: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &)fmm_globals_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelTree &)fmm_tree_proxy;
    implP|edge_length;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &)fmm_globals_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelTree &)fmm_tree_proxy;
    implP|edge_length;
  }
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_VanillaFMMWorkerT_marshall1,onPE);
}

/* DEFS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_inherit_lpole_from_parent_FMM_Msg,0+CK_MSG_INLINE);
}

/* DEFS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_multipole_contribution_from_child_FMM_Msg,0+CK_MSG_EXPEDITED);
}

/* DEFS: void evaluate(Eval_Message* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::evaluate(Eval_Message* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_evaluate_Eval_Message,0);
}

/* DEFS: threaded void asynchronous_check(int impl_noname_0);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::asynchronous_check(int impl_noname_0, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_asynchronous_check_marshall5,0);
}

/* DEFS: void check_waves_complete(int impl_noname_1);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::check_waves_complete(int impl_noname_1, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_1
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_1;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_1;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_check_waves_complete_marshall6,0);
}

/* DEFS: void request_collected_waves(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::request_collected_waves(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_request_collected_waves_void,0);
}

/* DEFS: void receive_incoming_wave(CkMessage* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_incoming_wave(CkMessage* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_incoming_wave_CkMessage,0+CK_MSG_INLINE);
}

/* DEFS: void solve(double impl_noname_2, const CkCallback &impl_noname_3);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::solve(double impl_noname_2, const CkCallback &impl_noname_3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double impl_noname_2, const CkCallback &impl_noname_3
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_2;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_3;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_2;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_3;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_solve_marshall9,0);
}

/* DEFS: void make_plane_waves(int dummy);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int dummy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|dummy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|dummy;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_make_plane_waves_marshall10,0);
}

/* DEFS: VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_VanillaFMMWorkerT_CkMigrateMessage=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_VanillaFMMWorkerT_CkMigrateMessage(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  new (impl_obj) VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ((CkMigrateMessage*)impl_msg);
}

/* DEFS: VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  CkArrayID CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::ckNew(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{ 
  //Marshall: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &)fmm_globals_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelTree &)fmm_tree_proxy;
    implP|edge_length;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &)fmm_globals_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelTree &)fmm_tree_proxy;
    implP|edge_length;
  }
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_VanillaFMMWorkerT_marshall1,opts);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_VanillaFMMWorkerT_marshall1=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_VanillaFMMWorkerT_marshall1(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > fmm_globals_proxy; implP|fmm_globals_proxy;
  CProxy_ParallelTree fmm_tree_proxy; implP|fmm_tree_proxy;
  double edge_length; implP|edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > (fmm_globals_proxy, fmm_tree_proxy, edge_length);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callmarshall_VanillaFMMWorkerT_marshall1(char* impl_buf,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj) {
  /*Unmarshall pup'd fields: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > fmm_globals_proxy; implP|fmm_globals_proxy;
  CProxy_ParallelTree fmm_tree_proxy; implP|fmm_tree_proxy;
  double edge_length; implP|edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > (fmm_globals_proxy, fmm_tree_proxy, edge_length);
  return implP.size();
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_VanillaFMMWorkerT_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > fmm_globals_proxy; implP|fmm_globals_proxy;
  CProxy_ParallelTree fmm_tree_proxy; implP|fmm_tree_proxy;
  double edge_length; implP|edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("fmm_globals_proxy");
  implDestP|fmm_globals_proxy;
  if (implDestP.hasComments()) implDestP.comment("fmm_tree_proxy");
  implDestP|fmm_tree_proxy;
  if (implDestP.hasComments()) implDestP.comment("edge_length");
  implDestP|edge_length;
}

/* DEFS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_inherit_lpole_from_parent_FMM_Msg,0+CK_MSG_INLINE);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_inherit_lpole_from_parent_FMM_Msg=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_inherit_lpole_from_parent_FMM_Msg(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  impl_obj->inherit_lpole_from_parent((FMM_Msg<MultipoleHolderT<NTERMS > >*)impl_msg);
}

/* DEFS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_multipole_contribution_from_child_FMM_Msg,0+CK_MSG_EXPEDITED);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_multipole_contribution_from_child_FMM_Msg=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_receive_multipole_contribution_from_child_FMM_Msg(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  impl_obj->receive_multipole_contribution_from_child((FMM_Msg<MultipoleHolderT<NTERMS > >*)impl_msg);
}

/* DEFS: void evaluate(Eval_Message* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::evaluate(Eval_Message* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_evaluate_Eval_Message,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_evaluate_Eval_Message=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_evaluate_Eval_Message(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  impl_obj->evaluate((Eval_Message*)impl_msg);
}

/* DEFS: threaded void asynchronous_check(int impl_noname_0);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::asynchronous_check(int impl_noname_0, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_asynchronous_check_marshall5,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_asynchronous_check_marshall5=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_asynchronous_check_marshall5(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CthThread tid = CthCreate((CthVoidFn)_callthr_asynchronous_check_marshall5, new CkThrCallArg(impl_msg,impl_obj), 0);
  ((Chare *)impl_obj)->CkAddThreadListeners(tid,impl_msg);
  CthAwaken(tid);
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callthr_asynchronous_check_marshall5(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  *impl_obj = (VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  *) impl_arg->obj;
  delete impl_arg;
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  int impl_noname_0; implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->asynchronous_check(impl_noname_0);
  delete impl_msg_typed;
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_asynchronous_check_marshall5(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  int impl_noname_0; implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
}

/* DEFS: void check_waves_complete(int impl_noname_1);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::check_waves_complete(int impl_noname_1, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_1
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_1;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_1;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_check_waves_complete_marshall6,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_check_waves_complete_marshall6=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_check_waves_complete_marshall6(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int impl_noname_1*/
  PUP::fromMem implP(impl_buf);
  int impl_noname_1; implP|impl_noname_1;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->check_waves_complete(impl_noname_1);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callmarshall_check_waves_complete_marshall6(char* impl_buf,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj) {
  /*Unmarshall pup'd fields: int impl_noname_1*/
  PUP::fromMem implP(impl_buf);
  int impl_noname_1; implP|impl_noname_1;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->check_waves_complete(impl_noname_1);
  return implP.size();
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_check_waves_complete_marshall6(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int impl_noname_1*/
  PUP::fromMem implP(impl_buf);
  int impl_noname_1; implP|impl_noname_1;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_1");
  implDestP|impl_noname_1;
}

/* DEFS: void request_collected_waves(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::request_collected_waves(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_request_collected_waves_void,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_request_collected_waves_void=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_request_collected_waves_void(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->request_collected_waves();
}

/* DEFS: void receive_incoming_wave(CkMessage* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_incoming_wave(CkMessage* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_incoming_wave_CkMessage,0+CK_MSG_INLINE);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_incoming_wave_CkMessage=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_receive_incoming_wave_CkMessage(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  impl_obj->receive_incoming_wave((CkMessage*)impl_msg);
}

/* DEFS: void solve(double impl_noname_2, const CkCallback &impl_noname_3);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::solve(double impl_noname_2, const CkCallback &impl_noname_3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double impl_noname_2, const CkCallback &impl_noname_3
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_2;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_3;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_2;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_3;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_solve_marshall9,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_solve_marshall9=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_solve_marshall9(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double impl_noname_2, const CkCallback &impl_noname_3*/
  PUP::fromMem implP(impl_buf);
  double impl_noname_2; implP|impl_noname_2;
  CkCallback impl_noname_3; implP|impl_noname_3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->solve(impl_noname_2, impl_noname_3);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callmarshall_solve_marshall9(char* impl_buf,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj) {
  /*Unmarshall pup'd fields: double impl_noname_2, const CkCallback &impl_noname_3*/
  PUP::fromMem implP(impl_buf);
  double impl_noname_2; implP|impl_noname_2;
  CkCallback impl_noname_3; implP|impl_noname_3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->solve(impl_noname_2, impl_noname_3);
  return implP.size();
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_solve_marshall9(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double impl_noname_2, const CkCallback &impl_noname_3*/
  PUP::fromMem implP(impl_buf);
  double impl_noname_2; implP|impl_noname_2;
  CkCallback impl_noname_3; implP|impl_noname_3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_2");
  implDestP|impl_noname_2;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_3");
  implDestP|impl_noname_3;
}

/* DEFS: void make_plane_waves(int dummy);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int dummy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|dummy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|dummy;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_make_plane_waves_marshall10,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_make_plane_waves_marshall10=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_make_plane_waves_marshall10(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int dummy*/
  PUP::fromMem implP(impl_buf);
  int dummy; implP|dummy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->make_plane_waves(dummy);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callmarshall_make_plane_waves_marshall10(char* impl_buf,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj) {
  /*Unmarshall pup'd fields: int dummy*/
  PUP::fromMem implP(impl_buf);
  int dummy; implP|dummy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->make_plane_waves(dummy);
  return implP.size();
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_make_plane_waves_marshall10(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int dummy*/
  PUP::fromMem implP(impl_buf);
  int dummy; implP|dummy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("dummy");
  implDestP|dummy;
}

/* DEFS: VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DEFS: VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
 */

/* DEFS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_inherit_lpole_from_parent_FMM_Msg,0+CK_MSG_INLINE);
}

/* DEFS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_multipole_contribution_from_child_FMM_Msg,0+CK_MSG_EXPEDITED);
}

/* DEFS: void evaluate(Eval_Message* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::evaluate(Eval_Message* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_evaluate_Eval_Message,0);
}

/* DEFS: threaded void asynchronous_check(int impl_noname_0);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::asynchronous_check(int impl_noname_0, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_asynchronous_check_marshall5,0);
}

/* DEFS: void check_waves_complete(int impl_noname_1);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::check_waves_complete(int impl_noname_1, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_1
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_1;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_1;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_check_waves_complete_marshall6,0);
}

/* DEFS: void request_collected_waves(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::request_collected_waves(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_request_collected_waves_void,0);
}

/* DEFS: void receive_incoming_wave(CkMessage* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_incoming_wave(CkMessage* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_incoming_wave_CkMessage,0+CK_MSG_INLINE);
}

/* DEFS: void solve(double impl_noname_2, const CkCallback &impl_noname_3);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::solve(double impl_noname_2, const CkCallback &impl_noname_3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double impl_noname_2, const CkCallback &impl_noname_3
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_2;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_3;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_2;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_3;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_solve_marshall9,0);
}

/* DEFS: void make_plane_waves(int dummy);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int dummy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|dummy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|dummy;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_make_plane_waves_marshall10,0);
}

#endif /*CK_TEMPLATES_ONLY*/
#ifdef CK_TEMPLATES_ONLY
template < int NTERMS, int NLAMBS, int NWAVES > void CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
// REG: VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
  __idx_VanillaFMMWorkerT_CkMigrateMessage = CkRegisterEp("VanillaFMMWorkerT(CkMigrateMessage* impl_msg)",
     (CkCallFnPtr)_call_VanillaFMMWorkerT_CkMigrateMessage, 0, __idx, 0);
  CkRegisterMigCtor(__idx, __idx_VanillaFMMWorkerT_CkMigrateMessage);

// REG: VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
  __idx_VanillaFMMWorkerT_marshall1 = CkRegisterEp("VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length)",
     (CkCallFnPtr)_call_VanillaFMMWorkerT_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_VanillaFMMWorkerT_marshall1,(CkMarshallUnpackFn)_callmarshall_VanillaFMMWorkerT_marshall1);
  CkRegisterMessagePupFn(__idx_VanillaFMMWorkerT_marshall1,(CkMessagePupFn)_marshallmessagepup_VanillaFMMWorkerT_marshall1);

// REG: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
  __idx_inherit_lpole_from_parent_FMM_Msg = CkRegisterEp("inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg)",
     (CkCallFnPtr)_call_inherit_lpole_from_parent_FMM_Msg, CMessage_FMM_Msg<MultipoleHolderT<NTERMS > >::__idx, __idx, 0);

// REG: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
  __idx_receive_multipole_contribution_from_child_FMM_Msg = CkRegisterEp("receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg)",
     (CkCallFnPtr)_call_receive_multipole_contribution_from_child_FMM_Msg, CMessage_FMM_Msg<MultipoleHolderT<NTERMS > >::__idx, __idx, 0);

// REG: void evaluate(Eval_Message* impl_msg);
  __idx_evaluate_Eval_Message = CkRegisterEp("evaluate(Eval_Message* impl_msg)",
     (CkCallFnPtr)_call_evaluate_Eval_Message, CMessage_Eval_Message::__idx, __idx, 0);

// REG: threaded void asynchronous_check(int impl_noname_0);
  __idx_asynchronous_check_marshall5 = CkRegisterEp("asynchronous_check(int impl_noname_0)",
     (CkCallFnPtr)_call_asynchronous_check_marshall5, CkMarshallMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(__idx_asynchronous_check_marshall5,(CkMessagePupFn)_marshallmessagepup_asynchronous_check_marshall5);

// REG: void check_waves_complete(int impl_noname_1);
  __idx_check_waves_complete_marshall6 = CkRegisterEp("check_waves_complete(int impl_noname_1)",
     (CkCallFnPtr)_call_check_waves_complete_marshall6, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_check_waves_complete_marshall6,(CkMarshallUnpackFn)_callmarshall_check_waves_complete_marshall6);
  CkRegisterMessagePupFn(__idx_check_waves_complete_marshall6,(CkMessagePupFn)_marshallmessagepup_check_waves_complete_marshall6);

// REG: void request_collected_waves(void);
  __idx_request_collected_waves_void = CkRegisterEp("request_collected_waves(void)",
     (CkCallFnPtr)_call_request_collected_waves_void, 0, __idx, 0);

// REG: void receive_incoming_wave(CkMessage* impl_msg);
  __idx_receive_incoming_wave_CkMessage = CkRegisterEp("receive_incoming_wave(CkMessage* impl_msg)",
     (CkCallFnPtr)_call_receive_incoming_wave_CkMessage, CMessage_CkMessage::__idx, __idx, 0);

// REG: void solve(double impl_noname_2, const CkCallback &impl_noname_3);
  __idx_solve_marshall9 = CkRegisterEp("solve(double impl_noname_2, const CkCallback &impl_noname_3)",
     (CkCallFnPtr)_call_solve_marshall9, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_solve_marshall9,(CkMarshallUnpackFn)_callmarshall_solve_marshall9);
  CkRegisterMessagePupFn(__idx_solve_marshall9,(CkMessagePupFn)_marshallmessagepup_solve_marshall9);

// REG: void make_plane_waves(int dummy);
  __idx_make_plane_waves_marshall10 = CkRegisterEp("make_plane_waves(int dummy)",
     (CkCallFnPtr)_call_make_plane_waves_marshall10, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_make_plane_waves_marshall10,(CkMarshallUnpackFn)_callmarshall_make_plane_waves_marshall10);
  CkRegisterMessagePupFn(__idx_make_plane_waves_marshall10,(CkMessagePupFn)_marshallmessagepup_make_plane_waves_marshall10);

}
#endif

/* DEFS: array VanillaFMMWorkerT<18,18,300 >: ArrayElement;
 */

/* DEFS: message FMM_Msg<MultipoleHolderT<18 > >;
 */
#ifndef CK_TEMPLATES_ONLY
#endif

/* DEFS: message FMM_Msg<MultiHolder<6,PlaneWaveHolderT<18,300 > > >;
 */
#ifndef CK_TEMPLATES_ONLY
#endif


} // namespace vanilla_fmm

#ifndef CK_TEMPLATES_ONLY
void _registervanilla_fmm_worker(void)
{
  static int _done = 0; if(_done) return; _done = 1;
      _registervanilla_fmm_evals();

using namespace vanilla_fmm;


/* REG: array VanillaFMMWorkerT<18,18,300 >: ArrayElement;
*/
  CkIndex_VanillaFMMWorkerT<18,18,300 >::__register("VanillaFMMWorkerT<18,18,300 >", sizeof(VanillaFMMWorkerT<18,18,300 >));

/* REG: message FMM_Msg<MultipoleHolderT<18 > >;
*/
CMessage_FMM_Msg<MultipoleHolderT<18 > >::__register("FMM_Msg<MultipoleHolderT<18 > >", sizeof(FMM_Msg<MultipoleHolderT<18 > >),(CkPackFnPtr) FMM_Msg<MultipoleHolderT<18 > >::pack,(CkUnpackFnPtr) FMM_Msg<MultipoleHolderT<18 > >::unpack);

/* REG: message FMM_Msg<MultiHolder<6,PlaneWaveHolderT<18,300 > > >;
*/
CMessage_FMM_Msg<MultiHolder<6,PlaneWaveHolderT<18,300 > > >::__register("FMM_Msg<MultiHolder<6,PlaneWaveHolderT<18,300 > > >", sizeof(FMM_Msg<MultiHolder<6,PlaneWaveHolderT<18,300 > > >),(CkPackFnPtr) FMM_Msg<MultiHolder<6,PlaneWaveHolderT<18,300 > > >::pack,(CkUnpackFnPtr) FMM_Msg<MultiHolder<6,PlaneWaveHolderT<18,300 > > >::unpack);



}
#endif
