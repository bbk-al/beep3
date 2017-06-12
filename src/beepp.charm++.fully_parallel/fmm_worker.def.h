


namespace beepp {
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

/* DEFS: template < int NTERMS, int NLAMBS, int NWAVES > array FMMWorkerT: ArrayElement{
FMMWorkerT(CkMigrateMessage* impl_msg);
FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
void form_multipoles(void);
void make_plane_waves(int dummy);
void pass_multipole_upwards(void);
void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
void evaluate(void);
void debug_chk(void);
void check_waves_complete(int impl_noname_2);
void request_collected_waves(void);
void receive_incoming_wave(CkMessage* impl_msg);
};
 */
#ifdef CK_TEMPLATES_ONLY
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx=0;
#endif
#ifdef CK_TEMPLATES_ONLY
/* DEFS: FMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DEFS: FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::insert(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length, int onPE, const CkEntryOptions *impl_e_opts)
{ 
  //Marshall: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &)fmm_globals_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &)fmm_tree_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FH_Values_NodeGroup &)fh_proxy;
    implP|edge_length;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &)fmm_globals_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &)fmm_tree_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FH_Values_NodeGroup &)fh_proxy;
    implP|edge_length;
  }
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_FMMWorkerT_marshall1,onPE);
}

/* DEFS: void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double impl_noname_0, double impl_noname_1, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    implP|impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
    implP|impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_solve_marshall2,0);
}

/* DEFS: void form_multipoles(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::form_multipoles(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_form_multipoles_void,0);
}

/* DEFS: void make_plane_waves(int dummy);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts) 
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
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_make_plane_waves_marshall4,0+CK_MSG_INLINE);
}

/* DEFS: void pass_multipole_upwards(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::pass_multipole_upwards(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_pass_multipole_upwards_void,0);
}

/* DEFS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_inherit_lpole_from_parent_FMM_Msg,0+CK_MSG_INLINE);
}

/* DEFS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_multipole_contribution_from_child_FMM_Msg,0+CK_MSG_EXPEDITED);
}

/* DEFS: void evaluate(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::evaluate(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_evaluate_void,0);
}

/* DEFS: void debug_chk(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::debug_chk(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_debug_chk_void,0);
}

/* DEFS: void check_waves_complete(int impl_noname_2);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::check_waves_complete(int impl_noname_2, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_2
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_2;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_2;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_check_waves_complete_marshall10,0);
}

/* DEFS: void request_collected_waves(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::request_collected_waves(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_request_collected_waves_void,0);
}

/* DEFS: void receive_incoming_wave(CkMessage* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_incoming_wave(CkMessage* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_incoming_wave_CkMessage,0+CK_MSG_INLINE);
}

/* DEFS: FMMWorkerT(CkMigrateMessage* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_FMMWorkerT_CkMigrateMessage=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_FMMWorkerT_CkMigrateMessage(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  new (impl_obj) FMMWorkerT < NTERMS, NLAMBS, NWAVES > ((CkMigrateMessage*)impl_msg);
}

/* DEFS: FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  CkArrayID CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::ckNew(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{ 
  //Marshall: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &)fmm_globals_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &)fmm_tree_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FH_Values_NodeGroup &)fh_proxy;
    implP|edge_length;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &)fmm_globals_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &)fmm_tree_proxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FH_Values_NodeGroup &)fh_proxy;
    implP|edge_length;
  }
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_FMMWorkerT_marshall1,opts);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_FMMWorkerT_marshall1=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_FMMWorkerT_marshall1(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > fmm_globals_proxy; implP|fmm_globals_proxy;
  CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > fmm_tree_proxy; implP|fmm_tree_proxy;
  CProxy_FH_Values_NodeGroup fh_proxy; implP|fh_proxy;
  double edge_length; implP|edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) FMMWorkerT < NTERMS, NLAMBS, NWAVES > (fmm_globals_proxy, fmm_tree_proxy, fh_proxy, edge_length);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callmarshall_FMMWorkerT_marshall1(char* impl_buf,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj) {
  /*Unmarshall pup'd fields: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > fmm_globals_proxy; implP|fmm_globals_proxy;
  CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > fmm_tree_proxy; implP|fmm_tree_proxy;
  CProxy_FH_Values_NodeGroup fh_proxy; implP|fh_proxy;
  double edge_length; implP|edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) FMMWorkerT < NTERMS, NLAMBS, NWAVES > (fmm_globals_proxy, fmm_tree_proxy, fh_proxy, edge_length);
  return implP.size();
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_FMMWorkerT_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > fmm_globals_proxy; implP|fmm_globals_proxy;
  CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > fmm_tree_proxy; implP|fmm_tree_proxy;
  CProxy_FH_Values_NodeGroup fh_proxy; implP|fh_proxy;
  double edge_length; implP|edge_length;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("fmm_globals_proxy");
  implDestP|fmm_globals_proxy;
  if (implDestP.hasComments()) implDestP.comment("fmm_tree_proxy");
  implDestP|fmm_tree_proxy;
  if (implDestP.hasComments()) implDestP.comment("fh_proxy");
  implDestP|fh_proxy;
  if (implDestP.hasComments()) implDestP.comment("edge_length");
  implDestP|edge_length;
}

/* DEFS: void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double impl_noname_0, double impl_noname_1, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    implP|impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
    implP|impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_solve_marshall2,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_solve_marshall2=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_solve_marshall2(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double impl_noname_0, double impl_noname_1, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  double impl_noname_0; implP|impl_noname_0;
  double impl_noname_1; implP|impl_noname_1;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->solve(impl_noname_0, impl_noname_1, cb);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callmarshall_solve_marshall2(char* impl_buf,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj) {
  /*Unmarshall pup'd fields: double impl_noname_0, double impl_noname_1, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  double impl_noname_0; implP|impl_noname_0;
  double impl_noname_1; implP|impl_noname_1;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->solve(impl_noname_0, impl_noname_1, cb);
  return implP.size();
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_solve_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double impl_noname_0, double impl_noname_1, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  double impl_noname_0; implP|impl_noname_0;
  double impl_noname_1; implP|impl_noname_1;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_1");
  implDestP|impl_noname_1;
  if (implDestP.hasComments()) implDestP.comment("cb");
  implDestP|cb;
}

/* DEFS: void form_multipoles(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::form_multipoles(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_form_multipoles_void,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_form_multipoles_void=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_form_multipoles_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->form_multipoles();
}

/* DEFS: void make_plane_waves(int dummy);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts) 
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
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_make_plane_waves_marshall4,0+CK_MSG_INLINE);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_make_plane_waves_marshall4=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_make_plane_waves_marshall4(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
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
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callmarshall_make_plane_waves_marshall4(char* impl_buf,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj) {
  /*Unmarshall pup'd fields: int dummy*/
  PUP::fromMem implP(impl_buf);
  int dummy; implP|dummy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->make_plane_waves(dummy);
  return implP.size();
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_make_plane_waves_marshall4(PUP::er &implDestP,void *impl_msg) {
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

/* DEFS: void pass_multipole_upwards(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::pass_multipole_upwards(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_pass_multipole_upwards_void,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_pass_multipole_upwards_void=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_pass_multipole_upwards_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->pass_multipole_upwards();
}

/* DEFS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_inherit_lpole_from_parent_FMM_Msg,0+CK_MSG_INLINE);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_inherit_lpole_from_parent_FMM_Msg=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_inherit_lpole_from_parent_FMM_Msg(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  impl_obj->inherit_lpole_from_parent((FMM_Msg<MultipoleHolderT<NTERMS > >*)impl_msg);
}

/* DEFS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_multipole_contribution_from_child_FMM_Msg,0+CK_MSG_EXPEDITED);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_multipole_contribution_from_child_FMM_Msg=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_receive_multipole_contribution_from_child_FMM_Msg(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  impl_obj->receive_multipole_contribution_from_child((FMM_Msg<MultipoleHolderT<NTERMS > >*)impl_msg);
}

/* DEFS: void evaluate(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::evaluate(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_evaluate_void,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_evaluate_void=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_evaluate_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->evaluate();
}

/* DEFS: void debug_chk(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::debug_chk(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_debug_chk_void,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_debug_chk_void=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_debug_chk_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->debug_chk();
}

/* DEFS: void check_waves_complete(int impl_noname_2);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::check_waves_complete(int impl_noname_2, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_2
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_2;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_2;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_check_waves_complete_marshall10,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_check_waves_complete_marshall10=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_check_waves_complete_marshall10(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int impl_noname_2*/
  PUP::fromMem implP(impl_buf);
  int impl_noname_2; implP|impl_noname_2;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->check_waves_complete(impl_noname_2);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_callmarshall_check_waves_complete_marshall10(char* impl_buf,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj) {
  /*Unmarshall pup'd fields: int impl_noname_2*/
  PUP::fromMem implP(impl_buf);
  int impl_noname_2; implP|impl_noname_2;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->check_waves_complete(impl_noname_2);
  return implP.size();
}
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_marshallmessagepup_check_waves_complete_marshall10(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: int impl_noname_2*/
  PUP::fromMem implP(impl_buf);
  int impl_noname_2; implP|impl_noname_2;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_2");
  implDestP|impl_noname_2;
}

/* DEFS: void request_collected_waves(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::request_collected_waves(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_request_collected_waves_void,0);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_request_collected_waves_void=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_request_collected_waves_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->request_collected_waves();
}

/* DEFS: void receive_incoming_wave(CkMessage* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_incoming_wave(CkMessage* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_incoming_wave_CkMessage,0+CK_MSG_INLINE);
}
template < int NTERMS, int NLAMBS, int NWAVES >  int CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_incoming_wave_CkMessage=0;
template < int NTERMS, int NLAMBS, int NWAVES >  void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::_call_receive_incoming_wave_CkMessage(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES >  * impl_obj)
{
  impl_obj->receive_incoming_wave((CkMessage*)impl_msg);
}

/* DEFS: FMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DEFS: FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
 */

/* DEFS: void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double impl_noname_0, double impl_noname_1, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_0;
    implP|impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_0;
    implP|impl_noname_1;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_solve_marshall2,0);
}

/* DEFS: void form_multipoles(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::form_multipoles(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_form_multipoles_void,0);
}

/* DEFS: void make_plane_waves(int dummy);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts) 
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
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_make_plane_waves_marshall4,0+CK_MSG_INLINE);
}

/* DEFS: void pass_multipole_upwards(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::pass_multipole_upwards(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_pass_multipole_upwards_void,0);
}

/* DEFS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_inherit_lpole_from_parent_FMM_Msg,0+CK_MSG_INLINE);
}

/* DEFS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_multipole_contribution_from_child_FMM_Msg,0+CK_MSG_EXPEDITED);
}

/* DEFS: void evaluate(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::evaluate(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_evaluate_void,0);
}

/* DEFS: void debug_chk(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::debug_chk(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_debug_chk_void,0);
}

/* DEFS: void check_waves_complete(int impl_noname_2);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::check_waves_complete(int impl_noname_2, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int impl_noname_2
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_2;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_2;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_check_waves_complete_marshall10,0);
}

/* DEFS: void request_collected_waves(void);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::request_collected_waves(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_request_collected_waves_void,0);
}

/* DEFS: void receive_incoming_wave(CkMessage* impl_msg);
 */
template < int NTERMS, int NLAMBS, int NWAVES >  void CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::receive_incoming_wave(CkMessage* impl_msg) 
{
  ckCheck();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__idx_receive_incoming_wave_CkMessage,0+CK_MSG_INLINE);
}

#endif /*CK_TEMPLATES_ONLY*/
#ifdef CK_TEMPLATES_ONLY
template < int NTERMS, int NLAMBS, int NWAVES > void CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES > ::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
// REG: FMMWorkerT(CkMigrateMessage* impl_msg);
  __idx_FMMWorkerT_CkMigrateMessage = CkRegisterEp("FMMWorkerT(CkMigrateMessage* impl_msg)",
     (CkCallFnPtr)_call_FMMWorkerT_CkMigrateMessage, 0, __idx, 0);
  CkRegisterMigCtor(__idx, __idx_FMMWorkerT_CkMigrateMessage);

// REG: FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
  __idx_FMMWorkerT_marshall1 = CkRegisterEp("FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length)",
     (CkCallFnPtr)_call_FMMWorkerT_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_FMMWorkerT_marshall1,(CkMarshallUnpackFn)_callmarshall_FMMWorkerT_marshall1);
  CkRegisterMessagePupFn(__idx_FMMWorkerT_marshall1,(CkMessagePupFn)_marshallmessagepup_FMMWorkerT_marshall1);

// REG: void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
  __idx_solve_marshall2 = CkRegisterEp("solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb)",
     (CkCallFnPtr)_call_solve_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_solve_marshall2,(CkMarshallUnpackFn)_callmarshall_solve_marshall2);
  CkRegisterMessagePupFn(__idx_solve_marshall2,(CkMessagePupFn)_marshallmessagepup_solve_marshall2);

// REG: void form_multipoles(void);
  __idx_form_multipoles_void = CkRegisterEp("form_multipoles(void)",
     (CkCallFnPtr)_call_form_multipoles_void, 0, __idx, 0);

// REG: void make_plane_waves(int dummy);
  __idx_make_plane_waves_marshall4 = CkRegisterEp("make_plane_waves(int dummy)",
     (CkCallFnPtr)_call_make_plane_waves_marshall4, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_make_plane_waves_marshall4,(CkMarshallUnpackFn)_callmarshall_make_plane_waves_marshall4);
  CkRegisterMessagePupFn(__idx_make_plane_waves_marshall4,(CkMessagePupFn)_marshallmessagepup_make_plane_waves_marshall4);

// REG: void pass_multipole_upwards(void);
  __idx_pass_multipole_upwards_void = CkRegisterEp("pass_multipole_upwards(void)",
     (CkCallFnPtr)_call_pass_multipole_upwards_void, 0, __idx, 0);

// REG: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
  __idx_inherit_lpole_from_parent_FMM_Msg = CkRegisterEp("inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg)",
     (CkCallFnPtr)_call_inherit_lpole_from_parent_FMM_Msg, CMessage_FMM_Msg<MultipoleHolderT<NTERMS > >::__idx, __idx, 0);

// REG: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
  __idx_receive_multipole_contribution_from_child_FMM_Msg = CkRegisterEp("receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg)",
     (CkCallFnPtr)_call_receive_multipole_contribution_from_child_FMM_Msg, CMessage_FMM_Msg<MultipoleHolderT<NTERMS > >::__idx, __idx, 0);

// REG: void evaluate(void);
  __idx_evaluate_void = CkRegisterEp("evaluate(void)",
     (CkCallFnPtr)_call_evaluate_void, 0, __idx, 0);

// REG: void debug_chk(void);
  __idx_debug_chk_void = CkRegisterEp("debug_chk(void)",
     (CkCallFnPtr)_call_debug_chk_void, 0, __idx, 0);

// REG: void check_waves_complete(int impl_noname_2);
  __idx_check_waves_complete_marshall10 = CkRegisterEp("check_waves_complete(int impl_noname_2)",
     (CkCallFnPtr)_call_check_waves_complete_marshall10, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_check_waves_complete_marshall10,(CkMarshallUnpackFn)_callmarshall_check_waves_complete_marshall10);
  CkRegisterMessagePupFn(__idx_check_waves_complete_marshall10,(CkMessagePupFn)_marshallmessagepup_check_waves_complete_marshall10);

// REG: void request_collected_waves(void);
  __idx_request_collected_waves_void = CkRegisterEp("request_collected_waves(void)",
     (CkCallFnPtr)_call_request_collected_waves_void, 0, __idx, 0);

// REG: void receive_incoming_wave(CkMessage* impl_msg);
  __idx_receive_incoming_wave_CkMessage = CkRegisterEp("receive_incoming_wave(CkMessage* impl_msg)",
     (CkCallFnPtr)_call_receive_incoming_wave_CkMessage, CMessage_CkMessage::__idx, __idx, 0);

}
#endif

/* DEFS: array FMMWorkerT<9,9,67 >: ArrayElement;
 */

/* DEFS: message FMM_Msg<MultipoleHolderT<9 > >;
 */
#ifndef CK_TEMPLATES_ONLY
#endif

/* DEFS: message FMM_Msg<MultiHolder<6,PlaneWaveHolderT<9,67 > > >;
 */
#ifndef CK_TEMPLATES_ONLY
#endif

} // namespace beepp

#ifndef CK_TEMPLATES_ONLY
void _registerfmm_worker(void)
{
  static int _done = 0; if(_done) return; _done = 1;
      _registerfh_values_nodegroup();



using namespace beepp;


/* REG: array FMMWorkerT<9,9,67 >: ArrayElement;
*/
  CkIndex_FMMWorkerT<9,9,67 >::__register("FMMWorkerT<9,9,67 >", sizeof(FMMWorkerT<9,9,67 >));

/* REG: message FMM_Msg<MultipoleHolderT<9 > >;
*/
CMessage_FMM_Msg<MultipoleHolderT<9 > >::__register("FMM_Msg<MultipoleHolderT<9 > >", sizeof(FMM_Msg<MultipoleHolderT<9 > >),(CkPackFnPtr) FMM_Msg<MultipoleHolderT<9 > >::pack,(CkUnpackFnPtr) FMM_Msg<MultipoleHolderT<9 > >::unpack);

/* REG: message FMM_Msg<MultiHolder<6,PlaneWaveHolderT<9,67 > > >;
*/
CMessage_FMM_Msg<MultiHolder<6,PlaneWaveHolderT<9,67 > > >::__register("FMM_Msg<MultiHolder<6,PlaneWaveHolderT<9,67 > > >", sizeof(FMM_Msg<MultiHolder<6,PlaneWaveHolderT<9,67 > > >),(CkPackFnPtr) FMM_Msg<MultiHolder<6,PlaneWaveHolderT<9,67 > > >::pack,(CkUnpackFnPtr) FMM_Msg<MultiHolder<6,PlaneWaveHolderT<9,67 > > >::unpack);


}
#endif
