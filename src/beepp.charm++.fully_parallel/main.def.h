













/* DEFS: nodegroup ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >: NodeGroup;
 */

/* DEFS: nodegroup ParallelFMMOctree<CharmNodePatch,beepp::CProxy_FMMWorker >: NodeGroup;
 */

/* DEFS: readonly CProxy_Main MainProxy;
 */
extern CProxy_Main MainProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MainProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|MainProxy;
}
#endif

/* DEFS: readonly CProxy_MeshWorkingCopy MeshWorkingCopyProxy;
 */
extern CProxy_MeshWorkingCopy MeshWorkingCopyProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MeshWorkingCopyProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|MeshWorkingCopyProxy;
}
#endif

/* DEFS: readonly CProxy_MeshLibrary MeshLibraryProxy;
 */
extern CProxy_MeshLibrary MeshLibraryProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MeshLibraryProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|MeshLibraryProxy;
}
#endif

/* DEFS: readonly beepp::CProxy_FMM_Globals_NodeGroup FMM_Globals_Proxy;
 */
extern beepp::CProxy_FMM_Globals_NodeGroup FMM_Globals_Proxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_FMM_Globals_Proxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|FMM_Globals_Proxy;
}
#endif

/* DEFS: readonly beepp::CProxy_ParallelTree ParallelFMMOctreeProxy;
 */
extern beepp::CProxy_ParallelTree ParallelFMMOctreeProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_ParallelFMMOctreeProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|ParallelFMMOctreeProxy;
}
#endif

/* DEFS: readonly CProxy_FMMWorker FMMWorkerProxy;
 */
extern CProxy_FMMWorker FMMWorkerProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_FMMWorkerProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|FMMWorkerProxy;
}
#endif

/* DEFS: readonly CProxy_ExplicitWorker ExplicitWorkerProxy;
 */
extern CProxy_ExplicitWorker ExplicitWorkerProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_ExplicitWorkerProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|ExplicitWorkerProxy;
}
#endif

/* DEFS: readonly CProxy_RHS_Handler RHS_HandlerProxy;
 */
extern CProxy_RHS_Handler RHS_HandlerProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_RHS_HandlerProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|RHS_HandlerProxy;
}
#endif

/* DEFS: readonly CProxy_GMRES GMRESProxy;
 */
extern CProxy_GMRES GMRESProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_GMRESProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|GMRESProxy;
}
#endif

/* DEFS: readonly CProxy_IterationHandler IterationHandlerProxy;
 */
extern CProxy_IterationHandler IterationHandlerProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_IterationHandlerProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|IterationHandlerProxy;
}
#endif

/* DEFS: readonly CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;
 */
extern CProxy_FH_Values_NodeGroup FH_Values_NodeGroupProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_FH_Values_NodeGroupProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|FH_Values_NodeGroupProxy;
}
#endif

/* DEFS: readonly CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
 */
extern CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_OpenCL_NodeGroupProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|OpenCL_NodeGroupProxy;
}
#endif

/* DEFS: readonly ComlibInstanceHandle streaming_strat;
 */
extern ComlibInstanceHandle streaming_strat;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_streaming_strat(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|streaming_strat;
}
#endif

/* DEFS: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
threaded void create_workers(const CkCallback &impl_noname_0);
threaded void get_rhs(void);
void quiescenceHandler(void);
  initcall void initManualLB(void);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Main::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Main(CkArgMsg* impl_msg);
 */
CkChareID CProxy_Main::ckNew(CkArgMsg* impl_msg, int impl_onPE)
{
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::__idx_Main_CkArgMsg, impl_msg, &impl_ret, impl_onPE);
  return impl_ret;
}
void CProxy_Main::ckNew(CkArgMsg* impl_msg, CkChareID* pcid, int impl_onPE)
{
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::__idx_Main_CkArgMsg, impl_msg, pcid, impl_onPE);
}
  CProxy_Main::CProxy_Main(CkArgMsg* impl_msg, int impl_onPE)
{
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::__idx_Main_CkArgMsg, impl_msg, &impl_ret, impl_onPE);
  ckSetChareID(impl_ret);
}
 int CkIndex_Main::__idx_Main_CkArgMsg=0;
void CkIndex_Main::_call_Main_CkArgMsg(void* impl_msg,Main * impl_obj)
{
  new (impl_obj) Main((CkArgMsg*)impl_msg);
}

/* DEFS: threaded void create_workers(const CkCallback &impl_noname_0);
 */
void CProxy_Main::create_workers(const CkCallback &impl_noname_0, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkCallback &impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_0;
  }
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_create_workers_marshall2, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::__idx_create_workers_marshall2, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_create_workers_marshall2, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_Main::__idx_create_workers_marshall2=0;
void CkIndex_Main::_call_create_workers_marshall2(void* impl_msg,Main * impl_obj)
{
  CthThread tid = CthCreate((CthVoidFn)_callthr_create_workers_marshall2, new CkThrCallArg(impl_msg,impl_obj), 0);
  ((Chare *)impl_obj)->CkAddThreadListeners(tid,impl_msg);
  CthAwaken(tid);
}
void CkIndex_Main::_callthr_create_workers_marshall2(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  Main *impl_obj = (Main *) impl_arg->obj;
  delete impl_arg;
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  CkCallback impl_noname_0; implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->create_workers(impl_noname_0);
  delete impl_msg_typed;
}
void CkIndex_Main::_marshallmessagepup_create_workers_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  CkCallback impl_noname_0; implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
}

/* DEFS: threaded void get_rhs(void);
 */
void CProxy_Main::get_rhs(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_get_rhs_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::__idx_get_rhs_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_get_rhs_void, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_Main::__idx_get_rhs_void=0;
void CkIndex_Main::_call_get_rhs_void(void* impl_msg,Main * impl_obj)
{
  CthThread tid = CthCreate((CthVoidFn)_callthr_get_rhs_void, new CkThrCallArg(impl_msg,impl_obj), 0);
  ((Chare *)impl_obj)->CkAddThreadListeners(tid,impl_msg);
  CthAwaken(tid);
}
void CkIndex_Main::_callthr_get_rhs_void(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  Main *impl_obj = (Main *) impl_arg->obj;
  delete impl_arg;
  CkFreeSysMsg(impl_msg);
  impl_obj->get_rhs();
}

/* DEFS: void quiescenceHandler(void);
 */
void CProxy_Main::quiescenceHandler(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_quiescenceHandler_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::__idx_quiescenceHandler_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_quiescenceHandler_void, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_Main::__idx_quiescenceHandler_void=0;
void CkIndex_Main::_call_quiescenceHandler_void(void* impl_msg,Main * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->quiescenceHandler();
}


#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Main::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeMainChare);
  CkRegisterBase(__idx, CkIndex_Chare::__idx);
// REG: Main(CkArgMsg* impl_msg);
  __idx_Main_CkArgMsg = CkRegisterEp("Main(CkArgMsg* impl_msg)",
     (CkCallFnPtr)_call_Main_CkArgMsg, CMessage_CkArgMsg::__idx, __idx, 0);
  CkRegisterMainChare(__idx, __idx_Main_CkArgMsg);

// REG: threaded void create_workers(const CkCallback &impl_noname_0);
  __idx_create_workers_marshall2 = CkRegisterEp("create_workers(const CkCallback &impl_noname_0)",
     (CkCallFnPtr)_call_create_workers_marshall2, CkMarshallMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(__idx_create_workers_marshall2,(CkMessagePupFn)_marshallmessagepup_create_workers_marshall2);

// REG: threaded void get_rhs(void);
  __idx_get_rhs_void = CkRegisterEp("get_rhs(void)",
     (CkCallFnPtr)_call_get_rhs_void, 0, __idx, 0);

// REG: void quiescenceHandler(void);
  __idx_quiescenceHandler_void = CkRegisterEp("quiescenceHandler(void)",
     (CkCallFnPtr)_call_quiescenceHandler_void, 0, __idx, 0);

      _registerInitCall(Main::initManualLB,1);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registermain(void)
{
  static int _done = 0; if(_done) return; _done = 1;
      _registermesh_working_copy();

      _registerfmm_globals_nodegroup();

      _registermesh_library();

      _registerparallel_fmm_octree();

      _registerfmm_worker();

      _registerexplicit_worker();

      _registerrhs_handler();

      _registergmres();

      _registeriteration_handler();

      _registerfh_values_nodegroup();

      _registeropencl_nodegroup();

      _registercomlib();

      _registervanilla_fmm_worker();

      _registervanilla_fmm_evals();

/* REG: nodegroup ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >: NodeGroup;
*/
  CkIndex_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >::__register("ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >", sizeof(ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >));

/* REG: nodegroup ParallelFMMOctree<CharmNodePatch,beepp::CProxy_FMMWorker >: NodeGroup;
*/
  CkIndex_ParallelFMMOctree<CharmNodePatch,beepp::CProxy_FMMWorker >::__register("ParallelFMMOctree<CharmNodePatch,beepp::CProxy_FMMWorker >", sizeof(ParallelFMMOctree<CharmNodePatch,beepp::CProxy_FMMWorker >));

  CkRegisterReadonly("MainProxy","CProxy_Main",sizeof(MainProxy),(void *) &MainProxy,__xlater_roPup_MainProxy);

  CkRegisterReadonly("MeshWorkingCopyProxy","CProxy_MeshWorkingCopy",sizeof(MeshWorkingCopyProxy),(void *) &MeshWorkingCopyProxy,__xlater_roPup_MeshWorkingCopyProxy);

  CkRegisterReadonly("MeshLibraryProxy","CProxy_MeshLibrary",sizeof(MeshLibraryProxy),(void *) &MeshLibraryProxy,__xlater_roPup_MeshLibraryProxy);

  CkRegisterReadonly("FMM_Globals_Proxy","beepp::CProxy_FMM_Globals_NodeGroup",sizeof(FMM_Globals_Proxy),(void *) &FMM_Globals_Proxy,__xlater_roPup_FMM_Globals_Proxy);

  CkRegisterReadonly("ParallelFMMOctreeProxy","beepp::CProxy_ParallelTree",sizeof(ParallelFMMOctreeProxy),(void *) &ParallelFMMOctreeProxy,__xlater_roPup_ParallelFMMOctreeProxy);

  CkRegisterReadonly("FMMWorkerProxy","CProxy_FMMWorker",sizeof(FMMWorkerProxy),(void *) &FMMWorkerProxy,__xlater_roPup_FMMWorkerProxy);

  CkRegisterReadonly("ExplicitWorkerProxy","CProxy_ExplicitWorker",sizeof(ExplicitWorkerProxy),(void *) &ExplicitWorkerProxy,__xlater_roPup_ExplicitWorkerProxy);

  CkRegisterReadonly("RHS_HandlerProxy","CProxy_RHS_Handler",sizeof(RHS_HandlerProxy),(void *) &RHS_HandlerProxy,__xlater_roPup_RHS_HandlerProxy);

  CkRegisterReadonly("GMRESProxy","CProxy_GMRES",sizeof(GMRESProxy),(void *) &GMRESProxy,__xlater_roPup_GMRESProxy);

  CkRegisterReadonly("IterationHandlerProxy","CProxy_IterationHandler",sizeof(IterationHandlerProxy),(void *) &IterationHandlerProxy,__xlater_roPup_IterationHandlerProxy);

  CkRegisterReadonly("FH_Values_NodeGroupProxy","CProxy_FH_Values_NodeGroup",sizeof(FH_Values_NodeGroupProxy),(void *) &FH_Values_NodeGroupProxy,__xlater_roPup_FH_Values_NodeGroupProxy);

  CkRegisterReadonly("OpenCL_NodeGroupProxy","CProxy_OpenCL_NodeGroup",sizeof(OpenCL_NodeGroupProxy),(void *) &OpenCL_NodeGroupProxy,__xlater_roPup_OpenCL_NodeGroupProxy);

  CkRegisterReadonly("streaming_strat","ComlibInstanceHandle",sizeof(streaming_strat),(void *) &streaming_strat,__xlater_roPup_streaming_strat);

/* REG: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
threaded void create_workers(const CkCallback &impl_noname_0);
threaded void get_rhs(void);
void quiescenceHandler(void);
  initcall void initManualLB(void);
};
*/
  CkIndex_Main::__register("Main", sizeof(Main));

}
extern "C" void CkRegisterMainModule(void) {
  _registermain();
}
#endif
