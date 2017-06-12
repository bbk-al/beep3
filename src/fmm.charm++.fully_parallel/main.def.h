





/* DEFS: nodegroup ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >: NodeGroup;
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

/* DEFS: readonly CProxy_FMM_Globals_NodeGroup FMM_Globals_Proxy;
 */
extern CProxy_FMM_Globals_NodeGroup FMM_Globals_Proxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_FMM_Globals_Proxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|FMM_Globals_Proxy;
}
#endif

/* DEFS: readonly CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy;
 */
extern CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_ParallelFMMOctreeProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|ParallelFMMOctreeProxy;
}
#endif

/* DEFS: readonly vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy;
 */
extern vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_VanillaFMMWorkerProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|VanillaFMMWorkerProxy;
}
#endif

/* DEFS: readonly vanilla_fmm::CProxy_Vanilla_FMM_Evals Vanilla_FMM_Evals_Proxy;
 */
extern vanilla_fmm::CProxy_Vanilla_FMM_Evals Vanilla_FMM_Evals_Proxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_Vanilla_FMM_Evals_Proxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|Vanilla_FMM_Evals_Proxy;
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

/* DEFS: readonly CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
 */
extern CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_OpenCL_NodeGroupProxy(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|OpenCL_NodeGroupProxy;
}
#endif

/* DEFS: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
threaded void create_workers(void);
void fmm_worker_complete(void);
void completed(vanilla_fmm::Eval_Message* impl_msg);
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

/* DEFS: threaded void create_workers(void);
 */
void CProxy_Main::create_workers(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_create_workers_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::__idx_create_workers_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_create_workers_void, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_Main::__idx_create_workers_void=0;
void CkIndex_Main::_call_create_workers_void(void* impl_msg,Main * impl_obj)
{
  CthThread tid = CthCreate((CthVoidFn)_callthr_create_workers_void, new CkThrCallArg(impl_msg,impl_obj), 0);
  ((Chare *)impl_obj)->CkAddThreadListeners(tid,impl_msg);
  CthAwaken(tid);
}
void CkIndex_Main::_callthr_create_workers_void(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  Main *impl_obj = (Main *) impl_arg->obj;
  delete impl_arg;
  CkFreeSysMsg(impl_msg);
  impl_obj->create_workers();
}

/* DEFS: void fmm_worker_complete(void);
 */
void CProxy_Main::fmm_worker_complete(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_fmm_worker_complete_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::__idx_fmm_worker_complete_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_fmm_worker_complete_void, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_Main::__idx_fmm_worker_complete_void=0;
void CkIndex_Main::_call_fmm_worker_complete_void(void* impl_msg,Main * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->fmm_worker_complete();
}

/* DEFS: void completed(vanilla_fmm::Eval_Message* impl_msg);
 */
void CProxy_Main::completed(vanilla_fmm::Eval_Message* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::__idx_completed_Eval_Message, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::__idx_completed_Eval_Message, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Main::__idx_completed_Eval_Message, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_Main::__idx_completed_Eval_Message=0;
void CkIndex_Main::_call_completed_Eval_Message(void* impl_msg,Main * impl_obj)
{
  impl_obj->completed((vanilla_fmm::Eval_Message*)impl_msg);
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

// REG: threaded void create_workers(void);
  __idx_create_workers_void = CkRegisterEp("create_workers(void)",
     (CkCallFnPtr)_call_create_workers_void, 0, __idx, 0);

// REG: void fmm_worker_complete(void);
  __idx_fmm_worker_complete_void = CkRegisterEp("fmm_worker_complete(void)",
     (CkCallFnPtr)_call_fmm_worker_complete_void, 0, __idx, 0);

// REG: void completed(vanilla_fmm::Eval_Message* impl_msg);
  __idx_completed_Eval_Message = CkRegisterEp("completed(vanilla_fmm::Eval_Message* impl_msg)",
     (CkCallFnPtr)_call_completed_Eval_Message, vanilla_fmm::CMessage_Eval_Message::__idx, __idx, 0);

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
      _registerfmm_globals_nodegroup();

      _registerparallel_fmm_octree();

      _registeropencl_nodegroup();

      _registervanilla_fmm_worker();

      _registervanilla_fmm_evals();

      _registercomlib();

/* REG: nodegroup ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >: NodeGroup;
*/
  CkIndex_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >::__register("ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >", sizeof(ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >));

  CkRegisterReadonly("MainProxy","CProxy_Main",sizeof(MainProxy),(void *) &MainProxy,__xlater_roPup_MainProxy);

  CkRegisterReadonly("FMM_Globals_Proxy","CProxy_FMM_Globals_NodeGroup",sizeof(FMM_Globals_Proxy),(void *) &FMM_Globals_Proxy,__xlater_roPup_FMM_Globals_Proxy);

  CkRegisterReadonly("ParallelFMMOctreeProxy","CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker >",sizeof(ParallelFMMOctreeProxy),(void *) &ParallelFMMOctreeProxy,__xlater_roPup_ParallelFMMOctreeProxy);

  CkRegisterReadonly("VanillaFMMWorkerProxy","vanilla_fmm::CProxy_VanillaFMMWorker",sizeof(VanillaFMMWorkerProxy),(void *) &VanillaFMMWorkerProxy,__xlater_roPup_VanillaFMMWorkerProxy);

  CkRegisterReadonly("Vanilla_FMM_Evals_Proxy","vanilla_fmm::CProxy_Vanilla_FMM_Evals",sizeof(Vanilla_FMM_Evals_Proxy),(void *) &Vanilla_FMM_Evals_Proxy,__xlater_roPup_Vanilla_FMM_Evals_Proxy);

  CkRegisterReadonly("streaming_strat","ComlibInstanceHandle",sizeof(streaming_strat),(void *) &streaming_strat,__xlater_roPup_streaming_strat);

  CkRegisterReadonly("OpenCL_NodeGroupProxy","CProxy_OpenCL_NodeGroup",sizeof(OpenCL_NodeGroupProxy),(void *) &OpenCL_NodeGroupProxy,__xlater_roPup_OpenCL_NodeGroupProxy);

/* REG: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
threaded void create_workers(void);
void fmm_worker_complete(void);
void completed(vanilla_fmm::Eval_Message* impl_msg);
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
