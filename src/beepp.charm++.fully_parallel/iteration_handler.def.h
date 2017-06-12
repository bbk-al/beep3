
/* DEFS: chare IterationHandler: Chare{
IterationHandler(double _universe_edge_length, double _beta);
void do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals);
threaded void phase_two(void);
void done_evals(void);
void reduce_fh_results(void);
void done_reduction(CkReductionMsg* impl_msg);
sync void set_num_workers(unsigned int num_workers_in);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_IterationHandler::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: IterationHandler(double _universe_edge_length, double _beta);
 */
CkChareID CProxy_IterationHandler::ckNew(double _universe_edge_length, double _beta, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: double _universe_edge_length, double _beta
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|_universe_edge_length;
    implP|_beta;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|_universe_edge_length;
    implP|_beta;
  }
  CkChareID impl_ret;
  CkCreateChare(CkIndex_IterationHandler::__idx, CkIndex_IterationHandler::__idx_IterationHandler_marshall1, impl_msg, &impl_ret, impl_onPE);
  return impl_ret;
}
void CProxy_IterationHandler::ckNew(double _universe_edge_length, double _beta, CkChareID* pcid, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: double _universe_edge_length, double _beta
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|_universe_edge_length;
    implP|_beta;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|_universe_edge_length;
    implP|_beta;
  }
  CkCreateChare(CkIndex_IterationHandler::__idx, CkIndex_IterationHandler::__idx_IterationHandler_marshall1, impl_msg, pcid, impl_onPE);
}
  CProxy_IterationHandler::CProxy_IterationHandler(double _universe_edge_length, double _beta, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: double _universe_edge_length, double _beta
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|_universe_edge_length;
    implP|_beta;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|_universe_edge_length;
    implP|_beta;
  }
  CkChareID impl_ret;
  CkCreateChare(CkIndex_IterationHandler::__idx, CkIndex_IterationHandler::__idx_IterationHandler_marshall1, impl_msg, &impl_ret, impl_onPE);
  ckSetChareID(impl_ret);
}
 int CkIndex_IterationHandler::__idx_IterationHandler_marshall1=0;
void CkIndex_IterationHandler::_call_IterationHandler_marshall1(void* impl_msg,IterationHandler * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double _universe_edge_length, double _beta*/
  PUP::fromMem implP(impl_buf);
  double _universe_edge_length; implP|_universe_edge_length;
  double _beta; implP|_beta;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) IterationHandler(_universe_edge_length, _beta);
}
int CkIndex_IterationHandler::_callmarshall_IterationHandler_marshall1(char* impl_buf,IterationHandler * impl_obj) {
  /*Unmarshall pup'd fields: double _universe_edge_length, double _beta*/
  PUP::fromMem implP(impl_buf);
  double _universe_edge_length; implP|_universe_edge_length;
  double _beta; implP|_beta;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) IterationHandler(_universe_edge_length, _beta);
  return implP.size();
}
void CkIndex_IterationHandler::_marshallmessagepup_IterationHandler_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double _universe_edge_length, double _beta*/
  PUP::fromMem implP(impl_buf);
  double _universe_edge_length; implP|_universe_edge_length;
  double _beta; implP|_beta;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("_universe_edge_length");
  implDestP|_universe_edge_length;
  if (implDestP.hasComments()) implDestP.comment("_beta");
  implDestP|_beta;
}

/* DEFS: void do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals);
 */
void CProxy_IterationHandler::do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkCallback &cb, const FH_Values &fh_vals
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)fh_vals;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)fh_vals;
  }
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_IterationHandler::__idx_do_bem_fmm_iteration_marshall2, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_IterationHandler::__idx_do_bem_fmm_iteration_marshall2, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_IterationHandler::__idx_do_bem_fmm_iteration_marshall2, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_IterationHandler::__idx_do_bem_fmm_iteration_marshall2=0;
void CkIndex_IterationHandler::_call_do_bem_fmm_iteration_marshall2(void* impl_msg,IterationHandler * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &cb, const FH_Values &fh_vals*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  FH_Values fh_vals; implP|fh_vals;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->do_bem_fmm_iteration(cb, fh_vals);
}
int CkIndex_IterationHandler::_callmarshall_do_bem_fmm_iteration_marshall2(char* impl_buf,IterationHandler * impl_obj) {
  /*Unmarshall pup'd fields: const CkCallback &cb, const FH_Values &fh_vals*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  FH_Values fh_vals; implP|fh_vals;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->do_bem_fmm_iteration(cb, fh_vals);
  return implP.size();
}
void CkIndex_IterationHandler::_marshallmessagepup_do_bem_fmm_iteration_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &cb, const FH_Values &fh_vals*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  FH_Values fh_vals; implP|fh_vals;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("cb");
  implDestP|cb;
  if (implDestP.hasComments()) implDestP.comment("fh_vals");
  implDestP|fh_vals;
}

/* DEFS: threaded void phase_two(void);
 */
void CProxy_IterationHandler::phase_two(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_IterationHandler::__idx_phase_two_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_IterationHandler::__idx_phase_two_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_IterationHandler::__idx_phase_two_void, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_IterationHandler::__idx_phase_two_void=0;
void CkIndex_IterationHandler::_call_phase_two_void(void* impl_msg,IterationHandler * impl_obj)
{
  CthThread tid = CthCreate((CthVoidFn)_callthr_phase_two_void, new CkThrCallArg(impl_msg,impl_obj), 0);
  ((Chare *)impl_obj)->CkAddThreadListeners(tid,impl_msg);
  CthAwaken(tid);
}
void CkIndex_IterationHandler::_callthr_phase_two_void(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  IterationHandler *impl_obj = (IterationHandler *) impl_arg->obj;
  delete impl_arg;
  CkFreeSysMsg(impl_msg);
  impl_obj->phase_two();
}

/* DEFS: void done_evals(void);
 */
void CProxy_IterationHandler::done_evals(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_IterationHandler::__idx_done_evals_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_IterationHandler::__idx_done_evals_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_IterationHandler::__idx_done_evals_void, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_IterationHandler::__idx_done_evals_void=0;
void CkIndex_IterationHandler::_call_done_evals_void(void* impl_msg,IterationHandler * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->done_evals();
}

/* DEFS: void reduce_fh_results(void);
 */
void CProxy_IterationHandler::reduce_fh_results(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_IterationHandler::__idx_reduce_fh_results_void, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_IterationHandler::__idx_reduce_fh_results_void, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_IterationHandler::__idx_reduce_fh_results_void, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_IterationHandler::__idx_reduce_fh_results_void=0;
void CkIndex_IterationHandler::_call_reduce_fh_results_void(void* impl_msg,IterationHandler * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->reduce_fh_results();
}

/* DEFS: void done_reduction(CkReductionMsg* impl_msg);
 */
void CProxy_IterationHandler::done_reduction(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_IterationHandler::__idx_done_reduction_CkReductionMsg, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_IterationHandler::__idx_done_reduction_CkReductionMsg, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_IterationHandler::__idx_done_reduction_CkReductionMsg, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_IterationHandler::__idx_done_reduction_CkReductionMsg=0;
void CkIndex_IterationHandler::_call_done_reduction_CkReductionMsg(void* impl_msg,IterationHandler * impl_obj)
{
  impl_obj->done_reduction((CkReductionMsg*)impl_msg);
}

/* DEFS: sync void set_num_workers(unsigned int num_workers_in);
 */
void CProxy_IterationHandler::set_num_workers(unsigned int num_workers_in, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: unsigned int num_workers_in
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_workers_in;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_workers_in;
  }
  CkFreeSysMsg(CkRemoteCall(CkIndex_IterationHandler::__idx_set_num_workers_marshall7, impl_msg, &ckGetChareID()));
}
 int CkIndex_IterationHandler::__idx_set_num_workers_marshall7=0;
void CkIndex_IterationHandler::_call_set_num_workers_marshall7(void* impl_msg,IterationHandler * impl_obj)
{
  int impl_ref = CkGetRefNum(impl_msg), impl_src = CkGetSrcPe(impl_msg);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int num_workers_in*/
  PUP::fromMem implP(impl_buf);
  unsigned int num_workers_in; implP|num_workers_in;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  void *impl_retMsg=CkAllocSysMsg();
    impl_obj->set_num_workers(num_workers_in);
  CkSendToFutureID(impl_ref, impl_retMsg, impl_src);
}
void CkIndex_IterationHandler::_marshallmessagepup_set_num_workers_marshall7(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int num_workers_in*/
  PUP::fromMem implP(impl_buf);
  unsigned int num_workers_in; implP|num_workers_in;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("num_workers_in");
  implDestP|num_workers_in;
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_IterationHandler::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeChare);
  CkRegisterBase(__idx, CkIndex_Chare::__idx);
// REG: IterationHandler(double _universe_edge_length, double _beta);
  __idx_IterationHandler_marshall1 = CkRegisterEp("IterationHandler(double _universe_edge_length, double _beta)",
     (CkCallFnPtr)_call_IterationHandler_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_IterationHandler_marshall1,(CkMarshallUnpackFn)_callmarshall_IterationHandler_marshall1);
  CkRegisterMessagePupFn(__idx_IterationHandler_marshall1,(CkMessagePupFn)_marshallmessagepup_IterationHandler_marshall1);

// REG: void do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals);
  __idx_do_bem_fmm_iteration_marshall2 = CkRegisterEp("do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals)",
     (CkCallFnPtr)_call_do_bem_fmm_iteration_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_do_bem_fmm_iteration_marshall2,(CkMarshallUnpackFn)_callmarshall_do_bem_fmm_iteration_marshall2);
  CkRegisterMessagePupFn(__idx_do_bem_fmm_iteration_marshall2,(CkMessagePupFn)_marshallmessagepup_do_bem_fmm_iteration_marshall2);

// REG: threaded void phase_two(void);
  __idx_phase_two_void = CkRegisterEp("phase_two(void)",
     (CkCallFnPtr)_call_phase_two_void, 0, __idx, 0);

// REG: void done_evals(void);
  __idx_done_evals_void = CkRegisterEp("done_evals(void)",
     (CkCallFnPtr)_call_done_evals_void, 0, __idx, 0);

// REG: void reduce_fh_results(void);
  __idx_reduce_fh_results_void = CkRegisterEp("reduce_fh_results(void)",
     (CkCallFnPtr)_call_reduce_fh_results_void, 0, __idx, 0);

// REG: void done_reduction(CkReductionMsg* impl_msg);
  __idx_done_reduction_CkReductionMsg = CkRegisterEp("done_reduction(CkReductionMsg* impl_msg)",
     (CkCallFnPtr)_call_done_reduction_CkReductionMsg, CMessage_CkReductionMsg::__idx, __idx, 0);

// REG: sync void set_num_workers(unsigned int num_workers_in);
  __idx_set_num_workers_marshall7 = CkRegisterEp("set_num_workers(unsigned int num_workers_in)",
     (CkCallFnPtr)_call_set_num_workers_marshall7, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMessagePupFn(__idx_set_num_workers_marshall7,(CkMessagePupFn)_marshallmessagepup_set_num_workers_marshall7);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registeriteration_handler(void)
{
  static int _done = 0; if(_done) return; _done = 1;
      _registerfmm_worker();

/* REG: chare IterationHandler: Chare{
IterationHandler(double _universe_edge_length, double _beta);
void do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals);
threaded void phase_two(void);
void done_evals(void);
void reduce_fh_results(void);
void done_reduction(CkReductionMsg* impl_msg);
sync void set_num_workers(unsigned int num_workers_in);
};
*/
  CkIndex_IterationHandler::__register("IterationHandler", sizeof(IterationHandler));

}
#endif
