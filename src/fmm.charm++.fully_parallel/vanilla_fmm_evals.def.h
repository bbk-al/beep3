namespace vanilla_fmm {
/* DEFS: message Eval_Message{
EvalPt data[];
}
;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_Eval_Message::operator new(size_t s){
  return Eval_Message::alloc(__idx, s, 0, 0);
}
void *CMessage_Eval_Message::operator new(size_t s, int* sz){
  return Eval_Message::alloc(__idx, s, sz, 0);
}
void *CMessage_Eval_Message::operator new(size_t s, int* sz,const int pb){
  return Eval_Message::alloc(__idx, s, sz, pb);
}
void *CMessage_Eval_Message::operator new(size_t s, int sz0) {
  int sizes[1];
  sizes[0] = sz0;
  return Eval_Message::alloc(__idx, s, sizes, 0);
}
void *CMessage_Eval_Message::operator new(size_t s, int sz0, const int p) {
  int sizes[1];
  sizes[0] = sz0;
  return Eval_Message::alloc(__idx, s, sizes, p);
}
void* CMessage_Eval_Message::alloc(int msgnum, size_t sz, int *sizes, int pb) {
  CkpvAccess(_offsets)[0] = ALIGN8(sz);
  if(sizes==0)
    CkpvAccess(_offsets)[1] = CkpvAccess(_offsets)[0];
  else
    CkpvAccess(_offsets)[1] = CkpvAccess(_offsets)[0] + ALIGN8(sizeof(EvalPt)*sizes[0]);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[1], pb);
}
CMessage_Eval_Message::CMessage_Eval_Message() {
Eval_Message *newmsg = (Eval_Message *)this;
  newmsg->data = (EvalPt *) ((char *)newmsg + CkpvAccess(_offsets)[0]);
}
void CMessage_Eval_Message::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_Eval_Message::pack(Eval_Message *msg) {
  msg->data = (EvalPt *) ((char *)msg->data - (char *)msg);
  return (void *) msg;
}
Eval_Message* CMessage_Eval_Message::unpack(void* buf) {
  Eval_Message *msg = (Eval_Message *) buf;
  msg->data = (EvalPt *) ((size_t)msg->data + (char *)msg);
  return msg;
}
int CMessage_Eval_Message::__idx=0;
#endif

/* DEFS: chare Vanilla_FMM_Evals: Chare{
Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy);
void receive_eval_results(Eval_Message* impl_msg);
void evaluate(const CkCallback &cb);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Vanilla_FMM_Evals::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy);
 */
CkChareID CProxy_Vanilla_FMM_Evals::ckNew(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<Vector > &)pts;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<Vector > &)pts;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
  }
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Vanilla_FMM_Evals::__idx, CkIndex_Vanilla_FMM_Evals::__idx_Vanilla_FMM_Evals_marshall1, impl_msg, &impl_ret, impl_onPE);
  return impl_ret;
}
void CProxy_Vanilla_FMM_Evals::ckNew(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, CkChareID* pcid, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<Vector > &)pts;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<Vector > &)pts;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
  }
  CkCreateChare(CkIndex_Vanilla_FMM_Evals::__idx, CkIndex_Vanilla_FMM_Evals::__idx_Vanilla_FMM_Evals_marshall1, impl_msg, pcid, impl_onPE);
}
  CProxy_Vanilla_FMM_Evals::CProxy_Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, int impl_onPE, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<Vector > &)pts;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<Vector > &)pts;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &)ParallelFMMOctreeProxy;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_VanillaFMMWorker &)VanillaFMMWorkerProxy;
  }
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Vanilla_FMM_Evals::__idx, CkIndex_Vanilla_FMM_Evals::__idx_Vanilla_FMM_Evals_marshall1, impl_msg, &impl_ret, impl_onPE);
  ckSetChareID(impl_ret);
}
 int CkIndex_Vanilla_FMM_Evals::__idx_Vanilla_FMM_Evals_marshall1=0;
void CkIndex_Vanilla_FMM_Evals::_call_Vanilla_FMM_Evals_marshall1(void* impl_msg,Vanilla_FMM_Evals * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy*/
  PUP::fromMem implP(impl_buf);
  std::vector<Vector > pts; implP|pts;
  CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy; implP|ParallelFMMOctreeProxy;
  CProxy_VanillaFMMWorker VanillaFMMWorkerProxy; implP|VanillaFMMWorkerProxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) Vanilla_FMM_Evals(pts, ParallelFMMOctreeProxy, VanillaFMMWorkerProxy);
}
int CkIndex_Vanilla_FMM_Evals::_callmarshall_Vanilla_FMM_Evals_marshall1(char* impl_buf,Vanilla_FMM_Evals * impl_obj) {
  /*Unmarshall pup'd fields: const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy*/
  PUP::fromMem implP(impl_buf);
  std::vector<Vector > pts; implP|pts;
  CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy; implP|ParallelFMMOctreeProxy;
  CProxy_VanillaFMMWorker VanillaFMMWorkerProxy; implP|VanillaFMMWorkerProxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) Vanilla_FMM_Evals(pts, ParallelFMMOctreeProxy, VanillaFMMWorkerProxy);
  return implP.size();
}
void CkIndex_Vanilla_FMM_Evals::_marshallmessagepup_Vanilla_FMM_Evals_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy*/
  PUP::fromMem implP(impl_buf);
  std::vector<Vector > pts; implP|pts;
  CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy; implP|ParallelFMMOctreeProxy;
  CProxy_VanillaFMMWorker VanillaFMMWorkerProxy; implP|VanillaFMMWorkerProxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("pts");
  implDestP|pts;
  if (implDestP.hasComments()) implDestP.comment("ParallelFMMOctreeProxy");
  implDestP|ParallelFMMOctreeProxy;
  if (implDestP.hasComments()) implDestP.comment("VanillaFMMWorkerProxy");
  implDestP|VanillaFMMWorkerProxy;
}

/* DEFS: void receive_eval_results(Eval_Message* impl_msg);
 */
void CProxy_Vanilla_FMM_Evals::receive_eval_results(Eval_Message* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Vanilla_FMM_Evals::__idx_receive_eval_results_Eval_Message, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Vanilla_FMM_Evals::__idx_receive_eval_results_Eval_Message, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Vanilla_FMM_Evals::__idx_receive_eval_results_Eval_Message, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_Vanilla_FMM_Evals::__idx_receive_eval_results_Eval_Message=0;
void CkIndex_Vanilla_FMM_Evals::_call_receive_eval_results_Eval_Message(void* impl_msg,Vanilla_FMM_Evals * impl_obj)
{
  impl_obj->receive_eval_results((Eval_Message*)impl_msg);
}

/* DEFS: void evaluate(const CkCallback &cb);
 */
void CProxy_Vanilla_FMM_Evals::evaluate(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Vanilla_FMM_Evals::__idx_evaluate_marshall3, impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Vanilla_FMM_Evals::__idx_evaluate_marshall3, impl_msg, &ckGetChareID(),destPE);
  }
  else CkSendMsg(CkIndex_Vanilla_FMM_Evals::__idx_evaluate_marshall3, impl_msg, &ckGetChareID(),0);
}
 int CkIndex_Vanilla_FMM_Evals::__idx_evaluate_marshall3=0;
void CkIndex_Vanilla_FMM_Evals::_call_evaluate_marshall3(void* impl_msg,Vanilla_FMM_Evals * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->evaluate(cb);
}
int CkIndex_Vanilla_FMM_Evals::_callmarshall_evaluate_marshall3(char* impl_buf,Vanilla_FMM_Evals * impl_obj) {
  /*Unmarshall pup'd fields: const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->evaluate(cb);
  return implP.size();
}
void CkIndex_Vanilla_FMM_Evals::_marshallmessagepup_evaluate_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("cb");
  implDestP|cb;
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Vanilla_FMM_Evals::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeChare);
  CkRegisterBase(__idx, CkIndex_Chare::__idx);
// REG: Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy);
  __idx_Vanilla_FMM_Evals_marshall1 = CkRegisterEp("Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy)",
     (CkCallFnPtr)_call_Vanilla_FMM_Evals_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_Vanilla_FMM_Evals_marshall1,(CkMarshallUnpackFn)_callmarshall_Vanilla_FMM_Evals_marshall1);
  CkRegisterMessagePupFn(__idx_Vanilla_FMM_Evals_marshall1,(CkMessagePupFn)_marshallmessagepup_Vanilla_FMM_Evals_marshall1);

// REG: void receive_eval_results(Eval_Message* impl_msg);
  __idx_receive_eval_results_Eval_Message = CkRegisterEp("receive_eval_results(Eval_Message* impl_msg)",
     (CkCallFnPtr)_call_receive_eval_results_Eval_Message, CMessage_Eval_Message::__idx, __idx, 0);

// REG: void evaluate(const CkCallback &cb);
  __idx_evaluate_marshall3 = CkRegisterEp("evaluate(const CkCallback &cb)",
     (CkCallFnPtr)_call_evaluate_marshall3, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_evaluate_marshall3,(CkMarshallUnpackFn)_callmarshall_evaluate_marshall3);
  CkRegisterMessagePupFn(__idx_evaluate_marshall3,(CkMessagePupFn)_marshallmessagepup_evaluate_marshall3);

}
#endif

} // namespace vanilla_fmm

#ifndef CK_TEMPLATES_ONLY
void _registervanilla_fmm_evals(void)
{
  static int _done = 0; if(_done) return; _done = 1;
using namespace vanilla_fmm;
/* REG: message Eval_Message{
EvalPt data[];
}
;
*/
CMessage_Eval_Message::__register("Eval_Message", sizeof(Eval_Message),(CkPackFnPtr) Eval_Message::pack,(CkUnpackFnPtr) Eval_Message::unpack);

/* REG: chare Vanilla_FMM_Evals: Chare{
Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy);
void receive_eval_results(Eval_Message* impl_msg);
void evaluate(const CkCallback &cb);
};
*/
  CkIndex_Vanilla_FMM_Evals::__register("Vanilla_FMM_Evals", sizeof(Vanilla_FMM_Evals));


}
#endif
