/* DEFS: array ExplicitWorker: ArrayElement{
ExplicitWorker(CkMigrateMessage* impl_msg);
ExplicitWorker(void);
void precalc(double kappa, unsigned int total_num_patches);
void evaluate(double kappa);
void triggerLB(void);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_ExplicitWorker::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: ExplicitWorker(CkMigrateMessage* impl_msg);
 */

/* DEFS: ExplicitWorker(void);
 */
void CProxyElement_ExplicitWorker::insert(int onPE)
{ 
  void *impl_msg = CkAllocSysMsg();
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_ExplicitWorker::__idx_ExplicitWorker_void,onPE);
}

/* DEFS: void precalc(double kappa, unsigned int total_num_patches);
 */
void CProxyElement_ExplicitWorker::precalc(double kappa, unsigned int total_num_patches, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa, unsigned int total_num_patches
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    implP|total_num_patches;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
    implP|total_num_patches;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_ExplicitWorker::__idx_precalc_marshall2,0);
}

/* DEFS: void evaluate(double kappa);
 */
void CProxyElement_ExplicitWorker::evaluate(double kappa, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_ExplicitWorker::__idx_evaluate_marshall3,0);
}

/* DEFS: void triggerLB(void);
 */
void CProxyElement_ExplicitWorker::triggerLB(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_ExplicitWorker::__idx_triggerLB_void,0);
}

/* DEFS: ExplicitWorker(CkMigrateMessage* impl_msg);
 */
 int CkIndex_ExplicitWorker::__idx_ExplicitWorker_CkMigrateMessage=0;
void CkIndex_ExplicitWorker::_call_ExplicitWorker_CkMigrateMessage(void* impl_msg,ExplicitWorker * impl_obj)
{
  new (impl_obj) ExplicitWorker((CkMigrateMessage*)impl_msg);
}

/* DEFS: ExplicitWorker(void);
 */
CkArrayID CProxy_ExplicitWorker::ckNew(const CkArrayOptions &opts)
{ 
  void *impl_msg = CkAllocSysMsg();
   return ckCreateArray((CkArrayMessage *)impl_msg,CkIndex_ExplicitWorker::__idx_ExplicitWorker_void,opts);
}
 int CkIndex_ExplicitWorker::__idx_ExplicitWorker_void=0;
void CkIndex_ExplicitWorker::_call_ExplicitWorker_void(void* impl_msg,ExplicitWorker * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  new (impl_obj) ExplicitWorker();
}

/* DEFS: void precalc(double kappa, unsigned int total_num_patches);
 */
void CProxy_ExplicitWorker::precalc(double kappa, unsigned int total_num_patches, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa, unsigned int total_num_patches
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    implP|total_num_patches;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
    implP|total_num_patches;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_ExplicitWorker::__idx_precalc_marshall2,0);
}
 int CkIndex_ExplicitWorker::__idx_precalc_marshall2=0;
void CkIndex_ExplicitWorker::_call_precalc_marshall2(void* impl_msg,ExplicitWorker * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double kappa, unsigned int total_num_patches*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  unsigned int total_num_patches; implP|total_num_patches;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->precalc(kappa, total_num_patches);
}
int CkIndex_ExplicitWorker::_callmarshall_precalc_marshall2(char* impl_buf,ExplicitWorker * impl_obj) {
  /*Unmarshall pup'd fields: double kappa, unsigned int total_num_patches*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  unsigned int total_num_patches; implP|total_num_patches;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->precalc(kappa, total_num_patches);
  return implP.size();
}
void CkIndex_ExplicitWorker::_marshallmessagepup_precalc_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double kappa, unsigned int total_num_patches*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  unsigned int total_num_patches; implP|total_num_patches;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("kappa");
  implDestP|kappa;
  if (implDestP.hasComments()) implDestP.comment("total_num_patches");
  implDestP|total_num_patches;
}

/* DEFS: void evaluate(double kappa);
 */
void CProxy_ExplicitWorker::evaluate(double kappa, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_ExplicitWorker::__idx_evaluate_marshall3,0);
}
 int CkIndex_ExplicitWorker::__idx_evaluate_marshall3=0;
void CkIndex_ExplicitWorker::_call_evaluate_marshall3(void* impl_msg,ExplicitWorker * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double kappa*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->evaluate(kappa);
}
int CkIndex_ExplicitWorker::_callmarshall_evaluate_marshall3(char* impl_buf,ExplicitWorker * impl_obj) {
  /*Unmarshall pup'd fields: double kappa*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->evaluate(kappa);
  return implP.size();
}
void CkIndex_ExplicitWorker::_marshallmessagepup_evaluate_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: double kappa*/
  PUP::fromMem implP(impl_buf);
  double kappa; implP|kappa;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("kappa");
  implDestP|kappa;
}

/* DEFS: void triggerLB(void);
 */
void CProxy_ExplicitWorker::triggerLB(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_ExplicitWorker::__idx_triggerLB_void,0);
}
 int CkIndex_ExplicitWorker::__idx_triggerLB_void=0;
void CkIndex_ExplicitWorker::_call_triggerLB_void(void* impl_msg,ExplicitWorker * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->triggerLB();
}

/* DEFS: ExplicitWorker(CkMigrateMessage* impl_msg);
 */

/* DEFS: ExplicitWorker(void);
 */

/* DEFS: void precalc(double kappa, unsigned int total_num_patches);
 */
void CProxySection_ExplicitWorker::precalc(double kappa, unsigned int total_num_patches, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa, unsigned int total_num_patches
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    implP|total_num_patches;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
    implP|total_num_patches;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_ExplicitWorker::__idx_precalc_marshall2,0);
}

/* DEFS: void evaluate(double kappa);
 */
void CProxySection_ExplicitWorker::evaluate(double kappa, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|kappa;
  }
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_ExplicitWorker::__idx_evaluate_marshall3,0);
}

/* DEFS: void triggerLB(void);
 */
void CProxySection_ExplicitWorker::triggerLB(void) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_ExplicitWorker::__idx_triggerLB_void,0);
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_ExplicitWorker::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
// REG: ExplicitWorker(CkMigrateMessage* impl_msg);
  __idx_ExplicitWorker_CkMigrateMessage = CkRegisterEp("ExplicitWorker(CkMigrateMessage* impl_msg)",
     (CkCallFnPtr)_call_ExplicitWorker_CkMigrateMessage, 0, __idx, 0);
  CkRegisterMigCtor(__idx, __idx_ExplicitWorker_CkMigrateMessage);

// REG: ExplicitWorker(void);
  __idx_ExplicitWorker_void = CkRegisterEp("ExplicitWorker(void)",
     (CkCallFnPtr)_call_ExplicitWorker_void, 0, __idx, 0);
  CkRegisterDefaultCtor(__idx, __idx_ExplicitWorker_void);

// REG: void precalc(double kappa, unsigned int total_num_patches);
  __idx_precalc_marshall2 = CkRegisterEp("precalc(double kappa, unsigned int total_num_patches)",
     (CkCallFnPtr)_call_precalc_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_precalc_marshall2,(CkMarshallUnpackFn)_callmarshall_precalc_marshall2);
  CkRegisterMessagePupFn(__idx_precalc_marshall2,(CkMessagePupFn)_marshallmessagepup_precalc_marshall2);

// REG: void evaluate(double kappa);
  __idx_evaluate_marshall3 = CkRegisterEp("evaluate(double kappa)",
     (CkCallFnPtr)_call_evaluate_marshall3, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_evaluate_marshall3,(CkMarshallUnpackFn)_callmarshall_evaluate_marshall3);
  CkRegisterMessagePupFn(__idx_evaluate_marshall3,(CkMessagePupFn)_marshallmessagepup_evaluate_marshall3);

// REG: void triggerLB(void);
  __idx_triggerLB_void = CkRegisterEp("triggerLB(void)",
     (CkCallFnPtr)_call_triggerLB_void, 0, __idx, 0);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registerexplicit_worker(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: array ExplicitWorker: ArrayElement{
ExplicitWorker(CkMigrateMessage* impl_msg);
ExplicitWorker(void);
void precalc(double kappa, unsigned int total_num_patches);
void evaluate(double kappa);
void triggerLB(void);
};
*/
  CkIndex_ExplicitWorker::__register("ExplicitWorker", sizeof(ExplicitWorker));

}
#endif
