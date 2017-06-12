/* DEFS: group FH_Values_NodeGroup: IrrGroup{
FH_Values_NodeGroup(unsigned int impl_noname_0);
void set(const FH_Values &vals, const CkCallback &cb);
void reduce(const CkCallback &cb);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_FH_Values_NodeGroup::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: FH_Values_NodeGroup(unsigned int impl_noname_0);
 */

/* DEFS: void set(const FH_Values &vals, const CkCallback &cb);
 */
void CProxyElement_FH_Values_NodeGroup::set(const FH_Values &vals, const CkCallback &cb, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const FH_Values &vals, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgBranch(CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: void reduce(const CkCallback &cb);
 */
void CProxyElement_FH_Values_NodeGroup::reduce(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkGroupMsgPrep(CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgBranch(CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: FH_Values_NodeGroup(unsigned int impl_noname_0);
 */
CkGroupID CProxy_FH_Values_NodeGroup::ckNew(unsigned int impl_noname_0, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int impl_noname_0
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
  return CkCreateGroup(CkIndex_FH_Values_NodeGroup::__idx, CkIndex_FH_Values_NodeGroup::__idx_FH_Values_NodeGroup_marshall1, impl_msg);
}
  CProxy_FH_Values_NodeGroup::CProxy_FH_Values_NodeGroup(unsigned int impl_noname_0, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int impl_noname_0
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
  ckSetGroupID(CkCreateGroup(CkIndex_FH_Values_NodeGroup::__idx, CkIndex_FH_Values_NodeGroup::__idx_FH_Values_NodeGroup_marshall1, impl_msg));
}
 int CkIndex_FH_Values_NodeGroup::__idx_FH_Values_NodeGroup_marshall1=0;
void CkIndex_FH_Values_NodeGroup::_call_FH_Values_NodeGroup_marshall1(void* impl_msg,FH_Values_NodeGroup * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  unsigned int impl_noname_0; implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) FH_Values_NodeGroup(impl_noname_0);
}
int CkIndex_FH_Values_NodeGroup::_callmarshall_FH_Values_NodeGroup_marshall1(char* impl_buf,FH_Values_NodeGroup * impl_obj) {
  /*Unmarshall pup'd fields: unsigned int impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  unsigned int impl_noname_0; implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) FH_Values_NodeGroup(impl_noname_0);
  return implP.size();
}
void CkIndex_FH_Values_NodeGroup::_marshallmessagepup_FH_Values_NodeGroup_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  unsigned int impl_noname_0; implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
}

/* DEFS: void set(const FH_Values &vals, const CkCallback &cb);
 */
void CProxy_FH_Values_NodeGroup::set(const FH_Values &vals, const CkCallback &cb, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const FH_Values &vals, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupID(),0);
}
void CProxy_FH_Values_NodeGroup::set(const FH_Values &vals, const CkCallback &cb, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const FH_Values &vals, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkSendMsgBranchMulti(CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_FH_Values_NodeGroup::set(const FH_Values &vals, const CkCallback &cb, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const FH_Values &vals, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  CkSendMsgBranchGroup(CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupID(), grp,0);
}
 int CkIndex_FH_Values_NodeGroup::__idx_set_marshall2=0;
void CkIndex_FH_Values_NodeGroup::_call_set_marshall2(void* impl_msg,FH_Values_NodeGroup * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const FH_Values &vals, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  FH_Values vals; implP|vals;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->set(vals, cb);
}
int CkIndex_FH_Values_NodeGroup::_callmarshall_set_marshall2(char* impl_buf,FH_Values_NodeGroup * impl_obj) {
  /*Unmarshall pup'd fields: const FH_Values &vals, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  FH_Values vals; implP|vals;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->set(vals, cb);
  return implP.size();
}
void CkIndex_FH_Values_NodeGroup::_marshallmessagepup_set_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const FH_Values &vals, const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  FH_Values vals; implP|vals;
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("vals");
  implDestP|vals;
  if (implDestP.hasComments()) implDestP.comment("cb");
  implDestP|cb;
}

/* DEFS: void reduce(const CkCallback &cb);
 */
void CProxy_FH_Values_NodeGroup::reduce(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkGroupMsgPrep(CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupID(),0);
}
void CProxy_FH_Values_NodeGroup::reduce(const CkCallback &cb, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
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
  CkSendMsgBranchMulti(CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_FH_Values_NodeGroup::reduce(const CkCallback &cb, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
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
  CkSendMsgBranchGroup(CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupID(), grp,0);
}
 int CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3=0;
void CkIndex_FH_Values_NodeGroup::_call_reduce_marshall3(void* impl_msg,FH_Values_NodeGroup * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->reduce(cb);
}
int CkIndex_FH_Values_NodeGroup::_callmarshall_reduce_marshall3(char* impl_buf,FH_Values_NodeGroup * impl_obj) {
  /*Unmarshall pup'd fields: const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->reduce(cb);
  return implP.size();
}
void CkIndex_FH_Values_NodeGroup::_marshallmessagepup_reduce_marshall3(PUP::er &implDestP,void *impl_msg) {
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

/* DEFS: FH_Values_NodeGroup(unsigned int impl_noname_0);
 */

/* DEFS: void set(const FH_Values &vals, const CkCallback &cb);
 */
void CProxySection_FH_Values_NodeGroup::set(const FH_Values &vals, const CkCallback &cb, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const FH_Values &vals, const CkCallback &cb
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)vals;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)cb;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_FH_Values_NodeGroup::__idx_set_marshall2, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

/* DEFS: void reduce(const CkCallback &cb);
 */
void CProxySection_FH_Values_NodeGroup::reduce(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkGroupMsgPrep(CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_FH_Values_NodeGroup::__idx_reduce_marshall3, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_FH_Values_NodeGroup::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_IrrGroup::__idx);
   CkRegisterGroupIrr(__idx,FH_Values_NodeGroup::isIrreducible());
// REG: FH_Values_NodeGroup(unsigned int impl_noname_0);
  __idx_FH_Values_NodeGroup_marshall1 = CkRegisterEp("FH_Values_NodeGroup(unsigned int impl_noname_0)",
     (CkCallFnPtr)_call_FH_Values_NodeGroup_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_FH_Values_NodeGroup_marshall1,(CkMarshallUnpackFn)_callmarshall_FH_Values_NodeGroup_marshall1);
  CkRegisterMessagePupFn(__idx_FH_Values_NodeGroup_marshall1,(CkMessagePupFn)_marshallmessagepup_FH_Values_NodeGroup_marshall1);

// REG: void set(const FH_Values &vals, const CkCallback &cb);
  __idx_set_marshall2 = CkRegisterEp("set(const FH_Values &vals, const CkCallback &cb)",
     (CkCallFnPtr)_call_set_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_set_marshall2,(CkMarshallUnpackFn)_callmarshall_set_marshall2);
  CkRegisterMessagePupFn(__idx_set_marshall2,(CkMessagePupFn)_marshallmessagepup_set_marshall2);

// REG: void reduce(const CkCallback &cb);
  __idx_reduce_marshall3 = CkRegisterEp("reduce(const CkCallback &cb)",
     (CkCallFnPtr)_call_reduce_marshall3, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_reduce_marshall3,(CkMarshallUnpackFn)_callmarshall_reduce_marshall3);
  CkRegisterMessagePupFn(__idx_reduce_marshall3,(CkMessagePupFn)_marshallmessagepup_reduce_marshall3);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registerfh_values_nodegroup(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: group FH_Values_NodeGroup: IrrGroup{
FH_Values_NodeGroup(unsigned int impl_noname_0);
void set(const FH_Values &vals, const CkCallback &cb);
void reduce(const CkCallback &cb);
};
*/
  CkIndex_FH_Values_NodeGroup::__register("FH_Values_NodeGroup", sizeof(FH_Values_NodeGroup));

}
#endif
