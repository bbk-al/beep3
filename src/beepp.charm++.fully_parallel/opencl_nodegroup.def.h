/* DEFS: nodegroup OpenCL_NodeGroup: NodeGroup{
OpenCL_NodeGroup(unsigned int num_patches);
void run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
void collate_bem_results(const CkCallback &cb);
void precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_OpenCL_NodeGroup::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: OpenCL_NodeGroup(unsigned int num_patches);
 */

/* DEFS: void run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
 */
void CProxyElement_OpenCL_NodeGroup::run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkArrayIndexOctreeIndexer &idxer, double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSend(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: void collate_bem_results(const CkCallback &cb);
 */
void CProxyElement_OpenCL_NodeGroup::collate_bem_results(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSend(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: void precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
 */
void CProxyElement_OpenCL_NodeGroup::precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkArrayIndexOctreeIndexer &idxer, double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSend(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: OpenCL_NodeGroup(unsigned int num_patches);
 */
CkGroupID CProxy_OpenCL_NodeGroup::ckNew(unsigned int num_patches, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int num_patches
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_patches;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_patches;
  }
  return CkCreateNodeGroup(CkIndex_OpenCL_NodeGroup::__idx, CkIndex_OpenCL_NodeGroup::__idx_OpenCL_NodeGroup_marshall1, impl_msg);
}
  CProxy_OpenCL_NodeGroup::CProxy_OpenCL_NodeGroup(unsigned int num_patches, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int num_patches
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_patches;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_patches;
  }
  ckSetGroupID(CkCreateNodeGroup(CkIndex_OpenCL_NodeGroup::__idx, CkIndex_OpenCL_NodeGroup::__idx_OpenCL_NodeGroup_marshall1, impl_msg));
}
 int CkIndex_OpenCL_NodeGroup::__idx_OpenCL_NodeGroup_marshall1=0;
void CkIndex_OpenCL_NodeGroup::_call_OpenCL_NodeGroup_marshall1(void* impl_msg,OpenCL_NodeGroup * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int num_patches*/
  PUP::fromMem implP(impl_buf);
  unsigned int num_patches; implP|num_patches;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) OpenCL_NodeGroup(num_patches);
}
int CkIndex_OpenCL_NodeGroup::_callmarshall_OpenCL_NodeGroup_marshall1(char* impl_buf,OpenCL_NodeGroup * impl_obj) {
  /*Unmarshall pup'd fields: unsigned int num_patches*/
  PUP::fromMem implP(impl_buf);
  unsigned int num_patches; implP|num_patches;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) OpenCL_NodeGroup(num_patches);
  return implP.size();
}
void CkIndex_OpenCL_NodeGroup::_marshallmessagepup_OpenCL_NodeGroup_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int num_patches*/
  PUP::fromMem implP(impl_buf);
  unsigned int num_patches; implP|num_patches;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("num_patches");
  implDestP|num_patches;
}

/* DEFS: void run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
 */
void CProxy_OpenCL_NodeGroup::run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkArrayIndexOctreeIndexer &idxer, double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupBroadcast(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg, ckGetGroupID(),0);
}
 int CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2=0;
void CkIndex_OpenCL_NodeGroup::_call_run_bem_marshall2(void* impl_msg,OpenCL_NodeGroup * impl_obj)
{
  if(CmiTryLock(impl_obj->__nodelock)) {
    impl_msg = CkCopyMsg(&impl_msg);
    CkSendMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2,impl_msg,CkMyNode(),impl_obj->CkGetNodeGroupID());
    return;
  }
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkArrayIndexOctreeIndexer &idxer, double kappa*/
  PUP::fromMem implP(impl_buf);
  CkArrayIndexOctreeIndexer idxer; implP|idxer;
  double kappa; implP|kappa;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->run_bem(idxer, kappa);
  CmiUnlock(impl_obj->__nodelock);
}
void CkIndex_OpenCL_NodeGroup::_marshallmessagepup_run_bem_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkArrayIndexOctreeIndexer &idxer, double kappa*/
  PUP::fromMem implP(impl_buf);
  CkArrayIndexOctreeIndexer idxer; implP|idxer;
  double kappa; implP|kappa;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("idxer");
  implDestP|idxer;
  if (implDestP.hasComments()) implDestP.comment("kappa");
  implDestP|kappa;
}

/* DEFS: void collate_bem_results(const CkCallback &cb);
 */
void CProxy_OpenCL_NodeGroup::collate_bem_results(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupBroadcast(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg, ckGetGroupID(),0);
}
 int CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3=0;
void CkIndex_OpenCL_NodeGroup::_call_collate_bem_results_marshall3(void* impl_msg,OpenCL_NodeGroup * impl_obj)
{
  if(CmiTryLock(impl_obj->__nodelock)) {
    impl_msg = CkCopyMsg(&impl_msg);
    CkSendMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3,impl_msg,CkMyNode(),impl_obj->CkGetNodeGroupID());
    return;
  }
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->collate_bem_results(cb);
  CmiUnlock(impl_obj->__nodelock);
}
void CkIndex_OpenCL_NodeGroup::_marshallmessagepup_collate_bem_results_marshall3(PUP::er &implDestP,void *impl_msg) {
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

/* DEFS: void precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
 */
void CProxy_OpenCL_NodeGroup::precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkArrayIndexOctreeIndexer &idxer, double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupBroadcast(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg, ckGetGroupID(),0);
}
 int CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4=0;
void CkIndex_OpenCL_NodeGroup::_call_precalc_bem_marshall4(void* impl_msg,OpenCL_NodeGroup * impl_obj)
{
  if(CmiTryLock(impl_obj->__nodelock)) {
    impl_msg = CkCopyMsg(&impl_msg);
    CkSendMsgNodeBranch(CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4,impl_msg,CkMyNode(),impl_obj->CkGetNodeGroupID());
    return;
  }
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkArrayIndexOctreeIndexer &idxer, double kappa*/
  PUP::fromMem implP(impl_buf);
  CkArrayIndexOctreeIndexer idxer; implP|idxer;
  double kappa; implP|kappa;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->precalc_bem(idxer, kappa);
  CmiUnlock(impl_obj->__nodelock);
}
void CkIndex_OpenCL_NodeGroup::_marshallmessagepup_precalc_bem_marshall4(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkArrayIndexOctreeIndexer &idxer, double kappa*/
  PUP::fromMem implP(impl_buf);
  CkArrayIndexOctreeIndexer idxer; implP|idxer;
  double kappa; implP|kappa;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("idxer");
  implDestP|idxer;
  if (implDestP.hasComments()) implDestP.comment("kappa");
  implDestP|kappa;
}

/* DEFS: OpenCL_NodeGroup(unsigned int num_patches);
 */

/* DEFS: void run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
 */
void CProxySection_OpenCL_NodeGroup::run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkArrayIndexOctreeIndexer &idxer, double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSectionSend(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgNodeBranchMulti(CkIndex_OpenCL_NodeGroup::__idx_run_bem_marshall2, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

/* DEFS: void collate_bem_results(const CkCallback &cb);
 */
void CProxySection_OpenCL_NodeGroup::collate_bem_results(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSectionSend(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgNodeBranchMulti(CkIndex_OpenCL_NodeGroup::__idx_collate_bem_results_marshall3, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

/* DEFS: void precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
 */
void CProxySection_OpenCL_NodeGroup::precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CkArrayIndexOctreeIndexer &idxer, double kappa
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CkArrayIndexOctreeIndexer &)idxer;
    implP|kappa;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSectionSend(ckDelegatedPtr(),CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgNodeBranchMulti(CkIndex_OpenCL_NodeGroup::__idx_precalc_bem_marshall4, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_OpenCL_NodeGroup::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_NodeGroup::__idx);
   CkRegisterGroupIrr(__idx,OpenCL_NodeGroup::isIrreducible());
// REG: OpenCL_NodeGroup(unsigned int num_patches);
  __idx_OpenCL_NodeGroup_marshall1 = CkRegisterEp("OpenCL_NodeGroup(unsigned int num_patches)",
     (CkCallFnPtr)_call_OpenCL_NodeGroup_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_OpenCL_NodeGroup_marshall1,(CkMarshallUnpackFn)_callmarshall_OpenCL_NodeGroup_marshall1);
  CkRegisterMessagePupFn(__idx_OpenCL_NodeGroup_marshall1,(CkMessagePupFn)_marshallmessagepup_OpenCL_NodeGroup_marshall1);

// REG: void run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
  __idx_run_bem_marshall2 = CkRegisterEp("run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa)",
     (CkCallFnPtr)_call_run_bem_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMessagePupFn(__idx_run_bem_marshall2,(CkMessagePupFn)_marshallmessagepup_run_bem_marshall2);

// REG: void collate_bem_results(const CkCallback &cb);
  __idx_collate_bem_results_marshall3 = CkRegisterEp("collate_bem_results(const CkCallback &cb)",
     (CkCallFnPtr)_call_collate_bem_results_marshall3, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMessagePupFn(__idx_collate_bem_results_marshall3,(CkMessagePupFn)_marshallmessagepup_collate_bem_results_marshall3);

// REG: void precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
  __idx_precalc_bem_marshall4 = CkRegisterEp("precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa)",
     (CkCallFnPtr)_call_precalc_bem_marshall4, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMessagePupFn(__idx_precalc_bem_marshall4,(CkMessagePupFn)_marshallmessagepup_precalc_bem_marshall4);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registeropencl_nodegroup(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: nodegroup OpenCL_NodeGroup: NodeGroup{
OpenCL_NodeGroup(unsigned int num_patches);
void run_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
void collate_bem_results(const CkCallback &cb);
void precalc_bem(const CkArrayIndexOctreeIndexer &idxer, double kappa);
};
*/
  CkIndex_OpenCL_NodeGroup::__register("OpenCL_NodeGroup", sizeof(OpenCL_NodeGroup));

}
#endif
