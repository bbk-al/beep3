/* DEFS: template < typename CType, typename CProxy_FMMWorkerT > nodegroup ParallelFMMOctree: NodeGroup{
ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
void insert(const std::vector<CType > &items);
void insert(const CType &single);
void finalize(const CkCallback &cb);
void clear_waves(void);
void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
};
 */
#ifdef CK_TEMPLATES_ONLY
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx=0;
#endif
#ifdef CK_TEMPLATES_ONLY
/* DEFS: ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
 */

/* DEFS: void insert(const std::vector<CType > &items);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::insert(const std::vector<CType > &items, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const std::vector<CType > &items
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<CType > &)items;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<CType > &)items;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: void insert(const CType &single);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::insert(const CType &single, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CType &single
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CType &)single;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CType &)single;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: void finalize(const CkCallback &cb);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::finalize(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: void clear_waves(void);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::clear_waves(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
}

/* DEFS: void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::request_data(const CProxy_FMMWorkerT &FMMWorkerProxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_FMMWorkerT &FMMWorkerProxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMMWorkerT &)FMMWorkerProxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMMWorkerT &)FMMWorkerProxy;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else CkSendMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg, ckGetGroupPe(), ckGetGroupID(),0+CK_MSG_IMMEDIATE);
}

/* DEFS: ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
 */
template < typename CType, typename CProxy_FMMWorkerT >  CkGroupID CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::ckNew(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_1;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_1;
  }
  return CkCreateNodeGroup(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx, CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_ParallelFMMOctree_marshall1, impl_msg);
}
template < typename CType, typename CProxy_FMMWorkerT >    CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::CProxy_ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1, const CkEntryOptions *impl_e_opts)
{
  //Marshall: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_1;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|max_items_per_node;
    //Have to cast away const-ness to get pup routine
    implP|(Vector &)universe_centre;
    implP|universe_edge_length;
    implP|impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(CkCallback &)impl_noname_1;
  }
  ckSetGroupID(CkCreateNodeGroup(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx, CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_ParallelFMMOctree_marshall1, impl_msg));
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_ParallelFMMOctree_marshall1=0;
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_call_ParallelFMMOctree_marshall1(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1*/
  PUP::fromMem implP(impl_buf);
  unsigned int max_items_per_node; implP|max_items_per_node;
  Vector universe_centre; implP|universe_centre;
  double universe_edge_length; implP|universe_edge_length;
  unsigned int impl_noname_0; implP|impl_noname_0;
  CkCallback impl_noname_1; implP|impl_noname_1;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) ParallelFMMOctree < CType, CProxy_FMMWorkerT > (max_items_per_node, universe_centre, universe_edge_length, impl_noname_0, impl_noname_1);
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_callmarshall_ParallelFMMOctree_marshall1(char* impl_buf,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj) {
  /*Unmarshall pup'd fields: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1*/
  PUP::fromMem implP(impl_buf);
  unsigned int max_items_per_node; implP|max_items_per_node;
  Vector universe_centre; implP|universe_centre;
  double universe_edge_length; implP|universe_edge_length;
  unsigned int impl_noname_0; implP|impl_noname_0;
  CkCallback impl_noname_1; implP|impl_noname_1;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) ParallelFMMOctree < CType, CProxy_FMMWorkerT > (max_items_per_node, universe_centre, universe_edge_length, impl_noname_0, impl_noname_1);
  return implP.size();
}
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_marshallmessagepup_ParallelFMMOctree_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1*/
  PUP::fromMem implP(impl_buf);
  unsigned int max_items_per_node; implP|max_items_per_node;
  Vector universe_centre; implP|universe_centre;
  double universe_edge_length; implP|universe_edge_length;
  unsigned int impl_noname_0; implP|impl_noname_0;
  CkCallback impl_noname_1; implP|impl_noname_1;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("max_items_per_node");
  implDestP|max_items_per_node;
  if (implDestP.hasComments()) implDestP.comment("universe_centre");
  implDestP|universe_centre;
  if (implDestP.hasComments()) implDestP.comment("universe_edge_length");
  implDestP|universe_edge_length;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_1");
  implDestP|impl_noname_1;
}

/* DEFS: void insert(const std::vector<CType > &items);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::insert(const std::vector<CType > &items, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const std::vector<CType > &items
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<CType > &)items;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<CType > &)items;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupBroadcast(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg, ckGetGroupID(),0);
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2=0;
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_call_insert_marshall2(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj)
{
  if(CmiTryLock(impl_obj->__nodelock)) {
    impl_msg = CkCopyMsg(&impl_msg);
    CkSendMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2,impl_msg,CkMyNode(),impl_obj->CkGetNodeGroupID());
    return;
  }
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::vector<CType > &items*/
  PUP::fromMem implP(impl_buf);
  std::vector<CType > items; implP|items;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->insert(items);
  CmiUnlock(impl_obj->__nodelock);
}
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_marshallmessagepup_insert_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::vector<CType > &items*/
  PUP::fromMem implP(impl_buf);
  std::vector<CType > items; implP|items;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("items");
  implDestP|items;
}

/* DEFS: void insert(const CType &single);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::insert(const CType &single, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CType &single
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CType &)single;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CType &)single;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupBroadcast(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg, ckGetGroupID(),0);
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3=0;
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_call_insert_marshall3(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj)
{
  if(CmiTryLock(impl_obj->__nodelock)) {
    impl_msg = CkCopyMsg(&impl_msg);
    CkSendMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3,impl_msg,CkMyNode(),impl_obj->CkGetNodeGroupID());
    return;
  }
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CType &single*/
  PUP::fromMem implP(impl_buf);
  CType single; implP|single;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->insert(single);
  CmiUnlock(impl_obj->__nodelock);
}
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_marshallmessagepup_insert_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CType &single*/
  PUP::fromMem implP(impl_buf);
  CType single; implP|single;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("single");
  implDestP|single;
}

/* DEFS: void finalize(const CkCallback &cb);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::finalize(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupBroadcast(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg, ckGetGroupID(),0);
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4=0;
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_call_finalize_marshall4(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->finalize(cb);
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_callmarshall_finalize_marshall4(char* impl_buf,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj) {
  /*Unmarshall pup'd fields: const CkCallback &cb*/
  PUP::fromMem implP(impl_buf);
  CkCallback cb; implP|cb;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->finalize(cb);
  return implP.size();
}
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_marshallmessagepup_finalize_marshall4(PUP::er &implDestP,void *impl_msg) {
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

/* DEFS: void clear_waves(void);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::clear_waves(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupBroadcast(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg, ckGetGroupID(),0);
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void=0;
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_call_clear_waves_void(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  impl_obj->clear_waves();
}

/* DEFS: void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::request_data(const CProxy_FMMWorkerT &FMMWorkerProxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_FMMWorkerT &FMMWorkerProxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMMWorkerT &)FMMWorkerProxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMMWorkerT &)FMMWorkerProxy;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupBroadcast(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg, ckGetGroupID());
  } else CkBroadcastMsgNodeBranch(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg, ckGetGroupID(),0+CK_MSG_IMMEDIATE);
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6=0;
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_call_request_data_marshall6(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CProxy_FMMWorkerT &FMMWorkerProxy*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMMWorkerT FMMWorkerProxy; implP|FMMWorkerProxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->request_data(FMMWorkerProxy);
}
template < typename CType, typename CProxy_FMMWorkerT >  int CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_callmarshall_request_data_marshall6(char* impl_buf,ParallelFMMOctree < CType, CProxy_FMMWorkerT >  * impl_obj) {
  /*Unmarshall pup'd fields: const CProxy_FMMWorkerT &FMMWorkerProxy*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMMWorkerT FMMWorkerProxy; implP|FMMWorkerProxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->request_data(FMMWorkerProxy);
  return implP.size();
}
template < typename CType, typename CProxy_FMMWorkerT >  void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::_marshallmessagepup_request_data_marshall6(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const CProxy_FMMWorkerT &FMMWorkerProxy*/
  PUP::fromMem implP(impl_buf);
  CProxy_FMMWorkerT FMMWorkerProxy; implP|FMMWorkerProxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("FMMWorkerProxy");
  implDestP|FMMWorkerProxy;
}

/* DEFS: ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
 */

/* DEFS: void insert(const std::vector<CType > &items);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::insert(const std::vector<CType > &items, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const std::vector<CType > &items
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<CType > &)items;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<CType > &)items;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSectionSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgNodeBranchMulti(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall2, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

/* DEFS: void insert(const CType &single);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::insert(const CType &single, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CType &single
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CType &)single;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CType &)single;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSectionSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgNodeBranchMulti(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_insert_marshall3, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

/* DEFS: void finalize(const CkCallback &cb);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::finalize(const CkCallback &cb, const CkEntryOptions *impl_e_opts)
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
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSectionSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgNodeBranchMulti(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_finalize_marshall4, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

/* DEFS: void clear_waves(void);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::clear_waves(void)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg();
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSectionSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgNodeBranchMulti(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_clear_waves_void, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}

/* DEFS: void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
 */
template < typename CType, typename CProxy_FMMWorkerT >  void CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::request_data(const CProxy_FMMWorkerT &FMMWorkerProxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_FMMWorkerT &FMMWorkerProxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMMWorkerT &)FMMWorkerProxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(CProxy_FMMWorkerT &)FMMWorkerProxy;
  }
  if (ckIsDelegated()) {
     CkNodeGroupMsgPrep(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg, ckGetGroupID());
     ckDelegatedTo()->NodeGroupSectionSend(ckDelegatedPtr(),CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp = (ckGetNumSections()>1) ? CkCopyMsg((void **) &impl_msg) : impl_msg;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgNodeBranchMulti(CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__idx_request_data_marshall6, impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0+CK_MSG_IMMEDIATE);
    }
  }
}

#endif /*CK_TEMPLATES_ONLY*/
#ifdef CK_TEMPLATES_ONLY
template < typename CType, typename CProxy_FMMWorkerT > void CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_NodeGroup::__idx);
   CkRegisterGroupIrr(__idx,ParallelFMMOctree < CType, CProxy_FMMWorkerT > ::isIrreducible());
// REG: ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
  __idx_ParallelFMMOctree_marshall1 = CkRegisterEp("ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1)",
     (CkCallFnPtr)_call_ParallelFMMOctree_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_ParallelFMMOctree_marshall1,(CkMarshallUnpackFn)_callmarshall_ParallelFMMOctree_marshall1);
  CkRegisterMessagePupFn(__idx_ParallelFMMOctree_marshall1,(CkMessagePupFn)_marshallmessagepup_ParallelFMMOctree_marshall1);

// REG: void insert(const std::vector<CType > &items);
  __idx_insert_marshall2 = CkRegisterEp("insert(const std::vector<CType > &items)",
     (CkCallFnPtr)_call_insert_marshall2, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMessagePupFn(__idx_insert_marshall2,(CkMessagePupFn)_marshallmessagepup_insert_marshall2);

// REG: void insert(const CType &single);
  __idx_insert_marshall3 = CkRegisterEp("insert(const CType &single)",
     (CkCallFnPtr)_call_insert_marshall3, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMessagePupFn(__idx_insert_marshall3,(CkMessagePupFn)_marshallmessagepup_insert_marshall3);

// REG: void finalize(const CkCallback &cb);
  __idx_finalize_marshall4 = CkRegisterEp("finalize(const CkCallback &cb)",
     (CkCallFnPtr)_call_finalize_marshall4, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_finalize_marshall4,(CkMarshallUnpackFn)_callmarshall_finalize_marshall4);
  CkRegisterMessagePupFn(__idx_finalize_marshall4,(CkMessagePupFn)_marshallmessagepup_finalize_marshall4);

// REG: void clear_waves(void);
  __idx_clear_waves_void = CkRegisterEp("clear_waves(void)",
     (CkCallFnPtr)_call_clear_waves_void, 0, __idx, 0);

// REG: void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
  __idx_request_data_marshall6 = CkRegisterEp("request_data(const CProxy_FMMWorkerT &FMMWorkerProxy)",
     (CkCallFnPtr)_call_request_data_marshall6, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP+CK_EP_TRACEDISABLE);
  CkRegisterMarshallUnpackFn(__idx_request_data_marshall6,(CkMarshallUnpackFn)_callmarshall_request_data_marshall6);
  CkRegisterMessagePupFn(__idx_request_data_marshall6,(CkMessagePupFn)_marshallmessagepup_request_data_marshall6);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registerparallel_fmm_octree(void)
{
  static int _done = 0; if(_done) return; _done = 1;

}
#endif
