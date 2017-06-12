/* DEFS: nodegroup MeshLibrary: NodeGroup{
MeshLibrary(const std::vector<std::string > &filenames);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_MeshLibrary::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: MeshLibrary(const std::vector<std::string > &filenames);
 */

/* DEFS: MeshLibrary(const std::vector<std::string > &filenames);
 */
CkGroupID CProxy_MeshLibrary::ckNew(const std::vector<std::string > &filenames, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const std::vector<std::string > &filenames
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<std::string > &)filenames;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<std::string > &)filenames;
  }
  return CkCreateNodeGroup(CkIndex_MeshLibrary::__idx, CkIndex_MeshLibrary::__idx_MeshLibrary_marshall1, impl_msg);
}
  CProxy_MeshLibrary::CProxy_MeshLibrary(const std::vector<std::string > &filenames, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const std::vector<std::string > &filenames
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<std::string > &)filenames;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::vector<std::string > &)filenames;
  }
  ckSetGroupID(CkCreateNodeGroup(CkIndex_MeshLibrary::__idx, CkIndex_MeshLibrary::__idx_MeshLibrary_marshall1, impl_msg));
}
 int CkIndex_MeshLibrary::__idx_MeshLibrary_marshall1=0;
void CkIndex_MeshLibrary::_call_MeshLibrary_marshall1(void* impl_msg,MeshLibrary * impl_obj)
{
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::vector<std::string > &filenames*/
  PUP::fromMem implP(impl_buf);
  std::vector<std::string > filenames; implP|filenames;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) MeshLibrary(filenames);
}
int CkIndex_MeshLibrary::_callmarshall_MeshLibrary_marshall1(char* impl_buf,MeshLibrary * impl_obj) {
  /*Unmarshall pup'd fields: const std::vector<std::string > &filenames*/
  PUP::fromMem implP(impl_buf);
  std::vector<std::string > filenames; implP|filenames;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj) MeshLibrary(filenames);
  return implP.size();
}
void CkIndex_MeshLibrary::_marshallmessagepup_MeshLibrary_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::vector<std::string > &filenames*/
  PUP::fromMem implP(impl_buf);
  std::vector<std::string > filenames; implP|filenames;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("filenames");
  implDestP|filenames;
}

/* DEFS: MeshLibrary(const std::vector<std::string > &filenames);
 */

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_MeshLibrary::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_NodeGroup::__idx);
   CkRegisterGroupIrr(__idx,MeshLibrary::isIrreducible());
// REG: MeshLibrary(const std::vector<std::string > &filenames);
  __idx_MeshLibrary_marshall1 = CkRegisterEp("MeshLibrary(const std::vector<std::string > &filenames)",
     (CkCallFnPtr)_call_MeshLibrary_marshall1, CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(__idx_MeshLibrary_marshall1,(CkMarshallUnpackFn)_callmarshall_MeshLibrary_marshall1);
  CkRegisterMessagePupFn(__idx_MeshLibrary_marshall1,(CkMessagePupFn)_marshallmessagepup_MeshLibrary_marshall1);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registermesh_library(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: nodegroup MeshLibrary: NodeGroup{
MeshLibrary(const std::vector<std::string > &filenames);
};
*/
  CkIndex_MeshLibrary::__register("MeshLibrary", sizeof(MeshLibrary));

}
#endif
