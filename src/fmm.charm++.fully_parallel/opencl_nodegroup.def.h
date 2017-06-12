/* DEFS: nodegroup OpenCL_NodeGroup: NodeGroup{
OpenCL_NodeGroup(void);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_OpenCL_NodeGroup::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: OpenCL_NodeGroup(void);
 */

/* DEFS: OpenCL_NodeGroup(void);
 */
CkGroupID CProxy_OpenCL_NodeGroup::ckNew(void)
{
  void *impl_msg = CkAllocSysMsg();
  return CkCreateNodeGroup(CkIndex_OpenCL_NodeGroup::__idx, CkIndex_OpenCL_NodeGroup::__idx_OpenCL_NodeGroup_void, impl_msg);
}
 int CkIndex_OpenCL_NodeGroup::__idx_OpenCL_NodeGroup_void=0;
void CkIndex_OpenCL_NodeGroup::_call_OpenCL_NodeGroup_void(void* impl_msg,OpenCL_NodeGroup * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  new (impl_obj) OpenCL_NodeGroup();
}

/* DEFS: OpenCL_NodeGroup(void);
 */

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_OpenCL_NodeGroup::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_NodeGroup::__idx);
   CkRegisterGroupIrr(__idx,OpenCL_NodeGroup::isIrreducible());
// REG: OpenCL_NodeGroup(void);
  __idx_OpenCL_NodeGroup_void = CkRegisterEp("OpenCL_NodeGroup(void)",
     (CkCallFnPtr)_call_OpenCL_NodeGroup_void, 0, __idx, 0);
  CkRegisterDefaultCtor(__idx, __idx_OpenCL_NodeGroup_void);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registeropencl_nodegroup(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: nodegroup OpenCL_NodeGroup: NodeGroup{
OpenCL_NodeGroup(void);
};
*/
  CkIndex_OpenCL_NodeGroup::__register("OpenCL_NodeGroup", sizeof(OpenCL_NodeGroup));

}
#endif
