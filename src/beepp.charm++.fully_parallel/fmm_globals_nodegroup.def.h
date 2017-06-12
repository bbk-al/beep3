/* DEFS: template < int NTERMS, int NLAMBS > nodegroup FMM_Globals_NodeGroupT: NodeGroup{
FMM_Globals_NodeGroupT(void);
};
 */
#ifdef CK_TEMPLATES_ONLY
template < int NTERMS, int NLAMBS >  int CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ::__idx=0;
#endif
#ifdef CK_TEMPLATES_ONLY
/* DEFS: FMM_Globals_NodeGroupT(void);
 */

/* DEFS: FMM_Globals_NodeGroupT(void);
 */
template < int NTERMS, int NLAMBS >  CkGroupID CProxy_FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ::ckNew(void)
{
  void *impl_msg = CkAllocSysMsg();
  return CkCreateNodeGroup(CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ::__idx, CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ::__idx_FMM_Globals_NodeGroupT_void, impl_msg);
}
template < int NTERMS, int NLAMBS >  int CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ::__idx_FMM_Globals_NodeGroupT_void=0;
template < int NTERMS, int NLAMBS >  void CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ::_call_FMM_Globals_NodeGroupT_void(void* impl_msg,FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  new (impl_obj) FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ();
}

/* DEFS: FMM_Globals_NodeGroupT(void);
 */

#endif /*CK_TEMPLATES_ONLY*/
#ifdef CK_TEMPLATES_ONLY
template < int NTERMS, int NLAMBS > void CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_NodeGroup::__idx);
   CkRegisterGroupIrr(__idx,FMM_Globals_NodeGroupT < NTERMS, NLAMBS > ::isIrreducible());
// REG: FMM_Globals_NodeGroupT(void);
  __idx_FMM_Globals_NodeGroupT_void = CkRegisterEp("FMM_Globals_NodeGroupT(void)",
     (CkCallFnPtr)_call_FMM_Globals_NodeGroupT_void, 0, __idx, 0);
  CkRegisterDefaultCtor(__idx, __idx_FMM_Globals_NodeGroupT_void);

}
#endif

/* DEFS: nodegroup FMM_Globals_NodeGroupT<2,2 >: NodeGroup;
 */

/* DEFS: nodegroup FMM_Globals_NodeGroupT<4,4 >: NodeGroup;
 */

/* DEFS: nodegroup FMM_Globals_NodeGroupT<9,9 >: NodeGroup;
 */

/* DEFS: nodegroup FMM_Globals_NodeGroupT<18,18 >: NodeGroup;
 */

#ifndef CK_TEMPLATES_ONLY
void _registerfmm_globals_nodegroup(void)
{
  static int _done = 0; if(_done) return; _done = 1;

/* REG: nodegroup FMM_Globals_NodeGroupT<2,2 >: NodeGroup;
*/
  CkIndex_FMM_Globals_NodeGroupT<2,2 >::__register("FMM_Globals_NodeGroupT<2,2 >", sizeof(FMM_Globals_NodeGroupT<2,2 >));

/* REG: nodegroup FMM_Globals_NodeGroupT<4,4 >: NodeGroup;
*/
  CkIndex_FMM_Globals_NodeGroupT<4,4 >::__register("FMM_Globals_NodeGroupT<4,4 >", sizeof(FMM_Globals_NodeGroupT<4,4 >));

/* REG: nodegroup FMM_Globals_NodeGroupT<9,9 >: NodeGroup;
*/
  CkIndex_FMM_Globals_NodeGroupT<9,9 >::__register("FMM_Globals_NodeGroupT<9,9 >", sizeof(FMM_Globals_NodeGroupT<9,9 >));

/* REG: nodegroup FMM_Globals_NodeGroupT<18,18 >: NodeGroup;
*/
  CkIndex_FMM_Globals_NodeGroupT<18,18 >::__register("FMM_Globals_NodeGroupT<18,18 >", sizeof(FMM_Globals_NodeGroupT<18,18 >));

}
#endif
