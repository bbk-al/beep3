/* DEFS: chare GMRES: Chare{
GMRES(void);
threaded sync void solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_GMRES::__idx=0;
#endif
#ifndef CK_TEMPLATES_ONLY
/* DEFS: GMRES(void);
 */
CkChareID CProxy_GMRES::ckNew(int impl_onPE)
{
  void *impl_msg = CkAllocSysMsg();
  CkChareID impl_ret;
  CkCreateChare(CkIndex_GMRES::__idx, CkIndex_GMRES::__idx_GMRES_void, impl_msg, &impl_ret, impl_onPE);
  return impl_ret;
}
void CProxy_GMRES::ckNew(CkChareID* pcid, int impl_onPE)
{
  void *impl_msg = CkAllocSysMsg();
  CkCreateChare(CkIndex_GMRES::__idx, CkIndex_GMRES::__idx_GMRES_void, impl_msg, pcid, impl_onPE);
}
 int CkIndex_GMRES::__idx_GMRES_void=0;
void CkIndex_GMRES::_call_GMRES_void(void* impl_msg,GMRES * impl_obj)
{
  CkFreeSysMsg(impl_msg);
  new (impl_obj) GMRES();
}

/* DEFS: threaded sync void solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs);
 */
void CProxy_GMRES::solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)output_filename;
    implP|kappa;
    implP|Dsolvent;
    implP|num_patches;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)rhs;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(std::string &)output_filename;
    implP|kappa;
    implP|Dsolvent;
    implP|num_patches;
    //Have to cast away const-ness to get pup routine
    implP|(FH_Values &)rhs;
  }
  CkFreeSysMsg(CkRemoteCall(CkIndex_GMRES::__idx_solve_marshall2, impl_msg, &ckGetChareID()));
}
 int CkIndex_GMRES::__idx_solve_marshall2=0;
void CkIndex_GMRES::_call_solve_marshall2(void* impl_msg,GMRES * impl_obj)
{
  CthThread tid = CthCreate((CthVoidFn)_callthr_solve_marshall2, new CkThrCallArg(impl_msg,impl_obj), 0);
  ((Chare *)impl_obj)->CkAddThreadListeners(tid,impl_msg);
  CthAwaken(tid);
}
void CkIndex_GMRES::_callthr_solve_marshall2(CkThrCallArg *impl_arg)
{
  void *impl_msg = impl_arg->msg;
  GMRES *impl_obj = (GMRES *) impl_arg->obj;
  delete impl_arg;
  int impl_ref = CkGetRefNum(impl_msg), impl_src = CkGetSrcPe(impl_msg);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs*/
  PUP::fromMem implP(impl_buf);
  std::string output_filename; implP|output_filename;
  double kappa; implP|kappa;
  double Dsolvent; implP|Dsolvent;
  unsigned int num_patches; implP|num_patches;
  FH_Values rhs; implP|rhs;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  void *impl_retMsg=CkAllocSysMsg();
    impl_obj->solve(output_filename, kappa, Dsolvent, num_patches, rhs);
  CkSendToFutureID(impl_ref, impl_retMsg, impl_src);
  delete impl_msg_typed;
}
void CkIndex_GMRES::_marshallmessagepup_solve_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  /*Unmarshall pup'd fields: const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs*/
  PUP::fromMem implP(impl_buf);
  std::string output_filename; implP|output_filename;
  double kappa; implP|kappa;
  double Dsolvent; implP|Dsolvent;
  unsigned int num_patches; implP|num_patches;
  FH_Values rhs; implP|rhs;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("output_filename");
  implDestP|output_filename;
  if (implDestP.hasComments()) implDestP.comment("kappa");
  implDestP|kappa;
  if (implDestP.hasComments()) implDestP.comment("Dsolvent");
  implDestP|Dsolvent;
  if (implDestP.hasComments()) implDestP.comment("num_patches");
  implDestP|num_patches;
  if (implDestP.hasComments()) implDestP.comment("rhs");
  implDestP|rhs;
}

#endif /*CK_TEMPLATES_ONLY*/
#ifndef CK_TEMPLATES_ONLY
void CkIndex_GMRES::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeChare);
  CkRegisterBase(__idx, CkIndex_Chare::__idx);
// REG: GMRES(void);
  __idx_GMRES_void = CkRegisterEp("GMRES(void)",
     (CkCallFnPtr)_call_GMRES_void, 0, __idx, 0);
  CkRegisterDefaultCtor(__idx, __idx_GMRES_void);

// REG: threaded sync void solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs);
  __idx_solve_marshall2 = CkRegisterEp("solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs)",
     (CkCallFnPtr)_call_solve_marshall2, CkMarshallMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(__idx_solve_marshall2,(CkMessagePupFn)_marshallmessagepup_solve_marshall2);

}
#endif

#ifndef CK_TEMPLATES_ONLY
void _registergmres(void)
{
  static int _done = 0; if(_done) return; _done = 1;
/* REG: chare GMRES: Chare{
GMRES(void);
threaded sync void solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs);
};
*/
  CkIndex_GMRES::__register("GMRES", sizeof(GMRES));

}
#endif
