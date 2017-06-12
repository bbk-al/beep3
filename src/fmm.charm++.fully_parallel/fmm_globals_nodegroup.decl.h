#ifndef _DECL_fmm_globals_nodegroup_H_
#define _DECL_fmm_globals_nodegroup_H_
#include "charm++.h"
/* DECLS: template < int NTERMS, int NLAMBS > nodegroup FMM_Globals_NodeGroupT: NodeGroup{
FMM_Globals_NodeGroupT(void);
};
 */
template < int NTERMS, int NLAMBS >  class FMM_Globals_NodeGroupT;
template < int NTERMS, int NLAMBS >  class CkIndex_FMM_Globals_NodeGroupT;
template < int NTERMS, int NLAMBS >  class CProxy_FMM_Globals_NodeGroupT;
template < int NTERMS, int NLAMBS >  class CProxyElement_FMM_Globals_NodeGroupT;
template < int NTERMS, int NLAMBS >  class CProxySection_FMM_Globals_NodeGroupT;
/* --------------- index object ------------------ */
template < int NTERMS, int NLAMBS > class CkIndex_FMM_Globals_NodeGroupT:public CProxyElement_NodeGroup{
  public:
    typedef FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  local_t;
    typedef CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  index_t;
    typedef CProxy_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  proxy_t;
    typedef CProxyElement_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  element_t;
    typedef CProxySection_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: FMM_Globals_NodeGroupT(void);
 */
    static int __idx_FMM_Globals_NodeGroupT_void;
    static int ckNew(void) { return __idx_FMM_Globals_NodeGroupT_void; }
    static void _call_FMM_Globals_NodeGroupT_void(void* impl_msg,FMM_Globals_NodeGroupT < NTERMS, NLAMBS > * impl_obj);

};
/* --------------- element proxy ------------------ */
template < int NTERMS, int NLAMBS > class CProxyElement_FMM_Globals_NodeGroupT: public CProxyElement_NodeGroup{
  public:
    typedef FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  local_t;
    typedef CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  index_t;
    typedef CProxy_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  proxy_t;
    typedef CProxyElement_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  element_t;
    typedef CProxySection_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  section_t;

    CProxyElement_FMM_Globals_NodeGroupT(void) {}
    CProxyElement_FMM_Globals_NodeGroupT(const IrrGroup *g) : CProxyElement_NodeGroup(g){  }
    CProxyElement_FMM_Globals_NodeGroupT(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_NodeGroup(_gid,_onPE,CK_DELCTOR_ARGS){  }
    CProxyElement_FMM_Globals_NodeGroupT(CkGroupID _gid,int _onPE) : CProxyElement_NodeGroup(_gid,_onPE){  }
int ckIsDelegated(void) const {return CProxyElement_NodeGroup::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxyElement_NodeGroup::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxyElement_NodeGroup::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxyElement_NodeGroup::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxyElement_NodeGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxyElement_NodeGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxyElement_NodeGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxyElement_NodeGroup::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxyElement_NodeGroup::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxyElement_NodeGroup::ckSetReductionClient(cb); }
int ckGetGroupPe(void) const
{return CProxyElement_NodeGroup::ckGetGroupPe();}
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxyElement_NodeGroup::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxyElement_NodeGroup::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxyElement_NodeGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxyElement_NodeGroup::ckSetGroupID(g);
    }
    FMM_Globals_NodeGroupT < NTERMS, NLAMBS > * ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static FMM_Globals_NodeGroupT < NTERMS, NLAMBS > * ckLocalBranch(CkGroupID gID) {
      return (FMM_Globals_NodeGroupT < NTERMS, NLAMBS > *)CkLocalNodeBranch(gID);
    }
/* DECLS: FMM_Globals_NodeGroupT(void);
 */

};
/* ---------------- collective proxy -------------- */
template < int NTERMS, int NLAMBS > class CProxy_FMM_Globals_NodeGroupT: public CProxy_NodeGroup{
  public:
    typedef FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  local_t;
    typedef CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  index_t;
    typedef CProxy_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  proxy_t;
    typedef CProxyElement_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  element_t;
    typedef CProxySection_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  section_t;

    CProxy_FMM_Globals_NodeGroupT(void) {}
    CProxy_FMM_Globals_NodeGroupT(const IrrGroup *g) : CProxy_NodeGroup(g){  }
    CProxy_FMM_Globals_NodeGroupT(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_NodeGroup(_gid,CK_DELCTOR_ARGS){  }
    CProxy_FMM_Globals_NodeGroupT(CkGroupID _gid) : CProxy_NodeGroup(_gid){  }
    CProxyElement_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  operator[](int onPE) const
      {return CProxyElement_FMM_Globals_NodeGroupT < NTERMS, NLAMBS > (ckGetGroupID(),onPE,CK_DELCTOR_CALL);}
int ckIsDelegated(void) const {return CProxy_NodeGroup::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxy_NodeGroup::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxy_NodeGroup::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxy_NodeGroup::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxy_NodeGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxy_NodeGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxy_NodeGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxy_NodeGroup::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxy_NodeGroup::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxy_NodeGroup::ckSetReductionClient(cb); }
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxy_NodeGroup::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxy_NodeGroup::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_NodeGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxy_NodeGroup::ckSetGroupID(g);
    }
    FMM_Globals_NodeGroupT < NTERMS, NLAMBS > * ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static FMM_Globals_NodeGroupT < NTERMS, NLAMBS > * ckLocalBranch(CkGroupID gID) {
      return (FMM_Globals_NodeGroupT < NTERMS, NLAMBS > *)CkLocalNodeBranch(gID);
    }
/* DECLS: FMM_Globals_NodeGroupT(void);
 */
    static CkGroupID ckNew(void);

};
/* ---------------- section proxy -------------- */
template < int NTERMS, int NLAMBS > class CProxySection_FMM_Globals_NodeGroupT: public CProxySection_NodeGroup{
  public:
    typedef FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  local_t;
    typedef CkIndex_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  index_t;
    typedef CProxy_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  proxy_t;
    typedef CProxyElement_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  element_t;
    typedef CProxySection_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  section_t;

    CProxySection_FMM_Globals_NodeGroupT(void) {}
    CProxySection_FMM_Globals_NodeGroupT(const IrrGroup *g) : CProxySection_NodeGroup(g){  }
    CProxySection_FMM_Globals_NodeGroupT(const CkGroupID &_gid,const int *_pelist,int _npes,CK_DELCTOR_PARAM) : CProxySection_NodeGroup(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_FMM_Globals_NodeGroupT(const CkGroupID &_gid,const int *_pelist,int _npes) : CProxySection_NodeGroup(_gid,_pelist,_npes){  }
    CProxySection_FMM_Globals_NodeGroupT(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes) : CProxySection_NodeGroup(n,_gid,_pelist,_npes){  }
    CProxySection_FMM_Globals_NodeGroupT(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes,CK_DELCTOR_PARAM) : CProxySection_NodeGroup(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
int ckIsDelegated(void) const {return CProxySection_NodeGroup::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxySection_NodeGroup::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxySection_NodeGroup::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxySection_NodeGroup::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxySection_NodeGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxySection_NodeGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxySection_NodeGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxySection_NodeGroup::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxySection_NodeGroup::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxySection_NodeGroup::ckSetReductionClient(cb); }
inline int ckGetNumSections() const
{ return CProxySection_NodeGroup::ckGetNumSections(); }
inline CkSectionInfo &ckGetSectionInfo()
{ return CProxySection_NodeGroup::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
{ return CProxySection_NodeGroup::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
{ return CProxySection_NodeGroup::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
{ return CProxySection_NodeGroup::ckGetSectionID(i); }
inline CkGroupID ckGetGroupIDn(int i) const
{ return CProxySection_NodeGroup::ckGetGroupIDn(i); }
inline int *ckGetElements() const
{ return CProxySection_NodeGroup::ckGetElements(); }
inline int *ckGetElements(int i) const
{ return CProxySection_NodeGroup::ckGetElements(i); }
inline int ckGetNumElements() const
{ return CProxySection_NodeGroup::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{ return CProxySection_NodeGroup::ckGetNumElements(i); }
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxySection_NodeGroup::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxySection_NodeGroup::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxySection_NodeGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxySection_NodeGroup::ckSetGroupID(g);
    }
    FMM_Globals_NodeGroupT < NTERMS, NLAMBS > * ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static FMM_Globals_NodeGroupT < NTERMS, NLAMBS > * ckLocalBranch(CkGroupID gID) {
      return (FMM_Globals_NodeGroupT < NTERMS, NLAMBS > *)CkLocalNodeBranch(gID);
    }
/* DECLS: FMM_Globals_NodeGroupT(void);
 */

};
template < int NTERMS, int NLAMBS > 
class CBase_FMM_Globals_NodeGroupT : public CBaseT1<NodeGroup, CProxy_FMM_Globals_NodeGroupT < NTERMS, NLAMBS >  > { };





extern void _registerfmm_globals_nodegroup(void);
#endif
