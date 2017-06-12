#ifndef _DECL_opencl_nodegroup_H_
#define _DECL_opencl_nodegroup_H_
#include "charm++.h"
/* DECLS: nodegroup OpenCL_NodeGroup: NodeGroup{
OpenCL_NodeGroup(void);
};
 */
 class OpenCL_NodeGroup;
 class CkIndex_OpenCL_NodeGroup;
 class CProxy_OpenCL_NodeGroup;
 class CProxyElement_OpenCL_NodeGroup;
 class CProxySection_OpenCL_NodeGroup;
/* --------------- index object ------------------ */
class CkIndex_OpenCL_NodeGroup:public CProxyElement_NodeGroup{
  public:
    typedef OpenCL_NodeGroup local_t;
    typedef CkIndex_OpenCL_NodeGroup index_t;
    typedef CProxy_OpenCL_NodeGroup proxy_t;
    typedef CProxyElement_OpenCL_NodeGroup element_t;
    typedef CProxySection_OpenCL_NodeGroup section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: OpenCL_NodeGroup(void);
 */
    static int __idx_OpenCL_NodeGroup_void;
    static int ckNew(void) { return __idx_OpenCL_NodeGroup_void; }
    static void _call_OpenCL_NodeGroup_void(void* impl_msg,OpenCL_NodeGroup* impl_obj);

};
/* --------------- element proxy ------------------ */
class CProxyElement_OpenCL_NodeGroup: public CProxyElement_NodeGroup{
  public:
    typedef OpenCL_NodeGroup local_t;
    typedef CkIndex_OpenCL_NodeGroup index_t;
    typedef CProxy_OpenCL_NodeGroup proxy_t;
    typedef CProxyElement_OpenCL_NodeGroup element_t;
    typedef CProxySection_OpenCL_NodeGroup section_t;

    CProxyElement_OpenCL_NodeGroup(void) {}
    CProxyElement_OpenCL_NodeGroup(const IrrGroup *g) : CProxyElement_NodeGroup(g){  }
    CProxyElement_OpenCL_NodeGroup(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_NodeGroup(_gid,_onPE,CK_DELCTOR_ARGS){  }
    CProxyElement_OpenCL_NodeGroup(CkGroupID _gid,int _onPE) : CProxyElement_NodeGroup(_gid,_onPE){  }
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
    OpenCL_NodeGroup* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static OpenCL_NodeGroup* ckLocalBranch(CkGroupID gID) {
      return (OpenCL_NodeGroup*)CkLocalNodeBranch(gID);
    }
/* DECLS: OpenCL_NodeGroup(void);
 */

};
PUPmarshall(CProxyElement_OpenCL_NodeGroup)
/* ---------------- collective proxy -------------- */
class CProxy_OpenCL_NodeGroup: public CProxy_NodeGroup{
  public:
    typedef OpenCL_NodeGroup local_t;
    typedef CkIndex_OpenCL_NodeGroup index_t;
    typedef CProxy_OpenCL_NodeGroup proxy_t;
    typedef CProxyElement_OpenCL_NodeGroup element_t;
    typedef CProxySection_OpenCL_NodeGroup section_t;

    CProxy_OpenCL_NodeGroup(void) {}
    CProxy_OpenCL_NodeGroup(const IrrGroup *g) : CProxy_NodeGroup(g){  }
    CProxy_OpenCL_NodeGroup(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_NodeGroup(_gid,CK_DELCTOR_ARGS){  }
    CProxy_OpenCL_NodeGroup(CkGroupID _gid) : CProxy_NodeGroup(_gid){  }
    CProxyElement_OpenCL_NodeGroup operator[](int onPE) const
      {return CProxyElement_OpenCL_NodeGroup(ckGetGroupID(),onPE,CK_DELCTOR_CALL);}
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
    OpenCL_NodeGroup* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static OpenCL_NodeGroup* ckLocalBranch(CkGroupID gID) {
      return (OpenCL_NodeGroup*)CkLocalNodeBranch(gID);
    }
/* DECLS: OpenCL_NodeGroup(void);
 */
    static CkGroupID ckNew(void);

};
PUPmarshall(CProxy_OpenCL_NodeGroup)
/* ---------------- section proxy -------------- */
class CProxySection_OpenCL_NodeGroup: public CProxySection_NodeGroup{
  public:
    typedef OpenCL_NodeGroup local_t;
    typedef CkIndex_OpenCL_NodeGroup index_t;
    typedef CProxy_OpenCL_NodeGroup proxy_t;
    typedef CProxyElement_OpenCL_NodeGroup element_t;
    typedef CProxySection_OpenCL_NodeGroup section_t;

    CProxySection_OpenCL_NodeGroup(void) {}
    CProxySection_OpenCL_NodeGroup(const IrrGroup *g) : CProxySection_NodeGroup(g){  }
    CProxySection_OpenCL_NodeGroup(const CkGroupID &_gid,const int *_pelist,int _npes,CK_DELCTOR_PARAM) : CProxySection_NodeGroup(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_OpenCL_NodeGroup(const CkGroupID &_gid,const int *_pelist,int _npes) : CProxySection_NodeGroup(_gid,_pelist,_npes){  }
    CProxySection_OpenCL_NodeGroup(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes) : CProxySection_NodeGroup(n,_gid,_pelist,_npes){  }
    CProxySection_OpenCL_NodeGroup(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes,CK_DELCTOR_PARAM) : CProxySection_NodeGroup(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
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
    OpenCL_NodeGroup* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static OpenCL_NodeGroup* ckLocalBranch(CkGroupID gID) {
      return (OpenCL_NodeGroup*)CkLocalNodeBranch(gID);
    }
/* DECLS: OpenCL_NodeGroup(void);
 */

};
PUPmarshall(CProxySection_OpenCL_NodeGroup)
typedef CBaseT1<NodeGroup, CProxy_OpenCL_NodeGroup> CBase_OpenCL_NodeGroup;

extern void _registeropencl_nodegroup(void);
#endif
