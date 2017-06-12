#ifndef _DECL_mesh_library_H_
#define _DECL_mesh_library_H_
#include "charm++.h"
/* DECLS: nodegroup MeshLibrary: NodeGroup{
MeshLibrary(const std::vector<std::string > &filenames);
};
 */
 class MeshLibrary;
 class CkIndex_MeshLibrary;
 class CProxy_MeshLibrary;
 class CProxyElement_MeshLibrary;
 class CProxySection_MeshLibrary;
/* --------------- index object ------------------ */
class CkIndex_MeshLibrary:public CProxyElement_NodeGroup{
  public:
    typedef MeshLibrary local_t;
    typedef CkIndex_MeshLibrary index_t;
    typedef CProxy_MeshLibrary proxy_t;
    typedef CProxyElement_MeshLibrary element_t;
    typedef CProxySection_MeshLibrary section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: MeshLibrary(const std::vector<std::string > &filenames);
 */
    static int __idx_MeshLibrary_marshall1;
    static int ckNew(const std::vector<std::string > &filenames) { return __idx_MeshLibrary_marshall1; }
    static void _call_MeshLibrary_marshall1(void* impl_msg,MeshLibrary* impl_obj);
    static int _callmarshall_MeshLibrary_marshall1(char* impl_buf,MeshLibrary* impl_obj);
    static void _marshallmessagepup_MeshLibrary_marshall1(PUP::er &p,void *msg);

};
/* --------------- element proxy ------------------ */
class CProxyElement_MeshLibrary: public CProxyElement_NodeGroup{
  public:
    typedef MeshLibrary local_t;
    typedef CkIndex_MeshLibrary index_t;
    typedef CProxy_MeshLibrary proxy_t;
    typedef CProxyElement_MeshLibrary element_t;
    typedef CProxySection_MeshLibrary section_t;

    CProxyElement_MeshLibrary(void) {}
    CProxyElement_MeshLibrary(const IrrGroup *g) : CProxyElement_NodeGroup(g){  }
    CProxyElement_MeshLibrary(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_NodeGroup(_gid,_onPE,CK_DELCTOR_ARGS){  }
    CProxyElement_MeshLibrary(CkGroupID _gid,int _onPE) : CProxyElement_NodeGroup(_gid,_onPE){  }
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
    MeshLibrary* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MeshLibrary* ckLocalBranch(CkGroupID gID) {
      return (MeshLibrary*)CkLocalNodeBranch(gID);
    }
/* DECLS: MeshLibrary(const std::vector<std::string > &filenames);
 */

};
PUPmarshall(CProxyElement_MeshLibrary)
/* ---------------- collective proxy -------------- */
class CProxy_MeshLibrary: public CProxy_NodeGroup{
  public:
    typedef MeshLibrary local_t;
    typedef CkIndex_MeshLibrary index_t;
    typedef CProxy_MeshLibrary proxy_t;
    typedef CProxyElement_MeshLibrary element_t;
    typedef CProxySection_MeshLibrary section_t;

    CProxy_MeshLibrary(void) {}
    CProxy_MeshLibrary(const IrrGroup *g) : CProxy_NodeGroup(g){  }
    CProxy_MeshLibrary(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_NodeGroup(_gid,CK_DELCTOR_ARGS){  }
    CProxy_MeshLibrary(CkGroupID _gid) : CProxy_NodeGroup(_gid){  }
    CProxyElement_MeshLibrary operator[](int onPE) const
      {return CProxyElement_MeshLibrary(ckGetGroupID(),onPE,CK_DELCTOR_CALL);}
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
    MeshLibrary* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MeshLibrary* ckLocalBranch(CkGroupID gID) {
      return (MeshLibrary*)CkLocalNodeBranch(gID);
    }
/* DECLS: MeshLibrary(const std::vector<std::string > &filenames);
 */
    static CkGroupID ckNew(const std::vector<std::string > &filenames, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_MeshLibrary(const std::vector<std::string > &filenames, const CkEntryOptions *impl_e_opts=NULL);

};
PUPmarshall(CProxy_MeshLibrary)
/* ---------------- section proxy -------------- */
class CProxySection_MeshLibrary: public CProxySection_NodeGroup{
  public:
    typedef MeshLibrary local_t;
    typedef CkIndex_MeshLibrary index_t;
    typedef CProxy_MeshLibrary proxy_t;
    typedef CProxyElement_MeshLibrary element_t;
    typedef CProxySection_MeshLibrary section_t;

    CProxySection_MeshLibrary(void) {}
    CProxySection_MeshLibrary(const IrrGroup *g) : CProxySection_NodeGroup(g){  }
    CProxySection_MeshLibrary(const CkGroupID &_gid,const int *_pelist,int _npes,CK_DELCTOR_PARAM) : CProxySection_NodeGroup(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_MeshLibrary(const CkGroupID &_gid,const int *_pelist,int _npes) : CProxySection_NodeGroup(_gid,_pelist,_npes){  }
    CProxySection_MeshLibrary(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes) : CProxySection_NodeGroup(n,_gid,_pelist,_npes){  }
    CProxySection_MeshLibrary(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes,CK_DELCTOR_PARAM) : CProxySection_NodeGroup(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
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
    MeshLibrary* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MeshLibrary* ckLocalBranch(CkGroupID gID) {
      return (MeshLibrary*)CkLocalNodeBranch(gID);
    }
/* DECLS: MeshLibrary(const std::vector<std::string > &filenames);
 */

};
PUPmarshall(CProxySection_MeshLibrary)
typedef CBaseT1<NodeGroup, CProxy_MeshLibrary> CBase_MeshLibrary;

extern void _registermesh_library(void);
#endif
