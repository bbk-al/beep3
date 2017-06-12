#ifndef _DECL_fh_values_nodegroup_H_
#define _DECL_fh_values_nodegroup_H_
#include "charm++.h"
/* DECLS: group FH_Values_NodeGroup: IrrGroup{
FH_Values_NodeGroup(unsigned int impl_noname_0);
void set(const FH_Values &vals, const CkCallback &cb);
void reduce(const CkCallback &cb);
};
 */
 class FH_Values_NodeGroup;
 class CkIndex_FH_Values_NodeGroup;
 class CProxy_FH_Values_NodeGroup;
 class CProxyElement_FH_Values_NodeGroup;
 class CProxySection_FH_Values_NodeGroup;
/* --------------- index object ------------------ */
class CkIndex_FH_Values_NodeGroup:public CProxyElement_IrrGroup{
  public:
    typedef FH_Values_NodeGroup local_t;
    typedef CkIndex_FH_Values_NodeGroup index_t;
    typedef CProxy_FH_Values_NodeGroup proxy_t;
    typedef CProxyElement_FH_Values_NodeGroup element_t;
    typedef CProxySection_FH_Values_NodeGroup section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: FH_Values_NodeGroup(unsigned int impl_noname_0);
 */
    static int __idx_FH_Values_NodeGroup_marshall1;
    static int ckNew(unsigned int impl_noname_0) { return __idx_FH_Values_NodeGroup_marshall1; }
    static void _call_FH_Values_NodeGroup_marshall1(void* impl_msg,FH_Values_NodeGroup* impl_obj);
    static int _callmarshall_FH_Values_NodeGroup_marshall1(char* impl_buf,FH_Values_NodeGroup* impl_obj);
    static void _marshallmessagepup_FH_Values_NodeGroup_marshall1(PUP::er &p,void *msg);

/* DECLS: void set(const FH_Values &vals, const CkCallback &cb);
 */
    static int __idx_set_marshall2;
    static int set(const FH_Values &vals, const CkCallback &cb) { return __idx_set_marshall2; }
    static void _call_set_marshall2(void* impl_msg,FH_Values_NodeGroup* impl_obj);
    static int _callmarshall_set_marshall2(char* impl_buf,FH_Values_NodeGroup* impl_obj);
    static void _marshallmessagepup_set_marshall2(PUP::er &p,void *msg);

/* DECLS: void reduce(const CkCallback &cb);
 */
    static int __idx_reduce_marshall3;
    static int reduce(const CkCallback &cb) { return __idx_reduce_marshall3; }
    static void _call_reduce_marshall3(void* impl_msg,FH_Values_NodeGroup* impl_obj);
    static int _callmarshall_reduce_marshall3(char* impl_buf,FH_Values_NodeGroup* impl_obj);
    static void _marshallmessagepup_reduce_marshall3(PUP::er &p,void *msg);

};
/* --------------- element proxy ------------------ */
class CProxyElement_FH_Values_NodeGroup: public CProxyElement_IrrGroup{
  public:
    typedef FH_Values_NodeGroup local_t;
    typedef CkIndex_FH_Values_NodeGroup index_t;
    typedef CProxy_FH_Values_NodeGroup proxy_t;
    typedef CProxyElement_FH_Values_NodeGroup element_t;
    typedef CProxySection_FH_Values_NodeGroup section_t;

    CProxyElement_FH_Values_NodeGroup(void) {}
    CProxyElement_FH_Values_NodeGroup(const IrrGroup *g) : CProxyElement_IrrGroup(g){  }
    CProxyElement_FH_Values_NodeGroup(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_IrrGroup(_gid,_onPE,CK_DELCTOR_ARGS){  }
    CProxyElement_FH_Values_NodeGroup(CkGroupID _gid,int _onPE) : CProxyElement_IrrGroup(_gid,_onPE){  }
int ckIsDelegated(void) const {return CProxyElement_IrrGroup::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxyElement_IrrGroup::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxyElement_IrrGroup::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxyElement_IrrGroup::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxyElement_IrrGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxyElement_IrrGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxyElement_IrrGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxyElement_IrrGroup::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxyElement_IrrGroup::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxyElement_IrrGroup::ckSetReductionClient(cb); }
int ckGetGroupPe(void) const
{return CProxyElement_IrrGroup::ckGetGroupPe();}
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxyElement_IrrGroup::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxyElement_IrrGroup::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxyElement_IrrGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxyElement_IrrGroup::ckSetGroupID(g);
    }
    FH_Values_NodeGroup* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static FH_Values_NodeGroup* ckLocalBranch(CkGroupID gID) {
      return (FH_Values_NodeGroup*)CkLocalBranch(gID);
    }
/* DECLS: FH_Values_NodeGroup(unsigned int impl_noname_0);
 */

/* DECLS: void set(const FH_Values &vals, const CkCallback &cb);
 */
    void set(const FH_Values &vals, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void reduce(const CkCallback &cb);
 */
    void reduce(const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

};
PUPmarshall(CProxyElement_FH_Values_NodeGroup)
/* ---------------- collective proxy -------------- */
class CProxy_FH_Values_NodeGroup: public CProxy_IrrGroup{
  public:
    typedef FH_Values_NodeGroup local_t;
    typedef CkIndex_FH_Values_NodeGroup index_t;
    typedef CProxy_FH_Values_NodeGroup proxy_t;
    typedef CProxyElement_FH_Values_NodeGroup element_t;
    typedef CProxySection_FH_Values_NodeGroup section_t;

    CProxy_FH_Values_NodeGroup(void) {}
    CProxy_FH_Values_NodeGroup(const IrrGroup *g) : CProxy_IrrGroup(g){  }
    CProxy_FH_Values_NodeGroup(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_IrrGroup(_gid,CK_DELCTOR_ARGS){  }
    CProxy_FH_Values_NodeGroup(CkGroupID _gid) : CProxy_IrrGroup(_gid){  }
    CProxyElement_FH_Values_NodeGroup operator[](int onPE) const
      {return CProxyElement_FH_Values_NodeGroup(ckGetGroupID(),onPE,CK_DELCTOR_CALL);}
int ckIsDelegated(void) const {return CProxy_IrrGroup::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxy_IrrGroup::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxy_IrrGroup::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxy_IrrGroup::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxy_IrrGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxy_IrrGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxy_IrrGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxy_IrrGroup::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxy_IrrGroup::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxy_IrrGroup::ckSetReductionClient(cb); }
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxy_IrrGroup::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxy_IrrGroup::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_IrrGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxy_IrrGroup::ckSetGroupID(g);
    }
    FH_Values_NodeGroup* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static FH_Values_NodeGroup* ckLocalBranch(CkGroupID gID) {
      return (FH_Values_NodeGroup*)CkLocalBranch(gID);
    }
/* DECLS: FH_Values_NodeGroup(unsigned int impl_noname_0);
 */
    static CkGroupID ckNew(unsigned int impl_noname_0, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_FH_Values_NodeGroup(unsigned int impl_noname_0, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void set(const FH_Values &vals, const CkCallback &cb);
 */
    void set(const FH_Values &vals, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);
    void set(const FH_Values &vals, const CkCallback &cb, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    void set(const FH_Values &vals, const CkCallback &cb, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void reduce(const CkCallback &cb);
 */
    void reduce(const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);
    void reduce(const CkCallback &cb, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    void reduce(const CkCallback &cb, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

};
PUPmarshall(CProxy_FH_Values_NodeGroup)
/* ---------------- section proxy -------------- */
class CProxySection_FH_Values_NodeGroup: public CProxySection_IrrGroup{
  public:
    typedef FH_Values_NodeGroup local_t;
    typedef CkIndex_FH_Values_NodeGroup index_t;
    typedef CProxy_FH_Values_NodeGroup proxy_t;
    typedef CProxyElement_FH_Values_NodeGroup element_t;
    typedef CProxySection_FH_Values_NodeGroup section_t;

    CProxySection_FH_Values_NodeGroup(void) {}
    CProxySection_FH_Values_NodeGroup(const IrrGroup *g) : CProxySection_IrrGroup(g){  }
    CProxySection_FH_Values_NodeGroup(const CkGroupID &_gid,const int *_pelist,int _npes,CK_DELCTOR_PARAM) : CProxySection_IrrGroup(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_FH_Values_NodeGroup(const CkGroupID &_gid,const int *_pelist,int _npes) : CProxySection_IrrGroup(_gid,_pelist,_npes){  }
    CProxySection_FH_Values_NodeGroup(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes) : CProxySection_IrrGroup(n,_gid,_pelist,_npes){  }
    CProxySection_FH_Values_NodeGroup(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes,CK_DELCTOR_PARAM) : CProxySection_IrrGroup(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
int ckIsDelegated(void) const {return CProxySection_IrrGroup::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxySection_IrrGroup::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxySection_IrrGroup::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxySection_IrrGroup::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxySection_IrrGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxySection_IrrGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxySection_IrrGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxySection_IrrGroup::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxySection_IrrGroup::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxySection_IrrGroup::ckSetReductionClient(cb); }
inline int ckGetNumSections() const
{ return CProxySection_IrrGroup::ckGetNumSections(); }
inline CkSectionInfo &ckGetSectionInfo()
{ return CProxySection_IrrGroup::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
{ return CProxySection_IrrGroup::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
{ return CProxySection_IrrGroup::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
{ return CProxySection_IrrGroup::ckGetSectionID(i); }
inline CkGroupID ckGetGroupIDn(int i) const
{ return CProxySection_IrrGroup::ckGetGroupIDn(i); }
inline int *ckGetElements() const
{ return CProxySection_IrrGroup::ckGetElements(); }
inline int *ckGetElements(int i) const
{ return CProxySection_IrrGroup::ckGetElements(i); }
inline int ckGetNumElements() const
{ return CProxySection_IrrGroup::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{ return CProxySection_IrrGroup::ckGetNumElements(i); }
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxySection_IrrGroup::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxySection_IrrGroup::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxySection_IrrGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxySection_IrrGroup::ckSetGroupID(g);
    }
    FH_Values_NodeGroup* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static FH_Values_NodeGroup* ckLocalBranch(CkGroupID gID) {
      return (FH_Values_NodeGroup*)CkLocalBranch(gID);
    }
/* DECLS: FH_Values_NodeGroup(unsigned int impl_noname_0);
 */

/* DECLS: void set(const FH_Values &vals, const CkCallback &cb);
 */
    void set(const FH_Values &vals, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void reduce(const CkCallback &cb);
 */
    void reduce(const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

};
PUPmarshall(CProxySection_FH_Values_NodeGroup)
typedef CBaseT1<Group, CProxy_FH_Values_NodeGroup> CBase_FH_Values_NodeGroup;

extern void _registerfh_values_nodegroup(void);
#endif
