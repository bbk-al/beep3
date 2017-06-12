#ifndef _DECL_parallel_fmm_octree_H_
#define _DECL_parallel_fmm_octree_H_
#include "charm++.h"
/* DECLS: template < typename CType, typename CProxy_FMMWorkerT > nodegroup ParallelFMMOctree: NodeGroup{
ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
void insert(const std::vector<CType > &items);
void insert(const CType &single);
void finalize(const CkCallback &cb);
void clear_waves(void);
void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
};
 */
template < typename CType, typename CProxy_FMMWorkerT >  class ParallelFMMOctree;
template < typename CType, typename CProxy_FMMWorkerT >  class CkIndex_ParallelFMMOctree;
template < typename CType, typename CProxy_FMMWorkerT >  class CProxy_ParallelFMMOctree;
template < typename CType, typename CProxy_FMMWorkerT >  class CProxyElement_ParallelFMMOctree;
template < typename CType, typename CProxy_FMMWorkerT >  class CProxySection_ParallelFMMOctree;
/* --------------- index object ------------------ */
template < typename CType, typename CProxy_FMMWorkerT > class CkIndex_ParallelFMMOctree:public CProxyElement_NodeGroup{
  public:
    typedef ParallelFMMOctree < CType, CProxy_FMMWorkerT >  local_t;
    typedef CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  index_t;
    typedef CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  proxy_t;
    typedef CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  element_t;
    typedef CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
 */
    static int __idx_ParallelFMMOctree_marshall1;
    static int ckNew(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1) { return __idx_ParallelFMMOctree_marshall1; }
    static void _call_ParallelFMMOctree_marshall1(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);
    static int _callmarshall_ParallelFMMOctree_marshall1(char* impl_buf,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);
    static void _marshallmessagepup_ParallelFMMOctree_marshall1(PUP::er &p,void *msg);

/* DECLS: void insert(const std::vector<CType > &items);
 */
    static int __idx_insert_marshall2;
    static int insert(const std::vector<CType > &items) { return __idx_insert_marshall2; }
    static void _call_insert_marshall2(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);
    static void _marshallmessagepup_insert_marshall2(PUP::er &p,void *msg);

/* DECLS: void insert(const CType &single);
 */
    static int __idx_insert_marshall3;
    static int insert(const CType &single) { return __idx_insert_marshall3; }
    static void _call_insert_marshall3(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);
    static void _marshallmessagepup_insert_marshall3(PUP::er &p,void *msg);

/* DECLS: void finalize(const CkCallback &cb);
 */
    static int __idx_finalize_marshall4;
    static int finalize(const CkCallback &cb) { return __idx_finalize_marshall4; }
    static void _call_finalize_marshall4(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);
    static int _callmarshall_finalize_marshall4(char* impl_buf,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);
    static void _marshallmessagepup_finalize_marshall4(PUP::er &p,void *msg);

/* DECLS: void clear_waves(void);
 */
    static int __idx_clear_waves_void;
    static int clear_waves(void) { return __idx_clear_waves_void; }
    static void _call_clear_waves_void(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);

/* DECLS: void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
 */
    static int __idx_request_data_marshall6;
    static int request_data(const CProxy_FMMWorkerT &FMMWorkerProxy) { return __idx_request_data_marshall6; }
    static void _call_request_data_marshall6(void* impl_msg,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);
    static int _callmarshall_request_data_marshall6(char* impl_buf,ParallelFMMOctree < CType, CProxy_FMMWorkerT > * impl_obj);
    static void _marshallmessagepup_request_data_marshall6(PUP::er &p,void *msg);

};
/* --------------- element proxy ------------------ */
template < typename CType, typename CProxy_FMMWorkerT > class CProxyElement_ParallelFMMOctree: public CProxyElement_NodeGroup{
  public:
    typedef ParallelFMMOctree < CType, CProxy_FMMWorkerT >  local_t;
    typedef CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  index_t;
    typedef CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  proxy_t;
    typedef CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  element_t;
    typedef CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  section_t;

    CProxyElement_ParallelFMMOctree(void) {}
    CProxyElement_ParallelFMMOctree(const IrrGroup *g) : CProxyElement_NodeGroup(g){  }
    CProxyElement_ParallelFMMOctree(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_NodeGroup(_gid,_onPE,CK_DELCTOR_ARGS){  }
    CProxyElement_ParallelFMMOctree(CkGroupID _gid,int _onPE) : CProxyElement_NodeGroup(_gid,_onPE){  }
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
    ParallelFMMOctree < CType, CProxy_FMMWorkerT > * ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static ParallelFMMOctree < CType, CProxy_FMMWorkerT > * ckLocalBranch(CkGroupID gID) {
      return (ParallelFMMOctree < CType, CProxy_FMMWorkerT > *)CkLocalNodeBranch(gID);
    }
/* DECLS: ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
 */

/* DECLS: void insert(const std::vector<CType > &items);
 */
    void insert(const std::vector<CType > &items, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void insert(const CType &single);
 */
    void insert(const CType &single, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void finalize(const CkCallback &cb);
 */
    void finalize(const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void clear_waves(void);
 */
    void clear_waves(void);

/* DECLS: void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
 */
    void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy, const CkEntryOptions *impl_e_opts=NULL);

};
/* ---------------- collective proxy -------------- */
template < typename CType, typename CProxy_FMMWorkerT > class CProxy_ParallelFMMOctree: public CProxy_NodeGroup{
  public:
    typedef ParallelFMMOctree < CType, CProxy_FMMWorkerT >  local_t;
    typedef CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  index_t;
    typedef CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  proxy_t;
    typedef CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  element_t;
    typedef CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  section_t;

    CProxy_ParallelFMMOctree(void) {}
    CProxy_ParallelFMMOctree(const IrrGroup *g) : CProxy_NodeGroup(g){  }
    CProxy_ParallelFMMOctree(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_NodeGroup(_gid,CK_DELCTOR_ARGS){  }
    CProxy_ParallelFMMOctree(CkGroupID _gid) : CProxy_NodeGroup(_gid){  }
    CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  operator[](int onPE) const
      {return CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT > (ckGetGroupID(),onPE,CK_DELCTOR_CALL);}
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
    ParallelFMMOctree < CType, CProxy_FMMWorkerT > * ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static ParallelFMMOctree < CType, CProxy_FMMWorkerT > * ckLocalBranch(CkGroupID gID) {
      return (ParallelFMMOctree < CType, CProxy_FMMWorkerT > *)CkLocalNodeBranch(gID);
    }
/* DECLS: ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
 */
    static CkGroupID ckNew(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void insert(const std::vector<CType > &items);
 */
    void insert(const std::vector<CType > &items, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void insert(const CType &single);
 */
    void insert(const CType &single, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void finalize(const CkCallback &cb);
 */
    void finalize(const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void clear_waves(void);
 */
    void clear_waves(void);

/* DECLS: void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
 */
    void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy, const CkEntryOptions *impl_e_opts=NULL);

};
/* ---------------- section proxy -------------- */
template < typename CType, typename CProxy_FMMWorkerT > class CProxySection_ParallelFMMOctree: public CProxySection_NodeGroup{
  public:
    typedef ParallelFMMOctree < CType, CProxy_FMMWorkerT >  local_t;
    typedef CkIndex_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  index_t;
    typedef CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  proxy_t;
    typedef CProxyElement_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  element_t;
    typedef CProxySection_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  section_t;

    CProxySection_ParallelFMMOctree(void) {}
    CProxySection_ParallelFMMOctree(const IrrGroup *g) : CProxySection_NodeGroup(g){  }
    CProxySection_ParallelFMMOctree(const CkGroupID &_gid,const int *_pelist,int _npes,CK_DELCTOR_PARAM) : CProxySection_NodeGroup(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_ParallelFMMOctree(const CkGroupID &_gid,const int *_pelist,int _npes) : CProxySection_NodeGroup(_gid,_pelist,_npes){  }
    CProxySection_ParallelFMMOctree(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes) : CProxySection_NodeGroup(n,_gid,_pelist,_npes){  }
    CProxySection_ParallelFMMOctree(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes,CK_DELCTOR_PARAM) : CProxySection_NodeGroup(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
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
    ParallelFMMOctree < CType, CProxy_FMMWorkerT > * ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static ParallelFMMOctree < CType, CProxy_FMMWorkerT > * ckLocalBranch(CkGroupID gID) {
      return (ParallelFMMOctree < CType, CProxy_FMMWorkerT > *)CkLocalNodeBranch(gID);
    }
/* DECLS: ParallelFMMOctree(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, unsigned int impl_noname_0, const CkCallback &impl_noname_1);
 */

/* DECLS: void insert(const std::vector<CType > &items);
 */
    void insert(const std::vector<CType > &items, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void insert(const CType &single);
 */
    void insert(const CType &single, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void finalize(const CkCallback &cb);
 */
    void finalize(const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void clear_waves(void);
 */
    void clear_waves(void);

/* DECLS: void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy);
 */
    void request_data(const CProxy_FMMWorkerT &FMMWorkerProxy, const CkEntryOptions *impl_e_opts=NULL);

};
template < typename CType, typename CProxy_FMMWorkerT > 
class CBase_ParallelFMMOctree : public CBaseT1<NodeGroup, CProxy_ParallelFMMOctree < CType, CProxy_FMMWorkerT >  > { };

extern void _registerparallel_fmm_octree(void);
#endif
