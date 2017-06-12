#ifndef _DECL_fmm_worker_H_
#define _DECL_fmm_worker_H_
#include "charm++.h"
#include "fh_values_nodegroup.decl.h"

using beepp::MultipoleHolderT;

using beepp::PlaneWaveHolderT;

namespace beepp {
/* DECLS: template < class HType > message FMM_Msg;
 */
template < class HType > class FMM_Msg;
template < class HType > class CMessage_FMM_Msg:public CkMessage{
  public:
    static int __idx;
    void* operator new(size_t, void*p) { return p; }
    void* operator new(size_t);
    void* operator new(size_t, int*, const int);
    void* operator new(size_t, int*);
#if CMK_MULTIPLE_DELETE
    void operator delete(void*p, void*){dealloc(p);}
    void operator delete(void*p){dealloc(p);}
    void operator delete(void*p, int*, const int){dealloc(p);}
    void operator delete(void*p, int*){dealloc(p);}
#endif
    void operator delete(void*p, size_t){dealloc(p);}
    static void* alloc(int,size_t, int*, int);
    static void dealloc(void *p);
    CMessage_FMM_Msg < HType > ();
    static void *pack(FMM_Msg < HType >  *p);
    static FMM_Msg < HType > * unpack(void* p);
    void *operator new(size_t, const int);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, const int){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

/* DECLS: template < int NTERMS, int NLAMBS, int NWAVES > array FMMWorkerT: ArrayElement{
FMMWorkerT(CkMigrateMessage* impl_msg);
FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
void form_multipoles(void);
void make_plane_waves(int dummy);
void pass_multipole_upwards(void);
void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
void evaluate(void);
void debug_chk(void);
void check_waves_complete(int impl_noname_2);
void request_collected_waves(void);
void receive_incoming_wave(CkMessage* impl_msg);
};
 */
template < int NTERMS, int NLAMBS, int NWAVES >  class FMMWorkerT;
template < int NTERMS, int NLAMBS, int NWAVES >  class CkIndex_FMMWorkerT;
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxy_FMMWorkerT;
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxyElement_FMMWorkerT;
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxySection_FMMWorkerT;
/* --------------- index object ------------------ */
template < int NTERMS, int NLAMBS, int NWAVES > class CkIndex_FMMWorkerT:public CProxyElement_ArrayElement{
  public:
    typedef FMMWorkerT < NTERMS, NLAMBS, NWAVES >  local_t;
    typedef CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  index_t;
    typedef CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  proxy_t;
    typedef CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  element_t;
    typedef CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: FMMWorkerT(CkMigrateMessage* impl_msg);
 */
    static int __idx_FMMWorkerT_CkMigrateMessage;
    static int ckNew(CkMigrateMessage* impl_msg) { return __idx_FMMWorkerT_CkMigrateMessage; }
    static void _call_FMMWorkerT_CkMigrateMessage(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
 */
    static int __idx_FMMWorkerT_marshall1;
    static int ckNew(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length) { return __idx_FMMWorkerT_marshall1; }
    static void _call_FMMWorkerT_marshall1(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static int _callmarshall_FMMWorkerT_marshall1(char* impl_buf,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _marshallmessagepup_FMMWorkerT_marshall1(PUP::er &p,void *msg);

/* DECLS: void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
 */
    static int __idx_solve_marshall2;
    static int solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb) { return __idx_solve_marshall2; }
    static void _call_solve_marshall2(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static int _callmarshall_solve_marshall2(char* impl_buf,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _marshallmessagepup_solve_marshall2(PUP::er &p,void *msg);

/* DECLS: void form_multipoles(void);
 */
    static int __idx_form_multipoles_void;
    static int form_multipoles(void) { return __idx_form_multipoles_void; }
    static void _call_form_multipoles_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void make_plane_waves(int dummy);
 */
    static int __idx_make_plane_waves_marshall4;
    static int make_plane_waves(int dummy) { return __idx_make_plane_waves_marshall4; }
    static void _call_make_plane_waves_marshall4(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static int _callmarshall_make_plane_waves_marshall4(char* impl_buf,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _marshallmessagepup_make_plane_waves_marshall4(PUP::er &p,void *msg);

/* DECLS: void pass_multipole_upwards(void);
 */
    static int __idx_pass_multipole_upwards_void;
    static int pass_multipole_upwards(void) { return __idx_pass_multipole_upwards_void; }
    static void _call_pass_multipole_upwards_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    static int __idx_inherit_lpole_from_parent_FMM_Msg;
    static int inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) { return __idx_inherit_lpole_from_parent_FMM_Msg; }
    static void _call_inherit_lpole_from_parent_FMM_Msg(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    static int __idx_receive_multipole_contribution_from_child_FMM_Msg;
    static int receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) { return __idx_receive_multipole_contribution_from_child_FMM_Msg; }
    static void _call_receive_multipole_contribution_from_child_FMM_Msg(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void evaluate(void);
 */
    static int __idx_evaluate_void;
    static int evaluate(void) { return __idx_evaluate_void; }
    static void _call_evaluate_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void debug_chk(void);
 */
    static int __idx_debug_chk_void;
    static int debug_chk(void) { return __idx_debug_chk_void; }
    static void _call_debug_chk_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void check_waves_complete(int impl_noname_2);
 */
    static int __idx_check_waves_complete_marshall10;
    static int check_waves_complete(int impl_noname_2) { return __idx_check_waves_complete_marshall10; }
    static void _call_check_waves_complete_marshall10(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static int _callmarshall_check_waves_complete_marshall10(char* impl_buf,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _marshallmessagepup_check_waves_complete_marshall10(PUP::er &p,void *msg);

/* DECLS: void request_collected_waves(void);
 */
    static int __idx_request_collected_waves_void;
    static int request_collected_waves(void) { return __idx_request_collected_waves_void; }
    static void _call_request_collected_waves_void(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void receive_incoming_wave(CkMessage* impl_msg);
 */
    static int __idx_receive_incoming_wave_CkMessage;
    static int receive_incoming_wave(CkMessage* impl_msg) { return __idx_receive_incoming_wave_CkMessage; }
    static void _call_receive_incoming_wave_CkMessage(void* impl_msg,FMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

};
/* --------------- element proxy ------------------ */
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxyElement_FMMWorkerT : public CProxyElement_ArrayElement{
  public:
    typedef FMMWorkerT < NTERMS, NLAMBS, NWAVES >  local_t;
    typedef CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  index_t;
    typedef CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  proxy_t;
    typedef CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  element_t;
    typedef CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  section_t;

    CProxyElement_FMMWorkerT(void) {}
    CProxyElement_FMMWorkerT(const ArrayElement *e) : CProxyElement_ArrayElement(e){  }
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxyElement_ArrayElement::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxyElement_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxyElement_ArrayElement::pup(p);
    }
int ckIsDelegated(void) const {return CProxyElement_ArrayElement::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxyElement_ArrayElement::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxyElement_ArrayElement::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxyElement_ArrayElement::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxyElement_ArrayElement::ckCheck();}
inline operator CkArrayID () const {return ckGetArrayID();}
inline static CkArrayID ckCreateEmptyArray(void){ return CProxyElement_ArrayElement::ckCreateEmptyArray(); }
inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts){ return CProxyElement_ArrayElement::ckCreateArray(m,ctor,opts); }
inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx){ CProxyElement_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const{ CProxyElement_ArrayElement::ckBroadcast(m,ep,opts); }
inline CkArrayID ckGetArrayID(void) const{ return CProxyElement_ArrayElement::ckGetArrayID();}
inline CkArray *ckLocalBranch(void) const{ return CProxyElement_ArrayElement::ckLocalBranch(); }
inline CkLocMgr *ckLocMgr(void) const{ return CProxyElement_ArrayElement::ckLocMgr(); }
inline void doneInserting(void) { CProxyElement_ArrayElement::doneInserting(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxyElement_ArrayElement::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxyElement_ArrayElement::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxyElement_ArrayElement::ckSetReductionClient(cb); }
inline void ckInsert(CkArrayMessage *m,int ctor,int onPe)
  { CProxyElement_ArrayElement::ckInsert(m,ctor,onPe); }
inline void ckSend(CkArrayMessage *m, int ep, int opts = 0) const
  { CProxyElement_ArrayElement::ckSend(m,ep,opts); }
inline void *ckSendSync(CkArrayMessage *m, int ep) const
  { return CProxyElement_ArrayElement::ckSendSync(m,ep); }
inline const CkArrayIndex &ckGetIndex() const
  { return CProxyElement_ArrayElement::ckGetIndex(); }
    FMMWorkerT < NTERMS, NLAMBS, NWAVES >  *ckLocal(void) const
      { return (FMMWorkerT < NTERMS, NLAMBS, NWAVES >  *)CProxyElement_ArrayElement::ckLocal(); }
    CProxyElement_FMMWorkerT(const CkArrayID &aid,const CkArrayIndexOctreeIndexer &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS) {}
    CProxyElement_FMMWorkerT(const CkArrayID &aid,const CkArrayIndexOctreeIndexer &idx)
        :CProxyElement_ArrayElement(aid,idx) {}
/* DECLS: FMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DECLS: FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
 */
    void insert(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length, int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
 */
    void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void form_multipoles(void);
 */
    void form_multipoles(void) ;

/* DECLS: void make_plane_waves(int dummy);
 */
    void make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void pass_multipole_upwards(void);
 */
    void pass_multipole_upwards(void) ;

/* DECLS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void evaluate(void);
 */
    void evaluate(void) ;

/* DECLS: void debug_chk(void);
 */
    void debug_chk(void) ;

/* DECLS: void check_waves_complete(int impl_noname_2);
 */
    void check_waves_complete(int impl_noname_2, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void request_collected_waves(void);
 */
    void request_collected_waves(void) ;

/* DECLS: void receive_incoming_wave(CkMessage* impl_msg);
 */
    void receive_incoming_wave(CkMessage* impl_msg) ;

};
/* ---------------- collective proxy -------------- */
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxy_FMMWorkerT : public CProxy_ArrayElement{
  public:
    typedef FMMWorkerT < NTERMS, NLAMBS, NWAVES >  local_t;
    typedef CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  index_t;
    typedef CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  proxy_t;
    typedef CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  element_t;
    typedef CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  section_t;

    CProxy_FMMWorkerT(void) {}
    CProxy_FMMWorkerT(const ArrayElement *e) : CProxy_ArrayElement(e){  }
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxy_ArrayElement::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxy_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_ArrayElement::pup(p);
    }
int ckIsDelegated(void) const {return CProxy_ArrayElement::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxy_ArrayElement::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxy_ArrayElement::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxy_ArrayElement::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxy_ArrayElement::ckCheck();}
inline operator CkArrayID () const {return ckGetArrayID();}
inline static CkArrayID ckCreateEmptyArray(void){ return CProxy_ArrayElement::ckCreateEmptyArray(); }
inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts){ return CProxy_ArrayElement::ckCreateArray(m,ctor,opts); }
inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx){ CProxy_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const{ CProxy_ArrayElement::ckBroadcast(m,ep,opts); }
inline CkArrayID ckGetArrayID(void) const{ return CProxy_ArrayElement::ckGetArrayID();}
inline CkArray *ckLocalBranch(void) const{ return CProxy_ArrayElement::ckLocalBranch(); }
inline CkLocMgr *ckLocMgr(void) const{ return CProxy_ArrayElement::ckLocMgr(); }
inline void doneInserting(void) { CProxy_ArrayElement::doneInserting(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxy_ArrayElement::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxy_ArrayElement::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxy_ArrayElement::ckSetReductionClient(cb); }
    static CkArrayID ckNew(void) {return ckCreateEmptyArray();}
//Generalized array indexing:
    CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  operator [] (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > (ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  operator() (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > (ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxy_FMMWorkerT(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_FMMWorkerT(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: FMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DECLS: FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
 */
    static CkArrayID ckNew(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
 */
    void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void form_multipoles(void);
 */
    void form_multipoles(void) ;

/* DECLS: void make_plane_waves(int dummy);
 */
    void make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void pass_multipole_upwards(void);
 */
    void pass_multipole_upwards(void) ;

/* DECLS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void evaluate(void);
 */
    void evaluate(void) ;

/* DECLS: void debug_chk(void);
 */
    void debug_chk(void) ;

/* DECLS: void check_waves_complete(int impl_noname_2);
 */
    void check_waves_complete(int impl_noname_2, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void request_collected_waves(void);
 */
    void request_collected_waves(void) ;

/* DECLS: void receive_incoming_wave(CkMessage* impl_msg);
 */
    void receive_incoming_wave(CkMessage* impl_msg) ;

};
/* ---------------- section proxy -------------- */
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxySection_FMMWorkerT : public CProxySection_ArrayElement{
  public:
    typedef FMMWorkerT < NTERMS, NLAMBS, NWAVES >  local_t;
    typedef CkIndex_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  index_t;
    typedef CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  proxy_t;
    typedef CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  element_t;
    typedef CProxySection_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  section_t;

    CProxySection_FMMWorkerT(void) {}
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxySection_ArrayElement::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxySection_ArrayElement::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxySection_ArrayElement::pup(p);
    }
int ckIsDelegated(void) const {return CProxySection_ArrayElement::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxySection_ArrayElement::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxySection_ArrayElement::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxySection_ArrayElement::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxySection_ArrayElement::ckCheck();}
inline operator CkArrayID () const {return ckGetArrayID();}
inline static CkArrayID ckCreateEmptyArray(void){ return CProxySection_ArrayElement::ckCreateEmptyArray(); }
inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts){ return CProxySection_ArrayElement::ckCreateArray(m,ctor,opts); }
inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx){ CProxySection_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const{ CProxySection_ArrayElement::ckBroadcast(m,ep,opts); }
inline CkArrayID ckGetArrayID(void) const{ return CProxySection_ArrayElement::ckGetArrayID();}
inline CkArray *ckLocalBranch(void) const{ return CProxySection_ArrayElement::ckLocalBranch(); }
inline CkLocMgr *ckLocMgr(void) const{ return CProxySection_ArrayElement::ckLocMgr(); }
inline void doneInserting(void) { CProxySection_ArrayElement::doneInserting(); }
inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxySection_ArrayElement::setReductionClient(fn,param); }
inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
{ CProxySection_ArrayElement::ckSetReductionClient(fn,param); }
inline void ckSetReductionClient(CkCallback *cb) const
{ CProxySection_ArrayElement::ckSetReductionClient(cb); }
inline void ckSend(CkArrayMessage *m, int ep, int opts = 0)
 { CProxySection_ArrayElement::ckSend(m,ep,opts); }
inline CkSectionInfo &ckGetSectionInfo()
  { return CProxySection_ArrayElement::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
  { return CProxySection_ArrayElement::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
  { return CProxySection_ArrayElement::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
  { return CProxySection_ArrayElement::ckGetSectionID(i); }
inline CkArrayID ckGetArrayIDn(int i) const
{return CProxySection_ArrayElement::ckGetArrayIDn(i); } 
inline CkArrayIndex *ckGetArrayElements() const
  { return CProxySection_ArrayElement::ckGetArrayElements(); }
inline CkArrayIndex *ckGetArrayElements(int i) const
{return CProxySection_ArrayElement::ckGetArrayElements(i); }
inline int ckGetNumElements() const
  { return CProxySection_ArrayElement::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{return CProxySection_ArrayElement::ckGetNumElements(i); } 
//Generalized array indexing:
    CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  operator [] (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > (ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  operator() (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_FMMWorkerT < NTERMS, NLAMBS, NWAVES > (ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxySection_FMMWorkerT(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_FMMWorkerT(const CkArrayID &aid, CkArrayIndex *elems, int nElems) 
        :CProxySection_ArrayElement(aid,elems,nElems) {}
    CProxySection_FMMWorkerT(const CkSectionID &sid)       :CProxySection_ArrayElement(sid) {}
    CProxySection_FMMWorkerT(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_FMMWorkerT(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) {}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
/* DECLS: FMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DECLS: FMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelFMMOctree<CharmNodePatch,CProxy_FMMWorkerT<NTERMS,NLAMBS,NWAVES > > &fmm_tree_proxy, const CProxy_FH_Values_NodeGroup &fh_proxy, double edge_length);
 */

/* DECLS: void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb);
 */
    void solve(double impl_noname_0, double impl_noname_1, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void form_multipoles(void);
 */
    void form_multipoles(void) ;

/* DECLS: void make_plane_waves(int dummy);
 */
    void make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void pass_multipole_upwards(void);
 */
    void pass_multipole_upwards(void) ;

/* DECLS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void evaluate(void);
 */
    void evaluate(void) ;

/* DECLS: void debug_chk(void);
 */
    void debug_chk(void) ;

/* DECLS: void check_waves_complete(int impl_noname_2);
 */
    void check_waves_complete(int impl_noname_2, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void request_collected_waves(void);
 */
    void request_collected_waves(void) ;

/* DECLS: void receive_incoming_wave(CkMessage* impl_msg);
 */
    void receive_incoming_wave(CkMessage* impl_msg) ;

};
template < int NTERMS, int NLAMBS, int NWAVES > 
class CBase_FMMWorkerT : public CBaseT1<ArrayElementT<OctreeIndexer>, CProxy_FMMWorkerT < NTERMS, NLAMBS, NWAVES >  > { };




} // namespace beepp

extern void _registerfmm_worker(void);
#endif
