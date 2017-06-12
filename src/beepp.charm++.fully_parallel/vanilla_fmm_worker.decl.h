#ifndef _DECL_vanilla_fmm_worker_H_
#define _DECL_vanilla_fmm_worker_H_
#include "charm++.h"
#include "vanilla_fmm_evals.decl.h"

namespace vanilla_fmm {
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

/* DECLS: template < int NTERMS, int NLAMBS, int NWAVES > array VanillaFMMWorkerT: ArrayElement{
VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
threaded void evaluate(Eval_Message* impl_msg);
void check_waves_complete(int impl_noname_0);
void request_collected_waves(void);
void receive_incoming_wave(CkMessage* impl_msg);
void solve(double impl_noname_1, const CkCallback &impl_noname_2);
void make_plane_waves(int dummy);
};
 */
template < int NTERMS, int NLAMBS, int NWAVES >  class VanillaFMMWorkerT;
template < int NTERMS, int NLAMBS, int NWAVES >  class CkIndex_VanillaFMMWorkerT;
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxy_VanillaFMMWorkerT;
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxyElement_VanillaFMMWorkerT;
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxySection_VanillaFMMWorkerT;
/* --------------- index object ------------------ */
template < int NTERMS, int NLAMBS, int NWAVES > class CkIndex_VanillaFMMWorkerT:public CProxyElement_ArrayElement{
  public:
    typedef VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  local_t;
    typedef CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  index_t;
    typedef CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  proxy_t;
    typedef CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  element_t;
    typedef CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
 */
    static int __idx_VanillaFMMWorkerT_CkMigrateMessage;
    static int ckNew(CkMigrateMessage* impl_msg) { return __idx_VanillaFMMWorkerT_CkMigrateMessage; }
    static void _call_VanillaFMMWorkerT_CkMigrateMessage(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
 */
    static int __idx_VanillaFMMWorkerT_marshall1;
    static int ckNew(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length) { return __idx_VanillaFMMWorkerT_marshall1; }
    static void _call_VanillaFMMWorkerT_marshall1(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static int _callmarshall_VanillaFMMWorkerT_marshall1(char* impl_buf,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _marshallmessagepup_VanillaFMMWorkerT_marshall1(PUP::er &p,void *msg);

/* DECLS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    static int __idx_inherit_lpole_from_parent_FMM_Msg;
    static int inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) { return __idx_inherit_lpole_from_parent_FMM_Msg; }
    static void _call_inherit_lpole_from_parent_FMM_Msg(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    static int __idx_receive_multipole_contribution_from_child_FMM_Msg;
    static int receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) { return __idx_receive_multipole_contribution_from_child_FMM_Msg; }
    static void _call_receive_multipole_contribution_from_child_FMM_Msg(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: threaded void evaluate(Eval_Message* impl_msg);
 */
    static int __idx_evaluate_Eval_Message;
    static int evaluate(Eval_Message* impl_msg) { return __idx_evaluate_Eval_Message; }
    static void _call_evaluate_Eval_Message(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _callthr_evaluate_Eval_Message(CkThrCallArg *);

/* DECLS: void check_waves_complete(int impl_noname_0);
 */
    static int __idx_check_waves_complete_marshall5;
    static int check_waves_complete(int impl_noname_0) { return __idx_check_waves_complete_marshall5; }
    static void _call_check_waves_complete_marshall5(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static int _callmarshall_check_waves_complete_marshall5(char* impl_buf,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _marshallmessagepup_check_waves_complete_marshall5(PUP::er &p,void *msg);

/* DECLS: void request_collected_waves(void);
 */
    static int __idx_request_collected_waves_void;
    static int request_collected_waves(void) { return __idx_request_collected_waves_void; }
    static void _call_request_collected_waves_void(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void receive_incoming_wave(CkMessage* impl_msg);
 */
    static int __idx_receive_incoming_wave_CkMessage;
    static int receive_incoming_wave(CkMessage* impl_msg) { return __idx_receive_incoming_wave_CkMessage; }
    static void _call_receive_incoming_wave_CkMessage(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);

/* DECLS: void solve(double impl_noname_1, const CkCallback &impl_noname_2);
 */
    static int __idx_solve_marshall8;
    static int solve(double impl_noname_1, const CkCallback &impl_noname_2) { return __idx_solve_marshall8; }
    static void _call_solve_marshall8(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static int _callmarshall_solve_marshall8(char* impl_buf,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _marshallmessagepup_solve_marshall8(PUP::er &p,void *msg);

/* DECLS: void make_plane_waves(int dummy);
 */
    static int __idx_make_plane_waves_marshall9;
    static int make_plane_waves(int dummy) { return __idx_make_plane_waves_marshall9; }
    static void _call_make_plane_waves_marshall9(void* impl_msg,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static int _callmarshall_make_plane_waves_marshall9(char* impl_buf,VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > * impl_obj);
    static void _marshallmessagepup_make_plane_waves_marshall9(PUP::er &p,void *msg);

};
/* --------------- element proxy ------------------ */
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxyElement_VanillaFMMWorkerT : public CProxyElement_ArrayElement{
  public:
    typedef VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  local_t;
    typedef CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  index_t;
    typedef CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  proxy_t;
    typedef CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  element_t;
    typedef CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  section_t;

    CProxyElement_VanillaFMMWorkerT(void) {}
    CProxyElement_VanillaFMMWorkerT(const ArrayElement *e) : CProxyElement_ArrayElement(e){  }
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
    VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  *ckLocal(void) const
      { return (VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  *)CProxyElement_ArrayElement::ckLocal(); }
    CProxyElement_VanillaFMMWorkerT(const CkArrayID &aid,const CkArrayIndexOctreeIndexer &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS) {}
    CProxyElement_VanillaFMMWorkerT(const CkArrayID &aid,const CkArrayIndexOctreeIndexer &idx)
        :CProxyElement_ArrayElement(aid,idx) {}
/* DECLS: VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DECLS: VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
 */
    void insert(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length, int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: threaded void evaluate(Eval_Message* impl_msg);
 */
    void evaluate(Eval_Message* impl_msg) ;

/* DECLS: void check_waves_complete(int impl_noname_0);
 */
    void check_waves_complete(int impl_noname_0, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void request_collected_waves(void);
 */
    void request_collected_waves(void) ;

/* DECLS: void receive_incoming_wave(CkMessage* impl_msg);
 */
    void receive_incoming_wave(CkMessage* impl_msg) ;

/* DECLS: void solve(double impl_noname_1, const CkCallback &impl_noname_2);
 */
    void solve(double impl_noname_1, const CkCallback &impl_noname_2, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void make_plane_waves(int dummy);
 */
    void make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts=NULL) ;

};
/* ---------------- collective proxy -------------- */
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxy_VanillaFMMWorkerT : public CProxy_ArrayElement{
  public:
    typedef VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  local_t;
    typedef CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  index_t;
    typedef CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  proxy_t;
    typedef CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  element_t;
    typedef CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  section_t;

    CProxy_VanillaFMMWorkerT(void) {}
    CProxy_VanillaFMMWorkerT(const ArrayElement *e) : CProxy_ArrayElement(e){  }
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
    CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  operator [] (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > (ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  operator() (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > (ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxy_VanillaFMMWorkerT(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_VanillaFMMWorkerT(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DECLS: VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
 */
    static CkArrayID ckNew(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: threaded void evaluate(Eval_Message* impl_msg);
 */
    void evaluate(Eval_Message* impl_msg) ;

/* DECLS: void check_waves_complete(int impl_noname_0);
 */
    void check_waves_complete(int impl_noname_0, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void request_collected_waves(void);
 */
    void request_collected_waves(void) ;

/* DECLS: void receive_incoming_wave(CkMessage* impl_msg);
 */
    void receive_incoming_wave(CkMessage* impl_msg) ;

/* DECLS: void solve(double impl_noname_1, const CkCallback &impl_noname_2);
 */
    void solve(double impl_noname_1, const CkCallback &impl_noname_2, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void make_plane_waves(int dummy);
 */
    void make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts=NULL) ;

};
/* ---------------- section proxy -------------- */
template < int NTERMS, int NLAMBS, int NWAVES >  class CProxySection_VanillaFMMWorkerT : public CProxySection_ArrayElement{
  public:
    typedef VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  local_t;
    typedef CkIndex_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  index_t;
    typedef CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  proxy_t;
    typedef CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  element_t;
    typedef CProxySection_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  section_t;

    CProxySection_VanillaFMMWorkerT(void) {}
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
    CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  operator [] (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > (ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  operator() (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES > (ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxySection_VanillaFMMWorkerT(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_VanillaFMMWorkerT(const CkArrayID &aid, CkArrayIndex *elems, int nElems) 
        :CProxySection_ArrayElement(aid,elems,nElems) {}
    CProxySection_VanillaFMMWorkerT(const CkSectionID &sid)       :CProxySection_ArrayElement(sid) {}
    CProxySection_VanillaFMMWorkerT(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_VanillaFMMWorkerT(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) {}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
/* DECLS: VanillaFMMWorkerT(CkMigrateMessage* impl_msg);
 */

/* DECLS: VanillaFMMWorkerT(const CProxy_FMM_Globals_NodeGroupT<NTERMS,NLAMBS > &fmm_globals_proxy, const CProxy_ParallelTree &fmm_tree_proxy, double edge_length);
 */

/* DECLS: void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void inherit_lpole_from_parent(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg);
 */
    void receive_multipole_contribution_from_child(FMM_Msg<MultipoleHolderT<NTERMS > >* impl_msg) ;

/* DECLS: threaded void evaluate(Eval_Message* impl_msg);
 */
    void evaluate(Eval_Message* impl_msg) ;

/* DECLS: void check_waves_complete(int impl_noname_0);
 */
    void check_waves_complete(int impl_noname_0, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void request_collected_waves(void);
 */
    void request_collected_waves(void) ;

/* DECLS: void receive_incoming_wave(CkMessage* impl_msg);
 */
    void receive_incoming_wave(CkMessage* impl_msg) ;

/* DECLS: void solve(double impl_noname_1, const CkCallback &impl_noname_2);
 */
    void solve(double impl_noname_1, const CkCallback &impl_noname_2, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void make_plane_waves(int dummy);
 */
    void make_plane_waves(int dummy, const CkEntryOptions *impl_e_opts=NULL) ;

};
template < int NTERMS, int NLAMBS, int NWAVES > 
class CBase_VanillaFMMWorkerT : public CBaseT1<ArrayElementT<OctreeIndexer>, CProxy_VanillaFMMWorkerT < NTERMS, NLAMBS, NWAVES >  > { };

} // namespace vanilla_fmm





extern void _registervanilla_fmm_worker(void);
#endif
