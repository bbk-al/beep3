#ifndef _DECL_iteration_handler_H_
#define _DECL_iteration_handler_H_
#include "charm++.h"
#include "fmm_worker.decl.h"

/* DECLS: chare IterationHandler: Chare{
IterationHandler(double _universe_edge_length, double _beta);
void do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals);
threaded void phase_two(void);
void done_evals(void);
void reduce_fh_results(void);
void done_reduction(CkReductionMsg* impl_msg);
sync void set_num_workers(unsigned int num_workers_in);
};
 */
 class IterationHandler;
 class CkIndex_IterationHandler;
 class CProxy_IterationHandler;
/* --------------- index object ------------------ */
class CkIndex_IterationHandler:public CProxy_Chare{
  public:
    typedef IterationHandler local_t;
    typedef CkIndex_IterationHandler index_t;
    typedef CProxy_IterationHandler proxy_t;
    typedef CProxy_IterationHandler element_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: IterationHandler(double _universe_edge_length, double _beta);
 */
    static int __idx_IterationHandler_marshall1;
    static int ckNew(double _universe_edge_length, double _beta) { return __idx_IterationHandler_marshall1; }
    static void _call_IterationHandler_marshall1(void* impl_msg,IterationHandler* impl_obj);
    static int _callmarshall_IterationHandler_marshall1(char* impl_buf,IterationHandler* impl_obj);
    static void _marshallmessagepup_IterationHandler_marshall1(PUP::er &p,void *msg);

/* DECLS: void do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals);
 */
    static int __idx_do_bem_fmm_iteration_marshall2;
    static int do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals) { return __idx_do_bem_fmm_iteration_marshall2; }
    static void _call_do_bem_fmm_iteration_marshall2(void* impl_msg,IterationHandler* impl_obj);
    static int _callmarshall_do_bem_fmm_iteration_marshall2(char* impl_buf,IterationHandler* impl_obj);
    static void _marshallmessagepup_do_bem_fmm_iteration_marshall2(PUP::er &p,void *msg);

/* DECLS: threaded void phase_two(void);
 */
    static int __idx_phase_two_void;
    static int phase_two(void) { return __idx_phase_two_void; }
    static void _call_phase_two_void(void* impl_msg,IterationHandler* impl_obj);
    static void _callthr_phase_two_void(CkThrCallArg *);

/* DECLS: void done_evals(void);
 */
    static int __idx_done_evals_void;
    static int done_evals(void) { return __idx_done_evals_void; }
    static void _call_done_evals_void(void* impl_msg,IterationHandler* impl_obj);

/* DECLS: void reduce_fh_results(void);
 */
    static int __idx_reduce_fh_results_void;
    static int reduce_fh_results(void) { return __idx_reduce_fh_results_void; }
    static void _call_reduce_fh_results_void(void* impl_msg,IterationHandler* impl_obj);

/* DECLS: void done_reduction(CkReductionMsg* impl_msg);
 */
    static int __idx_done_reduction_CkReductionMsg;
    static int done_reduction(CkReductionMsg* impl_msg) { return __idx_done_reduction_CkReductionMsg; }
    static void _call_done_reduction_CkReductionMsg(void* impl_msg,IterationHandler* impl_obj);

/* DECLS: sync void set_num_workers(unsigned int num_workers_in);
 */
    static int __idx_set_num_workers_marshall7;
    static int set_num_workers(unsigned int num_workers_in) { return __idx_set_num_workers_marshall7; }
    static void _call_set_num_workers_marshall7(void* impl_msg,IterationHandler* impl_obj);
    static void _marshallmessagepup_set_num_workers_marshall7(PUP::er &p,void *msg);

};
/* --------------- element proxy ------------------ */
class CProxy_IterationHandler:public CProxy_Chare{
  public:
    typedef IterationHandler local_t;
    typedef CkIndex_IterationHandler index_t;
    typedef CProxy_IterationHandler proxy_t;
    typedef CProxy_IterationHandler element_t;

    CProxy_IterationHandler(void) {};
    CProxy_IterationHandler(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_IterationHandler(const Chare *c) : CProxy_Chare(c){  }
int ckIsDelegated(void) const {return CProxy_Chare::ckIsDelegated();}
inline CkDelegateMgr *ckDelegatedTo(void) const {return CProxy_Chare::ckDelegatedTo();}
inline CkDelegateData *ckDelegatedPtr(void) const {return CProxy_Chare::ckDelegatedPtr();}
CkGroupID ckDelegatedIdx(void) const {return CProxy_Chare::ckDelegatedIdx();}
inline void ckCheck(void) const {CProxy_Chare::ckCheck();}
const CkChareID &ckGetChareID(void) const
{ return CProxy_Chare::ckGetChareID(); }
operator const CkChareID &(void) const {return ckGetChareID();}
    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL) {
      CProxy_Chare::ckDelegate(dTo,dPtr);
    }
    void ckUndelegate(void) {
      CProxy_Chare::ckUndelegate();
    }
    void pup(PUP::er &p) {
      CProxy_Chare::pup(p);
    }
    void ckSetChareID(const CkChareID &c) {
      CProxy_Chare::ckSetChareID(c);
    }
    IterationHandler *ckLocal(void) const
     { return (IterationHandler *)CkLocalChare(&ckGetChareID()); }
/* DECLS: IterationHandler(double _universe_edge_length, double _beta);
 */
    static CkChareID ckNew(double _universe_edge_length, double _beta, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);
    static void ckNew(double _universe_edge_length, double _beta, CkChareID* pcid, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_IterationHandler(double _universe_edge_length, double _beta, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals);
 */
    void do_bem_fmm_iteration(const CkCallback &cb, const FH_Values &fh_vals, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: threaded void phase_two(void);
 */
    void phase_two(void);

/* DECLS: void done_evals(void);
 */
    void done_evals(void);

/* DECLS: void reduce_fh_results(void);
 */
    void reduce_fh_results(void);

/* DECLS: void done_reduction(CkReductionMsg* impl_msg);
 */
    void done_reduction(CkReductionMsg* impl_msg);

/* DECLS: sync void set_num_workers(unsigned int num_workers_in);
 */
    void set_num_workers(unsigned int num_workers_in, const CkEntryOptions *impl_e_opts=NULL);

};
PUPmarshall(CProxy_IterationHandler)
typedef CBaseT1<Chare, CProxy_IterationHandler> CBase_IterationHandler;

extern void _registeriteration_handler(void);
#endif
