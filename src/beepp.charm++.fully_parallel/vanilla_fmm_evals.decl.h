#ifndef _DECL_vanilla_fmm_evals_H_
#define _DECL_vanilla_fmm_evals_H_
#include "charm++.h"
namespace vanilla_fmm {
/* DECLS: message Eval_Message{
EvalPt data[];
}
;
 */
class Eval_Message;
class CMessage_Eval_Message:public CkMessage{
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
    CMessage_Eval_Message();
    static void *pack(Eval_Message *p);
    static Eval_Message* unpack(void* p);
    void *operator new(size_t, int);
    void *operator new(size_t, int, const int);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, int){dealloc(p);}
    void operator delete(void *p, int, const int){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

/* DECLS: chare Vanilla_FMM_Evals: Chare{
Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy);
void receive_eval_results(Eval_Message* impl_msg);
void evaluate(const CkCallback &cb);
};
 */
 class Vanilla_FMM_Evals;
 class CkIndex_Vanilla_FMM_Evals;
 class CProxy_Vanilla_FMM_Evals;
/* --------------- index object ------------------ */
class CkIndex_Vanilla_FMM_Evals:public CProxy_Chare{
  public:
    typedef Vanilla_FMM_Evals local_t;
    typedef CkIndex_Vanilla_FMM_Evals index_t;
    typedef CProxy_Vanilla_FMM_Evals proxy_t;
    typedef CProxy_Vanilla_FMM_Evals element_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy);
 */
    static int __idx_Vanilla_FMM_Evals_marshall1;
    static int ckNew(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy) { return __idx_Vanilla_FMM_Evals_marshall1; }
    static void _call_Vanilla_FMM_Evals_marshall1(void* impl_msg,Vanilla_FMM_Evals* impl_obj);
    static int _callmarshall_Vanilla_FMM_Evals_marshall1(char* impl_buf,Vanilla_FMM_Evals* impl_obj);
    static void _marshallmessagepup_Vanilla_FMM_Evals_marshall1(PUP::er &p,void *msg);

/* DECLS: void receive_eval_results(Eval_Message* impl_msg);
 */
    static int __idx_receive_eval_results_Eval_Message;
    static int receive_eval_results(Eval_Message* impl_msg) { return __idx_receive_eval_results_Eval_Message; }
    static void _call_receive_eval_results_Eval_Message(void* impl_msg,Vanilla_FMM_Evals* impl_obj);

/* DECLS: void evaluate(const CkCallback &cb);
 */
    static int __idx_evaluate_marshall3;
    static int evaluate(const CkCallback &cb) { return __idx_evaluate_marshall3; }
    static void _call_evaluate_marshall3(void* impl_msg,Vanilla_FMM_Evals* impl_obj);
    static int _callmarshall_evaluate_marshall3(char* impl_buf,Vanilla_FMM_Evals* impl_obj);
    static void _marshallmessagepup_evaluate_marshall3(PUP::er &p,void *msg);

};
/* --------------- element proxy ------------------ */
class CProxy_Vanilla_FMM_Evals:public CProxy_Chare{
  public:
    typedef Vanilla_FMM_Evals local_t;
    typedef CkIndex_Vanilla_FMM_Evals index_t;
    typedef CProxy_Vanilla_FMM_Evals proxy_t;
    typedef CProxy_Vanilla_FMM_Evals element_t;

    CProxy_Vanilla_FMM_Evals(void) {};
    CProxy_Vanilla_FMM_Evals(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_Vanilla_FMM_Evals(const Chare *c) : CProxy_Chare(c){  }
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
    Vanilla_FMM_Evals *ckLocal(void) const
     { return (Vanilla_FMM_Evals *)CkLocalChare(&ckGetChareID()); }
/* DECLS: Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy);
 */
    static CkChareID ckNew(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);
    static void ckNew(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, CkChareID* pcid, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_Vanilla_FMM_Evals(const std::vector<Vector > &pts, const CProxy_ParallelFMMOctree<Charge,CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void receive_eval_results(Eval_Message* impl_msg);
 */
    void receive_eval_results(Eval_Message* impl_msg);

/* DECLS: void evaluate(const CkCallback &cb);
 */
    void evaluate(const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

};
PUPmarshall(CProxy_Vanilla_FMM_Evals)
typedef CBaseT1<Chare, CProxy_Vanilla_FMM_Evals> CBase_Vanilla_FMM_Evals;

} // namespace vanilla_fmm

extern void _registervanilla_fmm_evals(void);
#endif
