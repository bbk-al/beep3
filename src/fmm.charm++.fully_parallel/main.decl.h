#ifndef _DECL_main_H_
#define _DECL_main_H_
#include "charm++.h"
#include "fmm_globals_nodegroup.decl.h"

#include "parallel_fmm_octree.decl.h"

#include "opencl_nodegroup.decl.h"

#include "vanilla_fmm_worker.decl.h"

#include "vanilla_fmm_evals.decl.h"

#include "comlib.decl.h"


/* DECLS: readonly CProxy_Main MainProxy;
 */

/* DECLS: readonly CProxy_FMM_Globals_NodeGroup FMM_Globals_Proxy;
 */

/* DECLS: readonly CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > ParallelFMMOctreeProxy;
 */

/* DECLS: readonly vanilla_fmm::CProxy_VanillaFMMWorker VanillaFMMWorkerProxy;
 */

/* DECLS: readonly vanilla_fmm::CProxy_Vanilla_FMM_Evals Vanilla_FMM_Evals_Proxy;
 */

/* DECLS: readonly ComlibInstanceHandle streaming_strat;
 */

/* DECLS: readonly CProxy_OpenCL_NodeGroup OpenCL_NodeGroupProxy;
 */

/* DECLS: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
threaded void create_workers(void);
void fmm_worker_complete(void);
void completed(vanilla_fmm::Eval_Message* impl_msg);
void quiescenceHandler(void);
  initcall void initManualLB(void);
};
 */
 class Main;
 class CkIndex_Main;
 class CProxy_Main;
/* --------------- index object ------------------ */
class CkIndex_Main:public CProxy_Chare{
  public:
    typedef Main local_t;
    typedef CkIndex_Main index_t;
    typedef CProxy_Main proxy_t;
    typedef CProxy_Main element_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: Main(CkArgMsg* impl_msg);
 */
    static int __idx_Main_CkArgMsg;
    static int ckNew(CkArgMsg* impl_msg) { return __idx_Main_CkArgMsg; }
    static void _call_Main_CkArgMsg(void* impl_msg,Main* impl_obj);

/* DECLS: threaded void create_workers(void);
 */
    static int __idx_create_workers_void;
    static int create_workers(void) { return __idx_create_workers_void; }
    static void _call_create_workers_void(void* impl_msg,Main* impl_obj);
    static void _callthr_create_workers_void(CkThrCallArg *);

/* DECLS: void fmm_worker_complete(void);
 */
    static int __idx_fmm_worker_complete_void;
    static int fmm_worker_complete(void) { return __idx_fmm_worker_complete_void; }
    static void _call_fmm_worker_complete_void(void* impl_msg,Main* impl_obj);

/* DECLS: void completed(vanilla_fmm::Eval_Message* impl_msg);
 */
    static int __idx_completed_Eval_Message;
    static int completed(vanilla_fmm::Eval_Message* impl_msg) { return __idx_completed_Eval_Message; }
    static void _call_completed_Eval_Message(void* impl_msg,Main* impl_obj);

/* DECLS: void quiescenceHandler(void);
 */
    static int __idx_quiescenceHandler_void;
    static int quiescenceHandler(void) { return __idx_quiescenceHandler_void; }
    static void _call_quiescenceHandler_void(void* impl_msg,Main* impl_obj);


};
/* --------------- element proxy ------------------ */
class CProxy_Main:public CProxy_Chare{
  public:
    typedef Main local_t;
    typedef CkIndex_Main index_t;
    typedef CProxy_Main proxy_t;
    typedef CProxy_Main element_t;

    CProxy_Main(void) {};
    CProxy_Main(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_Main(const Chare *c) : CProxy_Chare(c){  }
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
    Main *ckLocal(void) const
     { return (Main *)CkLocalChare(&ckGetChareID()); }
/* DECLS: Main(CkArgMsg* impl_msg);
 */
    static CkChareID ckNew(CkArgMsg* impl_msg, int onPE=CK_PE_ANY);
    static void ckNew(CkArgMsg* impl_msg, CkChareID* pcid, int onPE=CK_PE_ANY);
    CProxy_Main(CkArgMsg* impl_msg, int onPE=CK_PE_ANY);

/* DECLS: threaded void create_workers(void);
 */
    void create_workers(void);

/* DECLS: void fmm_worker_complete(void);
 */
    void fmm_worker_complete(void);

/* DECLS: void completed(vanilla_fmm::Eval_Message* impl_msg);
 */
    void completed(vanilla_fmm::Eval_Message* impl_msg);

/* DECLS: void quiescenceHandler(void);
 */
    void quiescenceHandler(void);


};
PUPmarshall(CProxy_Main)
typedef CBaseT1<Chare, CProxy_Main> CBase_Main;

extern void _registermain(void);
extern "C" void CkRegisterMainModule(void);
#endif
