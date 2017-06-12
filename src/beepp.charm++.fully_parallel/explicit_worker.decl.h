#ifndef _DECL_explicit_worker_H_
#define _DECL_explicit_worker_H_
#include "charm++.h"
/* DECLS: array ExplicitWorker: ArrayElement{
ExplicitWorker(CkMigrateMessage* impl_msg);
ExplicitWorker(void);
void precalc(double kappa, unsigned int total_num_patches);
void evaluate(double kappa);
void triggerLB(void);
};
 */
 class ExplicitWorker;
 class CkIndex_ExplicitWorker;
 class CProxy_ExplicitWorker;
 class CProxyElement_ExplicitWorker;
 class CProxySection_ExplicitWorker;
/* --------------- index object ------------------ */
class CkIndex_ExplicitWorker:public CProxyElement_ArrayElement{
  public:
    typedef ExplicitWorker local_t;
    typedef CkIndex_ExplicitWorker index_t;
    typedef CProxy_ExplicitWorker proxy_t;
    typedef CProxyElement_ExplicitWorker element_t;
    typedef CProxySection_ExplicitWorker section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: ExplicitWorker(CkMigrateMessage* impl_msg);
 */
    static int __idx_ExplicitWorker_CkMigrateMessage;
    static int ckNew(CkMigrateMessage* impl_msg) { return __idx_ExplicitWorker_CkMigrateMessage; }
    static void _call_ExplicitWorker_CkMigrateMessage(void* impl_msg,ExplicitWorker* impl_obj);

/* DECLS: ExplicitWorker(void);
 */
    static int __idx_ExplicitWorker_void;
    static int ckNew(void) { return __idx_ExplicitWorker_void; }
    static void _call_ExplicitWorker_void(void* impl_msg,ExplicitWorker* impl_obj);

/* DECLS: void precalc(double kappa, unsigned int total_num_patches);
 */
    static int __idx_precalc_marshall2;
    static int precalc(double kappa, unsigned int total_num_patches) { return __idx_precalc_marshall2; }
    static void _call_precalc_marshall2(void* impl_msg,ExplicitWorker* impl_obj);
    static int _callmarshall_precalc_marshall2(char* impl_buf,ExplicitWorker* impl_obj);
    static void _marshallmessagepup_precalc_marshall2(PUP::er &p,void *msg);

/* DECLS: void evaluate(double kappa);
 */
    static int __idx_evaluate_marshall3;
    static int evaluate(double kappa) { return __idx_evaluate_marshall3; }
    static void _call_evaluate_marshall3(void* impl_msg,ExplicitWorker* impl_obj);
    static int _callmarshall_evaluate_marshall3(char* impl_buf,ExplicitWorker* impl_obj);
    static void _marshallmessagepup_evaluate_marshall3(PUP::er &p,void *msg);

/* DECLS: void triggerLB(void);
 */
    static int __idx_triggerLB_void;
    static int triggerLB(void) { return __idx_triggerLB_void; }
    static void _call_triggerLB_void(void* impl_msg,ExplicitWorker* impl_obj);

};
/* --------------- element proxy ------------------ */
 class CProxyElement_ExplicitWorker : public CProxyElement_ArrayElement{
  public:
    typedef ExplicitWorker local_t;
    typedef CkIndex_ExplicitWorker index_t;
    typedef CProxy_ExplicitWorker proxy_t;
    typedef CProxyElement_ExplicitWorker element_t;
    typedef CProxySection_ExplicitWorker section_t;

    CProxyElement_ExplicitWorker(void) {}
    CProxyElement_ExplicitWorker(const ArrayElement *e) : CProxyElement_ArrayElement(e){  }
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
    ExplicitWorker *ckLocal(void) const
      { return (ExplicitWorker *)CProxyElement_ArrayElement::ckLocal(); }
    CProxyElement_ExplicitWorker(const CkArrayID &aid,const CkArrayIndexOctreeIndexer &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS) {}
    CProxyElement_ExplicitWorker(const CkArrayID &aid,const CkArrayIndexOctreeIndexer &idx)
        :CProxyElement_ArrayElement(aid,idx) {}
/* DECLS: ExplicitWorker(CkMigrateMessage* impl_msg);
 */

/* DECLS: ExplicitWorker(void);
 */
    void insert(int onPE=-1);
/* DECLS: void precalc(double kappa, unsigned int total_num_patches);
 */
    void precalc(double kappa, unsigned int total_num_patches, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void evaluate(double kappa);
 */
    void evaluate(double kappa, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void triggerLB(void);
 */
    void triggerLB(void) ;

};
PUPmarshall(CProxyElement_ExplicitWorker)
/* ---------------- collective proxy -------------- */
 class CProxy_ExplicitWorker : public CProxy_ArrayElement{
  public:
    typedef ExplicitWorker local_t;
    typedef CkIndex_ExplicitWorker index_t;
    typedef CProxy_ExplicitWorker proxy_t;
    typedef CProxyElement_ExplicitWorker element_t;
    typedef CProxySection_ExplicitWorker section_t;

    CProxy_ExplicitWorker(void) {}
    CProxy_ExplicitWorker(const ArrayElement *e) : CProxy_ArrayElement(e){  }
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
    CProxyElement_ExplicitWorker operator [] (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_ExplicitWorker(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_ExplicitWorker operator() (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_ExplicitWorker(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxy_ExplicitWorker(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_ExplicitWorker(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: ExplicitWorker(CkMigrateMessage* impl_msg);
 */

/* DECLS: ExplicitWorker(void);
 */
    static CkArrayID ckNew(const CkArrayOptions &opts);

/* DECLS: void precalc(double kappa, unsigned int total_num_patches);
 */
    void precalc(double kappa, unsigned int total_num_patches, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void evaluate(double kappa);
 */
    void evaluate(double kappa, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void triggerLB(void);
 */
    void triggerLB(void) ;

};
PUPmarshall(CProxy_ExplicitWorker)
/* ---------------- section proxy -------------- */
 class CProxySection_ExplicitWorker : public CProxySection_ArrayElement{
  public:
    typedef ExplicitWorker local_t;
    typedef CkIndex_ExplicitWorker index_t;
    typedef CProxy_ExplicitWorker proxy_t;
    typedef CProxyElement_ExplicitWorker element_t;
    typedef CProxySection_ExplicitWorker section_t;

    CProxySection_ExplicitWorker(void) {}
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
    CProxyElement_ExplicitWorker operator [] (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_ExplicitWorker(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_ExplicitWorker operator() (const CkArrayIndexOctreeIndexer &idx) const
        {return CProxyElement_ExplicitWorker(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxySection_ExplicitWorker(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_ExplicitWorker(const CkArrayID &aid, CkArrayIndex *elems, int nElems) 
        :CProxySection_ArrayElement(aid,elems,nElems) {}
    CProxySection_ExplicitWorker(const CkSectionID &sid)       :CProxySection_ArrayElement(sid) {}
    CProxySection_ExplicitWorker(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_ExplicitWorker(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) {}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
/* DECLS: ExplicitWorker(CkMigrateMessage* impl_msg);
 */

/* DECLS: ExplicitWorker(void);
 */

/* DECLS: void precalc(double kappa, unsigned int total_num_patches);
 */
    void precalc(double kappa, unsigned int total_num_patches, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void evaluate(double kappa);
 */
    void evaluate(double kappa, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void triggerLB(void);
 */
    void triggerLB(void) ;

};
PUPmarshall(CProxySection_ExplicitWorker)
typedef CBaseT1<ArrayElementT<OctreeIndexer>, CProxy_ExplicitWorker> CBase_ExplicitWorker;

extern void _registerexplicit_worker(void);
#endif
