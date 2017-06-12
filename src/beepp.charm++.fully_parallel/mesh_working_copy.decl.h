#ifndef _DECL_mesh_working_copy_H_
#define _DECL_mesh_working_copy_H_
#include "charm++.h"
/* DECLS: array MeshWorkingCopy: ArrayElement{
MeshWorkingCopy(CkMigrateMessage* impl_msg);
MeshWorkingCopy(void);
void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
void calculate_energy(double kappa, double Dsolvent);
void write_output(const std::string &output_filename);
};
 */
 class MeshWorkingCopy;
 class CkIndex_MeshWorkingCopy;
 class CProxy_MeshWorkingCopy;
 class CProxyElement_MeshWorkingCopy;
 class CProxySection_MeshWorkingCopy;
/* --------------- index object ------------------ */
class CkIndex_MeshWorkingCopy:public CProxyElement_ArrayElement{
  public:
    typedef MeshWorkingCopy local_t;
    typedef CkIndex_MeshWorkingCopy index_t;
    typedef CProxy_MeshWorkingCopy proxy_t;
    typedef CProxyElement_MeshWorkingCopy element_t;
    typedef CProxySection_MeshWorkingCopy section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: MeshWorkingCopy(CkMigrateMessage* impl_msg);
 */
    static int __idx_MeshWorkingCopy_CkMigrateMessage;
    static int ckNew(CkMigrateMessage* impl_msg) { return __idx_MeshWorkingCopy_CkMigrateMessage; }
    static void _call_MeshWorkingCopy_CkMigrateMessage(void* impl_msg,MeshWorkingCopy* impl_obj);

/* DECLS: MeshWorkingCopy(void);
 */
    static int __idx_MeshWorkingCopy_void;
    static int ckNew(void) { return __idx_MeshWorkingCopy_void; }
    static void _call_MeshWorkingCopy_void(void* impl_msg,MeshWorkingCopy* impl_obj);

/* DECLS: void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
 */
    static int __idx_init_marshall2;
    static int init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5) { return __idx_init_marshall2; }
    static void _call_init_marshall2(void* impl_msg,MeshWorkingCopy* impl_obj);
    static int _callmarshall_init_marshall2(char* impl_buf,MeshWorkingCopy* impl_obj);
    static void _marshallmessagepup_init_marshall2(PUP::er &p,void *msg);

/* DECLS: void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
 */
    static int __idx_calculate_rhs_marshall3;
    static int calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb) { return __idx_calculate_rhs_marshall3; }
    static void _call_calculate_rhs_marshall3(void* impl_msg,MeshWorkingCopy* impl_obj);
    static int _callmarshall_calculate_rhs_marshall3(char* impl_buf,MeshWorkingCopy* impl_obj);
    static void _marshallmessagepup_calculate_rhs_marshall3(PUP::er &p,void *msg);

/* DECLS: void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
 */
    static int __idx_process_returned_eval_pts_Eval_Message;
    static int process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg) { return __idx_process_returned_eval_pts_Eval_Message; }
    static void _call_process_returned_eval_pts_Eval_Message(void* impl_msg,MeshWorkingCopy* impl_obj);

/* DECLS: void calculate_energy(double kappa, double Dsolvent);
 */
    static int __idx_calculate_energy_marshall5;
    static int calculate_energy(double kappa, double Dsolvent) { return __idx_calculate_energy_marshall5; }
    static void _call_calculate_energy_marshall5(void* impl_msg,MeshWorkingCopy* impl_obj);
    static int _callmarshall_calculate_energy_marshall5(char* impl_buf,MeshWorkingCopy* impl_obj);
    static void _marshallmessagepup_calculate_energy_marshall5(PUP::er &p,void *msg);

/* DECLS: void write_output(const std::string &output_filename);
 */
    static int __idx_write_output_marshall6;
    static int write_output(const std::string &output_filename) { return __idx_write_output_marshall6; }
    static void _call_write_output_marshall6(void* impl_msg,MeshWorkingCopy* impl_obj);
    static int _callmarshall_write_output_marshall6(char* impl_buf,MeshWorkingCopy* impl_obj);
    static void _marshallmessagepup_write_output_marshall6(PUP::er &p,void *msg);

};
/* --------------- element proxy ------------------ */
 class CProxyElement_MeshWorkingCopy : public CProxyElement_ArrayElement{
  public:
    typedef MeshWorkingCopy local_t;
    typedef CkIndex_MeshWorkingCopy index_t;
    typedef CProxy_MeshWorkingCopy proxy_t;
    typedef CProxyElement_MeshWorkingCopy element_t;
    typedef CProxySection_MeshWorkingCopy section_t;

    CProxyElement_MeshWorkingCopy(void) {}
    CProxyElement_MeshWorkingCopy(const ArrayElement *e) : CProxyElement_ArrayElement(e){  }
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
    MeshWorkingCopy *ckLocal(void) const
      { return (MeshWorkingCopy *)CProxyElement_ArrayElement::ckLocal(); }
    CProxyElement_MeshWorkingCopy(const CkArrayID &aid,const CkArrayIndex1D &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS) {}
    CProxyElement_MeshWorkingCopy(const CkArrayID &aid,const CkArrayIndex1D &idx)
        :CProxyElement_ArrayElement(aid,idx) {}
/* DECLS: MeshWorkingCopy(CkMigrateMessage* impl_msg);
 */

/* DECLS: MeshWorkingCopy(void);
 */
    void insert(int onPE=-1);
/* DECLS: void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
 */
    void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
 */
    void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
 */
    void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg) ;

/* DECLS: void calculate_energy(double kappa, double Dsolvent);
 */
    void calculate_energy(double kappa, double Dsolvent, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void write_output(const std::string &output_filename);
 */
    void write_output(const std::string &output_filename, const CkEntryOptions *impl_e_opts=NULL) ;

};
PUPmarshall(CProxyElement_MeshWorkingCopy)
/* ---------------- collective proxy -------------- */
 class CProxy_MeshWorkingCopy : public CProxy_ArrayElement{
  public:
    typedef MeshWorkingCopy local_t;
    typedef CkIndex_MeshWorkingCopy index_t;
    typedef CProxy_MeshWorkingCopy proxy_t;
    typedef CProxyElement_MeshWorkingCopy element_t;
    typedef CProxySection_MeshWorkingCopy section_t;

    CProxy_MeshWorkingCopy(void) {}
    CProxy_MeshWorkingCopy(const ArrayElement *e) : CProxy_ArrayElement(e){  }
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
    CProxyElement_MeshWorkingCopy operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_MeshWorkingCopy(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_MeshWorkingCopy operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_MeshWorkingCopy(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_MeshWorkingCopy operator [] (int idx) const 
        {return CProxyElement_MeshWorkingCopy(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxyElement_MeshWorkingCopy operator () (int idx) const 
        {return CProxyElement_MeshWorkingCopy(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxy_MeshWorkingCopy(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_MeshWorkingCopy(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: MeshWorkingCopy(CkMigrateMessage* impl_msg);
 */

/* DECLS: MeshWorkingCopy(void);
 */
    static CkArrayID ckNew(const CkArrayOptions &opts);
    static CkArrayID ckNew(const int s1);

/* DECLS: void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
 */
    void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
 */
    void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
 */
    void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg) ;

/* DECLS: void calculate_energy(double kappa, double Dsolvent);
 */
    void calculate_energy(double kappa, double Dsolvent, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void write_output(const std::string &output_filename);
 */
    void write_output(const std::string &output_filename, const CkEntryOptions *impl_e_opts=NULL) ;

};
PUPmarshall(CProxy_MeshWorkingCopy)
/* ---------------- section proxy -------------- */
 class CProxySection_MeshWorkingCopy : public CProxySection_ArrayElement{
  public:
    typedef MeshWorkingCopy local_t;
    typedef CkIndex_MeshWorkingCopy index_t;
    typedef CProxy_MeshWorkingCopy proxy_t;
    typedef CProxyElement_MeshWorkingCopy element_t;
    typedef CProxySection_MeshWorkingCopy section_t;

    CProxySection_MeshWorkingCopy(void) {}
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
    CProxyElement_MeshWorkingCopy operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_MeshWorkingCopy(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_MeshWorkingCopy operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_MeshWorkingCopy(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_MeshWorkingCopy operator [] (int idx) const 
        {return CProxyElement_MeshWorkingCopy(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
    CProxyElement_MeshWorkingCopy operator () (int idx) const 
        {return CProxyElement_MeshWorkingCopy(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex1D *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, int l, int u, int s) {
      CkVec<CkArrayIndex1D> al;
      for (int i=l; i<=u; i+=s) al.push_back(CkArrayIndex1D(i));
      return CkSectionID(aid, al.getVec(), al.size());
    } 
    CProxySection_MeshWorkingCopy(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_MeshWorkingCopy(const CkArrayID &aid, CkArrayIndex *elems, int nElems) 
        :CProxySection_ArrayElement(aid,elems,nElems) {}
    CProxySection_MeshWorkingCopy(const CkSectionID &sid)       :CProxySection_ArrayElement(sid) {}
    CProxySection_MeshWorkingCopy(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_MeshWorkingCopy(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) {}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
/* DECLS: MeshWorkingCopy(CkMigrateMessage* impl_msg);
 */

/* DECLS: MeshWorkingCopy(void);
 */

/* DECLS: void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5);
 */
    void init(unsigned int impl_noname_0, const Quaternion &impl_noname_1, const Vector &impl_noname_2, unsigned int impl_noname_3, double impl_noname_4, double impl_noname_5, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb);
 */
    void calculate_rhs(const CProxy_ParallelFMMOctree<Charge,vanilla_fmm::CProxy_VanillaFMMWorker > &ParallelFMMOctreeProxy, const vanilla_fmm::CProxy_VanillaFMMWorker &VanillaFMMWorkerProxy, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg);
 */
    void process_returned_eval_pts(vanilla_fmm::Eval_Message* impl_msg) ;

/* DECLS: void calculate_energy(double kappa, double Dsolvent);
 */
    void calculate_energy(double kappa, double Dsolvent, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void write_output(const std::string &output_filename);
 */
    void write_output(const std::string &output_filename, const CkEntryOptions *impl_e_opts=NULL) ;

};
PUPmarshall(CProxySection_MeshWorkingCopy)
typedef CBaseT1<ArrayElementT<CkIndex1D>, CProxy_MeshWorkingCopy> CBase_MeshWorkingCopy;

extern void _registermesh_working_copy(void);
#endif
