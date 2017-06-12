#ifndef _DECL_rhs_handler_H_
#define _DECL_rhs_handler_H_
#include "charm++.h"
/* DECLS: chare RHS_Handler: Chare{
RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length);
void add_charges(const std::vector<Charge > &charges);
threaded void get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb);
void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg);
};
 */
 class RHS_Handler;
 class CkIndex_RHS_Handler;
 class CProxy_RHS_Handler;
/* --------------- index object ------------------ */
class CkIndex_RHS_Handler:public CProxy_Chare{
  public:
    typedef RHS_Handler local_t;
    typedef CkIndex_RHS_Handler index_t;
    typedef CProxy_RHS_Handler proxy_t;
    typedef CProxy_RHS_Handler element_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length);
 */
    static int __idx_RHS_Handler_marshall1;
    static int ckNew(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length) { return __idx_RHS_Handler_marshall1; }
    static void _call_RHS_Handler_marshall1(void* impl_msg,RHS_Handler* impl_obj);
    static int _callmarshall_RHS_Handler_marshall1(char* impl_buf,RHS_Handler* impl_obj);
    static void _marshallmessagepup_RHS_Handler_marshall1(PUP::er &p,void *msg);

/* DECLS: void add_charges(const std::vector<Charge > &charges);
 */
    static int __idx_add_charges_marshall2;
    static int add_charges(const std::vector<Charge > &charges) { return __idx_add_charges_marshall2; }
    static void _call_add_charges_marshall2(void* impl_msg,RHS_Handler* impl_obj);
    static int _callmarshall_add_charges_marshall2(char* impl_buf,RHS_Handler* impl_obj);
    static void _marshallmessagepup_add_charges_marshall2(PUP::er &p,void *msg);

/* DECLS: threaded void get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb);
 */
    static int __idx_get_rhs_marshall3;
    static int get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb) { return __idx_get_rhs_marshall3; }
    static void _call_get_rhs_marshall3(void* impl_msg,RHS_Handler* impl_obj);
    static void _callthr_get_rhs_marshall3(CkThrCallArg *);
    static void _marshallmessagepup_get_rhs_marshall3(PUP::er &p,void *msg);

/* DECLS: void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg);
 */
    static int __idx_process_rhs_results_from_working_mesh_Eval_Message;
    static int process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg) { return __idx_process_rhs_results_from_working_mesh_Eval_Message; }
    static void _call_process_rhs_results_from_working_mesh_Eval_Message(void* impl_msg,RHS_Handler* impl_obj);

};
/* --------------- element proxy ------------------ */
class CProxy_RHS_Handler:public CProxy_Chare{
  public:
    typedef RHS_Handler local_t;
    typedef CkIndex_RHS_Handler index_t;
    typedef CProxy_RHS_Handler proxy_t;
    typedef CProxy_RHS_Handler element_t;

    CProxy_RHS_Handler(void) {};
    CProxy_RHS_Handler(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_RHS_Handler(const Chare *c) : CProxy_Chare(c){  }
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
    RHS_Handler *ckLocal(void) const
     { return (RHS_Handler *)CkLocalChare(&ckGetChareID()); }
/* DECLS: RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length);
 */
    static CkChareID ckNew(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);
    static void ckNew(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, CkChareID* pcid, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_RHS_Handler(unsigned int max_items_per_node, const Vector &universe_centre, double universe_edge_length, int onPE=CK_PE_ANY, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void add_charges(const std::vector<Charge > &charges);
 */
    void add_charges(const std::vector<Charge > &charges, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: threaded void get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb);
 */
    void get_rhs(double Dsolvent, unsigned int impl_noname_0, const CkCallback &cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg);
 */
    void process_rhs_results_from_working_mesh(vanilla_fmm::Eval_Message* impl_msg);

};
PUPmarshall(CProxy_RHS_Handler)
typedef CBaseT1<Chare, CProxy_RHS_Handler> CBase_RHS_Handler;

/* DECLS: message RHS_Message{
double data[];
}
;
 */
class RHS_Message;
class CMessage_RHS_Message:public CkMessage{
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
    CMessage_RHS_Message();
    static void *pack(RHS_Message *p);
    static RHS_Message* unpack(void* p);
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

extern void _registerrhs_handler(void);
#endif
