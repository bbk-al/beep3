#ifndef _DECL_gmres_H_
#define _DECL_gmres_H_
#include "charm++.h"
/* DECLS: chare GMRES: Chare{
GMRES(void);
threaded sync void solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs);
};
 */
 class GMRES;
 class CkIndex_GMRES;
 class CProxy_GMRES;
/* --------------- index object ------------------ */
class CkIndex_GMRES:public CProxy_Chare{
  public:
    typedef GMRES local_t;
    typedef CkIndex_GMRES index_t;
    typedef CProxy_GMRES proxy_t;
    typedef CProxy_GMRES element_t;

    static int __idx;
    static void __register(const char *s, size_t size);
/* DECLS: GMRES(void);
 */
    static int __idx_GMRES_void;
    static int ckNew(void) { return __idx_GMRES_void; }
    static void _call_GMRES_void(void* impl_msg,GMRES* impl_obj);

/* DECLS: threaded sync void solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs);
 */
    static int __idx_solve_marshall2;
    static int solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs) { return __idx_solve_marshall2; }
    static void _call_solve_marshall2(void* impl_msg,GMRES* impl_obj);
    static void _callthr_solve_marshall2(CkThrCallArg *);
    static void _marshallmessagepup_solve_marshall2(PUP::er &p,void *msg);

};
/* --------------- element proxy ------------------ */
class CProxy_GMRES:public CProxy_Chare{
  public:
    typedef GMRES local_t;
    typedef CkIndex_GMRES index_t;
    typedef CProxy_GMRES proxy_t;
    typedef CProxy_GMRES element_t;

    CProxy_GMRES(void) {};
    CProxy_GMRES(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_GMRES(const Chare *c) : CProxy_Chare(c){  }
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
    GMRES *ckLocal(void) const
     { return (GMRES *)CkLocalChare(&ckGetChareID()); }
/* DECLS: GMRES(void);
 */
    static CkChareID ckNew(int onPE=CK_PE_ANY);
    static void ckNew(CkChareID* pcid, int onPE=CK_PE_ANY);

/* DECLS: threaded sync void solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs);
 */
    void solve(const std::string &output_filename, double kappa, double Dsolvent, unsigned int num_patches, const FH_Values &rhs, const CkEntryOptions *impl_e_opts=NULL);

};
PUPmarshall(CProxy_GMRES)
typedef CBaseT1<Chare, CProxy_GMRES> CBase_GMRES;

extern void _registergmres(void);
#endif
