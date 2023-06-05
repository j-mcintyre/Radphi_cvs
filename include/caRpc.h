
struct CA_BCNA {
	int b;
	int c;
	int n;
	int a;
};
typedef struct CA_BCNA CA_BCNA;
bool_t xdr_CA_BCNA();


struct CA_FXD {
	int f;
	int ext;
	int dat;
};
typedef struct CA_FXD CA_FXD;
bool_t xdr_CA_FXD();


struct CA_DQ {
	int dat;
	int q;
};
typedef struct CA_DQ CA_DQ;
bool_t xdr_CA_DQ();


struct CA_XL {
	int ext;
	int l;
};
typedef struct CA_XL CA_XL;
bool_t xdr_CA_XL();


struct CA_LL {
	int lam;
	int l;
};
typedef struct CA_LL CA_LL;
bool_t xdr_CA_LL();


#define CASRVR ((u_long)0x2c0da009)
#define CAVERS ((u_long)1)
#define CA_RPCCDREG ((u_long)1)
extern int *ca_rpccdreg_1();
#define CA_RPCCFSA ((u_long)2)
extern CA_DQ *ca_rpccfsa_1();
#define CA_RPCCCCZ ((u_long)3)
extern int *ca_rpccccz_1();
#define CA_RPCCCCC ((u_long)4)
extern int *ca_rpccccc_1();
#define CA_RPCCCCI ((u_long)5)
extern int *ca_rpcccci_1();
#define CA_RPCCTCI ((u_long)6)
extern int *ca_rpcctci_1();
#define CA_RPCCCCD ((u_long)7)
extern int *ca_rpccccd_1();
#define CA_RPCCTCD ((u_long)8)
extern int *ca_rpcctcd_1();
#define CA_RPCCTGL ((u_long)9)
extern int *ca_rpcctgl_1();
#define CA_RPCCDLAM ((u_long)10)
extern int *ca_rpccdlam_1();
#define CA_RPCCCLM ((u_long)11)
extern int *ca_rpccclm_1();
#define CA_RPCCCLC ((u_long)12)
extern int *ca_rpccclc_1();
#define CA_RPCCTLM ((u_long)13)
extern int *ca_rpcctlm_1();
#define CA_RPCCCLWT ((u_long)14)
extern int *ca_rpccclwt_1();
#define CA_RPCCCINIT ((u_long)15)
extern int *ca_rpcccinit_1();
#define CA_RPCCTSTAT ((u_long)16)
extern int *ca_rpcctstat_1();

