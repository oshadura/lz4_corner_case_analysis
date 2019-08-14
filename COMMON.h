#ifndef	__common__
#define	__common__

#define	GEM_N			22
#define TG_N            	22

#define TGR_COUNTER	( (TG_N << 9) | 0 )
#define TGR_COUNTER2	( (TG_N << 9) | 1 )
#define TGT_LAM		( (TG_N << 9) | 8 )
#define TGC_CLRLAM	( (TG_N << 9) | 10)
#define TGW_COUNTER	( (TG_N << 9) | 16)
#define TGW_STATUS	( (TG_N << 9) | 17)
#define TGC_DISLAM	( (TG_N << 9) | 24)
#define TGC_ENLAM	( (TG_N << 9) | 26)
#define TGC_ENABLE	( (TG_N << 9) | 25)
#define TAG_READ	TGR_COUNTER

#define	SINGLES_TIMER	20

#define	IP_CFG_FNAME		"/usr/local/data/crate2ip"


#define MS_CMD_PORT		2000
#define MS_DAT_PORT		2002
#define MS_ERR_PORT		2004
#define DS_CMD_PORT		2006
#define DS_DAT_PORT		2008
#define DS_ERR_PORT		2010

#define TRUE 			1
#define DONE			1
#define FALSE 			0

typedef char* 			charptr;
typedef	unsigned char 		byte;
typedef byte* 			byteptr;

typedef	char			String10[11];	
typedef	char			String20[21];	
typedef	char			String50[51];	
typedef	char			String100[101];	
typedef	char			String200[201];	
typedef	char			String300[301];	

typedef	unsigned short		int16;
typedef	int16*			int16ptr;
typedef	unsigned int		int32;
typedef int32 			boolean;
typedef	int32*			int32ptr;

#define	SADCSIZE		(1 << 13)	/* 13 bit ADC */
#define	SINGBUFSIZE		(SADCSIZE * sizeof (int32))

#define	MAXDS			3
#define DSBUFSIZE		32000	// Decided by the LP RAM size 
#define DS_MEM_SIZE		32512	// FPGA Limit (q-stop + 256 extra)
#define MSBUFSIZE		(DSBUFSIZE * MAXDS)
#define MS_MEM_SIZE		(DS_MEM_SIZE * MAXDS)

#define DS_MAXNAFLIST		100
#define DS_MAXCTLLIST		20
#define DS_MAXINILIST		50
#define	DS_MAXSINGLES		8

#define MAXNAFLIST		(DS_MAXNAFLIST * MAXDS)
#define MAXCTLLIST		(DS_MAXCTLLIST * MAXDS)
#define MAXINILIST		(DS_MAXINILIST * MAXDS)

#define MAXSINGLES		(DS_MAXSINGLES * MAXDS)
#define MAXSCALIST		20

#define NAMESIZE		10		
#define MAXUSERINFO		300
#define ERRMSGSIZE		64
#define	MAX_IP_LEN		64

typedef
    struct ds_info {
	int32		cmdchan, datachan, statchan;

// added for q-Stop START
	int32		words_this_fragment;
// added for q-Stop END
// - ETS
	int32		words_per_fragment;	/* in 2 byte words */
	int32		num_reads;			
	int16		rdnaflist[DS_MAXNAFLIST];  
	int16		rdatamask[DS_MAXNAFLIST];	
	int32 		num_ctls;			
	int16		ctlnaflist[DS_MAXCTLLIST];
	int32		num_singles;
	int16		single_station[DS_MAXSINGLES];

	boolean		waiting_for_data;	/* for individual DSs */
	boolean		received_error_pack;
	
	int16		data[DS_MEM_SIZE];
	int16ptr	dsbp;
	
	String200	status_text;
	char		ds_ip [MAX_IP_LEN];
} ds_info;

enum CamacModes {
	DefaultCamac, HitPatCamac, QstopCamac, MaxCamacModes
};

enum BitPatModes {
	ExecAlways, ExecIfBitSet, ExecIfBitAndInc, EndOfList
};

typedef
    struct {
	int32		command;
	int32		target,result;
	int32		cnaf,data;			

	int32	 	num_singles;	
	int32 		single_crate[MAXSINGLES];
	int32 		single_station[MAXSINGLES];

	int32 		events_per_block;		
	int32 		events_to_send;		

	int32 		eventsize;			
	int32		rdcnaflist[MAXNAFLIST];
	int16		rddatamask[MAXNAFLIST];	

	int32 		num_ctls;			
	int32  		ctlcnaflist[MAXCTLLIST];
	
	int32 		num_scalers;		
	unsigned 	sca_rd_cnaflist[MAXSCALIST];
	unsigned	sca_clr_cnaflist[MAXSCALIST];

	int32 		num_inits;			
	unsigned	init_cnaflist[MAXINILIST];
	unsigned	init_datalist[MAXINILIST];

	char		rd_names[MAXNAFLIST][NAMESIZE];
	char		sing_names[MAXSINGLES][NAMESIZE];
	char		ctl_names[MAXCTLLIST][NAMESIZE];
	char		sca_names[MAXSCALIST][NAMESIZE];
	char		init_names[MAXINILIST][NAMESIZE];

	String300	filename;
	int32		block_size;	/* used by show_status ?? */
	String300	text;		/* user info & error messages */	
} cmdpack,*cmdpackptr;

typedef
    struct {
	int32		command;
	String200	text;	
	int32		result;
	int32		cnaf,data;			

	int32 		events_per_block;		
	int32 		events_to_send;		
	int32 		eventsize;			

	int32 		num_ctls;			
	int32	 	num_singles;	
	int32 		num_scalers;		
	int32 		num_inits;			

	int16		rdnaflist[DS_MAXNAFLIST];
	int16  		ctlnaflist[DS_MAXCTLLIST];

	int16		sca_rd_naflist[MAXSCALIST];
	int16		sca_clr_naflist[MAXSCALIST];

	int16		init_naflist[MAXINILIST];
	int32		init_datalist[MAXINILIST];

	int16		single_station[MAXSINGLES];
} srcmdpack,*srcmdpackptr;


enum	blocktype { event, start, stop, Pause, resume,
hgram, scaler, user, names, calibs, show, error, end_of_file, datarate};

typedef
    struct pilot {
	enum 	blocktype block;		
	int32 	blocksize;			
} pilotpack, pilotpackptr;
	
typedef
    struct dataheader {
	enum 	blocktype block;
	int32 	unitsize;		
	int32 	number_of_units;
	int32 	compstatus;
	int32 	size_in_bytes;		/* these many to follow */
} dataheader,*dataheadptr;

typedef
    struct data {
	dataheader header;		/* order is important */
	int16	data[MSBUFSIZE];
	int32 	empty;
	int32 	size;
} datapack,*datapackptr;

typedef
    struct rate {
	dataheader header;		
	double	evn_during_thisblock;
	double	evn_since_start;
	double	time_for_thisblock;
	double	time_since_start;;
} ratepack,*ratepackptr;

#define OPEN_DATA_FILE		1
#define CLOSE_DATA_FILE		2
#define START_COLLECT		3
#define STOP_COLLECT		4
#define PAUSE_COLLECT		5
#define RESUME_COLLECT		6
#define ENABLE_DISK		7
#define DISABLE_DISK		8
#define SHOW_STATUS		9
#define FSTOP_COLLECT		99

#define TEST_SOCKET		10
#define LAST_REPLY		11	

#define REG_EVENT		20
#define UNREG_EVENT		21

#define	Q_STOP_IGN		25
#define	BITPAT_ENABLE		26

#define ENABLE_SINGLES		30
#define DISABLE_SINGLES		31
#define READ_SINGLES		32
#define CLEAR_SINGLES		33

#define READ_SCALERS		40

#define TERMINATE		50
#define CAMAC_IO		51
#define	ADCINFO			52
#define	TEXTINFO		53

#define ERR_REPORT		60
#define ERR_ACK			61

#define SKIP_RUNS               70
#define SKIP_BLOCKS             71
#define REWIND_FILE             72
#define NEXT_FILE               73
#define PREV_FILE               74

#define	DS_CLOSE_CMDCHAN	81
#define DS_CAMAC_IO		82
#define	READ_DATA		83
#define	REG_LP			84
#define	SWAP_LPBUF		85
#define	ENABLE_INTR		86
#define	DISABLE_INTR		87

#endif
