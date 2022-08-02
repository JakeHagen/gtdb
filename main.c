#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>

#include "rocksdb/c.h"
//#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/kstring.h"

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

typedef struct _args_t
{
	char *command;
	char *database;
	char *input;
	char *output;
	int threads;
	int compression;

	char **argv;
	int argc;
}
args_t;

static void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}

//make key in the format chromposrefalt
kstring_t mk_key(bcf1_t *line) {
	kstring_t key = {0,0,0};
	bcf_unpack(line, BCF_UN_STR);
        kputw(line->rid, &key);
        kputw(line->pos, &key);
        int i;
        for (i = 0; i < line->n_allele; i++) {
        	kputs(line->d.allele[i], &key);
	}
        return key;
}

kstring_t mk_key_nfmts(bcf1_t *line) {
	kstring_t key = {0,0,0};
        bcf_unpack(line, BCF_UN_STR);
        kputw(line->rid, &key);
        kputw(line->pos, &key);
        int i;
        for (i = 0; i < line->n_allele; i++) {
                kputs(line->d.allele[i], &key);
        }
	kputc('-', &key);
        return key;
}

int verify_fmt_keys(bcf_hdr_t *hdr, bcf_hdr_t *fhdr) {
          khint_t k;
          for (k=0; k<*fhdr->n; k++) {
		  if (strcmp(hdr->id[0][k].key, fhdr->id[0][k].key)) {
			  fprintf(stderr, "%s != %s\n", hdr->id[0][k].key, fhdr->id[0][k].key);
			  return -1;
		  }

          }
          return 0;
  }

kstring_t mk_hdr_path(char *db_path) {
	kstring_t path = {0,0,0};
	kputs(db_path, &path);
	kputs("/hdr.bcf", &path);
	return path;
}

int set_write_options(rocksdb_options_t *options, args_t *args) {
	rocksdb_options_optimize_level_style_compaction(options, 0);
        rocksdb_options_set_enable_blob_files(options, 1);
        rocksdb_options_set_create_if_missing(options, 1);
        rocksdb_options_increase_parallelism(options, args->threads);
        rocksdb_options_set_compression(options, args->compression);
        rocksdb_options_set_blob_compression_type(options, args->compression);
        rocksdb_options_prepare_for_bulk_load(options);
        rocksdb_options_set_write_buffer_size(options, 100000000);

	return 0;
}

uint32_t get_nfmt(bcf_hdr_t *hdr) {
	uint32_t c = 0;
        int i = 0;
	while ( i < hdr->nhrec ) {
		if ( hdr->hrec[i]->type == BCF_HL_FMT ) { c++; }
                i++;
	}
	return c;
}

int store(args_t *args) {
	rocksdb_t *db;
	rocksdb_options_t *options = rocksdb_options_create();
        char *err = NULL;
        
	//open dbq
	set_write_options(options, args);
	db = rocksdb_open(options, args->database, &err);
	assert(!err);

	htsFile *ihf = hts_open(args->input, "r");
	
	hts_set_threads(ihf, args->threads);
	bcf_hdr_t *ihdr = bcf_hdr_read(ihf);
	//uint32_t nfmt = get_nfmt(ihdr);

	//write hdr into db directory
	//open hdr htsfile and write
	kstring_t ffp = mk_hdr_path(args->database);
	htsFile *fhf = hts_open(ks_str(&ffp), "wb");
	if ( bcf_hdr_write(fhf, ihdr)!=0 ) error("Error: cannot write to %s\n", ffp);

	//open output
	htsFile *ohf = hts_open(args->output, "wb");
	hts_set_threads(ohf, args->threads);
	bcf_hdr_t *ohdr = bcf_hdr_subset(ihdr, 0, 0, 0);
	//bcf_hdr_remove(ohdr, BCF_HL_FMT, NULL);
	if ( bcf_hdr_write(ohf, ohdr)!=0 ) error("Error: cannot write output header\n");

	bcf1_t *line = bcf_init();

	rocksdb_writeoptions_t *writeoptions = rocksdb_writeoptions_create();
	
	while (bcf_read(ihf, ihdr, line) == 0) {
		kstring_t key = mk_key(line);

		rocksdb_put(db, writeoptions, ks_str(&key), ks_len(&key), line->indiv.s, line->indiv.l, &err);
		assert(!err);

		kstring_t n_fmt = {0,0,0};
		kputw(line->n_fmt, &n_fmt);
		kstring_t fkey = mk_key_nfmts(line);
		rocksdb_put(db, writeoptions, ks_str(&fkey), ks_len(&fkey), n_fmt.s, n_fmt.l, &err);
		assert(!err);

		//bcf_update_info(ohdr, line, "nfmt", line->n_fmt, 1, BCF_HT_INT);
		bcf_subset(ohdr, line, 0, 0);
		bcf_write(ohf, ohdr, line);
	}

	rocksdb_writeoptions_destroy(writeoptions);
	rocksdb_options_destroy(options);
	rocksdb_close(db);

	bcf_hdr_destroy(ihdr);
	bcf_hdr_destroy(ohdr);
	hts_close(ihf);
	hts_close(fhf);
	hts_close(ohf);

	bcf_destroy(line);

	ks_free(&ffp);

	return 0;
}

void get(args_t *args) {
	//open input and get header, input should not have any individual sample fields (format fields) or samples,
	//but it should have format keys in header unless it was kept in bcf format while detached from format fields
	htsFile *hf = hts_open(args->input, "r");
	hts_set_threads(hf, args->threads);
	bcf_hdr_t *hdr = bcf_hdr_read(hf);

	//open hdr file that contains format fields and samples
	kstring_t ffp = mk_hdr_path(args->database);
        htsFile *fhf = hts_open(ks_str(&ffp), "r");
	bcf_hdr_t *fhdr = bcf_hdr_read(fhf);

	//check that the link between info and format keys and ids has been preserved
	if ( verify_fmt_keys(hdr, fhdr)!=0 ) error("Error: position of format keys does not match\n");

	//uint32_t nfmt = get_nfmt(fhdr); 

	bcf_hdr_merge(fhdr, hdr);
        bcf_hdr_sync(fhdr);

	htsFile *out = hts_open(args->output, "wb");
	hts_set_threads(out, args->threads);
	bcf_hdr_write(out, fhdr);

	bcf1_t *line = bcf_init();

	rocksdb_t *db;
        rocksdb_options_t *options = rocksdb_options_create();
        char *err = NULL;
        rocksdb_options_optimize_level_style_compaction(options, 0);
        rocksdb_options_set_enable_blob_files(options, 1);
        rocksdb_options_set_create_if_missing(options, 1);
        rocksdb_options_increase_parallelism(options, args->threads);
        rocksdb_options_set_compression(options, args->compression);
        rocksdb_options_set_blob_compression_type(options, args->compression);
        db = rocksdb_open_for_read_only(options, args->database, 1, &err);
	assert(!err);

	rocksdb_readoptions_t *readoptions = rocksdb_readoptions_create();

	//char *returned_value;
	while (bcf_read(hf, hdr, line) == 0) {
		free(line->indiv.s);
		kstring_t key = mk_key(line);

		line->indiv.s = rocksdb_get(db, readoptions, ks_str(&key), ks_len(&key), &line->indiv.l, &err);
		assert(!err);

		kstring_t fkey = mk_key_nfmts(line);
		kstring_t nfmtsc = {0,0,0};
		nfmtsc.s = rocksdb_get(db, readoptions, ks_str(&fkey), ks_len(&fkey), &nfmtsc.l, &err);
		line->n_fmt = atoi(ks_str(&nfmtsc));
		//fprintf(stderr, "%s\n", ks_str(&nfmtsc));
		line->n_sample =  bcf_hdr_nsamples(fhdr);
		bcf_write(out, fhdr, line);
	}

	rocksdb_readoptions_destroy(readoptions);
        rocksdb_options_destroy(options);
        rocksdb_close(db);
	
	bcf_hdr_destroy(hdr);
	bcf_hdr_destroy(fhdr);
        
	hts_close(hf);
	hts_close(out);
        hts_close(fhf);
        bcf_destroy(line);

        ks_free(&ffp);
}

int main(int argc, char **argv) {

	int c;
	args_t *args = (args_t*) calloc(1, sizeof(args_t));
	args->argc = argc;
	args->argv = argv;
	args->output = "-";
	args->input = "-";
	args->threads = 1;
	args->compression = 7;

	static struct option loptions[] =
	{
		{"database", required_argument, NULL, 'd'},
		{"output", required_argument, NULL, 'o'},
		{"threads", required_argument, NULL, 't'},
		{"compression", required_argument, NULL, 'c'},
	};

	while ((c = getopt_long(argc, argv, "d:o:t:c:", loptions, NULL)) >= 0)
	{
		switch (c) {
			case 'd': args->database = optarg; break;
			case 'o': args->output = optarg; break;
			case 't': args->threads = strtol(optarg, 0, 0); break;
			case 'c': args->compression = strtol(optarg, 0, 0); break;
		}
	}

	if (!args->database) { error("must give database location\n"); }

	if (argc < 2) { error("supply args"); }

	if (!strcmp(argv[1], "store") && !strcmp(argv[1], "get")) { error("command mut be store or get\n"); }
	optind++;

	int count = 0;
	while ( optind < argc )
	{
		if ( count > 1 ) { error("too many positonal arguments, need subcommand and one file name\n"); }
		args->input = argv[optind];
		optind++;
		count++;
	}

	if (strcmp(argv[argc-2], "store") == 0)
	{
		store(args);
	}
	if (strcmp(argv[argc-2], "get") == 0)
        {
        	get(args);
        }



  	return 0;
}
