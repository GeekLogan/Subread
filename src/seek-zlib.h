#ifndef __SEEK_ZLIB_H_
#define __SEEK_ZLIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "subread.h"

// returns 0 if OK; returns 1 if the file is not indexable; returns -1 if file doesn't exist.
int seekgz_open(const char * fname, seekable_zfile_t * fp);

// returns length in bytes if OK (length includes the line break at the end); returns 0 if EOF
int seekgz_gets(seekable_zfile_t * fp, char * buf, int buf_size);

void seekgz_tell(seekable_zfile_t * fp, seekable_position_t * pos);

void seekgz_seek(seekable_zfile_t * fp, seekable_position_t * pos);

int seekgz_next_char(seekable_zfile_t * fp);

void seekgz_close(seekable_zfile_t * fp);
#endif
