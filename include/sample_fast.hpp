#ifndef SAMPLE_FAST_HPP
#define SAMPLE_FAST_HPP
#include <fstream>
#include <string>
#include <iostream>
#include "zlib.h"
#include <vector>
#include <algorithm>
#include "get_time.hpp"
#include "strip_split_join.hpp"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

using namespace std;


#ifndef __AC_KHASH_H
#define __AC_KHASH_H

/*!
  @header
  Generic hash table library.
 */

#define AC_VERSION_KHASH_H "0.2.6"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* compipler specific configuration */

#if UINT_MAX == 0xffffffffu
typedef unsigned int khint32_t;
#elif ULONG_MAX == 0xffffffffu
typedef unsigned long khint32_t;
#endif

#if ULONG_MAX == ULLONG_MAX
typedef unsigned long khint64_t;
#else
typedef unsigned long long khint64_t;
#endif

#ifdef _MSC_VER
#define inline __inline
#endif

typedef khint32_t khint_t;
typedef khint_t khiter_t;

#define __ac_isempty(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&2)
#define __ac_isdel(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&1)
#define __ac_iseither(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&3)
#define __ac_set_isdel_false(flag, i) (flag[i>>4]&=~(1ul<<((i&0xfU)<<1)))
#define __ac_set_isempty_false(flag, i) (flag[i>>4]&=~(2ul<<((i&0xfU)<<1)))
#define __ac_set_isboth_false(flag, i) (flag[i>>4]&=~(3ul<<((i&0xfU)<<1)))
#define __ac_set_isdel_true(flag, i) (flag[i>>4]|=1ul<<((i&0xfU)<<1))

#ifdef KHASH_LINEAR
#define __ac_inc(k, m) 1
#else
#define __ac_inc(k, m) (((k)>>3 ^ (k)<<3) | 1) & (m)
#endif

#define __ac_fsize(m) ((m) < 16? 1 : (m)>>4)

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static const double __ac_HASH_UPPER = 0.77;

#define __KHASH_TYPE(name, khkey_t, khval_t) \
	typedef struct { \
		khint_t n_buckets, size, n_occupied, upper_bound; \
		khint32_t *flags; \
		khkey_t *keys; \
		khval_t *vals; \
	} kh_##name##_t;

#define KHASH_DECLARE(name, khkey_t, khval_t)		 					\
	__KHASH_TYPE(name, khkey_t, khval_t) 								\
	extern kh_##name##_t *kh_init_##name();								\
	extern void kh_destroy_##name(kh_##name##_t *h);					\
	extern void kh_clear_##name(kh_##name##_t *h);						\
	extern khint_t kh_get_##name(const kh_##name##_t *h, khkey_t key); 	\
	extern void kh_resize_##name(kh_##name##_t *h, khint_t new_n_buckets); \
	extern khint_t kh_put_##name(kh_##name##_t *h, khkey_t key, int *ret); \
	extern void kh_del_##name(kh_##name##_t *h, khint_t x);

#define KHASH_INIT2(name, SCOPE, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal) \
	__KHASH_TYPE(name, khkey_t, khval_t) 								\
	SCOPE kh_##name##_t *kh_init_##name() {								\
		return (kh_##name##_t*)calloc(1, sizeof(kh_##name##_t));		\
	}																	\
	SCOPE void kh_destroy_##name(kh_##name##_t *h)						\
	{																	\
		if (h) {														\
			free(h->keys); free(h->flags);								\
			free(h->vals);												\
			free(h);													\
		}																\
	}																	\
	SCOPE void kh_clear_##name(kh_##name##_t *h)						\
	{																	\
		if (h && h->flags) {											\
			memset(h->flags, 0xaa, __ac_fsize(h->n_buckets) * sizeof(khint32_t)); \
			h->size = h->n_occupied = 0;								\
		}																\
	}																	\
	SCOPE khint_t kh_get_##name(const kh_##name##_t *h, khkey_t key) 	\
	{																	\
		if (h->n_buckets) {												\
			khint_t inc, k, i, last, mask;								\
			mask = h->n_buckets - 1;									\
			k = __hash_func(key); i = k & mask;							\
			inc = __ac_inc(k, mask); last = i; /* inc==1 for linear probing */ \
			while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
				i = (i + inc) & mask; 									\
				if (i == last) return h->n_buckets;						\
			}															\
			return __ac_iseither(h->flags, i)? h->n_buckets : i;		\
		} else return 0;												\
	}																	\
	SCOPE void kh_resize_##name(kh_##name##_t *h, khint_t new_n_buckets) \
	{ /* This function uses 0.25*n_bucktes bytes of working space instead of [sizeof(key_t+val_t)+.25]*n_buckets. */ \
		khint32_t *new_flags = 0;										\
		khint_t j = 1;													\
		{																\
			kroundup32(new_n_buckets); 									\
			if (new_n_buckets < 4) new_n_buckets = 4;					\
			if (h->size >= (khint_t)(new_n_buckets * __ac_HASH_UPPER + 0.5)) j = 0;	/* requested size is too small */ \
			else { /* hash table size to be changed (shrink or expand); rehash */ \
				new_flags = (khint32_t*)malloc(__ac_fsize(new_n_buckets) * sizeof(khint32_t));	\
				memset(new_flags, 0xaa, __ac_fsize(new_n_buckets) * sizeof(khint32_t)); \
				if (h->n_buckets < new_n_buckets) {	/* expand */		\
					h->keys = (khkey_t*)realloc(h->keys, new_n_buckets * sizeof(khkey_t)); \
					if (kh_is_map) h->vals = (khval_t*)realloc(h->vals, new_n_buckets * sizeof(khval_t)); \
				} /* otherwise shrink */								\
			}															\
		}																\
		if (j) { /* rehashing is needed */								\
			for (j = 0; j != h->n_buckets; ++j) {						\
				if (__ac_iseither(h->flags, j) == 0) {					\
					khkey_t key = h->keys[j];							\
					khval_t val;										\
					khint_t new_mask;									\
					new_mask = new_n_buckets - 1; 						\
					if (kh_is_map) val = h->vals[j];					\
					__ac_set_isdel_true(h->flags, j);					\
					while (1) { /* kick-out process; sort of like in Cuckoo hashing */ \
						khint_t inc, k, i;								\
						k = __hash_func(key);							\
						i = k & new_mask;								\
						inc = __ac_inc(k, new_mask);					\
						while (!__ac_isempty(new_flags, i)) i = (i + inc) & new_mask; \
						__ac_set_isempty_false(new_flags, i);			\
						if (i < h->n_buckets && __ac_iseither(h->flags, i) == 0) { /* kick out the existing element */ \
							{ khkey_t tmp = h->keys[i]; h->keys[i] = key; key = tmp; } \
							if (kh_is_map) { khval_t tmp = h->vals[i]; h->vals[i] = val; val = tmp; } \
							__ac_set_isdel_true(h->flags, i); /* mark it as deleted in the old hash table */ \
						} else { /* write the element and jump out of the loop */ \
							h->keys[i] = key;							\
							if (kh_is_map) h->vals[i] = val;			\
							break;										\
						}												\
					}													\
				}														\
			}															\
			if (h->n_buckets > new_n_buckets) { /* shrink the hash table */ \
				h->keys = (khkey_t*)realloc(h->keys, new_n_buckets * sizeof(khkey_t)); \
				if (kh_is_map) h->vals = (khval_t*)realloc(h->vals, new_n_buckets * sizeof(khval_t)); \
			}															\
			free(h->flags); /* free the working space */				\
			h->flags = new_flags;										\
			h->n_buckets = new_n_buckets;								\
			h->n_occupied = h->size;									\
			h->upper_bound = (khint_t)(h->n_buckets * __ac_HASH_UPPER + 0.5); \
		}																\
	}																	\
	SCOPE khint_t kh_put_##name(kh_##name##_t *h, khkey_t key, int *ret) \
	{																	\
		khint_t x;														\
		if (h->n_occupied >= h->upper_bound) { /* update the hash table */ \
			if (h->n_buckets > (h->size<<1)) kh_resize_##name(h, h->n_buckets - 1); /* clear "deleted" elements */ \
			else kh_resize_##name(h, h->n_buckets + 1); /* expand the hash table */ \
		} /* TODO: to implement automatically shrinking; resize() already support shrinking */ \
		{																\
			khint_t inc, k, i, site, last, mask = h->n_buckets - 1;		\
			x = site = h->n_buckets; k = __hash_func(key); i = k & mask; \
			if (__ac_isempty(h->flags, i)) x = i; /* for speed up */	\
			else {														\
				inc = __ac_inc(k, mask); last = i;						\
				while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
					if (__ac_isdel(h->flags, i)) site = i;				\
					i = (i + inc) & mask; 								\
					if (i == last) { x = site; break; }					\
				}														\
				if (x == h->n_buckets) {								\
					if (__ac_isempty(h->flags, i) && site != h->n_buckets) x = site; \
					else x = i;											\
				}														\
			}															\
		}																\
		if (__ac_isempty(h->flags, x)) { /* not present at all */		\
			h->keys[x] = key;											\
			__ac_set_isboth_false(h->flags, x);							\
			++h->size; ++h->n_occupied;									\
			*ret = 1;													\
		} else if (__ac_isdel(h->flags, x)) { /* deleted */				\
			h->keys[x] = key;											\
			__ac_set_isboth_false(h->flags, x);							\
			++h->size;													\
			*ret = 2;													\
		} else *ret = 0; /* Don't touch h->keys[x] if present and not deleted */ \
		return x;														\
	}																	\
	SCOPE void kh_del_##name(kh_##name##_t *h, khint_t x)				\
	{																	\
		if (x != h->n_buckets && !__ac_iseither(h->flags, x)) {			\
			__ac_set_isdel_true(h->flags, x);							\
			--h->size;													\
		}																\
	}

#define KHASH_INIT(name, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal) \
	KHASH_INIT2(name, static inline, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal)

/* --- BEGIN OF HASH FUNCTIONS --- */

/*! @function
  @abstract     Integer hash function
  @param  key   The integer [khint32_t]
  @return       The hash value [khint_t]
 */
#define kh_int_hash_func(key) (khint32_t)(key)
/*! @function
  @abstract     Integer comparison function
 */
#define kh_int_hash_equal(a, b) ((a) == (b))
/*! @function
  @abstract     64-bit integer hash function
  @param  key   The integer [khint64_t]
  @return       The hash value [khint_t]
 */
#define kh_int64_hash_func(key) (khint32_t)((key)>>33^(key)^(key)<<11)
/*! @function
  @abstract     64-bit integer comparison function
 */
#define kh_int64_hash_equal(a, b) ((a) == (b))
/*! @function
  @abstract     const char* hash function
  @param  s     Pointer to a null terminated string
  @return       The hash value
 */
static inline khint_t __ac_X31_hash_string(const char *s)
{
	khint_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}
/*! @function
  @abstract     Another interface to const char* hash function
  @param  key   Pointer to a null terminated string [const char*]
  @return       The hash value [khint_t]
 */
#define kh_str_hash_func(key) __ac_X31_hash_string(key)
/*! @function
  @abstract     Const char* comparison function
 */
#define kh_str_hash_equal(a, b) (strcmp(a, b) == 0)

static inline khint_t __ac_Wang_hash(khint_t key)
{
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}
#define kh_int_hash_func2(k) __ac_Wang_hash((khint_t)key)

/* --- END OF HASH FUNCTIONS --- */

/* Other convenient macros... */

/*!
  @abstract Type of the hash table.
  @param  name  Name of the hash table [symbol]
 */
#define khash_t(name) kh_##name##_t

/*! @function
  @abstract     Initiate a hash table.
  @param  name  Name of the hash table [symbol]
  @return       Pointer to the hash table [khash_t(name)*]
 */
#define kh_init(name) kh_init_##name()

/*! @function
  @abstract     Destroy a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
 */
#define kh_destroy(name, h) kh_destroy_##name(h)

/*! @function
  @abstract     Reset a hash table without deallocating memory.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
 */
#define kh_clear(name, h) kh_clear_##name(h)

/*! @function
  @abstract     Resize a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  s     New size [khint_t]
 */
#define kh_resize(name, h, s) kh_resize_##name(h, s)

/*! @function
  @abstract     Insert a key to the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Key [type of keys]
  @param  r     Extra return code: 0 if the key is present in the hash table;
                1 if the bucket is empty (never used); 2 if the element in
				the bucket has been deleted [int*]
  @return       Iterator to the inserted element [khint_t]
 */
#define kh_put(name, h, k, r) kh_put_##name(h, k, r)

/*! @function
  @abstract     Retrieve a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Key [type of keys]
  @return       Iterator to the found element, or kh_end(h) is the element is absent [khint_t]
 */
#define kh_get(name, h, k) kh_get_##name(h, k)

/*! @function
  @abstract     Remove a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Iterator to the element to be deleted [khint_t]
 */
#define kh_del(name, h, k) kh_del_##name(h, k)

/*! @function
  @abstract     Test whether a bucket contains data.
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       1 if containing data; 0 otherwise [int]
 */
#define kh_exist(h, x) (!__ac_iseither((h)->flags, (x)))

/*! @function
  @abstract     Get key given an iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       Key [type of keys]
 */
#define kh_key(h, x) ((h)->keys[x])

/*! @function
  @abstract     Get value given an iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       Value [type of values]
  @discussion   For hash sets, calling this results in segfault.
 */
#define kh_val(h, x) ((h)->vals[x])

/*! @function
  @abstract     Alias of kh_val()
 */
#define kh_value(h, x) ((h)->vals[x])

/*! @function
  @abstract     Get the start iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       The start iterator [khint_t]
 */
#define kh_begin(h) (khint_t)(0)

/*! @function
  @abstract     Get the end iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       The end iterator [khint_t]
 */
#define kh_end(h) ((h)->n_buckets)

/*! @function
  @abstract     Get the number of elements in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       Number of elements in the hash table [khint_t]
 */
#define kh_size(h) ((h)->size)

/*! @function
  @abstract     Get the number of buckets in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       Number of buckets in the hash table [khint_t]
 */
#define kh_n_buckets(h) ((h)->n_buckets)

/* More conenient interfaces */

/*! @function
  @abstract     Instantiate a hash set containing integer keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_INT(name)										\
	KHASH_INIT(name, khint32_t, char, 0, kh_int_hash_func, kh_int_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing integer keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_INT(name, khval_t)								\
	KHASH_INIT(name, khint32_t, khval_t, 1, kh_int_hash_func, kh_int_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_INT64(name)										\
	KHASH_INIT(name, khint64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_INT64(name, khval_t)								\
	KHASH_INIT(name, khint64_t, khval_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

typedef const char *kh_cstr_t;
/*! @function
  @abstract     Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_STR(name)										\
	KHASH_INIT(name, kh_cstr_t, char, 0, kh_str_hash_func, kh_str_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_STR(name, khval_t)								\
	KHASH_INIT(name, kh_cstr_t, khval_t, 1, kh_str_hash_func, kh_str_hash_equal)

#endif /* __AC_KHASH_H */

namespace sample_fast
{
    /* sample */

    typedef uint64_t krint64_t;

    typedef struct {
	int n, m;
	uint64_t *a;
	int8_t *rev;
    } reglist_t;

    KHASH_MAP_INIT_STR(reg, reglist_t)
    KHASH_SET_INIT_INT64(64)

    typedef kh_reg_t reghash_t;

    reghash_t *stk_reg_read_alt(const char *fn)
    {
        reghash_t *h = kh_init(reg);
        gzFile fp;
        kstream_t *ks;
        int dret;
        kstring_t str = {0,0,0};

        fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
        if (fp == 0) return 0;
        ks = ks_init(fp);
        while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
            int i, c, st = -1, en = -1, rev = 0;
            char *p, *q;
            reglist_t *r;
            khint_t k;
            for (i = 0, p = q = str.s;; ++p) {
                if (*p == '\t' || *p == '\0') {
                    c = *p;
                    *p = 0;
                    if (i == 1) {
                        if (isdigit(*q)) st = strtol(q, &q, 10);
                        if (q != p) st = -1;
                    } else if (i == 2) {
                        if (isdigit(*q)) en = strtol(q, &q, 10);
                        if (q != p) en = -1;
                    } else if (i == 5) {
                        if (*q == '+') rev = 1;
                        else if (*q == '-') rev = -1;
                    }
                    ++i, q = p + 1;
                    if (c == 0) break;
                }
            }
            if (i == 0) continue;
            k = kh_get(reg, h, str.s);
            if (k == kh_end(h)) {
                int ret;
                char *s = strdup(str.s);
                k = kh_put(reg, h, s, &ret);
                memset(&kh_val(h, k), 0, sizeof(reglist_t));
            }
            r = &kh_val(h, k);
            if (en < 0 && st > 0) en = st, st = st - 1; // if there is only one column
            if (st < 0) st = 0, en = INT_MAX;
            if (r->n == r->m) {
                r->m = r->m? r->m<<1 : 4;
                r->a = (uint64_t*)realloc(r->a, r->m * 8);
                r->rev = (int8_t*)realloc(r->rev, r->m);
            }
            r->a[r->n] = (uint64_t)st<<32 | en;
            r->rev[r->n++] = rev;
        }
        ks_destroy(ks);
        gzclose(fp);
        free(str.s);
        return h;
    }

    reghash_t *stk_reg_read(const char *fn)
    {
        reghash_t *h = kh_init(reg);
        gzFile fp;
        kstream_t *ks;
        int dret;
        kstring_t *str;
        // read the list
        fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
        if (fp == 0) return 0;
        ks = ks_init(fp);
        str = (kstring_t *)calloc(1, sizeof(kstring_t));
        while (ks_getuntil(ks, 0, str, &dret) >= 0) {
            int beg = -1, end = -1;
            reglist_t *p;
            khint_t k = kh_get(reg, h, str->s);
            if (k == kh_end(h)) {
                int ret;
                char *s = strdup(str->s);
                k = kh_put(reg, h, s, &ret);
                memset(&kh_val(h, k), 0, sizeof(reglist_t));
            }
            p = &kh_val(h, k);
            if (dret != '\n') {
                if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
                    beg = atoi(str->s);
                    if (dret != '\n') {
                        if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
                            end = atoi(str->s);
                            if (end < 0) end = -1;
                        }
                    }
                }
            }
            // skip the rest of the line
            if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
            if (end < 0 && beg > 0) end = beg, beg = beg - 1; // if there is only one column
            if (beg < 0) beg = 0, end = INT_MAX;
            if (p->n == p->m) {
                p->m = p->m? p->m<<1 : 4;
                p->a = (uint64_t *)realloc(p->a, p->m * 8);
            }
            p->a[p->n++] = (uint64_t)beg<<32 | end;
        }
        ks_destroy(ks);
        gzclose(fp);
        free(str->s); free(str);
        return h;
    }

    void stk_reg_destroy(reghash_t *h)
    {
        khint_t k;
        if (h == 0) return;
        for (k = 0; k < kh_end(h); ++k) {
            if (kh_exist(h, k)) {
                free(kh_val(h, k).a);
                free(kh_val(h, k).rev);
                free((char*)kh_key(h, k));
            }
        }
        kh_destroy(reg, h);
    }

    /* constant table */

    unsigned char seq_nt16_table[256] = {
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
        15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
        15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
        15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
        15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
    };

    unsigned char seq_nt6_table[256] = {
        0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
    };

    char *seq_nt16_rev_table = (char*)"XACMGRSVTWYHKDBN";
    unsigned char seq_nt16to4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
    unsigned char seq_nt16comp_table[] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
    int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
    char comp_tab[] = {
        0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
        16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
        32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
        48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
        64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
        64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
    };

    static void stk_printstr(FILE *fp, const kstring_t *s, unsigned line_len)
    {
        if (line_len != UINT_MAX && line_len != 0) {
            int i, rest = s->l;
            for (i = 0; i < s->l; i += line_len, rest -= line_len) {
                fputc('\n', fp);
                if (rest > line_len) fwrite(s->s + i, 1, line_len, fp);
                else fwrite(s->s + i, 1, rest, fp);
            }
            fputc('\n', fp);
        } else {
            fputc('\n', fp);
            fputs(s->s, fp);
            fputc('\n', fp);
        }
    }

    static inline void stk_printseq_renamed(FILE *fp, const kseq_t *s, int line_len, const char *prefix, int64_t n)
    {
        fputc(s->qual.l? '@' : '>', fp);
        if (n >= 0) {
            if (prefix) fputs(prefix, fp);
            fprintf(fp, "%lld", (long long)n);
        } else fputs(s->name.s, fp);
        if (s->comment.l) {
            fputc(' ', fp); fputs(s->comment.s, fp);
        }
        stk_printstr(fp, &s->seq, line_len);
        if (s->qual.l) {
            fputc('+', fp);
            stk_printstr(fp, &s->qual, line_len);
        }
    }

    static inline void stk_printseq(FILE *fp, const kseq_t *s, int line_len)
    {
        stk_printseq_renamed(fp, s, line_len, 0, -1);
    }

    /* 
    64-bit Mersenne Twister pseudorandom number generator. Adapted from:
        http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
    which was written by Takuji Nishimura and Makoto Matsumoto and released
    under the 3-clause BSD license.
    */

    struct _krand_t;
    typedef struct _krand_t krand_t;

    #define KR_NN 312
    #define KR_MM 156
    #define KR_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
    #define KR_LM 0x7FFFFFFFULL /* Least significant 31 bits */

    struct _krand_t {
        int mti;
        krint64_t mt[KR_NN];
    };

    static void kr_srand0(krint64_t seed, krand_t *kr)
    {
        kr->mt[0] = seed;
        for (kr->mti = 1; kr->mti < KR_NN; ++kr->mti) 
            kr->mt[kr->mti] = 6364136223846793005ULL * (kr->mt[kr->mti - 1] ^ (kr->mt[kr->mti - 1] >> 62)) + kr->mti;
    }

    krand_t *kr_srand(krint64_t seed)
    {
        krand_t *kr;
        kr = (krand_t *)malloc(sizeof(krand_t));
        kr_srand0(seed, kr);
        return kr;
    }

    krint64_t kr_rand(krand_t *kr)
    {
        krint64_t x;
        static const krint64_t mag01[2] = { 0, 0xB5026F5AA96619E9ULL };
        if (kr->mti >= KR_NN) {
            int i;
            if (kr->mti == KR_NN + 1) kr_srand0(5489ULL, kr);
            for (i = 0; i < KR_NN - KR_MM; ++i) {
                x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
                kr->mt[i] = kr->mt[i + KR_MM] ^ (x>>1) ^ mag01[(int)(x&1)];
            }
            for (; i < KR_NN - 1; ++i) {
                x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
                kr->mt[i] = kr->mt[i + (KR_MM - KR_NN)] ^ (x>>1) ^ mag01[(int)(x&1)];
            }
            x = (kr->mt[KR_NN - 1] & KR_UM) | (kr->mt[0] & KR_LM);
            kr->mt[KR_NN - 1] = kr->mt[KR_MM - 1] ^ (x>>1) ^ mag01[(int)(x&1)];
            kr->mti = 0;
        }
        x = kr->mt[kr->mti++];
        x ^= (x >> 29) & 0x5555555555555555ULL;
        x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
        x ^= (x << 37) & 0xFFF7EEE000000000ULL;
        x ^= (x >> 43);
        return x;
    }

    #define kr_drand(_kr) ((kr_rand(_kr) >> 11) * (1.0/9007199254740992.0))

    static void cpy_kstr(kstring_t *dst, const kstring_t *src)
    {
        if (src->l == 0) return;
        if (src->l + 1 > dst->m) {
            dst->m = src->l + 1;
            kroundup32(dst->m);
            dst->s = (char *)realloc(dst->s, dst->m);
        }
        dst->l = src->l;
        memcpy(dst->s, src->s, src->l + 1);
    }

    static void cpy_kseq(kseq_t *dst, const kseq_t *src)
    {
        cpy_kstr(&dst->name, &src->name);
        cpy_kstr(&dst->seq,  &src->seq);
        cpy_kstr(&dst->qual, &src->qual);
        cpy_kstr(&dst->comment, &src->comment);
    }

    // int stk_sample(int argc, char *argv[])
    int stk_sample(const string & inputFileName, 
                   double frac,
                   const long int & seedNum, 
                   const string & mode
    )
    {
        int c, twopass = 0;
        uint64_t i, num = 0, n_seqs = 0;
        gzFile fp;
        kseq_t *seq;
        krand_t *kr = 0;

        // 种子
        kr = kr_srand(seedNum);

        if (mode == "2")
        {
            twopass = 1;
        }

        // 判断frac是否大于1，大于1用number，小于1用比例
        if (frac >= 1.0) num = (uint64_t)(frac + .499), frac = 0.;
        else if (twopass) {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                 << "when sampling a fraction, option -2 is ignored.\n";
            twopass = 0;
        }

        if (!twopass) { // the streaming version
            kseq_t *buf = 0;
            if (num > 0) buf = (kseq_t *)calloc(num, sizeof(kseq_t));
            if (num > 0 && buf == NULL) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                     << "Could not allocate enough memory for " << num << " sequences. Exiting...\n";
                free(kr);
                return 1;
            }

            // 打开文件
            fp = gzopen(inputFileName.c_str(), "r");
            if (fp == 0) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                     << "failed to open the input file.\n";
                return 1;
            }
            seq = kseq_init(fp);
            n_seqs = 0;
            while (kseq_read(seq) >= 0) {
                double r = kr_drand(kr);
                ++n_seqs;
                if (num) {
                    uint64_t y = n_seqs - 1 < num? n_seqs - 1 : (uint64_t)(r * n_seqs);
                    if (y < num) cpy_kseq(&buf[y], seq);
                } else if (r < frac) stk_printseq(stdout, seq, UINT_MAX);
            }
            for (i = 0; i < num; ++i) {
                kseq_t *p = &buf[i];
                if (p->seq.l) stk_printseq(stdout, p, UINT_MAX);
                free(p->seq.s); free(p->qual.s); free(p->name.s);
            }
            if (buf != NULL) free(buf);
        } else {
            uint64_t *buf;
            khash_t(64) *hash;
            int absent;

            // 1st pass
            buf = (uint64_t *)malloc(num * 8);
            for (i = 0; i < num; ++i) buf[i] = UINT64_MAX;
            fp = gzopen(inputFileName.c_str(), "r");
            if (fp == 0) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                     << "failed to open the input file.\n";
                fprintf(stderr, "[E::%s] failed to open the input file.\n", __func__);
                return 1;
            }
            seq = kseq_init(fp);
            n_seqs = 0;
            while (kseq_read(seq) >= 0) {
                double r = kr_drand(kr);
                uint64_t y;
                ++n_seqs;
                y = n_seqs - 1 < num? n_seqs - 1 : (uint64_t)(r * n_seqs);
                if (y < num) buf[y] = n_seqs;
            }
            kseq_destroy(seq);
            gzclose(fp);
            hash = kh_init(64);
            for (i = 0; i < num; ++i) kh_put(64, hash, buf[i], &absent);
            free(buf);
            // 2nd pass
            fp = gzopen(inputFileName.c_str(), "r");
            seq = kseq_init(fp);
            n_seqs = 0;
            while (kseq_read(seq) >= 0)
                if (kh_get(64, hash, ++n_seqs) != kh_end(hash))
                    stk_printseq(stdout, seq, UINT_MAX);
            kh_destroy(64, hash);
        }

        kseq_destroy(seq);
        gzclose(fp);
        free(kr);
        return 0;
    }
}

#endif