/*
  CSE 5242 Project 2, Fall 2024
  Author: Benji Ofori
          Gabe Richner
          Johnathan Ross
          Test
  See class project handout for more extensive documentation.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <asm/unistd.h>
#include <immintrin.h>

/* uncomment out the following DEBUG line for debug info, for experiment comment the DEBUG line  */
// #define DEBUG

/* compare two int64_t values - for use with qsort */
static int compare(const void *p1, const void *p2)
{
  int a,b;
  a = *(int64_t *)p1;
  b = *(int64_t *)p2;
  if (a<b) return -1;
  if (a==b) return 0;
  return 1;
}

/* initialize searches and data - data is sorted and searches is a random permutation of data */
int init(int64_t* data, int64_t* searches, int count) {
  for(int64_t i = 0; i < count; i++){
    searches[i] = rand();
    data[i] = searches[i] + 1;
  }
  qsort(data,count,sizeof(int64_t),compare);
}

/* initialize outer probes of band join */
int band_init(int64_t* outer, int64_t size)
{
  for(int64_t i = 0; i < size; i++){
    outer[i] = rand();
  }
}

inline int64_t simple_binary_search(int64_t* data, int64_t size, int64_t target)
{
  int64_t left=0;
  int64_t right=size;
  int64_t mid;

  while(left<=right) {
    mid = (left + right)/2;   /* ignore possibility of overflow of left+right */
    if (data[mid]==target) return mid;
    if (data[mid]<target) left=mid+1;
    else right = mid-1;
  }
  return -1; /* no match */
}

inline int64_t low_bin_search(int64_t* data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does only one comparison in the inner loop
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left=0;
  int64_t right=size;
  int64_t mid;

  while(left<right) {
    mid = (left + right)/2;   /* ignore possibility of overflow of left+right */
    if (data[mid]>=target)
      right=mid;
    else
      left=mid+1;
  }
  return right;
}

inline int64_t low_bin_nb_arithmetic(int64_t* data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using multiplication and addition to convert control dependencies
         to data dependencies
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left = 0;
  int64_t right = size;
  int64_t mid;

  while(left < right) {

    mid = (left + right) >> 1; // use a right shift instead of division to improve performance  
    // Check if the middle element is less then the target.
    int64_t lessThan = (data[mid] < target);
    left = left + lessThan * (mid + 1 - left); //update left if mid is less than the target.
    right = right - (1 - lessThan) * (right - mid); //updates right if mid is greater or equal to the target.

  }
  return right;
}

inline int64_t low_bin_nb_mask(int64_t* data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using bit masking operations to convert control dependencies
         to data dependencies
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left=0;
  int64_t right=size;
  int64_t mid;

  while(left < right) {
     mid = left + ((right - left) >> 1);
     //set up a mask used to update the right and left values
     int64_t conditionMask = (target > data[mid]) - 1; //if target > data[mid] then the mask will be 0xFFFFFFFFFFFFFFFF, otherwise set to 0
     left = (mid + 1) & conditionMask | left & ~conditionMask; //left updates when target > data[mid]
     right = mid & ~conditionMask | right & conditionMask; //right moves down
    
  }


  return right;
}

inline void low_bin_nb_4x(int64_t* data, int64_t size, int64_t* targets, int64_t* right)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using bit masking operations instead
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
     (d) does 4 searches at the same time in an interleaved fashion, so that an out-of-order processor
         can make progress even when some instructions are still waiting for their operands to be ready.

     Note that we're using the result array as the "right" elements in the search so no need for a return statement
  */

    //intialize the arrays for left and right
     // Initialize the arrays for left and right
    int64_t left[4] = {0, 0, 0, 0};
    int64_t mask[4];
    int64_t mid[4];

    for (int64_t i = 0; i < 4; i++) {
        right[i] = size; // Initialized to the size
    }

    int flag = 0;
    while (!flag) {
        flag = 1;
        for (int64_t i = 0; i < 4; i++) {
            if (right[i] > left[i]) {
                flag = 0;
                mid[i] = left[i] + ((right[i] - left[i]) / 2);
                mask[i] = (data[mid[i]] - targets[i]) >> 63;
                left[i] = (mid[i] + 1) & mask[i] | left[i] & ~mask[i];
                right[i] = mid[i] & ~mask[i] | right[i] & mask[i];
            }
        }
    }
}


/* The following union type is handy to output the contents of AVX512 data types */
union int8x4 {
  __m256i a;
  int64_t b[4];
};

void printavx(char* name, __m256i v) {
  union int8x4 n;

  n.a=v;
  printf("Value in %s is [%ld %ld %ld %ld ]\n",name,n.b[0],n.b[1],n.b[2],n.b[3]);
}

/*
 * Optinal for using AVX-512

  union int8x8 {
    __m512i a;
    int64_t b[8];
  };

  void printavx512(char* name, __m512i v) {
    union int8x4 n;

    n.a=v;
    printf("Value in %s is [%ld %ld %ld %ld %ld %ld %ld %ld ]\n",name,n.b[0],n.b[1],n.b[2],n.b[3]);
  }

 */


inline void low_bin_nb_simd(int64_t* data, int64_t size, __m256i target, __m256i* result)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using bit masking operations instead
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
     (d) does 4 searches at the same time using AVX2 intrinsics

     See https://www.intel.com/content/www/us/en/docs/cpp-compiler/developer-guide-reference/2021-8/intrinsics-for-avx2.html
     for documentation of the AVX512 intrinsics

     Note that we're using the result array as the "right" elements in the search, and that searchkey is being passed
     as an __m256i value rather than via a pointer.
  */


  //initalize variables right set to 0 and left set to size.
  __m256i left = _mm256_set1_epi64x(0);
  __m256i right = _mm256_set1_epi64x(size);
  __m256i one = _mm256_set1_epi64x(1);
  //updates the indicies by creating a constant vector of 1.

  //loop will compare until left == right
 while (!_mm256_testc_si256(_mm256_cmpeq_epi64(left, right), _mm256_set1_epi64x(-1))) {
        __m256i mid = _mm256_add_epi64(left, _mm256_srli_epi64(_mm256_sub_epi64(right, left), 1));
        __m256i mid_values = _mm256_i64gather_epi64((const long long int*)data, mid, 8);
        __m256i cmp = _mm256_cmpgt_epi64(target, mid_values);

        left = _mm256_blendv_epi8(left, _mm256_add_epi64(mid, one), cmp);
        right = _mm256_blendv_epi8(mid, right, cmp);
    }
  
  *result = right;


}

void bulk_bin_search(int64_t* data, int64_t size, int64_t* searchkeys, int64_t numsearches, int64_t* results, int repeats)
{
  for(int j=0; j<repeats; j++) {
    /* Function to test a large number of binary searches

       we might need repeats>1 to make sure the events we're measuring are not dominated by various
       overheads, particularly for small values of size and/or numsearches

       we assume that we want exactly "size" searches, where "size" is the length if the searchkeys array
     */
    for(int64_t i=0;i<numsearches; i++) {
#ifdef DEBUG
      printf("Searching for %ld...\n",searchkeys[i]);
#endif

      // Uncomment one of the following to measure it
      results[i] = low_bin_search(data,size,searchkeys[i]);
      //results[i] = low_bin_nb_arithmetic(data,size,searchkeys[i]);
      //results[i] = low_bin_nb_mask(data,size,searchkeys[i]);

#ifdef DEBUG
      printf("Result is %ld\n",results[i]);
#endif
    }
  }
}

void bulk_bin_search_4x(int64_t* data, int64_t size, int64_t* searchkeys, int64_t numsearches, int64_t* results, int repeats)
{
  register __m256i searchkey_4x;

  for(int j=0; j<repeats; j++) {
    /* Function to test a large number of binary searches using one of the 8x routines

       we might need repeats>1 to make sure the events we're measuring are not dominated by various
       overheads, particularly for small values of size and/or numsearches

       we assume that we want exactly "size" searches, where "size" is the length if the searchkeys array
     */
    int64_t extras = numsearches % 4;
    for(int64_t i=0;i<numsearches-extras; i+=4) {
#ifdef DEBUG
      printf("Searching for %ld %ld %ld %ld  ...\n",
	     searchkeys[i],searchkeys[i+1],searchkeys[i+2],searchkeys[i+3]);
#endif

      // Uncomment one of the following depending on which routine you want to profile

      // Algorithm A
       low_bin_nb_4x(data,size,&searchkeys[i],&results[i]);

      // Algorithm B
      // searchkey_4x = _mm256_loadu_si256((__m256i *)&searchkeys[i]);
      // low_bin_nb_simd(data,size,searchkey_4x,(__m256i *)&results[i]);

#ifdef DEBUG
      printf("Result is %ld %ld %ld %ld  ...\n",
	     results[i],results[i+1],results[i+2],results[i+3]);
#endif
    }
    /* a little bit more work if numsearches is not a multiple of 8 */
    for(int64_t i=numsearches-extras;i<numsearches; i++) {

      results[i] = low_bin_nb_mask(data,size,searchkeys[i]);

    }

  }
}


int64_t band_join(int64_t* inner, int64_t inner_size, int64_t* outer, int64_t outer_size, int64_t* inner_results, int64_t* outer_results, int64_t result_size, int64_t bound)
{
  /* In a band join we want matches within a range of values.  If p is the probe value from the outer table, then all
     reccords in the inner table with a key in the range [p-bound,p+bound] inclusive should be part of the result.

     Results are returned via two arrays. outer_results stores the index of the outer table row that matches, and
     inner_results stores the index of the inner table row that matches.  result_size tells you the size of the
     output array that has been allocated. You should make sure that you don't exceed this size.  If there are
     more results than can fit in the result arrays, then return early with just a prefix of the results in the result
     arrays. The return value of the function should be the number of output results.

  */
    int64_t result = 0;
    int flag = 0; // Flag used to stop the loop when the result arrays are full
    int64_t i = 0;

    while (i < outer_size && !flag) {
        int64_t target = outer[i];
        int64_t targetValues[4] = {target - bound, target - bound, target + bound, target + bound};
        int64_t right[4];

        low_bin_nb_4x(inner, inner_size, targetValues, right); // Finds the range

        // Collect the results in the range
        int64_t j = right[0];
        while (j < inner_size && inner[j] <= target + bound && !flag) {
            if (target - bound <= inner[j]) {
                if (result < result_size) {
                    inner_results[result] = j;
                    outer_results[result] = i;
                    result++;
                } else {
                    flag = 1; // Once the array is full, set the flag to end the loop
                }
            }
            j++;
        }

        i++;
    }

    // Handle remaining records if outer_size is not a multiple of 4
    for (; i < outer_size; ++i) {
        int64_t remaining_target = outer[i];
        int64_t remainRight = low_bin_nb_mask(inner, inner_size, remaining_target - bound);
        while (remainRight < inner_size && inner[remainRight] <= remaining_target + bound && !flag) {
            if (remaining_target - bound <= inner[remainRight]) {
                if (result < result_size) {
                    inner_results[result] = remainRight;
                    outer_results[result] = i;
                    result++;
                } else {
                    flag = 1; // Set the flag to indicate the result array is full
                }
            }
            remainRight++;
        }
    }

    return result;
}

int64_t band_join_simd(int64_t* inner, int64_t inner_size, int64_t* outer, int64_t outer_size, int64_t* inner_results, int64_t* outer_results, int64_t result_size, int64_t bound)
{
  /* In a band join we want matches within a range of values.  If p is the probe value from the outer table, then all
     reccords in the inner table with a key in the range [p-bound,p+bound] inclusive should be part of the result.

     Results are returned via two arrays. outer_results stores the index of the outer table row that matches, and
     inner_results stores the index of the inner table row that matches.  result_size tells you the size of the
     output array that has been allocated. You should make sure that you don't exceed this size.  If there are
     more results than can fit in the result arrays, then return early with just a prefix of the results in the result
     arrays. The return value of the function should be the number of output results.

     To do the binary search, you could use the low_bin_nb_simd you just implemented to search for the lower bounds in parallel

     Once you've found the lower bounds, do the following for each of the 4 search keys in turn:
        scan along the sorted inner array, generating outputs for each match, and making sure not to exceed the output array bounds.

     This inner scanning code does not have to use SIMD.
  */

     int64_t result = 0;
    __m256i boundVec = _mm256_set1_epi64x(bound);
    int64_t i = 0;

    // Process main chunks with SIMD
    while (i + 4 <= outer_size) {
        __m256i outerVec = _mm256_loadu_si256((__m256i*)&outer[i]);
        __m256i lowerBound = _mm256_sub_epi64(outerVec, boundVec);
        __m256i upperBound = _mm256_add_epi64(outerVec, boundVec);

        int64_t lowerBounds[4];
        int64_t upperBounds[4];
        _mm256_storeu_si256((__m256i*)lowerBounds, lowerBound);
        _mm256_storeu_si256((__m256i*)upperBounds, upperBound);

        for (int j = 0; j < 4; ++j) {
            int64_t lower = lowerBounds[j];
            int64_t upper = upperBounds[j];
            int64_t outerIndex = i + j;

            // Use i64gather to load values from the inner array
            __m256i indices = _mm256_set_epi64x(3, 2, 1, 0);
            __m256i innerValues = _mm256_i64gather_epi64((const long long int*)inner, indices, sizeof(int64_t));

            int64_t innerVals[4];
            _mm256_storeu_si256((__m256i*)innerVals, innerValues);

            for (int l = 0; l < 4; ++l) {
                if (innerVals[l] >= lower && innerVals[l] <= upper) {
                    if (result < result_size) {
                        outer_results[result] = outerIndex;
                        inner_results[result] = l;
                        result++;
                    } else {
                        return result;
                    }
                }
            }
        }
        i += 4;
    }

    // Handle remaining outer records without SIMD
    for (; i < outer_size; ++i) {
        int64_t lower = outer[i] - bound;
        int64_t upper = outer[i] + bound;

        for (int64_t k = 0; k < inner_size; ++k) {
            if (inner[k] >= lower && inner[k] <= upper) {
                if (result < result_size) {
                    outer_results[result] = i;
                    inner_results[result] = k;
                    result++;
                } else {
                    return result;
                }
            }
        }
    }

    return result;
}

int
main(int argc, char *argv[])
{
	 long long counter;
	 int64_t arraysize, outer_size, result_size;
	 int64_t bound;
	 int64_t *data, *queries, *results;
	 int ret;
	 struct timeval before, after;
	 int repeats;
	 int64_t total_results;

	 // band-join arrays
	 int64_t *outer, *outer_results, *inner_results;


	 if (argc >= 5)
	   {
	     arraysize = atoi(argv[1]);
	     outer_size = atoi(argv[2]);
	     result_size = atoi(argv[3]);
	     bound = atoi(argv[4]);
	   }
	 else
	   {
	     fprintf(stderr, "Usage: db5242 inner_size outer_size result_size bound <repeats>\n");
	     exit(EXIT_FAILURE);
	   }

	 if (argc >= 6)
	   repeats = atoi(argv[5]);
	 else
	   {
	     repeats=1;
	   }

	 /* allocate the array and the queries for searching */
	 ret=posix_memalign((void**) &data,64,arraysize*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }
	 ret=posix_memalign((void**) &queries,64,arraysize*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }
	 ret=posix_memalign((void**) &results,64,arraysize*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }

	 /* allocate the outer array and output arrays for band-join */
	 ret=posix_memalign((void**) &outer,64,outer_size*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }
	 ret=posix_memalign((void**) &outer_results,64,result_size*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }
	 ret=posix_memalign((void**) &inner_results,64,result_size*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }


	   /* code to initialize data structures goes here so that it is not included in the timing measurement */
	   init(data,queries,arraysize);
	   band_init(outer,outer_size);

#ifdef DEBUG
	   /* show the arrays */
	   printf("data: ");
	   for(int64_t i=0;i<arraysize;i++) printf("%ld ",data[i]);
	   printf("\n");
	   printf("queries: ");
	   for(int64_t i=0;i<arraysize;i++) printf("%ld ",queries[i]);
	   printf("\n");
	   printf("outer: ");
	   for(int64_t i=0;i<outer_size;i++) printf("%ld ",outer[i]);
	   printf("\n");
#endif


	   /* now measure... */

	   gettimeofday(&before,NULL);

	   /* the code that you want to measure goes here; make a function call */
	   bulk_bin_search(data,arraysize,queries,arraysize,results, repeats);

	   gettimeofday(&after,NULL);
	   printf("Time in bulk_bin_search loop is %ld microseconds or %f microseconds per search\n", (after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec), 1.0*((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/arraysize/repeats);



	   gettimeofday(&before,NULL);

	   /* the code that you want to measure goes here; make a function call */
	   bulk_bin_search_4x(data,arraysize,queries,arraysize,results, repeats);

	   gettimeofday(&after,NULL);
	   printf("Time in bulk_bin_search_4x loop is %ld microseconds or %f microseconds per search\n", (after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec), 1.0*((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/arraysize/repeats);


	   gettimeofday(&before,NULL);

	   /* the code that you want to measure goes here; make a function call */
	   total_results=band_join(data, arraysize, outer, outer_size, inner_results, outer_results, result_size, bound);

	   gettimeofday(&after,NULL);
	   printf("Band join result size is %ld with an average of %f matches per output record\n",total_results, 1.0*total_results/(1.0+outer_results[total_results-1]));
	   printf("Time in band_join loop is %ld microseconds or %f microseconds per outer record\n", (after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec), 1.0*((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/outer_size);

//#ifdef DEBUG
	   /* show the band_join results */
	   printf("band_join results: ");
	   for(int64_t i=0;i<total_results;i++) printf("(%ld,%ld) ",outer_results[i],inner_results[i]);
	   printf("\n");
//#endif

}

