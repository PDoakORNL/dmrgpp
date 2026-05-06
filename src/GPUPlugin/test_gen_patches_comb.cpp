#include "analysis.h"
#include "dmrg_types.h"
#include <stdio.h>
#include <stdlib.h>

/*
 ---------------------------------------
 simple program to test gen_patches_comb()
 ---------------------------------------
*/

IntegerType main()
{

	const IntegerType left_size         = 8;
	const IntegerType right_size        = 10;
	const IntegerType target_up         = (left_size + right_size) / 2;
	const IntegerType target_down       = target_up;
	const IntegerType keep_left_states  = 5000;
	const IntegerType keep_right_states = 4 * keep_left_states;
	const IntegerType max_patches       = (1 + target_up) * (1 + target_down);

	VectorSizeType left_patch_size_(max_patches + 1);
	VectorSizeType left_patch_up_(max_patches + 1);
	VectorSizeType left_patch_down_(max_patches + 1);

	VectorSizeType right_patch_size_(max_patches + 1);
	VectorSizeType right_patch_up_(max_patches + 1);
	VectorSizeType right_patch_down_(max_patches + 1);

#define left_patch_size(i) left_patch_size_[(i) - 1]
#define right_patch_size(i) right_patch_size_[(i) - 1]

	IntegerType npatches = gen_patches_comb(left_size,
	                                        right_size,
	                                        target_up,
	                                        target_down,
	                                        keep_left_states,
	                                        keep_right_states,

	                                        left_patch_size_,
	                                        right_patch_size_,

	                                        left_patch_up_,
	                                        left_patch_down_,

	                                        right_patch_up_,
	                                        right_patch_down_);

	printf("left_size=%d, right_size=%d, target_up=%d, target_down=%d\n",
	       left_size,
	       right_size,
	       target_up,
	       target_down);
	printf("keep_left_states=%d, keep_right_states=%d\n", keep_left_states, keep_right_states);

	printf("npatches=%d\n", npatches);
	{
		for (IntegerType ipatch = 1; ipatch <= npatches; ipatch++) {
			printf("ipatch=%d, left_patch_size=%lu, right_patch_size=%lu\n",
			       ipatch,
			       left_patch_size(ipatch),
			       right_patch_size(ipatch));
		};
	}

	exit(0);
	return (0);
}
