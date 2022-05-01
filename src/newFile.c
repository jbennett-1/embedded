// new file testing memory
#include <stdint.h>
#include "arm_math.h"
#define ARMCM7

struct testing_args{
        arm_matrix_instance_f32* input_buffer;
        
	float32_t* data;
	uint16_t vec_len;
        uint16_t vec_num;
	
};

void initialize(struct testing_args *args, uint16_t vec_len, uint16_t vec_num, float32_t* data){
        args->vec_len=vec_len;
        args->vec_num=vec_num;
        args->data=data;
}

void test_functions(struct testing_args *args){
	arm_mat_init_f32(args->input_buffer, args->vec_len, args->vec_num, args->data);
}

int main(){
        while(1);
	struct testing_args *args;
	test_functions(args);
}

