#include "functions.h"

int main(int argc, char** argv) {
//	_controlfp(~(_MCW_EM & (~_EM_INEXACT) & (~_EM_UNDERFLOW)), _MCW_EM);
	init();

    t_glob = 0.0;
	step = 0;
    save_csv(step);
	while (t_glob < TMAX) {
        t_glob += TAU;
		step++;

        copy_to_half();

        calc_integral_rhs();
        calc_new();
        calc_limiters();

        calc_integral_rhs();
        calc_new();

        calc_new_1();

        calc_limiters();
        calc_integral_rhs();
        calc_new();

        calc_new_2();
        calc_limiters();

        if (step % STEP_OUTPUT == 0) {
			save_csv(step);
			printf("STEP: %d\t\tTIME: %25.16e\n", step, t_glob);
		}

	}

    printf("ERROR_L2: %e\n", calc_err_L2());

	mem_deallocate();
	return 0;
}
